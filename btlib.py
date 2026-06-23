"""
btlib.py -- Exact arithmetic for the Bruhat-Tits tree of U_3(Z[1/3, omega]).

Companion library for "Buildings for Synthesis with Clifford+R"
(Deaconu, Gargava, Kalra, Mosca, Yard).

Everything is EXACT: elements of F = Q(omega) are represented as
(a + b*omega)/d with a, b, d integers.  No floating point anywhere.

Conventions (matching the paper):
  omega = exp(2*pi*i/3),   omega^2 = -1 - omega
  chi   = 1 - omega,       chi * conj(chi) = 3,    conj(chi) = -omega^2 * chi
  pi    = chi * O_F  (the ramified prime above 3); residue field F_3,
          with residue map  a + b*omega  |->  (a + b) mod 3   (omega == 1 mod chi)
  v_pi  = chi-adic valuation ("sde" in synthesis literature); v_pi(3) = 2.
  <x,y> = sum_i x_i * conj(y_i)        (the Hermitian form)
  Gram of a basis matrix A is  A* A    (entrywise conj-transpose times A);
          a lattice A*O^3 is self-dual  iff  A* A in GL_3(O_pi)
          (entries v_pi >= 0 and det a unit).   [Lemma le:gramm_matrix]
  l(g)  = -2 * min_{ij} v_pi(g_ij)     [Eq. defi_of_l]
  d~(g,h) = (l(g^-1 h) + l(h^-1 g))/2  = graph distance between pure vertices.

Vertex conventions for the tree:
  * pure vertex:        lattice class {pi^i L},  L self-dual, normalized so
                        v_pi(det basis) == 0.
  * alternating vertex: class {pi^i L, pi^i L#} with L# < L (index 9);
                        we store the BIG representative L, normalized so
                        v_pi(det basis) == -1.   (Then det L# has v = +1.)
Canonical vertex keys come from the column Hermite normal form over O_pi.
"""

from __future__ import annotations
from fractions import Fraction
from math import gcd, inf
import itertools

# ----------------------------------------------------------------------------
# Z[omega]
# ----------------------------------------------------------------------------

class Zw:
    """a + b*omega with a, b integers."""
    __slots__ = ("a", "b")

    def __init__(self, a, b=0):
        self.a = int(a)
        self.b = int(b)

    def __add__(self, o):  return Zw(self.a + o.a, self.b + o.b)
    def __sub__(self, o):  return Zw(self.a - o.a, self.b - o.b)
    def __neg__(self):     return Zw(-self.a, -self.b)

    def __mul__(self, o):
        if isinstance(o, int):
            return Zw(self.a * o, self.b * o)
        # (a+bw)(c+dw) = ac + (ad+bc) w + bd w^2 ;  w^2 = -1-w
        a, b, c, d = self.a, self.b, o.a, o.b
        return Zw(a * c - b * d, a * d + b * c - b * d)

    __rmul__ = __mul__

    def conj(self):
        # conj(a + b w) = a + b w^2 = (a - b) - b w
        return Zw(self.a - self.b, -self.b)

    def norm(self):
        """N(a+bw) = a^2 - a b + b^2  (a non-negative integer)."""
        return self.a * self.a - self.a * self.b + self.b * self.b

    def is_zero(self):
        return self.a == 0 and self.b == 0

    def res3(self):
        """Residue in O_F/chi = F_3:  a + b mod 3  (omega == 1)."""
        return (self.a + self.b) % 3

    def div_chi(self):
        """Exact division by chi = 1-omega.  x/chi = x*(2+omega)/3.
        Requires res3 == 0."""
        a, b = self.a, self.b
        # (a+bw)(2+w) = (2a - b) + (a + b) w
        na, nb = 2 * a - b, a + b
        assert na % 3 == 0 and nb % 3 == 0, "not divisible by chi"
        return Zw(na // 3, nb // 3)

    def v_chi(self):
        """chi-adic valuation; inf for 0."""
        if self.is_zero():
            return inf
        x, v = self, 0
        while x.res3() == 0:
            x = x.div_chi()
            v += 1
        return v

    def __eq__(self, o):
        return isinstance(o, Zw) and self.a == o.a and self.b == o.b

    def __hash__(self):
        return hash((self.a, self.b))

    def __repr__(self):
        return f"({self.a}{self.b:+d}w)"


# ----------------------------------------------------------------------------
# Q(omega):  (a + b*omega)/den
# ----------------------------------------------------------------------------

def _v3(n):
    """3-adic valuation of a positive integer."""
    v = 0
    while n % 3 == 0:
        n //= 3
        v += 1
    return v


class Fw:
    """Element of Q(omega): num/den with num in Z[omega], den a positive int."""
    __slots__ = ("num", "den")

    def __init__(self, num, den=1, _reduce=True):
        if isinstance(num, int):
            num = Zw(num, 0)
        den = int(den)
        if den < 0:
            num, den = -num, -den
        assert den > 0
        if _reduce and den != 1:
            g = gcd(gcd(abs(num.a), abs(num.b)), den)
            if g > 1:
                num = Zw(num.a // g, num.b // g)
                den //= g
        self.num = num
        self.den = den

    # -- arithmetic --------------------------------------------------------
    def __add__(self, o):
        o = _coerce(o)
        return Fw(self.num * o.den + o.num * self.den, self.den * o.den)

    __radd__ = __add__

    def __sub__(self, o):
        o = _coerce(o)
        return Fw(self.num * o.den - o.num * self.den, self.den * o.den)

    def __rsub__(self, o):
        return _coerce(o) - self

    def __neg__(self):
        return Fw(-self.num, self.den, _reduce=False)

    def __mul__(self, o):
        o = _coerce(o)
        return Fw(self.num * o.num, self.den * o.den)

    __rmul__ = __mul__

    def inv(self):
        """1/x = den * conj(num) / N(num)."""
        n = self.num.norm()
        assert n != 0, "division by zero"
        return Fw(self.num.conj() * self.den, n)

    def __truediv__(self, o):
        return self * _coerce(o).inv()

    def conj(self):
        return Fw(self.num.conj(), self.den, _reduce=False)

    # -- predicates / valuation --------------------------------------------
    def is_zero(self):
        return self.num.is_zero()

    def v_pi(self):
        """pi-adic valuation: v_chi(num) - 2*v_3(den).
        (Primes other than 3 in den are units in O_pi.)"""
        if self.is_zero():
            return inf
        return self.num.v_chi() - 2 * _v3(self.den)

    def is_integral(self):
        """Is x in O_pi?  (v_pi >= 0)"""
        return self.v_pi() >= 0

    def is_unit(self):
        return self.v_pi() == 0

    def res3(self):
        """Residue in O_pi/pi = F_3 (requires v_pi >= 0).

        Write den = 3^m * d with gcd(d,3)=1.  Then
        1/3 = -omega * chi^-2 (since 3 = -omega^2 chi^2), so
        x = num * (-omega)^m * chi^(-2m) / d, and we divide num by chi^(2m)
        exactly, then take residues; (-omega) has residue -1 = 2.
        """
        assert self.v_pi() >= 0, "res3 of a non-integral element"
        m = _v3(self.den)
        d = self.den // (3 ** m)
        z = self.num
        for _ in range(2 * m):
            z = z.div_chi()
        r = z.res3() * pow(2, m, 3) * pow(d, -1, 3)
        return r % 3

    # -- misc ----------------------------------------------------------------
    def key(self):
        return (self.num.a, self.num.b, self.den)

    def __eq__(self, o):
        o = _coerce(o)
        return (self.num * o.den) == (o.num * self.den)

    def __hash__(self):
        return hash(self.key())

    def __repr__(self):
        if self.den == 1:
            return repr(self.num)
        return f"{self.num!r}/{self.den}"


def _coerce(x):
    if isinstance(x, Fw):
        return x
    if isinstance(x, Zw):
        return Fw(x, 1)
    if isinstance(x, int):
        return Fw(Zw(x, 0), 1)
    if isinstance(x, Fraction):
        return Fw(Zw(x.numerator, 0), x.denominator)
    raise TypeError(f"cannot coerce {x!r}")


ZERO  = Fw(0)
ONE   = Fw(1)
OMEGA = Fw(Zw(0, 1))
CHI   = Fw(Zw(1, -1))            # 1 - omega
CHI_INV = CHI.inv()              # (2+omega)/3
I_SQRT3 = Fw(Zw(1, 2), 3)        # i/sqrt(3) = (omega - omega^2)/3 = (1+2w)/3


def chi_pow(k):
    """chi^k as an Fw, any integer k."""
    if k >= 0:
        r = ONE
        for _ in range(k):
            r = r * CHI
        return r
    r = ONE
    for _ in range(-k):
        r = r * CHI_INV
    return r


# ----------------------------------------------------------------------------
# Matrices over Q(omega)
# ----------------------------------------------------------------------------

class Mat:
    """Dense matrix over Fw (any square size; the paper uses n=3)."""
    __slots__ = ("m", "n")

    def __init__(self, rows):
        self.m = [[_coerce(x) for x in row] for row in rows]
        self.n = len(self.m)
        assert all(len(r) == self.n for r in self.m), "square matrices only"

    @staticmethod
    def identity(n=3):
        return Mat([[ONE if i == j else ZERO for j in range(n)]
                    for i in range(n)])

    @staticmethod
    def diag(*entries):
        n = len(entries)
        return Mat([[_coerce(entries[i]) if i == j else ZERO
                     for j in range(n)] for i in range(n)])

    def __getitem__(self, ij):
        i, j = ij
        return self.m[i][j]

    def __mul__(self, o):
        if isinstance(o, Mat):
            assert self.n == o.n
            n = self.n
            return Mat([[sum((self.m[i][k] * o.m[k][j] for k in range(n)),
                             ZERO) for j in range(n)] for i in range(n)])
        s = _coerce(o)
        return Mat([[x * s for x in row] for row in self.m])

    def __rmul__(self, s):
        return self * s

    def __add__(self, o):
        return Mat([[self.m[i][j] + o.m[i][j] for j in range(self.n)]
                    for i in range(self.n)])

    def __sub__(self, o):
        return Mat([[self.m[i][j] - o.m[i][j] for j in range(self.n)]
                    for i in range(self.n)])

    def __neg__(self):
        return Mat([[-x for x in row] for row in self.m])

    def conjT(self):
        """Conjugate transpose A*."""
        return Mat([[self.m[j][i].conj() for j in range(self.n)]
                    for i in range(self.n)])

    def T(self):
        return Mat([[self.m[j][i] for j in range(self.n)]
                    for i in range(self.n)])

    def det(self):
        n = self.n
        if n == 1:
            return self.m[0][0]
        if n == 2:
            return self.m[0][0] * self.m[1][1] - self.m[0][1] * self.m[1][0]
        # Laplace along first row (n is small here)
        total = ZERO
        for j in range(n):
            if self.m[0][j].is_zero():
                continue
            minor = Mat([row[:j] + row[j + 1:] for row in self.m[1:]])
            term = self.m[0][j] * minor.det()
            total = total + (term if j % 2 == 0 else -term)
        return total

    def inv(self):
        n = self.n
        d = self.det()
        assert not d.is_zero(), "singular matrix"
        if n == 1:
            return Mat([[d.inv()]])
        cof = [[None] * n for _ in range(n)]
        for i in range(n):
            for j in range(n):
                minor = Mat([row[:j] + row[j + 1:]
                             for k, row in enumerate(self.m) if k != i])
                c = minor.det()
                cof[i][j] = c if (i + j) % 2 == 0 else -c
        dinv = d.inv()
        return Mat([[cof[j][i] * dinv for j in range(n)] for i in range(n)])

    def col(self, j):
        return [self.m[i][j] for i in range(self.n)]

    def with_col(self, j, vec):
        rows = [row[:] for row in self.m]
        for i in range(self.n):
            rows[i][j] = _coerce(vec[i])
        return Mat(rows)

    def __eq__(self, o):
        return (isinstance(o, Mat) and self.n == o.n and
                all(self.m[i][j] == o.m[i][j]
                    for i in range(self.n) for j in range(self.n)))

    def key(self):
        return tuple(x.key() for row in self.m for x in row)

    def __hash__(self):
        return hash(self.key())

    def __repr__(self):
        return "Mat(" + ", ".join(repr(r) for r in self.m) + ")"

    # -- O_pi / unitarity predicates ---------------------------------------
    def min_v(self):
        return min(x.v_pi() for row in self.m for x in row)

    def is_O(self):
        """All entries in O_pi."""
        return self.min_v() >= 0

    def in_GL_O(self):
        """In GL_n(O_pi): integral entries and unit determinant."""
        return self.is_O() and self.det().v_pi() == 0

    def is_unitary(self):
        return self.conjT() * self == Mat.identity(self.n)


def ell(g: Mat):
    """l(g) = -2 min v_pi(entries)   [Eq. defi_of_l]."""
    return -2 * g.min_v()


def d_tilde(g: Mat, h: Mat):
    """d~(g,h) = (l(g^-1 h) + l(h^-1 g)) / 2."""
    gi_h = g.inv() * h
    return (ell(gi_h) + ell(gi_h.inv())) // 2


def sde(g: Mat):
    """Smallest denominator exponent = -min v_pi = l(g)/2 (for g in A)."""
    return ell(g) // 2


# ----------------------------------------------------------------------------
# Gates and the finite group U_3(O_F)
# ----------------------------------------------------------------------------

W  = OMEGA
W2 = OMEGA * OMEGA

H_GATE = Mat([[I_SQRT3, I_SQRT3, I_SQRT3],
              [I_SQRT3, I_SQRT3 * W, I_SQRT3 * W2],
              [I_SQRT3, I_SQRT3 * W2, I_SQRT3 * W]])

S_GATE = Mat.diag(1, W, 1)
R_GATE = Mat.diag(1, 1, -1)

_UNITS = [ONE, W, W2, -ONE, -W, -W2]      # the six units of Z[omega]


def monomial_matrices():
    """All 1296 elements of U_3(O_F): permutation x unit-diagonal."""
    out = []
    for perm in itertools.permutations(range(3)):
        for us in itertools.product(_UNITS, repeat=3):
            rows = [[ZERO] * 3 for _ in range(3)]
            for i in range(3):
                rows[i][perm[i]] = us[i]
            out.append(Mat(rows))
    return out


def in_Gamma(g: Mat):
    """Is g in U_3(Z[omega, 1/3])?  (unitary, entries with 3-power denom)"""
    if not g.is_unitary():
        return False
    for row in g.m:
        for x in row:
            d = x.den
            while d % 3 == 0:
                d //= 3
            if d != 1:
                return False
    return True


# ----------------------------------------------------------------------------
# Lattices  (Lambda = g * O_pi^3, g an invertible basis matrix)
# ----------------------------------------------------------------------------

def lattice_contains(g: Mat, h: Mat):
    """Does g*O^3 contain h*O^3?  iff g^-1 h has integral entries."""
    return (g.inv() * h).is_O()


def lattice_eq(g: Mat, h: Mat):
    """g*O^3 == h*O^3  iff  g^-1 h in GL(O)."""
    return (g.inv() * h).in_GL_O()


def dual_basis(g: Mat):
    """Basis of the dual lattice:  (g*)^-1."""
    return g.conjT().inv()


def gram(g: Mat):
    return g.conjT() * g


def is_self_dual(g: Mat):
    """Lemma le:gramm_matrix."""
    return gram(g).in_GL_O()


def lattice_index_log3(g: Mat, h: Mat):
    """log_3 [gO^3 : hO^3]  (requires containment) = v_pi(det(g^-1 h))."""
    q = g.inv() * h
    assert q.is_O(), "not contained"
    return q.det().v_pi()


# -- canonical form (column HNF over the DVR O_pi) ---------------------------

def _reduce_mod_chi_k(x: Fw, k):
    """Canonical representative of x mod chi^k * O_pi:
    r = sum_{i=v}^{k-1} d_i chi^i with digits d_i in {0,1,2}."""
    r = ZERO
    y = x
    while True:
        v = y.v_pi()
        if v >= k:
            break
        d = (y * chi_pow(-v)).res3()
        assert d != 0
        t = chi_pow(v) * d
        r = r + t
        y = y - t
    return r


def hnf_local(g: Mat):
    """Column Hermite normal form over O_pi:  returns g' = g*U,
    U in GL_3(O_pi), with g' lower-triangular, diagonal chi^{a_i}, and
    entries below each diagonal reduced mod that row's pivot.
    This is a canonical representative of the coset g*GL_3(O_pi),
    i.e. of the lattice g*O^3."""
    n = g.n
    work = Mat([row[:] for row in g.m])

    for r in range(n):
        # pivot: column (>= r) whose entry in row r has minimal valuation
        best, best_v = None, None
        for j in range(r, n):
            v = work.m[r][j].v_pi()
            if v != inf and (best_v is None or v < best_v):
                best, best_v = j, v
        assert best is not None, "singular matrix in hnf_local"
        if best != r:
            for i in range(n):
                work.m[i][r], work.m[i][best] = work.m[i][best], work.m[i][r]
        # normalize pivot to exactly chi^a  (multiply column by unit)
        a = best_v
        u = chi_pow(a) / work.m[r][r]          # a unit of O_pi
        for i in range(n):
            work.m[i][r] = work.m[i][r] * u
        # eliminate row r in later columns
        for j in range(r + 1, n):
            c = work.m[r][j] * chi_pow(-a)     # in O_pi
            if c.is_zero():
                continue
            for i in range(n):
                work.m[i][j] = work.m[i][j] - c * work.m[i][r]

    # reduce below-diagonal entries:  entry (i, j), j < i, mod chi^{a_i}
    for i in range(1, n):
        a_i = work.m[i][i].v_pi()
        for j in range(i):
            x = work.m[i][j]
            r_x = _reduce_mod_chi_k(x, a_i)
            c = (x - r_x) * chi_pow(-a_i)      # in O_pi
            if c.is_zero():
                continue
            for k_ in range(n):
                work.m[k_][j] = work.m[k_][j] - c * work.m[k_][i]
    return work


# ----------------------------------------------------------------------------
# Tree vertices
# ----------------------------------------------------------------------------

class Vertex:
    """A vertex of the tree B.

    kind = 'P' (pure): basis g of a self-dual lattice, v_pi(det)=0.
    kind = 'A' (alternating): basis g of the BIG lattice L (L# < L, index 9),
           v_pi(det) = -1.
    """
    __slots__ = ("kind", "g", "_key")

    def __init__(self, kind, g):
        self.kind = kind
        # normalize the pi-class so v_pi(det) is 0 (pure) or -1 (alternating)
        vd = g.det().v_pi()
        target = 0 if kind == "P" else -1
        assert (vd - target) % 3 == 0, (kind, vd)
        shift = (target - vd) // 3
        if shift:
            g = g * chi_pow(shift)
        self.g = g
        self._key = (kind, hnf_local(g).key())

    def key(self):
        return self._key

    def __eq__(self, o):
        return isinstance(o, Vertex) and self._key == o._key

    def __hash__(self):
        return hash(self._key)


ORIGIN = None  # set below after neighbor machinery


# -- F_3 linear algebra -------------------------------------------------------

def f3_lines():
    """The 13 one-dimensional subspaces of F_3^3, as normalized generators."""
    seen, out = set(), []
    for v in itertools.product(range(3), repeat=3):
        if v == (0, 0, 0) or v in seen:
            continue
        seen.add(v)
        seen.add(tuple((2 * x) % 3 for x in v))
        out.append(v)
    return out


F3_LINES = f3_lines()


def residue_matrix(g: Mat):
    """Entry-wise residue mod chi of an integral matrix, as 3x3 ints."""
    return [[g.m[i][j].res3() for j in range(g.n)] for i in range(g.n)]


def _f3_bilinear(M, x, y):
    return sum(x[i] * M[i][j] * y[j] for i in range(3) for j in range(3)) % 3


def isotropic_lines(M):
    """Lines V with V subset V-perp for the form <x,y> = x^T M y."""
    return [v for v in F3_LINES if _f3_bilinear(M, v, v) == 0]


def f3_planes():
    """All 13 two-dimensional subspaces, as pairs of spanning vectors."""
    planes, seen = [], set()
    for v1 in F3_LINES:
        for v2 in F3_LINES:
            span = frozenset(
                tuple((a * v1[i] + b * v2[i]) % 3 for i in range(3))
                for a in range(3) for b in range(3))
            if len(span) == 9 and span not in seen:
                seen.add(span)
                planes.append((v1, v2, span))
    return planes


F3_PLANES = f3_planes()


def selfdual_planes(M):
    """2-dim subspaces V with V == V-perp for the form x^T M y (M antisym)."""
    out = []
    for v1, v2, span in F3_PLANES:
        if (_f3_bilinear(M, v1, v2) == 0 and _f3_bilinear(M, v2, v1) == 0
                and _f3_bilinear(M, v1, v1) == 0
                and _f3_bilinear(M, v2, v2) == 0):
            # V subset V-perp; dims forces V-perp == V when rank(M) == 2
            out.append((v1, v2, span))
    return out


def _complete_to_unimodular(*vecs):
    """Complete the given F_3 vectors (lifted to {0,1,2} ints) to a basis
    of O^3: returns a 3x3 integer Mat with the given vectors as first
    columns and unit determinant mod 3."""
    cols = [list(v) for v in vecs]
    for e in ([1, 0, 0], [0, 1, 0], [0, 0, 1]):
        if len(cols) == 3:
            break
        trial = cols + [e]
        d = Mat([[Fw(trial[j][i]) for j in range(len(trial))] +
                 [ZERO] * (3 - len(trial)) for i in range(3)])
        # check rank over F_3 by determinant of the filled-out square later;
        # simpler: test linear independence mod 3 via brute force
        if _f3_rank([c[:] for c in trial]) == len(trial):
            cols.append(e)
    assert len(cols) == 3
    m = Mat([[Fw(cols[j][i]) for j in range(3)] for i in range(3)])
    assert m.det().v_pi() == 0, "completion not unimodular"
    return m


def _f3_rank(cols):
    """Rank over F_3 of a list of length-3 integer column vectors."""
    rows = [[c[i] % 3 for c in cols] for i in range(3)]
    rank, pr = 0, 0
    for pc in range(len(cols)):
        piv = next((r for r in range(pr, 3) if rows[r][pc] % 3), None)
        if piv is None:
            continue
        rows[pr], rows[piv] = rows[piv], rows[pr]
        inv = pow(rows[pr][pc], -1, 3)
        rows[pr] = [(x * inv) % 3 for x in rows[pr]]
        for r in range(3):
            if r != pr and rows[r][pc] % 3:
                f = rows[r][pc]
                rows[r] = [(rows[r][k] - f * rows[pr][k]) % 3
                           for k in range(len(cols))]
        pr += 1
        rank += 1
    return rank


# -- neighbors ----------------------------------------------------------------

def neighbors_of_pure(v: Vertex):
    """The 4 alternating neighbors of a pure vertex.

    Residue Gram form on L/pi L is symmetric non-degenerate; each isotropic
    line V = span(u) lifts to the big lattice  L_1 = g*(O^3 + chi^-1 u O).
    """
    assert v.kind == "P"
    M = residue_matrix(gram(v.g))
    lines = isotropic_lines(M)
    out = []
    for u in lines:
        U = _complete_to_unimodular(u)
        g1 = v.g * U * Mat.diag(CHI_INV, 1, 1)
        out.append(Vertex("A", g1))
    return out


def neighbors_of_alternating(v: Vertex):
    """The 4 pure neighbors of an alternating vertex.

    Residue form chi*<,> mod pi on L/pi L is antisymmetric of rank 2; each
    self-dual plane V = span(w1,w2) lifts to the self-dual lattice
    L_1 = g*(w1 O + w2 O + chi w3 O).
    """
    assert v.kind == "A"
    M = residue_matrix(gram(v.g) * CHI)
    planes = selfdual_planes(M)
    out = []
    for w1, w2, _span in planes:
        U = _complete_to_unimodular(w1, w2)
        g1 = v.g * U * Mat.diag(1, 1, CHI)
        out.append(Vertex("P", g1))
    return out


def neighbors(v: Vertex):
    return neighbors_of_pure(v) if v.kind == "P" else neighbors_of_alternating(v)


ORIGIN = Vertex("P", Mat.identity())


def bfs_tree(radius, start=None):
    """BFS from the origin (or `start`).  Returns (levels, parent, cycle_edges)
    where levels[r] is the list of vertices at distance r, parent maps vertex
    key -> parent vertex key, and cycle_edges collects any non-tree edge
    (must be empty iff B is a tree, locally)."""
    if start is None:
        start = ORIGIN
    seen = {start.key(): 0}
    parent = {start.key(): None}
    levels = [[start]]
    cycle_edges = []
    frontier = [start]
    for r in range(1, radius + 1):
        nxt = []
        for v in frontier:
            for w in neighbors(v):
                wk = w.key()
                if wk == parent[v.key()]:
                    continue
                if wk in seen:
                    cycle_edges.append((v.key(), wk))
                    continue
                seen[wk] = r
                parent[wk] = v.key()
                nxt.append(w)
        levels.append(nxt)
        frontier = nxt
    return levels, parent, cycle_edges


# ----------------------------------------------------------------------------
# Exact synthesis  (walk to the origin; proof of Theorem th:ring_equality)
# ----------------------------------------------------------------------------

def s0_coset_reps():
    """12 monomial matrices m such that {m * H_GATE * O^3} = S_0, the set of
    pure vertices at distance 2 from the origin."""
    reps, seen = [], set()
    e1 = H_GATE
    for m in monomial_matrices():
        g = m * e1
        v = Vertex("P", g)
        if v.key() not in seen:
            seen.add(v.key())
            reps.append(m)
    assert len(reps) == 12, len(reps)
    return reps


_S0_REPS = None


def synthesis_steps(U: Mat):
    """Decompose U in Gamma as  (m_1 H)(m_2 H)...(m_n H) * M  with m_i
    monomial and M monomial.  Returns (list_of_m, M).  Exact."""
    global _S0_REPS
    if _S0_REPS is None:
        _S0_REPS = s0_coset_reps()
    assert in_Gamma(U), "input not in U_3(Z[omega,1/3])"
    word = []
    cur = U
    while ell(cur) > 0:
        target = ell(cur) - 2
        for m in _S0_REPS:
            cand = (m * H_GATE).inv() * cur
            if ell(cand) == target:
                word.append(m)
                cur = cand
                break
        else:
            raise RuntimeError("no descending step found -- theory violated?")
    return word, cur


def reconstruct(word, M):
    out = Mat.identity()
    for m in word:
        out = out * (m * H_GATE)
    return out * M
