#!/usr/bin/env python
"""
v05_cartan.py -- Verification of Cartan decomposition and the l-function
lemmas of "Buildings for Synthesis with Clifford+R".

Claims verified (paper lines ~324-334, 929-1005, appendix ~1364-1507):

  [0] Sanity of the underlying primitives (btlib is a tool, not ground truth):
      v_pi cross-checked against an independent norm-based valuation
      (v_pi(x) = v_3(N(num)) - 2 v_3(den), valid since chi is the unique,
      ramified, residue-degree-1 prime over 3), matrix algebra cross-checked
      against numpy complex floats, Gram/self-duality reimplemented from
      scratch.
  [1] Theorem th:cartan -- Smith normal form over the DVR O_pi, implemented
      from scratch: g = k a k', k,k' in GL_3(O_pi),
      a = diag(chi^{l1},chi^{l2},chi^{l3}), l1>=l2>=l3, with uniqueness
      cross-checked via minors (l3 = min v(entries), l3+l2 = min v(2x2
      minors), sum = v_pi(det)) and via invariance under GL_3(O) translation.
  [2] Lemma le:diagonals -- g in A = {g : g* g in GL_3(O_pi)} has Cartan
      exponents (n, 0, -n), n >= 0.
  [3] Lemma le:ispositive (all 5 parts), randomized.  Also: A is NOT closed
      under inverses (explicit example), which the lemma carefully avoids
      claiming.
  [4] Lemma le:hecke_nghbs both directions; the (4,0)-impossibility step and
      the congruence-subgroup argument of its proof.
  [5] Proposition pr:chain (Appendix E): Cartan basis Gram conditions and the
      interpolating self-dual lattices L_i = chi^i v1 O + v2 O + chi^-i v3 O.

Deterministic (seeded).  Exits nonzero on any FAIL.
"""

import sys
import itertools
import random
from math import inf

sys.path.insert(0, "/Users/markdeaconu/projects/qutrits_v2")

import numpy as np

from btlib import (Zw, Fw, Mat, ZERO, ONE, OMEGA, CHI, CHI_INV, chi_pow,
                   H_GATE, S_GATE, R_GATE, Vertex, neighbors,
                   gram as bt_gram, is_self_dual as bt_is_self_dual,
                   ell as bt_ell, d_tilde as bt_d_tilde)

RNG = random.Random(20260609)
np.random.seed(20260609)

FAILURES = []


def check(name, ok, detail=""):
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f"  -- {detail}" if detail else ""))
    if not ok:
        FAILURES.append(name)


# =============================================================================
# Independent primitives (NOT relying on btlib's v_chi / gram / ell)
# =============================================================================

def v3int(n):
    """3-adic valuation of a nonzero integer."""
    n = abs(int(n))
    assert n != 0
    v = 0
    while n % 3 == 0:
        n //= 3
        v += 1
    return v


def vpi(x: Fw):
    """Independent v_pi: chi is the unique prime of Z[omega] over 3, ramified
    with residue degree 1, so v_chi(z) = v_3(N(z)) for z in Z[omega] where
    N(a+bw) = a^2 - ab + b^2; and v_pi(3) = 2 handles the denominator."""
    if x.num.is_zero():
        return inf
    return v3int(x.num.norm()) - 2 * v3int(x.den)


def min_v(M: Mat):
    return min(vpi(M.m[i][j]) for i in range(M.n) for j in range(M.n))


def ell(M: Mat):
    """l(g) = -2 min_ij v_pi(g_ij)   [eq:defi_of_l], independent v_pi."""
    return -2 * min_v(M)


def is_O(M: Mat):
    return min_v(M) >= 0


def in_GLO(M: Mat):
    return is_O(M) and vpi(M.det()) == 0


def herm(x, y):
    """<x,y> = sum_i x_i conj(y_i)  (paper convention)."""
    return sum((xi * yi.conj() for xi, yi in zip(x, y)), ZERO)


def my_gram(M: Mat):
    """Gram matrix [<v_i, v_j>] of the COLUMNS v_j of M."""
    cols = [M.col(j) for j in range(M.n)]
    return Mat([[herm(cols[i], cols[j]) for j in range(M.n)]
                for i in range(M.n)])


def self_dual(M: Mat):
    """Lemma le:gramm_matrix, reimplemented."""
    return in_GLO(my_gram(M))


def d_til(g: Mat, h: Mat):
    m = g.inv() * h
    s = ell(m) + ell(m.inv())
    assert s % 2 == 0
    return s // 2


def lat_eq(g, h):
    return in_GLO(g.inv() * h)


def lat_contains(g, h):
    """g O^3 contains h O^3 ?"""
    return is_O(g.inv() * h)


def in_A(g):
    """g in A = {g : g* g in GL_3(O_pi)}  (Lemma le:Aset)."""
    return in_GLO(g.conjT() * g)


# complex embedding for float cross-checks
OMC = complex(-0.5, 3 ** 0.5 / 2)


def to_c(x: Fw):
    return (x.num.a + x.num.b * OMC) / x.den


def mat_c(M: Mat):
    return np.array([[to_c(M.m[i][j]) for j in range((M.n))] for i in range(M.n)])


# =============================================================================
# Random generators (seeded RNG passed explicitly)
# =============================================================================

def rand_zw(rng, B=6):
    return Zw(rng.randint(-B, B), rng.randint(-B, B))


def rand_O(rng):
    """Random element of O_pi (denominator prime to 3 allowed)."""
    return Fw(rand_zw(rng), rng.choice([1, 1, 1, 1, 2, 5]))


def rand_unitO(rng):
    """Random unit of O_pi (v_pi = 0)."""
    while True:
        z = rand_zw(rng, 4)
        if not z.is_zero() and z.res3() != 0:
            return Fw(z, rng.choice([1, 1, 1, 2]))


def rand_perm_mat(rng):
    p = list(range(3))
    rng.shuffle(p)
    return Mat([[ONE if p[i] == j else ZERO for j in range(3)] for i in range(3)])


def rand_GLO(rng):
    """Random element of GL_3(O_pi): lower-tri x upper-tri (unit diagonals,
    integral entries) x permutation."""
    L = Mat([[rand_unitO(rng) if i == j else (rand_O(rng) if i > j else ZERO)
              for j in range(3)] for i in range(3)])
    U = Mat([[rand_unitO(rng) if i == j else (rand_O(rng) if i < j else ZERO)
              for j in range(3)] for i in range(3)])
    M = L * U * rand_perm_mat(rng)
    assert in_GLO(M)
    return M


def rand_GLF(rng):
    """Random g in GL_3(Q(omega)), entries with denominators incl. powers of 3,
    overall chi-power twist."""
    while True:
        M = Mat([[Fw(rand_zw(rng, 9), rng.choice([1, 1, 3, 9, 2]))
                  for _ in range(3)] for _ in range(3)])
        if not M.det().is_zero():
            return M * chi_pow(rng.randint(-2, 2))


GATES = [H_GATE, S_GATE, R_GATE]


def rand_gamma(rng, lo=2, hi=10):
    """Random word in the gates H, S, R (an element of U_3(Z[1/3,omega]))."""
    g = Mat.identity()
    for _ in range(rng.randint(lo, hi)):
        g = g * rng.choice(GATES)
    return g


PURE_POOL = []   # bases of pure tree vertices, filled by the walks


def rand_A(rng):
    """Random element of A: (basis of a self-dual lattice) * GL_3(O)."""
    if PURE_POOL and rng.random() < 0.5:
        base = rng.choice(PURE_POOL)
    else:
        base = rand_gamma(rng)
    g = base * rand_GLO(rng)
    assert in_A(g)
    return g


def random_walk(rng, steps, start=None):
    """Non-backtracking walk on the tree, starting at a pure vertex.
    Checks 4-regularity and adjacency symmetry along the way."""
    v = start if start is not None else Vertex("P", rand_gamma(rng))
    path = [v]
    for _ in range(steps):
        nbrs = neighbors(path[-1])
        assert len(nbrs) == 4, f"vertex not 4-regular: {len(nbrs)}"
        if len(path) > 1:
            prev = path[-2].key()
            keys = [w.key() for w in nbrs]
            assert prev in keys, "adjacency not symmetric"
            cand = [w for w in nbrs if w.key() != prev]
            assert len(cand) == 3
        else:
            cand = nbrs
        path.append(rng.choice(cand))
    return path


# =============================================================================
# Smith normal form over the DVR O_pi (from scratch)
# =============================================================================

def snf_dvr(g: Mat):
    """g = K * diag(chi^e1, chi^e2, chi^e3) * Kp with K, Kp in GL_3(O_pi)
    and e1 <= e2 <= e3 (ascending).  Pivot = entry of minimal v_pi."""
    n = g.n
    A = Mat([row[:] for row in g.m])
    K = Mat.identity(n)
    Kp = Mat.identity(n)
    # invariant: g == K * A * Kp at all times
    for r in range(n):
        best = None
        for i in range(r, n):
            for j in range(r, n):
                v = vpi(A.m[i][j])
                if v != inf and (best is None or v < best[0]):
                    best = (v, i, j)
        assert best is not None, "singular matrix"
        v, bi, bj = best
        if bi != r:                      # row swap; K <- K * P
            A.m[r], A.m[bi] = A.m[bi], A.m[r]
            for row in K.m:
                row[r], row[bi] = row[bi], row[r]
        if bj != r:                      # col swap; Kp <- P * Kp
            for row in A.m:
                row[r], row[bj] = row[bj], row[r]
            Kp.m[r], Kp.m[bj] = Kp.m[bj], Kp.m[r]
        u = A.m[r][r] * chi_pow(-v)      # unit of O_pi
        ui = u.inv()
        for j in range(n):               # scale row r by u^-1; K col r by u
            A.m[r][j] = A.m[r][j] * ui
        for i in range(n):
            K.m[i][r] = K.m[i][r] * u
        for i in range(n):               # clear column r
            if i == r or A.m[i][r].is_zero():
                continue
            c = A.m[i][r] * chi_pow(-v)  # in O_pi by pivot minimality
            assert vpi(c) >= 0
            for j in range(n):
                A.m[i][j] = A.m[i][j] - c * A.m[r][j]
            for t in range(n):           # K <- K (I + c e_{ir})
                K.m[t][r] = K.m[t][r] + c * K.m[t][i]
        for j in range(n):               # clear row r
            if j == r or A.m[r][j].is_zero():
                continue
            c = A.m[r][j] * chi_pow(-v)
            assert vpi(c) >= 0
            for i in range(n):
                A.m[i][j] = A.m[i][j] - c * A.m[i][r]
            for t in range(n):           # Kp <- (I + c e_{rj}) Kp
                Kp.m[r][t] = Kp.m[r][t] + c * Kp.m[j][t]
    exps = [vpi(A.m[i][i]) for i in range(n)]
    return K, A, Kp, exps


P_REV = Mat([[0, 0, 1], [0, 1, 0], [1, 0, 0]])


def cartan(g: Mat):
    """Paper-order Cartan decomposition: g = k1 a k2 with
    a = diag(chi^{l1}, chi^{l2}, chi^{l3}), l1 >= l2 >= l3."""
    K, A, Kp, exps = snf_dvr(g)
    k1 = K * P_REV
    a = P_REV * A * P_REV
    k2 = P_REV * Kp
    return k1, a, k2, list(reversed(exps))


def min_minor2_v(g: Mat):
    best = inf
    for rows in itertools.combinations(range(3), 2):
        for cols in itertools.combinations(range(3), 2):
            d = (g.m[rows[0]][cols[0]] * g.m[rows[1]][cols[1]]
                 - g.m[rows[0]][cols[1]] * g.m[rows[1]][cols[0]])
            v = vpi(d)
            if v < best:
                best = v
    return best


# =============================================================================
# Section 0: sanity of primitives
# =============================================================================

def section0():
    print("\n--- Section 0: sanity of primitives (btlib vs independent code) ---")

    # specific values
    ok = (vpi(Fw(3)) == 2 and vpi(CHI) == 1 and vpi(CHI_INV) == -1
          and vpi(OMEGA) == 0 and vpi(Fw(Zw(1, 2), 3)) == -1)
    check("0.1 v_pi specific values (v(3)=2, v(chi)=1, v(1/chi)=-1, v(w)=0, v(i/sqrt3)=-1)", ok)

    # v_pi vs btlib on random elements
    n_ok = 0
    N = 2000
    for _ in range(N):
        x = Fw(rand_zw(RNG, 30), RNG.choice([1, 2, 3, 9, 27, 5, 6, 15]))
        if x.num.is_zero():
            x = x + ONE
        if vpi(x) == x.v_pi():
            n_ok += 1
        # float consistency of the norm
        cx = to_c(x)
        assert abs(abs(cx) ** 2 - x.num.norm() / x.den ** 2) < 1e-6 * max(1.0, abs(cx) ** 2)
    check("0.2 independent v_pi == btlib v_pi", n_ok == N, f"{n_ok}/{N} random elements")

    # matrix algebra vs numpy
    n_ok = 0
    N = 20
    for _ in range(N):
        g = rand_GLF(RNG)
        h = rand_GLF(RNG)
        good = (g * g.inv() == Mat.identity())
        good &= np.allclose(mat_c(g.inv()), np.linalg.inv(mat_c(g)), atol=1e-8)
        good &= np.allclose(mat_c(g * h), mat_c(g) @ mat_c(h), atol=1e-8)
        good &= np.allclose(mat_c(g.conjT()), mat_c(g).conj().T, atol=1e-12)
        good &= abs(np.linalg.det(mat_c(g)) - to_c(g.det())) < 1e-6 * max(1.0, abs(to_c(g.det())))
        if good:
            n_ok += 1
    check("0.3 Mat inverse/product/conjT/det vs numpy floats", n_ok == N, f"{n_ok}/{N}")

    # gram matrix: my_gram vs float and vs btlib's A* A
    n_ok = 0
    N = 50
    for _ in range(N):
        g = rand_GLF(RNG)
        G = my_gram(g)
        gc = mat_c(g)
        good = np.allclose(mat_c(G), gc.T @ gc.conj(), atol=1e-8)
        # btlib gram = A* A; paper Gram [<v_i,v_j>] is its transpose
        good &= (bt_gram(g).T() == G)
        good &= (self_dual(g) == bt_is_self_dual(g))
        good &= (ell(g) == bt_ell(g))
        if good:
            n_ok += 1
    check("0.4 Gram (= (A*A)^T), self-duality, ell vs btlib + numpy", n_ok == N, f"{n_ok}/{N}")

    # gates
    ok = (H_GATE.is_unitary() and S_GATE.is_unitary() and R_GATE.is_unitary()
          and np.allclose(mat_c(H_GATE) @ mat_c(H_GATE).conj().T, np.eye(3), atol=1e-12)
          and ell(H_GATE) == 2)
    check("0.5 H,S,R unitary (exact + float), l(H)=2", ok)

    # self-duality basics
    ok = (self_dual(Mat.identity()) and self_dual(H_GATE)
          and not self_dual(Mat.identity() * CHI))
    check("0.6 self-duality basics: O^3, H O^3 self-dual; chi O^3 not", ok)


# =============================================================================
# Section 1: Smith normal form / Cartan decomposition (Theorem th:cartan)
# =============================================================================

def section1():
    print("\n--- Section 1: Cartan decomposition (SNF over O_pi, from scratch) ---")
    N = 100
    n_recon = n_glo = n_diag = n_det = n_min = n_minor = n_uni = 0
    for t in range(N):
        g = rand_GLF(RNG)
        K, A, Kp, exps = snf_dvr(g)
        if K * A * Kp == g:
            n_recon += 1
        if in_GLO(K) and in_GLO(Kp):
            n_glo += 1
        diag_ok = all(A.m[i][j].is_zero() for i in range(3) for j in range(3) if i != j)
        diag_ok &= all(A.m[i][i] == chi_pow(exps[i]) for i in range(3))
        diag_ok &= exps[0] <= exps[1] <= exps[2]
        if diag_ok:
            n_diag += 1
        if sum(exps) == vpi(g.det()):
            n_det += 1
        if exps[0] == min_v(g):           # smallest exponent = min valuation
            n_min += 1
        if exps[0] + exps[1] == min_minor2_v(g):
            n_minor += 1
        # uniqueness: GL_3(O)-translates have the same exponents
        k, kp = rand_GLO(RNG), rand_GLO(RNG)
        _, _, _, exps2 = snf_dvr(k * g * kp)
        if exps2 == exps:
            n_uni += 1
    check("1.1 reconstruction g == K a K'", n_recon == N, f"{n_recon}/{N} random g in GL_3(Q(w))")
    check("1.2 K, K' in GL_3(O_pi)", n_glo == N, f"{n_glo}/{N}")
    check("1.3 a = diag(chi^e), exponents sorted", n_diag == N, f"{n_diag}/{N}")
    check("1.4 sum(exps) == v_pi(det g)", n_det == N, f"{n_det}/{N}")
    check("1.5 smallest exponent (paper l3) == min_ij v_pi(g_ij)", n_min == N, f"{n_min}/{N}")
    check("1.6 e1+e2 == min v_pi over 2x2 minors (elementary divisors)", n_minor == N, f"{n_minor}/{N}")
    check("1.7 uniqueness: exponents invariant under k g k', k,k' in GL_3(O)", n_uni == N, f"{n_uni}/{N}")


# =============================================================================
# Section 2: Lemma le:diagonals
# =============================================================================

def section2():
    print("\n--- Section 2: Lemma le:diagonals (Cartan exponents of A are (n,0,-n)) ---")
    N = 50
    n_ok = 0
    n_seen = {}
    samples = []
    for _ in range(N - 10):
        samples.append(rand_A(RNG))
    # also bases of random tree pure vertices
    for _ in range(10):
        path = random_walk(RNG, 2 * RNG.randint(1, 2))
        v = path[-1]
        assert v.kind == "P"
        samples.append(v.g * rand_GLO(RNG))
    for g in samples:
        assert in_A(g), "generator broken: not in A"
        k1, a, k2, exps = cartan(g)
        ok = (k1 * a * k2 == g) and (exps[1] == 0) and (exps[0] == -exps[2]) and exps[0] >= 0
        if ok:
            n_ok += 1
            n_seen[exps[0]] = n_seen.get(exps[0], 0) + 1
    check("2.1 g in A => Cartan a = diag(chi^n, 1, chi^-n), n>=0",
          n_ok == N, f"{n_ok}/{N}; n-histogram {dict(sorted(n_seen.items()))}")


# =============================================================================
# Section 3: Lemma le:ispositive (5 parts)
# =============================================================================

def section3():
    print("\n--- Section 3: Lemma le:ispositive ---")

    # (1) l(k g k') = l(g)
    N = 200
    n_ok = sum(1 for _ in range(N)
               if (lambda g, k, kp: ell(k * g * kp) == ell(g))
               (rand_GLF(RNG), rand_GLO(RNG), rand_GLO(RNG)))
    check("3.1 l(k g k') == l(g) for k,k' in GL_3(O)", n_ok == N, f"{n_ok}/{N}")

    # (2) subadditivity
    N = 200
    n_ok = 0
    for _ in range(N):
        g1, g2 = rand_GLF(RNG), rand_GLF(RNG)
        if ell(g1 * g2) <= ell(g1) + ell(g2):
            n_ok += 1
    check("3.2 l(g1 g2) <= l(g1) + l(g2)", n_ok == N, f"{n_ok}/{N}")

    # (3) g in A: l(g) >= 0 and l(g^-1) = l(g)
    N = 200
    n_ok = 0
    n_not_invclosed = 0
    for _ in range(N):
        g = rand_A(RNG)
        if ell(g) >= 0 and ell(g.inv()) == ell(g):
            n_ok += 1
        if not in_GLO(g * g.conjT()):
            n_not_invclosed += 1     # i.e. g^-1 not in A
    check("3.3 g in A: l(g) >= 0 and l(g^-1) == l(g)", n_ok == N, f"{n_ok}/{N}")

    # A is NOT inverse-closed: explicit example
    # g = basis of self-dual lattice  chi(1,0,0) O + (1,-1,0) O + chi^-1 (1,1,1) O
    g0 = Mat([[CHI, 1, CHI_INV],
              [ZERO, -ONE, CHI_INV],
              [ZERO, ZERO, CHI_INV]])
    in_a = in_A(g0) and self_dual(g0)
    gg_star = g0 * g0.conjT()
    not_glo = not in_GLO(gg_star)        # so (g0^-1)* g0^-1 = (g0 g0*)^-1 not in GL_3(O)
    inv_not_in_A = not in_A(g0.inv())
    l_match = (ell(g0.inv()) == ell(g0))
    check("3.3b explicit g0 in A with g0 g0* NOT in GL_3(O)  (A not inverse-closed)",
          in_a and not_glo and inv_not_in_A and l_match,
          f"g0 in A: {in_a}; min v_pi(g0 g0*) = {min_v(gg_star)}; "
          f"g0^-1 in A: {not inv_not_in_A}; still l(g0^-1)==l(g0)=={ell(g0)}; "
          f"random samples with g g* not in GL_3(O): {n_not_invclosed}/{N}")

    # (4) l(g^-1 h) >= 0 for g,h in A
    N = 200
    n_ok = 0
    for _ in range(N):
        g, h = rand_A(RNG), rand_A(RNG)
        if ell(g.inv() * h) >= 0:
            n_ok += 1
    check("3.4 l(g^-1 h) >= 0 for g,h in A", n_ok == N, f"{n_ok}/{N}")

    # (5) metric axioms on A / GL_3(O)
    N = 100
    n_sym = n_zero_self = n_iff = 0
    for _ in range(N):
        g, h = rand_A(RNG), rand_A(RNG)
        d1, d2 = d_til(g, h), d_til(h, g)
        if d1 == d2 and d1 == bt_d_tilde(g, h):
            n_sym += 1
        if d_til(g, g * rand_GLO(RNG)) == 0:
            n_zero_self += 1
        if (d_til(g, h) == 0) == lat_eq(g, h):
            n_iff += 1
    check("3.5a d~ symmetric (and == btlib d_tilde)", n_sym == N, f"{n_sym}/{N}")
    check("3.5b d~(g, g k) == 0 for k in GL_3(O)", n_zero_self == N, f"{n_zero_self}/{N}")
    check("3.5c d~(g,h) == 0  iff  same lattice", n_iff == N, f"{n_iff}/{N}")

    N = 200
    n_tri = 0
    n_nonneg = 0
    for _ in range(N):
        g, h, f = rand_A(RNG), rand_A(RNG), rand_A(RNG)
        dgh, dgf, dfh = d_til(g, h), d_til(g, f), d_til(f, h)
        if dgh <= dgf + dfh:
            n_tri += 1
        if dgh >= 0 and dgf >= 0 and dfh >= 0:
            n_nonneg += 1
    check("3.5d triangle inequality d~(g,h) <= d~(g,f) + d~(f,h)", n_tri == N, f"{n_tri}/{N} triples")
    check("3.5e d~ >= 0", n_nonneg == N, f"{n_nonneg}/{N}")


# =============================================================================
# Section 4: Lemma le:hecke_nghbs
# =============================================================================

def incl_in_pinv(g, h):
    """Lambda_g subseteq pi^-1 Lambda_h  iff  chi * h^-1 g integral."""
    return is_O((h.inv() * g) * CHI)


def dist2_pair(rng):
    """A pair (g, h, alt) of distinct self-dual lattice bases at tree
    distance 2 (pure -> alternating -> pure), plus the middle vertex."""
    while True:
        path = random_walk(rng, 2)
        p0, alt, p1 = path
        if p1.key() != p0.key():
            return p0.g, p1.g, alt


def section4():
    print("\n--- Section 4: Lemma le:hecke_nghbs ---")

    # forward: distance-2 pairs
    N = 50
    n_struct = n_d2 = n_incl = 0
    for _ in range(N):
        g, h, alt = dist2_pair(RNG)
        PURE_POOL.append(g)
        # independent structural certificate that the tree distance is 2:
        # both L_g and L_h are self-dual, distinct, and sandwich the SAME
        # alternating pair  Lam# subset L subset Lam  (each a 1-simplex)
        lam = alt.g                       # big lattice of the alternating vertex
        lam_dual = lam.conjT().inv()      # its dual (the small one)
        ok = self_dual(g) and self_dual(h) and not lat_eq(g, h)
        ok &= alt.kind == "A" and vpi(lam.det()) == -1
        ok &= lat_contains(lam, lam_dual) and vpi((lam.inv() * lam_dual).det()) == 2
        for x in (g, h):
            ok &= lat_contains(lam, x) and lat_contains(x, lam_dual)
            ok &= vpi((lam.inv() * x).det()) == 1     # index 3 each step
        if ok:
            n_struct += 1
        if d_til(g, h) == 2:
            n_d2 += 1
        if incl_in_pinv(g, h) and incl_in_pinv(h, g):
            n_incl += 1
    check("4.1 structural certificate: pure-alt-pure with proper inclusions, L_g != L_h",
          n_struct == N, f"{n_struct}/{N} pairs")
    check("4.2 d~(g,h) == 2 at tree distance 2", n_d2 == N, f"{n_d2}/{N}")
    check("4.3 inclusions L_g <= pi^-1 L_h and L_h <= pi^-1 L_g at distance 2",
          n_incl == N, f"{n_incl}/{N}")

    # converse: distance 4 / 6 pairs
    for dist, npairs in ((4, 20), (6, 20)):
        n_d = n_noincl = 0
        for _ in range(npairs):
            path = random_walk(RNG, dist)
            g, h = path[0].g, path[-1].g
            if d_til(g, h) == dist:
                n_d += 1
            if not (incl_in_pinv(g, h) and incl_in_pinv(h, g)):
                n_noincl += 1
        check(f"4.4 distance-{dist} pairs: d~ == {dist} (so != 2)",
              n_d == npairs, f"{n_d}/{npairs}")
        check(f"4.5 distance-{dist} pairs: inclusions fail",
              n_noincl == npairs, f"{n_noincl}/{npairs}")

    # (4,0)-impossibility, random search over A x A
    N = 500
    n_40 = 0
    n_uneq = 0
    lvals = {}
    for _ in range(N):
        g, h = rand_A(RNG), rand_A(RNG)
        l1, l2 = ell(g.inv() * h), ell(h.inv() * g)
        lvals[(l1, l2)] = lvals.get((l1, l2), 0) + 1
        if {l1, l2} == {4, 0}:
            n_40 += 1
        if l1 != l2:
            n_uneq += 1
    check("4.6 no g,h in A with {l(g^-1 h), l(h^-1 g)} == {4, 0}",
          n_40 == 0, f"0 hits in {N} random pairs; "
          f"pairs with l(g^-1h) != l(h^-1g): {n_uneq} (always equal, as forced by self-duality)")

    # targeted search at the level of the proof: x in GL_3(O), d=diag(chi^n,1,chi^-n)
    N = 500
    n_40 = 0
    n_l4 = 0
    for _ in range(N):
        x = rand_GLO(RNG)
        for n in (1, 2):
            d = Mat.diag(chi_pow(n), 1, chi_pow(-n))
            m = d.inv() * x * d
            if ell(m) == 4:
                n_l4 += 1
                if ell(m.inv()) == 0:
                    n_40 += 1
    check("4.7 no x in GL_3(O) with l(d^-1 x d)=4 and l(d^-1 x^-1 d)=0",
          n_40 == 0, f"0 hits / {N} random x (cases with l(d^-1 x d)=4: {n_l4})")

    # congruence-subgroup argument:
    # w in GL_3(O) with d^-1 w d integral  =>  d^-1 w^-1 d integral
    # (with x = w^-1 this is exactly: d^-1 x^-1 d integral => d^-1 x d integral)
    n_ok = 0
    n_tot = 0
    for n in (1, 2):
        d = Mat.diag(chi_pow(n), 1, chi_pow(-n))
        for _ in range(100):
            # w = L * D * U with U's strict upper entries in chi^n O, chi^2n O, chi^n O
            L = Mat([[ONE if i == j else (rand_O(RNG) if i > j else ZERO)
                      for j in range(3)] for i in range(3)])
            D = Mat.diag(rand_unitO(RNG), rand_unitO(RNG), rand_unitO(RNG))
            U = Mat([[ONE, rand_O(RNG) * chi_pow(n), rand_O(RNG) * chi_pow(2 * n)],
                     [ZERO, ONE, rand_O(RNG) * chi_pow(n)],
                     [ZERO, ZERO, ONE]])
            w = L * D * U
            assert in_GLO(w)
            assert is_O(d.inv() * w * d), "construction broken"
            n_tot += 1
            if is_O(d.inv() * w.inv() * d):
                n_ok += 1
    check("4.8 congruence subgroup closed under inverse: d^-1 x^-1 d in O => d^-1 x d in O",
          n_ok == n_tot, f"{n_ok}/{n_tot} (n=1 and n=2)")


# =============================================================================
# Section 5: Proposition pr:chain
# =============================================================================

def section5():
    print("\n--- Section 5: Proposition pr:chain (interpolating self-dual lattices) ---")
    cases = [(2, 12), (3, 10), (4, 8)]
    n_tot = 0
    n_cartan = n_gram = n_selfdual = n_consec = 0
    for n, npairs in cases:
        for _ in range(npairs):
            path = random_walk(RNG, 2 * n)
            g, h = path[0].g, path[-1].g
            assert self_dual(g) and self_dual(h)
            n_tot += 1
            # Cartan of g^-1 h must be diag(chi^n, 1, chi^-n)
            k1, a, k2, exps = cartan(g.inv() * h)
            ok_cart = (exps == [n, 0, -n]) and (k1 * a * k2 == g.inv() * h)
            ok_cart &= (d_til(g, h) == 2 * n)
            g1 = g * k1
            # Lambda_g = g1 O^3 and Lambda_h = g1 a O^3
            ok_cart &= lat_eq(g1, g) and lat_eq(g1 * a, h)
            if ok_cart:
                n_cartan += 1
            # Gram conditions from the proof
            v1, v2, v3 = (g1.col(j) for j in range(3))
            ok_gram = (vpi(herm(v2, v3)) >= n
                       and vpi(herm(v3, v3)) >= 2 * n
                       and vpi(herm(v1, v3)) == 0
                       and vpi(herm(v2, v2)) == 0)
            if ok_gram:
                n_gram += 1
            # interpolating lattices L_i, i = 0..n  (L_0 = Lam_g, L_n = Lam_h)
            Ls = [g1 * Mat.diag(chi_pow(i), 1, chi_pow(-i)) for i in range(n + 1)]
            if all(self_dual(L) for L in Ls):
                n_selfdual += 1
            ok_consec = True
            for i in range(n):
                ok_consec &= (not lat_eq(Ls[i], Ls[i + 1]))
                ok_consec &= (d_til(Ls[i], Ls[i + 1]) == 2)
                ok_consec &= incl_in_pinv(Ls[i], Ls[i + 1]) and incl_in_pinv(Ls[i + 1], Ls[i])
            if ok_consec:
                n_consec += 1
    check("5.1 Cartan of g^-1 h is (n,0,-n) with 2n = tree distance; Lam_h = (g k1) a O^3",
          n_cartan == n_tot, f"{n_cartan}/{n_tot} pairs (n=2,3,4)")
    check("5.2 Gram conditions: <v2,v3> in pi^n O, <v3,v3> in pi^2n O, <v1,v3>,<v2,v2> units",
          n_gram == n_tot, f"{n_gram}/{n_tot}")
    check("5.3 all interpolating L_i = chi^i v1 O + v2 O + chi^-i v3 O self-dual (i=0..n)",
          n_selfdual == n_tot, f"{n_selfdual}/{n_tot}")
    check("5.4 consecutive L_i, L_i+1 distinct, d~ == 2, mutual pi^-1-inclusions",
          n_consec == n_tot, f"{n_consec}/{n_tot}")


# =============================================================================

def main():
    section0()
    section1()
    section2()
    section3()
    section4()
    section5()
    print()
    if FAILURES:
        print(f"OVERALL: FAIL ({len(FAILURES)} failed checks)")
        for f in FAILURES:
            print("  failed:", f)
        sys.exit(1)
    print("OVERALL: PASS (all checks passed)")
    sys.exit(0)


if __name__ == "__main__":
    main()
