"""
v07_distance_audit.py -- Adversarial audit of the metric/distance results of
"Buildings for Synthesis with Clifford+R":

    Lemma  le:Aset          (paper.tex 929-940)
    Lemma  le:diagonals     (paper.tex 943-950)
    Lemma  le:ispositive    (paper.tex 960-975, Appendix E 1419-1476)
    Lemma  le:hecke_nghbs   (paper.tex 981-990, Appendix F 1477-1507)
    Prop   pr:chain         (paper.tex 729-748, Appendix D 1364-1418)
    Prop   prop:distance    (paper.tex 992-1003, Appendix G 1508-1547)

btlib.py is treated as a TOOL, not ground truth: every btlib primitive used
here is cross-checked against from-scratch reimplementations (norm-based
valuation, hand-rolled Smith normal form over the DVR, independent Gram /
self-duality / lattice-sum code) in Section A.

Deterministic (seeded).  Prints PASS/FAIL per claim; exits nonzero on FAIL.

Run:  /Users/markdeaconu/projects/qutrits_v2/.venv/bin/python v07_distance_audit.py
"""

import sys, os, random, itertools
from math import inf

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
from btlib import (Zw, Fw, Mat, ZERO, ONE, OMEGA, CHI, CHI_INV, chi_pow,
                   H_GATE, S_GATE, R_GATE, Vertex, neighbors, hnf_local,
                   ORIGIN)

RNG = random.Random(20260610)
np.random.seed(20260610)

FAILURES = []


def check(name, ok, detail=""):
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f"  -- {detail}" if detail else ""))
    if not ok:
        FAILURES.append(name)


# =============================================================================
# Independent primitives (do NOT trust btlib's v_chi / ell / gram / hnf)
# =============================================================================

def v3int(n):
    n = abs(int(n)); assert n
    v = 0
    while n % 3 == 0:
        n //= 3; v += 1
    return v


def vpi(x: Fw):
    """v_pi via the norm: chi is the unique (ramified) prime over 3 with
    residue degree 1, so v_chi(z) = v_3(N(z)) for z in Z[omega],
    N(a+bw)=a^2-ab+b^2, and v_pi(integer d) = 2 v_3(d)."""
    if x.num.is_zero():
        return inf
    return v3int(x.num.norm()) - 2 * v3int(x.den)


def min_v(M):  return min(vpi(M.m[i][j]) for i in range(M.n) for j in range(M.n))
def ell(M):    return -2 * min_v(M)                       # eq:defi_of_l
def is_O(M):   return min_v(M) >= 0
def in_GLO(M): return is_O(M) and vpi(M.det()) == 0


def herm(x, y):
    """<x,y> = sum x_i conj(y_i)  (paper convention)."""
    return sum((a * b.conj() for a, b in zip(x, y)), ZERO)


def my_gram(M):
    cols = [M.col(j) for j in range(3)]
    return Mat([[herm(cols[i], cols[j]) for j in range(3)] for i in range(3)])


def self_dual(M):   return in_GLO(my_gram(M))             # le:gramm_matrix
def dualb(M):       return M.conjT().inv()                # basis of dual lattice
def lat_contains(g, h):  return is_O(g.inv() * h)         # gO^3 >= hO^3
def lat_eq(g, h):        return in_GLO(g.inv() * h)
def in_A(g):             return in_GLO(g.conjT() * g)     # le:Aset definition


def d_til(g, h):
    m = g.inv() * h
    s = ell(m) + ell(m.inv())
    return s / 2


# ---- Smith normal form over the DVR O_pi, from scratch ----------------------

def snf_dvr(g):
    """g = K diag(chi^e1..) Kp, K,Kp in GL_3(O), e ascending."""
    n = g.n
    A = Mat([row[:] for row in g.m]); K = Mat.identity(n); Kp = Mat.identity(n)
    for r in range(n):
        best = None
        for i in range(r, n):
            for j in range(r, n):
                v = vpi(A.m[i][j])
                if v != inf and (best is None or v < best[0]):
                    best = (v, i, j)
        v, bi, bj = best
        if bi != r:
            A.m[r], A.m[bi] = A.m[bi], A.m[r]
            for row in K.m: row[r], row[bi] = row[bi], row[r]
        if bj != r:
            for row in A.m: row[r], row[bj] = row[bj], row[r]
            Kp.m[r], Kp.m[bj] = Kp.m[bj], Kp.m[r]
        u = A.m[r][r] * chi_pow(-v); ui = u.inv()
        for j in range(n): A.m[r][j] = A.m[r][j] * ui
        for i in range(n): K.m[i][r] = K.m[i][r] * u
        for i in range(n):
            if i == r or A.m[i][r].is_zero(): continue
            c = A.m[i][r] * chi_pow(-v)
            for j in range(n): A.m[i][j] = A.m[i][j] - c * A.m[r][j]
            for t in range(n): K.m[t][r] = K.m[t][r] + c * K.m[t][i]
        for j in range(n):
            if j == r or A.m[r][j].is_zero(): continue
            c = A.m[r][j] * chi_pow(-v)
            for i in range(n): A.m[i][j] = A.m[i][j] - c * A.m[i][r]
            for t in range(n): Kp.m[r][t] = Kp.m[r][t] + c * Kp.m[j][t]
    return K, A, Kp, [vpi(A.m[i][i]) for i in range(n)]


P_REV = Mat([[0, 0, 1], [0, 1, 0], [1, 0, 0]])


def cartan(g):
    """g = k1 a k2, a = diag(chi^{l1},chi^{l2},chi^{l3}), l1>=l2>=l3."""
    K, A, Kp, e = snf_dvr(g)
    return K * P_REV, P_REV * A * P_REV, P_REV * Kp, list(reversed(e))


def cartan_exps(g):
    return cartan(g)[3]


# ---- independent column-HNF / lattice sum -----------------------------------

def span_basis(cols):
    """Triangular O_pi-basis of the span of >=3 columns (from scratch)."""
    cols = [list(c) for c in cols]
    n = len(cols)
    for r in range(3):
        best, bv = None, None
        for j in range(r, n):
            v = vpi(cols[j][r])
            if v != inf and (bv is None or v < bv):
                best, bv = j, v
        cols[r], cols[best] = cols[best], cols[r]
        u = chi_pow(bv) / cols[r][r]
        cols[r] = [x * u for x in cols[r]]
        for j in range(r + 1, n):
            c = cols[j][r] * chi_pow(-bv)
            if not c.is_zero():
                cols[j] = [cols[j][i] - c * cols[r][i] for i in range(3)]
    return Mat([[cols[j][i] for j in range(3)] for i in range(3)])


def mat_cols(g):       return [g.col(j) for j in range(3)]
def lat_sum(g, h):     return span_basis(mat_cols(g) + mat_cols(h))
def lat_cap(g, h):     return dualb(lat_sum(dualb(g), dualb(h)))


# ---- random generators -------------------------------------------------------

def rand_zw(B=6):   return Zw(RNG.randint(-B, B), RNG.randint(-B, B))
def rand_O():       return Fw(rand_zw(), RNG.choice([1, 1, 1, 2, 5]))


def rand_unitO():
    while True:
        z = rand_zw(4)
        if not z.is_zero() and z.res3() != 0:
            return Fw(z, RNG.choice([1, 1, 2]))


def rand_perm():
    p = list(range(3)); RNG.shuffle(p)
    return Mat([[ONE if p[i] == j else ZERO for j in range(3)] for i in range(3)])


def rand_GLO():
    L = Mat([[rand_unitO() if i == j else (rand_O() if i > j else ZERO)
              for j in range(3)] for i in range(3)])
    U = Mat([[rand_unitO() if i == j else (rand_O() if i < j else ZERO)
              for j in range(3)] for i in range(3)])
    M = L * U * rand_perm()
    assert in_GLO(M)
    return M


def rand_GLF():
    while True:
        M = Mat([[Fw(rand_zw(9), RNG.choice([1, 1, 3, 9, 2]))
                  for _ in range(3)] for _ in range(3)])
        if not M.det().is_zero():
            return M * chi_pow(RNG.randint(-2, 2))


GATES = [H_GATE, S_GATE, R_GATE]


def rand_gamma(lo=2, hi=9):
    g = Mat.identity()
    for _ in range(RNG.randint(lo, hi)):
        g = g * RNG.choice(GATES)
    return g


PURE_POOL = []


def rand_A():
    base = RNG.choice(PURE_POOL) if (PURE_POOL and RNG.random() < .5) else rand_gamma()
    g = base * rand_GLO()
    assert in_A(g)
    return g


def canon(v):
    return Vertex(v.kind, hnf_local(v.g))


def random_walk(steps, start=None):
    """Non-backtracking walk on the tree from a pure vertex.  In a tree this
    ends at graph distance == steps (independently re-certified in sec G)."""
    v = start if start is not None else canon(Vertex("P", rand_gamma()))
    path = [v]
    for _ in range(steps):
        nbrs = [canon(w) for w in neighbors(path[-1])]
        assert len(nbrs) == 4 and len({w.key() for w in nbrs}) == 4
        if len(path) > 1:
            prev = path[-2].key()
            assert prev in [w.key() for w in nbrs], "adjacency not symmetric"
            nbrs = [w for w in nbrs if w.key() != prev]
            assert len(nbrs) == 3
        path.append(RNG.choice(nbrs))
    return path


def graph_dist(v, w, maxr=3):
    """Exact graph distance (<= 2*maxr) by bidirectional BFS using btlib
    neighbors (whose correctness against brute force is v06's job; local
    4-regularity and symmetry are re-asserted in random_walk)."""
    lv = [{v.key()}]; lw = [{w.key()}]
    pv = {v.key(): v}; pw = {w.key(): w}
    seenv = set(lv[0]); seenw = set(lw[0])
    fv, fw = [v], [w]
    for _ in range(maxr):
        nxt = []
        for u in fv:
            for nb in neighbors(u):
                nb = canon(nb)
                if nb.key() not in seenv:
                    seenv.add(nb.key()); nxt.append(nb)
        lv.append({u.key() for u in nxt}); fv = nxt
        nxt = []
        for u in fw:
            for nb in neighbors(u):
                nb = canon(nb)
                if nb.key() not in seenw:
                    seenw.add(nb.key()); nxt.append(nb)
        lw.append({u.key() for u in nxt}); fw = nxt
    best = inf
    for i, si in enumerate(lv):
        for j, sj in enumerate(lw):
            if si & sj:
                best = min(best, i + j)
    return best


# =============================================================================
print("=" * 78)
print("SECTION A: sanity of primitives (btlib vs from-scratch vs floats)")
print("=" * 78)

samples = [CHI, Fw(3), OMEGA, ONE, CHI_INV, CHI * CHI, Fw(Zw(1, 2), 3),
           Fw(Zw(-2, 5), 9)] + [Fw(rand_zw(9), RNG.choice([1, 2, 3, 9]))
                                for _ in range(400)]
samples = [x for x in samples if not x.is_zero()]
check("A1 v_pi (norm-based) == btlib v_pi on 400+ samples",
      all(vpi(x) == x.v_pi() for x in samples)
      and vpi(CHI) == 1 and vpi(Fw(3)) == 2 and vpi(CHI_INV) == -1)

cc = CHI.conj() + OMEGA * OMEGA * CHI           # conj(chi) - (-w^2 chi) == 0
check("A2 conj(chi) == -omega^2 * chi (unit times chi) and v(conj x)==v(x)",
      cc.is_zero() and all(vpi(x.conj()) == vpi(x) for x in samples))

ok = True
for _ in range(40):
    e = sorted(RNG.randint(-3, 3) for _ in range(3))
    g = rand_GLO() * Mat.diag(*[chi_pow(k) for k in e]) * rand_GLO()
    K, A, Kp, exps = snf_dvr(g)
    ok &= (exps == e) and in_GLO(K) and in_GLO(Kp) and (K * A * Kp == g)
    k1, a, k2, dexp = cartan(g)
    ok &= dexp == list(reversed(e)) and (k1 * a * k2 == g)
check("A3 from-scratch SNF/Cartan: reconstructs g, K,Kp in GL3(O), exact exponents", ok)

ok = True
for _ in range(40):
    g = rand_GLF()
    h = g * rand_GLO()
    ok &= hnf_local(g).key() == hnf_local(h).key()          # same lattice
    g2 = rand_GLF()
    same = in_GLO(g.inv() * g2)
    ok &= (hnf_local(g).key() == hnf_local(g2).key()) == same
check("A4 btlib hnf_local is a faithful lattice invariant (vs g^-1 h in GL3(O))", ok)

ok = True
for M in (H_GATE, S_GATE, R_GATE):
    ok &= M.conjT() * M == Mat.identity(3) and in_A(M)
for _ in range(20):
    g = rand_GLF()
    gi = g.inv()
    ok &= g * gi == Mat.identity(3)
    ok &= np.allclose(
        np.array([[complex(x.num.a + x.num.b * complex(-.5, 3**.5/2)) / x.den
                   for x in row] for row in (g * gi).m]), np.eye(3))
check("A5 gates unitary & in A; exact inverse matches numpy floats", ok)

# build a pool of pure-vertex bases for later random generation
for _ in range(6):
    for v in random_walk(4):
        if v.kind == "P":
            PURE_POOL.append(v.g)

# =============================================================================
print("=" * 78)
print("SECTION B: Lemma le:Aset")
print("=" * 78)

# B1: the unstated step:  M O^3 = O^3  <=>  M in GL_3(O)
ok = True
for _ in range(40):
    M = rand_GLO()
    ok &= is_O(M) and is_O(M.inv())                  # GL_3(O) => both integral
    ok &= hnf_local(M).key() == hnf_local(Mat.identity(3)).key()
    N = M * Mat.diag(1, 1, CHI)                      # integral, det not unit
    ok &= is_O(N) and not in_GLO(N)
    ok &= hnf_local(N).key() != hnf_local(Mat.identity(3)).key()
check("B1 M O^3 == O^3  <=>  M in GL_3(O)   (the implicit final step)", ok)

# B2: self-dual lattice basis  =>  g* g in GL_3(O)  (g in A), and conversely
ok = True
for _ in range(40):
    g = RNG.choice(PURE_POOL) * rand_GLO()           # basis of a self-dual lattice
    ok &= lat_eq(dualb(g), g)                        # self-duality, directly
    ok &= in_A(g)                                    # the lemma's conclusion
for _ in range(40):
    g = rand_A()                                     # generic element of A
    ok &= lat_eq(dualb(g), g)                        # converse (le:gramm_matrix)
check("B2 self-dual <=> Gram in GL_3(O) <=> g in A  (40+40 random bases)", ok)

# B3: witness that l(M)=0 alone does NOT put M in GL_3(O) (det can be non-unit);
#     relevant to how le:ispositive(5)/le:hecke_nghbs must (and do) argue.
M = Mat.diag(1, CHI, CHI)
check("B3 witness: l(M)==0 but M not in GL_3(O) (so l=0 alone is not unimodularity)",
      ell(M) == 0 and not in_GLO(M) and ell(M.inv()) == 2)

# =============================================================================
print("=" * 78)
print("SECTION C: Lemma le:diagonals")
print("=" * 78)

ok = True
for _ in range(80):
    g = rand_A()
    e = cartan_exps(g)
    n = e[0]
    ok &= (e == [n, 0, -n]) and n >= 0 and ell(g) == 2 * n
check("C1 g in A => Cartan diagonal is (chi^n, 1, chi^-n), n = l(g)/2 >= 0", ok)

# C2: the two unstated ingredients of the uniqueness argument
ok = True
for _ in range(40):
    g = rand_GLF()
    e = cartan_exps(g)
    ok &= cartan_exps(g.conjT()) == e                       # conj-transpose keeps exps
    ok &= cartan_exps(rand_GLO() * g.inv()) == [-x for x in reversed(e)]
check("C2 exps(g*) == exps(g);  exps(gamma g^-1) == reversed-negated exps(g)", ok)

g = Mat.diag(CHI, 1, 1)
check("C3 hypothesis g in A is needed: diag(chi,1,1) has exps (1,0,0), not (n,0,-n)",
      (not in_A(g)) and cartan_exps(g) == [1, 0, 0])

# =============================================================================
print("=" * 78)
print("SECTION D: Lemma le:ispositive (Appendix E)")
print("=" * 78)

ok = True
for _ in range(80):
    g = rand_GLF()
    ok &= ell(rand_GLO() * g * rand_GLO()) == ell(g)
check("D1 part(1): l(k g k') == l(g)", ok)

ok = True
for _ in range(80):
    a, b = rand_GLF(), rand_GLF()
    ok &= ell(a * b) <= ell(a) + ell(b)
check("D2 part(2): l(g1 g2) <= l(g1) + l(g2)", ok)

ok = True
for _ in range(80):
    g = rand_A()
    ok &= vpi(g.det()) == 0                       # det-valuation step
    ok &= ell(g) >= 0 and ell(g.inv()) == ell(g)
check("D3 part(3): v(det g)==0, l(g) >= 0, l(g^-1)==l(g) on A", ok)

ok = True
for _ in range(80):
    g, h = rand_A(), rand_A()
    m = g.inv() * h
    ok &= ell(m) >= 0 and ell(m.inv()) >= 0
check("D4 part(4): l(g^-1 h) >= 0 for g,h in A", ok)

ok = True
for _ in range(60):
    x = rand_GLF()
    e = cartan_exps(x)
    ok &= ell(x) == -2 * e[2] and ell(x.inv()) == 2 * e[0]
check("D5 l(x) == -2*(smallest Cartan exp); l(x^-1) == 2*(largest)  [all x]", ok)

ok = True
for _ in range(50):
    g = rand_A()
    ok &= d_til(g, g * rand_GLO()) == 0           # indiscernible on cosets
    h = rand_A()
    if not lat_eq(g, h):
        ok &= d_til(g, h) > 0                     # positivity
    m = rand_A()
    ok &= d_til(g, h) <= d_til(g, m) + d_til(m, h)         # triangle
    ok &= d_til(g, h) == d_til(h, g)                       # symmetry
    ok &= d_til(g * rand_GLO(), h * rand_GLO()) == d_til(g, h)  # well-defined
check("D6 part(5): d~ is a well-defined metric on A/GL_3(O) (50 random checks)", ok)

# =============================================================================
print("=" * 78)
print("SECTION E: Lemma le:hecke_nghbs (Appendix F)")
print("=" * 78)


def incl(g, h):
    """Lambda_g <= pi^-1 Lambda_h  iff  chi * h^-1 g integral."""
    return is_O((h.inv() * g) * CHI)


# E1: exhaustive at the origin: full distance-2 shell
shell1 = [canon(w) for w in neighbors(ORIGIN)]
shell2 = {}
for a in shell1:
    for w in neighbors(a):
        w = canon(w)
        if w.key() != ORIGIN.key():
            shell2[w.key()] = w
I3 = Mat.identity(3)
ok = len(shell2) == 12
for w in shell2.values():
    h = w.g
    l1, l2 = ell(h), ell(h.inv())
    ok &= (l1, l2) == (2, 2) and d_til(I3, h) == 2
    ok &= incl(I3, h) and incl(h, I3) and not lat_eq(I3, h)
    e = cartan_exps(h)
    ok &= e == [1, 0, -1]
check("E1 exhaustive shell: all 12 distance-2 vertices have (l,l)==(2,2), d~==2, "
      "inclusions, exps (1,0,-1)", ok, f"|shell2|={len(shell2)}")

# E2: forward direction at random spots:  tree-distance-2  =>  d~ == 2 & inclusions
ok, n_pairs = True, 0
while n_pairs < 40:
    p0, alt, p1 = random_walk(2)
    if p1.key() == p0.key():
        continue
    n_pairs += 1
    g, h = p0.g, p1.g
    PURE_POOL.append(g)
    # independent structural certificate of distance 2 (shared 1-simplices):
    lam = alt.g
    ok &= self_dual(g) and self_dual(h) and not lat_eq(g, h)
    ok &= lat_contains(lam, g) and lat_contains(g, dualb(lam))
    ok &= lat_contains(lam, h) and lat_contains(h, dualb(lam))
    ok &= d_til(g, h) == 2 and incl(g, h) and incl(h, g)
check("E2 distance-2 pairs (40 random): d~==2 and both inclusions", ok)

# E3: the case split: (4,0) is impossible; in fact l(g^-1 h)==l(h^-1 g) always on A
ok, pat = True, {}
for _ in range(400):
    g, h = rand_A(), rand_A()
    m = g.inv() * h
    l1, l2 = ell(m), ell(m.inv())
    pat[(l1, l2)] = pat.get((l1, l2), 0) + 1
    ok &= l1 == l2                       # forced by self-duality (n2 == 0)
    ok &= l1 % 2 == 0                    # parity: d~ is an even integer on A
check("E3 400 random A-pairs: l(g^-1h) == l(h^-1g) always (so (4,0) never occurs; "
      "d~ even)", ok, f"l-patterns seen: {sorted(pat)}")

# E4: the congruence-subgroup step:  x in GL3(O), d^-1 x^-1 d integral
#                                    =>  d^-1 x d integral
ok, n_tot = True, 0
for n in (1, 2):
    d = Mat.diag(chi_pow(n), 1, chi_pow(-n))
    di = d.inv()
    for _ in range(120):
        # generic w in GL_3(O) such that d^-1 w d is integral:
        L = Mat([[ONE if i == j else (rand_O() if i > j else ZERO)
                  for j in range(3)] for i in range(3)])
        U = Mat([[ONE, rand_O() * chi_pow(n), rand_O() * chi_pow(2 * n)],
                 [ZERO, ONE, rand_O() * chi_pow(n)],
                 [ZERO, ZERO, ONE]])
        w = L * Mat.diag(rand_unitO(), rand_unitO(), rand_unitO()) * U
        assert in_GLO(w) and is_O(di * w * d)
        n_tot += 1
        m = di * w.inv() * d
        ok &= is_O(m) and in_GLO(di * w * d) and in_GLO(m)
check("E4 group closure: d^-1 x^-1 d integral => d^-1 x d integral (and both in "
      "GL3(O); det stays a unit)", ok, f"{n_tot} cases, n=1,2")

# E5: converse direction, adversarially: pairs at distances 0,2,4 -- the
#     inclusions hold iff distance <= 2, and with L_g != L_h they force d~ == 2.
ok = True
for dist in (0, 2, 4):
    for _ in range(15):
        path = random_walk(dist)
        g, h = path[0].g, path[-1].g
        both = incl(g, h) and incl(h, g)
        if lat_eq(g, h):
            ok &= d_til(g, h) == 0
        elif both:
            ok &= d_til(g, h) == 2          # the lemma's converse
        if dist == 4 and d_til(g, h) == 4:
            ok &= not both                  # inclusions must fail at distance 4
check("E5 converse: inclusions + distinct lattices => d~ == 2; inclusions fail "
      "at distance 4", ok)

# E6: also certify the n1==1, n3==-1 step for (2,2)-pairs
ok = True
for _ in range(20):
    path = random_walk(2)
    if path[2].key() == path[0].key():
        continue
    e = cartan_exps(path[0].g.inv() * path[2].g)
    ok &= e[0] == 1 and e[2] == -1 and e[1] == 0
check("E6 (2,2) case: Cartan exps of g^-1 h are exactly (1,0,-1)", ok)

# =============================================================================
print("=" * 78)
print("SECTION F: Proposition pr:chain (Appendix D)")
print("=" * 78)

ok_cart = ok_gram = ok_sd = ok_consec = ok_dualid = True
n_tot = 0
for n, reps in ((2, 8), (3, 6)):
    for _ in range(reps):
        path = random_walk(2 * n)
        g, h = path[0].g, path[-1].g
        assert self_dual(g) and self_dual(h)
        n_tot += 1
        gh = g.inv() * h
        k1, a, k2, e = cartan(gh)
        ok_cart &= (e == [n, 0, -n]) and (k1 * a * k2 == gh)
        g1 = g * k1
        ok_cart &= lat_eq(g1, g) and lat_eq(g1 * a, h)      # Lam_h = g1 a O^3
        # minimality of n: inclusions hold at n, fail at n-1
        ok_cart &= lat_contains(g, h * chi_pow(n)) and lat_contains(h, g * chi_pow(n))
        ok_cart &= not (lat_contains(g, h * chi_pow(n - 1))
                        and lat_contains(h, g * chi_pow(n - 1)))
        # duality bridge (corrected version of the displayed identity):
        #   ((g1^-1 h)*)^-1 == l1^-1 (g1^-1 h) l2,  l1=(g1* g1)^-1, l2=(h* h)^-1
        l1m = (g1.conjT() * g1).inv()
        l2m = (h.conjT() * h).inv()
        m = g1.inv() * h
        ok_dualid &= (m.conjT().inv() == l1m.inv() * m * l2m) \
                     and in_GLO(l1m) and in_GLO(l2m)
        # Gram facts used for self-duality of the interpolants
        v1, v2, v3 = (g1.col(j) for j in range(3))
        ok_gram &= (vpi(herm(v2, v3)) >= n and vpi(herm(v3, v3)) >= 2 * n
                    and vpi(herm(v1, v3)) == 0 and vpi(herm(v2, v2)) == 0)
        # interpolating lattices  L_i = g1 diag(chi^i,1,chi^-i) O^3
        Ls = [g1 * Mat.diag(chi_pow(i), 1, chi_pow(-i)) for i in range(n + 1)]
        ok_sd &= all(self_dual(L) for L in Ls)
        ok_sd &= lat_eq(Ls[0], g) and lat_eq(Ls[n], h)
        for i in range(n):
            ok_consec &= (not lat_eq(Ls[i], Ls[i + 1]))
            ok_consec &= incl(Ls[i], Ls[i + 1]) and incl(Ls[i + 1], Ls[i])
            ok_consec &= d_til(Ls[i], Ls[i + 1]) == 2
check(f"F1 Cartan of g^-1 h == (n,0,-n) with n minimal; Lam_h == (g k1) a O^3 "
      f"({n_tot} pairs, n=2,3)", ok_cart)
check("F2 duality bridge ((g1^-1 h)*)^-1 == l1^-1 (g1^-1 h) l2 holds exactly "
      "(corrected form of the displayed equation)", ok_dualid)
check("F3 Gram facts: v<v2,v3> >= n, v<v3,v3> >= 2n, <v1,v3>,<v2,v2> units", ok_gram)
check("F4 interpolants L_i all self-dual; L_0=Lam_g, L_n=Lam_h", ok_sd)
check("F5 consecutive L_i,L_{i+1}: distinct, mutual pi^-1 inclusions, d~ == 2", ok_consec)

# the paper's displayed identity, literally:  ((g1^-1 h)*)^-1 = l1 ((k1 a k2)*)^-1 l2
g, h = PURE_POOL[0], PURE_POOL[1] if len(PURE_POOL) > 1 else rand_gamma()
# (re-derive a fresh distance-4 pair for the literal check)
path = random_walk(4)
g, h = path[0].g, path[-1].g
gh = g.inv() * h
k1, a, k2, e = cartan(gh)
g1 = g * k1
l1m = (g1.conjT() * g1).inv()
l2m = (h.conjT() * h).inv()
lhs = (g1.inv() * h).conjT().inv()
rhs_literal = l1m * gh.conjT().inv() * l2m
check("F6 (expository) the displayed identity is garbled: literal RHS != LHS, "
      "but both lie in GL3(O) (k1 a k2)-double-coset so the conclusion stands",
      not (lhs == rhs_literal) and cartan_exps(lhs) == cartan_exps(rhs_literal))

# =============================================================================
print("=" * 78)
print("SECTION G: Proposition prop:distance (Appendix G)")
print("=" * 78)

# G1: d~ equals the true graph distance (independent bidirectional BFS)
ok = True
for steps, reps in ((2, 6), (4, 4), (6, 3)):
    for _ in range(reps):
        path = random_walk(steps)
        g, h = path[0], path[-1]
        dg = graph_dist(g, h)
        ok &= dg == steps                       # walk really is geodesic (tree)
        ok &= d_til(g.g, h.g) == dg
check("G1 d~(g,h) == graph distance (BFS-certified pairs at distance 2,4,6)", ok)

# G2: the Lambda = Lam_g + Lam_h construction at distance-2 pairs
ok_chain = ok_alt = ok_same = True
n_pairs = 0
while n_pairs < 25:
    p0, alt, p1 = random_walk(2)
    if p1.key() == p0.key():
        continue
    n_pairs += 1
    g, h = p0.g, p1.g
    S = lat_sum(g, h)                            # Lam_g + Lam_h
    C = lat_cap(g, h)                            # Lam_g cap Lam_h
    # dual(L_g + L_h) == L_g cap L_h  and  dual(L_g cap L_h) == L_g + L_h
    ok_same &= lat_eq(dualb(S), C) and lat_eq(dualb(C), S)
    # the displayed inclusion chain, line by line:
    ok_chain &= lat_contains(C, g * CHI) and lat_contains(C, h * CHI)
    ok_chain &= lat_contains(dualb(S), C) and lat_contains(C, dualb(S))  # equality
    ok_chain &= lat_contains(g, dualb(S)) and lat_contains(h, dualb(S))
    ok_chain &= lat_contains(S, g) and lat_contains(S, h)
    ok_chain &= lat_contains(dualb(C), S)
    ok_chain &= lat_contains(g * chi_pow(-1), dualb(C))
    ok_chain &= lat_contains(h * chi_pow(-1), dualb(C))
    # the required sandwich for BOTH pure vertices (eq. before line 1541)
    for x in (g, h):
        ok_chain &= lat_contains(dualb(S), x * CHI)
        ok_chain &= lat_contains(x, dualb(S)) and lat_contains(S, x)
        ok_chain &= lat_contains(x * chi_pow(-1), S)
    # Lambda gives a genuine alternating vertex adjacent to both pure vertices
    ok_alt &= not lat_eq(S, C)
    ok_alt &= vpi((S.inv() * C).det()) == 2      # index [Lam : Lam#] = 9
    va = Vertex("A", S)
    nb0 = {w.key() for w in neighbors(p0)}
    nb1 = {w.key() for w in neighbors(p1)}
    ok_alt &= va.key() in nb0 and va.key() in nb1
check("G2 dual(Lg+Lh) == Lg cap Lh == dual; the 'or' candidates coincide", ok_same)
check("G3 displayed inclusion chain (lines 1542-1546) verified, all 25 pairs", ok_chain)
check("G4 Lambda = Lg+Lh is an alternating vertex adjacent to BOTH pure "
      "vertices => d(g v0, h v0) == 2", ok_alt)

# G5: reverse normalization step: distance 2n => pi^n Lam_h <= Lam_g <= pi^-n Lam_h
ok = True
for n in (1, 2, 3):
    for _ in range(4):
        path = random_walk(2 * n)
        g, h = path[0].g, path[-1].g
        ok &= lat_contains(g, h * chi_pow(n)) and lat_contains(h, g * chi_pow(n))
        ok &= ell(g.inv() * h) <= 2 * n and ell(h.inv() * g) <= 2 * n
check("G5 chain-of-1-simplices step: d==2n => pi^n L_h <= L_g <= pi^-n L_h "
      "=> both l's <= 2n  (n=1,2,3)", ok)

# G6: the MISSING direction d <= d~ (not in the appendix text): repaired via
#     pr:chain -- the interpolants give an explicit path of length 2n.
ok = True
for n, reps in ((2, 4), (3, 3)):
    for _ in range(reps):
        path = random_walk(2 * n)
        g, h = path[0].g, path[-1].g
        assert d_til(g, h) == 2 * n
        gh = g.inv() * h
        k1, a, k2, e = cartan(gh)
        g1 = g * k1
        Ls = [g1 * Mat.diag(chi_pow(i), 1, chi_pow(-i)) for i in range(n + 1)]
        # each consecutive pair shares the alternating vertex L_i + L_{i+1}:
        for i in range(n):
            S = lat_sum(Ls[i], Ls[i + 1])
            va = Vertex("A", S)
            ki = {w.key() for w in neighbors(Vertex("P", Ls[i]))}
            kj = {w.key() for w in neighbors(Vertex("P", Ls[i + 1]))}
            ok &= va.key() in ki and va.key() in kj
        # hence an explicit path of length 2n exists: d(g v0, h v0) <= d~(g,h)
check("G6 repair of the missing direction d <= d~: pr:chain interpolants give "
      "an explicit path of length d~(g,h)", ok)

# =============================================================================
print("=" * 78)
if FAILURES:
    print(f"OVERALL: FAIL ({len(FAILURES)}):")
    for f in FAILURES:
        print("   -", f)
    sys.exit(1)
print("OVERALL: PASS -- all computational claims of le:Aset, le:diagonals,")
print("le:ispositive, le:hecke_nghbs, pr:chain, prop:distance verified.")
print()
print("NOTE (textual, not computational): the appendix proof of prop:distance")
print("(lines 1516-1535) derives d~ <= d twice and never derives d <= d~;")
print("section G6 above verifies that the missing direction follows from")
print("pr:chain + le:hecke_nghbs + the d~=2 <=> d=2 claim, so the gap is fixable.")
