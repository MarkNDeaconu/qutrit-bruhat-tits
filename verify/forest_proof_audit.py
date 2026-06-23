#!/usr/bin/env python
"""
forest_proof_audit.py -- adversarial audit of Appendix B (Proposition pr:forest,
"B is a forest") of paper.tex.

WHAT THIS SCRIPT DOES
=====================
Part 0: sanity-checks btlib.py primitives INDEPENDENTLY (floats, brute force,
        semantic re-derivations) before relying on them.
Part 1: demonstrates the flaw in the Appendix-B proof:
        the chain relations (B5)/(B6) and the 5-term nested display
            pi^n A_n  <  A_0  <  S_0  <  A_0#  <  (pi^n A_n)#
        are derived WITHOUT using the closing edge of the loop, hence they hold
        along any geodesic PATH in the (true) tree.  We build such paths
        (n = 1 and n = 2) and verify:
          (a) every displayed row relation used by the paper holds  [premises OK]
          (b) the 5-term display is properly nested with 5 distinct pi-classes
          (c) the pi-closure of the display is NOT totally ordered, so it is NOT
              a lattice chain in the sense of Definition de:lattice_chain
              (which requires stability under L -> pi L);  in particular the
              claimed contradiction with pr:simplex_classification (rank > 3
              chain = 2-simplex) does not arise.  The same inference applied to
              a path would "prove" geodesics of length >= 3 do not exist.
Part 2: verifies every step of the candidate repair:
        (i)   two pure vertices at distance 2 have a UNIQUE common alternating
              neighbour, with big lattice S + S' and small lattice S /\ S'.
        (ii)  STRONG forced first step: for pure v at distance 2k from the
              origin (k>=1), exactly one alternating neighbour x of v has
              small(x) >= pi^k O^3; it is the unique neighbour at distance
              2k-1, and big(x) = Lambda_v + pi^(k-1) O^3.
        (iii) Lemma C: alternating x at distance 2n-1 => small(x) >= pi^n O^3.
        (iv)  Case A of the cycle argument: an alternating vertex NEVER has two
              pure neighbours at distance d(x)-1 (neighbour-distance multiset
              is {d-1 once, d+1 thrice}).
        (v)   local cycle-freeness of the BFS ball (consistency).

Exits nonzero on any FAIL.  Deterministic (seeded RNG).
"""

import sys, os, random, itertools, cmath

sys.path.insert(0, "/Users/markdeaconu/projects/qutrits_v2")
from btlib import (Zw, Fw, Mat, ZERO, ONE, OMEGA, CHI, CHI_INV, chi_pow,
                   H_GATE, S_GATE, R_GATE, monomial_matrices, in_Gamma,
                   lattice_contains, lattice_eq, dual_basis, gram,
                   is_self_dual, lattice_index_log3, hnf_local, Vertex,
                   ORIGIN, neighbors, neighbors_of_pure,
                   neighbors_of_alternating, ell, d_tilde,
                   F3_PLANES, _complete_to_unimodular, _v3)

random.seed(12345)

FAILURES = []
def check(name, ok, extra=""):
    line = f"{'PASS' if ok else 'FAIL'}: {name}" + (f"  [{extra}]" if extra else "")
    print(line)
    if not ok:
        FAILURES.append(name)

# ---------------------------------------------------------------------------
# helpers (independent of btlib's tree machinery)
# ---------------------------------------------------------------------------

W3 = cmath.exp(2j * cmath.pi / 3)

def fw_to_c(x: Fw):
    return (x.num.a + x.num.b * W3) / x.den

def mat_to_c(M: Mat):
    return [[fw_to_c(M.m[i][j]) for j in range(M.n)] for i in range(M.n)]

def vdet(M: Mat):
    return M.det().v_pi()

def lattice_sum(A: Mat, B: Mat) -> Mat:
    """Basis of A*O^3 + B*O^3 (column reduction over the DVR O_pi)."""
    cols = [A.col(j) for j in range(3)] + [B.col(j) for j in range(3)]
    basis = []
    rem = cols
    for r in range(3):
        best, bv = None, None
        for c in rem:
            v = c[r].v_pi()
            if not c[r].is_zero() and (bv is None or v < bv):
                best, bv = c, v
        assert best is not None, "rank defect in lattice_sum"
        new_rem = []
        for c in rem:
            if c is best:
                continue
            if not c[r].is_zero():
                f = c[r] / best[r]           # in O_pi by pivot minimality
                assert f.v_pi() >= 0
                c = [c[i] - f * best[i] for i in range(3)]
            new_rem.append(c)
        basis.append(best)
        rem = new_rem
    for c in rem:
        assert all(x.is_zero() for x in c), "non-zero residual column"
    return Mat([[basis[j][i] for j in range(3)] for i in range(3)])

def lattice_meet(A: Mat, B: Mat) -> Mat:
    """A /\\ B = (A# + B#)#."""
    return dual_basis(lattice_sum(dual_basis(A), dual_basis(B)))

def pi_equivalent(A: Mat, B: Mat):
    dv = vdet(A) - vdet(B)
    if dv % 3:
        return False
    return lattice_eq(A, B * chi_pow(dv // 3))

def comparable(A: Mat, B: Mat):
    return lattice_contains(A, B) or lattice_contains(B, A)

def proper_chain(*Ms):
    """Each M_{i} properly contained in M_{i+1} (as lattices)."""
    for X, Y in zip(Ms, Ms[1:]):
        if not lattice_contains(Y, X):
            return False
        if lattice_eq(X, Y):
            return False
    return True

# ===========================================================================
print("=" * 76)
print("PART 0: independent sanity checks of btlib primitives")
print("=" * 76)

# 0a. Fw arithmetic vs complex floats -----------------------------------------
ok = True
for _ in range(300):
    a = Fw(Zw(random.randint(-9, 9), random.randint(-9, 9)),
           random.randint(1, 12))
    b = Fw(Zw(random.randint(-9, 9), random.randint(-9, 9)),
           random.randint(1, 12))
    if b.is_zero():
        continue
    for got, want in [(a + b, fw_to_c(a) + fw_to_c(b)),
                      (a * b, fw_to_c(a) * fw_to_c(b)),
                      (a - b, fw_to_c(a) - fw_to_c(b)),
                      (a / b, fw_to_c(a) / fw_to_c(b)),
                      (a.conj(), fw_to_c(a).conjugate())]:
        if abs(fw_to_c(got) - want) > 1e-9:
            ok = False
check("0a field arithmetic Fw matches complex floats", ok)

# 0b. v_chi(x) == v_3(Norm(x)) (chi is ramified, N(chi)=3) --------------------
ok = True
for _ in range(300):
    z = Zw(random.randint(-40, 40), random.randint(-40, 40))
    if z.is_zero():
        continue
    if z.v_chi() != _v3(z.norm()):
        ok = False
check("0b v_chi(x) == v_3(N(x)) on random Z[w] elements", ok)
check("0b conventions: v_pi(3)=2, v_pi(chi)=1, conj(chi)=-w^2 chi, w==1 mod chi",
      Fw(3).v_pi() == 2 and CHI.v_pi() == 1
      and CHI.conj() == -(OMEGA * OMEGA) * CHI
      and (OMEGA - ONE).v_pi() >= 1)

# 0c. dual_basis: semantic check of the duality pairing -----------------------
def rand_basis():
    while True:
        M = Mat([[Fw(Zw(random.randint(-3, 3), random.randint(-3, 3)),
                     random.choice([1, 1, 3])) for _ in range(3)]
                 for _ in range(3)])
        if not M.det().is_zero():
            return M

ok = True
for _ in range(25):
    g = rand_basis()
    D = dual_basis(g)
    # <x,y> = sum x_i conj(y_i) must be integral for x in D*O^3, y in g*O^3
    for _ in range(10):
        c = [random.randint(-2, 2) for _ in range(3)]
        cp = [random.randint(-2, 2) for _ in range(3)]
        x = [sum((D.m[i][j] * c[j] for j in range(3)), ZERO) for i in range(3)]
        y = [sum((g.m[i][j] * cp[j] for j in range(3)), ZERO) for i in range(3)]
        pair = sum((x[i] * y[i].conj() for i in range(3)), ZERO)
        if not pair.is_zero() and pair.v_pi() < 0:
            ok = False
    # maximality: D_j / chi pairs non-integrally with some basis vector of g
    for j in range(3):
        x = [D.m[i][j] * CHI_INV for i in range(3)]
        good = False
        for jj in range(3):
            pair = sum((x[i] * g.m[i][jj].conj() for i in range(3)), ZERO)
            if (not pair.is_zero()) and pair.v_pi() < 0:
                good = True
        if not good:
            ok = False
check("0c dual_basis satisfies (and maximises) the duality pairing", ok)

# 0d. is_self_dual (Gram criterion) == lattice equality with the dual ---------
ok = True
tests = [Mat.identity(), H_GATE, Mat.identity() * CHI,
         Mat.diag(CHI, 1, CHI_INV)]   # diag(chi, 1, chi^-1): NOT self-dual
for _ in range(40):
    tests.append(rand_basis())
for g in tests:
    if is_self_dual(g) != lattice_eq(g, dual_basis(g)):
        ok = False
check("0d is_self_dual(g) <=> g O^3 == dual lattice, on 44 lattices", ok)
check("0d H is unitary and in Gamma; O^3 self-dual; chi O^3 not self-dual",
      H_GATE.is_unitary() and in_Gamma(H_GATE)
      and is_self_dual(Mat.identity()) and not is_self_dual(Mat.identity() * CHI))

# 0e. hnf_local canonical under GL_3(O) ---------------------------------------
def rand_unimodular():
    """Random element of GL_3(O_pi) (product of elementary matrices/units)."""
    U = Mat.identity()
    units = [ONE, -ONE, OMEGA, -OMEGA, OMEGA * OMEGA]
    for _ in range(6):
        i, j = random.sample(range(3), 2)
        E = Mat.identity()
        E.m[i][j] = random.choice([ONE, -ONE, OMEGA, CHI, Fw(2)])
        U = U * E
    D = Mat.diag(*[random.choice(units) for _ in range(3)])
    U = U * D
    assert U.in_GL_O()
    return U

ok = True
for _ in range(30):
    g = rand_basis()
    U = rand_unimodular()
    if hnf_local(g).key() != hnf_local(g * U).key():
        ok = False
    if not lattice_eq(hnf_local(g), g):
        ok = False
check("0e hnf_local is a canonical form on lattices g*GL_3(O)", ok)

# 0f. neighbor functions vs BRUTE FORCE over all 13 index-3 sublattices -------
def index3_sublattices(g: Mat):
    out = []
    for w1, w2, _span in F3_PLANES:
        U = _complete_to_unimodular(w1, w2)
        out.append(g * U * Mat.diag(1, 1, CHI))
    return out

def brute_alt_neighbors(v: Vertex):
    """Alternating vertices {A, A#} with A < Lambda_v < A#, pi A# < A."""
    assert v.kind == "P"
    keys = set()
    for A in index3_sublattices(v.g):
        big = dual_basis(A)
        assert lattice_contains(v.g, A) and lattice_index_log3(v.g, A) == 1
        if not lattice_contains(big, v.g):
            continue
        if not lattice_contains(A, big * CHI):    # pi A# subset A
            continue
        keys.add(Vertex("A", big).key())
    return keys

def brute_pure_neighbors(v: Vertex):
    """Self-dual lattices S with A < S < A#  (index-3 sublattices of big)."""
    assert v.kind == "A"
    keys = set()
    for S in index3_sublattices(v.g):
        if lattice_eq(S, dual_basis(S)):
            keys.add(Vertex("P", S).key())
    return keys

ok = True
test_pure = [ORIGIN, Vertex("P", H_GATE)]
for v in test_pure:
    nb = {x.key() for x in neighbors_of_pure(v)}
    bf = brute_alt_neighbors(v)
    if not (len(nb) == 4 and nb == bf):
        ok = False
check("0f neighbors_of_pure == brute force over all index-3 sublattices", ok)

ok = True
test_alt = neighbors_of_pure(ORIGIN)[:2] + neighbors_of_pure(Vertex("P", H_GATE))[:1]
for v in test_alt:
    nb = {x.key() for x in neighbors_of_alternating(v)}
    bf = brute_pure_neighbors(v)
    if not (len(nb) == 4 and nb == bf):
        ok = False
check("0f neighbors_of_alternating == brute force (self-dual index-3 sublattices)",
      ok)

# 0g. BFS ball; d_tilde (metric from l) == BFS graph distance -----------------
RADIUS = 6
dist = {ORIGIN.key(): 0}
vert = {ORIGIN.key(): ORIGIN}
adj = {}
frontier = [ORIGIN]
for r in range(1, RADIUS + 1):
    nxt = []
    for v in frontier:
        nbs = neighbors(v)
        adj[v.key()] = [w.key() for w in nbs]
        for w in nbs:
            wk = w.key()
            if wk not in dist:
                dist[wk] = r
                vert[wk] = w
                nxt.append(w)
    frontier = nxt
sizes = [sum(1 for k in dist if dist[k] == r) for r in range(RADIUS + 1)]
check("0g BFS ball sizes 1,4,12,36,108,324,972 ((4,4)-biregular tree counts)",
      sizes == [1, 4, 12, 36, 108, 324, 972], f"sizes={sizes}")

ok = True
for k, d in dist.items():
    v = vert[k]
    if v.kind == "P" and d <= RADIUS - 1:
        if d_tilde(Mat.identity(), v.g) != d:
            ok = False
check("0g d_tilde(I,g) == BFS distance for all pure vertices to radius 5", ok)

# every vertex of the interior has exactly 4 neighbours; edges go level +-1;
# no non-tree edge (local cycle-freeness of B)
ok = True
tree_ok = True
for k, nbs in adj.items():
    if len(nbs) != 4 or len(set(nbs)) != 4:
        ok = False
    d = dist[k]
    closer = [w for w in nbs if dist.get(w, 10 ** 9) == d - 1]
    same = [w for w in nbs if dist.get(w, 10 ** 9) == d]
    if same:
        tree_ok = False
    if d > 0 and len(closer) != 1:
        tree_ok = False
check("0g every interior vertex has 4 distinct neighbours", ok)
check("0g ball radius 6 is a tree locally: unique closer neighbour, no cross edges",
      tree_ok)

# ===========================================================================
print()
print("=" * 76)
print("PART 1: the Appendix-B proof of pr:forest -- demonstration of the flaw")
print("=" * 76)
print("""
The proof derives, from the loop rows, only:
   (B5)  pi^n A_n  subset  A_0          (left diagonal,  rows 2..2n+1)
   (B6)  A_0#      subset  pi^-n A_n#   (right diagonal, dual of B5)
   plus  A_0 subset S_0 subset A_0#     (row 1)
and concludes the display  pi^n A_n < A_0 < S_0 < A_0# < (pi^n A_n)#
"is a self-dual lattice chain of rank greater than 3".
The closing edge (row 2n+2) is NEVER used.  So the same premises hold for a
geodesic PATH  S_0 - A_0 - S_1 - ... - S_n - A_n  in the true tree.  We build
such paths and show the premises hold but the conclusion object violates
Definition de:lattice_chain (no pi-stable total order exists).
""")

# build a geodesic path S0 - A0 - S1 - A1 - S2 - A2 from the origin ----------
S0v = ORIGIN
A0v = neighbors(S0v)[0]
S1v = [w for w in neighbors(A0v) if w != S0v][0]
A1v = [w for w in neighbors(S1v) if w != A0v][0]
S2v = [w for w in neighbors(A1v) if w != S1v][0]
A2v = [w for w in neighbors(S2v) if w != A1v][0]
path = [S0v, A0v, S1v, A1v, S2v, A2v]
check("1a path vertices pairwise distinct and at BFS distances 0..5",
      len({v.key() for v in path}) == 6
      and [dist[v.key()] for v in path] == [0, 1, 2, 3, 4, 5])

I3 = Mat.identity()
L0, a0 = A0v.g, dual_basis(A0v.g)       # A_0# and A_0
s1 = S1v.g
L1, a1 = A1v.g, dual_basis(A1v.g)       # A_1# and A_1
s2 = S2v.g
L2, a2 = A2v.g, dual_basis(A2v.g)       # A_2# and A_2

# (a) all row relations used by the paper hold with these representatives ----
rows_ok = (
    proper_chain(a0, I3, L0)            # row 1:  A_0 < S_0 < A_0#
    and proper_chain(a0, s1, L0)        # row 2:  A_0 < S_1 < A_0#
    and proper_chain(a1, s1, L1)        # row 3:  A_1 < S_1 < A_1#
    and proper_chain(a1, s2, L1)        # row 4:  A_1 < S_2 < A_1#
    and proper_chain(a2, s2, L2)        # row 5:  A_2 < S_2 < A_2#
    and lattice_contains(a0, L0 * CHI)  # pi A_0# < A_0  (0-simplex chain)
    and lattice_contains(a1, L1 * CHI)
    and lattice_contains(a0, s1 * CHI)  # pi S_1 < A_0   (row-2 diagonal)
    and lattice_contains(a1, s2 * CHI)
)
check("1a every displayed row relation holds on the geodesic path", rows_ok)

tel_ok = (lattice_contains(a0, a1 * CHI)                  # pi A_1 < A_0
          and lattice_contains(a1, a2 * CHI)
          and lattice_contains(a0, a2 * chi_pow(2))       # pi^2 A_2 < A_0
          and lattice_contains(dual_basis(a1) * CHI_INV, L0))   # A_0# < pi^-1 A_1#
check("1a telescoped relations (B5)/(B6) hold:  pi^n A_n < A_0,  A_0# < pi^-n A_n#",
      tel_ok)

# (b) the 5-term display is properly nested, with 5 distinct pi-classes ------
for n, (an, Ln) in [(1, (a1, L1)), (2, (a2, L2))]:
    disp = [an * chi_pow(n), a0, I3, L0, Ln * chi_pow(-n)]
    names = [f"pi^{n}A_{n}", "A_0", "S_0", "A_0#", f"pi^-{n}A_{n}#"]
    check(f"1b n={n}: display properly nested: " + " < ".join(names),
          proper_chain(*disp))
    distinct = all(not pi_equivalent(disp[i], disp[j])
                   for i in range(5) for j in range(i + 1, 5))
    check(f"1b n={n}: the 5 displayed lattices lie in 5 DISTINCT pi-classes",
          distinct)
    # self-dual display: dual of pi^n A_n is pi^-n A_n# (up to units)
    check(f"1b n={n}: display is 'self-dual' as a finite nested sequence",
          lattice_eq(dual_basis(disp[0]), disp[4])
          and lattice_eq(dual_basis(disp[1]), disp[3])
          and lattice_eq(dual_basis(disp[2]), disp[2]))

    # (c) ... but it is NOT a lattice chain per Definition de:lattice_chain:
    # the pi-closure is not totally ordered.
    # (c1) the wrap-around inclusion pi * (pi^n A_n)# subset pi^n A_n, i.e.
    #      A_n# subset pi^(2n-1) A_n, is impossible:
    wrap = lattice_contains(an * chi_pow(2 * n - 1), Ln)
    check(f"1c n={n}: wrap-around A_{n}# subset pi^(2n-1)A_{n} FAILS (as predicted)",
          not wrap,
          f"v(det) -1 vs {6 * n - 2} makes it impossible")
    # (c2) A_0 and A_n are pi-translates of displayed lattices with EQUAL
    #      determinant valuation; total order would force A_0 == A_n:
    check(f"1c n={n}: v(det A_0) == v(det A_{n}) == 1 but A_0 != A_{n}",
          vdet(a0) == 1 and vdet(an) == 1 and not lattice_eq(a0, an))
    check(f"1c n={n}: A_0 and A_{n} are INCOMPARABLE -> pi-closure of the "
          f"display is not totally ordered -> NOT a lattice chain "
          f"(Def. de:lattice_chain)",
          not comparable(a0, an))
    # (c3) exhaustive comparability scan over pi-translates in a window
    translates = []
    for nm, M in zip(names, disp):
        for i in range(-2, 3):
            translates.append((f"pi^{i}*{nm}", M * chi_pow(i)))
    bad = []
    for (n1, M1), (n2, M2) in itertools.combinations(translates, 2):
        if pi_equivalent(M1, M2):
            continue
        if not comparable(M1, M2):
            bad.append((n1, n2))
    check(f"1c n={n}: incomparable pi-translate pairs exist ({len(bad)} found, "
          f"e.g. {bad[0] if bad else None})", len(bad) > 0)

print("""
=> The premises of the final inference of Appendix B hold along plain geodesic
   paths of the true tree, while the conclusion ("a self-dual lattice chain of
   rank > 3 exists") is false there.  The inference is therefore INVALID:
   nestedness of the 5 displayed lattices does not produce an (admissible,
   i.e. pi-stable) lattice chain, and no contradiction with
   pr:simplex_classification is obtained.  The proof of pr:forest is broken.
""")

# ===========================================================================
print("=" * 76)
print("PART 2: verification of the candidate repair")
print("=" * 76)

# (i) unique common alternating neighbour of pure vertices at distance 2;
#     big = S + S', small = S /\ S' ------------------------------------------
alt_keys = [k for k in adj if vert[k].kind == "A"]          # dist <= 5
ok_gap, ok_sum = True, True
pair_to_alts = {}
for k in alt_keys:
    x = vert[k]
    big, small = x.g, dual_basis(x.g)
    pure_nbs = [vert[w] for w in adj[k]]
    for S in pure_nbs:                       # canonical gap: A < S < A#
        if not (lattice_contains(big, S.g) and lattice_contains(S.g, small)):
            ok_gap = False
    for Su, Sv in itertools.combinations(pure_nbs, 2):
        key = tuple(sorted([Su.key(), Sv.key()]))
        pair_to_alts.setdefault(key, []).append(k)
        if not lattice_eq(lattice_sum(Su.g, Sv.g), big):
            ok_sum = False
        if not lattice_eq(lattice_meet(Su.g, Sv.g), small):
            ok_sum = False
check("2(i) all pure neighbours of an alternating vertex lie in its canonical "
      "gap (A, A#)", ok_gap)
check("2(i) for every pair S != S' of pure neighbours: S+S' == A#  and  "
      "S/\\S' == A  (checked for all alternating vertices to radius 5)", ok_sum)
ok = True
n_pairs = 0
for (ku, kv), alts in pair_to_alts.items():
    if dist[ku] <= RADIUS - 2 and dist[kv] <= RADIUS - 2:
        n_pairs += 1
        if len(alts) != 1:
            ok = False
check("2(i) every pure pair at distance 2 has a UNIQUE common alternating "
      f"neighbour ({n_pairs} pairs checked)", ok)

# (ii) strong forced first step  ----------------------------------------------
ok_unique, ok_dist, ok_formula = True, True, True
n_checked = 0
for k, d in dist.items():
    v = vert[k]
    if v.kind != "P" or d == 0 or d > RADIUS - 1 or d % 2:
        continue
    kk = d // 2
    n_checked += 1
    hits = []
    closer = []
    for w in adj[k]:
        x = vert[w]
        small = dual_basis(x.g)
        if lattice_contains(small, I3 * chi_pow(kk)):     # pi^k O^3 <= A_x
            hits.append(x)
        if dist[w] == d - 1:
            closer.append(x)
    if len(hits) != 1:
        ok_unique = False
    if len(closer) != 1 or hits[0].key() != closer[0].key():
        ok_dist = False
    if not lattice_eq(hits[0].g, lattice_sum(v.g, I3 * chi_pow(kk - 1))):
        ok_formula = False
check(f"2(ii) pure v at distance 2k: EXACTLY ONE alternating neighbour x with "
      f"pi^k O^3 <= A_x ({n_checked} vertices, k=1,2)", ok_unique)
check("2(ii) that x is precisely the unique neighbour at distance 2k-1",
      ok_dist)
check("2(ii) big(x) == Lambda_v + pi^(k-1) O^3   (forced first-step formula)",
      ok_formula)

# long-range spot checks (k = 3, 4) using only the d~ metric ------------------
def random_pure_at(k_target, tries=400):
    reps = monomial_matrices()
    for _ in range(tries):
        g = Mat.identity()
        for _ in range(k_target):
            g = g * (random.choice(reps) * H_GATE)
        if d_tilde(I3, g) == 2 * k_target and is_self_dual(g):
            return Vertex("P", g)
    return None

for kk in (3, 4):
    v = random_pure_at(kk)
    if v is None:
        check(f"2(ii) long-range k={kk}: vertex construction", False)
        continue
    nbs = neighbors(v)
    hits = [x for x in nbs
            if lattice_contains(dual_basis(x.g), I3 * chi_pow(kk))]
    okf = (len(hits) == 1
           and lattice_eq(hits[0].g, lattice_sum(v.g, I3 * chi_pow(kk - 1))))
    # the hit has a pure neighbour at distance 2k-2; the other three have all
    # pure neighbours at distance >= 2k
    descend = any(d_tilde(I3, T.g) == 2 * kk - 2
                  for T in neighbors(hits[0]))
    others_ok = all(min(d_tilde(I3, T.g) for T in neighbors(x)) >= 2 * kk
                    for x in nbs if x.key() != hits[0].key())
    check(f"2(ii) long-range k={kk}: unique hit, formula, descent, and the "
          f"3 other neighbours do not descend", okf and descend and others_ok)

# (iii) Lemma C ---------------------------------------------------------------
ok = True
n_checked = 0
for k, d in dist.items():
    v = vert[k]
    if v.kind != "A":
        continue
    nn = (d + 1) // 2
    n_checked += 1
    if not lattice_contains(dual_basis(v.g), I3 * chi_pow(nn)):
        ok = False
check(f"2(iii) alternating x at distance 2n-1: pi^n O^3 <= A_x "
      f"({n_checked} vertices)", ok)

# (iv) Case A: an alternating vertex never has two closer pure neighbours ----
ok = True
n_checked = 0
for k in alt_keys:
    d = dist[k]
    ds = sorted(dist[w] for w in adj[k])
    n_checked += 1
    if ds != [d - 1] + [d + 1] * 3:
        ok = False
check(f"2(iv) every alternating vertex (<= radius 5) has pure-neighbour "
      f"distance multiset {{d-1, d+1, d+1, d+1}} ({n_checked} checked)", ok)

# (v) bipartite structure ------------------------------------------------------
ok = all(vert[w].kind != vert[k].kind for k in adj for w in adj[k])
check("2(v) graph is bipartite: every edge joins a pure and an alternating "
      "vertex", ok)

# ===========================================================================
print()
print("=" * 76)
if FAILURES:
    print(f"RESULT: {len(FAILURES)} FAILURE(S):")
    for f in FAILURES:
        print("  -", f)
    sys.exit(1)
print("RESULT: ALL CHECKS PASSED")
print("""Summary:
 * Part 1 demonstrates the Appendix-B proof of pr:forest is INVALID as written:
   its final inference (5 nested lattices => self-dual lattice chain of rank>3)
   is false -- the exhibited object is never a lattice chain in the sense of
   Definition de:lattice_chain, and the same premises hold on plain geodesics.
 * Part 2 verifies every computational ingredient of the proposed repair
   (forced first step / unique closer neighbour), which yields a correct proof.
""")
sys.exit(0)
