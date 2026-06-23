#!/usr/bin/env python
"""
v06_tree.py -- Verification of the global tree structure of the Bruhat-Tits
building B for U_3(Z[1/3, omega])  ("Buildings for Synthesis with Clifford+R").

Claims verified (paper Prop pr:simplex_classification, pr:forest, connectedness,
Prop prop:distance, Lemma le:hecke_nghbs neighborhood / Figure 1):

  PART 0  btlib primitives sanity-checked independently
          (v_pi vs 3-adic valuation of the field norm; numeric complex
           embedding vs numpy; residue map; hnf_local canonicity & column
           equivalence; self-duality via Gram vs via L == L^dual;
           neighbors() re-derived by brute-force lattice enumeration).
  PART 1  BFS from origin to radius 7:
          (a) every vertex (radius <= 6) has exactly 4 distinct neighbors;
              pure <-> alternating; symmetry w in N(v) <=> v in N(w) on
              ~200 sampled edges;
          (b) bipartite: BFS levels strictly alternate P/A;
          (c) sphere sizes exactly [1,4,12,36,108,324,972,2916];
          (d) zero non-tree edges encountered (empirical tree/forest);
          (e) canonical-key integrity under pi-scaling (chi, chi^3) and
              random right GL_3(O) changes of basis (~100 samples).
  PART 2  uniqueness of middle vertex for ~50 distance-2 pure pairs; the
          middle big lattice equals L_g + L_h, its dual the intersection
          (paper lines ~1536-1546).
  PART 3  d~(g,h) == BFS graph distance on ~100 random pure pairs (ball r<=6).
  PART 4  U_3(O_F) (1296 monomials) fixes the origin, permutes S_0 (the 12
          pure distance-2 vertices) and the 4 alternating distance-1 vertices,
          acts transitively on S_0; e_1 = H e_0 in S_0, d(e_0, e_1) = 2.
  PART 5  Gamma-action well-defined: ~20 random words in {H,S,R} acting on
          ~10 ball vertices give valid vertices of the same kind, commute with
          neighbors() (equivariance), and preserve d~ on sampled pure pairs.

Deterministic (seeded).  Exits nonzero on any failure.
"""
import sys, random, itertools
from math import inf

sys.path.insert(0, "/Users/markdeaconu/projects/qutrits_v2")
import numpy as np
import btlib
from btlib import (Zw, Fw, Mat, ZERO, ONE, OMEGA, CHI, CHI_INV, I_SQRT3,
                   chi_pow, ell, d_tilde, H_GATE, S_GATE, R_GATE,
                   monomial_matrices, in_Gamma, lattice_contains, lattice_eq,
                   dual_basis, gram, is_self_dual, lattice_index_log3,
                   hnf_local, Vertex, ORIGIN, neighbors, F3_LINES, F3_PLANES,
                   _complete_to_unimodular)

random.seed(20260610)

FAILURES = []
def report(name, ok, detail=""):
    tag = "PASS" if ok else "FAIL"
    print(f"[{tag}] {name}" + (f"  -- {detail}" if detail else ""))
    if not ok:
        FAILURES.append(name)

# ============================================================================
# Independent primitives (NOT relying on btlib internals)
# ============================================================================

def v3_int(n):
    n = abs(int(n)); assert n > 0
    v = 0
    while n % 3 == 0:
        n //= 3; v += 1
    return v

def v_pi_indep(x: Fw):
    """Independent pi-adic valuation: chi is the unique prime over 3 in
    Z[omega] and N(chi)=3, so v_chi(z) = v_3(N(z)) for z in Z[omega],
    N(a+bw) = a^2 - ab + b^2.  Then v_pi(num/den) = v_3(N(num)) - 2 v_3(den)."""
    if x.num.is_zero():
        return inf
    n = x.num.a * x.num.a - x.num.a * x.num.b + x.num.b * x.num.b
    return v3_int(n) - 2 * v3_int(x.den) if x.den % 3 == 0 else v3_int(n)

OMEGA_C = np.exp(2j * np.pi / 3)
def fw_c(x: Fw):
    return (x.num.a + x.num.b * OMEGA_C) / x.den
def mat_c(M: Mat):
    return np.array([[fw_c(M.m[i][j]) for j in range(M.n)] for i in range(M.n)])

def rand_fw(maxc=9, dens=(1, 1, 1, 3, 9)):
    return Fw(Zw(random.randint(-maxc, maxc), random.randint(-maxc, maxc)),
              random.choice(dens))

def rand_gl3_O(maxc=2):
    """Random element of GL_3(O_pi) with small Z[omega] entries."""
    while True:
        M = Mat([[Fw(Zw(random.randint(-maxc, maxc),
                        random.randint(-maxc, maxc))) for _ in range(3)]
                 for _ in range(3)])
        d = M.det()
        if not d.is_zero() and d.v_pi() == 0:
            return M

def in_GL_O_indep(M: Mat):
    """Independent GL_3(O_pi) membership using v_pi_indep only."""
    if any(v_pi_indep(M.m[i][j]) < 0 for i in range(3) for j in range(3)):
        return False
    return v_pi_indep(M.det()) == 0

def is_self_dual_indep(g: Mat):
    """Lemma le:gramm_matrix, recomputed with independent valuations."""
    return in_GL_O_indep(g.conjT() * g)

def span_basis(cols):
    """Independent column-HNF over the DVR O_pi for a 3 x n column list.
    Returns a Mat whose columns form a triangular basis of the O_pi-span.
    (Re-implemented from scratch; does NOT call btlib.hnf_local.)"""
    cols = [list(c) for c in cols]
    n = len(cols)
    assert n >= 3
    for r in range(3):
        best, bv = None, None
        for j in range(r, n):
            v = cols[j][r].v_pi()
            if v != inf and (bv is None or v < bv):
                best, bv = j, v
        assert best is not None, "rank < 3 in span_basis"
        cols[r], cols[best] = cols[best], cols[r]
        u = chi_pow(bv) / cols[r][r]              # unit of O_pi
        cols[r] = [x * u for x in cols[r]]
        for j in range(r + 1, n):
            c = cols[j][r] * chi_pow(-bv)         # in O_pi by minimality
            if c.is_zero():
                continue
            cols[j] = [cols[j][i] - c * cols[r][i] for i in range(3)]
    for j in range(3, n):                          # tail must be zero columns
        assert all(x.is_zero() for x in cols[j]), "span_basis leftover"
    return Mat([[cols[j][i] for j in range(3)] for i in range(3)])

def mat_cols(g: Mat):
    return [g.col(j) for j in range(3)]

def lattice_sum(g: Mat, h: Mat):
    return span_basis(mat_cols(g) + mat_cols(h))

def lattice_intersection(g: Mat, h: Mat):
    """L_g cap L_h = dual(dual L_g + dual L_h)."""
    return dual_basis(lattice_sum(dual_basis(g), dual_basis(h)))

def canon(v: Vertex):
    """Re-wrap a vertex on its (small, triangular) HNF basis to keep
    coefficient growth bounded during BFS."""
    return Vertex(v.kind, hnf_local(v.g))

# ============================================================================
# PART 0: primitive sanity checks
# ============================================================================
print("=" * 76)
print("PART 0: independent sanity checks of btlib primitives")
print("=" * 76)

# 0.1 v_pi vs norm-based valuation
samples = [CHI, Fw(3), I_SQRT3, OMEGA, ONE, CHI_INV, CHI * CHI,
           Fw(Zw(1, 2), 3)] + [rand_fw() for _ in range(500)]
bad = [x for x in samples if not x.is_zero() and x.v_pi() != v_pi_indep(x)]
report("0.1 v_pi == v_3(field norm) - 2 v_3(den)  (508 samples)",
       not bad and Fw(0).v_pi() == inf,
       f"{len(samples)} sampled; mismatches: {len(bad)}; "
       f"v(chi)=1, v(3)=2, v(i/sqrt3)=-1 confirmed")

# 0.2 numeric embedding vs numpy
ok = np.linalg.norm(mat_c(H_GATE).conj().T @ mat_c(H_GATE) - np.eye(3)) < 1e-12
gates = [H_GATE, S_GATE, R_GATE]
nm = 0
for _ in range(40):
    word = [random.choice(gates) for _ in range(random.randint(2, 8))]
    exact = Mat.identity()
    numer = np.eye(3, dtype=complex)
    for w in word:
        exact = exact * w
        numer = numer @ mat_c(w)
    if np.linalg.norm(mat_c(exact) - numer) > 1e-9:
        nm += 1
norm_ok = all(abs(abs(fw_c(x)) ** 2 - (x.num.norm() / x.den ** 2)) < 1e-9
              for x in samples[:200] if True)
report("0.2 exact arithmetic == numpy complex (H unitary, 40 gate words, norms)",
       ok and nm == 0 and norm_ok,
       f"H unitary numerically; {nm}/40 word mismatches; |x|^2==N(x) on 200")

# 0.3 residue map O_pi -> F_3
bad = 0
ints = [x for x in samples if not x.is_zero() and x.v_pi() >= 0][:200]
for x in ints:
    r = x.res3()
    if r not in (0, 1, 2) or (x - Fw(r)).v_pi() < 1:
        bad += 1
for _ in range(100):
    x, y = random.choice(ints), random.choice(ints)
    if (x + y).res3() != (x.res3() + y.res3()) % 3: bad += 1
    if (x * y).res3() != (x.res3() * y.res3()) % 3: bad += 1
report("0.3 res3 is the ring map O_pi -> F_3 (omega == 1 mod chi)",
       bad == 0 and OMEGA.res3() == 1 and CHI.res3() == 0,
       f"{len(ints)} representative + 100 hom checks; bad={bad}")

# 0.4 hnf_local: column equivalence, triangularity, canonicity
bad_eq = bad_tri = bad_canon = 0
for t in range(40):
    while True:
        g = Mat([[rand_fw(4) for _ in range(3)] for _ in range(3)])
        if not g.det().is_zero():
            break
    h = hnf_local(g)
    q = g.inv() * h
    if not (in_GL_O_indep(q) and in_GL_O_indep(q.inv())):
        bad_eq += 1
    if not all(h.m[i][j].is_zero() for i in range(3) for j in range(i + 1, 3)):
        bad_tri += 1
    if any(h.m[i][i] != chi_pow(h.m[i][i].v_pi()) for i in range(3)):
        bad_tri += 1
    k0 = h.key()
    for _ in range(2):
        if hnf_local(g * rand_gl3_O()).key() != k0:
            bad_canon += 1
    # cross-check against the independent span_basis (same lattice)
    s = span_basis(mat_cols(g))
    if not in_GL_O_indep(s.inv() * h) or not in_GL_O_indep(h.inv() * s):
        bad_eq += 1
report("0.4 hnf_local: lattice-equivalent, lower-tri chi-power diag, canonical",
       bad_eq + bad_tri + bad_canon == 0,
       f"40 random g, 80 GL_3(O) canonicity probes; eq_bad={bad_eq}, "
       f"tri_bad={bad_tri}, canon_bad={bad_canon}")

# 0.5 is_self_dual (Gram test) vs direct  L == L^dual
bad = 0
sd_seen = nsd_seen = 0
trials = []
trials.append(Mat.identity())
trials.append(H_GATE)
for _ in range(60):
    while True:
        g = Mat([[rand_fw(3) for _ in range(3)] for _ in range(3)])
        if not g.det().is_zero():
            break
    trials.append(g)
for m in monomial_matrices()[:30]:
    trials.append(m)
for g in trials:
    a = is_self_dual(g)
    b = lattice_eq(g, dual_basis(g))
    c = is_self_dual_indep(g)
    if not (a == b == c):
        bad += 1
    if a: sd_seen += 1
    else: nsd_seen += 1
report("0.5 self-duality: Gram in GL_3(O) <=> L == L^dual (3 implementations)",
       bad == 0 and sd_seen >= 2 and nsd_seen >= 2,
       f"{len(trials)} lattices ({sd_seen} self-dual, {nsd_seen} not); "
       f"disagreements: {bad}")

# ============================================================================
# PART 1: BFS to radius 7
# ============================================================================
print("=" * 76)
print("PART 1: BFS from origin to radius 7")
print("=" * 76)

RADIUS = 7
origin = canon(ORIGIN)
vert = {origin.key(): origin}           # key -> Vertex
level = {origin.key(): 0}
parent = {origin.key(): None}
nbr_keys = {}                           # key -> list of neighbor keys (r<=6)
levels = [[origin.key()]]
non_tree_edges = []
frontier = [origin]
for r in range(1, RADIUS + 1):
    nxt = []
    for v in frontier:
        ns = [canon(w) for w in neighbors(v)]
        nbr_keys[v.key()] = [w.key() for w in ns]
        for w in ns:
            wk = w.key()
            if wk == parent[v.key()]:
                continue
            if wk in level:
                non_tree_edges.append((v.key(), wk))
                continue
            level[wk] = r
            parent[wk] = v.key()
            vert[wk] = w
            nxt.append(w)
    levels.append([w.key() for w in nxt])
    frontier = nxt
print(f"  built ball: {len(vert)} vertices, sphere sizes "
      f"{[len(l) for l in levels]}")

# (c) sphere sizes
expected = [1, 4, 12, 36, 108, 324, 972, 2916]
got = [len(l) for l in levels]
report("1c  sphere sizes |S_r| = [1,4,12,36,108,324,972,2916] = 4*3^(r-1)",
       got == expected, f"got {got}; total {sum(got)} vertices")

# (d) no non-tree edges
report("1d  zero non-tree edges during BFS (B is a tree out to radius 7)",
       len(non_tree_edges) == 0, f"non-tree edges found: {len(non_tree_edges)}")

# (b) bipartite by level
bad_lvl = []
for r, lv in enumerate(levels):
    want = "P" if r % 2 == 0 else "A"
    if any(vert[k].kind != want for k in lv):
        bad_lvl.append(r)
report("1b  bipartite: even levels pure, odd levels alternating (all 4373)",
       not bad_lvl, f"levels violating: {bad_lvl or 'none'}")

# (a) degree 4 with kinds, all radius <= 6
bad_deg = 0
bad_kind = 0
n_checked = 0
for k, nks in nbr_keys.items():
    n_checked += 1
    if len(set(nks)) != 4 or len(nks) != 4:
        bad_deg += 1
    vkind = vert[k].kind
    for nk in nks:
        w = vert.get(nk)
        if w is not None and w.kind == vkind:
            bad_kind += 1
report("1a-i  every vertex (r<=6) has exactly 4 distinct neighbors, opposite kind",
       bad_deg == 0 and bad_kind == 0,
       f"{n_checked} vertices checked (4 pure-degree, 4 alt-degree); "
       f"bad degree: {bad_deg}, kind violations: {bad_kind}")

# (a) symmetry spot check on ~200 edges: w in N(v) iff v in N(w)
edge_pool = [(parent[k], k) for k in level if parent[k] is not None]
sample_edges = random.sample(edge_pool, 200)
bad_sym = 0
for pk, ck in sample_edges:
    child = vert[ck]
    back = [canon(w).key() for w in neighbors(child)]
    if pk not in back:
        bad_sym += 1
    if ck not in nbr_keys[pk]:
        bad_sym += 1
report("1a-ii  edge symmetry: w in N(v) <=> v in N(w) on 200 sampled edges",
       bad_sym == 0, f"200 edges, both directions; violations: {bad_sym}")

# (e) canonical-key integrity under chi-scaling and right GL_3(O)
bad_key = 0
n_samp = 100
keys_all = list(vert.keys())
for _ in range(n_samp):
    v = vert[random.choice(keys_all)]
    k0 = v.key()
    g = v.g
    if Vertex(v.kind, g * CHI).key() != k0:        bad_key += 1   # pi L
    if Vertex(v.kind, g * chi_pow(3)).key() != k0: bad_key += 1   # pi^3 L
    if Vertex(v.kind, g * rand_gl3_O()).key() != k0: bad_key += 1 # basis change
    u = random.choice([OMEGA, -ONE, OMEGA * OMEGA])
    if Vertex(v.kind, (g * rand_gl3_O()) * (u * CHI)).key() != k0: bad_key += 1
report("1e  canonical key invariant under chi^1/chi^3 scaling + random GL_3(O)",
       bad_key == 0, f"{n_samp} vertices x 4 rescalings; key changes: {bad_key}")

# extra: brute-force re-derivation of neighbors() (independent of residue-form
# machinery): enumerate ALL 13 index-3 super/sub-lattices and test the
# self-dual chain conditions of Prop pr:simplex_classification directly.
def brute_neighbors_of_pure(v):
    """All lattices L1 with L = v.g O^3 c L1, [L1:L]=3, L1# c L1 with
    [L1:L1#]=9 and pi L1 c L1# (chain eq:chain_alternating)."""
    out = []
    cands = []
    for u in F3_LINES:
        U = _complete_to_unimodular(u)
        g1 = v.g * U * Mat.diag(CHI_INV, 1, 1)
        cands.append(g1)
    # the 13 candidates must be pairwise distinct lattices, all containing L
    for i in range(13):
        for j in range(i + 1, 13):
            assert not lattice_eq(cands[i], cands[j]), "candidate collision"
    for g1 in cands:
        assert lattice_contains(g1, v.g) and lattice_index_log3(g1, v.g) == 1
        d = dual_basis(g1)
        if (lattice_contains(g1, d) and lattice_index_log3(g1, d) == 2
                and lattice_contains(d, g1 * CHI)):
            out.append(Vertex("A", g1))
    return out

def brute_neighbors_of_alt(v):
    """All self-dual L0 with L1# c L0 c L1 = v.g O^3, [L1:L0]=3."""
    out = []
    cands = []
    for w1, w2, _ in F3_PLANES:
        U = _complete_to_unimodular(w1, w2)
        cands.append(v.g * U * Mat.diag(1, 1, CHI))
    for i in range(13):
        for j in range(i + 1, 13):
            assert not lattice_eq(cands[i], cands[j]), "candidate collision"
    dv = dual_basis(v.g)
    for g0 in cands:
        assert lattice_contains(v.g, g0) and lattice_index_log3(v.g, g0) == 1
        if is_self_dual_indep(g0) and lattice_contains(g0, dv):
            out.append(Vertex("P", g0))
    return out

bad_bf = 0
n_p = n_a = 0
pool_P = [k for k in keys_all if vert[k].kind == "P" and level[k] <= 5]
pool_A = [k for k in keys_all if vert[k].kind == "A" and level[k] <= 5]
for k in random.sample(pool_P, 20):
    v = vert[k]
    bf = sorted(w.key() for w in brute_neighbors_of_pure(v))
    lib = sorted(w.key() for w in neighbors(v))
    n_p += 1
    if bf != lib or len(bf) != 4:
        bad_bf += 1
for k in random.sample(pool_A, 20):
    v = vert[k]
    bf = sorted(w.key() for w in brute_neighbors_of_alt(v))
    lib = sorted(w.key() for w in neighbors(v))
    n_a += 1
    if bf != lib or len(bf) != 4:
        bad_bf += 1
report("1x  neighbors() == brute-force enumeration of all 13 index-3 "
       "super/sub-lattices + chain conditions",
       bad_bf == 0,
       f"{n_p} pure + {n_a} alternating vertices; each: 13 candidates, "
       f"exactly 4 valid; mismatches: {bad_bf}")

# ============================================================================
# PART 2: uniqueness of the middle vertex for distance-2 pure pairs
# ============================================================================
print("=" * 76)
print("PART 2: middle vertex uniqueness + sum/intersection lattices")
print("=" * 76)

n_pairs = 50
bad_unique = bad_sum = bad_int = bad_chain = bad_dt = 0
pairs_done = 0
pool_P4 = [k for k in keys_all if vert[k].kind == "P" and level[k] <= 4]
while pairs_done < n_pairs:
    vk = random.choice(pool_P4)
    v = vert[vk]
    a = vert[random.choice(nbr_keys[vk])]
    wks = [x for x in nbr_keys[a.key()] if x != vk]
    w = vert[random.choice(wks)]
    pairs_done += 1
    g, h = v.g, w.g
    # exactly one common alternating neighbor, and it is a
    common = set(nbr_keys[vk]) & set(nbr_keys[w.key()])
    if common != {a.key()}:
        bad_unique += 1
    if d_tilde(g, h) != 2:
        bad_dt += 1
    # sum and intersection (independent span_basis)
    S = lattice_sum(g, h)
    I = lattice_intersection(g, h)
    # middle big lattice == L_g + L_h ; middle dual == L_g cap L_h
    if not lattice_eq(a.g, S):
        bad_sum += 1
    if not lattice_eq(dual_basis(a.g), I):
        bad_int += 1
    # consistency: dual(S) == I  (self-duality of L_g, L_h)
    if not lattice_eq(dual_basis(S), I):
        bad_int += 1
    # rigorous "I really is the intersection": I c L_g, I c L_h, [L_g:I]=3,
    # and L_g not c L_h  => any sublattice of L_g containing I and contained
    # in L_h equals I (index 3 is prime).
    okI = (lattice_contains(g, I) and lattice_contains(h, I)
           and lattice_index_log3(g, I) == 1
           and lattice_index_log3(h, I) == 1
           and not lattice_contains(h, g) and not lattice_eq(g, h))
    if not okI:
        bad_int += 1
    # chain eq (1538)/(1544):  pi L_g c I c L_g c S c pi^-1 L_g  (same for h)
    for x in (g, h):
        if not (lattice_contains(I, x * CHI) and lattice_contains(x, I)
                and lattice_contains(S, x)
                and lattice_contains(x * CHI_INV, S)):
            bad_chain += 1
    # index sanity: [S : L_g] = 3, [S : I] = 9
    if lattice_index_log3(S, g) != 1 or lattice_index_log3(S, I) != 2:
        bad_chain += 1
report("2a  exactly ONE common alternating neighbor for distance-2 pure pairs",
       bad_unique == 0, f"{pairs_done} pairs; non-unique/missing: {bad_unique}")
report("2b  middle big lattice == L_g + L_h (independent column-HNF span)",
       bad_sum == 0, f"{pairs_done} pairs; mismatches: {bad_sum}")
report("2c  middle dual lattice == L_g cap L_h (== dual(L_g + L_h))",
       bad_int == 0, f"{pairs_done} pairs x 3 cross-checks; mismatches: {bad_int}")
report("2d  sandwich chain pi L c I c L c S c pi^-1 L for both L_g, L_h; "
       "indices [S:L]=3, [S:I]=9",
       bad_chain == 0, f"{pairs_done} pairs; violations: {bad_chain}")
report("2e  d~(g,h) == 2 for all constructed pairs",
       bad_dt == 0, f"{pairs_done} pairs; violations: {bad_dt}")

# ============================================================================
# PART 3: d~ == BFS graph distance
# ============================================================================
print("=" * 76)
print("PART 3: d~ == graph distance (Prop prop:distance)")
print("=" * 76)

def tree_dist(k1, k2):
    a, b, d = k1, k2, 0
    while level[a] > level[b]:
        a = parent[a]; d += 1
    while level[b] > level[a]:
        b = parent[b]; d += 1
    while a != b:
        a = parent[a]; b = parent[b]; d += 2
    return d

pool_P6 = [k for k in keys_all if vert[k].kind == "P" and level[k] <= 6]
n_dist = 100
bad_d = 0
dist_hist = {}
for _ in range(n_dist):
    k1, k2 = random.choice(pool_P6), random.choice(pool_P6)
    dg = tree_dist(k1, k2)
    dt = d_tilde(vert[k1].g, vert[k2].g)
    dist_hist[dg] = dist_hist.get(dg, 0) + 1
    if dg != dt:
        bad_d += 1
        print(f"    MISMATCH: graph {dg} vs d~ {dt}")
report("3   d~(g,h) == BFS graph distance on 100 random pure pairs (r<=6)",
       bad_d == 0,
       f"100 pairs, graph-distance histogram {dict(sorted(dist_hist.items()))}; "
       f"mismatches: {bad_d}")

# ============================================================================
# PART 4: transitivity of U_3(O_F) on S_0
# ============================================================================
print("=" * 76)
print("PART 4: U_3(O_F) on the radius-2 ball; S_0 transitivity; e_1 = H e_0")
print("=" * 76)

mons = monomial_matrices()
report("4a  #U_3(O_F) = 3! * 6^3 = 1296 monomial matrices, all unitary, in Gamma",
       len(mons) == 1296 and all(m.is_unitary() for m in mons)
       and all(in_Gamma(m) for m in mons[:100]),
       f"{len(mons)} matrices; unitarity all, Gamma membership on 100")

S0_keys = set(levels[2])
A1_keys = set(levels[1])
e1 = canon(Vertex("P", H_GATE))
report("4b  e_1 = H e_0 is a pure vertex with d(e_0, e_1) = 2, e_1 in S_0",
       e1.key() in S0_keys and tree_dist(origin.key(), e1.key()) == 2
       and d_tilde(Mat.identity(), H_GATE) == 2 and ell(H_GATE) == 2,
       f"l(H)={ell(H_GATE)}, d~(I,H)={d_tilde(Mat.identity(), H_GATE)}, "
       f"e_1 in S_0: {e1.key() in S0_keys}")

bad_fix = sum(1 for m in mons if Vertex("P", m).key() != origin.key())
report("4c  all 1296 monomials fix the origin vertex",
       bad_fix == 0, f"1296 checked; moved origin: {bad_fix}")

orbit = {Vertex("P", m * H_GATE).key() for m in mons}
report("4d  U_3(O_F)-orbit of e_1 == S_0 exactly (transitive on 12 vertices)",
       orbit == S0_keys, f"|orbit| = {len(orbit)}, |S_0| = {len(S0_keys)}, "
       f"equal: {orbit == S0_keys}")

bad_perm = 0
S0_verts = [vert[k] for k in levels[2]]
A1_verts = [vert[k] for k in levels[1]]
for m in mons:
    im_S0 = {Vertex("P", m * s.g).key() for s in S0_verts}
    if im_S0 != S0_keys:
        bad_perm += 1
        continue
    im_A1 = {Vertex("A", m * a.g).key() for a in A1_verts}
    if im_A1 != A1_keys:
        bad_perm += 1
report("4e  every monomial permutes S_0 (12) and the 4 alternating r=1 vertices",
       bad_perm == 0, f"1296 monomials x 16 ball vertices; non-permutations: "
       f"{bad_perm}")

# ============================================================================
# PART 5: Gamma-action well-defined and isometric
# ============================================================================
print("=" * 76)
print("PART 5: Gamma = <H,S,R> action: validity, equivariance, isometry")
print("=" * 76)

def random_gamma():
    word = [random.choice("HSR") for _ in range(random.randint(8, 20))]
    g = Mat.identity()
    for c in word:
        g = g * {"H": H_GATE, "S": S_GATE, "R": R_GATE}[c]
    return g, "".join(word)

gammas = [random_gamma() for _ in range(20)]
bad_gamma = sum(1 for g, _ in gammas if not in_Gamma(g))

test_keys = []
for lv in (0, 1, 2, 3, 3, 4, 5, 5, 6, 6):
    test_keys.append(random.choice(levels[lv]))
test_verts = [vert[k] for k in test_keys]

bad_valid = bad_equiv = 0
n_act = 0
for gm, word in gammas:
    for v in test_verts:
        n_act += 1
        gv_basis = gm * v.g
        if v.kind == "P":
            if not is_self_dual_indep(gv_basis):
                bad_valid += 1
                continue
        else:
            d = dual_basis(gv_basis)
            if not (lattice_contains(gv_basis, d)
                    and lattice_index_log3(gv_basis, d) == 2
                    and lattice_contains(d, gv_basis * CHI)):
                bad_valid += 1
                continue
        gv = canon(Vertex(v.kind, gv_basis))
        # equivariance: gamma . N(v) == N(gamma . v)
        lhs = sorted(Vertex(n.kind, gm * n.g).key() for n in neighbors(v))
        rhs = sorted(n.key() for n in neighbors(gv))
        if lhs != rhs:
            bad_equiv += 1
report("5a  20 random words in {H,S,R} are in Gamma = U_3(Z[1/3,omega])",
       bad_gamma == 0, f"word lengths 8-20; failures: {bad_gamma}")
report("5b  gamma.v is a valid vertex of the same kind (self-dual -> self-dual,"
       " alternating chain -> alternating chain)",
       bad_valid == 0, f"{n_act} (gamma, vertex) actions; invalid images: "
       f"{bad_valid}")
report("5c  gamma commutes with neighbors(): gamma.N(v) == N(gamma.v)",
       bad_equiv == 0, f"{n_act} actions x 4 neighbors; mismatches: {bad_equiv}")

bad_iso = 0
n_iso = 0
for gm, word in gammas:
    for _ in range(3):
        k1, k2 = random.choice(pool_P4), random.choice(pool_P4)
        d0 = d_tilde(vert[k1].g, vert[k2].g)
        d1 = d_tilde(gm * vert[k1].g, gm * vert[k2].g)
        n_iso += 1
        if d0 != d1:
            bad_iso += 1
report("5d  d~(gamma x, gamma y) == d~(x, y) (graph distance preserved)",
       bad_iso == 0, f"{n_iso} (gamma, pair) samples; violations: {bad_iso}")

# ============================================================================
print("=" * 76)
if FAILURES:
    print(f"OVERALL: FAIL ({len(FAILURES)} failed): {FAILURES}")
    sys.exit(1)
print("OVERALL: PASS (all tree-structure claims verified)")
sys.exit(0)
