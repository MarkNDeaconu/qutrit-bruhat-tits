#!/usr/bin/env python
"""
v04_residue_forms.py -- Verify the residue-form lemmas of
"Buildings for Synthesis with Clifford+R" (paper.tex ~lines 567-709 and the
appendix proofs ~1245-1363) on actual lattices, with exact arithmetic.

Sections:
  0. Independent sanity checks of the btlib primitives this script relies on
     (numpy-float crosschecks, norm-based v_pi, search-based res3, exact
     inverse/dual identities, HNF canonicality, brute-force coset counting).
  1. Lemma le:nondegen on >= 30 self-dual lattices (origin, H*O^3, random
     unitary words in H,S,R, random pure tree vertices from walks <= 8):
     residue Gram form well-defined / symmetric / non-degenerate over F_3.
  2. ~30 alternating vertices: [L:dual(L)] = 9 (image of dual is a line);
     Lemma le:nondegen2: chi*<,> mod pi is antisymmetric, rank exactly 2,
     radical = image of dual(L); explicit GL_3(F_3) change of basis to the
     Lemma `antisym` shape [[0,a,b],[-a,0,0],[-b,0,0]].
  3. Proof details of le:nondegen2: conj(chi)/chi == -1 mod pi (and the
     intermediate identity conj(chi)/chi = -1 + conj(chi)), and
     conj(u) == u mod pi for (units and all elements) u in O_pi.
  4. Pure-vertex proposition at the origin: enumerate ALL 13 lines of F_3^3;
     exactly the 4 isotropic ones give valid alternating chains
     pi L < pi L1 < dual(L1) < L < L1 (all inclusions strict + exact index
     arithmetic); dual(L1)/pi L == V^perp and pi L1 / pi L == V.
  5. pr:alternatings_has_pure_neighbours on 10 alternating vertices:
     enumerate ALL 13 planes; exactly the 4 self-dual planes (w.r.t. the
     antisymmetric residue form) lift to genuinely self-dual lattices L1
     satisfying the sandwich chain (eq:chain_sandiwich), each plane
     containing the radical line.

Deterministic (seeded RNGs).  Prints PASS/FAIL per claim, exits nonzero on
any failure.
"""
import sys
import os
import itertools
import random
from math import inf

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))

import numpy as np
from btlib import (Zw, Fw, Mat, ZERO, ONE, OMEGA, CHI, CHI_INV, I_SQRT3,
                   chi_pow, H_GATE, S_GATE, R_GATE, gram, dual_basis,
                   is_self_dual, lattice_contains, lattice_eq,
                   lattice_index_log3, hnf_local, residue_matrix,
                   isotropic_lines, selfdual_planes, F3_LINES, F3_PLANES,
                   Vertex, ORIGIN, neighbors, neighbors_of_pure,
                   neighbors_of_alternating, _complete_to_unimodular,
                   _f3_rank)

FAILURES = []


def check(name, cond, detail=""):
    status = "PASS" if cond else "FAIL"
    line = f"[{status}] {name}"
    if detail:
        line += f" -- {detail}"
    print(line)
    if not cond:
        FAILURES.append(name)
    return cond


def section(title):
    print(f"\n=== {title} ===")


# ----------------------------------------------------------------------------
# helpers (independent of btlib internals where it matters)
# ----------------------------------------------------------------------------

W_C = complex(-0.5, np.sqrt(3.0) / 2.0)  # numeric omega


def to_c(x: Fw) -> complex:
    return (x.num.a + x.num.b * W_C) / x.den


def mat_to_np(g: Mat):
    return np.array([[to_c(g[i, j]) for j in range(3)] for i in range(3)])


def v3int(n: int) -> int:
    n = abs(n)
    v = 0
    while n != 0 and n % 3 == 0:
        n //= 3
        v += 1
    return v


def v_pi_indep(x: Fw):
    """Independent v_pi: v_3(N(num)) - 2 v_3(den), N(a+bw)=a^2-ab+b^2.
    Uses only integer arithmetic (chi is the unique prime above 3, N(chi)=3)."""
    if x.num.is_zero():
        return inf
    n = x.num.a * x.num.a - x.num.a * x.num.b + x.num.b * x.num.b
    return v3int(n) - 2 * v3int(x.den)


def res3_indep(x: Fw) -> int:
    """Independent residue: the unique r in {0,1,2} with v_pi(x - r) >= 1."""
    rs = [r for r in range(3) if v_pi_indep(x - Fw(r)) >= 1]
    assert len(rs) == 1, (x, rs)
    return rs[0]


def mat_vec(g: Mat, x):
    return [sum((g[i, k] * x[k] for k in range(3)), ZERO) for i in range(3)]


def herm(u, v):
    """<u, v> = sum_k u_k conj(v_k)   (paper convention)."""
    return sum((u[k] * v[k].conj() for k in range(3)), ZERO)


def herm_pair(g: Mat, i: int, j: int) -> Fw:
    """<g e_i, g e_j> computed directly (independent of btlib.gram)."""
    return sum((g[k, i] * g[k, j].conj() for k in range(3)), ZERO)


def form_matrix(g: Mat, scale: Fw):
    """B[i][j] = res3(scale * <g e_i, g e_j>); requires integrality.
    Then form(x, y) = x^T B y over F_3.  Returns (B, all_integral)."""
    B = [[None] * 3 for _ in range(3)]
    ok = True
    for i in range(3):
        for j in range(3):
            val = scale * herm_pair(g, i, j)
            if val.v_pi() < 0:
                ok = False
                B[i][j] = None
            else:
                B[i][j] = val.res3()
    return B, ok


def f3_bil(B, x, y):
    return sum(x[i] * B[i][j] * y[j] for i in range(3) for j in range(3)) % 3


def f3_radical(B):
    """{x : form(x, y) = 0 for all y} = ker(B^T), as the full set incl. 0."""
    return [x for x in itertools.product(range(3), repeat=3)
            if all(sum(x[i] * B[i][j] for i in range(3)) % 3 == 0
                   for j in range(3))]


def f3_colspace(Mint):
    """Set of all F_3 combinations of the columns of an int 3x3 matrix."""
    cols = [tuple(Mint[i][j] % 3 for i in range(3)) for j in range(3)]
    out = set()
    for a, b, c in itertools.product(range(3), repeat=3):
        out.add(tuple((a * cols[0][i] + b * cols[1][i] + c * cols[2][i]) % 3
                      for i in range(3)))
    return out


def f3_span(*vecs):
    out = set()
    for coeffs in itertools.product(range(3), repeat=len(vecs)):
        out.add(tuple(sum(c * v[i] for c, v in zip(coeffs, vecs)) % 3
                      for i in range(3)))
    return out


def res_int_matrix(g: Mat):
    """Entrywise res3 of an integral Fw matrix (independent of btlib's
    residue_matrix only in call-site; res3 itself is sanity-checked)."""
    return [[g[i, j].res3() for j in range(3)] for i in range(3)]


def strictly_contains(g: Mat, h: Mat) -> bool:
    return lattice_contains(g, h) and not lattice_eq(g, h)


def self_dual_two_ways(g: Mat):
    """(gram criterion, direct dual-lattice equality) -- must agree."""
    r1 = is_self_dual(g)                       # Gram in GL_3(O_pi)
    r2 = lattice_eq(g, dual_basis(g))          # L == L^dual as lattices
    return r1, r2


def random_walk_vertex(length, rng):
    cur, prev = ORIGIN, None
    for _ in range(length):
        nbrs = neighbors(cur)
        if prev is not None:
            filt = [w for w in nbrs if w.key() != prev.key()]
            if filt:
                nbrs = filt
        nxt = rng.choice(nbrs)
        prev, cur = cur, nxt
    return cur


def count_cosets(C: Mat):
    """Number of cosets of C*O^3 in O^3, brute force over {0,1,2}^3 reps.
    Valid when chi*O^3 (hence 3*O^3) is inside C*O^3."""
    Ci = C.inv()
    classes = []
    for x in itertools.product(range(3), repeat=3):
        placed = False
        for cl in classes:
            y = cl[0]
            d = [Fw(x[i] - y[i]) for i in range(3)]
            if all(t.v_pi() >= 0 for t in mat_vec(Ci, d)):
                cl.append(x)
                placed = True
                break
        if not placed:
            classes.append([x])
    return classes


# ============================================================================
# Section 0: sanity-check the primitives
# ============================================================================
section("Section 0: independent sanity checks of btlib primitives")

rng0 = random.Random(40904)


def rand_fw(rng, dens=(1, 1, 1, 2, 3, 5, 9, 27)):
    return Fw(Zw(rng.randrange(-30, 31), rng.randrange(-30, 31)),
              rng.choice(dens))


# 0a: field arithmetic vs numpy complex
ok = True
n_arith = 0
for _ in range(300):
    x, y = rand_fw(rng0), rand_fw(rng0)
    for expr, num in ((x * y, to_c(x) * to_c(y)),
                      (x + y, to_c(x) + to_c(y)),
                      (x - y, to_c(x) - to_c(y)),
                      (x.conj(), to_c(x).conjugate())):
        ok &= abs(to_c(expr) - num) < 1e-9
        n_arith += 1
    if not x.is_zero():
        ok &= (x * x.inv() == ONE)
        n_arith += 1
check("0a field arithmetic (Fw vs numpy complex, exact inverses)", ok,
      f"{n_arith} expressions")

# 0b: v_pi vs independent norm-based valuation
ok = True
samples = [CHI, CHI_INV, I_SQRT3, Fw(3), Fw(1, 3), OMEGA, CHI * CHI,
           Fw(Zw(1, 2), 3) * CHI]
for _ in range(500):
    samples.append(rand_fw(rng0))
for x in samples:
    if not x.is_zero():
        ok &= (x.v_pi() == v_pi_indep(x))
check("0b v_pi vs norm-based independent implementation", ok,
      f"{len(samples)} samples")

# 0c: res3 vs search-based definition
ok = True
n_res = 0
for _ in range(400):
    x = rand_fw(rng0)
    v = x.v_pi()
    if v is inf:
        continue
    if v < 0:
        x = x * chi_pow(-v)   # make integral, often with ugly denominators
    ok &= (x.res3() == res3_indep(x))
    n_res += 1
check("0c res3 vs search-based definition (unique r with v_pi(x-r)>=1)", ok,
      f"{n_res} integral samples")

# 0d: paper convention identities
ok = (Fw(3) == -(OMEGA * OMEGA) * CHI * CHI)
ok &= (Fw(3).v_pi() == 2)
ok &= (I_SQRT3 == (OMEGA - OMEGA * OMEGA) / Fw(3))
ok &= abs(to_c(I_SQRT3) - 1j / np.sqrt(3)) < 1e-12
ok &= (CHI.conj() == Fw(Zw(2, 1)))        # conj(1-w) = 2+w
ok &= (CHI * CHI.conj() == Fw(3))
check("0d ring identities: 3 = -w^2 chi^2, v_pi(3)=2, i/sqrt3=(1+2w)/3", ok)

# 0e: gates unitary (exact + numeric), H matches its definition
Hnp = (1j / np.sqrt(3)) * np.array([[1, 1, 1],
                                    [1, np.exp(2j * np.pi / 3), np.exp(4j * np.pi / 3)],
                                    [1, np.exp(4j * np.pi / 3), np.exp(2j * np.pi / 3)]])
ok = H_GATE.is_unitary() and S_GATE.is_unitary() and R_GATE.is_unitary()
ok &= np.max(np.abs(mat_to_np(H_GATE) - Hnp)) < 1e-12
check("0e H,S,R exactly unitary; H matches (i/sqrt3)[[1,1,1],[1,w,w2],[1,w2,w]]", ok)

# 0f: matrix inverse / dual-basis defining identities (exact), gram vs numpy
ok = True
n_mat = 0
for _ in range(30):
    while True:
        g = Mat([[rand_fw(rng0, dens=(1, 1, 2, 3)) for _ in range(3)]
                 for _ in range(3)])
        if not g.det().is_zero():
            break
    ok &= (g * g.inv() == Mat.identity())
    D = dual_basis(g)
    ok &= (D.conjT() * g == Mat.identity())      # defining pairing of dual
    A = mat_to_np(g)
    ok &= np.max(np.abs(mat_to_np(gram(g)) - A.conj().T @ A)) < 1e-8
    for i in range(3):
        for j in range(3):
            ok &= (gram(g)[j, i] == herm_pair(g, i, j))  # gram == <ge_i,ge_j>
    n_mat += 1
check("0f exact inverse, dual pairing D*g=I, gram vs numpy & direct <,>", ok,
      f"{n_mat} random matrices")

# 0g: hnf_local canonicality:  hnf(g k) == hnf(g) for k in GL_3(O_pi)
ok = True
n_hnf = 0
for _ in range(20):
    while True:
        g = Mat([[rand_fw(rng0, dens=(1, 1, 3)) for _ in range(3)]
                 for _ in range(3)])
        if not g.det().is_zero():
            break
    while True:
        k = Mat([[Fw(Zw(rng0.randrange(-3, 4), rng0.randrange(-3, 4)),
                     rng0.choice([1, 1, 2])) for _ in range(3)]
                 for _ in range(3)])
        if not k.det().is_zero() and k.det().v_pi() == 0 and k.is_O():
            break
    h = hnf_local(g)
    ok &= (hnf_local(g * k).key() == h.key())
    ok &= lattice_eq(h, g)
    n_hnf += 1
check("0g hnf_local: GL_3(O)-invariant key and same lattice", ok,
      f"{n_hnf} random (g,k) pairs")

# 0h: teeth check -- a non-self-dual lattice must fail the tests
g_bad = Mat.diag(CHI, 1, 1)
r1, r2 = self_dual_two_ways(g_bad)
Bbad, intg = form_matrix(g_bad, ONE)
detbad = round(np.linalg.det(np.array(Bbad))) % 3 if intg else None
check("0h negative control: diag(chi,1,1)*O^3 is NOT self-dual and its "
      "residue Gram is degenerate", (not r1) and (not r2) and intg and detbad == 0,
      f"det residue Gram = {detbad} mod 3")

# ============================================================================
# Section 1: Lemma le:nondegen on >= 30 self-dual lattices
# ============================================================================
section("Section 1: Lemma le:nondegen (residue Gram symmetric non-degenerate)")

rng1 = random.Random(112233)
GATES = {"H": H_GATE, "S": S_GATE, "R": R_GATE}

lattices = [("origin O^3", Mat.identity()), ("H*O^3", H_GATE)]
# random unitary words gamma in H,S,R  ->  gamma O^3 self-dual
for t in range(10):
    word = "".join(rng1.choice("HSR") for _ in range(rng1.randrange(4, 13)))
    gmat = Mat.identity()
    for ch in word:
        gmat = gmat * GATES[ch]
    # crosscheck the word product numerically
    Anp = np.eye(3, dtype=complex)
    for ch in word:
        Anp = Anp @ mat_to_np(GATES[ch])
    assert np.max(np.abs(mat_to_np(gmat) - Anp)) < 1e-6
    lattices.append((f"word {word}", gmat))
# random pure tree vertices from walks of even length <= 8
seen_keys = set(Vertex("P", g).key() for _, g in lattices)
attempts = 0
while len(lattices) < 32 and attempts < 120:
    attempts += 1
    ln = rng1.choice([2, 4, 6, 8])
    v = random_walk_vertex(ln, rng1)
    assert v.kind == "P", "even-length walk must end pure (bipartite)"
    if v.key() not in seen_keys:
        seen_keys.add(v.key())
        lattices.append((f"walk{ln} #{attempts}", v.g))

n_distinct = len({Vertex('P', g).key() for _, g in lattices})
check("1.0 collected >= 30 self-dual lattices (>= 25 distinct vertices)",
      len(lattices) >= 30 and n_distinct >= 25,
      f"{len(lattices)} lattices, {n_distinct} distinct vertex classes, "
      f"{attempts} walk attempts")

ok_sd = ok_wd = ok_sym = ok_nondeg = ok_lift = ok_iso4 = ok_agree = True
for label, g in lattices:
    r1, r2 = self_dual_two_ways(g)
    ok_agree &= (r1 == r2)
    ok_sd &= (r1 and r2)
    B, intg = form_matrix(g, ONE)             # res3(<g e_i, g e_j>)
    ok_wd &= intg
    ok_wd &= all(gram(g)[i, j].v_pi() >= 0 for i in range(3) for j in range(3))
    if not intg:
        continue
    # consistency with btlib residue_matrix (transpose relation)
    Mb = residue_matrix(gram(g))
    ok_agree &= all(B[i][j] == Mb[j][i] for i in range(3) for j in range(3))
    ok_sym &= all(B[i][j] == B[j][i] for i in range(3) for j in range(3))
    ok_nondeg &= (round(np.linalg.det(np.array(B, dtype=float))) % 3 != 0)
    ok_nondeg &= (_f3_rank([[B[i][j] for i in range(3)] for j in range(3)]) == 3)
    # well-definedness: form value depends only on lifts mod chi
    for _ in range(2):
        x = [Fw(rng1.randrange(3)) for _ in range(3)]
        y = [Fw(rng1.randrange(3)) for _ in range(3)]
        z = [Fw(Zw(rng1.randrange(-2, 3), rng1.randrange(-2, 3))) for _ in range(3)]
        x2 = [x[i] + CHI * z[i] for i in range(3)]
        v1 = herm(mat_vec(g, x), mat_vec(g, y))
        v2 = herm(mat_vec(g, x2), mat_vec(g, y))
        ok_lift &= ((v1 - v2).v_pi() >= 1)
    # Lemma le:mod3_codes instance: exactly 4 isotropic lines
    ok_iso4 &= (len([u for u in F3_LINES if f3_bil(B, u, u) == 0]) == 4)

n_lat = len(lattices)
check("1.1 all lattices self-dual (Gram criterion AND direct dual equality)",
      ok_sd and ok_agree, f"{n_lat} lattices, two independent routes agree")
check("1.2 le:nondegen (a) residue Gram well-defined: all <g e_i,g e_j> in "
      "O_pi + lift-invariance mod chi", ok_wd and ok_lift,
      f"{n_lat} lattices, 9 entries each + 2 random lift tests each")
check("1.3 le:nondegen (b) residue form symmetric mod pi", ok_sym,
      f"{n_lat} lattices")
check("1.4 le:nondegen (c) residue form non-degenerate over F_3 "
      "(det != 0 and rank 3)", ok_nondeg, f"{n_lat} lattices")
check("1.5 bonus (le:mod3_codes on data): exactly 4 isotropic lines each",
      ok_iso4, f"{n_lat} lattices x 13 lines")

# ============================================================================
# Section 2: alternating vertices, Lemma le:nondegen2
# ============================================================================
section("Section 2: alternating vertices and Lemma le:nondegen2")

rng2 = random.Random(556677)
alts = {}
attempts = 0
while len(alts) < 30 and attempts < 300:
    attempts += 1
    ln = rng2.choice([1, 3, 5, 7])
    v = random_walk_vertex(ln, rng2)
    assert v.kind == "A", "odd-length walk must end alternating (bipartite)"
    if v.key() not in alts:
        alts[v.key()] = v
alt_list = list(alts.values())
check("2.0 collected 30 distinct alternating vertices", len(alt_list) == 30,
      f"{len(alt_list)} vertices from {attempts} walks of length 1..7")

ok_idx = ok_img1 = ok_chain = ok_anti = ok_rank2 = ok_rad = ok_wd2 = True
ok_basis = ok_cosets = ok_negsym = True
n_coset_checked = 0
for idx_v, v in enumerate(alt_list):
    g = v.g                       # big representative L, v_pi(det) = -1
    D = dual_basis(g)
    # (chain_alternating) pi L < dual(L) < L, strictly
    ok_chain &= strictly_contains(g, D)
    ok_chain &= strictly_contains(D, g * CHI)
    # (a) index [L : dual(L)] = 9 exactly via v_pi(det)
    C = g.inv() * D
    ok_idx &= C.is_O()
    ok_idx &= (lattice_index_log3(g, D) == 2)
    ok_idx &= (C.det().v_pi() == 2)
    ok_idx &= (g.det().v_pi() == -1 and D.det().v_pi() == 1)
    # brute-force coset count for the first 3 vertices (independent of det)
    if idx_v < 3:
        classes = count_cosets(C)
        ok_cosets &= (len(classes) == 9 and all(len(c) == 3 for c in classes))
        n_coset_checked += 1
    # image of dual(L) in L/pi L is 1-dimensional
    Cbar = res_int_matrix(C)
    ok_img1 &= (_f3_rank([[Cbar[i][j] for i in range(3)] for j in range(3)]) == 1)
    dual_img = f3_colspace(Cbar)              # the line as a set of 3 vectors
    ok_img1 &= (len(dual_img) == 3)
    # the plain form <,> mod pi is NOT integral (needs the chi twist):
    ok_wd2 &= (min(gram(g)[i, j].v_pi() for i in range(3) for j in range(3)) == -1)
    # (b) chi*<,> mod pi: well-defined, antisymmetric, rank 2, radical = image
    B, intg = form_matrix(g, CHI)
    ok_wd2 &= intg
    if not intg:
        continue
    Mb = residue_matrix(gram(g) * CHI)        # btlib's version: B == Mb^T
    ok_negsym &= all(B[i][j] == Mb[j][i] for i in range(3) for j in range(3))
    ok_anti &= all(B[i][i] == 0 for i in range(3))
    ok_anti &= all((B[i][j] + B[j][i]) % 3 == 0 for i in range(3) for j in range(3))
    ok_rank2 &= (_f3_rank([[B[i][j] for i in range(3)] for j in range(3)]) == 2)
    rad = f3_radical(B)
    ok_rad &= (len(rad) == 3)                 # 1-dimensional radical
    ok_rad &= (set(rad) == dual_img)          # radical == image of dual(L)
    # (c) explicit change of basis to the Lemma `antisym` shape
    r = next(x for x in rad if x != (0, 0, 0))
    p3 = r
    p2 = next(u for u in itertools.product(range(3), repeat=3)
              if u != (0, 0, 0) and len(f3_span(u, r)) == 9)
    p1 = next(u for u in itertools.product(range(3), repeat=3)
              if f3_bil(B, u, p2) != 0)
    P = [p1, p2, p3]   # columns
    Pm = [[P[j][i] for j in range(3)] for i in range(3)]
    ok_basis &= (_f3_rank(P) == 3)            # P in GL_3(F_3)
    N = [[f3_bil(B, P[i], P[j]) for j in range(3)] for i in range(3)]
    a, b = N[0][1], N[0][2]
    shape_ok = (N[0][0] == N[1][1] == N[2][2] == 0
                and N[1][2] == 0 and N[2][1] == 0
                and N[1][0] == (-a) % 3 and N[2][0] == (-b) % 3
                and (a, b) != (0, 0))
    ok_basis &= shape_ok
    # radical of the model form is span((0,b,-a)) (Lemma antisym statement),
    # and P maps it back onto the radical of B = image of dual(L)
    vmodel = (0, b % 3, (-a) % 3)
    ok_basis &= (set(f3_radical(N)) == f3_span(vmodel))
    Pv = tuple(sum(Pm[i][k] * vmodel[k] for k in range(3)) % 3 for i in range(3))
    ok_basis &= (Pv in dual_img and Pv != (0, 0, 0))

n_alt = len(alt_list)
check("2.1 eq:chain_alternating: pi L < dual(L) < L strictly", ok_chain,
      f"{n_alt} vertices")
check("2.2 (a) [L : dual(L)] = 9 exactly via v_pi(det); image of dual is a "
      "line", ok_idx and ok_img1, f"{n_alt} vertices, det-valuation = 2 each")
check("2.3 (a') brute-force coset count [L:dual(L)] = 9 (independent of det "
      "formula)", ok_cosets, f"{n_coset_checked} vertices, 27 reps -> 9 classes of size 3")
check("2.4 le:nondegen2: plain <,> NOT integral (min v_pi = -1) but chi*<,> "
      "well-defined", ok_wd2, f"{n_alt} vertices")
check("2.5 le:nondegen2: chi*<,> mod pi antisymmetric", ok_anti,
      f"{n_alt} vertices")
check("2.6 le:nondegen2: rank exactly 2", ok_rank2, f"{n_alt} vertices")
check("2.7 le:nondegen2: radical == image of dual(L) in L/pi L", ok_rad,
      f"{n_alt} vertices")
check("2.8 (c) explicit P in GL_3(F_3): P^T B P = [[0,a,b],[-a,0,0],[-b,0,0]],"
      " radical (0,b,-a) maps onto image of dual(L)", ok_basis,
      f"{n_alt} explicit changes of basis")
check("2.9 internal consistency: my form matrix == transpose of btlib's "
      "residue_matrix(chi*Gram)", ok_negsym, f"{n_alt} vertices")

# ============================================================================
# Section 3: proof details of le:nondegen2 (paper lines ~1290-1299)
# ============================================================================
section("Section 3: proof details (conj(chi)/chi = -1 mod pi; conj(u) = u mod pi)")

u1 = CHI.conj() / CHI
ok = (u1 == -(OMEGA * OMEGA))                       # conj(chi)/chi = -w^2
ok &= (u1 == Fw(3) / CHI - ONE)                     # = (3-chi)/chi
ok &= (u1 == CHI.conj() - ONE)                      # = -1 + conj(chi)  (paper)
ok &= (u1.v_pi() == 0)                              # a unit of O_pi
ok &= ((u1 + ONE).v_pi() >= 1)                      # == -1 mod pi
ok &= (u1.res3() == 2)
check("3.1 conj(chi)/chi = -w^2 = -1 + conj(chi), a unit, == -1 mod pi", ok,
      f"v_pi(conj(chi)/chi + 1) = {(u1 + ONE).v_pi()}")

# structural: conj(a+bw) - (a+bw) = -b(1+2w) and v_pi(1+2w) = 1
ok = (Fw(Zw(1, 2)).v_pi() == 1)
rng3 = random.Random(778899)
for _ in range(50):
    z = Zw(rng3.randrange(-50, 51), rng3.randrange(-50, 51))
    ok &= (Fw(z.conj()) - Fw(z) == Fw(z.conj() - z))
    ok &= (z.conj() - z == Zw(0, 0) + (-z.b) * Zw(1, 2))
check("3.2 structural identity: conj(x) - x = -b*(1+2w), v_pi(1+2w)=1", ok,
      "50 random Z[w] elements")

ok = True
n_units = n_ints = 0
trials = 0
while (n_units < 100 or n_ints < 200) and trials < 5000:
    trials += 1
    x = rand_fw(rng3)
    v = x.v_pi()
    if v is inf or v < 0:
        continue
    diff_v = (x.conj() - x).v_pi()
    ok &= (diff_v >= 1)
    n_ints += 1
    if v == 0:
        ok &= (x.conj().res3() == x.res3())
        n_units += 1
check("3.3 conj(u) == u mod pi for units (and all integral elements) of O_pi",
      ok and n_units >= 100, f"{n_units} units, {n_ints} integral elements")

# ============================================================================
# Section 4: pure-vertex bijection at the origin (13-line enumeration)
# ============================================================================
section("Section 4: pure-vertex proposition at the origin (4 isotropic lines"
        " <-> 4 alternating chains)")

g0 = Mat.identity()
B0, intg = form_matrix(g0, ONE)
assert intg
iso = [u for u in F3_LINES if f3_bil(B0, u, u) == 0]
check("4.1 residue Gram at origin is the identity; exactly 4 isotropic lines",
      B0 == [[1, 0, 0], [0, 1, 0], [0, 0, 1]] and len(iso) == 4,
      f"isotropic lines: {sorted(iso)}")

passing, vertex_of_line, details_fail = [], {}, []
ok_idx4 = ok_perp = ok_V = True
for u in F3_LINES:
    U = _complete_to_unimodular(u)
    g1 = U * Mat.diag(CHI_INV, 1, 1)          # L1 = O^3 + chi^-1 u O
    D1 = dual_basis(g1)
    conds = {
        "L<L1": strictly_contains(g1, g0),
        "dual(L1)<L": strictly_contains(g0, D1),
        "piL1<dual(L1)": strictly_contains(D1, g1 * CHI),
        "piL<piL1": strictly_contains(g1 * CHI, g0 * CHI),
    }
    good = all(conds.values())
    if good:
        # exact index arithmetic along the chain
        ok_idx4 &= (lattice_index_log3(g1, g0) == 1)          # [L1:L]=3
        ok_idx4 &= (lattice_index_log3(g0, D1) == 1)          # [L:dual(L1)]=3
        ok_idx4 &= (lattice_index_log3(D1, g1 * CHI) == 1)    # [dual(L1):piL1]=3
        ok_idx4 &= (lattice_index_log3(g1 * CHI, g0 * CHI) == 1)  # [piL1:piL]=3
        ok_idx4 &= (lattice_index_log3(g1, D1) == 2)          # [L1:dual(L1)]=9
        ok_idx4 &= (lattice_index_log3(g0, g0 * CHI) == 3)    # [L:piL]=27
        # dual(L1)/pi L == V^perp  inside L/pi L = F_3^3
        Vperp = {x for x in itertools.product(range(3), repeat=3)
                 if f3_bil(B0, u, x) == 0}
        img_dual = f3_colspace(res_int_matrix(D1))
        ok_perp &= (img_dual == Vperp and len(img_dual) == 9)
        # pi L1 / pi L == V
        img_piL1 = f3_colspace(res_int_matrix(g1 * CHI))
        ok_V &= (img_piL1 == f3_span(u))
        ok_V &= (f3_span(u) <= Vperp)         # V subset V^perp
        passing.append(u)
        vertex_of_line[u] = Vertex("A", g1)
    else:
        details_fail.append((u, [k for k, val in conds.items() if not val]))

check("4.2 of all 13 lines, EXACTLY the 4 isotropic ones give valid chains "
      "pi L < pi L1 < dual(L1) < L < L1", sorted(passing) == sorted(iso),
      f"passing={sorted(passing)}; non-isotropic failures: {details_fail}")
check("4.3 exact index arithmetic: [L1:L]=[L:dualL1]=[dualL1:piL1]=[piL1:piL]"
      "=3, [L1:dualL1]=9, [L:piL]=27", ok_idx4, "4 chains")
check("4.4 dual(L1)/pi L == V^perp (as 9-element subspaces of F_3^3)", ok_perp,
      "4 chains")
check("4.5 pi L1/pi L == V and V subset V^perp", ok_V, "4 chains")

keys_constructed = {v.key() for v in vertex_of_line.values()}
keys_btlib = {v.key() for v in neighbors_of_pure(ORIGIN)}
check("4.6 bijection: 4 distinct alternating vertices == btlib's 4 neighbors "
      "of the origin", len(keys_constructed) == 4
      and keys_constructed == keys_btlib and len(keys_btlib) == 4,
      f"{len(keys_constructed)} constructed vs {len(keys_btlib)} from btlib")

# ============================================================================
# Section 5: pr:alternatings_has_pure_neighbours (13-plane enumeration)
# ============================================================================
section("Section 5: alternating vertices have exactly 4 pure neighbours "
        "(sandwich chain eq:chain_sandiwich)")

ten = alt_list[:10]
ok_count = ok_match = ok_sd5 = ok_chain5 = ok_idx5 = ok_radline = ok_dist = True
total_planes_pass = 0
for v in ten:
    g = v.g
    D = dual_basis(g)
    B, _ = form_matrix(g, CHI)
    rad = f3_radical(B)
    good_spans, good_vertices = [], []
    for (w1, w2, span) in F3_PLANES:
        U = _complete_to_unimodular(w1, w2)
        g1 = U_lat = g * U * Mat.diag(1, 1, CHI)   # L1 = g(w1 O + w2 O + chi O^3)
        r1, r2 = self_dual_two_ways(g1)
        sd = r1 and r2
        ok_sd5 &= (r1 == r2)
        if not sd:
            continue
        # sandwich chain: pi L1 < pi L < dual(L) < L1 < L  (all strict)
        c = (strictly_contains(g, g1)
             and strictly_contains(g1, D)
             and strictly_contains(D, g * CHI)
             and strictly_contains(g * CHI, g1 * CHI))
        ok_chain5 &= c
        if c:
            ok_idx5 &= (lattice_index_log3(g, g1) == 1)        # [L:L1]=3
            ok_idx5 &= (lattice_index_log3(g1, D) == 1)        # [L1:dualL]=3
            ok_idx5 &= (lattice_index_log3(D, g * CHI) == 1)   # [dualL:piL]=3
            ok_idx5 &= (lattice_index_log3(g * CHI, g1 * CHI) == 1)
        good_spans.append(frozenset(span))
        good_vertices.append(Vertex("P", g1))
        # the plane contains the radical line (image of dual(L))
        ok_radline &= (set(rad) <= span)
    total_planes_pass += len(good_spans)
    ok_count &= (len(good_spans) == 4 and len(set(good_spans)) == 4)
    # exactly the self-dual planes w.r.t. the antisymmetric residue form
    Mb = residue_matrix(gram(g) * CHI)
    btlib_spans = {frozenset(s) for _, _, s in selfdual_planes(Mb)}
    my_spans = {frozenset(s) for (x1, x2, s) in F3_PLANES
                if all(f3_bil(B, p, q) == 0
                       for p in (x1, x2) for q in (x1, x2))}
    ok_match &= (set(good_spans) == btlib_spans == my_spans)
    # the 4 lifted vertices are distinct and equal btlib's neighbors
    kc = {w.key() for w in good_vertices}
    kb = {w.key() for w in neighbors_of_alternating(v)}
    ok_dist &= (len(kc) == 4 and kc == kb)

check("5.1 exactly 4 of the 13 planes lift to self-dual lattices, per vertex",
      ok_count, f"10 vertices x 13 planes, {total_planes_pass} total passes "
      f"(expected 40)")
check("5.2 the passing planes == self-dual planes of the antisymmetric form "
      "(mine == btlib's)", ok_match, "10 vertices")
check("5.3 each lifted L1 genuinely self-dual (Gram criterion AND direct "
      "dual equality agree)", ok_sd5, "40 lattices, 2 routes each")
check("5.4 sandwich chain pi L1 < pi L < dual(L) < L1 < L, all strict",
      ok_chain5, "40 chains")
check("5.5 exact indices along sandwich: [L:L1]=[L1:dualL]=[dualL:piL]="
      "[piL:piL1]=3", ok_idx5, "40 chains")
check("5.6 each self-dual plane contains the radical line (image of dual(L))",
      ok_radline, "40 planes")
check("5.7 4 distinct pure vertices per alternating vertex, matching btlib "
      "neighbors", ok_dist, "10 vertices")

# ============================================================================
section("Summary")
if FAILURES:
    print(f"FAILED: {len(FAILURES)} check(s):")
    for f in FAILURES:
        print(f"  - {f}")
    sys.exit(1)
print("ALL CHECKS PASSED")
sys.exit(0)
