#!/usr/bin/env python
"""v03_finite_field.py -- exhaustive verification of the F_3 lemmas in
"Buildings for Synthesis with Clifford+R" (paper.tex lines ~491-563,
pure-vertex proposition ~592-611, Appendix C ~1244-1281).

Claims verified (all enumerations are EXHAUSTIVE and from scratch; btlib
is only cross-checked at the end, never used as ground truth):

 1. Lemma le:mod3_codes: for every invertible symmetric A in M_3(F_3)
    (subset of all 3^6 = 729 symmetric matrices), the form
    <x,y>_A = x^T A y has exactly 8 nonzero isotropic vectors and exactly
    4 one-dimensional subspaces V with V subset V^perp, independent of
    whether det A is a quadratic residue.
 2. Lemma antisym (Appendix C): for all 8 pairs (a,b) != (0,0) and
    A = [[0,a,b],[-a,0,0],[-b,0,0]]:
      (i)  v = (0,b,-a) spans the radical (unique up to scaling);
      (ii) exactly 4 two-dimensional V with F_3 v subset V and V^perp = V;
      (iii) the quotient argument: the form descends to F_3^2 as
            c*(x1 y2 - x2 y1) (c = a or b), all 4 lines of F_3^2 are
            self-dual, and their lifts are exactly the 4 planes of (ii).
 3. Pure-vertex proposition ingredients: 13 lines and 13 planes in F_3^3;
    for the identity Gram matrix the 4 isotropic lines are the spans of
    (1,1,1),(1,1,-1),(1,-1,1),(-1,1,1).
 4. Cross-check btlib's F3 helpers (F3_LINES, F3_PLANES, isotropic_lines,
    selfdual_planes) against the from-scratch enumeration, plus an
    independent sympy/numpy sanity check of determinants and isotropic
    vector counts.

Deterministic; exits nonzero on any failure.
"""

import itertools
import random
import sys

sys.path.insert(0, "/Users/markdeaconu/projects/qutrits_v2")

FAILURES = []


def check(name, cond, detail=""):
    status = "PASS" if cond else "FAIL"
    line = f"{status}: {name}"
    if detail:
        line += f" -- {detail}"
    print(line)
    if not cond:
        FAILURES.append(name)


# ---------------------------------------------------------------------------
# From-scratch F_3 linear algebra (independent of btlib)
# ---------------------------------------------------------------------------

VECS = list(itertools.product(range(3), repeat=3))
NONZERO = [v for v in VECS if v != (0, 0, 0)]


def smul(c, v):
    return tuple((c * x) % 3 for x in v)


def vadd(u, v):
    return tuple((a + b) % 3 for a, b in zip(u, v))


def line_of(v):
    """The 1-dim subspace spanned by nonzero v, as a frozenset incl. 0."""
    return frozenset(smul(c, v) for c in range(3))


def span2(v1, v2):
    """Set of all F_3-combinations of v1, v2."""
    return frozenset(vadd(smul(a, v1), smul(b, v2))
                     for a in range(3) for b in range(3))


ALL_LINES = sorted({line_of(v) for v in NONZERO}, key=sorted)
ALL_PLANES = sorted({span2(v1, v2) for v1 in NONZERO for v2 in NONZERO
                     if len(span2(v1, v2)) == 9}, key=sorted)


def bil(A, x, y):
    """<x,y>_A = x^T A y mod 3 (matches the paper's explicit expansion
    in Appendix C, eq. line ~1269)."""
    return sum(x[i] * A[i][j] * y[j] for i in range(3) for j in range(3)) % 3


def bil_Axy(A, x, y):
    """<Ax, y> = x^T A^T y mod 3 (the paper's *stated* convention
    <x,y>_A = <Ax,y>).  Equals bil() for symmetric A; differs by a global
    sign for antisymmetric A."""
    return sum(A[i][j] * x[j] * y[i] for i in range(3) for j in range(3)) % 3


def perp(A, S, form=bil):
    """V^perp = {w : <u,w>_A = 0 for all u in S} as a frozenset."""
    return frozenset(w for w in VECS if all(form(A, u, w) == 0 for u in S))


def det3(A):
    """Determinant mod 3 by explicit cofactor formula."""
    return (A[0][0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1])
            - A[0][1] * (A[1][0] * A[2][2] - A[1][2] * A[2][0])
            + A[0][2] * (A[1][0] * A[2][1] - A[1][1] * A[2][0])) % 3


# ---------------------------------------------------------------------------
# Part 3 first (basic counts used everywhere)
# ---------------------------------------------------------------------------

print("== Part 3: lines/planes of F_3^3 and the identity Gram matrix ==")

check("13 one-dimensional subspaces of F_3^3",
      len(ALL_LINES) == 13, f"found {len(ALL_LINES)} (formula (3^3-1)/(3-1)=13)")
check("13 two-dimensional subspaces of F_3^3",
      len(ALL_PLANES) == 13, f"found {len(ALL_PLANES)}")

ID3 = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
iso_lines_id = {L for L in ALL_LINES
                if L <= perp(ID3, L)}
expected_id = {line_of((1, 1, 1)), line_of((1, 1, 2)),
               line_of((1, 2, 1)), line_of((2, 1, 1))}
check("identity Gram: isotropic lines = spans of (1,1,1),(1,1,-1),(1,-1,1),(-1,1,1)",
      iso_lines_id == expected_id and len(iso_lines_id) == 4,
      f"found {len(iso_lines_id)} isotropic lines; sets equal: "
      f"{iso_lines_id == expected_id}")

# ---------------------------------------------------------------------------
# Part 1: Lemma le:mod3_codes, exhaustively over all 729 symmetric matrices
# ---------------------------------------------------------------------------

print("\n== Part 1: Lemma le:mod3_codes (symmetric invertible forms) ==")

sym_matrices = []
for a11, a22, a33, a12, a13, a23 in itertools.product(range(3), repeat=6):
    sym_matrices.append([[a11, a12, a13], [a12, a22, a23], [a13, a23, a33]])
check("enumerated all 3^6 symmetric matrices", len(sym_matrices) == 729,
      f"count = {len(sym_matrices)}")

invertible = [A for A in sym_matrices if det3(A) != 0]
n_inv = len(invertible)
by_det = {1: 0, 2: 0}
for A in invertible:
    by_det[det3(A)] += 1
print(f"INFO: invertible symmetric matrices: {n_inv} of 729 "
      f"(det=1: {by_det[1]}, det=2: {by_det[2]})")
check("invertible symmetric count matches orbit computation 2*|GL_3(F_3)|/48 = 468",
      n_inv == 468 and by_det[1] == 234 and by_det[2] == 234,
      f"total {n_inv}, det=1: {by_det[1]}, det=2: {by_det[2]}")

bad_veccount = []
bad_linecount = []
bad_convention = []
iso_vec_counts = set()
iso_line_counts = set()
counts_by_disc = {1: set(), 2: set()}
for A in invertible:
    nvec = sum(1 for x in NONZERO if bil(A, x, x) == 0)
    iso_vec_counts.add(nvec)
    counts_by_disc[det3(A)].add(nvec)
    if nvec != 8:
        bad_veccount.append(A)
    # full subspace condition V subset V^perp (all 9 pairs annihilate)
    good = [L for L in ALL_LINES if L <= perp(A, L)]
    iso_line_counts.add(len(good))
    if len(good) != 4:
        bad_linecount.append(A)
    # both stated conventions <Ax,y> and x^T A y must agree (A symmetric)
    if any(bil(A, x, y) != bil_Axy(A, x, y)
           for x in NONZERO[:4] for y in NONZERO):
        bad_convention.append(A)

check("all 468 invertible symmetric A: exactly 8 nonzero isotropic vectors",
      not bad_veccount,
      f"counts seen: {sorted(iso_vec_counts)}; failures: {len(bad_veccount)}")
check("8 isotropic vectors independent of discriminant class",
      counts_by_disc[1] == {8} and counts_by_disc[2] == {8},
      f"det=1 counts: {sorted(counts_by_disc[1])}, "
      f"det=2 counts: {sorted(counts_by_disc[2])}")
check("all 468 invertible symmetric A: exactly 4 lines V with V subset V^perp",
      not bad_linecount,
      f"counts seen: {sorted(iso_line_counts)}; failures: {len(bad_linecount)}")
check("conventions <Ax,y> and x^T A y agree for symmetric A",
      not bad_convention, f"failures: {len(bad_convention)}")

# ---------------------------------------------------------------------------
# Part 2: Lemma antisym (Appendix C), all 8 pairs (a,b) != (0,0)
# ---------------------------------------------------------------------------

print("\n== Part 2: Lemma antisym (Appendix C) ==")

pairs = [(a, b) for a in range(3) for b in range(3) if (a, b) != (0, 0)]
check("8 pairs (a,b) in F_3^2 \\ {(0,0)}", len(pairs) == 8, f"count = {len(pairs)}")

LINES_F32 = sorted({frozenset(((c * u) % 3, (c * w) % 3) for c in range(3))
                    for u in range(3) for w in range(3) if (u, w) != (0, 0)},
                   key=sorted)
check("F_3^2 has exactly 4 lines", len(LINES_F32) == 4,
      f"count = {len(LINES_F32)} (formula (9-1)/2 = 4)")

all_rad_ok, all_plane_ok, all_quot_ok, all_lift_ok = True, True, True, True
all_sign_note, literal_id_ok_pairs, literal_id_bad_pairs = True, [], []
plane_counts = set()

for (a, b) in pairs:
    A = [[0, a, b], [(-a) % 3, 0, 0], [(-b) % 3, 0, 0]]
    v = (0, b % 3, (-a) % 3)

    # -- sanity: paper's explicit expansion (line ~1269) vs conventions
    # paper claims <x,y>_A = a(x1y2-x2y1) + b(x1y3-x3y1); check which
    # convention realizes it.
    def paper_form(x, y):
        return (a * (x[0] * y[1] - x[1] * y[0])
                + b * (x[0] * y[2] - x[2] * y[0])) % 3
    matches_xAy = all(bil(A, x, y) == paper_form(x, y)
                      for x in VECS for y in VECS)
    matches_Axy = all(bil_Axy(A, x, y) == paper_form(x, y)
                      for x in VECS for y in VECS)
    if not (matches_xAy and not matches_Axy):
        all_sign_note = False
    # both conventions give the same radical / self-dual sets (global sign):
    sign_irrelevant = all(
        perp(A, S, bil) == perp(A, S, bil_Axy) for S in ALL_PLANES)

    # (i) radical = F_3 v, unique line
    radical = frozenset(x for x in VECS
                        if all(bil(A, x, y) == 0 for y in VECS))
    rad_ok = (radical == line_of(v) and len(radical) == 3)
    all_rad_ok &= rad_ok and sign_irrelevant

    # (ii) self-dual planes containing v
    sd_planes = [P for P in ALL_PLANES if perp(A, P) == P]
    sd_with_v = [P for P in sd_planes if v in P]
    plane_counts.add(len(sd_with_v))
    # v in P is automatic when P^perp = P (radical subset P^perp):
    plane_ok = (len(sd_with_v) == 4 and len(sd_planes) == 4)
    all_plane_ok &= plane_ok

    # (iii) quotient argument.
    #  Well-definedness on cosets of F_3 v:
    descends = all(bil(A, vadd(x, smul(t, v)), vadd(y, smul(s, v)))
                   == bil(A, x, y)
                   for x in VECS for y in VECS
                   for t in range(3) for s in range(3))
    #  Identification psi: F_3^3 / F_3 v -> F_3^2.  The paper's literal map
    #  x -> (x1, x2) (for a != 0) only kills v when b == 0; record this, and
    #  use the corrected map x -> (x1, x2 + (b/a) x3) (inverse of a is a
    #  itself in F_3), which reduces to the paper's when b == 0.
    if a != 0:
        literal_kills_v = (v[0] % 3 == 0 and v[1] % 3 == 0)  # psi(v)=(0,b)
        (literal_id_ok_pairs if literal_kills_v
         else literal_id_bad_pairs).append((a, b))
        inv_a = a  # 1*1=1, 2*2=4=1 mod 3
        def psi(x):
            return (x[0] % 3, (x[1] + b * inv_a * x[2]) % 3)
        c = a
    else:  # a == 0, b != 0: identify via (x1, x3), form b(x1y3 - x3y1)
        def psi(x):
            return (x[0] % 3, x[2] % 3)
        c = b
    psi_kills_v = (psi(v) == (0, 0))
    # psi is surjective with kernel exactly F_3 v -> bijection on quotient
    psi_fibers_ok = (len({psi(x) for x in VECS}) == 9 and
                     all((psi(x) == psi(y)) ==
                         (tuple((xi - yi) % 3 for xi, yi in zip(x, y))
                          in line_of(v))
                         for x in VECS for y in VECS))
    # descended form is exactly c*(u1 w2 - u2 w1) on F_3^2:
    form_matches = all(
        bil(A, x, y) == (c * (psi(x)[0] * psi(y)[1]
                              - psi(x)[1] * psi(y)[0])) % 3
        for x in VECS for y in VECS)
    # all 4 lines of F_3^2 self-dual under c*standard-symplectic:
    def symp(u, w):
        return (c * (u[0] * w[1] - u[1] * w[0])) % 3
    sd_lines_f32 = [L for L in LINES_F32
                    if frozenset(w for w in itertools.product(range(3),
                                                              repeat=2)
                                 if all(symp(u, w) == 0 for u in L)) == L]
    quot_ok = (descends and psi_kills_v and psi_fibers_ok and form_matches
               and len(sd_lines_f32) == 4)
    all_quot_ok &= quot_ok

    # lifts F_3 v + F_3 w of the 4 self-dual lines equal the planes of (ii)
    lifted = set()
    for L in sd_lines_f32:
        w2 = next(u for u in L if u != (0, 0))
        wlift = next(x for x in VECS if psi(x) == w2)
        lifted.add(span2(v, wlift))
    all_lift_ok &= (lifted == set(sd_with_v))

check("(i) radical is exactly F_3*(0,b,-a) for all 8 pairs (both sign conventions)",
      all_rad_ok)
check("(ii) exactly 4 self-dual planes V (= V^perp) containing v, all 8 pairs",
      all_plane_ok, f"counts seen: {sorted(plane_counts)}; "
      "self-dual planes automatically contain v")
check("(iii) form descends to F_3^2 as c*(x1y2-x2y1), all 4 lines of F_3^2 "
      "self-dual, all 8 pairs", all_quot_ok)
check("(iii) lifts F_3 v + F_3 w of the 4 quotient lines = the 4 planes of (ii)",
      all_lift_ok)
check("paper's explicit expansion (App. C) equals x^T A y = <x,Ay>, i.e. is "
      "the NEGATIVE of the stated convention <Ax,y> (inconsequential sign typo)",
      all_sign_note,
      "global sign does not change radical/perp/self-dual sets (verified)")
check("paper's literal identification x -> (x1,x2) kills v exactly when b == 0",
      sorted(literal_id_ok_pairs) == [(1, 0), (2, 0)]
      and sorted(literal_id_bad_pairs) == [(1, 1), (1, 2), (2, 1), (2, 2)],
      f"well-defined for (a,b) in {sorted(literal_id_ok_pairs)}; NOT "
      f"well-defined for {sorted(literal_id_bad_pairs)} (corrected map "
      "x -> (x1, x2 + (b/a)x3) used; lemma conclusion unaffected)")

# ---------------------------------------------------------------------------
# Part 4: cross-check btlib helpers + independent sympy/numpy sanity checks
# ---------------------------------------------------------------------------

print("\n== Part 4: btlib cross-checks and independent sanity checks ==")

import btlib  # noqa: E402

check("btlib.F3_LINES has 13 entries spanning the 13 distinct lines",
      len(btlib.F3_LINES) == 13 and
      {line_of(v) for v in btlib.F3_LINES} == set(ALL_LINES),
      f"len = {len(btlib.F3_LINES)}")
check("btlib.F3_PLANES has 13 entries whose spans are the 13 distinct planes",
      len(btlib.F3_PLANES) == 13 and
      {frozenset(p[2]) for p in btlib.F3_PLANES} == set(ALL_PLANES),
      f"len = {len(btlib.F3_PLANES)}")

mismatch_iso = 0
for A in invertible:
    mine = {L for L in ALL_LINES if L <= perp(A, L)}
    theirs = {line_of(v) for v in btlib.isotropic_lines(A)}
    if mine != theirs:
        mismatch_iso += 1
check("btlib.isotropic_lines agrees with from-scratch enumeration on all 468 "
      "invertible symmetric matrices", mismatch_iso == 0,
      f"mismatches: {mismatch_iso}")

mismatch_sd = 0
for (a, b) in pairs:
    A = [[0, a, b], [(-a) % 3, 0, 0], [(-b) % 3, 0, 0]]
    mine = {P for P in ALL_PLANES if perp(A, P) == P}
    theirs = {frozenset(p[2]) for p in btlib.selfdual_planes(A)}
    if mine != theirs:
        mismatch_sd += 1
check("btlib.selfdual_planes agrees with from-scratch enumeration on all 8 "
      "antisymmetric matrices", mismatch_sd == 0, f"mismatches: {mismatch_sd}")

# independent determinant check (sympy) over all 729 symmetric matrices
import sympy  # noqa: E402

sym_det_ok = all((int(sympy.Matrix(A).det()) % 3) == det3(A)
                 for A in sym_matrices)
check("from-scratch det mod 3 matches sympy on all 729 symmetric matrices",
      sym_det_ok)

# independent numpy isotropic-vector count on a seeded random sample
import numpy as np  # noqa: E402

rng = random.Random(12345)
sample = rng.sample(invertible, 25)
X = np.array(NONZERO)  # 26 x 3
np_ok = True
for A in sample:
    M = np.array(A)
    q = np.einsum('ni,ij,nj->n', X, M, X) % 3
    if int((q == 0).sum()) != 8:
        np_ok = False
check("numpy quadratic-form evaluation: 8 nonzero isotropic vectors on a "
      "seeded sample of 25 invertible symmetric matrices (seed 12345)", np_ok)

# ---------------------------------------------------------------------------

print()
if FAILURES:
    print(f"OVERALL: FAIL ({len(FAILURES)} failed checks): {FAILURES}")
    sys.exit(1)
print("OVERALL: PASS (all checks passed)")
