#!/usr/bin/env python
"""
v08_synthesis.py -- End-to-end verification of the exact synthesis algorithm
for the qutrit Clifford+R gate set (paper: "Buildings for Synthesis with
Clifford+R"; algorithmic content of Theorem th:ring_equality_p_is_3 and the
proof in Section "Proof of main theorem" / ss:proof).

Checks (one PASS/FAIL line per claim, nonzero exit on any FAIL):

  0. Independent sanity of btlib primitives used here (btlib is a tool, not
     ground truth): Fw/Mat arithmetic vs numpy complex floats, v_pi vs a
     from-scratch norm-based valuation, monomials = U_3(O_F) count 1296,
     local-HNF canonical-form invariance under GL_3(O_pi).
  1. Round trip: 100 seeded random words over {H,S,R} of lengths 1..40;
     synthesis_steps(U) reconstructs U EXACTLY; word length == sde(U)
     (= l(U)/2, computed from scratch); trailing M is monomial.
  2. Descent invariant: at every synthesis step l drops by exactly 2 and
     EXACTLY 1 of the 12 coset reps descends; full distribution of
     candidate l-values is {l-2:1, l:2, l+2:9} at every step.
  3. Geodesic/minimality: word length == d(e_0, U e_0)/2 with d the BFS
     graph distance in the tree (cross-checked on examples with sde 0..3).
  4. Normal-form counting: #S_0 = 12, #U_3(O_F) = 1296, stabilizer of e_1
     inside U_3(O_F) has exactly 108 elements (paper: #G = 108);
     #{U : sde(U)=1} = 12*1296 = 15552 (full enumeration, all distinct);
     #sphere(4) = 108 = 4*3^3 via BFS; #{U : sde(U)=2} = 108*1296 = 139968
     (full enumeration m1 H m2 H m3 over 12 x 9 x 1296, all distinct);
     3000 random sde-2 unitaries each occur in the enumeration (uniqueness
     of the normal form).
  5. Edge cases: I, omega*I, -I, S, R, H, H^2, SHSHS, HSH, RHSHR.
  6. Performance probe: time synthesis for sde = 5, 10, 20, 40.

Deterministic: all RNGs seeded.
"""

import itertools
import os
import random
import sys
import time
from collections import Counter
from math import inf

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))

import numpy as np

from btlib import (
    Zw, Fw, Mat, ZERO, ONE, OMEGA, CHI, chi_pow,
    H_GATE, S_GATE, R_GATE, I_SQRT3, _UNITS,
    monomial_matrices, in_Gamma, ell, sde, d_tilde,
    Vertex, ORIGIN, bfs_tree, lattice_eq, hnf_local,
    s0_coset_reps, synthesis_steps, reconstruct,
)

SEED = 20260609
FAILURES = []


def check(name, ok, detail=""):
    tag = "PASS" if ok else "FAIL"
    print(f"[{tag}] {name}" + (f"  -- {detail}" if detail else ""))
    if not ok:
        FAILURES.append(name)
    return ok


# ---------------------------------------------------------------------------
# From-scratch reimplementations (independent of btlib internals)
# ---------------------------------------------------------------------------

OMEGA_C = complex(-0.5, np.sqrt(3.0) / 2.0)   # exp(2*pi*i/3)


def fw_to_c(x: Fw) -> complex:
    return (x.num.a + x.num.b * OMEGA_C) / x.den


def mat_to_np(M: Mat) -> np.ndarray:
    return np.array([[fw_to_c(x) for x in row] for row in M.m], dtype=complex)


def v3_int(n: int) -> int:
    """3-adic valuation of a nonzero integer (from scratch)."""
    n = abs(n)
    assert n != 0
    v = 0
    while n % 3 == 0:
        n //= 3
        v += 1
    return v


def vpi_scratch(x: Fw):
    """v_pi from scratch via the field norm:
    for z = a + b*omega, N(z) = a^2 - a*b + b^2 = z * conj(z), so
    v_pi(z) = v_3(N(z)) (since v_pi(3) = 2 and v_pi(conj z) = v_pi(z));
    v_pi(num/den) = v_3(N(num)) - 2*v_3(den)."""
    if x.num.is_zero():
        return inf
    nrm = x.num.a * x.num.a - x.num.a * x.num.b + x.num.b * x.num.b
    return v3_int(nrm) - 2 * (v3_int(x.den) if x.den % 3 == 0 else 0)


def ell_scratch(M: Mat) -> int:
    return -2 * min(vpi_scratch(x) for row in M.m for x in row)


def sde_scratch(M: Mat) -> int:
    e = ell_scratch(M)
    assert e % 2 == 0
    return e // 2


UNIT_KEYS = {(1, 0), (0, 1), (-1, -1), (-1, 0), (0, -1), (1, 1)}
#            1       w       w^2=-1-w   -1       -w       -w^2=1+w


def is_monomial_scratch(M: Mat) -> bool:
    """One nonzero entry per row and column, each a unit of Z[omega]."""
    rowc, colc = [0, 0, 0], [0, 0, 0]
    for i in range(3):
        for j in range(3):
            x = M.m[i][j]
            if x.is_zero():
                continue
            if x.den != 1 or (x.num.a, x.num.b) not in UNIT_KEYS:
                return False
            rowc[i] += 1
            colc[j] += 1
    return rowc == [1, 1, 1] and colc == [1, 1, 1]


def mono_from(perm, uidx):
    """Monomial matrix with rows[i][perm[i]] = UNITS[uidx[i]] (btlib order)."""
    rows = [[ZERO] * 3 for _ in range(3)]
    for i in range(3):
        rows[i][perm[i]] = _UNITS[uidx[i]]
    return Mat(rows)


# ===========================================================================
# Section 0: independent sanity checks of btlib primitives
# ===========================================================================
print("=" * 76)
print("Section 0: sanity-check btlib primitives independently")
print("=" * 76)

rng = random.Random(SEED)

# 0a. Gate matrices match the paper, numerically and exactly unitary.
H_np = (1j / np.sqrt(3.0)) * np.array(
    [[1, 1, 1], [1, OMEGA_C, OMEGA_C**2], [1, OMEGA_C**2, OMEGA_C]], dtype=complex)
S_np = np.diag([1, OMEGA_C, 1]).astype(complex)
R_np = np.diag([1, 1, -1]).astype(complex)
ok = (np.max(np.abs(mat_to_np(H_GATE) - H_np)) < 1e-12
      and np.max(np.abs(mat_to_np(S_GATE) - S_np)) < 1e-12
      and np.max(np.abs(mat_to_np(R_GATE) - R_np)) < 1e-12)
check("0a gates H,S,R match paper definition (numpy floats)", ok,
      f"i/sqrt3 = (1+2w)/3 = {fw_to_c(I_SQRT3):.6f}")
check("0b gates exactly unitary (exact arithmetic)",
      H_GATE.is_unitary() and S_GATE.is_unitary() and R_GATE.is_unitary())

# 0c. v_chi on Z[omega] vs from-scratch norm-based valuation (2000 samples).
n_ok, n_tot = 0, 0
chiZ = Zw(1, -1)
for _ in range(2000):
    a, b = rng.randint(-200, 200), rng.randint(-200, 200)
    if a == 0 and b == 0:
        continue
    x = Zw(a, b)
    for _ in range(rng.randint(0, 5)):   # mix in high valuations
        x = x * chiZ
    n_tot += 1
    if x.v_chi() == v3_int(x.a * x.a - x.a * x.b + x.b * x.b):
        n_ok += 1
check("0c Zw.v_chi == v_3(norm) on random Z[omega] elements",
      n_ok == n_tot, f"{n_ok}/{n_tot} samples agree")

# 0d. Fw.v_pi vs scratch (with denominators).
n_ok, n_tot = 0, 0
for _ in range(2000):
    a, b = rng.randint(-300, 300), rng.randint(-300, 300)
    if a == 0 and b == 0:
        continue
    d = rng.choice([1, 2, 3, 9, 27, 81, 6, 12, 243])
    x = Fw(Zw(a, b), d)
    if x.is_zero():
        continue
    n_tot += 1
    if x.v_pi() == vpi_scratch(x):
        n_ok += 1
check("0d Fw.v_pi == norm-based v_pi (independent reimplementation)",
      n_ok == n_tot, f"{n_ok}/{n_tot} samples agree")

# 0e. Mat multiply / inverse / det vs numpy.
n_ok, n_tot = 0, 0
for _ in range(50):
    A = Mat([[Fw(Zw(rng.randint(-9, 9), rng.randint(-9, 9)),
                 rng.choice([1, 1, 3])) for _ in range(3)] for _ in range(3)])
    B = Mat([[Fw(Zw(rng.randint(-9, 9), rng.randint(-9, 9)),
                 rng.choice([1, 1, 3])) for _ in range(3)] for _ in range(3)])
    if A.det().is_zero():
        continue
    n_tot += 1
    An, Bn = mat_to_np(A), mat_to_np(B)
    ok = (np.max(np.abs(mat_to_np(A * B) - An @ Bn)) < 1e-9
          and np.max(np.abs(mat_to_np(A.inv()) - np.linalg.inv(An))) < 1e-6
          and abs(fw_to_c(A.det()) - np.linalg.det(An)) < 1e-6)
    n_ok += ok
check("0e Mat mul/inv/det agree with numpy on random matrices",
      n_ok == n_tot, f"{n_ok}/{n_tot}")

# 0f. monomial_matrices: 1296 = 3! * 6^3 distinct, all unitary, all in Gamma.
MONOS = monomial_matrices()
keys = {m.key() for m in MONOS}
ok = (len(MONOS) == 1296 and len(keys) == 1296
      and all(m.is_unitary() for m in MONOS)
      and all(is_monomial_scratch(m) for m in MONOS)
      and all(in_Gamma(m) for m in MONOS))
check("0f monomial_matrices(): exactly 1296 = 3!*6^3 distinct unitary monomials",
      ok, f"count={len(MONOS)}, distinct={len(keys)}")

# 0g. hnf_local canonical: hnf(g*U) == hnf(g) for random U in GL_3(O_pi),
#     and hnf(g) spans the same lattice as g.
def random_gl_O(rng):
    """Random element of GL_3(O_pi): product of monomials and integer
    unipotents (all have unit determinant at pi)."""
    g = MONOS[rng.randrange(1296)]
    for _ in range(3):
        i, j = rng.sample(range(3), 2)
        rows = [[ONE if r == c else ZERO for c in range(3)] for r in range(3)]
        rows[i][j] = Fw(Zw(rng.randint(-3, 3), rng.randint(-3, 3)))
        g = g * Mat(rows) * MONOS[rng.randrange(1296)]
    return g

n_ok, n_tot = 0, 0
for _ in range(40):
    g = Mat([[Fw(Zw(rng.randint(-6, 6), rng.randint(-6, 6)),
                 rng.choice([1, 3])) for _ in range(3)] for _ in range(3)])
    if g.det().is_zero():
        continue
    U_gl = random_gl_O(rng)
    n_tot += 1
    n_ok += (hnf_local(g * U_gl).key() == hnf_local(g).key()
             and lattice_eq(g, hnf_local(g)))
check("0g hnf_local canonical on GL_3(O_pi)-cosets (vertex keys well-defined)",
      n_ok == n_tot, f"{n_ok}/{n_tot}")

# 0h. d_tilde(g,h) == l(g^-1 h) for unitary g,h  (Lemma le:ispositive(3)).
n_ok = 0
for _ in range(10):
    w1 = ''.join(rng.choice('HSR') for _ in range(rng.randint(1, 8)))
    w2 = ''.join(rng.choice('HSR') for _ in range(rng.randint(1, 8)))
    G = {'H': H_GATE, 'S': S_GATE, 'R': R_GATE}
    g = Mat.identity()
    for c in w1:
        g = g * G[c]
    h = Mat.identity()
    for c in w2:
        h = h * G[c]
    n_ok += (d_tilde(g, h) == ell(g.inv() * h) == ell(h.inv() * g))
check("0h d~(g,h) == l(g^-1 h) == l(h^-1 g) for unitaries (l(x)=l(x^-1))",
      n_ok == 10, f"{n_ok}/10")


def random_unitary_with_sde(target, rng):
    """Random reduced word: multiply by random (m*H) keeping ascending steps."""
    cur = Mat.identity()
    while ell(cur) < 2 * target:
        cand = cur * (MONOS[rng.randrange(1296)] * H_GATE)
        if ell(cand) == ell(cur) + 2:
            cur = cand
    return cur


# ===========================================================================
# Section 1: round trip on 100 random words (lengths 1..40)
# ===========================================================================
print()
print("=" * 76)
print("Section 1: round-trip synthesis on 100 random {H,S,R} words")
print("=" * 76)

rng = random.Random(SEED + 1)
GATES = {'H': H_GATE, 'S': S_GATE, 'R': R_GATE}

word_strs, unitaries = [], []
for _ in range(100):
    L = rng.randint(1, 40)
    w = ''.join(rng.choice('HSR') for _ in range(L))
    word_strs.append(w)
    U = Mat.identity()
    for ch in w:
        U = U * GATES[ch]
    unitaries.append(U)

t0 = time.perf_counter()
synth_results = []
n_recon = n_len = n_mono = n_num = n_gamma = 0
sdes = []
for w, U in zip(word_strs, unitaries):
    n_gamma += in_Gamma(U)
    word, M = synthesis_steps(U)
    synth_results.append((word, M))
    R = reconstruct(word, M)
    n_recon += (R == U)                                   # EXACT equality
    s = sde_scratch(U)
    sdes.append(s)
    n_len += (len(word) == s and ell(U) == ell_scratch(U))
    n_mono += is_monomial_scratch(M)
    n_num += (np.max(np.abs(mat_to_np(R) - mat_to_np(U))) < 1e-9)
t1 = time.perf_counter()

check("1a all 100 inputs lie in U_3(Z[omega,1/3]) (unitary, 3-power denoms)",
      n_gamma == 100, f"{n_gamma}/100")
check("1b reconstruct(word, M) == U EXACTLY for all 100 words",
      n_recon == 100, f"{n_recon}/100 exact matrix equalities over Q(omega)")
check("1c word length == sde(U) = l(U)/2 (sde from scratch) for all 100",
      n_len == 100,
      f"{n_len}/100; sde range [{min(sdes)},{max(sdes)}], "
      f"mean {sum(sdes)/len(sdes):.2f}")
check("1d trailing M is a monomial (from-scratch test) for all 100",
      n_mono == 100, f"{n_mono}/100")
check("1e numpy float cross-check |reconstruct - U| < 1e-9 for all 100",
      n_num == 100, f"{n_num}/100  (total section time {t1-t0:.1f}s)")

# Random {H,S,R} words have low sde (synthesis finds much shorter circuits),
# so additionally round-trip 30 seeded random REDUCED words with sde 1..30.
rng_hi = random.Random(SEED + 6)
hi_pairs = []
n_ok, hi_sdes = 0, []
for _ in range(30):
    target = rng_hi.randint(1, 30)
    U = random_unitary_with_sde(target, rng_hi)
    word, M = synthesis_steps(U)
    hi_pairs.append((U, (word, M)))
    s = sde_scratch(U)
    hi_sdes.append(s)
    n_ok += (s == target and reconstruct(word, M) == U and len(word) == s
             and is_monomial_scratch(M) and in_Gamma(U)
             and np.max(np.abs(mat_to_np(reconstruct(word, M))
                               - mat_to_np(U))) < 1e-8)
check("1f exact round trip on 30 random reduced words with sde 1..30",
      n_ok == 30,
      f"{n_ok}/30; sde range [{min(hi_sdes)},{max(hi_sdes)}], "
      f"mean {sum(hi_sdes)/len(hi_sdes):.1f}")


# ===========================================================================
# Section 2: descent invariant -- exactly 1 of 12 reps descends each step
# ===========================================================================
print()
print("=" * 76)
print("Section 2: descent invariant (instrumented synthesis)")
print("=" * 76)

REPS = s0_coset_reps()
check("2a s0_coset_reps(): exactly 12 monomial coset reps (#S_0 = 12)",
      len(REPS) == 12 and all(is_monomial_scratch(m) for m in REPS),
      f"count={len(REPS)}")

MH = [m * H_GATE for m in REPS]
MHinv = [x.inv() for x in MH]

desc_counter = Counter()        # number of descending reps per step
pattern_counter = Counter()     # distribution of l(cand) - l(cur) per step
steps_total = 0
all_match = True
all_drop2 = True

descent_inputs = list(zip(unitaries, synth_results)) + hi_pairs

for U, (lib_word, lib_M) in descent_inputs:
    cur = U
    my_word = []
    while ell(cur) > 0:
        l = ell(cur)
        cands = [MHinv[i] * cur for i in range(12)]
        ells = [ell(c) for c in cands]
        rel = tuple(sorted(Counter(e - l for e in ells).items()))
        pattern_counter[rel] += 1
        desc_idx = [i for i, e in enumerate(ells) if e == l - 2]
        desc_counter[len(desc_idx)] += 1
        steps_total += 1
        if len(desc_idx) == 0:
            all_drop2 = False
            break
        i = desc_idx[0]
        my_word.append(REPS[i])
        if ells[i] != l - 2:
            all_drop2 = False
        cur = cands[i]
    # compare with library run (descent is unique => words must coincide)
    same = (len(my_word) == len(lib_word)
            and all(a.key() == b.key() for a, b in zip(my_word, lib_word))
            and cur == lib_M)
    all_match = all_match and same

check("2b l decreases by EXACTLY 2 at every synthesis step",
      all_drop2, f"{steps_total} steps over 130 words (incl. 30 high-sde)")
check("2c EXACTLY 1 of the 12 coset reps descends at every step",
      set(desc_counter) == {1},
      f"distribution of #descending-reps per step: {dict(desc_counter)}")
check("2d candidate l-distribution is {l-2:1, l:2, l+2:9} at every step",
      set(pattern_counter) == {((-2, 1), (0, 2), (2, 9))},
      f"observed patterns: {dict(pattern_counter)}")
check("2e instrumented walk reproduces synthesis_steps word and M exactly",
      all_match)


# ===========================================================================
# Section 3: geodesic / minimality vs BFS tree distance
# ===========================================================================
print()
print("=" * 76)
print("Section 3: word length == tree distance / 2 (BFS cross-check)")
print("=" * 76)

t0 = time.perf_counter()
levels, parent, cycles = bfs_tree(6)
t1 = time.perf_counter()
sizes = [len(l) for l in levels]
check("3a BFS to radius 6: level sizes 1,4,12,36,108,324,972 (=4*3^(r-1))",
      sizes == [1, 4, 12, 36, 108, 324, 972],
      f"sizes={sizes}  ({t1-t0:.1f}s)")
check("3b no cycle edges up to radius 6 (B is a tree locally)",
      cycles == [], f"{len(cycles)} cycle edges")

distmap = {}
for r, lev in enumerate(levels):
    for v in lev:
        distmap[v.key()] = r

# collect examples with sde 0..3 from seeded random words
rng = random.Random(SEED + 2)
examples = []
need = Counter()
while len(examples) < 24:
    w = ''.join(rng.choice('HSR') for _ in range(rng.randint(0, 10)))
    U = Mat.identity()
    for ch in w:
        U = U * GATES[ch]
    s = sde_scratch(U)
    if s <= 3 and need[s] < 8:
        need[s] += 1
        examples.append((w, U, s))

n_ok = 0
cover = Counter()
for w, U, s in examples:
    vk = Vertex("P", U).key()
    word, M = synthesis_steps(U)
    ok = (vk in distmap and distmap[vk] == ell_scratch(U) == 2 * len(word) == 2 * s)
    n_ok += ok
    cover[s] += 1
check("3c word length == d(e_0, U e_0)/2 == BFS distance/2 on 24 examples",
      n_ok == len(examples),
      f"{n_ok}/{len(examples)} examples, sde coverage {dict(sorted(cover.items()))}")

# the 12 coset-rep vertices are exactly the BFS sphere of radius 2 (S_0)
s0_keys_bfs = {v.key() for v in levels[2]}
s0_keys_reps = {Vertex("P", m * H_GATE).key() for m in REPS}
check("3d {m H e_0 : m in 12 reps} == BFS sphere(2) == S_0",
      s0_keys_reps == s0_keys_bfs and len(s0_keys_reps) == 12,
      f"|S_0| = {len(s0_keys_bfs)}")


# ===========================================================================
# Section 4: normal form & counting
# ===========================================================================
print()
print("=" * 76)
print("Section 4: normal form and exact counts")
print("=" * 76)

# 4a. stabilizer of e_1 inside U_3(O_F): the paper claims #G = 108 and
#     #U_3(O_F) e_1 = 1296/108 = 12 by orbit-stabilizer.
t0 = time.perf_counter()
e1_key = Vertex("P", H_GATE).key()
fiber = Counter(Vertex("P", m * H_GATE).key() for m in MONOS)
t1 = time.perf_counter()
check("4a stabilizer of e_1 in U_3(O_F) has exactly 108 elements (#G = 108)",
      fiber[e1_key] == 108, f"#G = {fiber[e1_key]}  ({t1-t0:.1f}s)")
check("4b orbit U_3(O_F).e_1 has exactly 12 vertices, all fibers = 108 = 1296/12",
      len(fiber) == 12 and set(fiber.values()) == {108}
      and set(fiber) == s0_keys_bfs,
      f"orbit size {len(fiber)}, fiber sizes {sorted(set(fiber.values()))}")

# 4c. l(H m2 H) distribution over the 12 reps: {0:1, 2:2, 4:9}
#     => exactly 9 'suitable' middle reps for sde-2 normal forms.
hmh = [(m, H_GATE * m * H_GATE) for m in REPS]
dist_hmh = Counter(ell(x) for _, x in hmh)
SUITABLE = [m for m, x in hmh if ell(x) == 4]
check("4c l(H m2 H) over the 12 reps is {0:1, 2:2, 4:9}; 9 suitable middles",
      dict(dist_hmh) == {0: 1, 2: 2, 4: 9},
      f"distribution {dict(sorted(dist_hmh.items()))}")

# --- fast enumeration helper: key of core*m3 without a full matmul --------
PERMS = list(itertools.permutations(range(3)))
PERMINVS = [tuple(p.index(c) for c in range(3)) for p in PERMS]


def all_product_keys(cores):
    """Keys of core * m3 over all 1296 monomials m3, for each core.
    Right-multiplying by m3 (rows[i][perm[i]] = u_i) permutes and unit-scales
    columns, so keys are assembled from precomputed scaled columns."""
    keyset = set()
    total = 0
    for core in cores:
        scaled = [[[(core.m[r][k] * u).key() for r in range(3)]
                   for u in _UNITS] for k in range(3)]
        for pinv in PERMINVS:
            for us in itertools.product(range(6), repeat=3):
                kt = tuple(scaled[pinv[c]][us[pinv[c]]][r]
                           for r in range(3) for c in range(3))
                keyset.add(kt)
                total += 1
    return keyset, total


# validate the fast-key trick against exact matrix products
rng = random.Random(SEED + 3)
n_ok = 0
for _ in range(20):
    core = REPS[rng.randrange(12)] * H_GATE
    pi_i = rng.randrange(6)
    us = tuple(rng.randrange(6) for _ in range(3))
    m3 = mono_from(PERMS[pi_i], us)
    pinv = PERMINVS[pi_i]
    scaled = [[[(core.m[r][k] * u).key() for r in range(3)]
               for u in _UNITS] for k in range(3)]
    fast = tuple(scaled[pinv[c]][us[pinv[c]]][r]
                 for r in range(3) for c in range(3))
    n_ok += (fast == (core * m3).key())
check("4d fast column-permute keys == exact matrix-product keys (20 samples)",
      n_ok == 20, f"{n_ok}/20")

# 4e. sde = 1: enumerate all m1 H m3, m1 in 12 reps, m3 in 1296 monomials.
cores1 = [m * H_GATE for m in REPS]
ok_core1 = all(ell(c) == 2 == ell_scratch(c) for c in cores1)
t0 = time.perf_counter()
keys1, total1 = all_product_keys(cores1)
t1 = time.perf_counter()
check("4e sde=1 normal forms m1*H*m3: 12*1296 = 15552 products, ALL distinct",
      ok_core1 and total1 == 15552 and len(keys1) == 15552,
      f"products={total1}, distinct={len(keys1)}  ({t1-t0:.1f}s)")

# 4f. sphere(4) via BFS == 108 == 4*3^3, and the 108 sde-2 cores hit exactly
#     those vertices, bijectively.
cores2 = []
for m1 in REPS:
    for m2 in SUITABLE:
        cores2.append(m1 * (H_GATE * m2 * H_GATE))
ok_core2 = all(ell(c) == 4 == ell_scratch(c) for c in cores2)
core2_vkeys = [Vertex("P", c).key() for c in cores2]
sphere4_keys = {v.key() for v in levels[4]}
check("4f #sphere(4) = 108 = 4*3^3 (BFS) and the 108 cores m1*H*m2*H hit "
      "each sphere-4 vertex exactly once",
      ok_core2 and len(sphere4_keys) == 108
      and len(cores2) == 108 and len(set(core2_vkeys)) == 108
      and set(core2_vkeys) == sphere4_keys,
      f"|sphere(4)|={len(sphere4_keys)}, cores={len(cores2)}, "
      f"distinct core vertices={len(set(core2_vkeys))}")

# 4g. sde = 2: full enumeration m1 H m2 H m3 over 12 x 9 x 1296.
t0 = time.perf_counter()
keys2, total2 = all_product_keys(cores2)
t1 = time.perf_counter()
check("4g sde=2 normal forms m1*H*m2*H*m3: 4*3^3*1296 = 139968 products, "
      "ALL distinct (normal form unique)",
      total2 == 139968 and len(keys2) == 139968 and 139968 == 108 * 1296,
      f"products={total2}, distinct={len(keys2)}  ({t1-t0:.1f}s)")

# 4h. random sde-1 products M0*H*M1 land in the sde=1 enumeration.
rng = random.Random(SEED + 4)
n_ok = 0
N1 = 1000
for _ in range(N1):
    U = MONOS[rng.randrange(1296)] * H_GATE * MONOS[rng.randrange(1296)]
    n_ok += (ell_scratch(U) == 2 and U.key() in keys1)
check("4h 1000 random M0*H*M1 all have sde=1 and lie in the 15552-element set",
      n_ok == N1, f"{n_ok}/{N1}")

# 4i. 3000 random sde-2 unitaries each occur in the 139968-element set
#     (occurrence is exactly once since the enumeration has no collisions).
t0 = time.perf_counter()
n_ok, n_drawn, n_acc = 0, 0, 0
while n_acc < 3000:
    U = (MONOS[rng.randrange(1296)] * H_GATE * MONOS[rng.randrange(1296)]
         * H_GATE * MONOS[rng.randrange(1296)])
    n_drawn += 1
    if ell_scratch(U) != 4:
        continue
    n_acc += 1
    n_ok += (U.key() in keys2)
t1 = time.perf_counter()
check("4i 3000 random sde-2 unitaries all occur (exactly once) in the "
      "normal-form enumeration",
      n_ok == 3000,
      f"{n_ok}/3000 found; acceptance {n_acc}/{n_drawn} "
      f"~ {n_acc/n_drawn:.3f} (theory 972/1296 = 0.750)  ({t1-t0:.1f}s)")


# ===========================================================================
# Section 5: edge cases
# ===========================================================================
print()
print("=" * 76)
print("Section 5: edge cases")
print("=" * 76)

I3 = Mat.identity()
cases = [
    ("I",        I3,                                            0),
    ("omega*I",  I3 * OMEGA,                                    0),
    ("-I",       -I3,                                           0),
    ("S",        S_GATE,                                        0),
    ("R",        R_GATE,                                        0),
    ("H",        H_GATE,                                        1),
    ("H^2",      H_GATE * H_GATE,                               0),
    ("SHSHS",    S_GATE * H_GATE * S_GATE * H_GATE * S_GATE,    None),
    ("HSH",      H_GATE * S_GATE * H_GATE,                      None),
    ("RHSHR",    R_GATE * H_GATE * S_GATE * H_GATE * R_GATE,    None),
]

all_ok = True
details = []
for name, U, expected_sde in cases:
    word, M = synthesis_steps(U)
    s = sde_scratch(U)
    ok = (reconstruct(word, M) == U and len(word) == s
          and is_monomial_scratch(M) and in_Gamma(U))
    if expected_sde is not None:
        ok = ok and (s == expected_sde)
    if name == "I":
        ok = ok and word == [] and M == I3
    if name == "H^2":
        ok = ok and is_monomial_scratch(U)   # H^2 is itself monomial
    all_ok = all_ok and ok
    details.append(f"{name}:sde={s}{'' if ok else '!FAIL'}")
check("5a edge cases synthesize exactly (I, wI, -I, S, R, H, H^2, SHSHS, ...)",
      all_ok, "; ".join(details))


# ===========================================================================
# Section 6: performance probe
# ===========================================================================
print()
print("=" * 76)
print("Section 6: performance probe (synthesis time vs sde)")
print("=" * 76)

rng = random.Random(SEED + 5)

perf_ok = True
perf_lines = []
for s in [5, 10, 20, 40]:
    U = random_unitary_with_sde(s, rng)
    t0 = time.perf_counter()
    word, M = synthesis_steps(U)
    t1 = time.perf_counter()
    ms = (t1 - t0) * 1000.0
    ok = (len(word) == s and reconstruct(word, M) == U
          and is_monomial_scratch(M))
    perf_ok = perf_ok and ok
    perf_lines.append(f"sde={s}: {ms:8.1f} ms total, {ms/s:6.1f} ms/step")
    print(f"    sde={s:3d}: synthesis {ms:8.1f} ms  ({ms/s:6.2f} ms per step), "
          f"round-trip {'OK' if ok else 'FAIL'}")
check("6a synthesis correct at sde = 5, 10, 20, 40 (timings above)",
      perf_ok, " | ".join(perf_lines))


# ===========================================================================
print()
print("=" * 76)
if FAILURES:
    print(f"RESULT: {len(FAILURES)} FAILURE(S): {FAILURES}")
    sys.exit(1)
print("RESULT: ALL CHECKS PASSED")
sys.exit(0)
