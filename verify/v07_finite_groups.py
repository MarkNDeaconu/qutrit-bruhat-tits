#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
v07_finite_groups.py -- verify the finite-group computations in the proof of
the main theorem of "Buildings for Synthesis with Clifford+R"
(paper.tex, subsection "Proof of main theorem", lines ~1007-1053).

Claims verified
  1. U_3(O_F) is exactly the 1296 monomial matrices:
     Diophantine step  m^2+n^2-mn <= 1  has the 7 listed solutions,
     #U_3(O_F) = 3! * 6^3 = 1296, and a *complete* enumeration shows any
     A in GL_3(O_F) with A*A = I is monomial.
  2. #G = 108 where G = U_3(O_F) ∩ H U_3(O_F) H^{-1}
     (= stabiliser of e_1 = H e_0 inside the stabiliser of the origin),
     computed two independent ways; brief group-structure report.
  3. Orbit-stabiliser: #(U_3(O_F) e_1) = 1296/108 = 12 and the orbit equals
     S_0 = {pure vertices at tree distance 2 from origin} (independent BFS).
  4. Stabiliser in U_3(O_F) of each alternating neighbour of the origin has
     order 324 = 1296/4 (transitivity on the 4 isotropic lines (+-1,+-1,+-1));
     amalgam data |G_P| = |G_A| = 1296, |G_E| = 324 (H fixes the midpoint
     alternating vertex and, with the 324 monomials, acts transitively on its
     4 pure neighbours).
  5. U_3(O_F) is contained in <H, S, R>, with explicit words.
  6. #{U in Gamma : l(U) = 2} = 1296^2/108 = 12 * 1296 = 15552 (full
     enumeration + collision-subgroup correspondence on samples).

btlib is treated as a tool, not ground truth: Section 0 independently
sanity-checks every primitive relied on (numpy floats for arithmetic and
unitarity, norm-based valuation for v_pi, invariance/faithfulness of the
canonical lattice key, and a from-scratch re-derivation of the local tree
structure at the origin and at an alternating vertex).

Deterministic (seeded RNG).  Prints PASS/FAIL per claim; exits non-zero on
any failure.
"""

import os
import sys
import time
import random
import itertools
from math import inf
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))

import numpy as np

from btlib import (Zw, Fw, Mat, ZERO, ONE, OMEGA, CHI, CHI_INV, I_SQRT3,
                   chi_pow, ell, d_tilde, H_GATE, S_GATE, R_GATE,
                   monomial_matrices, in_Gamma, lattice_contains, lattice_eq,
                   dual_basis, gram, is_self_dual, lattice_index_log3,
                   hnf_local, Vertex, ORIGIN, neighbors_of_pure,
                   neighbors_of_alternating, bfs_tree, isotropic_lines,
                   F3_LINES, _complete_to_unimodular,
                   synthesis_steps, reconstruct)

RNG = random.Random(20260610)
FAILURES = []


def check(name, cond, detail=""):
    status = "PASS" if cond else "FAIL"
    print(f"[{status}] {name}" + (f"  -- {detail}" if detail else ""))
    if not cond:
        FAILURES.append(name)
    return cond


W = OMEGA
W2 = OMEGA * OMEGA
I3 = Mat.identity()
Hinv = H_GATE.conjT()        # H is unitary, so H^-1 = H*
Sinv = S_GATE.conjT()

# ============================================================================
# Section 0: independent sanity checks of the btlib primitives we rely on
# ============================================================================
print("=" * 78)
print("Section 0: independent sanity checks of btlib primitives")
print("=" * 78)
t0 = time.time()

OMEGA_C = np.exp(2j * np.pi / 3)


def fw_to_c(x: Fw) -> complex:
    return (x.num.a + x.num.b * OMEGA_C) / x.den


def mat_to_c(g: Mat) -> np.ndarray:
    return np.array([[fw_to_c(x) for x in row] for row in g.m])


# 0a. gate definitions match the paper, numerically unitary
H_np = (1j / np.sqrt(3)) * np.array([[1, 1, 1],
                                     [1, OMEGA_C, OMEGA_C ** 2],
                                     [1, OMEGA_C ** 2, OMEGA_C]])
check("0a: H_GATE == (i/sqrt3)[[1,1,1],[1,w,w2],[1,w2,w]] numerically",
      np.allclose(mat_to_c(H_GATE), H_np, atol=1e-12))
check("0a: H,S,R numerically unitary and exactly unitary",
      all(np.allclose(mat_to_c(g).conj().T @ mat_to_c(g), np.eye(3),
                      atol=1e-12) and g.is_unitary()
          for g in (H_GATE, S_GATE, R_GATE)))
check("0a: Hinv = H* satisfies H*Hinv == I exactly", H_GATE * Hinv == I3)
check("0a: i/sqrt3 == (1+2w)/3 numerically",
      abs(fw_to_c(I_SQRT3) - 1j / np.sqrt(3)) < 1e-14)

# 0b. exact matrix products match numpy on random samples
MONOS = monomial_matrices()
ok = True
for _ in range(60):
    a = RNG.choice(MONOS)
    b = RNG.choice(MONOS)
    p = a * H_GATE * b
    if not np.allclose(mat_to_c(p), mat_to_c(a) @ H_np @ mat_to_c(b),
                       atol=1e-10):
        ok = False
check("0b: 60 random exact products m1*H*m2 match numpy float products", ok)


# 0c. v_pi cross-checked against the norm formula v_pi(x) = v_3(N(x)),
#     N(x) = x*conj(x) in Q  (valid since chi is the only prime above 3 and
#     v_pi restricted to Q equals 2*v_3).
def v3_fraction(q: Fraction):
    if q == 0:
        return inf
    num, den, v = q.numerator, q.denominator, 0
    while num % 3 == 0:
        num //= 3
        v += 1
    while den % 3 == 0:
        den //= 3
        v -= 1
    return v


def v_pi_indep(x: Fw):
    if x.is_zero():
        return inf
    return v3_fraction(Fraction(x.num.norm(), x.den ** 2))


ok = True
for _ in range(400):
    x = Fw(Zw(RNG.randint(-40, 40), RNG.randint(-40, 40)),
           RNG.choice([1, 1, 2, 3, 9, 27, 5, 15]))
    if x.v_pi() != v_pi_indep(x):
        ok = False
    if not x.is_zero() and x.conj().v_pi() != x.v_pi():
        ok = False
check("0c: v_pi == v_3(norm) on 400 random elements; v_pi(conj x)==v_pi(x)", ok)
check("0c: v_pi(3) == 2, v_pi(chi) == 1, v_pi(i/sqrt3) == -1",
      Fw(3).v_pi() == 2 and CHI.v_pi() == 1 and I_SQRT3.v_pi() == -1)


def ell_indep(g: Mat):
    return -2 * min(v_pi_indep(x) for row in g.m for x in row)


ok = True
for _ in range(200):
    p = RNG.choice(MONOS) * H_GATE * RNG.choice(MONOS)
    if ell(p) != ell_indep(p):
        ok = False
check("0c: ell == independent norm-based ell on 200 random m1*H*m2", ok)


# 0d. canonical lattice key (hnf_local): invariance under GL_3(O_pi) change
#     of basis, faithfulness (different lattice => different key), and
#     hnf_local(g) spans the same lattice as g.
def random_unimodular():
    """Random element of GL_3(O_pi) (product of elementary/permutation/unit
    matrices, occasionally with denominator 2 -- a unit in O_pi)."""
    U = Mat.identity()
    for _ in range(6):
        kind = RNG.randrange(3)
        if kind == 0:
            i, j = RNG.sample(range(3), 2)
            c = Fw(Zw(RNG.randint(-2, 2), RNG.randint(-2, 2)),
                   RNG.choice([1, 1, 1, 2]))
            E = Mat.identity()
            E.m[i][j] = c
            U = U * E
        elif kind == 1:
            perm = list(range(3))
            RNG.shuffle(perm)
            rows = [[ZERO] * 3 for _ in range(3)]
            for i in range(3):
                rows[i][perm[i]] = ONE
            U = U * Mat(rows)
        else:
            us = [RNG.choice([ONE, W, W2, -ONE, -W, -W2]) for _ in range(3)]
            U = U * Mat.diag(*us)
    return U


def random_basis():
    while True:
        g = Mat([[Zw(RNG.randint(-3, 3), RNG.randint(-3, 3))
                  for _ in range(3)] for _ in range(3)])
        if not g.det().is_zero():
            return g * chi_pow(RNG.randint(-1, 1))


ok_inv, ok_span, ok_neg = True, True, True
for _ in range(40):
    g = random_basis()
    U = random_unimodular()
    if hnf_local(g).key() != hnf_local(g * U).key():
        ok_inv = False
    if not lattice_eq(g, hnf_local(g)):
        ok_span = False
    g2 = g * Mat.diag(CHI, 1, 1)        # strictly smaller lattice
    if hnf_local(g).key() == hnf_local(g2).key() or lattice_eq(g, g2):
        ok_neg = False
check("0d: hnf key invariant under 40 random GL_3(O_pi) basis changes", ok_inv)
check("0d: hnf_local(g) spans the same lattice as g (40 random)", ok_span)
check("0d: distinct lattices get distinct keys (40 random negative tests)",
      ok_neg)
# faithfulness on the family used later: key equality <=> lattice equality
reps_by_key = {}
for m in RNG.sample(MONOS, 80):
    reps_by_key.setdefault(hnf_local(m * H_GATE).key(), []).append(m * H_GATE)
keys_list = list(reps_by_key)
ok = True
for k in keys_list:
    gs = reps_by_key[k]
    if len(gs) >= 2 and not lattice_eq(gs[0], gs[1]):
        ok = False
for k1, k2 in itertools.combinations(keys_list[:6], 2):
    if lattice_eq(reps_by_key[k1][0], reps_by_key[k2][0]):
        ok = False
check("0d: on the m*H family, key equality <=> lattice equality (sampled)", ok)

print(f"(section 0 took {time.time()-t0:.1f}s)")

# ============================================================================
# Section 1: U_3(O_F) = monomial matrices  (paper Proposition, lines 1019-1032)
# ============================================================================
print("=" * 78)
print("Section 1: U_3(O_F) = monomial matrices, #U_3(O_F) = 1296")
print("=" * 78)
t0 = time.time()

# 1a. Diophantine step (paper line 1025), exhaustive over |m|,|n| <= 3.
sols = [(m, n) for m in range(-3, 4) for n in range(-3, 4)
        if m * m + n * n - m * n <= 1]
expected = {(0, 0), (1, 0), (-1, 0), (0, 1), (0, -1), (1, 1), (-1, -1)}
check("1a: m^2+n^2-mn <= 1 has exactly the 7 listed solutions (|m|,|n|<=3)",
      set(sols) == expected and len(sols) == 7,
      f"found {sorted(sols)}")
# completeness beyond the box: m^2+n^2-mn = ((2m-n)^2+3n^2)/4 >= (3/4)max^2,
# so any solution has |m|,|n| <= 1; the box check below is consistent.
big = [(m, n) for m in range(-3, 4) for n in range(-3, 4)
       if max(abs(m), abs(n)) >= 2 and m * m + n * n - m * n <= 1]
check("1a: no solutions with max(|m|,|n|) >= 2 (so the box is exhaustive)",
      big == [])

# 1b. alpha in Z[omega] with N(alpha) <= 1:  exactly 0 and the 6 units.
alphas = [Zw(a, b) for a in range(-4, 5) for b in range(-4, 5)
          if Zw(a, b).norm() <= 1]
unit_set = {Zw(1, 0), Zw(-1, 0), Zw(0, 1), Zw(0, -1), Zw(-1, -1), Zw(1, 1)}
# Zw(0,1)=w, Zw(-1,-1)=w^2, Zw(1,1)=-w^2
check("1b: exactly 7 alpha with N(alpha)<=1: 0 and the 6 units +-1,+-w,+-w^2",
      len(alphas) == 7 and set(alphas) == unit_set | {Zw(0, 0)}
      and sum(1 for a in alphas if a.norm() == 1) == 6)

# 1c. unit-norm columns: a in O_F^3 with <a,a>=1.  Each entry has norm <= 1
#     (the norms are non-negative integers summing to 1), so entries range
#     over the 7-element set; the enumeration below is therefore COMPLETE.
cols = [c for c in itertools.product(alphas, repeat=3)
        if sum(x.norm() for x in c) == 1]
ok_shape = all(sum(1 for x in c if not x.is_zero()) == 1 and
               max(x.norm() for x in c) == 1 for c in cols)
check("1c: exactly 18 = 3*6 unit-norm columns, each = (unit)*e_i",
      len(cols) == 18 and ok_shape, f"#cols = {len(cols)}")

# 1d. COMPLETE enumeration of {A in GL_3(O_F) : A*A = I}: the diagonal of
#     A*A = I says each column is unit-norm, hence comes from the 18-set.
unitary_integral = []
for c1, c2, c3 in itertools.product(cols, repeat=3):
    A = Mat([[Fw(c1[i]), Fw(c2[i]), Fw(c3[i])] for i in range(3)])
    if A.conjT() * A == I3:
        unitary_integral.append(A)
mono_keys = {m.key() for m in MONOS}
ui_keys = {A.key() for A in unitary_integral}
check("1d: enumeration of all unitary integral matrices yields exactly 1296",
      len(unitary_integral) == 1296 and len(ui_keys) == 1296,
      f"found {len(unitary_integral)} from 18^3 = {18**3} column triples")
check("1d: they coincide with the 1296 = 3!*6^3 monomial matrices",
      ui_keys == mono_keys and len(mono_keys) == 1296
      and len(MONOS) == 6 * 6 * 6 * 6,
      "1296 == 3! * 6^3")
check("1d: every monomial matrix is exactly unitary and lies in Gamma",
      all(m.is_unitary() for m in MONOS)
      and all(in_Gamma(m) for m in RNG.sample(MONOS, 60)))
check("1d: every monomial fixes the origin vertex (is in GL_3(O_pi))",
      all(m.in_GL_O() for m in MONOS))
print("   (note: Gamma ∩ Stab(origin) = Gamma ∩ GL_3(O_pi) consists of")
print("    unitary matrices with entries in Z[1/3,w] ∩ O_pi = Z[w], hence")
print("    Stab(e_0, Gamma) = U_3(O_F) = the 1296 monomials, by 1d.)")
print(f"(section 1 took {time.time()-t0:.1f}s)")

# ============================================================================
# Section 2: the 108 claim
# ============================================================================
print("=" * 78)
print("Section 2: #G = #(U_3(O_F) ∩ H U_3(O_F) H^-1) = 108")
print("=" * 78)
t0 = time.time()

# count 1: m preserves the lattice H O^3, i.e. H^-1 m H in GL_3(O_pi)
G_list = [m for m in MONOS if (Hinv * m * H_GATE).in_GL_O()]
check("2a: #{m monomial : H^-1 m H in GL_3(O_pi)} == 108",
      len(G_list) == 108, f"count = {len(G_list)}")

# count 2: same thing via the vertex machinery (m fixes Vertex of H O^3)
e1v = Vertex("P", H_GATE)
G_keys = {m.key() for m in G_list}
G2_keys = {m.key() for m in MONOS if Vertex("P", m * H_GATE).key() == e1v.key()}
check("2b: #{m : m fixes the vertex H*O^3} == 108 and equals the 2a set",
      len(G2_keys) == 108 and G2_keys == G_keys)

# count 3: H^-1 m H is itself MONOMIAL (the literal intersection
# U_3(O_F) ∩ H U_3(O_F) H^-1); must agree since H^-1 m H is unitary, in
# Gamma, and integral => monomial by Section 1.
K_keys = {m.key() for m in MONOS if (Hinv * m * H_GATE).key() in mono_keys}
check("2c: #{m : H^-1 m H monomial} == 108 and equals the 2a set",
      len(K_keys) == 108 and K_keys == G_keys)

check("2d: orbit-stabilizer arithmetic 1296/108 == 12",
      1296 % 108 == 0 and 1296 // 108 == 12)

# group-structure report for G (order 108)
Ginv = {m.key(): m.conjT() for m in G_list}
ok_closed = True
comm_keys = {}
for a in G_list:
    ai = Ginv[a.key()]
    for b in G_list:
        p = a * b
        if p.key() not in G_keys:
            ok_closed = False
        c = ai * Ginv[b.key()] * p          # a^-1 b^-1 a b ... careful below
        # commutator [a,b] = a^-1 b^-1 a b:
        # ai * Ginv[b] * a * b = ai*(b^-1)*(a*b)?  p = a*b so c = ai*b^-1*a*b
        comm_keys[c.key()] = c
check("2e: G is closed under multiplication (a genuine subgroup)", ok_closed)


def order_of(m):
    p, k = m, 1
    while not p == I3:
        p = p * m
        k += 1
        assert k <= 200
    return k


orders = {}
for m in G_list:
    orders[order_of(m)] = orders.get(order_of(m), 0) + 1
center = [z for z in G_list
          if all((z * h).key() == (h * z).key() for h in G_list)]
scalar_keys = {(Mat.diag(u, u, u)).key() for u in (ONE, W, W2, -ONE, -W, -W2)}

# derived subgroup: close the commutator set under multiplication
derived = dict(comm_keys)
frontier = list(derived.values())
while frontier:
    new = []
    for x in frontier:
        for gkey, g in list(comm_keys.items()):
            p = x * g
            if p.key() not in derived:
                derived[p.key()] = p
                new.append(p)
    frontier = new
diag_in_G = [m for m in G_list
             if all(m.m[i][j].is_zero() for i in range(3) for j in range(3)
                    if i != j)]
exponent = 1
for o in orders:
    from math import lcm
    exponent = lcm(exponent, o)
check("2f: structure report computed", True)
print(f"    G: |G| = {len(G_list)} = 2^2*3^3, exponent {exponent}")
print(f"    order histogram: {dict(sorted(orders.items()))}")
print(f"    |Z(G)| = {len(center)}, Z(G) == scalar subgroup <wI,-I>: "
      f"{ {z.key() for z in center} == scalar_keys }")
print(f"    |[G,G]| = {len(derived)}, |G ∩ diagonal| = {len(diag_in_G)}")
print(f"(section 2 took {time.time()-t0:.1f}s)")

# ============================================================================
# Section 3: orbit U_3(O_F)*e_1 has size 12 and equals S_0
# ============================================================================
print("=" * 78)
print("Section 3: orbit of e_1 = H*O^3 under U_3(O_F) is S_0, #S_0 = 12")
print("=" * 78)
t0 = time.time()

orbit = {}
for m in MONOS:
    orbit.setdefault(Vertex("P", m * H_GATE).key(), []).append(m)
check("3a: orbit of vertex H*O^3 under the 1296 monomials has exactly 12 "
      "elements", len(orbit) == 12, f"#orbit = {len(orbit)}")
check("3b: each orbit fibre has size 108 (= #G, consistent with 2)",
      sorted(len(v) for v in orbit.values()) == [108] * 12)

levels, parent, cyc = bfs_tree(2)
S0 = levels[2]
check("3c: BFS: 4 alternating at distance 1, 12 pure at distance 2, no "
      "cycles", len(levels[1]) == 4 and all(v.kind == "A" for v in levels[1])
      and len(S0) == 12 and all(v.kind == "P" for v in S0) and cyc == [],
      f"|level1| = {len(levels[1])}, |level2| = {len(S0)}, cycles = {len(cyc)}")
kids = {}
for v in S0:
    kids[parent[v.key()]] = kids.get(parent[v.key()], 0) + 1
check("3d: each of the 4 alternating neighbours contributes exactly 3 new "
      "pure vertices (12 = 4 x 3)", sorted(kids.values()) == [3, 3, 3, 3])
check("3e: orbit key set == S_0 key set (and e_1 is in S_0)",
      set(orbit) == {v.key() for v in S0} and e1v.key() in set(orbit))
check("3f: d~(I, m*H) == 2 for one rep m per orbit element; d(e_0,e_1) == 2",
      all(d_tilde(I3, ms[0] * H_GATE) == 2 for ms in orbit.values())
      and d_tilde(I3, H_GATE) == 2)
ok = True
for v in S0:
    if not (is_self_dual(v.g) and lattice_eq(v.g, dual_basis(v.g))):
        ok = False
check("3g: every S_0 basis is self-dual (Gram test AND independent dual "
      "test)", ok)
print(f"(section 3 took {time.time()-t0:.1f}s)")

# ============================================================================
# Section 4: edge stabilisers, local transitivity, amalgam data
# ============================================================================
print("=" * 78)
print("Section 4: alternating-vertex stabilisers and amalgam structure")
print("=" * 78)
t0 = time.time()

alts = neighbors_of_pure(ORIGIN)
alt_keys = [a.key() for a in alts]
check("4a: origin has exactly 4 alternating neighbours (distinct keys)",
      len(alts) == 4 and len(set(alt_keys)) == 4)
iso = isotropic_lines([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
check("4b: isotropic lines of the residue form are the 4 lines "
      "(+-1,+-1,+-1)",
      set(iso) == {(1, 1, 1), (1, 1, 2), (1, 2, 1), (1, 2, 2)},
      f"lines = {sorted(iso)}")

# independent re-derivation of the alternating neighbours at the origin:
# every alternating neighbour is a lattice L with O^3 ⊂ L ⊂ chi^-1 O^3 and
# [L:O^3] = 3; such L correspond to the 13 lines of F_3^3.  The vertex
# condition (paper's simplex classification) is the CHAIN condition
#   pi*L ⊆ L# ⊊ L
# (note L# ⊊ L is automatic here, since L ⊇ O^3 forces L# ⊆ O^3 ⊆ L).
indep_alt_keys = set()
n_candidates, n_weak = 0, 0
for u in F3_LINES:
    B = _complete_to_unimodular(u) * Mat.diag(CHI_INV, 1, 1)
    assert lattice_contains(B, I3) and lattice_index_log3(B, I3) == 1
    assert lattice_contains(Mat.diag(CHI_INV, CHI_INV, CHI_INV), B)
    Bd = dual_basis(B)
    if lattice_contains(B, Bd) and not lattice_eq(B, Bd):
        n_weak += 1                               # holds for all 13
        assert lattice_index_log3(B, Bd) == 2     # [L : L#] = 9
    if lattice_contains(Bd, B * CHI):             # chain: pi*L ⊆ L#
        n_candidates += 1
        indep_alt_keys.add(Vertex("A", B).key())
check("4c: from-scratch enumeration over the 13 lines: all 13 satisfy "
      "L# ⊊ L but exactly 4 satisfy the chain condition pi*L ⊆ L#, "
      "matching btlib's 4 alternating neighbours",
      n_weak == 13 and n_candidates == 4 and indep_alt_keys == set(alt_keys),
      f"{n_candidates} of 13 satisfy pi*L ⊆ L# (L# ⊊ L holds for {n_weak})")

stabs = []
for A in alts:
    stab = [m for m in MONOS if Vertex("A", m * A.g).key() == A.key()]
    stabs.append(stab)
check("4d: the monomial stabiliser of EACH alternating neighbour has order "
      "324 = 1296/4", [len(s) for s in stabs] == [324] * 4,
      f"orders = {[len(s) for s in stabs]}")
orb_alt = {Vertex("A", m * alts[0].g).key() for m in MONOS}
check("4e: U_3(O_F) is transitive on the 4 alternating neighbours",
      orb_alt == set(alt_keys), f"#orbit = {len(orb_alt)}")

# the midpoint alternating vertex A* between e_0 and e_1 is fixed by H
A_star = None
for A in alts:
    pn = neighbors_of_alternating(A)
    if e1v.key() in {p.key() for p in pn}:
        A_star = A
        pn_star = pn
        break
check("4f: exactly one alternating neighbour A* of the origin is adjacent "
      "to e_1", A_star is not None)
check("4g: H fixes A* (so Stab(A*,Gamma) is strictly bigger than its 324 "
      "monomials)", Vertex("A", H_GATE * A_star.g).key() == A_star.key())

# independent re-derivation of the pure neighbours of A*: middle lattices
# A*# ⊂ L ⊂ A* (index 3 each way), via the 13 lines in chi^-1 A*# / A*#.
Ad = dual_basis(A_star.g)
indep_pn_keys = set()
n_mid = 0
for u in F3_LINES:
    L = Ad * _complete_to_unimodular(u) * Mat.diag(CHI_INV, 1, 1)
    if not lattice_contains(A_star.g, L):
        continue
    n_mid += 1
    assert lattice_contains(L, Ad) and lattice_index_log3(L, Ad) == 1
    sd = is_self_dual(L) and lattice_eq(L, dual_basis(L))
    assert sd, "middle lattice not self-dual?!"
    indep_pn_keys.add(Vertex("P", L).key())
check("4h: from-scratch middle-lattice enumeration gives exactly the 4 pure "
      "neighbours of A* (all self-dual), matching btlib",
      n_mid == 4 and indep_pn_keys == {p.key() for p in pn_star}
      and ORIGIN.key() in indep_pn_keys and e1v.key() in indep_pn_keys)

# orbit of the origin under <monomial stab of A*, H> inside pure nbrs of A*
stab_star = stabs[alts.index(A_star)]
pn_keys = {p.key() for p in pn_star}
verts = {ORIGIN.key(): ORIGIN}
frontier = [ORIGIN]
while frontier:
    new = []
    for v in frontier:
        for g in stab_star + [H_GATE]:
            w = Vertex("P", g * v.g)
            if w.key() not in verts:
                verts[w.key()] = w
                new.append(w)
    frontier = new
check("4i: <Stab_monomial(A*), H> acts transitively on the 4 pure "
      "neighbours of A*", set(verts) == pn_keys and len(verts) == 4,
      f"orbit size = {len(verts)}")
print("    amalgam data confirmed: |G_P| = |Stab(e_0)| = 1296,")
print("    |G_E| = |Stab(edge)| = |Stab(e_0) ∩ Stab(A*)| = 324,")
print(f"    |G_A| = |Stab(A*)| = (orbit 4) x (edge stab 324) = "
      f"{4*324} = 1296;  indices 4 and 4 <-> (4,4)-biregular tree.")
print(f"(section 4 took {time.time()-t0:.1f}s)")

# ============================================================================
# Section 5: U_3(O_F) ⊆ <H, S, R>
# ============================================================================
print("=" * 78)
print("Section 5: all 1296 monomials lie in <H, S, R> (with explicit words)")
print("=" * 78)
t0 = time.time()

TOKEN = {"H": H_GATE, "h": Hinv, "S": S_GATE, "s": Sinv, "R": R_GATE}


def eval_word(word):
    out = I3
    for t in word:
        out = out * TOKEN[t]
    return out


P23 = Mat([[1, 0, 0], [0, 0, 1], [0, 1, 0]])
check("5a: H^2 == -P23 (hence H has order 4 and H^-1 = H^3)",
      H_GATE * H_GATE == -P23 and eval_word("HHHH") == I3)

# BFS over words, collecting monomial values
GENS5 = [(H_GATE, "H"), (Hinv, "h"), (S_GATE, "S"), (Sinv, "s"),
         (R_GATE, "R")]
seen = {I3.key(): ""}
frontier = [(I3, "")]
mono_found = {I3.key(): (I3, "")}
level_sizes = []
MAXDEPTH = 6
for depth in range(1, MAXDEPTH + 1):
    nf = []
    for (m, w) in frontier:
        for (g, t) in GENS5:
            p = m * g
            k = p.key()
            if k in seen:
                continue
            seen[k] = w + t
            nf.append((p, w + t))
            if k in mono_keys:
                mono_found[k] = (p, w + t)
    level_sizes.append(len(nf))
    frontier = nf
print(f"    BFS over words: level sizes {level_sizes}, "
      f"{len(seen)} distinct matrices, {len(mono_found)} monomials found")


def close_words(gens):
    elems = {I3.key(): (I3, "")}
    frontier = [(I3, "")]
    while frontier:
        new = []
        for (m, w) in frontier:
            for (g, gw) in gens:
                p = m * g
                k = p.key()
                if k not in elems:
                    elems[k] = (p, w + gw)
                    new.append((p, w + gw))
        frontier = new
    return elems


# greedy generating set from the found monomials (shortest words first)
gens_sel = []
closed = {I3.key(): (I3, "")}
for (m, w) in sorted(mono_found.values(), key=lambda t: len(t[1])):
    if w and m.key() not in closed:
        gens_sel.append((m, w))
        closed = close_words(gens_sel)
print(f"    after closure of BFS monomials: subgroup of order {len(closed)}")

# augment with H-conjugates that are monomial (still words in H,S,R)
rounds = 0
while len(closed) < 1296 and rounds < 8:
    rounds += 1
    new = []
    for (m, w) in list(closed.values()):
        for (c, cw) in ((H_GATE * m * Hinv, "H" + w + "h"),
                        (Hinv * m * H_GATE, "h" + w + "H")):
            if c.key() in mono_keys and c.key() not in closed:
                new.append((c, cw))
    if not new:
        break
    for (m, w) in sorted(new, key=lambda t: len(t[1])):
        if m.key() not in closed:
            gens_sel.append((m, w))
            closed = close_words(gens_sel)
    print(f"    augmentation round {rounds}: subgroup order now {len(closed)}")

check("5b: the monomial subgroup generated by <H,S,R>-words has order "
      "exactly 1296 (i.e. U_3(O_F) ⊆ <H,S,R>)",
      len(closed) == 1296 and set(closed) == mono_keys,
      f"order = {len(closed)}")
print("    generating monomial words used:")
for (m, w) in gens_sel:
    print(f"      length {len(w):3d}: {w}")

targets = {
    "-I": Mat.diag(-1, -1, -1),
    "omega*I": Mat.diag(W, W, W),
    "3-cycle (e0->e1->e2->e0)": Mat([[0, 0, 1], [1, 0, 0], [0, 1, 0]]),
    "P12": Mat([[0, 1, 0], [1, 0, 0], [0, 0, 1]]),
    "diag(w,1,1)": Mat.diag(W, 1, 1),
}
ok = True
for name, tgt in targets.items():
    k = tgt.key()
    if k not in closed or eval_word(closed[k][1]) != tgt:
        ok = False
        print(f"    MISSING/WRONG word for {name}")
    else:
        w = closed[k][1]
        print(f"    {name}: word of length {len(w)}"
              + (f": {w}" if len(w) <= 60 else f" (first 60: {w[:60]}...)"))
check("5c: explicit verified words exist for -I, wI, the 3-cycle, P12, "
      "diag(w,1,1)", ok)

sample = RNG.sample(list(closed.values()), 25)
check("5d: 25 random closure words re-evaluate to their matrices exactly",
      all(eval_word(w) == m for (m, w) in sample))
print(f"(section 5 took {time.time()-t0:.1f}s)")

# ============================================================================
# Section 6: # of sde-1 unitaries in Gamma is 12 * 1296 = 15552
# ============================================================================
print("=" * 78)
print("Section 6: #{U in Gamma : l(U) = 2} = 1296^2/108 = 15552")
print("=" * 78)
t0 = time.time()

check("6a: arithmetic 1296^2 / 108 == 12 * 1296 == 15552",
      1296 ** 2 // 108 == 15552 and 1296 ** 2 % 108 == 0)

reps = [orbit[k][0] for k in orbit]            # 12 coset reps of G
# disjointness of the 12 sets r_i H M  <=>  H^-1 r_j^-1 r_i H not in GL(O_pi)
bad = sum(1 for i in range(12) for j in range(12) if i != j
          and (Hinv * reps[j].conjT() * reps[i] * H_GATE).in_GL_O())
check("6b: the 12 coset sets r_i H U_3(O_F) are pairwise disjoint "
      "(132 checks)", bad == 0, f"violations = {bad}")

prod_keys = set()
ell_ok = True
for r in reps:
    rH = r * H_GATE
    for m2 in MONOS:
        p = rH * m2
        prod_keys.add(p.key())
        if ell(p) != 2:
            ell_ok = False
check("6c: FULL enumeration: exactly 15552 distinct products m1*H*m2, "
      "all with l = 2", len(prod_keys) == 15552 and ell_ok,
      f"#distinct = {len(prod_keys)}")

ok_sample = True
for _ in range(2000):
    m1, m2 = RNG.choice(MONOS), RNG.choice(MONOS)
    p = m1 * H_GATE * m2
    if ell(p) != 2 or p.key() not in prod_keys:
        ok_sample = False
check("6d: 2000 random products m1*H*m2 all have sde 1 and land in the "
      "enumerated set", ok_sample)
check("6d': 30 random products are exactly unitary and in Gamma",
      all(in_Gamma(RNG.choice(MONOS) * H_GATE * RNG.choice(MONOS))
          for _ in range(30)))

# collision-subgroup correspondence:  m3 H m4 = U  <=>  m3 in m1*G
ok_coll = True
for trial in range(3):
    m1, m2 = RNG.choice(MONOS), RNG.choice(MONOS)
    U = m1 * H_GATE * m2
    direct = []
    for m3 in MONOS:
        m4 = Hinv * m3.conjT() * U
        if m4.key() in mono_keys:
            direct.append((m3, m4))
    predicted_m3 = {(m1 * g.conjT()).key() for g in G_list}
    if (len(direct) != 108
            or {m3.key() for (m3, _) in direct} != predicted_m3
            or any(not (m3 * H_GATE * m4 == U)
                   for (m3, m4) in direct[:6])):
        ok_coll = False
check("6e: for 3 random U = m1 H m2, exactly 108 collision pairs "
      "(m3,m4), m3 ranging over the coset m1*G", ok_coll)

# completeness: every U in Gamma with l(U)=2 is some m1 H m2.
# Reason (machine-checked ingredients): for unitary U, l(U^-1) = l(U*) = l(U)
# (0c: v_pi(conj) = v_pi), so d~(I,U) = l(U) = 2, so U*v0 is in S_0, which
# equals the monomial orbit of e_1 (3e); then (mH)^-1 U fixes v0 and lies in
# Gamma, hence is monomial (1d).
ok = True
for _ in range(200):
    U = RNG.choice(MONOS) * H_GATE * RNG.choice(MONOS)
    if ell(U.conjT()) != ell(U) or d_tilde(I3, U) != 2:
        ok = False
check("6f: l(U*) == l(U) and d~(I,U) == 2 on 200 sampled sde-1 unitaries "
      "(completeness ingredients; with 3e and 1d => count is exactly 15552)",
      ok)

ok_syn = True
for _ in range(20):
    U = RNG.choice(MONOS) * H_GATE * RNG.choice(MONOS)
    word, M = synthesis_steps(U)
    if len(word) != 1 or not (reconstruct(word, M) == U):
        ok_syn = False
check("6g: exact synthesis writes 20 sampled sde-1 unitaries as (m H) * "
      "monomial (word length 1)", ok_syn)
print(f"(section 6 took {time.time()-t0:.1f}s)")

# ============================================================================
print("=" * 78)
if FAILURES:
    print(f"OVERALL: FAIL ({len(FAILURES)} failed checks):")
    for f in FAILURES:
        print(f"  - {f}")
    sys.exit(1)
print("OVERALL: PASS (all checks)")
sys.exit(0)
