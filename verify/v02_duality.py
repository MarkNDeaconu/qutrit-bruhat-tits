#!/usr/bin/env python
"""
v02_duality.py -- Verification of dual-lattice facts in
"Buildings for Synthesis with Clifford+R", paper Sections 2.5-2.6
(paper.tex lines ~272-323).

Claims verified
---------------
 0. Sanity of btlib primitives used here (independent v_pi, residue facts,
    exact matrix inverse vs numpy floats, hnf invariance, lattice_eq).
 1. Dual basis = (A*)^{-1}; double dual == original lattice (as lattices);
    membership criterion x in dual(L) <=> <x, basis col> in O_pi  (line 282).
 2. Proposition le:dual_equivalence, bullet 1: L1 < L2 => dual(L2) < dual(L1).
 3. Proposition le:dual_equivalence, bullet 2: pi-equivalence preserved by
    dualization (with the exact exponent dual(chi^k L) = chi^{-k} dual(L)).
 4. det identity det(dual L) = conj(det L)^{-1} as O_pi-modules (line 315),
    plus v_pi(conj x) = v_pi(x) (line 317).
 5. Lemma le:dual_scaling (rank 3): dual(L) = pi^i L  =>  i even and
    pi^{i/2} L self-dual.  Constructed examples + random scan + alternating
    vertices are never pi-equivalent to their dual.
 6. Remark after le:dual_scaling ("fails for even dimensions"), investigated:
    6a. -1 is NOT a norm of F_pi over Q_3 (brute force mod 3, 9, 27);
        norms of units are == 1 mod 3.
    6b. The standard 2-dim Hermitian form x conj(x) + y conj(y) is
        ANISOTROPIC over F_pi (exhaustive primitive residue scan).
    6c. Dim 2, standard form: odd exponent i in dual(L) = pi^i L is
        IMPOSSIBLE.  Proof chain verified computationally:
        (i)  every rank-2 lattice is k diag(chi^a, chi^b) O^2 with
             k in GL_2(O)  (own Cartan decomposition, verified exactly);
        (ii) scaling shifts i by even amounts, so wlog L = k diag(chi^s,1)O^2
             with i = -s; for s odd the diagonal Gram entry
             (chi^{-s} Gram)_{22} (s>0) or (chi^{-s} Gram)_{11} (s<0) has
             negative valuation because Q(primitive column of k) is a UNIT
             (this is exactly anisotropy mod pi from 6b)  => no odd example.
        (iii) skew-unimodular residue lemma: N in GL_2(O), N* = -w^2 N
             forces N == antidiag(b,-b) mod pi and det N == 1 mod 3,
             while the Gram/norm side forces det == 2 mod 3 (independent
             confirmation of the obstruction).
    6d. Dim 2 with the HYPERBOLIC form x1 conj(y2) + x2 conj(y1): odd
        example EXISTS (L = diag(1, chi^-1) O^2 has dual_J(L) = pi L).
    6e. Dim 4, standard form: explicit ODD example exists:
        L = O^4 + chi^{-1} v1 O + chi^{-1} v2 O,  v1=(1,1,1,0), v2=(0,1,-1,1)
        (a totally isotropic plane mod 3) has dual(L) = pi L, i = 1 odd,
        and no pi-rescaling of L is self-dual.
    => the Remark's "fails for even dimensions" is TRUE for dim 4 with the
       standard form, but FALSE for dim 2 with the standard form (where the
       lemma's conclusion still holds); in dim 2 failure requires a
       different (isotropic) form.

Deterministic: seeded RNG.  Exit code 0 iff every check passes.
"""

import sys
import random
from math import inf

sys.path.insert(0, "/Users/markdeaconu/projects/qutrits_v2")

import numpy as np

from btlib import (Zw, Fw, Mat, ZERO, ONE, OMEGA, CHI, CHI_INV, chi_pow,
                   dual_basis, gram, lattice_eq, lattice_contains,
                   lattice_index_log3, is_self_dual, hnf_local,
                   H_GATE, S_GATE, R_GATE,
                   Vertex, ORIGIN, neighbors)

SEED = 20260609
rng = random.Random(SEED)

RESULTS = []


def check(name, ok, detail=""):
    RESULTS.append((name, bool(ok), detail))
    line = f"{'PASS' if ok else 'FAIL'}: {name}"
    if detail:
        line += f"  [{detail}]"
    print(line)
    return bool(ok)


# ---------------------------------------------------------------------------
# Independent (from-scratch) primitives.  These do NOT call btlib's v_pi.
# ---------------------------------------------------------------------------

def my_v3(n):
    """3-adic valuation of a nonzero integer."""
    n = abs(int(n))
    assert n != 0
    v = 0
    while n % 3 == 0:
        n //= 3
        v += 1
    return v


def my_vpi(x: Fw):
    """Independent v_pi:  v_pi(x) = v_3(N(x)) for the ramified prime over 3,
    where N((a+bw)/d) = (a^2 - ab + b^2)/d^2  (integer norm computation)."""
    a, b, d = x.num.a, x.num.b, x.den
    if a == 0 and b == 0:
        return inf
    return my_v3(a * a - a * b + b * b) - 2 * my_v3(d)


def my_integral(x: Fw):
    return my_vpi(x) >= 0


def inner(x, y):
    """<x,y> = sum x_i conj(y_i), exact."""
    return sum((xi * yi.conj() for xi, yi in zip(x, y)), ZERO)


def matvec(M: Mat, v):
    return [sum((M.m[i][k] * v[k] for k in range(M.n)), ZERO)
            for i in range(M.n)]


def det_int(rows):
    """Exact integer determinant by minor expansion (small matrices)."""
    n = len(rows)
    if n == 1:
        return rows[0][0]
    tot = 0
    for j in range(n):
        if rows[0][j] == 0:
            continue
        minor = [r[:j] + r[j + 1:] for r in rows[1:]]
        tot += (-1) ** j * rows[0][j] * det_int(minor)
    return tot


# ---------------------------------------------------------------------------
# Random generators (deterministic via rng)
# ---------------------------------------------------------------------------

DENS = [1, 1, 2, 3, 5, 6, 7, 9, 15]   # mixed denominators incl. non-3 primes


def rand_Fw(lim=6):
    return Fw(Zw(rng.randint(-lim, lim), rng.randint(-lim, lim)),
              rng.choice(DENS))


def rand_basis(n=3):
    while True:
        A = Mat([[rand_Fw() for _ in range(n)] for _ in range(n)])
        if not A.det().is_zero():
            return A


def rand_GLO(n=3):
    """Random element of GL_n(O_pi) with entries in Z[omega]."""
    while True:
        M = Mat([[Fw(Zw(rng.randint(-2, 2), rng.randint(-2, 2)))
                  for _ in range(n)] for _ in range(n)])
        d = M.det()
        if (not d.is_zero()) and d.v_pi() == 0:
            return M


def rand_int_mat(n=3):
    """Random integral (Z[omega]-entry) nonsingular matrix (any det)."""
    while True:
        M = Mat([[Fw(Zw(rng.randint(-3, 3), rng.randint(-3, 3)))
                  for _ in range(n)] for _ in range(n)])
        if not M.det().is_zero():
            return M


def rand_int_vec(n, lim=4):
    return [Fw(Zw(rng.randint(-lim, lim), rng.randint(-lim, lim)))
            for _ in range(n)]


W_C = complex(-0.5, 3 ** 0.5 / 2)     # omega as a float


def to_c(M: Mat):
    return np.array([[(e.num.a + e.num.b * W_C) / e.den for e in row]
                     for row in M.m])


# ===========================================================================
# Part 0: sanity of btlib primitives (btlib is a tool, not ground truth)
# ===========================================================================

def part0():
    print("\n--- Part 0: btlib primitive sanity checks ---")

    # v_pi against independent norm-based implementation
    n_tested = 0
    ok = True
    for _ in range(500):
        x = rand_Fw(lim=20)
        if x.is_zero():
            continue
        n_tested += 1
        if x.v_pi() != my_vpi(x):
            ok = False
            print("  mismatch:", x, x.v_pi(), my_vpi(x))
    # targeted values
    for val, expect in [(CHI, 1), (Fw(3), 2), (Fw(1, 3), -2), (OMEGA, 0),
                        (Fw(1, 2), 0), (Fw(1, 5), 0), (Fw(Zw(1, -1), 9), -3),
                        (I := Fw(Zw(1, 2), 3), -1)]:
        n_tested += 1
        if not (val.v_pi() == my_vpi(val) == expect):
            ok = False
            print("  targeted mismatch:", val, val.v_pi(), my_vpi(val), expect)
    check("P0.1 btlib v_pi == independent norm-based v_pi",
          ok, f"{n_tested} samples incl. targeted values; v_pi(3)=2, v_pi(chi)=1")

    # v_pi(conj x) == v_pi(x)  (used in proof of le:dual_scaling, line 317)
    ok = all((x := rand_Fw(lim=15)).is_zero() or
             my_vpi(x.conj()) == my_vpi(x) for _ in range(300))
    check("P0.2 v_pi(conj x) == v_pi(x) (conj fixes pi)", ok, "300 samples")

    # residue map: ring hom, omega == 1, conj trivial on residue field
    ok = True
    for _ in range(300):
        x, y = rand_Fw(), rand_Fw()
        if not (x.is_integral() and y.is_integral()):
            continue
        ok &= ((x + y).res3() == (x.res3() + y.res3()) % 3)
        ok &= ((x * y).res3() == (x.res3() * y.res3()) % 3)
        ok &= (x.conj().res3() == x.res3())
        ok &= ((x.res3() == 0) == (my_vpi(x) >= 1 or x.is_zero()))
    ok &= (OMEGA.res3() == 1) and (CHI.res3() == 0)
    check("P0.3 res3 is a ring hom to F_3, omega==1, conj trivial, "
          "ker = pi", ok, "300 samples + omega/chi")

    # exact matrix algebra vs numpy floats
    ok = True
    for _ in range(20):
        for n in (2, 3, 4):
            A = rand_basis(n)
            Ac = to_c(A)
            ok &= np.allclose(to_c(A.inv()), np.linalg.inv(Ac), atol=1e-8)
            dc = A.det()
            ok &= abs((dc.num.a + dc.num.b * W_C) / dc.den
                      - np.linalg.det(Ac)) < 1e-6 * max(1, abs(np.linalg.det(Ac)))
            ok &= np.allclose(to_c(A.conjT()), Ac.conj().T, atol=1e-12)
            ok &= np.allclose(to_c(dual_basis(A)),
                              np.linalg.inv(Ac.conj().T), atol=1e-7)
            ok &= (A * A.inv() == Mat.identity(n))
    check("P0.4 Mat inv/det/conjT/dual_basis agree with numpy; "
          "A*A^-1 == I exactly", ok, "20 trials x dims 2,3,4")

    # lattice_eq positive and negative controls
    okp = okn = True
    for _ in range(50):
        A = rand_basis(3)
        U = rand_GLO(3)
        okp &= lattice_eq(A, A * U)
        okn &= not lattice_eq(A, A * Mat.diag(CHI, 1, 1))
        okn &= not lattice_eq(A, A * CHI)
    check("P0.5 lattice_eq: A*O^3 == A*U*O^3 (U in GL(O)); "
          "differs under chi-scalings", okp and okn, "50 trials each control")

    # lattice_contains spot check via explicit points
    ok = True
    for _ in range(30):
        B2 = rand_basis(3)
        T = rand_int_mat(3)
        B1 = B2 * T
        ok &= lattice_contains(B2, B1)
        x = matvec(B1, rand_int_vec(3))
        coords = matvec(B2.inv(), x)
        ok &= all(my_integral(c) or c.is_zero() for c in coords)
    check("P0.6 lattice_contains consistent with explicit point coordinates",
          ok, "30 trials, 1 random point each")

    # hnf_local: canonical and lattice-preserving
    ok = True
    for _ in range(30):
        A = rand_basis(3)
        U = rand_GLO(3)
        ok &= lattice_eq(hnf_local(A), A)
        ok &= hnf_local(A * U).key() == hnf_local(A).key()
    check("P0.7 hnf_local preserves the lattice and is GL(O)-invariant",
          ok, "30 trials")


# ===========================================================================
# Part 1: dual basis, double dual, membership criterion (100 random bases)
# ===========================================================================

def part1():
    print("\n--- Part 1: dual basis (A*)^-1, double dual, membership ---")
    N = 100
    ok_dd = ok_wd = ok_mem_pos = ok_mem_neg = ok_lat_pt = True
    for t in range(N):
        A = rand_basis(3)
        U1, U2 = rand_GLO(3), rand_GLO(3)
        B = A * U1                      # same lattice, scrambled basis
        D = dual_basis(B)               # (B*)^-1

        # well-definedness: dual computed from different bases agrees
        ok_wd &= lattice_eq(D, dual_basis(A))

        # double dual == original lattice (scramble dual basis first so the
        # test is lattice equality, not trivial matrix equality)
        Dscr = D * U2
        ok_dd &= lattice_eq(dual_basis(Dscr), A)
        ok_dd &= not (dual_basis(Dscr) == A)  # really a lattice-level test

        # membership, positive direction: x = D c  (c integral)
        c = rand_int_vec(3)
        x = matvec(D, c)
        for j in range(3):
            val = inner(x, A.col(j))
            ok_mem_pos &= (val.is_zero() or my_integral(val))
        # ... and against a random lattice point y = A c'
        y = matvec(A, rand_int_vec(3))
        v2 = inner(x, y)
        ok_lat_pt &= (v2.is_zero() or my_integral(v2))

        # membership, negative direction: w has one non-integral coordinate
        j0 = rng.randrange(3)
        w = rand_int_vec(3)
        w[j0] = w[j0] + CHI_INV                      # v_pi = -1
        xbad = matvec(dual_basis(A), w)
        valbad = inner(xbad, A.col(j0))
        # the defining identity gives <xbad, a_j> = w_j exactly
        ok_mem_neg &= (valbad == w[j0]) and (my_vpi(valbad) < 0)
    check("P1.1 dual well-defined: dual of A*O^3 independent of basis",
          ok_wd, f"{N} random bases, mixed denoms {sorted(set(DENS))}")
    check("P1.2 double dual == original lattice (lattice equality)",
          ok_dd, f"{N} bases, dual basis scrambled by GL(O) first")
    check("P1.3 x = (A*)^-1 c (c integral) pairs integrally with all "
          "basis columns", ok_mem_pos, f"{N} bases x 3 columns")
    check("P1.4 ... and with random lattice points", ok_lat_pt, f"{N} points")
    check("P1.5 x outside dual has a witness column: <x,a_j> = w_j "
          "non-integral", ok_mem_neg, f"{N} bases")


# ===========================================================================
# Part 2: inclusion reversal (Prop le:dual_equivalence bullet 1)
# ===========================================================================

def part2():
    print("\n--- Part 2: inclusion reversal under dualization ---")
    N = 60
    ok_rev = ok_strict = ok_idx = True
    n_strict = 0
    for t in range(N):
        A = rand_basis(3)
        T = rand_int_mat(3)
        B2 = A * rand_GLO(3)            # L2
        B1 = A * T * rand_GLO(3)        # L1 = L2 * T  (subset of L2)
        assert lattice_contains(B2, B1)
        D1, D2 = dual_basis(B1), dual_basis(B2)
        ok_rev &= lattice_contains(D1, D2)           # dual(L1) >= dual(L2)
        idx = lattice_index_log3(B2, B1)             # log_3 [L2:L1]
        ok_idx &= (lattice_index_log3(D1, D2) == idx)
        if idx > 0:                                  # strict inclusion
            n_strict += 1
            ok_strict &= not lattice_contains(D2, D1)
    check("P2.1 L1 < L2  =>  dual(L2) < dual(L1)", ok_rev, f"{N} pairs")
    check("P2.2 strict inclusion stays strict after dualization",
          ok_strict, f"{n_strict}/{N} pairs were strict")
    check("P2.3 index preserved: [L2:L1] == [dual L1 : dual L2]",
          ok_idx, f"{N} pairs")


# ===========================================================================
# Part 3: pi-equivalence preserved by dualization (bullet 2)
# ===========================================================================

def part3():
    print("\n--- Part 3: pi-equivalence preserved by dualization ---")
    N = 60
    ok_eq = ok_exp = True
    for t in range(N):
        A = rand_basis(3)
        k = rng.choice([-3, -2, -1, 1, 2, 3])
        B1 = A * rand_GLO(3)
        B2 = (A * chi_pow(k)) * rand_GLO(3)          # L2 = pi^k L1
        D1, D2 = dual_basis(B1), dual_basis(B2)
        # pi-equivalent, with the exact exponent dual(pi^k L) = pi^-k dual(L)
        hits = [j for j in range(-5, 6) if lattice_eq(D2, D1 * chi_pow(j))]
        ok_eq &= (len(hits) == 1)
        ok_exp &= (hits == [-k])
    check("P3.1 duals of pi-equivalent lattices are pi-equivalent "
          "(unique exponent in [-5,5])", ok_eq, f"{N} pairs")
    check("P3.2 exact exponent: dual(pi^k L) = pi^(-k) dual(L)",
          ok_exp, f"{N} pairs, k in -3..3 nonzero")


# ===========================================================================
# Part 4: determinant identity (proof of le:dual_scaling)
# ===========================================================================

def part4():
    print("\n--- Part 4: det(dual L) = conj(det L)^-1 as O_pi-modules ---")
    N = 100
    ok_exact = ok_mod = ok_scal = True
    for t in range(N):
        A = rand_basis(3)
        # exact identity for the distinguished basis (A*)^-1:
        ok_exact &= (dual_basis(A).det() * A.det().conj() == ONE)
        # module-level: any basis of the dual has det of valuation -v(det A),
        # and the unit ambiguity is exactly GL(O) determinants
        B = A * rand_GLO(3)
        D = dual_basis(B) * rand_GLO(3)
        ok_mod &= (my_vpi(D.det()) == -my_vpi(A.det()))
        ratio = D.det() * B.det().conj()
        ok_mod &= (my_vpi(ratio) == 0)               # unit, i.e. same module
        # det(pi^i L) = pi^(3i) det(L)  (rank 3; also used in the proof)
        i = rng.randint(-3, 3)
        ok_scal &= (my_vpi((B * chi_pow(i)).det()) == my_vpi(B.det()) + 3 * i)
    check("P4.1 exact: det((A*)^-1) * conj(det A) == 1", ok_exact,
          f"{N} bases")
    check("P4.2 module: v_pi(det dual) == -v_pi(det L); ratio to "
          "conj(det L)^-1 is a unit", ok_mod, f"{N} scrambled bases")
    check("P4.3 det(pi^i L) = pi^(3i) det(L) (rank 3)", ok_scal,
          f"{N} samples, i in -3..3")


# ===========================================================================
# Part 5: Lemma le:dual_scaling in rank 3
# ===========================================================================

def rand_unitary3():
    """Random product of the paper's gates (unitary over F)."""
    gates = [H_GATE, S_GATE, R_GATE]
    U = Mat.identity(3)
    for _ in range(rng.randint(3, 8)):
        U = U * rng.choice(gates)
    return U


def part5():
    print("\n--- Part 5: Lemma le:dual_scaling (rank 3) ---")
    cases = []           # (basis, source-tag)
    # (a) unitary * GL(O) -> self-dual lattice, then scale by chi^m
    for _ in range(25):
        U = rand_unitary3()
        cases.append((U * rand_GLO(3), "unitary*GLO"))
    # (b) tree machinery: pure vertices from random walks are self-dual
    tree_pure = []
    for _ in range(8):
        cur = ORIGIN
        for _ in range(rng.randint(2, 6)):
            cur = rng.choice(neighbors(cur))
        if cur.kind == "P":
            tree_pure.append(cur.g)
        else:
            nb = rng.choice(neighbors(cur))   # neighbors of A-vertex are pure
            tree_pure.append(nb.g)
    for g in tree_pure:
        cases.append((g, "tree-pure"))
    # (c) trivial self-dual O^3 scrambled
    for _ in range(10):
        cases.append((Mat.identity(3), "identity"))

    n_cases = 0
    ok_selfdual = ok_int = ok_eq = ok_even = ok_half = True
    seen_i = set()
    for g0, tag in cases:
        # confirm base lattice is self-dual *via the definition*
        ok_selfdual &= lattice_eq(dual_basis(g0), g0)
        ok_selfdual &= is_self_dual(g0)     # gram criterion agrees
        m = rng.randint(-3, 3)
        B = (g0 * chi_pow(m)) * rand_GLO(3)        # Lambda = pi^m L0
        n_cases += 1
        # i forced by determinants: 3i = -2 v(det B)
        k = my_vpi(B.det())
        ok_int &= (2 * k) % 3 == 0
        i = (-2 * k) // 3
        seen_i.add(i)
        ok_eq &= lattice_eq(dual_basis(B), B * chi_pow(i))
        ok_even &= (i % 2 == 0)
        Bs = (B * chi_pow(i // 2)) * rand_GLO(3)
        ok_half &= lattice_eq(dual_basis(Bs), Bs)
    check("P5.1 constructed base lattices are self-dual (definition == "
          "gram criterion)", ok_selfdual, f"{len(cases)} lattices")
    check("P5.2 dual(L) = pi^i L holds with i = -2 v_pi(det)/3",
          ok_int and ok_eq, f"{n_cases} scaled lattices, i values {sorted(seen_i)}")
    check("P5.3 the exponent i is always EVEN", ok_even,
          f"{n_cases} lattices")
    check("P5.4 pi^(i/2) L is self-dual", ok_half, f"{n_cases} lattices")

    # random scan: arbitrary lattices; any that are pi-equiv to dual must
    # have even i (det obstruction: 3 | 2k forces 3 | k, i = -2(k/3) even)
    n_scan, n_equiv, n_det_blocked = 300, 0, 0
    ok_scan = True
    for _ in range(n_scan):
        A = rand_basis(3)
        k = my_vpi(A.det())
        if (2 * k) % 3 != 0:
            n_det_blocked += 1
            # det obstruction: cannot be pi-equivalent to its dual at all
            ok_scan &= not any(lattice_eq(dual_basis(A), A * chi_pow(j))
                               for j in range(-4, 5))
            continue
        i = (-2 * k) // 3
        if lattice_eq(dual_basis(A), A * chi_pow(i)):
            n_equiv += 1
            ok_scan &= (i % 2 == 0)
    check("P5.5 random scan: every lattice pi-equiv to its dual has even i; "
          "det-blocked ones are never equivalent", ok_scan,
          f"{n_scan} random lattices: {n_equiv} equivalent, "
          f"{n_det_blocked} det-blocked (3 not| v_pi(det))")

    # alternating tree vertices: dual strictly inside, index 9, and NOT
    # pi-equivalent to the dual (v(det) = -1 => 3 not| 2k)
    ok_alt = True
    n_alt = 0
    for _ in range(6):
        cur = rng.choice(neighbors(ORIGIN))     # alternating vertex
        for _ in range(rng.randint(0, 3)):
            nxt = rng.choice(neighbors(cur))
            cur = rng.choice(neighbors(nxt))    # stay on alternating side
        if cur.kind != "A":
            continue
        n_alt += 1
        g = cur.g
        D = dual_basis(g)
        ok_alt &= lattice_contains(g, D) and not lattice_eq(g, D)
        ok_alt &= (lattice_index_log3(g, D) == 2)          # index 9
        ok_alt &= not any(lattice_eq(D, g * chi_pow(j)) for j in range(-4, 5))
    check("P5.6 alternating vertices: dual(L) < L of index 9, never "
          "pi-equivalent to L", ok_alt, f"{n_alt} alternating vertices")


# ===========================================================================
# Part 6: the Remark -- even dimensions
# ===========================================================================

def cartan2(g: Mat):
    """Own Cartan/Smith decomposition over O_pi for 2x2:
    returns (k, (a1,a2), kp) with g = k diag(chi^a1, chi^a2) kp,
    k, kp in GL_2(O_pi), a1 <= a2.  Correctness is verified by the caller."""
    E = Mat([[ZERO, ONE], [ONE, ZERO]])
    W = Mat([row[:] for row in g.m])
    P, Q = Mat.identity(2), Mat.identity(2)      # W = P g Q
    vals = [(W.m[i][j].v_pi(), i, j) for i in range(2) for j in range(2)]
    v0, i0, j0 = min(vals, key=lambda t: (t[0], t[1], t[2]))
    if i0 == 1:
        W, P = E * W, E * P
    if j0 == 1:
        W, Q = W * E, Q * E
    a1 = W.m[0][0].v_pi()
    u = chi_pow(a1) / W.m[0][0]                  # unit
    L = Mat([[u, ZERO], [ZERO, ONE]])
    W, P = L * W, L * P
    c = W.m[1][0] * chi_pow(-a1)                 # in O_pi
    L = Mat([[ONE, ZERO], [-c, ONE]])
    W, P = L * W, L * P
    d = W.m[0][1] * chi_pow(-a1)                 # in O_pi
    R = Mat([[ONE, -d], [ZERO, ONE]])
    W, Q = W * R, Q * R
    a2 = W.m[1][1].v_pi()
    u2 = chi_pow(a2) / W.m[1][1]
    R = Mat([[ONE, ZERO], [ZERO, u2]])
    W, Q = W * R, Q * R
    return P.inv(), (a1, a2), Q.inv()


def part6():
    print("\n--- Part 6: the Remark (even dimensions) ---")

    # ---- 6a: -1 is not a norm of F_pi / Q_3 --------------------------------
    counts = []
    for k in (1, 2, 3):
        mod = 3 ** k
        sols = sum(1 for a in range(mod) for b in range(mod)
                   if (a * a - a * b + b * b) % mod == (mod - 1) % mod)
        counts.append(sols)
    check("P6a.1 x*conj(x) = -1 has NO solution mod 3, 9, 27 "
          "(so -1 is not a norm; unit norms need residue 1)",
          counts == [0, 0, 0], f"solution counts mod 3/9/27 = {counts}")

    unit_norm_res = sorted({(a * a - a * b + b * b) % 3
                            for a in range(3) for b in range(3)
                            if (a + b) % 3 != 0})
    check("P6a.2 norms of units of O_pi are == 1 mod 3 (exhaustive mod 3)",
          unit_norm_res == [1], f"residues found: {unit_norm_res}")

    n_minus2 = sum(1 for a in range(27) for b in range(27)
                   if (a * a - a * b + b * b) % 27 == 25)
    check("P6a.3 (info) x*conj(x) = -2 IS solvable mod 27 (-2 == 1 mod 3, "
          "consistent with unit-norm criterion)", n_minus2 > 0,
          f"{n_minus2} solutions mod 27")

    # ---- 6b: standard 2-dim form is anisotropic over F_pi ------------------
    # primitive (a,b) in (O/3)^2 means at least one coordinate is a unit.
    bad = []
    for a1 in range(3):
        for a2 in range(3):
            for b1 in range(3):
                for b2 in range(3):
                    unit_a = (a1 + a2) % 3 != 0
                    unit_b = (b1 + b2) % 3 != 0
                    if not (unit_a or unit_b):
                        continue
                    Na = (a1 * a1 - a1 * a2 + a2 * a2) % 3
                    Nb = (b1 * b1 - b1 * b2 + b2 * b2) % 3
                    if (Na + Nb) % 3 == 0:
                        bad.append((a1, a2, b1, b2))
    check("P6b.1 Q(x) = N(x1)+N(x2) is a UNIT for every primitive x "
          "(exhaustive 81-pair scan mod 3) => standard 2-dim form "
          "anisotropic over F_pi", not bad,
          f"{len(bad)} primitive isotropic residue pairs (expected 0)")

    # ---- 6c: dim 2, standard form: odd i impossible -------------------------
    # (i) own Cartan decomposition verified exactly on random 2x2 bases
    ok_cartan = True
    n_cart = 200
    for _ in range(n_cart):
        g = rand_basis(2)
        k, (a1, a2), kp = cartan2(g)
        ok_cartan &= k.in_GL_O() and kp.in_GL_O() and a1 <= a2
        ok_cartan &= (k * Mat.diag(chi_pow(a1), chi_pow(a2)) * kp == g)
        ok_cartan &= lattice_eq(g, k * Mat.diag(chi_pow(a1), chi_pow(a2)))
    check("P6c.1 every rank-2 lattice = k diag(chi^a1, chi^a2) O^2, "
          "k in GL_2(O) (own Cartan, exact)", ok_cartan,
          f"{n_cart} random bases")

    # (ii) for ALL k in GL_2(O) and odd s, the candidate dual(L) = pi^{-s} L
    # fails at a diagonal Gram entry, because Q(column of k) is a unit (6b).
    ok_diag = ok_fail = ok_consist = True
    n_k = 500
    for _ in range(n_k):
        k = rand_GLO(2)
        M = gram(k)                                  # k* k in GL_2(O)
        ok_diag &= (my_vpi(M.m[0][0]) == 0 and my_vpi(M.m[1][1]) == 0)
        for s in (1, -1, 3, -3):
            B = k * Mat.diag(chi_pow(s), ONE)        # L = k diag(chi^s,1) O^2
            G = gram(B)
            cond_gram = (G * chi_pow(-s)).in_GL_O()  # dual = pi^{-s} L ?
            cond_lat = lattice_eq(dual_basis(B), B * chi_pow(-s))
            ok_consist &= (cond_gram == cond_lat)
            ok_fail &= not cond_lat
            # the failing entry is the predicted diagonal one
            bad_entry = (G.m[1][1] * chi_pow(-s) if s > 0
                         else G.m[0][0] * chi_pow(-s))
            ok_fail &= my_vpi(bad_entry) < 0
    check("P6c.2 Gram diagonal of k*k is unit for k in GL_2(O) "
          "(= anisotropy mod pi)", ok_diag, f"{n_k} samples")
    check("P6c.3 dual(L) = pi^i L with i ODD fails for ALL "
          "L = k diag(chi^s,1) O^2, s in {+-1,+-3} -- failing witness is "
          "the predicted diagonal entry", ok_fail,
          f"{n_k} x 4 candidates, all fail")
    check("P6c.4 gram criterion chi^i G in GL_2(O) <=> lattice_eq dual test",
          ok_consist, f"{4 * n_k} comparisons")

    # combine: random search over generic 2x2 lattices
    n_search, n_equiv2, odd_found = 4000, 0, 0
    ok_even2 = True
    for _ in range(n_search):
        A = rand_basis(2)
        i = -my_vpi(A.det())            # n=2: det forces i = -v(det)
        if lattice_eq(dual_basis(A), A * chi_pow(i)):
            n_equiv2 += 1
            if i % 2:
                odd_found += 1
                ok_even2 = False
    # seeded even examples so the search is non-vacuous
    n_seed = 0
    for _ in range(30):
        m = rng.randint(-2, 2)
        B = (Mat.identity(2) * chi_pow(m)) * rand_GLO(2)   # pi^m O^2
        i = -my_vpi(B.det())
        if lattice_eq(dual_basis(B), B * chi_pow(i)):
            n_seed += 1
            ok_even2 &= (i % 2 == 0) and (i == -2 * m)
    check("P6c.5 dim-2 random + seeded search: NO odd-exponent self-dual-"
          "up-to-pi lattice for the standard form", ok_even2,
          f"{n_search} random (hits: {n_equiv2}, odd: {odd_found}); "
          f"{n_seed}/30 seeded even examples verified")

    # (iii) skew-unimodular residue lemma: N = chi*H in GL_2(O), H Hermitian
    # => N* = -w^2 N, N == antidiag(b,-b) mod pi, det N == 1 mod 3.
    # The Gram side would force det == 2 mod 3: incompatible.
    MW2 = -(OMEGA * OMEGA)
    n_skew = 0
    ok_skew = True
    while n_skew < 30:
        h11 = Fw(rng.randint(-6, 6), rng.choice([1, 2, 5]))
        h22 = Fw(rng.randint(-6, 6), rng.choice([1, 2, 5]))
        z = Fw(Zw(rng.randint(-4, 4), rng.randint(-4, 4)),
               rng.choice([1, 2, 5]))
        if z.is_zero() or my_vpi(z) != 0:
            continue
        h12 = z * CHI_INV                       # v_pi = -1 exactly
        H = Mat([[h11, h12], [h12.conj(), h22]])
        N = H * CHI
        if not N.in_GL_O():
            continue
        n_skew += 1
        ok_skew &= (N.conjT() == N * MW2)       # N* = -w^2 N
        r11, r22 = N.m[0][0].res3(), N.m[1][1].res3()
        r12, r21 = N.m[0][1].res3(), N.m[1][0].res3()
        ok_skew &= (r11 == 0 and r22 == 0)
        ok_skew &= (r21 == (-r12) % 3) and (r12 != 0)
        ok_skew &= (N.det().res3() == 1)
    check("P6c.6 skew-unimodular N (N* = -w^2 N) has zero diagonal mod pi "
          "and det == 1 mod 3; norm/Gram side needs det == 2 mod 3 "
          "=> obstruction", ok_skew, f"{n_skew} random skew-unimodular N")

    # ---- 6d: dim 2 with the HYPERBOLIC form: odd example exists -------------
    J2 = Mat([[ZERO, ONE], [ONE, ZERO]])

    def dual_basis_J(A, J):
        # dual_J = {x : A* J^T x in O^n}; J real symmetric here
        return (A.conjT() * J).inv()

    A2 = Mat.diag(ONE, CHI_INV)
    DJ = dual_basis_J(A2, J2)
    okJ = lattice_eq(DJ, A2 * CHI)              # dual_J(L) = pi L, i = 1 odd
    # defining-property test for the J-dual
    ok_def = True
    for _ in range(40):
        x = matvec(DJ, rand_int_vec(2))
        y = matvec(A2, rand_int_vec(2))
        val = sum((x[i] * J2.m[i][j].num.a * y[j].conj()
                   for i in range(2) for j in range(2)), ZERO)
        ok_def &= (val.is_zero() or my_integral(val))
    w = [ONE, CHI_INV]                          # non-integral coordinate 2
    xb = matvec(DJ, w)
    valb = sum((xb[i] * J2.m[i][j].num.a * A2.col(1)[j].conj()
                for i in range(2) for j in range(2)), ZERO)
    ok_neg = (my_vpi(valb) < 0)
    check("P6d.1 dim 2, hyperbolic form J: L = diag(1, chi^-1) O^2 has "
          "dual_J(L) = pi L  (i = 1 ODD)", okJ and ok_def and ok_neg,
          "exact; J-dual defining property checked on 40 points + witness")

    # ---- 6e: dim 4, standard form: explicit ODD example --------------------
    v1 = [1, 1, 1, 0]
    v2 = [0, 1, -1, 1]
    # totally isotropic plane mod 3 for the residue form sum x_i y_i:
    q1 = sum(x * x for x in v1) % 3
    q2 = sum(x * x for x in v2) % 3
    b12 = sum(x * y for x, y in zip(v1, v2)) % 3
    check("P6e.1 span{(1,1,1,0),(0,1,-1,1)} is a totally isotropic plane "
          "of the residue form on F_3^4", q1 == 0 and q2 == 0 and b12 == 0,
          f"Q(v1)={q1}, Q(v2)={q2}, B(v1,v2)={b12}")

    cols = [[Fw(x) * CHI_INV for x in v1],
            [Fw(x) * CHI_INV for x in v2],
            [ONE, ZERO, ZERO, ZERO],
            [ZERO, ONE, ZERO, ZERO]]
    A4 = Mat([[cols[j][i] for j in range(4)] for i in range(4)])
    # unimodularity of the underlying integer matrix [v1 v2 e1 e2]
    d_int = det_int([[v1[i], v2[i], 1 if i == 0 else 0, 1 if i == 1 else 0]
                     for i in range(4)])
    check("P6e.2 [v1 v2 e1 e2] is unimodular over Z (so A4 is a basis and "
          "L = O^4 + chi^-1 v1 O + chi^-1 v2 O)", d_int in (1, -1),
          f"integer det = {d_int}; v_pi(det A4) = {my_vpi(A4.det())} "
          "(must be -2)")

    G4 = gram(A4)
    N4 = G4 * CHI
    ok_int4 = N4.is_O() and N4.det().v_pi() == 0
    ok_dual4 = lattice_eq(dual_basis(A4), A4 * CHI)
    check("P6e.3 chi * Gram(A4) in GL_4(O)  <=>  dual(L) = pi L: "
          "ODD exponent i = 1 in dimension 4, standard form",
          ok_int4 and ok_dual4,
          f"det(chi G) v_pi = {N4.det().v_pi()}, lattice_eq confirms")

    # residue structure: chi*G antisymmetric, zero diagonal, full rank mod 3
    R4 = [[N4.m[i][j].res3() for j in range(4)] for i in range(4)]
    ok_res = all(R4[i][i] == 0 for i in range(4))
    ok_res &= all(R4[j][i] == (-R4[i][j]) % 3 for i in range(4)
                  for j in range(4))
    ok_res &= det_int(R4) % 3 != 0
    check("P6e.4 chi*Gram mod pi is antisymmetric, zero diagonal, "
          "nondegenerate (isotropic-plane structure)", ok_res,
          f"residue matrix {R4}, det mod 3 = {det_int(R4) % 3}")

    # defining property of the dim-4 dual (dual_basis correctness in dim 4)
    D4 = dual_basis(A4)
    ok_mem4 = True
    for _ in range(40):
        x = matvec(D4, rand_int_vec(4))
        for j in range(4):
            val = inner(x, A4.col(j))
            ok_mem4 &= (val.is_zero() or my_integral(val))
        y = matvec(A4, rand_int_vec(4))
        val = inner(x, y)
        ok_mem4 &= (val.is_zero() or my_integral(val))
    w = rand_int_vec(4)
    w[2] = w[2] + CHI_INV
    xb = matvec(D4, w)
    ok_mem4 &= my_vpi(inner(xb, A4.col(2))) < 0
    check("P6e.5 dim-4 dual_basis satisfies the defining property "
          "(40 points + non-membership witness)", ok_mem4, "exact")

    # no pi-rescaling of L is self-dual (lemma's conclusion genuinely fails)
    ok_nosd = True
    for m in range(-3, 4):
        Bm = A4 * chi_pow(m)
        ok_nosd &= not lattice_eq(dual_basis(Bm), Bm)
        # exponent pattern: dual(pi^m L) = pi^(1-2m) L  -- always odd
        ok_nosd &= lattice_eq(dual_basis(Bm), Bm * chi_pow(1 - 2 * m))
    check("P6e.6 no pi-rescaling pi^m L (m in -3..3) is self-dual; "
          "dual(pi^m L) = pi^(1-2m) pi^m L stays odd", ok_nosd, "7 scalings")

    print("\n  VERDICT on the Remark (line 320-321): 'fails for even "
          "dimensions' is")
    print("   * TRUE  in dim 4 with the paper's standard form "
          "(explicit lattice with dual(L) = pi L, i = 1 odd);")
    print("   * FALSE in dim 2 with the paper's standard form: the form is "
          "anisotropic over F_pi,")
    print("     and dual(L) = pi^i L forces i EVEN (verified impossibility "
          "of odd i);")
    print("   * in dim 2 odd exponents require a different (isotropic, "
          "e.g. hyperbolic) Hermitian form.")


# ===========================================================================

def main():
    print(f"v02_duality.py  (seed {SEED})")
    part0()
    part1()
    part2()
    part3()
    part4()
    part5()
    part6()

    n_fail = sum(1 for _, ok, _ in RESULTS if not ok)
    print(f"\n=== SUMMARY: {len(RESULTS) - n_fail}/{len(RESULTS)} checks "
          f"passed ===")
    if n_fail:
        for name, ok, det in RESULTS:
            if not ok:
                print(f"  FAILED: {name}  [{det}]")
        sys.exit(1)
    sys.exit(0)


if __name__ == "__main__":
    main()
