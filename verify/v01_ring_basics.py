#!/usr/bin/env python
"""
v01_ring_basics.py -- Verify arithmetic foundations of "Buildings for Synthesis
with Clifford+R", paper Section 2 (lines ~136-230) + Lemma le:padic_expansion
(line 214) + Lemma le:techincal_lemma (line 681).

Claims verified:
  1. 3 O_F = chi^2 O_F :  3 = -omega^2 chi^2 exactly, chi*conj(chi) = 3,
     conj(chi) = -omega^2 chi (unit multiple).
  2. Residue field O_F/chi = F_3, omega == 1 mod chi; residue map
     a+b*omega -> (a+b) mod 3 is a ring homomorphism with kernel chi*O_F.
  3. Z[1/3, omega] = Z[chi^{-1}] :  1/chi = (2+omega)/3  and  1/3 = -omega*chi^{-2}.
  4. Lemma le:padic_expansion: greedy digit expansion in {0,1,2} with unique
     digits; every x in F is chi^i * (integral) with i <= 0.
  5. Lemma le:techincal_lemma (contrapositive, constructive):
     min(v(x),v(y)) = n < 0  ==>  v(conj(x)x + conj(y)y) = 2n < 0,
     in both the v(x) < v(y) case and the v(x) == v(y) case (unit residues
     are +-1, so unit parts contribute 1 + 1 = 2 != 0 mod 3).
     Plus: v_pi restricted to Q is always even (= 2*v_3).
  6. v_pi properties: v(xy) = v(x)+v(y); v(x+y) >= min, with equality when
     v(x) != v(y); v_pi(3) = 2 so |3|_pi = 3^{-2/2} = 1/3.
  7. Gates: all 9 entries of H have v_pi = -1, hence l(H) = 2; H, S, R are
     exactly unitary; H^2 = -P23; det(H) = -1.

btlib.py is treated as a TOOL, not ground truth: its v_pi and res3 are
cross-checked against an INDEPENDENT implementation (valuation via the
3-adic valuation of the field norm N(a+b w) = a^2 - a b + b^2; residue via
"the unique r in {0,1,2} with v_pi(x - r) >= 1"), and Zw/Fw arithmetic is
cross-checked against numpy complex floats.

Deterministic (fixed seed). Prints PASS/FAIL per claim; exits nonzero on
any failure.
"""

import sys, os, random
from fractions import Fraction
from math import inf

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np

from btlib import (Zw, Fw, Mat, ZERO, ONE, OMEGA, CHI, CHI_INV, I_SQRT3,
                   chi_pow, H_GATE, S_GATE, R_GATE, ell)

random.seed(20260610)

FAILURES = []


def report(name, ok, detail=""):
    tag = "PASS" if ok else "FAIL"
    print(f"[{tag}] {name}" + (f" -- {detail}" if detail else ""))
    if not ok:
        FAILURES.append(name)


# ---------------------------------------------------------------------------
# Independent primitives (NOT from btlib)
# ---------------------------------------------------------------------------

W_C = np.exp(2j * np.pi / 3)          # float omega


def to_c(x):
    """Fw or Zw -> complex float (independent embedding)."""
    if isinstance(x, Zw):
        return x.a + x.b * W_C
    return (x.num.a + x.num.b * W_C) / x.den


def v3_int(n):
    """3-adic valuation of a nonzero integer."""
    n = abs(int(n))
    assert n != 0
    v = 0
    while n % 3 == 0:
        n //= 3
        v += 1
    return v


def v_pi_indep(x):
    """Independent v_pi: for x = (a+b w)/d,
       v_pi = v_3(N(a+b w)) - 2*v_3(d),  N(a+b w) = a^2 - a b + b^2.
    Justification: v_pi(z) = v_3(N(z)) for z in Z[omega] because the norm is
    multiplicative, N(chi) = 3, primes other than chi have norm coprime to 3,
    and units have norm 1; and v_pi(3) = 2 (ramification)."""
    if isinstance(x, Zw):
        x = Fw(x, 1)
    if x.num.is_zero():
        return inf
    a, b = x.num.a, x.num.b
    return v3_int(a * a - a * b + b * b) - 2 * v3_int(x.den)


def res_indep(x):
    """Independent residue in O_pi/chi = F_3 of an integral x:
    the unique r in {0,1,2} with v_pi(x - r) >= 1."""
    assert v_pi_indep(x) >= 0
    hits = [r for r in range(3) if v_pi_indep(x - Fw(r)) >= 1]
    assert len(hits) == 1, f"residue not unique/existing: {x!r} -> {hits}"
    return hits[0]


CHI_Z = Zw(1, -1)                      # chi = 1 - omega as a Zw
CHI_CONJ_Z = Zw(2, 1)                  # conj(chi) = 2 + omega


def chi_pow_indep(k):
    """chi^k as Fw, built without btlib's chi_pow:
    chi^{-1} = conj(chi)/3 = (2+omega)/3."""
    if k >= 0:
        z = Zw(1, 0)
        for _ in range(k):
            z = z * CHI_Z
        return Fw(z, 1)
    z = Zw(1, 0)
    for _ in range(-k):
        z = z * CHI_CONJ_Z
    return Fw(z, 3 ** (-k))


# -- random generators --------------------------------------------------------

def rand_Zw(lo=-100, hi=100, nonzero=False):
    while True:
        z = Zw(random.randint(lo, hi), random.randint(lo, hi))
        if not (nonzero and z.is_zero()):
            return z


def rand_Fw(nonzero=False):
    """Random element of F with denominator 3^m * (coprime part)."""
    z = rand_Zw(nonzero=nonzero)
    den = (3 ** random.randint(0, 3)) * random.choice([1, 1, 2, 5, 7])
    return Fw(z, den)


def rand_unit():
    """Random unit of O_pi (v_pi = 0), with an occasional 3-coprime
    denominator (which is a unit of O_pi)."""
    while True:
        z = rand_Zw(nonzero=True)
        if (z.a + z.b) % 3 != 0:
            d = random.choice([1, 1, 1, 2, 5, 7])
            return Fw(z, d)


# ---------------------------------------------------------------------------
# Sanity check 0: btlib's arithmetic & primitives vs independent ones
# ---------------------------------------------------------------------------

def check0_sanity():
    n = 400
    ok_arith = True
    for _ in range(n):
        x, y = rand_Fw(), rand_Fw(nonzero=True)
        cx, cy = to_c(x), to_c(y)
        for got, want in [(x + y, cx + cy), (x - y, cx - cy),
                          (x * y, cx * cy), (x / y, cx / cy),
                          (x.conj(), np.conj(cx)), (-x, -cx)]:
            if abs(to_c(got) - want) > 1e-7 * (1 + abs(want)):
                ok_arith = False
    report("sanity: Fw +,-,*,/,conj,neg vs numpy complex",
           ok_arith, f"{n} random pairs, tol 1e-7 relative")

    # norm formula N(a+bw) = a^2-ab+b^2 = |a+bw|^2 (floats)
    ok_norm = all(abs(rz.norm() - abs(to_c(rz)) ** 2) < 1e-6
                  for rz in (rand_Zw(nonzero=True) for _ in range(n)))
    report("sanity: Zw.norm == |.|^2 (float)", ok_norm, f"{n} samples")

    # btlib v_pi vs independent norm-based v_pi
    bad_v = 0
    for _ in range(n):
        x = rand_Fw(nonzero=True)
        if x.v_pi() != v_pi_indep(x):
            bad_v += 1
    report("sanity: btlib Fw.v_pi == independent norm-based v_pi",
           bad_v == 0, f"{n} samples, {bad_v} mismatches")

    # btlib res3 vs independent residue (on integral elements)
    bad_r = 0
    cnt = 0
    for _ in range(n):
        x = rand_Fw(nonzero=True)
        v = v_pi_indep(x)
        if v < 0:
            x = x * chi_pow_indep(-v)     # make integral
        cnt += 1
        if x.res3() != res_indep(x):
            bad_r += 1
    report("sanity: btlib Fw.res3 == independent residue",
           bad_r == 0, f"{cnt} integral samples, {bad_r} mismatches")

    # Zw multiplication closed-form vs floats (used by residue-hom check)
    ok_m = True
    for _ in range(n):
        z, w = rand_Zw(), rand_Zw()
        a, b, c, d = z.a, z.b, w.a, w.b
        prod = Zw(a * c - b * d, a * d + b * c - b * d)   # independent formula
        if prod != z * w or abs(to_c(prod) - to_c(z) * to_c(w)) > 1e-5:
            ok_m = False
    report("sanity: Zw multiplication formula vs floats", ok_m, f"{n} samples")


# ---------------------------------------------------------------------------
# 1. 3 O_F = chi^2 O_F
# ---------------------------------------------------------------------------

def check1_chi_squared():
    w = Zw(0, 1)
    chi = Zw(1, 0) - w                              # 1 - omega
    w2 = w * w
    ok_a = (-(w2 * (chi * chi)) == Zw(3, 0))
    report("1a: 3 == -omega^2 * chi^2 (exact in Z[omega])", ok_a,
           f"-w^2*chi^2 = {-(w2 * (chi * chi))!r}")

    ok_b = (chi * chi.conj() == Zw(3, 0))
    report("1b: chi * conj(chi) == 3", ok_b, f"chi*conj(chi) = {chi * chi.conj()!r}")

    # conj(chi) = -omega^2 * chi  (unit multiple => same ideal)
    ok_c = (chi.conj() == -(w2 * chi))
    report("1c: conj(chi) == -omega^2 * chi (unit multiple)", ok_c)

    # chi^2 = -3*omega, so 3 | chi^2 and chi^2 | 3 in Z[omega] (equal ideals)
    ok_d = (chi * chi == Zw(0, -3))
    # 3 / chi^2 = -omega^2 in Z[omega]; chi^2 / 3 = -omega in Z[omega]... wait
    q1 = Fw(3) / Fw(chi * chi)                      # should be -omega^2 = 1+omega
    ok_e = (q1 == Fw(Zw(1, 1)))                     # -w^2 = 1 + w
    report("1d: chi^2 == -3*omega and 3/chi^2 == -omega^2 in Z[omega] "
           "(=> 3 O_F = chi^2 O_F)", ok_d and ok_e,
           f"chi^2 = {chi * chi!r}, 3/chi^2 = {q1!r}")

    ok_f = abs(abs(to_c(Fw(chi))) ** 2 - 3) < 1e-12
    report("1e: |chi|^2 == 3 numerically", ok_f)


# ---------------------------------------------------------------------------
# 2. Residue field O_F/chi = F_3, omega == 1 mod chi, residue map is ring hom
# ---------------------------------------------------------------------------

def res_map(z: Zw):
    return (z.a + z.b) % 3


def chi_divides(z: Zw):
    """Independent: chi | z iff z*conj(chi)/3 = ((2a-b)+(a+b)w)/3 in Z[omega]."""
    return (2 * z.a - z.b) % 3 == 0 and (z.a + z.b) % 3 == 0


def check2_residue():
    n = 1000
    ok_add = ok_mul = True
    for _ in range(n):
        x, y = rand_Zw(), rand_Zw()
        if res_map(x + y) != (res_map(x) + res_map(y)) % 3:
            ok_add = False
        if res_map(x * y) != (res_map(x) * res_map(y)) % 3:
            ok_mul = False
    ok_one = (res_map(Zw(1, 0)) == 1)
    report("2a: residue map a+b*omega -> (a+b) mod 3 is a ring hom",
           ok_add and ok_mul and ok_one,
           f"{n} random pairs: additive ok={ok_add}, multiplicative ok={ok_mul}, f(1)=1 ok={ok_one}")

    # kernel == chi O_F  (both directions, random sample)
    bad_ker = 0
    for _ in range(n):
        z = rand_Zw()
        if (res_map(z) == 0) != chi_divides(z):
            bad_ker += 1
    report("2b: ker(residue map) == chi*O_F", bad_ker == 0,
           f"{n} samples, {bad_ker} mismatches (res==0 iff chi|z)")

    # omega == 1 mod chi:  omega - 1 = -chi
    w = Zw(0, 1)
    ok_w = (res_map(w) == 1) and (w - Zw(1, 0) == -CHI_Z) and chi_divides(w - Zw(1, 0))
    report("2c: omega == 1 mod chi (omega - 1 = -chi)", ok_w)

    # image is all of F_3 -> residue field has exactly 3 elements
    ok_img = sorted({res_map(Zw(r, 0)) for r in range(3)}) == [0, 1, 2]
    report("2d: residue map surjects onto F_3 = {0,1,2}", ok_img)


# ---------------------------------------------------------------------------
# 3. Z[1/3, omega] = Z[chi^{-1}]
# ---------------------------------------------------------------------------

def check3_ring_equality():
    # 1/chi = (2+omega)/3  (element of Z[omega, 1/3])
    inv_chi = Fw(Zw(2, 1), 3)
    ok_a = (Fw(CHI_Z) * inv_chi == ONE) and (CHI_INV == inv_chi)
    report("3a: 1/chi == (2+omega)/3 in Z[omega,1/3]", ok_a,
           f"chi * (2+omega)/3 = {Fw(CHI_Z) * inv_chi!r}")

    # 1/3 = -omega * chi^{-2}  (element of Z[omega, 1/chi])
    lhs = -(OMEGA * chi_pow_indep(-2))
    ok_b = (lhs == Fw(1, 3))
    # equivalent integral statement: chi^2 == -3*omega
    ok_c = (CHI_Z * CHI_Z == Zw(0, -3))
    report("3b: 1/3 == -omega * chi^{-2} in Z[omega,1/chi]", ok_b and ok_c,
           f"-omega*chi^-2 = {lhs!r}; chi^2 = {CHI_Z * CHI_Z!r} = -3w")

    ok_f = abs(to_c(inv_chi) - 1 / to_c(Fw(CHI_Z))) < 1e-12
    report("3c: 1/chi == (2+omega)/3 numerically", ok_f)


# ---------------------------------------------------------------------------
# 4. Lemma le:padic_expansion
# ---------------------------------------------------------------------------

def check4_padic_expansion():
    n, K = 200, 12
    bad = 0
    bad_unique = 0
    checked = 0
    for _ in range(n):
        # random x in O_pi as z/3^m
        z = rand_Zw(nonzero=True)
        m = random.randint(0, 3)
        x = Fw(z, 3 ** m)
        v = v_pi_indep(x)
        if v < 0:
            x = x * chi_pow_indep(-v)        # still of the form z'/3^m
        assert v_pi_indep(x) >= 0
        checked += 1
        # greedy digits: at step i, unique d in {0,1,2} with v(y - d chi^i) >= i+1
        y = x
        digits = []
        ok = True
        for i in range(K):
            if v_pi_indep(y) < i:
                ok = False
                break
            hits = [d for d in range(3)
                    if v_pi_indep(y - Fw(d) * chi_pow_indep(i)) >= i + 1]
            if len(hits) != 1:
                bad_unique += 1
                ok = False
                break
            d = hits[0]
            digits.append(d)
            y = y - Fw(d) * chi_pow_indep(i)
        if not ok:
            bad += 1
            continue
        # verify the tail: x - sum_{i<K} d_i chi^i has v_pi >= K
        tail = x
        for i, d in enumerate(digits):
            tail = tail - Fw(d) * chi_pow_indep(i)
        if not (v_pi_indep(tail) >= K and all(d in (0, 1, 2) for d in digits)):
            bad += 1
    report("4a: Lemma le:padic_expansion digit expansion "
           f"(digits in {{0,1,2}}, v_pi(x - sum_(i<{K}) d_i chi^i) >= {K})",
           bad == 0 and bad_unique == 0,
           f"{checked} random x in O_pi (z/3^m form), K={K}; "
           f"{bad} failures, {bad_unique} digit non-uniqueness events")

    # second part: every x in F_pi is chi^i * (integral), i in Z_{<=0}
    n2, bad2 = 200, 0
    for _ in range(n2):
        x = rand_Fw(nonzero=True)
        i = min(0, v_pi_indep(x))
        if not (i <= 0 and v_pi_indep(x * chi_pow_indep(-i)) >= 0):
            bad2 += 1
    report("4b: every x in F is chi^i * (element of O_pi) with i <= 0",
           bad2 == 0, f"{n2} random x in F, {bad2} failures")


# ---------------------------------------------------------------------------
# 5. Lemma le:techincal_lemma
# ---------------------------------------------------------------------------

def check5_technical_lemma():
    # conj preserves v_pi (used implicitly: v(conj(x)x) = 2 v(x))
    n0, bad0 = 300, 0
    for _ in range(n0):
        x = rand_Fw(nonzero=True)
        if v_pi_indep(x.conj()) != v_pi_indep(x):
            bad0 += 1
    report("5a: v_pi(conj(x)) == v_pi(x)", bad0 == 0,
           f"{n0} samples, {bad0} failures")

    # units of O_pi have residue +-1 (i.e. in {1,2}) and res(conj(u)*u) == 1
    n1, bad1 = 300, 0
    for _ in range(n1):
        u = rand_unit()
        r = res_indep(u)
        if r not in (1, 2) or res_indep(u.conj() * u) != 1:
            bad1 += 1
    report("5b: units u of O_pi have residue +-1 and res(conj(u)u) == 1",
           bad1 == 0, f"{n1} random units, {bad1} failures")

    # Case v(x) == v(y) == n < 0 (the interesting case):
    # v(conj(x)x + conj(y)y) == 2n exactly, because the unit parts give
    # residue 1 + 1 = 2 != 0 mod 3.
    nB, badB, badB_res = 250, 0, 0
    for _ in range(nB):
        n = random.randint(-6, -1)
        u1, u2 = rand_unit(), rand_unit()
        x = chi_pow_indep(n) * u1
        y = chi_pow_indep(n) * u2
        s = x.conj() * x + y.conj() * y
        if v_pi_indep(s) != 2 * n:
            badB += 1
        # the unit-part sum must have residue 2
        if res_indep(u1.conj() * u1 + u2.conj() * u2) != 2:
            badB_res += 1
    report("5c: v(x)==v(y)==n<0  =>  v(conj(x)x+conj(y)y) == 2n < 0",
           badB == 0 and badB_res == 0,
           f"{nB} pairs (n in [-6,-1]); {badB} valuation failures, "
           f"{badB_res} unit-residue(!=2) failures")

    # Case v(x) < v(y), n = v(x) < 0:
    nA, badA = 250, 0
    for _ in range(nA):
        n = random.randint(-6, -1)
        m = random.randint(n + 1, 4)
        x = chi_pow_indep(n) * rand_unit()
        y = chi_pow_indep(m) * rand_unit()
        if random.random() < 0.5:
            x, y = y, x          # also exercise the symmetric branch
        s = x.conj() * x + y.conj() * y
        if v_pi_indep(s) != 2 * n:
            badA += 1
    report("5d: min(v(x),v(y))=n<0, v(x)!=v(y)  =>  v(conj(x)x+conj(y)y) == 2n",
           badA == 0, f"{nA} pairs, {badA} failures")

    # Direct check of the lemma's implication on mixed random pairs
    nC, badC, hitC = 400, 0, 0
    for _ in range(nC):
        x, y = rand_Fw(nonzero=True), rand_Fw(nonzero=True)
        s = x.conj() * x + y.conj() * y
        if (s.is_zero() or v_pi_indep(s) >= 0):
            hitC += 1
            if not (v_pi_indep(x) >= 0 and v_pi_indep(y) >= 0):
                badC += 1
    report("5e: conj(x)x + conj(y)y in O_pi  =>  x, y in O_pi (direct)",
           badC == 0,
           f"{nC} mixed random pairs, {hitC} had integral sum, {badC} violations")

    # v_pi on Q (subset of Q_3) is always even, = 2 * v_3
    nQ, badQ = 100, 0
    for _ in range(nQ):
        p = random.randint(-200, 200)
        q = random.randint(1, 200)
        if p == 0:
            p = 9
        x = Fw(Zw(p, 0), q)
        v = v_pi_indep(x)
        v3 = v3_int(p) - v3_int(q)
        if v % 2 != 0 or v != 2 * v3 or x.v_pi() != v:
            badQ += 1
    report("5f: v_pi restricted to Q is even (v_pi(p/q) == 2*v_3(p/q))",
           badQ == 0, f"{nQ} random rationals, {badQ} failures")


# ---------------------------------------------------------------------------
# 6. v_pi properties
# ---------------------------------------------------------------------------

def check6_valuation():
    n = 500
    bad_mul = 0
    for _ in range(n):
        x, y = rand_Fw(nonzero=True), rand_Fw(nonzero=True)
        if v_pi_indep(x * y) != v_pi_indep(x) + v_pi_indep(y):
            bad_mul += 1
        if (x * y).v_pi() != x.v_pi() + y.v_pi():     # btlib agrees
            bad_mul += 1
    report("6a: v(xy) == v(x) + v(y)", bad_mul == 0,
           f"{n} pairs (independent + btlib), {bad_mul} failures")

    bad_add, bad_eq, neq_cnt = 0, 0, 0
    for _ in range(n):
        x, y = rand_Fw(nonzero=True), rand_Fw(nonzero=True)
        s = x + y
        vs = v_pi_indep(s) if not s.is_zero() else inf
        mn = min(v_pi_indep(x), v_pi_indep(y))
        if vs < mn:
            bad_add += 1
        if v_pi_indep(x) != v_pi_indep(y):
            neq_cnt += 1
            if vs != mn:
                bad_eq += 1
    report("6b: v(x+y) >= min(v(x),v(y)), equality when v(x) != v(y)",
           bad_add == 0 and bad_eq == 0,
           f"{n} pairs ({neq_cnt} with v(x)!=v(y)); "
           f"{bad_add} ultrametric failures, {bad_eq} equality failures")

    v3_ = Fw(3).v_pi()
    ok_v3 = (v3_ == 2) and (v_pi_indep(Fw(3)) == 2)
    abs3 = Fraction(3) ** Fraction(-v3_, 2) if v3_ % 2 == 0 else None
    ok_abs = (abs3 == Fraction(1, 3))
    report("6c: v_pi(3) == 2 and |3|_pi = 3^(-v/2) = 1/3", ok_v3 and ok_abs,
           f"v_pi(3) = {v3_}, |3|_pi = {abs3}")


# ---------------------------------------------------------------------------
# 7. Gates H, S, R
# ---------------------------------------------------------------------------

def check7_gates():
    # i/sqrt(3) == (1+2*omega)/3 numerically and exactly as stored
    ok_is3 = (abs(to_c(I_SQRT3) - 1j / np.sqrt(3)) < 1e-12 and
              I_SQRT3 == Fw(Zw(1, 2), 3))
    report("7a: i/sqrt(3) == (1+2*omega)/3", ok_is3,
           f"|float diff| = {abs(to_c(I_SQRT3) - 1j / np.sqrt(3)):.2e}")

    # build H independently and compare to btlib's H_GATE + numpy
    w = Fw(Zw(0, 1))
    Hind = Mat([[I_SQRT3, I_SQRT3, I_SQRT3],
                [I_SQRT3, I_SQRT3 * w, I_SQRT3 * w * w],
                [I_SQRT3, I_SQRT3 * w * w, I_SQRT3 * w]])
    Hnp = (1j / np.sqrt(3)) * np.array([[1, 1, 1],
                                        [1, W_C, W_C ** 2],
                                        [1, W_C ** 2, W_C]])
    ok_H = (Hind == H_GATE) and all(
        abs(to_c(H_GATE.m[i][j]) - Hnp[i][j]) < 1e-12
        for i in range(3) for j in range(3))
    report("7b: H_GATE matches paper definition (exact + numpy)", ok_H)

    # entries all have v_pi == -1  =>  l(H) = 2
    vs = sorted({v_pi_indep(H_GATE.m[i][j]) for i in range(3) for j in range(3)})
    ok_v = (vs == [-1]) and (ell(H_GATE) == 2)
    report("7c: all 9 entries of H have v_pi == -1, hence l(H) == 2",
           ok_v, f"entry valuations = {vs}, l(H) = {ell(H_GATE)}")

    # unitarity: exact (A* A == I) and numerical
    Idm = Mat.identity(3)
    ok_u = True
    for name, g, gnp in [("H", H_GATE, Hnp),
                         ("S", S_GATE, np.diag([1, W_C, 1])),
                         ("R", R_GATE, np.diag([1, 1, -1]))]:
        exact = (g.conjT() * g == Idm)
        numer = np.allclose(gnp.conj().T @ gnp, np.eye(3), atol=1e-12)
        if not (exact and numer):
            ok_u = False
            print(f"       unitarity failed for {name}: exact={exact}, num={numer}")
    report("7d: H, S, R are unitary (exact A*A == I, + numpy)", ok_u)

    # H^2 == -P23
    P23 = Mat([[1, 0, 0], [0, 0, 1], [0, 1, 0]])
    ok_h2 = (H_GATE * H_GATE == -P23) and np.allclose(Hnp @ Hnp,
              -np.array([[1, 0, 0], [0, 0, 1], [0, 1, 0]], dtype=complex),
              atol=1e-12)
    report("7e: H^2 == -P23 (exact + numpy)", ok_h2)

    # det(H) == -1
    d = H_GATE.det()
    ok_det = (d == -ONE) and abs(np.linalg.det(Hnp) - (-1)) < 1e-12
    report("7f: det(H) == -1 (exact + numpy)", ok_det, f"det(H) = {d!r}")


# ---------------------------------------------------------------------------

def main():
    print("=== v01_ring_basics.py : Section 2 arithmetic foundations ===")
    check0_sanity()
    check1_chi_squared()
    check2_residue()
    check3_ring_equality()
    check4_padic_expansion()
    check5_technical_lemma()
    check6_valuation()
    check7_gates()
    print()
    if FAILURES:
        print(f"OVERALL: FAIL ({len(FAILURES)} failed checks): {FAILURES}")
        sys.exit(1)
    print("OVERALL: PASS (all checks)")
    sys.exit(0)


if __name__ == "__main__":
    main()
