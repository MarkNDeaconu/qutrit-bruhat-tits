# Adversarial audit: Prop. pr:simplex_classification (App. A), Lemmas selfdual_oradjacent, le:dual_scaling, le:lattice_chains_are_finite

Auditor scope: paper.tex lines 305-322 (le:dual_scaling + remark), 363-420
(le:lattice_chains_are_finite, selfdual_oradjacent), 1064-1172 (Appendix A).
Sanity computations: exact arithmetic via btlib.py (primitives independently
re-checked: v_pi vs 3-adic valuation of the norm, dual_basis vs the defining
property, det identities); plus the pre-existing verify/v02_duality.py
(43/43 PASS) which independently reaches the same conclusion on the remark.

## Verdict: minor-gap
All *statements* audited are true. The proofs contain a cluster of fixable
gaps centered on chain periodicity, plus one factually wrong (non-load-bearing)
remark about even dimensions.

---

## 1. Lemma le:dual_scaling (lines 305-319) - SOUND

dual(L) = pi^i L (rank 3) ==> i even and pi^{i/2}L self-dual.

* det(dual L) = conj(det L)^{-1}: follows from line 282 (dual basis is
  (A^*)^{-1}); det((A^*)^{-1}) = conj(det A)^{-1}. As O_pi-modules. Verified
  exactly on 100 random bases. "One can check" is fair.
* Parity: det(L) = pi^k as fractional ideals; conjugation fixes the *ideal*
  pi (conj chi = -w^2 chi, a unit multiple; verified). So
  pi^{2k} = pi^{-3i} ==> 2k = -3i ==> i even. Sound. Note the equation lives
  at the ideal level, where it is unambiguous.
* "The second statement trivially follows": dual(pi^{i/2}L) =
  pi^{-i/2} dual(L) = pi^{i/2} L. Fine.
* Computational: 400-lattice random scan, every pi-equivalent dual had even
  exponent; constructed even examples verified; alternating lattices are
  never pi-equivalent to their dual (v_pi(det) obstruction mod 3).

### Remark lines 320-322 - ERROR (non-load-bearing)
"can be replaced by any other odd number" - TRUE (the parity argument is
2k = -ni, n odd forces i even).
"but fails for even dimensions" - FALSE for n = 2 with the paper's standard
form sum x_i conj(y_i): no rank-2 lattice satisfies dual(L) = pi^{odd} L.
Proof: wlog i = 1, so G = A^*A is Hermitian with chi*G in GL_2(O). Diagonal
entries of G lie in Q_3 cap chi^{-1}O = Z_3 (elements of Q_3 have even
v_pi). Then 3 det G = 3 G11 G22 - N(chi G12) == -N(unit) == 2 (mod 3),
while det G = N(det A) forces unit part == 1 (mod 3) (unit norms of the
ramified extension are == 1 mod 3; -1 is not a norm). Contradiction.
Verified: 2500-lattice search (0 hits) here, and independently in
verify/v02_duality.py P6c.1-P6c.6 (including an exact structural proof via
its own Cartan decomposition).
The remark IS true for n = 4: L = O^4 + chi^{-1}v1 O + chi^{-1}v2 O,
v1 = (1,1,1,0), v2 = (0,1,-1,1) (totally isotropic plane mod 3) satisfies
dual(L) = pi L exactly and no rescaling is self-dual (verified exactly).
General standard-form pattern (math, not load-bearing): odd-modular lattices
decompose into binary hyperbolic-type blocks of det class -1/3, so they
exist in the standard space iff (-1)^{n/2} == 1 mod 3, i.e. iff 4 | n.
Suggested wording: "for even n the parity argument breaks down and the
statement can genuinely fail (e.g. n = 4 with the standard form, or n = 2
for an isotropic Hermitian form); for the standard form it still holds when
n == 2 (mod 4)."

## 2. Lemma selfdual_oradjacent (403-420) - SOUND (informal at one point)

* Existence/rigidity of d: the proof asserts the iff with d(i-1) = d(i)+1
  without justification. Correct argument: (i) by definition of self-dual
  chain, dualization maps the chain set onto itself; it is injective
  (biduality) and order-reversing (Prop. le:dual_equivalence), and the
  indexing i -> Lambda_i is an order isomorphism with Z (proper inclusions);
  (ii) hence d (defined by Lambda_i = dual(Lambda_{d(i)})) is a strictly
  decreasing bijection of Z; (iii) every strictly decreasing bijection of Z
  is i -> c - i (order-automorphisms of Z are translations). d is in fact an
  involution. This is standard; the gap is cosmetic.
* Case split: the unique x with d(x) = 0 is either 2j or 2j+1 - exhaustive.
  The notational jump "d(j) = d(2j - j) = j" parses via
  d(i) = d(2j) + (2j - i) = 2j - i evaluated at i = j; correct. Likewise the
  odd case gives d(j) = j + 1, so Lambda_j = dual(Lambda_{j+1}). Correct.
* "at least one" - in fact exactly one, by parity of c. Fine as stated.
* Computational: for the explicit period-3 chain at the origin, verified
  d(i) = -i exactly (dual(Lambda_i) = Lambda_{-i} for i in [-3,3]).

### Anchor definition (1066-1072) - minor gap
"the unique half-integer": uniqueness is later used ("at most one anchor
point", line 1153) but never proven. One-line fix: c = d^{-1}(0)... more
precisely anchor = c/2 where d(i) = c - i; c is unique, so the anchor is;
the two cases are distinguished by the parity of c, so they are mutually
exclusive. Gap-fixable.

## 3. CENTRAL GAP: periodicity of chains (Def. 353; used at 1154, 1161, 1166)

Definition de:lattice_chain requires the chain to be "preserved under the
action Lambda -> varpi Lambda". Two readings:
 (a) weak: varpi S subset S;
 (b) strong: S stable under the Z-action, varpi S = S.
The appendix repeatedly uses full periodicity varpi Lambda_i = Lambda_{i-r}
for all i, with r = rank ("a rank-2 chain has period 2"; "alternating ...
is forced"; the rank-3 identification of the chain with
{pi^i L0, pi^i L1, pi^i L1^sharp}).

Under reading (a), periodicity is FALSE in general. Explicit counterexample
(verified exactly): S = {pi^t A : t in Z} u {pi^t B : t >= 1}, A = O^3,
B = O + O + chi^{-1} O, ordered
... pi^2 A < pi^2 B < pi A < pi B < A < pi^{-1} A < pi^{-2} A < ...
This is a properly nested Z-sequence with varpi S subset S, rank 2, and the
pi-shift on indices is 1 in one place and 2 in another (not a translation).

Fixes (any one suffices):
 1. Adopt the standard Bruhat-Tits / Abramenko-Nebe definition: there is
    r >= 1 with Lambda_{i-r} = varpi Lambda_i for ALL i.
 2. Read (b): varpi acts bijectively on S; a strictly order-preserving
    bijection of a chain isomorphic to Z is a translation; periodicity
    follows with some r >= 1.
 3. Note self-dual chains are automatically periodic even under (a): for
    each pi-class, T = {t : pi^t L in S} is upward closed; duality carries
    it to a class with scaling set -T, which must also be upward closed;
    hence T = Z for every class, varpi acts bijectively, and 2 applies.
    (So the building B is unaffected; verified that the counterexample S
    above is indeed not self-dual.)

With periodicity, "period = rank" needs the (easy, omitted) observation
that the window elements Lambda_{-r+1}, ..., Lambda_0 in (varpi L0, L0] are
pairwise inequivalent: if Lambda_a = pi^k Lambda_b, k >= 1, with
varpi L0 < Lambda_a <= Lambda_b <= L0, then Lambda_a <= varpi Lambda_b <=
varpi L0, contradiction. Hence r = #classes = rank, which is exactly
"rank-2 chain has period 2" and the rank-3 analogue. (Verified for the
explicit 1-simplex chain: the three window classes are pairwise
inequivalent, v_pi(det) == 0,1,2 mod 3.)

## 4. Lemma le:lattice_chains_are_finite (363-373) - TRUE, proof gappy

(i) "intermediate lattices form a flag of subspaces of k^n": needs
    varpi L0 <= L <= L0 ==> L / varpi L0 is a k-subspace (an R-submodule
    annihilated by varpi), the correspondence is injective and inclusion-
    preserving, and chain elements in the window are totally ordered; hence
    at most n-1 strictly between. Standard; "one can show" acceptable.
    (Verified on the explicit chain: dims 1 and 2, nested.)
(ii) Real gap: bounding the window does not bound the number of classes
    unless every class of the chain has a representative that is a CHAIN
    ELEMENT in (varpi L0, L0]. For Lambda_j > L0 this follows from weak
    closure (take minimal m >= 1 with varpi^m Lambda_j <= L0; then
    varpi^{m-1} Lambda_j is a chain element not contained in L0, hence
    contains it, so varpi L0 < varpi^m Lambda_j <= L0). For
    Lambda_j < varpi L0 the symmetric argument needs negative powers in the
    chain, i.e. periodicity (Section 3 fixes). With the periodic definition
    the lemma is immediate: the chain is {varpi^m Lambda_t : 0 <= t < r},
    so #classes <= r, and the window has r elements with r-1 strictly
    between, so r <= n. Gap-fixable.

## 5. Proposition pr:simplex_classification (Appendix A) - conclusions TRUE

(a) rank 1 (1136-1144): "To maintain rank 1, the only possibility is k=1" -
    justification mislabeled: the actual reason is that for k >= 2,
    pi Lambda_1 = pi^{1-k} L0 is a chain element strictly between the
    consecutive elements L0 < Lambda_1. Conclusion correct. Cosmetic.
    Anchor-1/2 elimination via le:dual_scaling (i = -1 odd): correct.
    "cannot contain any subchains": a proper subchain would have to be
    pi-closed with order type Z inside {pi^i L0}, forcing the whole chain;
    informal but correct.
(b) rank 2, anchor 0 (1150-1156): final contradiction
    dual(Lambda_1) = Lambda_{-1} (from d(i) = -i) = pi Lambda_1 (period 2)
    contradicts le:dual_scaling - correct GIVEN periodicity (central gap).
    The auxiliary sentence "Lambda_1 could not be pi-equivalent to a self
    dual lattice ..." additionally uses (silently) that the self-dual
    scaling would itself be a chain element (periodicity again) and anchor
    uniqueness (Section 2 note). All fixable; no hidden case.
(c) rank 2, anchor 1/2 (1158-1161): "Lambda_0 and Lambda_1 cannot be
    pi-equivalent (le:dual_scaling)" - needs one connecting word: if
    Lambda_1 = pi^{-k} L0 then consecutiveness forces k = 1, so
    dual(L0) = pi^{-1} L0, odd, contradiction. (Alternatively: k even would
    produce a self-dual chain element and an integer anchor.) "Alternating
    between scalings is forced" is exactly the periodicity claim (central
    gap). Otherwise correct.
(d) rank 3 (1163-1168):
    * "Because duality is an involution, at least one of the three classes
      contains a self-dual lattice": compressed. Precise: duality descends
      to the 3 classes (Prop. le:dual_equivalence) as an involution; any
      involution of an odd-size set has a fixed point; a fixed class means
      dual(L) = pi^k L, and by le:dual_scaling (NOT cited here, but
      essential) k is even and pi^{k/2} L is self-dual; it lies in the
      chain by periodicity. Without le:dual_scaling a fixed class need NOT
      contain a self-dual lattice - exactly what happens in the n = 4
      example above. Gap-fixable.
    * "In particular dual(Lambda_1) must be Lambda_{-1}": forced by
      d(i) = -i (anchor at 0): duality sends the successor of Lambda_0 to
      the predecessor of dual(Lambda_0) = Lambda_0. The "in particular"
      linkage to the preceding sentence is a non sequitur as written, but
      the fact is correct. Cosmetic.
    * "is a 1-simplex since it contains the two types of 0-simplices":
      should also note (i) the two contained chains are themselves
      self-dual chains (immediate), and (ii) the other two rank-2
      subchains, on class pairs {[L0],[L1]} and {[L0],[L1^sharp]}, are NOT
      self-dual (duality swaps [L1] and [L1^sharp]), so the rank-3 chain
      strictly contains only 0-simplices, making it a 1-simplex and nothing
      higher. Minor omission, fixable.
(e) rank >= 4 impossible (flag bound) and hence no 2-simplices: correct
    (a 2-simplex would strictly contain a rank-3 chain, hence have rank
    >= 4; a rank-3 chain has no proper rank-3 sub- or superchains since all
    classes are full).

Computational cross-checks (all PASS): the explicit 1-simplex chain at the
origin has the exact stated arrangement with proper inclusions; d(i) = -i
holds; the window contains exactly n-1 = 2 lattices strictly between
pi L0 and L0; alternating-vertex lattices satisfy pi L1 < L1^sharp < L1
with dual not pi-equivalent to L1.

## Issue list (by severity)

1. [ERROR, lines 320-322] Remark "fails for even dimensions" is false for
   n = 2 (standard form); true for n = 4. Non-load-bearing.
2. [GAP-FIXABLE, line 353 + 1154/1161/1166] Periodicity of chains used but
   not established; definition ambiguous; non-periodic weak chains exist.
   Three independent fixes given; building B unaffected (self-dual chains
   are automatically periodic).
3. [GAP-FIXABLE, lines 367-373] le:lattice_chains_are_finite: missing step
   "every class has a chain representative in the window".
4. [GAP-FIXABLE, lines 1163-1166] Rank-3 "involution ==> self-dual class"
   compressed; le:dual_scaling must be invoked (else false, cf. n = 4).
5. [GAP-FIXABLE, line 1067/1153] Anchor uniqueness asserted, not proven
   (one line via d(i) = c - i).
6. [COSMETIC] d's existence/rigidity in selfdual_oradjacent asserted
   (standard); "to maintain rank 1" mislabeled; "in particular" non
   sequitur in rank-3; 1-simplex completeness (no other subchains) omitted;
   flag claim "one can show".

Overall: minor-gap; all results stand, with the remark needing a factual
correction and the periodicity infrastructure needing to be made explicit.
