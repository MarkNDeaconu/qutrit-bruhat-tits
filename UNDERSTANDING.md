# Buildings for Synthesis with Clifford+R — Complete Technical Understanding

*Working notes for the qutrits_v2 project. Companion to `paper.tex` (Deaconu,
Gargava, Kalra, Mosca, Yard) and `btlib.py`. Verification results in
`VERIFICATION.md`; lemma-by-lemma scripts in `verify/`; proof audits in `audit/`.*

---

## 0. One-paragraph summary

The qutrit Clifford+R group ⟨H, S, R⟩ equals U₃(ℤ[1/3, ω]) (arithmeticity, first
proven by Kalra et al. 2025). This paper re-proves it geometrically: it constructs
the Bruhat–Tits building **B** of the unitary group U₃ over the ramified local
field F_π = ℚ₃(ω) as a complex of *self-dual lattice chains*, proves **B is a
(4,4)-biregular bipartite tree**, and shows the gate set acts transitively enough
on the tree (descent from any vertex toward the origin by monomial·H steps) to
generate the whole S-arithmetic group. Exact synthesis = walking a path to the
tree's origin; circuit length = graph distance = 2·sde.

## 1. The arithmetic substrate

- F = ℚ(ω), ω = e^{2πi/3}; O_F = ℤ[ω] (Eisenstein integers).
- χ = 1 − ω. **3 = χ·χ̄ = −ω²·χ²** — the prime 3 is *ramified*: π = χO_F,
  π² = 3O_F. Residue field O_F/π = 𝔽₃ with **ω ≡ 1 (mod χ)**; residue map is
  a + bω ↦ (a+b) mod 3.
- ℤ[1/3, ω] = ℤ[χ⁻¹] since 1/χ = (2+ω)/3 and 1/3 = −ω·χ⁻².
- v_π = χ-adic valuation; v_π(3) = 2; on ℚ₃ ⊂ F_π, v_π is always **even** (ramification
  index 2). v_π = "smallest denominator exponent" (sde) of the synthesis literature.
- Local data: F_π = ℚ₃(ω), O_π = ℤ₃[ω]; uniformizer χ; |x|_π = 3^{−v_π(x)/2}.
- Involution: ω ↦ ω̄ = ω²; fixed rings ℚ, ℤ, ℚ₃, ℤ₃. **Trivial on the residue
  field** (ω ≡ 1). Units of O_π reduce to ±1 mod π; hence for a unit u,
  ū·u ≡ u² ≡ 1 (mod π) — load-bearing fact (Lemma "technical").
- Hermitian form ⟨x, y⟩ = Σᵢ xᵢȳᵢ on F_π³.
- Gates (exact representations): i/√3 = (ω − ω²)/3 = (1+2ω)/3, so every entry of
  H has v_π = −1 exactly, and l(H) = 2. H² = −P₂₃; det H = −1; S = diag(1, ω, 1),
  R = diag(1, 1, −1) are monomial.

## 2. Lattice formalism (the data structures)

A lattice Λ = g·O_π³, g ∈ GL₃(F_π). Everything is decidable in exact arithmetic
over ℚ(ω):

| concept | computational test |
|---|---|
| Λ_g ⊆ Λ_h | h⁻¹g has all entries with v_π ≥ 0 |
| Λ_g = Λ_h | h⁻¹g ∈ GL₃(O_π): integral + det a unit |
| dual Λ♯ | basis (g*)⁻¹ |
| self-dual | **Gram test**: g*g ∈ GL₃(O_π) (Lemma le:gramm_matrix) |
| log₃ index [Λ_g : Λ_h] | v_π(det(g⁻¹h)) |
| canonical vertex key | column HNF over the DVR O_π (digits {0,1,2} per χ-power) |

Dualization: order-reversing involution, det(Λ♯) = (conj det Λ)⁻¹.

**Lemma le:dual_scaling** (parity): if Λ♯ = π^i Λ in odd rank n, then
i = −2k/n where det Λ = π^k, forcing i even, and π^{i/2}Λ is self-dual.
*Mechanism*: det(Λ♯) = π^{−k} vs det(π^iΛ) = π^{ni+k}; π fixed by conjugation.
In even rank this parity argument dies (i = −k in rank 2). Subtlety we verified:
for the *standard* rank-2 form xx̄ + yȳ odd self-duality is still impossible
(the form is anisotropic over F_π because −1 is not a norm: norms are 3^ℤ·(1+3ℤ₃));
genuine odd examples need rank 4 (the standard rank-4 form is isotropic:
1+1+4 = 6 and −6 = 3·(−2) is a norm) or a split form.

## 3. The building B

**Lattice chains** (Abramenko–Nebe model): nested {Λᵢ}, stable under Λ ↦ πΛ
("admissible"). A chain in dimension n has ≤ n π-equivalence classes ("rank"),
because the intermediate lattices between πΛ₀ and Λ₀ form a flag in
Λ₀/πΛ₀ ≅ 𝔽₃³. **B** = all *self-dual* chains (chain fixed setwise by dualization)
of rank ≤ 3. Simplicial structure: 0-simplices = minimal chains; a rank-k chain
is a (k−1)-simplex.

**Anchor** (Appendix A): the dual chain is the original re-indexed by an
order-reversing bijection d(i) = c − i. c even ⇒ some Λⱼ self-dual (anchor j);
c odd ⇒ Λⱼ = Λ♯ⱼ₊₁ (anchor j + ½). Exactly one of the two occurs.

**Classification (Prop pr:simplex_classification)**:
- **Pure vertex** P: {πⁱΛ}, Λ self-dual. (Rank 1; anchor at integer. Anchor at
  half-integer would mean Λ♯₀ = Λ₁ = π⁻¹Λ₀ — odd self-duality, killed by parity.)
- **Alternating vertex** A: {πⁱΛ, πⁱΛ♯} with Λ♯ ⊊ Λ, [Λ : Λ♯] = 9.
  (Rank 2; anchor at half-integer. Rank-2 with anchor at 0 dies: period-2 forces
  Λ♯₁ = πΛ₁, odd self-duality again.)
- **Edge** (1-simplex): {πⁱΛ₁, πⁱΛ₀, πⁱΛ♯₁}: a pure Λ₀ sandwiched
  Λ♯₁ ⊂ Λ₀ ⊂ Λ₁ with indices 3, 3. (Rank 3: one class must be self-dual since
  dualization is an order-reversing involution on 3 classes.)
- **No 2-simplices** (rank ≤ 3 in dimension 3) ⇒ B is a graph.

The graph is **bipartite** (every edge joins P to A) with **origin**
e₀ = {πⁱO_π³}.

### 3.1 Local structure: degree 4 + 4 via 𝔽₃ geometry

**At a pure vertex** Λ: the form ⟨,⟩ mod π on Λ/πΛ ≅ 𝔽₃³ is symmetric (involution
trivial mod π) and non-degenerate (by self-duality, Lemma le:nondegen).
Neighbors = chains πΛ ⊂ πΛ₁ ⊂ Λ♯₁ ⊂ Λ ↔ isotropic lines V ⊂ V^⊥ ⊂ 𝔽₃³
(V = πΛ₁/πΛ, V^⊥ = Λ♯₁/πΛ; the lift is Λ₁ = Λ + χ⁻¹·ṽ·O_π).
**Counting isotropic lines (Lemma le:mod3_codes)**: any non-degenerate ternary
quadratic form over 𝔽₃ has exactly q² − 1 = 8 nonzero isotropic vectors = **4
lines**, *independent of discriminant* (x²+y²+z² and x²+y²−z² both give 8).
For the origin the 4 lines are span(±1,±1,±1).

**At an alternating vertex** Λ (big rep): [Λ : Λ♯] = 9, so Λ♯/πΛ is a line.
The right residue form is **χ⟨,⟩ mod π** — *antisymmetric* because
χ̄/χ ≡ −1 (mod π) (Lemma le:nondegen2; the only place the ramification really
shows its teeth). Via a lower-triangular (HNF) basis one computes χ·Gram mod π =
[[0,a,b],[−a,0,0],[−b,0,0]], (a,b) ≠ (0,0) (else Λ would be self-dual), with
radical = image of Λ♯ (Lemma le:techincal_lemma supplies the entry-by-entry
integrality bootstrap: x̄x + ȳy ∈ O_π ⇒ x, y ∈ O_π, using residues ±1 ⇒ unit
parts sum to 2 ≢ 0 mod 3). Neighbors = self-dual lattices Λ♯ ⊂ Λ₁ ⊂ Λ
↔ Lagrangian-like planes V (radical ⊂ V = V^⊥): the form descends to a
non-degenerate symplectic form on 𝔽₃², where **all** q+1 = **4** lines are
self-dual (Lemma antisym). So alternating vertices also have degree **4**.

Sphere sizes from e₀: 4, 12, 36, 108, 324, … = 4·3^{r−1}. Pure vertices live at
even distance; |pure sphere at 2n| = 4·3^{2n−1}.

### 3.2 Metric structure

𝒜 = {α ∈ GL₃(F_π) : α*α ∈ GL₃(O_π)} = bases of self-dual lattices; pure
vertices = 𝒜/GL₃(O_π). For g ∈ 𝒜, the Cartan decomposition is forced symmetric:
**a = diag(χⁿ, 1, χ⁻ⁿ)** (Lemma le:diagonals; via g* = (g*g)·g⁻¹ and Cartan
uniqueness). Define l(g) = −2·min v_π(gᵢⱼ); then l(kak') = l(a) = 2n,
and **d̃(g,h) = ½[l(g⁻¹h) + l(h⁻¹g)] equals the graph distance** between pure
vertices (Prop prop:distance). Notes:
- l(g⁻¹) = l(g) for g ∈ 𝒜 — but 𝒜 is *not* inverse-closed (gg* ∈ GL₃(O_π) can
  fail), which is exactly why d̃ is defined symmetrized.
- d̃ = 2 ⇔ πΛ_h ⊂ Λ_g ⊂ π⁻¹Λ_h (Lemma le:hecke_nghbs). The hypothetical
  asymmetric case l = (4,0) is killed by a congruence-subgroup trick:
  d⁻¹GL₃(O)d ∩ GL₃(O) is a group, hence inverse-closed.
- Distance from origin: **d(e₀, U·e₀) = l(U) = 2·sde(U)** for U ∈ Γ. This is the
  bridge between tree geometry and circuit complexity.

**Connectedness** (Prop pr:chain): for self-dual Λ_g, Λ_h with minimal n such that
πⁿΛ_g ⊆ Λ_h ⊆ π⁻ⁿΛ_g, a Cartan basis gives Λ_g = ⟨v₁,v₂,v₃⟩,
Λ_h = ⟨χⁿv₁, v₂, χ⁻ⁿv₃⟩; Gram constraints (⟨v₂,v₃⟩ ∈ πⁿO, ⟨v₃,v₃⟩ ∈ π²ⁿO,
⟨v₁,v₃⟩, ⟨v₂,v₂⟩ units) make every interpolant
Lᵢ = ⟨χⁱv₁, v₂, χ⁻ⁱv₃⟩ self-dual — a geodesic of pure vertices.

**Treeness** (Prop pr:forest): B has no cycles. ⚠️ **Audit-confirmed flaw** (see
`audit/forest-proof-audit.md`) — the appendix-B proof constructs a 5-term nested
family and calls it a rank>3 chain, but the family is **not π-stable** (the
wrap-around inclusion A_n♯ ⊆ π^{2n−1}A_n it would need is false: det valuations
−1 vs 6n−2), so nestedness alone contradicts nothing. Decisive soundness check
found by the auditor: the proof's premises never use the loop's *closing* edge,
so the same argument run on a plain geodesic path would "prove" no paths of
length ≥ 3 exist. The *theorem is true* (general Bruhat–Tits theory: U₃ here has
relative rank 1; and our BFS to radius 7 — 4373 vertices — finds zero non-tree
edges). Clean elementary repair, machine-verified in
`verify/forest_proof_audit.py`, using only the paper's toolkit:
1. Two pure vertices at distance 2 have a **unique** common alternating
   neighbor: its big lattice is forced to be S + S′ (index-3 squeeze), with
   small lattice S ∩ S′ = (S + S′)♯.
2. More generally the first alternating vertex on *any* geodesic from Λ_g toward
   Λ_h (distance 2n) is forced: big lattice = Λ_g + π^{n−1}Λ_h. (Cartan basis:
   that sum is ⟨v₁, v₂, χ⁻¹v₃⟩; any geodesic first-step L₁ ⊇ Λ_g with index 3
   and L₁ ⊇ π^{n−1}Λ_h, while minimality of n gives π^{n−1}Λ_h ⊄ Λ_g, forcing
   equality.)
3. Forced first steps ⇒ unique geodesics ⇒ (with bipartiteness) no cycles.

## 4. The main theorem and the synthesis algorithm

Γ = U₃(ℤ[χ⁻¹]) ⊇ ℋ = ⟨H, S, R⟩. Stabilizer of e₀ in Γ is
U₃(O_F) = **monomial matrices** {permutation × diag(±ω^k)}, order 3!·6³ = **1296**
(proof: αᾱ = m²+n²−mn is a positive-definite integer form; unit rows force one
unit entry per row). Two finite computations anchor the proof:
- 𝒢 = Stab(e₀) ∩ Stab(e₁) where e₁ = H·e₀: **#𝒢 = 108** (monomials m with
  H⁻¹mH integral).
- Orbit-stabilizer: #(U₃(O_F)·e₁) = 1296/108 = **12** = #S₀ (pure sphere at
  distance 2: 4 alternating neighbors × 3). So **U₃(O_F) is transitive on S₀**.
  (Robustness note: orbit ⊆ S₀, #orbit = 12, #S₀ ≤ 12 by walk-counting, so
  equality holds *without* invoking treeness.)

**Descent induction**: for any pure v ≠ e₀ at distance 2n, the geodesic's
distance-2 vertex e ∈ S₀ equals m·e₁ for some monomial m; then (mH)⁻¹v is at
distance 2n−2. Induct: v = (m₁H)(m₂H)⋯(m_nH)·e₀. Since pure vertices = Γ·e₀
and stabilizers are monomial ⊆ ℋ, **Γ = ℋ**. ∎
(The proof needs U₃(O_F) ⊆ ⟨H,S,R⟩ — asserted in the paper without proof, and
citing Kalra et al. Cor 5.8 would be circular since that is the theorem being
re-proven. Verified and *repaired* here: ⟨H², S, R⟩ alone has order only 72
(every element fixes the first coordinate axis), but the length-7 palindrome
**T = HSHHHSH = diag(ω,1,ω)·P₁₃** supplies the missing transposition, and
⟨H², T, S, R⟩ = all 1296 monomials. Other certificates: −I = (HSHSH)³,
ωI as a length-21 word. See verify/v07_finite_groups.py, v07_main_theorem.py.)

**Robustness note** (auditor): the main theorem's proof needs *connectedness*
(pr:chain) and *bipartiteness* (classification) but **not** pr:forest — #S₀ = 12
follows from orbit-stabilizer (≥ 12) plus walk-counting (≤ 12) without treeness.
So the Appendix B flaw does not propagate to the main result.

**Algorithm (exact synthesis)**, implemented in `btlib.synthesis_steps`:
```
while l(U) > 0:                  # l(U) = 2·sde(U) = distance to origin
    find m among 12 coset reps with l((mH)⁻¹U) = l(U) − 2   # exactly 1 works
    U ← (mH)⁻¹U; emit (m, H)
emit final monomial U            # l = 0 ⇒ monomial
```
Cost: O(sde) steps × 12 trials × O(1) exact 3×3 algebra. Output length is
**optimal**: number of H-blocks = sde(U) = d(e₀, Ue₀)/2 (geodesic).

## 5. Structural corollaries (the "rich structure" to exploit)

These follow from the paper's results and are the foundation of the project:

1. **Cosets ↔ vertices.** Pure vertices ↔ Γ/U₃(O_F); the tree is the coset graph.
   Alternating vertices ↔ Γ/Stab(A); Stab_Γ(alternating vertex) also has order
   1296 (it surjects onto the symplectic 𝔽₃-data; its monomial part has order 324).
2. **Amalgam / Bass–Serre.** Γ acts on the tree edge-transitively with quotient a
   single edge ⇒ **Γ ≅ G_P ∗_{G_E} G_A** with |G_P| = |G_A| = 1296, |G_E| = 324.
   So the Clifford+R group is *virtually free*, with Euler characteristic
   χ(Γ) = 2/1296 − 1/324 = −1/648: any torsion-free finite-index subgroup of
   index N is free of rank 1 + N/648. Explicit finite presentation extractable.
3. **Normal form & exact counting.** Unique normal form
   U = (m₁H)(m₂H)⋯(m_nH)·M with m₁ ∈ 12 reps, subsequent mᵢ ∈ 9 non-backtracking
   reps, M ∈ U₃(O_F). Hence
   **#{U ∈ Γ : sde(U) = n} = 4·3^{2n−1}·1296 = 12·9^{n−1}·1296** for n ≥ 1
   (15552 for n = 1, 139968 for n = 2, …). Exact "R-count"-style census of the
   gate set, for free, from the tree.
4. **Membership test.** U ∈ Γ ⇔ unitary + entries in ℤ[ω, 1/3] — O(1) exact test;
   synthesis then constructs the witness word (constructive membership =
   arithmeticity in action; contrast with thin groups where this is hopeless).
5. **Geodesic visualization.** Synthesis is literally a walk e₀-ward in a
   4-regular tree; radial embeddings (as in the paper's Fig. 1) animate it.

## 6. Project-facing notes

- `btlib.py` is the exact-arithmetic core (no dependencies; numpy/sympy/networkx/
  matplotlib are installed in `.venv` for cross-checks and future visualization).
- Canonical vertex hashing via HNF over O_π makes BFS, orbit computations, and
  uniqueness checks exact and fast (ball of radius 7 ≈ 4400 vertices in seconds).
- Everything generalizes mechanically to: other vertex stabilizers, Hecke-operator
  experiments (adjacency = Hecke at π), spectral gap / Ramanujan-ness questions,
  covering-rate statistics for *approximate* synthesis, amalgam presentations.
- Floating-point is never needed; ℚ(ω) arithmetic with 3-power denominators is
  closed under all operations used.

## 7. Errata for paper.tex (all machine-verified; none affect the main theorem)

1. **Appendix B (pr:forest): proof invalid as written** — see §3.2 above;
   complete repair in `audit/forest-proof-audit.md`.
2. **Line 663**: "[Λ:Λ♯] = N(det g)⁻¹" is false (evaluates to 3 where brute
   force gives 9). Correct: [Λ:Λ♯] = 3^{−2v_π(det g)}, automatically an even
   power of 3; downstream conclusion (index 9) unaffected.
3. **Lines 1019–1021**: U₃(O_F) ⊆ ⟨H,S,R⟩ asserted, not proven, load-bearing,
   and circular if cited from Kalra Cor 5.8. Short self-contained proof
   available via T = HSHHHSH (§4 above).
4. **Remark after Lemma 2.11 (lines 320–322)**: "fails for even dimensions" is
   wrong for dim 2 with the paper's own standard form (which is *anisotropic*
   over F_π since −1 is not a norm; odd self-duality impossible). Failure
   genuinely starts at dim 4: Λ = O⁴ + χ⁻¹(1,1,1,0)O + χ⁻¹(0,1,−1,1)O has
   Λ♯ = πΛ. Suggested wording: "fails for even dimensions ≥ 4 (and for dim 2
   with isotropic Hermitian forms)".
5. **Lemma le:diagonals (line 943)**: missing hypothesis g ∈ 𝒜 (false otherwise:
   diag(χ,1,1)); proof also needs the stated fact that conj/transpose preserve
   Cartan exponents (χ̄ = unit·χ).
6. **prop:distance (Appendix G)**: proves d̃ ≤ d twice; the direction d ≤ d̃ is
   missing (fix: pr:chain interpolants + le:hecke_nghbs; text in audit report).
7. **pr:chain (line 1396)**: displayed duality identity garbled (false on random
   instances as printed); correct bridge: ((g₁⁻¹h)*)⁻¹ = l₁⁻¹(g₁⁻¹h)l₂.
8. **Appendix C**: sign-convention mismatch (⟨Ax,y⟩ vs displayed expansion =
   xᵀAy; global sign, harmless) and the quotient map (x₁,x₂,x₃) ↦ (x₁,x₂) is
   well-defined only when b = 0; use (x₁, x₂ + (b/a)x₃) in general.
9. **Definition 3.1 (lattice chain)**: Appendix A uses periodicity
   ϖΛᵢ = Λᵢ₋ᵣ ("period 2", etc.), which the stated definition does not supply
   under a weak reading; adopt the standard Abramenko–Nebe periodic form, or
   note self-dual chains are automatically periodic (audit report has the
   one-paragraph lemma).
10. Cosmetic: "le:techincal_lemma" label typo; le:ispositive proof says
    U₃(O_π) where GL₃(O_π) is meant (line 1443); monomial-proposition proof
    conflates rows/columns of A*A = I (line 1031).
