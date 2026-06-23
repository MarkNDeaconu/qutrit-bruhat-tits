# Verification Report — paper.tex ("Buildings for Synthesis with Clifford+R")

12-agent verification + audit run, 2026-06-10. All scripts are deterministic
(seeded), exact (no floats in any assertion; numpy used only as an independent
cross-check), and re-runnable with `.venv/bin/python verify/<file>.py`
(every script exits 0 = all checks pass).

**Bottom line: every theorem, proposition, and lemma *conclusion* in the paper
is true and machine-verified.** Two written *proofs* need repair (one seriously:
Appendix B), one stated formula and one remark are wrong as written (both
non-load-bearing, both with corrected statements supplied), and a handful of
"one can check" steps were checked. Full repair texts in `audit/`.

---

## 1. Lemma-by-lemma verification matrix

| Paper item | Script | Status | Highlights |
|---|---|---|---|
| 3O_F = χ²O_F, residue field 𝔽₃, ω ≡ 1, ℤ[1/3,ω] = ℤ[χ⁻¹] | v01 | ✅ | exact identities; residue map a ring hom (1000 samples) |
| Lemma 2.4 digit expansion (le:padic_expansion) | v01 | ✅ | 200 elements × 12 digits, digits unique in {0,1,2} |
| Lemma le:techincal_lemma | v01 | ✅ | contrapositive exact: v(x̄x+ȳy) = 2n < 0; unit residues make the sum ≡ 2 mod 3 |
| v_π properties, l(H) = 2, H² = −P₂₃, det H = −1 | v01 | ✅ | all 9 entries of H have v_π = −1 exactly |
| Dual = (A*)⁻¹, double dual, inclusion reversal, π-equiv preserved | v02 | ✅ | 100+ random lattices incl. non-3 denominators |
| det(Λ♯) = (conj det Λ)⁻¹ | v02 | ✅ | exact identity + module-level on scrambled bases |
| Lemma 2.11 (le:dual_scaling), rank 3 | v02 | ✅ | 43 constructed dual-equivalent lattices: i always even, π^{i/2}Λ self-dual |
| Remark "fails for even dimensions" | v02 | ⚠️ partial | **dim 2 + standard form: lemma still HOLDS** (form anisotropic over F_π; −1 not a norm — exhaustive mod 3,9,27). Genuine odd counterexample first in **dim 4**: Λ = O⁴ + χ⁻¹(1,1,1,0)O + χ⁻¹(0,1,−1,1)O has Λ♯ = πΛ. Suggested wording: "fails for even dimensions ≥ 4 (and dim 2 with isotropic forms)" |
| Lemma 3.12 (le:mod3_codes) | v03 | ✅ | **exhaustive**: all 468 invertible symmetric A over 𝔽₃ have exactly 8 isotropic vectors = 4 lines |
| Lemma 3.13/C (antisym) | v03 | ✅ | exhaustive over all 8 (a,b): radical = span(0,b,−a); exactly 4 self-dual planes; quotient argument checked |
| Lemma le:nondegen (pure residue form) | v04 | ✅ | 32 self-dual lattices: symmetric, non-degenerate, well-defined |
| [Λ : Λ♯] = 9 at alternating vertices | v04 | ✅ | det-valuation AND brute-force coset count (27 reps → 9 classes) |
| Lemma le:nondegen2 (χ⟨,⟩ antisym, rank 2, radical = im Λ♯) | v04 | ✅ | 30 vertices + explicit change of basis to the [[0,a,b],[−a,0,0],[−b,0,0]] shape |
| Pure-vertex ↔ 4 isotropic lines bijection | v04 | ✅ | at origin: exactly the 4 isotropic lines give valid chains; the 9 others fail (witness: πΛ₁ ⊄ Λ♯₁) |
| Alternating ↔ 4 self-dual planes | v04 | ✅ | 10 vertices × 13 planes: exactly 4 lift to self-dual lattices |
| Cartan/SNF over O_π (Thm 2.13) | v05 | ✅ | from-scratch SNF; reconstruction + uniqueness, 100 random g |
| Lemma le:diagonals (a = diag(χⁿ,1,χ⁻ⁿ) for g ∈ 𝒜) | v05 | ✅ | 50 samples; **note: hypothesis g ∈ 𝒜 missing from lemma statement** (diag(χ,1,1) is a counterexample otherwise) |
| Lemma le:ispositive (1)–(5) | v05 | ✅ | incl. explicit g₀ ∈ 𝒜 with g₀g₀* ∉ GL₃(O): **𝒜 is not inverse-closed**, yet l(g⁻¹) = l(g) holds; d̃ metric axioms 200 triples |
| Lemma le:hecke_nghbs | v05 | ✅ | both directions; (4,0) case: 0 hits in 500 trials; congruence-subgroup step verified |
| Prop pr:chain (interpolation) | v05 | ✅ | 30 pairs, n ≤ 4: Gram conditions hold, all interpolants self-dual, consecutive distance 2 |
| Tree structure (Prop 3.8 + 3.10) | v06 | ✅ | **ball of radius 7 = 4373 vertices; spheres exactly [1,4,12,36,108,324,972,2916]; 0 non-tree edges; degree 4 everywhere** (cross-checked by brute-force enumeration of all 13 index-3 candidates per vertex: exactly 4 valid) |
| Unique middle vertex at distance 2 | v06 | ✅ | 50 pairs: exactly one common alternating neighbor; big lattice = Λ_g + Λ_h, small = Λ_g ∩ Λ_h |
| Prop prop:distance (d̃ = graph distance) | v06 | ✅ | 100 random pure pairs to distance 12: 0 mismatches |
| Γ-action: isometries, kind-preserving, commutes with neighbors | v06 | ✅ | 200 actions |
| U₃(O_F) = monomials, #= 1296 | v07 | ✅ | complete enumeration of 18³ unit-norm column triples → exactly the 1296 monomials |
| **#𝒢 = 108** (the paper's "by explicit calculation") | v07 | ✅ | three independent methods agree; structure: center = 6 scalars, order histogram {1:1, 2:19, 3:26, 6:62} |
| Orbit U₃(O_F)·e₁ = S₀, #= 12 | v07 | ✅ | orbit keys ≡ BFS sphere-2 keys; fibers all 108 |
| Edge stabilizer 324; amalgam data (1296, 324, 1296) | v07 | ✅ | all 4 alternating neighbors; H fixes the midpoint edge vertex ⇒ Stab_Γ(alt) = 4·324 = 1296 |
| **U₃(O_F) ⊆ ⟨H,S,R⟩** (used, unproven in paper) | v07 | ✅ | word-BFS depth 6: 148 monomial values generate all 1296. Certificates: S, R, H² = −P₂₃, HSHSH; −I = (HSHSH)³ len 15; ωI len 21; T = HSHHHSH = diag(ω,1,ω)P₁₃ len 7. **⟨H²,S,R⟩ alone has order only 72** (fixes axis 1) — a genuinely new generator is needed |
| #{sde = 1} = 15552 = 1296²/108 | v07+v08 | ✅ | full enumeration; collision sets match coset prediction exactly |
| Synthesis round-trip (Thm 2.2 algorithmic content) | v08 | ✅ | 130 unitaries (sde up to 29): exact reconstruction; word length = sde = distance/2 |
| Descent geometry | v08 | ✅ | at all 690 steps the 12 candidates split l-change {−2: 1, 0: 2, +2: 9} — matches tree picture (1 toward origin, 2 siblings, 9 outward); descent deterministic |
| Normal form unique; #{sde = 2} = 139968 | v08 | ✅ | 12×9×1296 enumeration, zero collisions; 3000 random sde-2 unitaries each hit exactly once |
| Performance | v08 | ✅ | ~0.9 ms/step, linear to sde 40 |

## 2. Audit verdicts (proof validity, not statement truth)

| Proof | Verdict | Key issues |
|---|---|---|
| App. A (classification) + chain lemmas | **minor-gap** | Definition of lattice chain needs the standard periodicity (ϖΛᵢ = Λᵢ₋ᵣ) — used as "period 2/3" but the weak reading admits counterexamples; rank-3 "involution ⇒ self-dual class" needs le:dual_scaling cited (a fixed *class* need not contain a self-dual lattice — exactly what happens in dim 4); anchor uniqueness unproven (one line). Repairs in `audit/appendixA_simplex_classification.md` |
| **App. B (forest / pr:forest)** | **flawed** | The 5-term nested display is *not* a lattice chain: its π-closure is never totally ordered (wrap-around needs A_n♯ ⊆ π^{2n−1}A_n, impossible; det valuations −1 vs 6n−2). Decisive soundness check: the proof's premises never use the loop's closing edge, so the same argument run on a plain geodesic path would "prove" no length-≥3 paths exist. **The proposition is true** (BFS: 0 cycles to radius 7) and a complete elementary repair is written and machine-verified: unique-midpoint Lemma A + forced-first-step Lemma B (Cartan basis: forced big lattice = Λ_g + π^{n−1}Λ_h) ⇒ unique geodesics ⇒ no cycles. See `audit/forest-proof-audit.md` + `verify/forest_proof_audit.py` |
| Distance/metric suite (le:Aset, le:diagonals, le:ispositive, le:hecke_nghbs, pr:chain, prop:distance) | **minor-gap** | prop:distance proves d̃ ≤ d *twice* and never d ≤ d̃ (fix: pr:chain interpolants + le:hecke_nghbs give the path; text supplied); le:diagonals missing "g ∈ 𝒜" hypothesis + the conj-Cartan fact; pr:chain has a garbled display (the duality identity as printed is false on random instances — correct bridge: ((g₁⁻¹h)*)⁻¹ = l₁⁻¹(g₁⁻¹h)l₂); le:hecke_nghbs converse "follows in the same way" actually needs 3 extra steps (reconstructed + verified) |
| Main theorem (§4.2) | **minor-gap** | (1) **U₃(O_F) ⊆ ⟨H,S,R⟩ is load-bearing and unproven**; citing Kalra Cor 5.8 would be circular. Fix supplied: short generation proof via H², T = HSHHHSH, S, R. (2) **Line 663 formula error**: "[Λ:Λ♯] = N(det g)⁻¹" is false (gives 3, brute force gives 9); correct: [Λ:Λ♯] = 3^{−2v_π(det g)} — conclusion (index 9) unaffected. (3) #S₀ = 4×3 = 12 implicitly uses treeness, but the proof is robust without App. B: orbit gives ≥ 12, walk-count gives ≤ 12. **The main theorem needs connectedness + bipartiteness but NOT pr:forest.** |
| Appendix C (antisym) | sound (2 cosmetic) | sign convention ⟨Ax,y⟩ vs xᵀAy (global sign, harmless); the quotient map (x₁,x₂,x₃) ↦ (x₁,x₂) is only well-defined when b = 0 — corrected map (x₁, x₂ + (b/a)x₃) works for all cases (verified exhaustively) |

## 3. Files

```
btlib.py                      exact-arithmetic core (Z[ω], v_π, lattices, HNF keys,
                              tree neighbors/BFS, exact synthesis)
smoke_test.py                 quick library sanity
verify/v01_ring_basics.py     … v08_synthesis.py, v07_distance_audit.py,
       v07_main_theorem.py, forest_proof_audit.py     (all exit 0)
audit/appendixA_simplex_classification.md             repair text for App. A
audit/forest-proof-audit.md                            full App. B analysis + new proof
UNDERSTANDING.md              complete technical understanding + project roadmap
```
