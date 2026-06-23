# Adversarial audit: Appendix B, Proposition `pr:forest` ("B is a forest")

**Target:** `paper.tex`, lines 1173–1242 (Appendix `se:forest_proof`), proving
Proposition `pr:forest` (stated at line 488): the 0- and 1-simplices of the
building `B` form a forest.

**Verification script:** `verify/forest_proof_audit.py`
(run with `.venv/bin/python`; 41 PASS, 0 FAIL, exit 0; deterministic).

**Verdict: the proposition is true, but the proof in Appendix B is invalid
(fatal error in the final inference).** A cheap repair within the paper's
toolkit exists; it is spelled out and machine-verified below.

---

## 1. What the proof does

Notation (paper's, from the displayed rows at lines 1181–1189): a loop based at
a pure vertex `{pi^i S_0}` traverses alternating vertices `{A_j, A_j#}` and pure
vertices `{S_j}`, with edges (1-simplex chains)

```
row 2j+1:  ... pi A_j < pi S_j     < A_j < S_j     < A_j# < pi^-1 S_j ...
row 2j+2:  ... pi A_j < pi S_{j+1} < A_j < S_{j+1} < A_j# < ...
row 2n+2:  ... pi A_n < pi S_0     < A_n < S_0     < A_n# < ...   (closing edge)
```

So the convention here is `A_j < S < A_j#` with `[A#:A] = 9`, `[S:A] = [A#:S] = 3`
(verified: line 663 of the paper gives the index 9; script checks indices).

From the rows the proof telescopes the "left diagonal"

```
pi A_{j+1} < pi S_{j+1} < A_j   ==>   pi^n A_n ⊂ A_0          (B5, line 1231)
```

and dually `A_0# ⊂ pi^-n A_n#` (B6, line 1232). It then asserts (lines
1235–1241) that

```
... ⊂ (pi^n A_n) ⊂ A_0 ⊂ S_0 ⊂ A_0# ⊂ (pi^n A_n)# ⊂ ...
```

"is a self-dual lattice chain" of "rank greater than 3 (a 2-simplex),
contradicting Proposition `pr:simplex_classification`."

### 1.1 The chain relations themselves are fine (audit item 1)

Tracing the inclusions: row 2j+3 gives `pi A_{j+1} < pi S_{j+1}`, row 2j+2 gives
`pi S_{j+1} < pi A_j# < A_j` (using `pi A# < A` from the 0-simplex chain). So
`pi^{j+1} A_{j+1} ⊂ pi^j A_j`, and telescoping over rows **2 .. 2n+1 only**
yields (B5); (B6) is its dual. Row 1 gives `A_0 ⊂ S_0 ⊂ A_0#`. The script
reproduces every row relation and the telescoping exactly on lattices of the
actual tree (checks `1a`). Two minor caveats:

* **Representative normalization (implicit).** Each row reuses representatives
  fixed by previous rows. For rows 1..2n+1 (a path) this is justifiable: each
  pure vertex has a *unique* self-dual representative (`(pi^i L)# = pi^-i L#`,
  so `i = 0` is forced), and each alternating vertex chain has a *unique*
  index-9 gap `(A, A#)` whose endpoints are a dual pair; any self-dual lattice
  in the edge chain must lie in *that* gap (dualizing a gap containing a
  self-dual lattice returns a gap containing the same lattice; gaps containing
  a given lattice are unique). The paper says none of this, but it is true and
  cheap. (Script check `2(i)`-gap.)
* **The closing row 2n+2 is NOT normalizable for free**: with `A_n` already
  fixed by row 2n+1 and `S_0` by row 1, the closing edge only gives
  `pi^c A_n < S_0 < pi^-c A_n#` for *some* `c ∈ Z`; the displayed `c = 0` is an
  unjustified assumption. Harmless, because row 2n+2 is **never used** in
  deriving (B5)/(B6) — which is itself the root of the disaster below.

### 1.2 The fatal error: the displayed object is not a lattice chain (audit item 2)

Definition `de:lattice_chain` (lines 340–354) requires a lattice chain to be a
**totally ordered** nested family that is **preserved under `L -> pi L`**
(admissibility is built into the paper's definition; see Remark
`re:lattice_chains_admissible`). Both `pr:simplex_classification` ("no
2-simplices", via Lemma `le:lattice_chains_are_finite`, whose proof at line
368 *uses* `pi L_0` being in the chain) and the notion of "rank" only apply to
such chains. The exhibited 5-term display fails this, **provably and
unavoidably**:

* **Wrap-around impossibility.** A pi-stable totally ordered chain containing
  the period `pi^n A_n ⊂ ... ⊂ (pi^n A_n)# = pi^-n A_n#` must satisfy
  `pi · (pi^n A_n)# ⊆ pi^n A_n`, i.e. `A_n# ⊆ pi^{2n-1} A_n`. Since
  `A_n ⊊ A_n#` with `[A_n# : A_n] = 9` (det-valuations: `v(det A_n#) = -1`
  vs. `v(det pi^{2n-1}A_n) = 6n-2`), this is impossible for every `n >= 1`.
* **Total order forces `A_0 = A_n`.** The pi-closure of the display contains
  both `A_0` and `A_n` (= `pi^{-n}·(pi^n A_n)`), which have *equal* determinant
  valuation (`v = 1` when `v(det S) = 0`). Comparable lattices with equal
  determinant valuation are equal. So a total order would force `A_0 = A_n`,
  contradicting the proof's own hypothesis that the loop uses two *distinct*
  alternating vertices. Hence **no pi-stable totally ordered family containing
  the five displayed pi-classes exists, ever** — the "self-dual lattice chain
  of rank greater than 3" cannot be formed, and no contradiction with
  `pr:simplex_classification` is obtained.
* **The proof proves too much.** Since (B5), (B6) and row 1 are derived without
  the closing edge, the identical premises hold along any geodesic **path**
  `S_0 – A_0 – S_1 – ... – S_n – A_n` in the true tree. If the final inference
  were valid, it would equally show that the building has no geodesics of
  length >= 3 — i.e. that the infinite (4,4)-biregular tree the paper itself
  constructs is empty. The script builds such paths for `n = 1, 2` and
  machine-verifies: all premises hold (checks `1a`, `1b`: properly nested,
  five distinct pi-classes, self-dual as a finite sequence), while the
  conclusion object fails Definition `de:lattice_chain` (checks `1c`:
  wrap-around fails; `A_0` vs `A_n` incomparable; 28 resp. 50 incomparable
  pi-translate pairs found in a window).
* Under the charitable non-admissible reading ("lattice chain" = mere nested
  sequence, à la Abramenko–Nebe), the display *is* a nested rank-5 sequence —
  but then it contradicts nothing: `pr:simplex_classification` and
  `le:lattice_chains_are_finite` only constrain admissible chains, and
  non-admissible nested sequences of arbitrary length exist trivially (e.g.
  along any geodesic, as just demonstrated).

There is also no quick salvage using the closing-edge relations: the
total-order obstruction above (`A_0 = A_n`) is valuation-theoretic and
independent of any additional inclusions a loop would supply.

**Conclusion: the proof is flawed, not merely under-justified.** The claim is
still true (rank-1 Bruhat–Tits theory), and a self-contained repair is below.

---

## 2. The repair (audit item 3) — verified

All ingredients are in the paper already: `pr:chain`'s Cartan basis (Appendix,
lines 1383–1418), the isotropic-line neighbor classification (lines 592–611 +
`le:mod3_codes`), `le:hecke_nghbs`, `prop:distance`, `le:dual_scaling`,
bipartiteness (line 483). Distances below are graph distances `d` (=
`d~ = (l(g^-1 h)+l(h^-1 g))/2` on pure vertices by `prop:distance`; the script
re-verifies `d~ == BFS distance` on the radius-5 ball, check `0g`).

**Lemma A (canonical representatives).** Each pure vertex has a unique
self-dual lattice `S` (if `pi^i S` were also self-dual then `i = 0`). Each
alternating vertex `{pi^i A, pi^i A#}` has a unique index-9 gap whose endpoints
are a dual pair `(A, A#)`, and every self-dual lattice lying in the vertex's
chain-gaps lies in *that* gap: if `B ⊂ S ⊂ T` with `(B,T)` a gap, dualizing
gives the gap `(T#, B#)` containing `S`, and gaps containing `S` are unique, so
`T = B#`. Consequently if pure `S != S'` are both adjacent to `{A, A#}`, then
`A ⊆ S ∩ S' ⊊ S` with `[S : A] = 3` forces `A = S ∩ S'`, and dually
`A# = S + S'`. In particular **two pure vertices at distance 2 have exactly one
common alternating neighbor** (the candidate fix's step (i): index argument,
`Lambda_1 = S + S'`). *Machine-verified* for all alternating vertices of the
radius-5 ball and all 240 pure pairs at distance 2 (checks `2(i)`).

**Lemma B (strong forced first step; fix's step (ii), strengthened).** Let
`Lambda_g, Lambda_S` be self-dual with `d = 2k >= 2`, and take the Cartan basis
of `pr:chain`: `Lambda_S = <v1,v2,v3>`, `Lambda_g = <chi^k v1, v2, chi^-k v3>`,
with Gram facts `<v2,v3> ∈ pi^k O`, `<v3,v3> ∈ pi^{2k} O`, `<v1,v3>, <v2,v2>`
units (proved at line 1408). Then **any** alternating neighbor `{A, A#}` of `S`
(canonical reps, `A ⊂ Lambda_S ⊂ A#`) with `pi^k Lambda_g ⊆ A` equals the
vertex with

```
A = <chi v1, v2, v3>,      A# = <v1, v2, chi^-1 v3> = Lambda_S + pi^{k-1} Lambda_g .
```

*Proof.* Neighbors of `S` correspond to isotropic lines `V` in
`Lambda_S/pi Lambda_S` (paper, lines 592–611), with `A` = preimage of
`V^perp`. The image of `pi^k Lambda_g = <chi^{2k}v1, chi^k v2, v3>` mod
`pi Lambda_S` is the line `<v3>` (nonzero because the Cartan exponent is
exactly `k`, i.e. the distance is exactly `2k`; this needs `k >= 1`). So
`pi^k Lambda_g ⊆ A` iff `v3 ∈ V^perp` iff `V ⊆ v3^perp = span(v2, v3)` (row 3
of the reduced Gram is `(unit, 0, 0)`). On that plane the form is
`(b v2 + c v3) · (b v2 + c v3) = b^2 <v2,v2> != 0` unless `b = 0`; the only
isotropic line is `<v3>`. Hence `V = <v3>` is forced, giving the displayed `A`,
`A#`. Moreover `T = <chi v1, v2, chi^-1 v3>` is self-dual (Gram entries
integral using `k >= 1`, determinant a unit; Lemma `le:gramm_matrix`) and is a
neighbor of this vertex with `pi^{k-1} T ⊆ Lambda_g ⊆ pi^{-(k-1)} T`, so
`d(g, T) <= 2k-2` and `d(g, {A,A#}) <= 2k-1`. ∎

*Machine-verified*: for **all** pure vertices at distances 2 and 4 in the ball
(120 vertices) and random vertices at distances 6 and 8, exactly one
alternating neighbor satisfies `pi^k O^3 ⊆ A`; it is exactly the unique
neighbor at distance `2k-1`, and its big lattice equals
`Lambda_v + pi^{k-1} O^3` (checks `2(ii)`).

**Lemma C.** If `x = {A, A#}` is alternating with `d(g, x) = 2n-1` (`g` pure,
canonical reps), then `pi^n Lambda_g ⊆ A`. *Proof:* `x` has a pure neighbor `T`
at distance `2n-2`, so `pi^{n-1} Lambda_g ⊆ T` (Cartan sandwich, easy direction
of `pr:chain`/`prop:distance`), and `T ⊆ A#` with `pi A# ⊆ A` gives
`pi^n Lambda_g ⊆ pi T ⊆ pi A# ⊆ A`. ∎ (Both `T` reps lie in the same gap by
Lemma A.) *Machine-verified* on all 364 alternating vertices of the ball
(check `2(iii)`).

**Theorem (no cycles).** Suppose `C` is a simple cycle. `B` is bipartite
(every 1-simplex joins one pure and one alternating 0-simplex, by
`pr:simplex_classification`; and there are no double edges, since a 1-simplex
chain is the union of its two 0-simplex chains, so an edge is determined by its
endpoints). Fix a pure vertex `g` on `C` and let `f` be a vertex of `C` at
maximal distance `D = d(g, f) >= 1`. Its two cycle-neighbors `u != w` satisfy
`d(g, ·) <= D` and, by bipartite parity, `d(g, u) = d(g, w) = D - 1`.

* **Case B: `f` pure, `D = 2n` (`n >= 1`).** Then `u, w` are alternating
  neighbors of `f` at distance `2n-1` from `g`. By Lemma C,
  `pi^n Lambda_g ⊆ A_u` and `pi^n Lambda_g ⊆ A_w`; by Lemma B (applied to the
  pair `(f, g)` at distance `2n`) both equal the *one* forced vertex, so
  `u = w` — contradicting simplicity of `C`.
* **Case A: `f = {A, A#}` alternating, `D = 2k+1`.** Then `u, w` are pure,
  distinct, at distance `2k`. If `k = 0` then `u = w = g`, contradiction; so
  `k >= 1`. By the distance-`2k` sandwich, `pi^k Lambda_g ⊆ u ∩ w`, and by
  Lemma A `u ∩ w = A`. Thus `f` is an alternating neighbor of `u` with
  `pi^k Lambda_g ⊆ A`, so by Lemma B `f` is the forced vertex and
  `d(g, f) <= 2k - 1 < D` — contradicting maximality of `D`. ∎

This settles the question raised in the assignment about "the pure vertex
between two alternating vertices": it does **not** need to be forced directly
(there are indeed 3 non-backtracking pure choices below an alternating vertex,
all legitimate); the farthest-vertex argument only ever needs the *closer*
neighbor to be unique, and the alternating-farthest case (Case A) reduces to
Lemma B via the intersection identity `A = u ∩ w` of Lemma A. The
"two distinct neighbors of the farthest vertex both reduce distance" argument
closes cleanly in both parities. *Machine-verified consistency*: for every
alternating vertex in the radius-5 ball the neighbor-distance multiset is
exactly `{d-1, d+1, d+1, d+1}` — the Case-A "bad configuration" never occurs
(check `2(iv)`); the whole radius-6 ball is a tree with unique closer
neighbors (check `0g`).

The repair needs only one glue fact not displayed verbatim in the paper:
"`d(g,h) = 2n` between pure vertices iff the minimal sandwich exponent /
Cartan invariant is `(n, 0, -n)`", which is exactly the content assembled from
`prop:distance` + the first paragraph of the `pr:chain` proof (lines
1384–1397) + `le:ispositive`; the script re-verifies it numerically
(`d~ == BFS distance`, check `0g`).

So the **candidate fix is verified correct** (steps (i), (ii), (iii)), with the
one refinement that step (ii) must be stated in the strong form of Lemma B
("any alternating neighbor whose small lattice contains `pi^k Lambda_g` is the
forced one"), which is what both Case A and Case B actually consume.

---

## 3. Findings summary

| # | Location | Severity | Finding |
|---|----------|----------|---------|
| 1 | lines 1235–1241 (final inference of the `pr:forest` proof) | **error** | The 5-term nested display is asserted to be "a self-dual lattice chain of rank greater than 3". Per Definition `de:lattice_chain` a chain must be pi-stable and totally ordered; the display's pi-closure is never totally ordered for distinct `A_0, A_n` (wrap-around needs `A_n# ⊆ pi^{2n-1}A_n`, impossible for `n>=1`; total order forces `A_0 = A_n` by determinant valuation). No contradiction with `pr:simplex_classification` arises. The premises hold on plain geodesic paths, so the argument would also "prove" the tree has no paths of length >= 3. The proof does not establish the proposition. |
| 2 | lines 1181–1189 (displayed rows) | gap-fixable | Representative normalization across rows is implicitly assumed. Justifiable for rows 1..2n+1 via the unique-self-dual-gap lemma (not stated in the paper); the closing row 2n+2 (`c = 0` scaling) is not justifiable, but is unused. |
| 3 | lines 1181–1241 | cosmetic | Properness of the inclusions (needed even for the rank count) is used but only follows from distinctness of the loop's vertices; stated only obliquely ("two distinct alternating vertices"). |

Proposition `pr:forest` itself is **true** (and locally confirmed by exhaustive
BFS to radius 6: 1457 vertices, no cycles), and the repair in §2 proves it with
the paper's own toolkit.

**Script:** `verify/forest_proof_audit.py` — 41 PASS / 0 FAIL, exit 0,
runtime ~2 s, deterministic (seed 12345). btlib primitives were independently
cross-checked against complex floats, norm-based valuations, semantic duality
pairings, and brute-force neighbor enumeration before being relied upon
(checks `0a`–`0g`).
