# Data Contract — engine ⇆ stage

All exact arithmetic lives in Python (btlib, frozen + vizgen). The stage (browser)
receives only JSON. No arithmetic is ever re-implemented in TypeScript.

## Vertex addresses (the spine of everything)

The tree is homogeneous: the origin has 4 neighbors; every other vertex has 1
parent + 3 children. Therefore every vertex is named by a digit string:

- `""`  = origin (pure, e₀)
- first digit ∈ {0,1,2,3}, all later digits ∈ {0,1,2}
- parent(addr) = addr[:-1]; depth = len(addr); kind = depth even ? "P" : "A"

Child ordering is canonical: children = btlib neighbors minus parent, sorted by
the canonical HNF vertex key (deterministic across runs/machines — verified in
v09). The stage renders the skeleton procedurally from addresses alone.

Layout (pure function of address, computed stage-side):
- root interval I("") = [0, 2π); child c of v splits I(v) into k equal parts
  (k = 4 at root, else 3), takes the c-th; θ(addr) = midpoint of I(addr)
- hyperbolic radius ρ = depth · λ (λ = ln 3); Poincaré coordinate
  z = tanh(ρ/2)·e^{iθ}

## REST endpoints (engine, localhost)

### GET /api/health
`{ "ok": true, "btlib": "frozen", "ball_check": [1,4,12,36] }`

### POST /api/walk   body: `{ "word": "HSHHHSH" }` or `{ "matrix": [[ "...", ...]] }`
Word letters ∈ {H,S,R} (case-insensitive). Matrix entries are strings parsed as
exact elements of ℤ[ω,1/3]: integers, `w` for ω, `/3^k` denominators, e.g.
`"(1+2w)/3"`. Engine validates unitarity + ring membership exactly.

Response:
```jsonc
{
  "ok": true,
  "word": "HSHHHSH",
  "sde": 0,                       // = l(U)/2 = optimal H-count = distance/2
  "steps": [                      // one per letter, in order
    { "letter": "H", "type": "move",
      "from": "",  "mid": "2", "to": "21" },     // addresses
    { "letter": "S", "type": "fix", "at": "21" } // S/R never move the vertex
  ],
  "trail": ["", "2", "21", "2", ""],  // TREE-level vertex addresses (midpoints
                                      // included): 2 entries appended per H step
  "edgeOwner": [0, 0, 1, 1],          // per trail edge: ordinal of the H that made it
  "geodesic": ["", "2", "21"],        // ALL vertices on the geodesic e₀→Ue₀ (P,A,P,A,…)
  "reduction": [                      // straightening schedule, replayed IN ORDER on
                                      // the evolving trail copy:
    { "index": 2,                     // remove trail[i] and trail[i+1]
      "removed": ["21", "2"] }        // (the two vertex addrs removed, for display)
  ],
  // Invariants (engine-asserted): applying all reduction events to trail yields
  // exactly `geodesic`; an H token dies when both its edges are removed; tokens
  // can half-cancel pairwise (lateral moves) — remaining full tokens = sde.
  "normalForm": {                 // U = (m₁H)(m₂H)…(m_nH)·M
    "blocks": [ { "monomial": "diag(1,ω,1)·P₁₃", "h": true }, ... ],
    "tailMonomial": "−P₂₃"        // M (or "I")
  },
  "matrix": { "entries": [[{ "re": 0.33, "im": 0.66, "str": "(1+2w)/3", "vpi": -1 }, ...]],
              "sde": 0 },
  "error": null                   // or human-readable validation error (ok:false)
}
```
Invariant (verified): len(trail) − 1 = number of H letters; after applying all
`reduction` events to `trail`, the remaining pure vertices are exactly the pure
entries of `geodesic`; sde = (len(geodesic) − 1)/2.

### POST /api/synthesis   body: same as /api/walk
Adds the descent oracle:
```jsonc
{ ...same as walk...,
  "oracle": [   // one entry per geodesic descent step, from Ue₀ back to origin
    { "at": "21..", "candidates": [ { "rep": 3, "lChange": -2, "chosen": true },
                                    ... 12 total, lChange ∈ {-2,0,2} ... ] }
  ] }
```
Invariant (verified, v08): every step has exactly one lChange = −2 (chosen),
two 0, nine +2.

### GET /api/vertex/{addr}
Inspector card (exact data, display-ready):
```jsonc
{ "addr": "21", "kind": "P", "depth": 2,
  "basis": { "entries": [[...fw display...]] },        // canonical HNF basis
  "gram":  { "entries": [[...]] },
  "residue": [[0,1,2],[...],[...]],                    // Gram (P) or χ·Gram (A) mod χ
  "branches": [                                        // edge ↔ 𝔽₃ data
    { "child": null, "toward": "2", "line": [1,1,1] }, // parent edge
    { "child": 0, "toward": "210", "line": [1,1,2] }, ...
  ],                                                   // lines (P) / planes (A)
  "unitary": { "entries": [[...]] } | null             // representative U: Ue₀ = v (pure only)
}
```

### GET /api/random_word?length=40&seed=7  → `{ "word": "HSRH..." }`
### POST /api/equal  body `{ "a": <word|matrix>, "b": ... }` →
`{ "equal": true|false, "dTilde": 4, "vertexA": "21", "vertexB": "21" }`

## Notes for the stage
- Walks are returned whole; the stage owns all animation timing.
- Engine caps word length (default 200) and matrix sde (default 40); errors are
  `{ok:false, error}` — display, never crash.
- All numbers in display payloads are floats/ints/strings; BigInt never crosses
  the wire.
