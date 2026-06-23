# Running "The Tree Behind the Gates" (MVP — M0+M1)

## Quick start

```bash
cd ~/projects/qutrits_v2
make demo        # builds stage if needed, starts engine, opens the browser
```

Manual equivalent: `make stage` (once), then `make engine`, then open
http://127.0.0.1:8137. Dev loop with hot reload: `make engine` in one terminal,
`make stage-dev` in another (Vite on :5173 proxies /api to the engine).

Everything runs locally in `.venv` + `stage/node_modules`; no network needed
after the initial `npm install`.

## What's in the MVP

- **Poincaré-disk tree browser**: the ball of radius 8 (13,121 vertices) drawn
  procedurally — vermilion pure / azure alternating vertices, gold origin halo,
  boundary rim = the Cantor set of ends. Drag to pan, wheel to zoom,
  double-click a vertex to focus, `0` resets, `L` toggles light theme
  (projector-friendly).
- **The flagship walk** (Play / Append / Synthesize):
  - **Play** starts fresh from the origin: snaps the camera home, animates the
    typed word as a walk (H = pulse across two edges through the midpoint; S/R =
    compass spin in place), and marks the vertex it lands on with a **pulsating
    green node** (= U·e₀ for the typed word).
  - **Append** (green button, key `A` or Shift+Enter) continues from the green
    node, adding whatever is in the text box as a further segment; the green
    node moves to the new endpoint. It extends *whatever path is on screen* —
    including a freshly **synthesized** geodesic — so synthesize-then-append
    grows the clean path from its end. Everything else (Play, examples, Random)
    starts over from the origin — only Append continues.
  - **Synthesize** snaps the camera back to the origin, collapses the messy
    trail to the **geodesic** (lobes burning off, gate tokens annihilating),
    then **plays a pulse tracing the clean path** from the origin to the green
    node. Same destination, optimal route: H-count = distance/2 = sde
    (optimality banner). Works on the full accumulated word (all appends).
  - The text box is the *scratch input* and clears after each Play/Append; the
    accumulated word lives in the green node + scrubber tokens.
  - You can also paste an exact 3×3 matrix over ℤ[ω,1/3] (Matrix… button);
    Append is disabled in matrix mode (a pasted matrix has no word suffix).
- **Inspector**: click any vertex → exact HNF basis, Gram matrix, residue form
  mod χ, and (pure vertices) a representative unitary. All exact strings with
  valuation badges, computed live by the engine.
- **Engine** (FastAPI, 127.0.0.1:8137): every number on screen comes from the
  frozen verified `btlib.py` via `vizgen/` — faithfulness machine-checked by
  `verify/v09_vizgen.py` (31 checks).

## Visual QA checklist (needs human eyes — please walk through once)

1. Load: tree fades in by depth (~600 ms); 60 fps when idle; origin halo
   "breathes" gently.
2. **Scenic example "collapses to 4 edges"** → Play: a long walk sprawls out in
   several directions and back over itself; a **green node** lights up and
   pulsates where it lands. Synthesize: camera snaps to the origin, the tangle
   collapses to a clean 4-edge geodesic running straight to the green node;
   banner reads sde = 4. This is the headline contrast (crazy → clean).
3. **Append flow**: Play a short word (e.g. `HSH`); green node appears. Type
   `RHH`, click **Append** (or press `A`): the walk continues from the green
   node along the new segment, green node moves. Repeat a couple times. Then
   **Synthesize**: it snaps to the origin and draws the shortest path to the
   final green node for the *whole* accumulated word.
4. **Play resets, Append continues**: after building a path with appends,
   clicking Play (or an example) wipes it and starts fresh from the origin;
   only Append continues from the green node.
5. **T = HSHHHSH** → Play: 5 H-pulses wander; green node lands back at the
   origin. Synthesize: all tokens annihilate, sde = 0, tail monomial
   `diag(w,1,w)·P(321)`. (The length-7 palindrome that closes the
   monomial-generation gap — good talking point.)
6. Click vertices at several depths: inspector cards render; origin card works;
   Esc closes.
7. `L`: light theme readable; green node + trail still visible.
8. Disconnect the engine (Ctrl-C) with the page open: toast/status shows
   engine unavailable; no crash.

## Curated talking points (all badges true)

- "Every pixel is exact arithmetic — no floats anywhere in the math path."
- The 1–2–9 split: at every vertex exactly one of the 12 step choices descends
  (machine-verified at 690/690 steps + 37/37 in v09).
- Tree distance counts **H** — the *denominator* gate, which is Clifford. This
  inverts the Clifford+T folklore; never say "T-count."
- Census: #{U : sde = n} = 12·9^(n−1)·1296 exactly (15,552 at n=1 — fully
  enumerated).

## Roadmap next (per approved plan)

M2: descent-oracle overlay + 𝔽₃ compass + equality tester UI.
M3: two-worlds split view; universality toggle (certify ⟨H,S⟩ first via
`vizgen/orbits.py` — not yet written); elliptic/hyperbolic axes; identity vault.
M4: census dashboard, random-walk lab, boundary/amalgam overlays, presenter
keys 1–9. M5: approximation sandbox, exact naive-sampler drift chain.
