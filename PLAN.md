# Plan: WebMM Playground v1 — interactive physics-toy page (site/playground.html)

## Context
Decision: build the "fun & cool" teaching sandbox (#1) as a standalone page in
this repo, sharing the existing `pkg/webmm.js` build and Pages deploy. RDKit.js
exposes no dynamics, so this is the pure-WebMM showcase. Anchor experience: a
live-MD molecule the user can grab, pull, and heat up, with force-glow coloring.
The current `MDLive` handle is read-only (step + getters only) — three small
additive engine exports make it perturbable. Verified against `src/md/mod.rs`:
`MDRunner` owns `coords`/`vel`/`grad`; MB-init + rescale machinery already exists
in `new()`; forces are the negated cached gradient.

## Phase 1 — Engine: 3 additive MDRunner methods + MDLive WASM wrappers
- `MDRunner::set_atom_position(i, pos)`: set `coords[i]`, zero `vel[i]`
  (kinematic drag), then refresh `pe`/`grad` via one `energy_and_gradient`
  call — the cached gradient corresponds to the old position, and the next
  integrator step would otherwise consume a stale force.
- `MDRunner::rescale_temperature(t)`: factor the MB-init + COM-removal +
  exact-rescale block of `new()` into a private `init_velocities(t)` helper;
  call it from both `new()` and `rescale_temperature` (handles T≈0 restart);
  update `self.temperature_k` so the Langevin thermostat targets the new T.
- `MDRunner::force_magnitudes() -> Vec<f64>`: per-atom `|-grad|` (kcal/mol/Å)
  from the cached gradient.
- WASM wrappers on `MDLive` only: `set_atom_position(i, x, y, z)`,
  `rescale_temperature(t)`, `force_magnitudes()`. No changes to existing
  exports; batch APIs and `MetaDLive` untouched (metad perturbation is v1.5).
- Native tests: (a) pe after set_atom_position equals a fresh eval at the same
  coords; (b) temperature() after rescale_temperature(t) ≈ t (exact rescale);
  (c) force_magnitudes match MMFF analytical gradient on identical coords.
- Gates: `cargo test` green, `cargo clippy --all-targets` 0 warnings,
  `python3 scripts/benchmark_mmff.py --no-speed` 230/230 (MMFF untouched),
  wasm-pack build + `pkg/webmm.d.ts` reflects the new exports.

## Phase 2 — Spike: 3Dmol 2.4 atom picking + screen→world drag (headless CDP)
Riskiest site piece — verify before building the page:
- Atom pick: 3Dmol per-atom click callback (`clickable` + `callback`) fires on
  the rebuilt models the live loop creates.
- Drag mapping: mouse deltas → displacement in the screen-parallel plane
  through the grabbed atom, unprojected via the view quaternion
  (`viewer.getView()`); confirm stable under rotation/zoom.
- Harness: existing CDP setup (python http.server over `pkg/` + headless
  Chrome with SwiftShader).
- If picking proves unreliable in 3Dmol 2.4: fallback = drag the atom nearest
  to the mousedown point projected onto all atom screen positions (manual
  project via the same quaternion math).

## Phase 3 — site/playground.html v1
- Dark glow theme; 3Dmol viewer; energy/temperature chart; temperature slider
  (0–1500 K, live via `rescale_temperature`); run/pause; speed pacing pattern
  copied from index.html (`setTimeout(1000/fps)`, ~120 frames/run chunking).
- Grab & pull: pointerdown near atom → per-frame `set_atom_position` +
  `step(chunk)` + render; release restores normal dynamics.
- Force glow: per-frame atom colors from `force_magnitudes()` (white→red ramp).
- Preset gallery (~8 recognizable molecules as embedded, pre-verified SDFs —
  caffeine, capsaicin, adrenaline, serotonin, dopamine, sucrose, vanillin,
  menthol; each must embed + minimize cleanly; no documented ETKDG outliers,
  no P(=O) compounds).
- Metadynamics section with live FES: ported from index.html as-is
  (MetaDLive batch + FES canvas), re-skinned only.
- Shareable URL: preset id + temperature in the query string.
- Input scope: presets + paste SDF/molblock only (no SMILES — no parser).
- Cross-link: header link back to the demo (index.html); placeholder
  "Open in Molecule Clipboard" link for later.

## Verification (CDP headless)
- Spike passes: atom pick + drag mapping move the intended atom.
- Full page: preset loads → MD runs live; slider mid-run change moves T;
  dragging an atom spikes PE; force-glow colors change frame-to-frame;
  metad FES draws; zero console errors.
- All Phase 1 gates green. Commit, push.

## Docs per workflow
- This file overwrites PLAN.md. After completion: update `CODE_STATUS.md`
  (new exports documented under Recently Completed; Current Focus rewritten);
  README API section gains the three `MDLive` methods + playground link.

## Out of scope
- SMILES parsing, 2D sketcher, descriptors (molecule-clipboard's domain).
- Beat-the-optimizer game, manual dihedral twisting, MetaDLive perturbation
  exports (v1.5 candidates).
- No changes to existing WASM export signatures or index.html behavior.
- No physics changes (dt, integrators, thermostat algorithm, force field).
