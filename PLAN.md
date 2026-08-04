# Plan: Live trajectory animation (MD / metadynamics synced with simulation steps)

## Context
The demo currently runs the whole MD/metadynamics simulation synchronously
(`run_md_from_sdf`/`run_metadynamics_from_sdf`) and only then animates the
recorded trajectory. The user wants the trajectory shown **lively**, i.e. the
viewer updating while the simulation runs. This requires the JS side to step
the simulation in small chunks (a synchronous WASM call blocks the browser,
so frames can't paint mid-run) — hence a stateful WASM handle that JS drives
with `requestAnimationFrame`.

## Phase 1 — Engine ownership refactor (minimal, needed for a self-contained handle)
- `src/md/mod.rs`: `MDRunner<'a> { ff: &'a dyn ForceField }` →
  `MDRunner { ff: Rc<dyn ForceField> }`; `new`/`from_molecule` take
  `Rc<dyn ForceField>`. Update call sites: `src/lib.rs` (×2), `src/md/mod.rs`
  tests, `src/solvation/mod.rs` test.
- `src/metad/mod.rs`: `MetaDynamics<'a> { ff: &'a dyn ForceField }` →
  `MetaDynamics { ff: Box<dyn ForceField> }`; `new` takes
  `Box<dyn ForceField>`. Update call sites (lib.rs, metad tests).
- Rationale: the live handle must own its force field chain (MMFF →
  MetaDynamics → MDRunner). Rc (not Box) for the runner's FF lets
  `MetaDLive` share the same `MetaDynamics` object it queries (last_cv,
  hill_count, FES).

## Phase 2 — New WASM exports (additive; existing exports unchanged)
`src/lib.rs`:
- `MDLive` (constructor `(sdf, MDOptions)`): `success()`, `error()`,
  `n_atoms()`, `step(n)`, `coords()`, `potential_energy()`, `temperature()`,
  `time_fs()`, `steps_done()`.
- `MetaDLive` (constructor `(sdf, MetaDOptions)`): same + `last_cv()`,
  `hill_count()`, `hill_centers()`, `fes_s(grid)`, `fes_f(grid)`.
- Tests (native, deterministic seed): MDLive steps advance coords/energy
  finite; MetaDLive steps deposit hills, FES arrays match grid size.

## Phase 3 — Live site loop (site/index.html)
- `runMD`: create `MDLive`, push the initial frame, then loop
  `state.step(chunk); pushFrame(); await requestAnimationFrame` (chunk =
  nSteps / ~80–300 target frames). Each frame updates the 3Dmol viewer +
  chart + readouts + status. On completion: enable playback of the recorded
  trajectory.
- `runMetad`: same with `MetaDLive` + per-frame CV value; after the loop,
  compute + draw the FES, show hill count.
- Guard: `runId` token incremented by `loadMol` (and each run start) so a
  preset click mid-run stops the loop cleanly (no viewer model mismatch).
- Unit-aware CV readout (dihedral ° / distance Å).
- Existing synchronous exports stay for API stability (site no longer uses
  them).

## Phase 4 — Docs per workflow
- Overwrite `PLAN.md`; update `CODE_STATUS.md`; README API section: note the
  live handles.

## Verification
- `cargo test` (new + existing, incl. determinism — live stepping must equal
  the synchronous trajectory), `cargo clippy --all-targets` 0 warnings.
- `wasm-pack build --target web --release`, stage `site/index.html`.
- CDP headless Chrome: run MD live — frames advance during the run, slider
  has ~80–300 frames, playback works, zero console errors; run MetaD live —
  same + FES drawn with hill dots. Preset click mid-run stops cleanly.
- Commit, push (Pages redeploys).

## Out of scope
- Changing existing WASM export signatures (kept stable).
- Any changes to integrator/FF numerics (identical step sequence → identical
  trajectory).
