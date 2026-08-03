# Plan: Fix demo page WASM API bugs (Embed / Optimize / Metadynamics buttons)

## Context
User reports three runtime errors on the deployed demo (site/index.html):
1. **ETKDG Embed**: `Error: res.success is not a function`
2. **Run Metadynamics**: `MetaD error: Cannot read properties of undefined (reading 'toFixed')`
3. **Optimize (MMFF94s)**: `Error: res.coordinates is not a function`

Root-caused by comparing site calls against the generated bindings
(pkg/webmm.d.ts) and REPRODUCED in node against the real wasm build:

- `ETKDGResult` exposes `get_success()` / `get_coordinates()` / `get_n_atoms()` /
  `get_error()` — the site calls `success()` / `coordinates()` / `n_atoms()`.
  Node repro: `EMBED ERROR (site bug): res.success is not a function` ✓.
- `OptimizationResult` exposes only properties (`n_atoms`, `final_energy`,
  `converged`, `iterations`, `message`) + methods (`get_final_energy()`,
  `get_coord(atom, dim)`, …) — **no coordinates accessor at all**. The site
  calls `coordinates()`, `final_energy()`, `steps()`. Node repro:
  `OPTIMIZE ERROR (site bug): res.coordinates is not a function` ✓.
- `MetaDResult` calls in the site all match the bindings (success(), n_atoms(),
  n_frames(), energies(), times_fs(), cv_values(), fes_s(), fes_f(),
  hill_centers(), n_hills()); node repro of the full metad flow with the site's
  exact defaults (caffeine, dihedral 8,2,1,0, 5000 steps, snap 50) passes with
  consistent arrays (101 frames / 72 FES points / 100 hills). The reported
  toFixed crash is therefore either a stale cached webmm_bg.wasm in the
  browser or an empty-array edge case — the site's `updateReadouts` /
  FES-span code assumes non-empty arrays and calls `.toFixed` unconditionally.
- MD button: all calls match the bindings; node repro passes (no fix needed).

## Phase 1 — Fix confirmed API mismatches (site/index.html)
- `runEmbed`: use `get_success()` / `get_coordinates()` / `get_n_atoms()` /
  `get_error()`.
- `runOptimize`: build the flat coordinate array via `get_coord(a, d)` over
  `res.n_atoms`; read `final_energy` / `iterations` as properties; drop the
  nonexistent `steps()` call (use `iterations`).

## Phase 2 — Harden metad/MD readout code against empty/undefined arrays
- `updateReadouts`: bail out when `traj.energies` is missing or `frame` is out
  of range; guard `times`/`cvValues`/`temps` reads with fallbacks.
- `runMetad` FES span: guard `fes_f().length` (show "—" instead of computing
  on an empty array).

## Phase 3 — Docs per workflow
- Overwrite `PLAN.md` (this file).
- Update `CODE_STATUS.md` (Recently Completed).

## Verification
- Node repro script (repro_demo.mjs, temporary) re-run with the site's exact
  call patterns: embed/optimize succeed; metad readouts + FES span print
  without errors.
- `cargo test` 199/199, clippy 0 (no Rust change; site-only).
- Stage `site/index.html` into local `pkg/` and commit; push; Pages deploys.
- Advise the user to hard-refresh the deployed page (bust cached wasm) and
  re-test; if metad still errors, capture the browser console stack.

## Out of scope
- Changing the Rust/WASM API (bindings are the contract; site was stale).
- Cache-busting schemes (GitHub Pages 10-min cache; hard refresh suffices).
- Other buttons (MD works; viewer/playback untouched).
