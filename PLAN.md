# Plan: README refresh — align with 0.6.0 codebase

## Context
README.md is stale: last touched in the GitHub Pages commit (d3b43db), predating
the MD engine, metadynamics, GBSA implicit solvation, the MMFF 230/230 RDKit
validation, the site/ + gitignored pkg/ layout, and the pinned wasm-pack
devDependency. It also contains a broken rustup URL and recommends
`cargo install wasm-pack`, which AGENTS.md forbids. Docs-only pass — no code
changes.

## Phase 1 — README.md refresh (keep existing structure/format)
- **Intro**: mention MMFF94/MMFF94s + ETKDG v3 + L-BFGS + MD + metadynamics + GBSA.
- **Project structure**: reflect the actual tree — `src/forces.rs` (ForceField
  trait), `src/md/`, `src/metad/`, `src/solvation/`, full `src/mmff/` module list,
  `site/index.html` (committed demo), `pkg/` (gitignored build output),
  `scripts/` (RDKit validation/benchmark tooling), `examples/`.
- **Progress**: add MD (velocity-Verlet/BAOAB Langevin), metadynamics
  (well-tempered, FES), GBSA (OBC2 + LCPO SA), MMFF validation 230/230 vs RDKit;
  update test count 86 → 199.
- **Build instructions**: fix the rustup URL typo; replace `cargo install
  wasm-pack` with the pinned devDependency (`npm install` + `npm run build`,
  matching package.json); keep the manual wasm-bindgen option with a version-
  match note; note `wasm32-simd128` in `.cargo/config.toml` (scoped to wasm
  target, native unaffected).
- **Demo usage**: add the ⚙ Optimize / 🔥 MD / 🔬 Metadynamics / 🌐 Embed 3D
  buttons and trajectory playback.
- **WASM API**: document `init()`, `generate_initial_coordinates_wasm`,
  `optimize_from_sdf` (+ `_direct`), `run_md_from_sdf`, `run_metadynamics_from_sdf`
  with their options/results (from src/lib.rs).
- **Validation**: brief section pointing at `scripts/benchmark_mmff.py`
  (230/230 gate) and `scripts/validate_etkdg.py` (multi-seed harness).
- **Resources**: add links to `docs/` (atom-type-coverage.md,
  validation-energy-analysis.md).

## Phase 2 — Docs per workflow
- Overwrite `PLAN.md` (this file).
- Update `CODE_STATUS.md`: Recently Completed entry.

## Verification
- README renders (markdown structure intact); every path/command in the README
  matches the repo (verified against Cargo.toml, package.json, src tree,
  site/index.html button ids, lib.rs exports).
- `cargo test` still 199/199 (no code change, sanity only).
- Commit the pass.

## Out of scope
- License section (TBD — not decided).
- Rewriting the historical Progress list into a changelog.
- Updating `docs/` contents or `site/index.html` copy.
