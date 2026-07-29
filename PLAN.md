# Plan: Fix 3 confirmed ETKDG-vs-RDKit discrepancies (B1, B3, D1)

## Final Status: DONE — shipped B3 + D1; reverted B1 (regresses; pipeline co-adapted)

## Context
Systematic review of `src/etkdg/mod.rs` vs RDKit 2025.09.3 source
(`Code/DistGeom/TriangleSmooth.cpp`, `Code/DistGeom/DistGeomUtils.cpp`,
`Code/GraphMol/DistGeomHelpers/BoundsMatrixBuilder.cpp` + `Embedder.cpp`).
`ETKDG_MMFF_REVIEW.md` marks all 3 issues "FIXED" but the code shows otherwise.
Baseline (current HEAD, 130 mols, seed 42): r=0.8565, RMSD=24.82.

## Completed
- [x] **B3. `is_larger_sp2_atom` ring check** (`src/etkdg/mod.rs:623`)
  - Added the `numAtomRings > 0` condition RDKit's `isLargerSP2Atom` requires
    (added `rings: &[Vec<usize>]` param + 2 call sites in `build_distance_bounds`).
  - SHIPPED — correct, net neutral-to-positive (r 0.8565→0.8568).

- [x] **D1. `set_ring_angle` SP3D2** (`src/etkdg/mod.rs:487`)
  - 135° → 90° to match RDKit `_setRingAngle` (octahedral cis ligands).
  - SHIPPED — correct; no-op on the 130-mol val set (no Sp3D2 ring atoms).

## Reverted (with reason)
- [x] **B1. Triangle-smoothing lower-bound formula** (`src/etkdg/mod.rs:206`) — REVERTED
  - The formula `|lower[ik] − lower[kj]|` is genuinely wrong vs RDKit's
    `max(L_ik−U_kj, L_kj−U_ik)` and can yield infeasible lower bounds.
  - Fixed the formula, but it REGRESSED the harness in every variant:
    formula-only → r=0.8520; formula + single Floyd–Warshall sweep (matching
    RDKit's `triangleSmoothBounds`, removing the `while changed` fixed point)
    → r=0.8498. Dominant mover: cyclopentene +43.7 (a known hard 5-ring).
  - Root cause: WebMM's downstream pipeline (L-BFGS + torsion-snap + H-trilateration
    + H-only relaxation) is co-adapted to the old (over-loose) bound matrix.
    RDKit-faithful smoothing in isolation makes the embedding worse. B1 needs a
    bundled pipeline re-tuning effort, not a one-line fix. Reverted to original.

## Validation
- `cargo build --release --example dump_etkdg_geom` → regenerate
  `scripts/etkdg_val/webmm_ref.json` → `scripts/validate_etkdg.py`.
- r 0.8565→0.8568, RMSD 24.82→24.77; 129/130 embed OK; RDKit ref unchanged.
- `cargo test`: 191 passed / 0 failed; `cargo clippy --all-targets`: 0 warnings.

## Non-goals (deferred, tracked in review)
B2 (MT19937 RNG — only affects bit-exact RDKit reproducibility), D2–D12
(check_and_set widen-vs-tighten, sp2 123/114 heuristic, 3D long-range
filter/double-count, 4D BFGS, full torsion SMARTS table, macrocycle handling),
and a bundled B1+pipeline re-tuning.
