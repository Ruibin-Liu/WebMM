# Plan: M2 — full Phase 2+3+4 bundle (RDKit-faithful embedding path)

## Goal
Test the bundle thesis: when bounds + minimizer + 3D-FF change *together* (and the
co-adapted workarounds drop), does the pipeline re-tune and lift r past 0.8603?
(3/3 individual changes regressed; only the bundle can.)

## Architecture (no function duplication)
- `static EXP_RDKIT_ALL: AtomicBool` (default false).
- Extract the embed body to `embed_impl(mol, config)` (shared).
- `generate_initial_coords_default`: `EXP_RDKIT_ALL=false`; `embed_impl`.
- `generate_initial_coords_rdkit`: `EXP_RDKIT_ALL=true`; `embed_impl`.
- Each concern is gated on `EXP_RDKIT_ALL` inside the shared functions:
  - **Bounds:** B1 (smoothing formula `max(L_ik−U_kj,…)` + single sweep),
    D2 (`check_and_set` widen), D3 (sp² distribute/flat-120), D14 (`in_ring` snapshot).
  - **Minimizer:** D7 (4D L-BFGS via `lbfgs_minimize` helper).
  - **3D FF:** D4 (long-range: all non-(1-2/1-3/1-4) pairs, no filter, no
    double-count), D5 (pin 1-2/1-3 to current dist ± tol).
  - **Workarounds:** torsion-snap / bond-snap / H-trilateration / H-relax SKIPPED.
- Default path: all off → byte-identical r=0.8603.

## Deferred (next turns)
- D6 (UFF `InversionContrib` Fourier planarity) — keep current harmonic planarity
  for the first bundle measurement; swap if the bundle shows legs.
- D8 (`computeInitialCoords` retry), D9 (tetra/chiral retry gates).

## Sequence (measure cumulatively; default must stay 0.8603)
1. Toggle + `embed_impl` extraction + `_default`/`_rdkit` setters. Verify 0.8603 both.
2. + B1/D2/D3/D14 (bounds). Measure rdkit path.
3. + D7 (4D L-BFGS). Measure.
4. + D4/D5 (3D FF). Measure.
5. + workarounds-off. Measure = **full bundle (minus D6/retry)**.
6. If r > 0.8603 → add D6/retry, iterate. If r < 0.8603 → ablate (granular toggles)
   to find the culprit; report whether the bundle thesis holds.

## Gate
`cargo test` 191; my code clippy-clean; default r=0.8603 unchanged; deterministic
(run1==run2). Measure rdkit path via `WEBMM_RDKIT_FAITHFUL=1`.
