# Add regression tests + symmetry invariant + expand validation set

**Status:** ✅ COMPLETE

## Phase 1: Regression tests + symmetry invariant (DONE)
- `src/lib.rs` `regression_tests` module: 6 tests pinning the 3 silent bugs
  (OOP cyclic, DFSB skip, sb_type order) + bond-param symmetry + per-term
  breakdown + purine end-to-end.
- `src/prop_tests.rs`: energy_equals_breakdown_sum invariant.
- 172 tests pass (was 165), 0 warnings.

## Phase 2: Validation set expansion (DONE)
- 21 new molecules added (130 total). Atom types 0.00% mismatch on all.
- Bond-param fixes for 3-ring (CR3R), cumulated (C=C=C), sulfonate chemistry.
- Fixed: allene, cyclopropene, p_toluene_sulfonic_acid, ethylene_oxide.

## Metrics
- 130-mol set: r=1.0000, RMSD=0.281, 1 outlier >1.0 (aziridine)
- Original 109-mol set: unchanged (RMSD=0.244, 0 outliers >1.0)
- Known gap: 3-ring torsion V3 (aziridine +1.24, cyclopropane +0.71) —
  tor_type 0 CR3R-CR3R central bond, RDKit V3=0.236. Future torsion-table work.
- 172 tests pass, 0 warnings; WASM rebuilt.
