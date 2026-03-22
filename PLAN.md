# Plan: MMFF94 Tasks 11-16 — ETKDG Improvements & Integration Tests

## Task 11: ETKDG 1-3/1-4 Distance Bounds
- File: src/etkdg/mod.rs `build_distance_bounds`
- Add angle-derived (1-3) bounds using law of cosines
- Add torsion-derived (1-4) bounds using trans/eclipsed extremes
- Add ring closure tightening for upper bounds

## Task 12: ETKDG Proper 4D-to-3D Projection
- File: src/etkdg/mod.rs `project_to_3d`
- Replace trivial drop-4th-coordinate with eigenvector projection via Gram matrix + power iteration

## Task 13: ETKDG FF-Based Refinement
- File: src/etkdg/mod.rs
- Replace `dgff_minimization` with `refine_with_ff` using actual MMFF94 + L-BFGS optimizer

## Task 14: ETKDG Multi-Conformer
- File: src/etkdg/mod.rs `generate_initial_coords_with_config`
- Generate multiple conformers, return lowest-energy one

## Task 15: Integration Tests
- File: src/lib.rs tests module
- Add ring detection test for benzene from SDF
- Add charges non-zero test for water
- Add full optimization test for ethanol

## Task 16: Update Documentation
- Files: README.md, CODE_STATUS.md
- Update completed features list, remove partial/TODO items
