# Code Status

## Project Summary
WebMM is a WASM-based molecular geometry optimizer using MMFF94/MMFF94s force field and L-BFGS optimization.

## Current Focus
Phase 25: ETKDG v3 Full Exact Port — Phase 7 COMPLETE: Fill All Remaining Gaps

Phase 1 (Distance Bounds) completed:
- `ComputedData` struct, 1-5 geometry helpers, `set_15_bounds_helper`, H-bond detection, macrocycle 14-config

Phase 2 (Embedding Pipeline) completed:
- Fixed `find_chiral_centers` volume bounds, basic chirality detection, `minimize_fourth_dimension`, `check_tetrahedral`

Phase 3 (Torsion Preferences) completed:
- Full 6-term Fourier series, 15 pattern categories, basic knowledge ring torsions

Phase 4 (Force Field) completed:
- Replaced O(n³) finite-difference gradient with full analytical gradients
- Distance bounds: analytical dE/dr for all pairs
- Chiral volumes: analytical gradient via scalar triple product derivatives
- Dihedral gradient: standard MD formula (verified by translational invariance)
- Torsion preferences: analytical dE/dφ chained through dihedral gradient
- Planarity: analytical for ring/exocyclic torsions, numerical for impropers (negligible cost)
- All 122 tests pass, 0 clippy warnings, 4× speedup (19s → 4.4s)

Phase 5 (Aromaticity) completed:
- Rewrote `is_aromatic` with RDKit-style candidate check: all ring atoms must be candidates
- Added `is_aromatic_candidate`, `count_pi_electrons`, `has_multiple_bond_in_ring`, `ring_has_heteroatom`, `estimate_total_neighbors`
- Correctly handles 5-membered heteroaromatics in Kekulé form: furan, thiophene, imidazole, pyrrole
- Distinguishes pyrrole-like N (2e⁻) from pyridine-like N (1e⁻) via neighbor count
- Correctly rejects non-aromatic rings: cyclohexane, cyclohexene, 2,5-dihydrofuran, cyclopentadiene
- Added 3 unit tests: `test_is_aromatic_furan_kekule`, `test_is_aromatic_imidazole_kekule`, `test_is_aromatic_2_5_dihydrofuran`
- All 125 tests pass, 0 clippy errors

Phase 6 (Conformer Validation) completed:
- Added `has_vdw_clash` function: rejects conformers with non-bonded atoms closer than 60% of VDW sum
- Integrated clash check into conformer selection loop in `generate_initial_coords_with_config`
- Conformer must pass planarity, double-bond geometry, chirality, stereo, AND be clash-free to be accepted
- All 125 tests pass, WASM build verified

Phase 7 (Fill All Remaining Gaps) completed:
- Amide/ester trans preference in 1-4 bounds: non-ring 1-4 paths through amide (C=O–N) or ester (C=O–O) bonds now enforce trans geometry via `compute_14_dist_trans` with `GEN_DIST_TOL` window
- Added `is_ester_bond` helper: detects ester single-bond C–O adjacent to carbonyl, excludes anhydrides
- Wired up `force_trans_amides` config: default changed from `false` to `true` (RDKit ETKDGv3 default); bounds builder now reads `config.force_trans_amides` for non-ring 1-4 paths
- Bond length validation: added `bond_lengths_reasonable` function checking all bond distances against expected covalent-radius sums (with bond-type scaling); tolerance ±0.30 Å
- Integrated bond validation into conformer selection loop alongside existing checks
- All 125 tests pass, 0 clippy errors, WASM build verified

## Completed
- Phase 1: Fixed molecule layer
  - Parser: atom symbol from correct columns (31-34), molecule name from correct line, removed fake second-bond parsing, V2000 charge encoding, bond block boundary check
  - Hybridization: now considers bond order (pi bonds) for correct sp1/sp2/sp3 assignment
  - Adjacency: cached in Molecule struct, built during parsing, get_neighbors returns &[usize] from cache
  - Tests: added water parser test and charge encoding test
- Phase 2: Fixed MMFF gradients and energy formulas
  - Bond gradient: fixed double divide-by-r (1/r applied twice), added r < 1e-10 guard, tests for direction/equilibrium/numerical
  - Angle gradient: added missing factor of 2 in dE/dtheta, added sin_theta collinear guard
  - Torsion gradient: fixed wrong imports, replaced completely wrong analytical gradient with correct numerical finite-difference gradient, all 4 atoms now receive gradients
  - OOP gradient: fixed wrong imports, replaced sin(energy) bug with numerical finite-difference gradient
  - VDW energy: fixed beta-term cancellation bug (repulsive/attractive terms swapped), implemented proper buffered 14-7 with attractive well, gradient via numerical differentiation
  - Electrostatics: fixed reversed gradient signs (missing negative sign), removed unused import
  - mod.rs: added H atom type to MMFFAtomType enum, added atomic_number==1 match, fixed torsion gradient destructuring (was discarding g1/g4 and misapplying g2/g3)
  - Bond params: added H-C_3, H-N_3, H-O_3, H-S_3, C_2-C_3, C_2-C_2, C_2-O_2, C_2-O_R, C_2-N_2 (symmetric H matching)
  - Angle params: added C_3-C_3-H, H-C_3-H, C_3-N_3-H, C_3-O_3-H, C_3-C_2-H, H-C_2-H, C_2-C_3-C_3
- Phase 3: Fixed L-BFGS optimizer
  - Rewrote two-loop recursion: correct backward loop (s^T*q, alpha per step, q -= alpha*y), H0 gamma scaling, forward loop (y^T*r, beta, r += (alpha-beta)*s)
  - Fixed iteration count: was `converged as usize` (always 1), now tracks `iter + 1`
  - Added energy change convergence check (prev_energy tracking)
- Phase 4: Fixed MMFF orchestrator
  - Single-pass energy+gradient via `calculate_energy_and_gradient`; `calculate_energy` and `calculate_gradient` delegate to it
  - Cached topology (angles, torsions, oops) computed once in `new()` instead of recomputed every call
  - Added 1-3 pair exclusion for VDW (angle endpoints excluded via `excluded_pairs` HashSet)
  - O(1) VDW exclusion lookup replacing O(N*B) bond scan
  - Safe fallback for unsupported atom types (no more panic)
  - Deduplicated MMFFVariant: mmff/mod.rs re-exports from lib.rs
- Phase 5: Fixed WASM API and ETKDG
  - OptimizationResult now stores and exposes coordinates via `get_coord(atom_idx, coord_idx)`
  - OptimizationOptions has WASM constructor
  - Removed dead InternalOptimizationResult struct
  - MMFFVariant uses derive Default
  - ETKDG: replaced js_sys::Math::random() with getrandom-based random_f64() (works on non-WASM)
  - ETKDG: removed unused Bond import, redundant as usize casts
- Phase 6: Embedded MMFF parameters from JSON
  - `load_mmff_params()` now uses `include_str!` to embed data/mmff94_sample_parameters.json at compile time
  - Added `get_bond_params_from_json()` helper for bond parameter JSON fallback lookup
  - Bond params `get_bond_params` falls back to JSON lookup when no hardcoded match exists
  - Fixed JSON file missing closing brace for root object
- Phase 7: Fixed tests
  - Removed 2 ignored parser tests (test_parse_simple_sdf, test_parse_benzene_sdf) — redundant with existing tests
  - Added end-to-end integration test: parse SDF -> embed coords -> create FF -> optimize
  - Added unit tests for `load_mmff_params()` and `get_bond_params_from_json()`
- Phase 8: Fixed all clippy warnings (58 -> 0)
  - `#[allow(non_camel_case_types)]` on MMFFAtomType enum
  - Renamed `dE_dr` to `d_e_dr` in bond.rs
  - Removed duplicate unreachable patterns in angle.rs
  - Removed duplicate `grad_j` in vdw.rs, removed unused `mut` in test
  - Prefixed unused variables with `_` (mol, atom, iteration, z, atom_types)
  - Fixed unreachable patterns and overlapping ranges in atom type assignment
  - Collapsed nested `if` in graph.rs
  - Used range pattern `5..=7` in parser.rs
  - Used `.clamp()` and `.is_empty()` in etkdg
  - `#[allow(clippy::too_many_arguments)]` on armijo_line_search
  - Removed needless borrows in mmff/mod.rs
  - Cleaned up redundant re-exports in molecule/mod.rs
  - Removed wasm-pack from dev-dependencies (CLI tool, not a crate)
- Phase 9: Updated documentation
  - Fixed pkg/index.html: `result.converged()` -> `result.get_converged()`
- Phase 10: MMFF94 ring detection, aromaticity, atom typing, property table
  - SSSR ring detection via BFS shortest-path tree + canonical deduplication (graph.rs)
  - Aromaticity detection via ring membership + Huckel rule with heteroatom handling (graph.rs)
  - Context-sensitive atom typing: formal charge, neighbor C=O, ether O, amide N (mmff/mod.rs)
  - Atom type property table from Halgren 1996 Table II (mmff/atom_types.rs)
  - 10 new tests: 4 ring detection, 4 aromaticity, 2 property table

- Phase 13: Test suite expansion
  - Angle: 6 new tests (equilibrium energy, straight line, numerical gradient, coincident atoms, linear atoms)
  - Electrostatics: 8 new tests (like charges, magnitude, dielectric, 3D geometry, zero charge, coincident, numerical gradient)
  - Estimation: 7 new tests (double/triple/aromatic bond scaling, heteroatom bonds, linear/trigonal angles, symmetric types)
  - Charges: 3 new tests (ammonia, methane symmetry, single atom)
  - Integration: 11 new tests (single atom, H2, acetylene, dimethyl sulfide, atom type assignment for methane/formaldehyde/hydroxide/ether, optimizer convergence, MMFF94/MMFF94s variants, WASM API)
  - V3000 parser: 2 new tests (water, methane)
  - Property-based tests: 6 proptest invariants (energy finiteness, bond equilibrium, convexity, VDW well, parser robustness, gradient consistency)

## In Progress
- Phase 14: Fix GitHub Pages 3D display
  - Fixed missing JS wrappers for OptimizationOptions setter methods (wasm-bindgen codegen bug)
  - Fixed ETKDGResult.n_atoms -> ETKDGResult.get_n_atoms in site HTML
  - Changed site/index.html to use options.convergence.xxx property setters

## Completed
- Phase 12: BCI charges, parameter estimation, expanded tables, per-type VDW/OOP
  - Task 5: Bond Charge Increment method (src/mmff/charges.rs) — FBCI initialization, 37-entry BCI table, electronegativity-based fallback estimation, charge neutralization
  - Task 6: Parameter estimation rules (src/mmff/estimation.rs) — Halgren eqs 10-13 for bond params (kb, r0 with bond order corrections), eqs 14-16 for angle params (ktheta, theta0 by central geometry)
  - Task 7: Expanded bond parameters from ~15 to 80+ entries: C-C, C-N, C-O, C-S, C-H, N-H, O-H, S-H, halogen, N-N, O-O, P bonds with symmetric matching
  - Task 8: Expanded angle parameters from ~10 to 80+ entries: C-C-C, C=C-C, aromatic, C-N-C, aromatic N, C-O-C, carbonyl, C-S-C, C-P-C, halogen angles with symmetric matching
  - Task 9: Expanded torsion parameters from 3 to 17 central-bond types: C_3-C_3, C_3-C_2, C_2=C_2, C_3-C_AR, C_AR-C_AR, C_3-N_3, C_3-N_PL3, C_3-N_AM, C_AR-N_AR, C_3-O_3, C_3-O_R, C_AR-O_R, C_3-S_3, C_3-P_3, N_3-C_2. Default zero-barrier fallback for unknown torsions.
  - Task 10: Per-type VDW (r0, epsilon, alpha) and OOP (k_oop) from atom type property table — replaced grouped per-element match with individual per-type values
  - Bond/angle parameter lookup now chains: hardcoded table -> estimation fallback -> JSON fallback
  - 6 new tests: 3 charges (ethanol, water, neutralization), 3 estimation (bond, angle, ions)
- Phase 11: ETKDG improvements and integration tests
  - Task 11: Added 1-3 (angle-derived) and 1-4 (torsion-derived) distance bounds + ring closure tightening to `build_distance_bounds`
  - Task 12: Replaced trivial 4D-to-3D projection with eigenvector projection via Gram matrix + power iteration
  - Task 13: Replaced custom DGFF minimization with FF-based refinement using actual MMFF94 + L-BFGS optimizer
  - Task 14: Multi-conformer generation with lowest-energy selection
  - Task 15: 3 new integration tests (ring detection from SDF, benzene embedding, ethanol optimization)
  - Task 16: Updated README.md (removed Partial/TODO, listed all completed features) and CODE_STATUS.md

## Upcoming
- Fix remaining test molecule SDFs in site/index.html with RDKit-generated versions
- Deploy to GitHub Pages

- 2026-04-14 — Replaced BCI charge system with RDKit-compatible model: 31 verified BCI entries, PBCI fallback, directional sign convention
- 2026-04-14 — Fixed H subtype type mapping (H_ONC→21, H_OAR→29, H_OH→31, H_N3→23, H_COOH→24, H_NAM→28)
- 2026-03-28 — Implemented stretch-bend coupling term with 30+ parameter entries from RDKit
- 2026-03-28 — Added H subtype classification (H_OH, H_ONC, H_COOH, H_OAR, H_N3, H_NAM) with `base_type()` normalization
- 2026-03-27 — Fixed MMFF94 angle energy constant from 0.000043945 to 71.9662 (= C_bn/2)
- 2026-03-27 — Calibrated bond and angle parameters against RDKit MMFF94 reference values
- 2026-03-27 — Verified angle formula E = 71.9662 * Z_IJK * dtheta^2 via RDKit numerical differentiation (water gives C_eff = 71.9662 exactly)
- 2026-03-21 — Added adjacency cache to Molecule struct to avoid O(N+E) rebuilds per call
- 2026-03-21 — Used numerical finite-difference gradients for torsion, OOP, and VDW (reliable over complex analytical derivations)
- 2026-03-21 — Swapped VDW repulsive/attractive terms in buffered 14-7 formula (user spec had them inverted)
- 2026-03-21 — Used getrandom with js feature for cross-target random numbers (replaces js_sys::Math::random)
- 2026-03-21 — Cached topology and excluded pairs in MMFFForceField struct (computed once, used everywhere)
- 2026-03-21 — Single-pass energy+gradient in MMFFForceField (halves force field evaluation cost)
- 2026-03-21 — Used include_str! to embed MMFF94 JSON parameters at compile time (works for both native and WASM targets)
- 2026-03-21 — JSON parameter fallback in get_bond_params for types not in hardcoded match table
- 2026-03-22 — Three-tier parameter lookup: hardcoded table -> estimation fallback -> JSON fallback
- 2026-03-22 — Per-type VDW/OOP from atom type property table (replaces grouped per-element values)
- 2026-03-22 — Zero-barrier default for unknown torsions (flat landscape) instead of skipping entirely


## Constraints & Assumptions

## Known Risks / Issues
- RDKit distinguishes subtypes of O_3, C_AR that our simplified typing does not
- CO2 linear molecule test is flaky (ETKDG embedding occasionally produces degenerate geometry)
- Electrostatic energy differs from RDKit total energy due to different Eel formulation (charges now match)
