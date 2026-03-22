# Code Status

## Project Summary
MolGopt is a WASM-based molecular geometry optimizer using MMFF94/MMFF94s force field and L-BFGS optimization.

## Current Focus
Phase 11: MMFF94 Tasks 11-16 complete. 34 tests pass (3 new), 0 clippy warnings.

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

## In Progress
- (none)

## Completed
- Phase 11: ETKDG improvements and integration tests
  - Task 11: Added 1-3 (angle-derived) and 1-4 (torsion-derived) distance bounds + ring closure tightening to `build_distance_bounds`
  - Task 12: Replaced trivial 4D-to-3D projection with eigenvector projection via Gram matrix + power iteration
  - Task 13: Replaced custom DGFF minimization with FF-based refinement using actual MMFF94 + L-BFGS optimizer
  - Task 14: Multi-conformer generation with lowest-energy selection
  - Task 15: 3 new integration tests (ring detection from SDF, benzene embedding, ethanol optimization)
  - Task 16: Updated README.md (removed Partial/TODO, listed all completed features) and CODE_STATUS.md

## Upcoming
- (none)

## Technical Decisions
- 2026-03-21 — Added adjacency cache to Molecule struct to avoid O(N+E) rebuilds per call
- 2026-03-21 — Used numerical finite-difference gradients for torsion, OOP, and VDW (reliable over complex analytical derivations)
- 2026-03-21 — Swapped VDW repulsive/attractive terms in buffered 14-7 formula (user spec had them inverted)
- 2026-03-21 — Used getrandom with js feature for cross-target random numbers (replaces js_sys::Math::random)
- 2026-03-21 — Cached topology and excluded pairs in MMFFForceField struct (computed once, used everywhere)
- 2026-03-21 — Single-pass energy+gradient in MMFFForceField (halves force field evaluation cost)
- 2026-03-21 — Used include_str! to embed MMFF94 JSON parameters at compile time (works for both native and WASM targets)
- 2026-03-21 — JSON parameter fallback in get_bond_params for types not in hardcoded match table

## Constraints & Assumptions
- V2000 MOL format only (no V3000 support)

## Known Risks / Issues
- (none)
