# Code Status

## Project Summary
WebMM is a WASM-based molecular geometry optimizer using MMFF94/MMFF94s force field and L-BFGS optimization.

## Current Focus
Phase 18: Full RDKit-compatible MMFF94s implementation.
- VDW: Added aTerm^7 damping factor matching RDKit's buffered 14-7 formula
  - Old: E = ε * α^7 * (α^7 - 2) where α = 1.12 / (ρ + 0.12)
  - New: E = ε * a^7 * b where a = 1.07*R*/(R+0.07*R*), b = 1.12*R*^7/(R^7+0.12*R*^7) - 2.0
- Electrostatics: Added buffered distance (r + 0.05) matching RDKit
  - Old: E = 332.06 * q_i * q_j / r²
  - New: E = 332.07 * q_i * q_j / (r + 0.05)
- 1-4 interaction scaling: VDW and electrostatics scaled by 0.75 for torsion-related pairs
  - Added `one_four_pairs` HashSet to MMFFForceField, populated from torsion endpoints
- BCI charges: Re-enabled with corrected charge transfer direction
  - More electronegative atoms (more negative fbci) gain negative charge
  - BCI table updated with correct MMFF94 values (C=O: 0.42, C-O: 0.35, O-H: 0.40, C-H: 0.05)
- SDF parser: Fixed coordinate field width from 9 to 10 chars (V2000 standard)
- Acetic acid test: E2E test verifying optimizer preserves RDKit reference geometry
- All 116 tests passing, WASM rebuilt

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
- Fix BCI charge calculation accuracy (electrostatic energy off for acetic acid)
- Fix remaining test molecule SDFs in site/index.html with RDKit-generated versions
- Deploy to GitHub Pages

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
- BCI charge calculation produces inaccurate charges for some molecules (acetic acid electrostatic energy ~60 kcal/mol off from RDKit)
- RDKit distinguishes subtypes of O_3, C_AR that our simplified typing does not
- CO2 linear molecule test is flaky (ETKDG embedding occasionally produces degenerate geometry)
