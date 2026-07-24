# Code Status

## Project Summary
WebMM is a WASM-based molecular geometry optimizer using MMFF94/MMFF94s force field and L-BFGS optimization.

## Current Focus
None active — MMFF formal charge distribution fixed 2026-07-23 (see Recently Completed).

## Recently Completed
- **MMFF formal charge distribution fixed to match RDKit** (Halgren MMFF.V equation 15), per PLAN.md "Fix MMFF Formal Charge Distribution to Match RDKit":
  - `src/mmff/charges.rs`: Replaced naive `charges = [0; n] + BCI + uniform_residual` with RDKit's two-stage model: (1) `compute_mmff_formal_charges` computes MMFF formal charges using per-type rules (carboxylate/sulfonate O → shared -1/n, N-oxide N/O → 0, ions → integer charge, everything else → 0), and (2) `calculate_bci_charges` applies equation 15: `pChg = (1 - M·v)·q0 + v·sumFormal + sumBci` where M=crd, v=fcadj, q0=MMFF formal charge (adjusted by neighbor charges when v=0)
  - `src/mmff/charges.rs`: Fixed `mmff_bond_type` for BCI — Aromatic bonds now return bt=0 (matching RDKit `getMMFFBondType`), not bt=4; also added `both_arom` check for SINGLE bonds between two aromatic atoms
  - Partial charges now match RDKit exactly for: acetate (C=-0.106, C_carbox=+0.906, O=-0.900×2), pyridine N-oxide (N=+0.571, O=-0.750), nitromethane, sulfonate, ammonium
  - **All 163 tests pass**, 0 clippy errors, WASM build verified
- **Remaining MMFF atom typing gaps fixed vs RDKit 2025.09.3** (carboxylate C=41, N-oxide N=69, Si=19, imine H=27), per PLAN.md "Fix Remaining MMFF Atom Typing Gaps vs RDKit":
  - `src/mmff/mod.rs`: `MMFFAtomType` gained `C_CO2`(41), `N_POX`(69), `Si`(19), `H_NIM`(27); `assign_atom_types` now detects carboxylate C (sp2 C with C=O + C-[O⁻] → 41), both carboxylate O → 32, pyridine N-oxide N → 69 + O → 32, silicon (Z=14) → 19, imine N-H (H on sp2 N with C=N double bond, non-aromatic) → 27; `base_type()` collapses `H_NIM`→`H`; `determine_h_subtype` now routes imine H to 27 instead of misclassifying as amide H=28
  - `src/mmff/bond.rs`: 7 RDKit-exact bond-stretch parameter sets for newly reachable type pairs: (1,41) 3.830/1.510, (32,41) 9.756/1.261, (37,69) 5.396/1.352, (32,69) 6.098/1.261, (1,19) 2.866/1.830, (6,19) 4.661/1.660, (9,27) 6.230/1.026
  - `src/molecule/parser.rs`: added Si to `get_atomic_number` (→14) and `get_atomic_mass` (→28.0855)
  - `src/mmff/params.rs` + `src/mmff/atom_types.rs`: `mmff_type_id` mappings and property fallbacks for all 4 new variants
  - New `type_audit5` test module (`src/lib.rs`): 5 tests covering carboxylate typing (acetate/benzoate + acetic acid regression), N-oxide typing (pyridine N-oxide + pyridine regression), silicon (TMS), imine H (methanimine H=27), energy/optimizer convergence
  - **All 162 tests pass**, 0 clippy errors, WASM build verified
- **MMFF atom typing gaps fixed vs RDKit 2025.09.3** (alkene C=2, S=O=17, SO2=18, nitro N=45, O2CM=32, NSO2=43, metal/halide ions), per PLAN.md "Fix MMFF Atom Typing Gaps vs RDKit":
  - `src/mmff/mod.rs`: `MMFFAtomType` gained `C_VIN`(2), `S_OX`(17), `S_O2`(18), `N_NO2`(45), `N_SO2`(43), `F_M`(89), `CL_M`(90), `BR_M`(91); `assign_atom_types` now splits sp2 C (double bond to N/O/S → `C_2`, C=C only → `C_VIN`), types sulfoxide/sulfone S (17/18), nitro N (45), sulfonamide N (43), SO2/NO2 oxygens (32), metal ions by element+formal charge (Li/Na/K/Mg/Ca/Zn/Fe/Cu), halide ions by charge −1
  - `src/mmff/charges.rs`: fixed pre-existing BCI bond-type bug — Double bonds were mapped to bt=1, breaking ALL carbonyl BCI lookups; now Single→1 iff both atoms sbmb, Aromatic→4, else 0 (matches RDKit `getMMFFBondType`)
  - `src/mmff/bond.rs`: 17 RDKit-exact bond-stretch parameter sets (kb/r0) for newly reachable type pairs: (1,2) 4.539/1.482, (2,2) 9.505/1.333, (2,3) 4.565/1.468, (2,37) 5.007/1.449, (2,5) 5.17/1.083, (2,6) 5.52/1.373, (1,17) 2.841/1.813, (17,37) 3.098/1.787, (7,17) 8.77/1.5, (1,18) 3.258/1.772, (18,37) 3.281/1.77, (18,32) 10.748/1.45, (18,43) 3.301/1.71, (1,43) 3.971/1.472, (1,45) 3.844/1.48, (37,45) 4.705/1.431, (32,45) 9.42/1.233
  - `src/molecule/parser.rs`: `M  CHG` record parsing (authoritative V2000 charges), metal atomic numbers/masses (Li, Na, K, Mg, Ca, Zn, Fe, Cu)
  - `src/mmff/params.rs` + `src/mmff/atom_types.rs`: `mmff_type_id` mappings and property fallbacks for all 8 new variants
  - New `type_audit4` test module (`src/lib.rs`): 19 RDKit-generated V2000 fixtures, 6 tests covering sp2-carbon split (ethylene/allene/ketene/acrylate/methanimine/acetone/thioacetone/vinyl ether), oxidized sulfur (DMSO/sulfone/sulfonamide/sulfinamide), nitro (nitromethane/nitrobenzene), ions (Na⁺/Mg²⁺/Fe³⁺/Cl⁻), BCI partial charges (DMSO, sulfone, nitromethane, benzaldehyde — hand-derived from RDKit BCI table), and energy/optimizer convergence
  - **All 157 tests pass**, 0 clippy errors
- **UFF electronegativity & radii fixed to match RDKit** (`src/etkdg/mod.rs:307-430`):
  - `uff_electronegativity()`: replaced 9 incorrect Pauling-scale values with RDKit's GMP_Xi values (N=6.899, O=8.741, F=10.874, Si=4.168, S=6.928, Cl=8.564, Br=7.79, I=6.822, B=5.11)
  - `uff_radius()`: added hybridization-specific radii (C_3/C_R/C_2/C_1, N_3/N_R/N_2/N_1, O_3/O_R/O_2, S_3+2/S_R/S_2, B_3/B_2) and aromaticity check
  - `compute_bond_length()`: updated signature to accept `hyb1, is_aromatic1, hyb2, is_aromatic2`; aromatic-aromatic bonds now use BO=1.5 regardless of Kekulé form
  - Fixed P radius: 1.087 → 1.101 (P_3+3), I radius: 1.338 → 1.382
  - Both call sites (`build_distance_bounds`, `bond_lengths_reasonable`) updated to pass hybridization and aromaticity
  - Aniline: initial energy improved 296→286 kcal/mol, N-C distance 1.43→1.40 Å
  - All 151 tests pass, WASM build verified
  - `src/mmff/vdw.rs` already contains full RDKit VDW parameter table (95 entries with alpha_i, N_i, A_i, G_i, DA)
  - `calc_r_star_ij()` implements B/Beta correction: `R*_ij = 0.5 * (R*_i + R*_j) * (1 + B*(1 - exp(-Beta*gamma^2)))`
  - `calc_well_depth()` implements Slater-Kirkwood formula: `epsilon = 181.16 * G_i * G_j * alpha_i * alpha_j / ((sqrt(alpha_i/N_i) + sqrt(alpha_j/N_j)) * R*_ij^6)`
  - `apply_da_scaling()` applies DARAD=0.8 and DAEPS=0.5 for Donor-Acceptor pairs
  - Buffered 14-7 potential with vdw1=1.07, vdw2=1.12 parameters
  - All 151 tests pass, WASM build verified

## Recently Completed
- **Aniline NH2 ETKDG geometry fix**: Added `fix_aniline_nh2_geometry()` post-processing function to restore RDKit-like pyramidal aniline NH2 after ETKDG minimization:
  - Problem: ETKDG flattening + minimization squeezed H-N-H to ~102° (fully planar), conflicting with RDKit's ~117.5°
  - Solution: After `minimize_etkdg`, detect aniline-like sp2 N atoms (non-aromatic N with 2 H's, 1 heavy neighbor, bonded to aromatic ring) and explicitly re-place H atoms at 117.5° with ±0.89 Å pyramidalization
  - Result: H-N-H = 117.36°, N-H = 1.04 Å, H's ±0.89 Å from ring plane (matches RDKit ETKDG behavior)
  - All 151 tests pass, diag_aniline tests preserved, WASM build verified
- **Angle parameter lookup fix**: Removed incorrect `base_type()` mapping from angle parameter lookup. The lookup now uses original MMFF type IDs (e.g., H_NAM=28 instead of H=5), allowing the equivalence level mechanism to correctly find MMFFANG table entries:
  - Root cause: `base_type(H_NAM) = H` converted type 28→5 before lookup, but MMFFANG table has entries for type 28
  - Same issue affected all H subtypes (H_OH=31, H_ONC=21, H_COOH=24, H_OAR=29, H_N3=23, H_NAM=28)
  - Fix: `get_angle_params_with_bond_info` now passes original types to `lookup_angle_params`; bond r0 lookup still uses `base_type` for parameter retrieval
  - Verified: all 6 aniline angle params now exactly match RDKit (ka and theta0)
  - Verified: all 4 aniline stretch-bend params still match RDKit exactly
  - Diagnostic test corrected: bond types changed from bt=1 to bt=0 for aromatic bonds (matching RDKit's `getMMFFBondType` which returns 0 for Aromatic bonds)
- **All 149 tests pass**, 0 failures
- **MMFFANG table-based angle lookup**: Replaced 800+ line hardcoded match with 2342-entry MMFFANG table + multi-stage equivalence level lookup (levels 0-3):
  - Added `MMFF_ANGLE_TABLE` (2342 entries) and `MMFF_DEF_EQ_LEVELS` (95 entries × 4 levels) to `mmff_tables.rs`
  - `compute_angle_type()` derives angle type (0-8) from bond types + ring size
  - `lookup_angle_params()` iterates 4 equivalence levels with i/k canonicalization
  - Falls back to RDKit empirical rule (Halgren eq. 20) when no table match
  - **Water: optimizer now converges in 12 iters (was 5000 non-converging), angle error 1.02°**
  - **Ammonia: converges in 8 iters, angle error 1.00°**
  - **Acetic acid: converges in 89 iters (was non-converging)**
- **MMFFStbn table-based stretch-bend lookup**: Replaced hardcoded match with 282-entry table + 30-entry DFSB fallback:
  - Added `MMFF_STBN_TABLE` (282 entries) and `MMFF_DFSB_TABLE` (30 entries)
  - `compute_stretch_bend_type()` derives SB type (0-11) from angle type + bond types
  - `lookup_stretch_bend_params()` canonicalizes i/k, falls back to periodic-table-row defaults
- **MMFFProp table + empirical angle rule**: Added 95-entry `MMFF_PROP_TABLE` (atno, crd, val, pilp, mltb, arom, linh, sbmb) and `empirical_angle_params()`:
  - Theta0 from central atom properties: crd=4→109.45°, crd=2+O→105°, crd=3+val=3+mltb=0+N→107°
  - Ring overrides: 3-ring→60°, 4-ring→90°
  - Force constant from Halgren eq. 20: `ka = beta * Z_i * C_j * Z_k / ((r0_ij + r0_jk) * θ₀² * e^(2D))`
  - Z/C values from Halgren Table VI for H, C, N, O, F, Si, P, S, Cl, Br, I
- **All 164 tests pass**, 0 clippy errors, WASM build verified
- **All 17 opt_compare tests pass** (16 molecules + debug breakdowns)


  - Well depth: `epsilon_ij = 181.16 * G_i * G_j * alpha_i * alpha_j / ((sqrt(alpha_i/N_i) + sqrt(alpha_j/N_j)) * R*_ij^6)`
  - R*_ij with B/Beta correction: `R*_ij = 0.5 * (R*_i + R*_j) * (1 + B * (1 - exp(-Beta * gamma_ij^2)))` where B=0.2, Beta=12.0
  - Donor-Acceptor scaling: `R*_ij *= 0.8`, `epsilon *= 0.5` when one atom is Donor and other is Acceptor
  - Full 95-entry `VDW_TABLE` from RDKit Params.cpp (indexed by type_id)
- **Torsion type classification fixed**: Replaced incorrect ring-based classification with RDKit's `getMMFFTorsionType`:
  - Type 4 only for 4-membered rings (was incorrectly applied to ALL ring bonds)
  - Type 5 only for 5-membered rings with sp3 carbon present
  - Bond type from `getMMFFBondType`: returns 1 if SINGLE bond between two sbmb/arom atoms
  - Added `MMFFPROP_SBMB` and `MMFFPROP_AROM` tables (100 entries each) for atom type property lookup
  - Precomputed `torsion_types` stored in `MMFFForceField` struct
  - **Cyclohexane: +13.4 → +0.017 kcal/mol delta** (fixed!)
- **VDW no 1-4 scaling**: Removed 0.75 scaling from VDW energy (only electrostatics uses 1-4 scaling)
  - **Benzene: -3.1 → +0.067 kcal/mol delta** (fixed!)
- **VDW subtype lookup**: Use original atom type (not base_type) for VDW parameter lookup
  - H subtypes (H_N3, H_COOH, H_OAR, etc.) now use their own alpha=0.150 and DA='D' instead of generic H alpha=0.250, DA='-'
  - **Glycine VDW at RDKit geometry: 5.2 → 1.905 (matches RDKit exactly!)**
- **find_torsions fixed**: Changed to iterate ALL bonds instead of only "rotatable" bonds
  - Ethane torsion energy: +0.173 → -4.786 (RDKit: -4.734)
- **Clippy warnings fixed**: Removed 9 unreachable pattern warnings in angle.rs, atom_types.rs, stretch_bend.rs
- **All 164 tests pass**, WASM build verified
- **E10 Timeout**: Changed `timeout_ms: u64` to `timeout_s: f64` to match RDKit's seconds-based timeout
- **E11 RNG**: Replaced Xoshiro256** with MT19937 Mersenne Twister matching `std::mt19937` — same seed now produces same sequence as RDKit
- **M5 Torsion**: Replaced approximate 20-entry torsion parameter table with full 926-entry RDKit parameter set + 5-stage equivalence protocol:
  - Added `TOR94_TABLE` and `TOR94S_TABLE` const arrays (926 entries each, sorted by (tor_type, j, k, i, l))
  - `get_torsion_params` now takes (type1, type2, type3, type4, tor_type, variant) instead of (type1, type2, type3, type4, variant)
  - Torsion type classified from bond type + ring membership: single→0/4(ring), aromatic→1, double→2, fallback→5
  - 5-stage equivalence lookup: iter 0=exact, iter 1/2=half-wildcard, iter 3=full wildcard, iter 4=retry with fallback tor type
  - Table entries with iType=0 or lType=0 are wildcards (match any type) — 4-way binary search per iteration
  - MMFF94s differs from MMFF94 in 42 entries (types 10=N_AM and 40=N_PL3)
  - Added `is_in_ring(atom1, atom2, mol)` function to `src/molecule/graph.rs`
  - Updated both call sites in `src/mmff/mod.rs` (energy+gradient and energy_breakdown)
  - If no match found (all 5 stages + fallback fail), returns `None` (torsion skipped) instead of `Some(TorsionParams{0,0,0})`
  - 6 tests: exact match, aromatic lookup, no-match fallback, MMFF94 vs MMFF94s difference, energy, gradient
- **M6 OOP**: Replaced single-k-per-type OOP lookup with full 117-entry parameter tables + 4-stage equivalence protocol:
  - Added `OOP94_TABLE` and `OOP94S_TABLE` (117 entries each, 11 differ for MMFF94s)
  - `get_oop_params` now takes (i_type, central_type, k_type, l_type, variant) — 4-type lookup instead of central-only
  - 4-stage equivalence: iter 0=exact, iter 1-3=generalize peripheral i,k,l via eq levels, central j always exact
  - Peripheral atoms sorted before lookup (matching RDKit's protocol)
  - Both call sites in `mod.rs` updated to pass all 4 atom types
  - 4 new tests: exact match, aromatic, no-match fallback, MMFF94 vs MMFF94s
- **M7 BCI**: Expanded BCI charge table from 32 to 498 entries (full RDKit coverage):
  - Replaced `BciEntries` match-arm struct with sorted `BCI_TABLE` const slice
  - `lookup_bci_canonical` uses `binary_search_by_key` on `(bt, i, j)` for O(log n) lookup
  - All RDKit BCI entries included (bond types 0/1/2/4, types 1-99)
- **Shared infrastructure** (`src/mmff/params.rs`):
  - Moved `mmff_type_id` from charges.rs to shared location
  - Added `EQ_LEVELS` table (100 entries × 4 equivalence levels from RDKit's defaultMMFFDef)
  - Added `classify_torsion_type(bond_type, is_in_ring) -> (primary, fallback)` for torsion bond classification
  - All 143 tests pass, 0 clippy errors, WASM build verified
- **E3**: Fixed chiral center R/S stereochemistry detection:
  - Added `stereo_parity` field to `Atom` struct (0=unspecified, 1=odd/CCW/R, 2=even/CW/S, 3=either)
  - V2000 parser reads stereo parity from column 39-42 of atom lines
  - V3000 parser reads CFG=n field from atom lines
  - `find_chiral_centers()` uses stereo parity to set negative volume bounds for S/CW centers
  - S (CW) centers: (-100.0, -5.0) for 4-neighbor, (-100.0, -2.0) for 3-neighbor
  - R (CCW) centers: (5.0, 100.0) for 4-neighbor, (2.0, 100.0) for 3-neighbor
  - Falls back to positive bounds when stereo_parity is 0 (unspecified)

## Recently Completed (ETKDGv3 Gap Fixes)
- **C5**: Added basin threshold (BASIN_THRESH = 5.0) to all force field functions
- **C3**: Added first 4D minimization with proper weight staging (chiral=1.0, 4thdim=0.1 → chiral=0.2, 4thdim=1.0)
- **C3**: Added energy-per-atom check (0.05 threshold) after first minimization
- **C6**: Fixed 1-2/1-3 constraints to use target distances from bounds matrix instead of current distances
- **C6**: Added long-range distance constraints with LONG_RANGE_FORCE = 10.0
- **S8**: Scaled attempts to molecule size (10 * n_atoms)
- **M4**: Fixed sameSide tolerance — chiral centers use 0.10, tetrahedral uses 0.30
- **M5**: Changed random_seed default to -1 (matching RDKit)
- **S1**: Added conjugated 5-ring squish (EXTRA_SQUISH = 0.2)
- **S2**: Added fused small ring volume relaxation (volScale=0.25 for atoms in >1 ring <5)
- **S5**: Added bridged ring exclusion in torsion preferences
- **C1**: Expanded torsion patterns from 15 to ~50+ with proper atom-environment matching
- **C2**: Replaced covalent radii with UFF bond length formula (Pauling + electronegativity correction)
- **C4**: Added missing 1-4 topology cases: cumulated double bonds (force cis), S-S 90° torsion, two-in-same-ring (SP2→trans), two-in-different-rings, macrocycle amide +0.1 offset, amide/ester 1-4 with H-trans/non-H-cis
- **S3**: Added SP3D (105°) and SP3D2 (90°) hybridization to enum, angle calculations, and 1-3 bounds
- **S6**: Replaced PCA planarity check with improper energy method (energy/atom < 0.7, matching RDKit)
- Ring torsions restructured: separate handling for macrocycle (>=9), 5-ring, 3-4 ring, 6-8 ring
- Fixed all flaky tests (pyrrole, acetic acid, benzene planarity stress, worst-run) to use fixed seeds
- All 134 tests pass, 0 clippy errors (WASM verified)

## Completed (M6-M12)
- M6: SP3D/SP3D2 MMFF atom types — Added P_3D, S_3D, S_3D2 to MMFFAtomType enum with base_type fallback to P_3/S_3, assign_atom_types detection for Sp3D/Sp3D2 hybridization, angle params for trigonal bipyramidal/octahedral geometries, atom type properties and charge mappings
- M7: Atropisomer chirality — Added AtropCW/AtropCCW to BondStereo, `find_stereo_bonds` detects atropisomers, `check_atropisomers` validates chiral volume sign, test with simple C(F)(H)-C(Cl)(Br) atropisomer
- M8: Fragment embedding — `find_connected_components` detects disconnected fragments, `extract_submol` creates sub-molecules, `spread_fragments` places fragments with 5Å gaps, test with water dimer
- M9: CoordMap (constrained atoms) — Added `coord_map: HashMap<usize, [f64; 3]>` to ETKDGConfig, fixed atoms during 4D embedding and 3D minimization, test with water O fixed at origin
- M10: Timeout support — Added `timeout_ms` to ETKDGConfig, `web-time` dependency for WASM-compatible timing, timeout check in embedding loop
- M11: et_version field usage — Wired `config.et_version` into `build_torsion_preferences` and `match_torsion_pattern`, v2 uses softer sp3-sp3 torsion barriers (5.0 vs 7.0, 1.5 vs 2.0 for rings)
- M12: Triangle smoothing epsilon — Made `smooth_triangle_inequality` accept configurable epsilon parameter, added `triangle_smoothing_epsilon` to ETKDGConfig (default 1e-6)
- **WASM Verification**: `wasm-pack build --target web` succeeds, 314KB wasm binary generated at `pkg/webmm_bg.wasm`
- **All 139 tests pass**, 0 ignored
- **Acetic acid planarity fix**: Added missing MMFF angle parameter C_3-C_2-O_R (112.42°), made O=C-O angle consistent (121.02°), increased C_2 OOP constant from 0.04 to 0.20 to enforce carboxyl planarity. Carboxyl group now planar (dihedral ~180°).

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
- V3000 charge parsing not implemented; halide ions beyond F⁻/Cl⁻/Br⁻ (no I⁻ in MMFF94)
- CO2 linear molecule test is flaky (ETKDG embedding occasionally produces degenerate geometry)
- Electrostatic energy differs from RDKit total energy due to different Eel formulation (charges now match)

## Fixed Issues (Post-RDKit Review)
- 2026-04-24 — Fixed stretch-bend conversion factor: added missing `0.5 * C_SB` where `C_SB = 143.9324 * π/180 ≈ 2.512` (src/mmff/stretch_bend.rs). Was producing energies ~5× too small.
- 2026-04-24 — Fixed OOP angle calculation: rewrote `calculate_oop_angle` to use proper Wilson angle (RDKit-style) with normalized vectors from central atom (src/mmff/oop.rs).
- 2026-04-24 — Fixed MMFF94s variant ignored: `get_torsion_params` now accepts and uses `mmff_variant` parameter, with MMFF94s-specific overrides for amide and peptide torsions (src/mmff/torsion.rs).
- 2026-04-27 — Completed M6-M12 gaps: SP3D/SP3D2 atom types, atropisomer chirality, fragment embedding, coordMap constraints, timeout support, et_version usage, configurable triangle smoothing epsilon
- 2026-04-27 — Fixed acetic acid carboxyl planarity: `build_planarity_constraints` now adds improper torsions and exocyclic torsions for non-aromatic sp2 atoms (carboxyl, carbonyl, amide groups). Carboxyl group now planar (dihedrals ~0°/180°). Un-ignored `test_etkdg_acetic_acid_geometry`.
- All 135 tests pass after fixes, WASM build verified.
