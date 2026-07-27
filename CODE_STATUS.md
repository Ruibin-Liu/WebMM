# Code Status

## Project Summary
WebMM is a WASM-based molecular geometry optimizer using MMFF94/MMFF94s force field and L-BFGS optimization.

## Current Focus
None active — atom typing matches RDKit 100% (0.00% mismatch) on a **130-molecule** validation set. Energy accuracy r=1.0000 / RMSD=0.170 kcal/mol (**0 outliers >1.0**).

## Recently Completed
- **Verified equivalence vs RDKit 2026.03.4 + scoped the one 2026 ETKDG change** (per request "re-run validation harness against 2026.03.x"):
  - **Harness re-run** (`gen_validation_set.py` SMILES→ETKDGv3→MMFF94s, 109 mols / 1324 atoms): atom types **0.00% mismatch**, charges mean |Δ| **0.0000**, energy Pearson **r=1.0000**, RMSD **0.10 kcal/mol** (worst purine +0.59). Geometry-matched RDKit-vs-RDKit (2025.09.3 vs 2026.03.4 on identical SDFs, 130 mols / 1550 atoms): **0 differences** in types/charges/energy — MMFF94s is byte-identical between versions.
  - **No MMFF changes** in any 2026.03.x release (grep of full release notes + all patch bodies = 0 `mmff`/`forcefield` hits). Our 2025.09.3 calibration carries over unchanged.
  - **The only ETKDGv3 code change** is #9143 / PR #9195 ("ETKDGv3 generating 90° twisted amides" in macrocycles). Corrected characterization: the +373/−369 was mostly a tabs→spaces reformat of `torsionPreferences_macrocycles.in`; the **functional change is ~12 SMARTS patterns** whose amide-N ring constraint was relaxed from `;r{9-}` (must be in ≥9-ring) → `;r` (any ring), so the planarizing V2=8.0 torsion preference fires for ring amides even when the N's ring-size match failed. (The 3 `[C:1][C;r{9-}:2](=O)...[NX3;r:3]` "experimental ring amide" patterns at the top exist in *both* versions — not the fix.)
  - **Gap scoping (WebMM vs #9143) — EMPIRICALLY CONFIRMED WebMM has a (worse) macrocycle-amide issue:** WebMM does NOT port RDKit's SMARTS torsion-preference table; its `match_torsion_pattern` (`src/etkdg/mod.rs`) does carry an `is_amide_bond`→V2=8.0 planarizing term for ≥9-ring central bonds, which by code reading *looked* sufficient. But running WebMM ETKDG on the actual #9143 reporter macrocycle (`O=C1NCC=CC[C@@H](NC(=O)c2c[nH]cc2Br)CC(=O)N2CCCC2C2CCN(C2)C(=O)c2csc1n2`, 66 atoms) across seeds 0xf00d/1/7/42/100/12345 produces **multiple ~90–127° twisted amides on every seed** (e.g. seed 0xf00d: C1-N2=127°, C18-N20=124°, C30-N28=105°; only C9-N8 planar). For comparison on the same molecule/seeds: RDKit **2025.09.3** (pre-fix) twists only the C30-N28 amide intermittently (3/6 seeds, 60–64°); RDKit **2026.03.4** (fixed) is planar on all amides/all seeds (C30-N28 4–20°). So WebMM's macrocycle-amide handling is **more broken than RDKit's pre-#9143 state**, not equivalent to the fixed state. Diagnostic: `examples/macrocycle_amide.rs` (reads `/tmp/macrocycle.sdf`). This is a real gap, scoped to macrocycles — does not affect the 130-mol validation set (no macrocycles).
  - **Root cause CONFIRMED (diagnosis)**: not SSSR detection — `find_rings` *does* find the macrocycle (returns rings of size 18/19), and `match_torsion_pattern`'s macrocycle branch *would* return the V2=8.0 planarizing term for the in-ring amide bonds. The term is never applied because `build_torsion_preferences` (`src/etkdg/mod.rs`) skips any central bond with `num_ring_bonds > 3` — and each macrocycle amide C-N bond is counted in **4** rings (find_rings returns the macrocycle as 4 overlapping rings: 18,18,19,19), so `num_ring_bonds=4 > 3` → `continue` → `match_torsion_pattern` is never called for it → no V2=8.0 term. Confirmed via temporary instrumentation: amide bonds 1-2, 18-20, 28-30 all logged SKIPPED (num_ring_bonds=4); these are exactly the 3 amides that came out twisted. The exocyclic amide 9-8 is not skipped and is planar. **Fix implemented (PARTIAL)**: changed `num_ring_bonds` in `build_torsion_preferences` to count only small rings (`r.len() < MIN_MACROCYCLE_SIZE`). Verified by instrumentation: the 3 in-macrocycle amide bonds (1-2, 18-20, 28-30) are no longer skipped and `match_torsion_pattern` now returns the V2=8.0 planarizer `[0,8,0,0,0,0]` for each (exocyclic amide 9-8 unaffected, stays planar). **However this is necessary-but-not-sufficient**: the macrocycle embedding STILL produces twisted amides (127°/124°/105°), geometry byte-identical to pre-fix.
  - **Root cause #2 — Phase 2 FIXED (gates), Phase 5 remains (embedding quality):** `generate_initial_coords_with_config` had three hard pre-refinement gates; every embedding attempt failed GATE2 (`check_tetrahedral`) on the `[C@@H]` stereocenter (atom 7, maximally flat `min_tp=0.0000`), aborting before `minimize_etkdg`. **Fixed (kept):** GATE2/GATE3 are now advisory (no `continue`) — `minimize_etkdg` now runs; the existing post-refinement acceptance gate (planarity + `check_chiral_centers` + `center_in_volume_tol` + stereo + clashes + bond-lengths + energy) judges the result. `check_tetrahedral`/`center_in_volume` retained `#[allow(dead_code)]`. **190 tests pass; 130-mol benchmark unchanged (0% type mismatch, r=1.0000)** — Phase 2 is ETKDG-only. **BUT amides still twist:** Phase 3 experiment (minimize_etkdg 300→3000 iters, step 0.1→0.3) showed some attempts now converge yet land in E≈1000–2000 local minima (torsion-pref energy TP≈540–756) — refinement budget/step is NOT the problem. **The real blocker is Phase 5: the 4D distance-geometry embedding produces a catastrophically bad starting point (E₀≈1575–5419 kcal/mol) for this 66-atom macrocycle; 3D refinement can't escape that basin.** Needs core-DG work (metric-matrix embedding / distance-bound tightening for large rings / `find_rings` over-enumeration fix) — substantial, separate effort. Also removed 2 shadowed duplicate match arms in `src/mmff/bond.rs` (`N_PL3-C_VIN`, `C_2-Cl`) that caused `unreachable pattern` warnings. `cargo build`+`cargo clippy --all-targets`: 0 warnings. Scoped to macrocycles; validation set unaffected.
  - **Fresh-plan execution (Phase 5 investigation):** instrumented the embedder — the 4D metric embedding is NOT broken (top-4 eigenvalues strongly positive 234–739 on all attempts; the `1.0-2.0*rng()` negative-eigenvalue fallback is not the cause). E₀≈5000 is **distributed** across long-range bounds (~782, ~245 violated pairs of ~2145), torsion prefs (~560–730), and 1-2/1-3/chiral/planarity — a globally poor 4D→3D projection for a strained 66-atom ring. Tried and reverted: (a) `minimize_etkdg` 300→3000 iters + step 0.1→0.3 — lands in E≈1000–2000 local minima; (b) embed from smoothed upper-bound matrix (`t=hi`) — E0 unchanged, amides still twisted, and killed per-attempt randomness. **None of the accessible levers work** — the macrocycle embedding is a fundamental DG-quality problem needing deeper work (4D minimization/collapse pipeline, proper SSSR `find_rings`, or a non-classical-MDS embedding strategy), plus real ETKDG regression coverage (the 130-mol benchmark uses RDKit geometry, not WebMM ETKDG, so it can't catch ETKDG regressions). All instrumentation reverted; kept state unchanged (0 warnings, 190 tests, val 0%/r=1.0000).
  - **Incidental (FIXED)**: removed the shadowed duplicate `(N_PL3, C_VIN, Single)` match arm in `src/mmff/bond.rs` (5.952/1.376, dead code) that caused the `unreachable pattern` build warning. Live arm (6.110/1.370) unchanged. Plus #8807 (Python `maxAttempts` kwarg) and #8977 (embed-params JSON) are API/serialization, not algorithm — N/A to the Rust reimplementation.
  - Setup note: ran RDKit 2026.03.4 via `PYTHONPATH` shadow of the PyPI wheel (uv's PyPI simple-index access was blocked in this env; `curl` worked). No repo files changed; `scripts/val_set/` restored (git clean).
- **Audited + fixed the worst remaining |Δ| cases** (per request "check the other cases with big |Δ|"):
  - `src/mmff/bond.rs`: N_4-H (ammonium) 1.012/5.5→1.028/6.163; C_AR-C_2 Single 1.484/5.0→1.457/4.488; added C_2-C_2 Single 1.489/4.418 (oxalic acid); added O_3/O_R-H_COOH 0.981/7.403 (carboxylic acid O-H); C-F 1.38/6.0→1.360/6.011 (sign-flipped stretch-bend, same class of bug as P-C/C_AR-N_AR).
  - All RDKit-verbose verified. Fixed ammonium (+0.99→~0), trifluoromethane (−0.94→fixed), oxalic_acid/benzoic_acid/acetic_acid (carboxylic family).
  - **Result**: 130-set r=1.0000, RMSD 0.248→**0.170**, max|Δ| 0.99→0.95, outliers >0.5: 9→4. 0 outliers >1.0. 190 tests pass, 0 warnings. WASM rebuilt.
  - **Remaining worst**: caffeine (+0.95) — sb gap spread across ~14 amide/carbonyl angles (O-C-N type 7,3,10); all bond params verified correct, gap is in stretch-bend kba table entries for amide angles (needs systematic STBN audit). carbon_disulfide (+0.61, C=S), benzoate (+0.61), cyclobutane (−0.59, 4-ring) — all <1.0, diminishing returns.

## Recently Completed
- **All clippy lints fixed → 0 warnings** (per PLAN.md "Fix all clippy lints"):
  - `cargo clippy --all-targets`: **70 warnings → 0**; `cargo build`/release clean; 190 tests pass; WASM rebuilt.
  - **Mechanical-equivalence fixes** (clippy-verified, no behavior change): `manual_rotate`→`rotate_left` (etkdg Xoshiro tempering); `iter_cloned_collect`→`to_vec`; `len_zero`/`length_comparison_to_one`→`is_empty`; `manual_range_contains`×8→`(lo..=hi).contains`; `nonminimal_bool` simplifications (graph.rs ring test, etkdg); `mul_by_neg_one`→unary `-` (electrostatics); `unnecessary_to_vec`×4 (opt_compare; `optimize` takes `&[[f64;3]]`); `clamp_like_pattern`→`.clamp` (torsion cos_phi); OR-pattern→range (charges halides); `field_assignment outside Default`→struct-update (lib.rs config); `assert_eq!` literal bool→`assert!`; `function call inside expect`×2→`unwrap_or_else`; removed duplicate `#[cfg(test)]` (lib.rs); `needless_range_loop`→`iter().enumerate()`/`zip()` on clean cases (oop, torsion invariance, lib.rs print loops/RMSD/worst_coords).
  - **Scoped `#[allow]`** (structural/pervasive, no refactor): etkdg/mod.rs module allow for `needless_range_loop`+`too_many_arguments`+`type_complexity` (distance-geometry matrix code); extended mmff_tables allow with `large_const_arrays`+`type_complexity`, torsion with `type_complexity` (lookup tables); `too_many_arguments` on `get_mmff_torsion_type`/`get_stretch_bend_params`; `unnecessary_cast` on `estimate_implicit_h` (removing the `as i32` casts triggers a *different* false positive — the `.max(0)` hypervalent-atom guard is load-bearing).
  - `cargo clippy --fix` was unusable: it is transactional and its `unnecessary_cast` suggestion (graph.rs) is a false positive causing E0308, rolling back the whole batch each run; all fixes applied manually with the 190-test suite as the regression net.
- **Fixed aziridine 3-ring torsion bug + regression test** (per PLAN.md "Fix aziridine 3-ring torsion gap"):
  - `src/molecule/graph.rs` (`find_torsions`): was creating degenerate torsions where atom1 == atom4 in 3-membered rings (the ring wraps around to the same atom, e.g. C2-C0-C1-C2). Added `l != k` filter. RDKit excludes these; WebMM didn't, inflating torsion energy (cyclopropane +0.71, aziridine +1.24).
  - `src/lib.rs` (`regression_tests`): `no_degenerate_torsions_in_3_ring` test asserts cyclopropane has 0 torsions with atom1==atom4 and all torsion atoms are distinct.
  - **Result**: cyclopropane now exact (Δ=0.00), aziridine +0.41 (was +1.24). 130-set: r=1.0000, RMSD 0.281→**0.248**, 0 outliers >1.0 (was 1). 190 tests pass (was 172), 0 warnings. WASM rebuilt.

## Recently Completed
- **Test suite hardened: print-only tests converted to real assertions** (per PLAN.md "Harden the test suite: convert print-only tests to real assertions"):
  - **Discovered `src/opt_compare.rs` was orphaned dead code** — never declared as a module, so its 17 tests never compiled or ran. Wired it in via `#[cfg(test)] mod opt_compare;` in `src/lib.rs` and stripped the file's self-wrapping `#[cfg(test)] mod opt_compare { ... }` so it is a normal module body.
  - `src/opt_compare.rs`: added `assert!(result.converged)` + `assert!((webmm_e - rdkit_e).abs() < 0.1)` to all 16 `opt_*` tests and 4 `debug_breakdowns` molecule blocks — previously pure `println!` diagnostics (0 asserts in the whole file). Observed max |Δ| = 0.019 kcal/mol (acetic_acid); most exact.
  - `src/lib.rs`: added asserts to 9 more print-only tests — `test_rdkit_energy_comparison` (|ΔE|<0.01, was ratio-print only), `test_aniline_sdf_compare` (total<0.01 vs RDKit 8.145), `test_diagnostic_ethane_etkdg_steps` (finite energy + C-C∈[1.45,1.55]), `test_worst_run` (benzene planarity <0.01 Å over 10 seeds), `check_webmm_etkdg_aniline` (H-N-H∈[110,125]° + NH₂ pyramidalization >0.3 Å), and the 4 type audits (`test_aniline_mmff_types`/`test_5ring_mmff_types`/`test_simple_molecule_types`/`test_thiophene_imidazole_types`) now assert exact atom MMFF types via `format!("{:?}")` (typing is 0.00% vs RDKit → current values encoded as regression baseline).
  - Renamed 4 duplicate test names for clean `cargo test` filtering: `type_audit4/5::test_new_type_energy_and_convergence` → `test_sp2_sulfur_nitro_energy_and_convergence` / `test_carboxylate_noxide_silicon_imine_energy_and_convergence`; `estimation/atom_types::test_ion_types_return_none` → `test_ion_bond_angle_estimation_returns_none` / `test_ion_atom_type_props_return_none`.
  - **189 tests pass** (was 172; +17 from newly-compiled opt_compare), 0 failures, 0 build warnings. All changes test-gated (`#[cfg(test)]`); no production code changed, WASM unaffected.
- **Regression tests + symmetry invariant + validation-set expansion** (per PLAN.md "Add regression tests + symmetry invariant + expand validation set"):
  - `src/lib.rs` (new `regression_tests` module): 6 tests pinning the 3 silent algorithmic bugs — OOP cyclic permutations (guanidinium→3 terms), DFSB not-skipped (P-O-H→Some), sb_type order-invariance (purine N-C-C), bond-param symmetry (21 curated heteroatom pairs), pyrimidine per-term breakdown vs RDKit, purine end-to-end tolerance.
  - `src/prop_tests.rs`: `energy_equals_breakdown_sum` invariant (catches dropped/double-counted terms).
  - `scripts/add_val_molecules.py` + 21 new molecules in `scripts/val_set/` (formamidine, acetamidine, biguanide, phosphoric/sulfuric/p-toluenesulfonic acid, benzimidazole, benzothiazole, xanthine, CF3CH/CH2Br2/CH2I2/vinyl_chloride/acetyl_chloride, tetramethylsilane, cyclopropene/allene/ethylene_oxide/aziridine, glycine_zwitterion, ibuprofen). Atom types 0.00% mismatch on all new molecules.
  - `src/mmff/bond.rs`: 3-ring (CR3R-CR3R/CR3R-O/CR3R-H/C_2-CR3R/C_VIN-CR3R) + cumulated (C_2=C_1/C_VIN=C_1 Double) + sulfonate (S_O2-O_3) bond params from RDKit verbose. Fixed allene, cyclopropene, p_toluene_sulfonic_acid, ethylene_oxide.
  - **Result**: 172 tests pass (was 165), 0 warnings. Original 109 set unchanged (RMSD 0.244, 0 outliers >1.0). Expanded 130-set: r=1.0000, RMSD 0.281, 1 outlier >1.0.
  - **Known gap**: aziridine (+1.24) + cyclopropane (+0.71) expose a 3-ring torsion V3 gap (tor_type 0, CR3R-CR3R central bond; RDKit V3=0.236). Needs torsion-table investigation.

## Recently Completed
- **All 7 energy outliers fixed → r=1.0000, RMSD=0.238 kcal/mol, 0 outliers >1.0** (per PLAN.md "Fix remaining 7 energy outliers"):
  - `src/mmff/bond.rs`: S_3-S_3 (added 2.050/2.531), C_3-S_3 (1.805/2.893), C-Cl (1.773/2.974), N_AR-N_AR Aromatic (1.246/5.002), C_3-P_3 (1.830/2.790), C_3-P_4 (1.810/2.980), P_4-O_CO2 Double (1.510/8.296), C5A-OFUR (1.360/5.787), NPYL-N5A (1.339/5.513), N5A-C5B (1.335/8.258), C_AR-I (2.075/1.781), C_2=N_2 Double (1.290/10.077), C_2-N_2 Single (1.360/6.385) — all from RDKit verbose.
  - `src/mmff/mmff_tables.rs` (`lookup_stretch_bend_params`): **DFSB row canonicalization bug** — the default stretch-bend lookup used raw row_i/row_k instead of canonical (min,max). For asymmetric angles like P-O-H (rows 2,1,0), the table stores (0,1,2) but the query checked (2,1,0) → no match → stretch-bend skipped entirely. Fixed to canonicalize. Root cause for trimethylphosphine + methylphosphonic_acid.
  - `src/mmff/stretch_bend.rs` (`get_stretch_bend_params`): **sb_type canonicalization** — RDKit computes the stretch-bend type using bond_type_1 = bond to the LOWER-type peripheral atom (canonical i where type_i ≤ type_k). WebMM passed raw bond order, so angles like N-C-C (9,3,2) got sb_type 2 instead of 1, using DFSB default kba (0.3) instead of specific STBN entries (0.61/0.227). Fixed by swapping bt_ij/bt_jk when ti > tk. Root cause for purine.
  - **Validation (109 mols)**: energy Pearson r→**1.0000**, RMSD 0.85→**0.238 kcal/mol**, max|Δ| 3.87→**0.99**, outliers >1.0: 7→**0**. Atom types 0.00%; charges mean|Δ|=0.0003. 165 tests pass, 0 warnings. WASM rebuilt.
  - Fixed all original outliers: guanidinine, pyrimidine, dioxane, acetonitrile, propyne, thiophene, thiazole, dimethyl_disulfide, CCl4, pyridazine, furan, pyrazole, trimethylphosphine, methylphosphonic_acid, iodobenzene, purine. Worst remaining: ammonium (+0.99), caffeine (+0.95) — both <1.0.
  - `src/mmff/bond.rs`: S_3-S_3 (added 2.050/2.531), C_3-S_3 (1.805/2.893), C-Cl (1.773/2.974), N_AR-N_AR Aromatic (1.246/5.002), C_3-P_3 (1.830/2.790), C_3-P_4 (1.810/2.980), P_4-O_CO2 Double (1.510/8.296), C5A-OFUR (1.360/5.787), NPYL-N5A (1.339/5.513), N5A-C5B (1.335/8.258), C_AR-I (2.075/1.781), C_2=N_2 Double (1.290/10.077), C_2-N_2 Single (1.360/6.385) — all from RDKit verbose.
  - `src/mmff/mmff_tables.rs` (`lookup_stretch_bend_params`): **DFSB row canonicalization bug** — the default stretch-bend lookup used raw row_i/row_k instead of canonical (min,max). For asymmetric angles like P-O-H (rows 2,1,0), the table stores (0,1,2) but the query checked (2,1,0) → no match → stretch-bend skipped entirely. Fixed to canonicalize. Root cause for trimethylphosphine + methylphosphonic_acid.
  - **Validation (109 mols)**: energy Pearson r→**1.0000**, RMSD 0.441→**0.251 kcal/mol**, max|Δ| 2.17→1.12, outliers >1.0: **7→1**. Atom types 0.00%; charges mean|Δ|=0.0003. 165 tests pass, 0 warnings. WASM rebuilt.
  - Fixed: dimethyl_disulfide, CCl4, pyridazine, furan, pyrazole, trimethylphosphine, methylphosphonic_acid (all Δ<0.3 except mpa +0.2), iodobenzene (+0.1).
  - Remaining: **purine (+1.12)** — stretch-bend sb_type classification for mixed Double/Single bond angles (compute_stretch_bend_type gives sb_type 2 where RDKit uses 1), causing kba to fall to DFSB default (0.3) instead of specific STBN entries (0.61/0.227). Subtle algorithmic issue needing RDKit source comparison.
  - `src/molecule/graph.rs` (`find_out_of_planes`): RDKit creates **3 cyclic OOP permutations** per 3-neighbor sp2 atom (each neighbor takes a turn as the out-of-plane atom1, other two define the reference plane). WebMM was creating only 1. Tripled the OOP terms, fixing guanidinium (OOP 1.81→5.57; Δ −3.9→−0.1).
  - `src/mmff/bond.rs`: C_3-O_R entry had stale params (1.43/5.0) instead of the RDKit-verified C_3-O_3 values (1.418/5.047). With 4 C-O bonds at 1.35 Å this caused +3.3 kcal/mol excess. Fixed dioxane (Δ +2.4→+0.01).
  - `src/mmff/bond.rs`: C_AR-N_AR Aromatic entry lacked the reverse `(N_AR, C_AR, Aromatic)` arm, so N→C-ordered aromatic ring bonds fell to wrong params (r0=1.368/4.0 instead of 1.333/5.737). This corrupted BOTH bond energy AND stretch-bend (which uses bond r0 for dr), flipping the stretch-bend sign. Corrected to RDKit-verbose values r0=1.333, kb=5.737. Fixed pyrimidine (Δ −2.4→0.0; all per-terms now match RDKit: bond/angle/sb/vdw/elec).
  - **Validation (109 mols)**: energy Pearson r 0.9997→**0.9999**, RMSD 0.85→**0.441 kcal/mol**, max|Δ| 3.87→2.17, outliers >1.0: **7**. Atom types still 0.00% mismatch; charges mean|Δ|=0.0003. 165 tests pass, 0 warnings. WASM rebuilt.
  - Additional fixes this session: C_3-C_1 single (1.459/4.707) + C_1-N_1 triple (1.160/16.582) from estimator/wrong values (fixed acetonitrile Δ+1.7→0.0, propyne Δ+1.8→0.0); 5-ring specific C5B-C5B (1.418/4.313), C5A-S_AR (1.717/3.589), C5A-H (1.080/5.531), C5B-H (1.080/5.506) replacing C_AR fallback (fixed thiophene Δ−1.5→0.0, thiazole Δ+1.9→0.0).
  - Remaining worst: purine (+2.17, fused-ring sb cascade), trimethylphosphine (−1.56, P-C), dimethyl_disulfide (−1.31, S-S), carbon_tetrachloride (+1.20, C-Cl), pyridazine (−1.13), furan (−1.06), methylphosphonic_acid (+1.05) — heteroatom-specific param gaps, each <2.2 kcal/mol.
- **Final 8 typing gaps closed → 0.00% mismatch** (per PLAN.md "Close final 8 typing gaps"):
  - `src/mmff/mod.rs`: gated `has_c_c_neighbor`/`has_c_n_neighbor` on `!n_owns_cn` (this N has no own C=X double bond) so the enamine/amidine → N_PL3 rules no longer catch the =N imine itself (→ N_2 instead). Fixed guanine/purine 6-ring N (4).
  - N_4 (quaternary ammonium) charge rule now excludes acyl/amidine/enamine Ns → charged guanidinium N → N_PL3. Fixed guanidinium (1).
  - H_NAM gated on non-aromatic N (aromatic-ring N-H → H_N3); acyl detection broadened to C=C (enamine). Fixed indole (1) + guanine NH2 H (2).
  - 5-ring dicoordinate N → N5A when adjacent to NPYL/OFUR/S_AR. Fixed pyrazole (1).
  - **Validation: atom types 8 → 0 (0.00%)** on 109-mol (1324 atoms) AND 43-file sets. Charges mean |Δ| 0.0092 → 0.0003. Energy r 0.967 → 0.9754, RMSD 8.65 → 7.53. Caffeine min −123.02 (RDKit −123.49). 165 tests pass, 0 warnings.
  - Atom typing now matches RDKit 2025.09.3 exactly on the full validation set. Remaining energy gaps are parameter-table issues (5-ring/strained-ring angles, P/S params), documented in docs/validation-energy-analysis.md.
- **Remaining typing gaps fixed + energy outliers documented** (per PLAN.md "Fix remaining typing gaps"):
  - `docs/validation-energy-analysis.md`: categorizes the 15 worst energy outliers by family with root causes.
  - New MMFF types: CR3R (22, 3-ring sp3 C), HNRP (36, ammonium H), S2CM (72, CS2 S), HOS (33, sulfonic-acid H).
  - 5-ring C5A/C5B restructured to two-pass (classify N/O/S first, then C by neighbor type: adjacent to NPYL/OFUR/S_AR or flanked by 2 heteroatoms → C5A, else C5B).
  - OH2 (water O) bug: now requires both neighbors be H (was firing for O-S/P-H).
  - H_NAM/H_N3 rule corrected (was backwards): H_NAM for amide/amidine/sulfonamide/aniline N-H; H_N3 for aromatic-ring N-H (pyrrole) + simple amines. Added amidine (C=N) + sulfonamide (SO2) detection; H on O-P → H_COOH.
  - Amidine/guanidinium N → N_PL3 (has_c_n_neighbor); HNRP only for sp3 ammonium.
  - **Validation: atom types 51 → 8 (0.60%); energy r 0.954 → 0.967, RMSD 11.04 → 8.65; caffeine minimized −115.5 → −123.02 (RDKit −123.49).** 165 tests pass, 0 warnings.
  - Remaining 8 (exotic): guanidinium amidine N (1), guanine/purine 6-ring N→N_2 (4), guanine NH2 H (2), pyrazole N5A (1) — RDKit-specific nuances, minimal energy impact.
- **#1 5-ring type_ids collapse fixed (energy-neutral)** (per PLAN.md "Fix #1: 5-ring type_ids collapse"):
  - `src/mmff/mod.rs`: `type_ids` now uses `mmff_type_id(at)` for all types (no base_type collapse). The collapse was a workaround for the since-fixed NPYL charge bug; with that fixed, un-collapsing is energy-neutral (caffeine breakdown byte-identical; 109-mol r=0.954/RMSD 11.04 unchanged).
  - Validation atom-type mismatch **91 → 51 (3.85%)** — the 56 dominant 5-ring cases (C5A/C5B/NPYL/N5B/OFUR) now report their real MMFF IDs. 165 tests pass, 0 warnings.
  - Remaining (cosmetic, no energy impact): C5A-vs-C5B alpha/beta swap (13 — RDKit keys on pyrrole-vs-pyridine-N position; both EQ-fall to C_2), plus H-subtype residuals (16), HNR+ charged-amine H (5), cyclopropane CR3R (3), water OH2-vs-O_3 (3), N_2-vs-N_PL3 (3), CS2 S2CM (2), H_COOH-vs-H_ONC (2).
- **Large-scale RDKit validation harness + H_NAM subtype fix** (per PLAN.md "Large-Scale RDKit Validation"):
  - Built `scripts/gen_validation_set.py` (109 curated explicit-H 3D molecules via RDKit ETKDG), `examples/dump_types_energy.rs` (WebMM types+energy), `scripts/validate.py` (compare). PubChem-drug fetch wired via curl (`--drugs`) since the sandbox blocks `urllib`/`subprocess` import.
  - **Results (109 mols, 1324 atoms)**: atom types **91 mismatches (6.87%)**; charges mean |Δ|=0.0107 (23 outliers: ammonium/purines/CS2); energy **Pearson r=0.952, RMSD 11.3 kcal/mol** with systematic high-energy bias for 5-ring heteroaromatics + strained rings.
  - `src/mmff/mod.rs` (`determine_h_subtype`): H on amide N (bonded to C=O) now → H_NAM(28), was H_N3(23). 165 tests pass; mismatches 105→91.
  - **Prioritized roadmap** (see PLAN.md): (1) 5-ring type collapse (56 cases, biggest — would also fix the largest energy outliers), (2) H-subtype residuals (16), (3) ammonium charge bug (Δ0.94), (4) cyclopropane CR3R(22) ×3, (5) CS2 S2CM(72) ×2, (6) N_2 vs N_PL3 ×3.
  - The earlier 43-file set's 0% under-sampled 5-ring heteroaromatics + amides; the 109-set is the meaningful benchmark.
- **Final 2 atom-type mismatches fixed → 0.0% mismatch rate** (per PLAN.md "Fix Final 2 Atom-Type Mismatches"):
  - **Enamine/vinylamine N → N_PL3** (`src/mmff/mod.rs`): RDKit types an N directly bonded to a C=C carbon as planar N_PL3(40). Verified via RDKit: vinylamine `C=C-N`→40, allylamine `C=CC-N`→8 (not adjacent), aniline→40, 2-azetine `C1=CNC1`→40. Added `has_c_c_neighbor` (mirror of `has_c_o_neighbor`) and rule `(7, _, false, _) if has_c_c_neighbor => N_PL3`. Fixed Nitrogen_004 4-ring N (was N_3).
  - **O=P → O_CO2** (`src/mmff/mod.rs`): RDKit types P=O oxygen as O_CO2(32) (reuses carbonyl-O type). Added rule for O double-bonded to P (Z=15). Same params as O_2 via EQ fallback, but type now matches. Fixed Sulfur_003 O=P.
  - **Result: type diff 2 → 0 mismatches (0.0%)** across 577 atoms in 43 valid molecules (`test.sdf` skipped as a no-explicit-H input artifact). 165 tests pass, 0 warnings; demo molecules converge unchanged (caffeine −115.47). WASM rebuilt + demo re-staged.
- **N_AM typing + 5-ring N (NPYL) classification + charge fix** (per PLAN.md "N_AM Typing + 5-Ring N (NPYL) Classification + Charge Fix"):
  - **Task 1 (test.sdf)**: confirmed the 10 `test.sdf` mismatches are an input artifact — it's a 2D mol with no explicit H. Adding explicit Hs via RDKit gives WebMM/RDKit **identical** types. No implicit-H inference in WebMM (by design, targets explicit-H 3D inputs).
  - **Task 2 (N_AM)**: the assumed "N_AM under-calibration" was a misdiagnosis — N_AM typing+charges were already correct (WebMM's N_AM partial charges match RDKit *exactly* for caffeine atom8/11 and for Drug-like_Fragments_006). The real caffeine energy gap was a **5-ring N mis-typing**: caffeine's N-Me pyrrole-type N was typed N5B with charge −1.03 (RDKit NPYL, +0.05).
  - `src/mmff/mod.rs`: re-enabled amide-N typing (non-aromatic N bonded to carbonyl C → N_AM(10), was N_PL3(40)).
  - `src/mmff/mod.rs` (5-ring N): replaced "H-count" rule (N with H → NPYL, else N5B) with **coordination** rule (tricoordinate → NPYL 39; dicoordinate → N5B 66), matching RDKit. Bond order can't be used because aromatic perception already replaced Kekulé double bonds.
  - **Caffeine energy vs TRUE (sanitized) RDKit**: at-coords −88.4 (was −18.8; RDKit −107.3); **minimized −115.5** (RDKit −123.5) — gap closed from ~66 → ~8 kcal/mol. (Earlier "Δ=0.11 vs −57.07" had compared against non-sanitized RDKit, the wrong reference.)
  - Type diff **15 → 2 real mismatches (0.3%)** (test.sdf skipped as a no-explicit-H input artifact — verified identical to RDKit once Hs are added): the 3 amide-N cases resolved. Remaining 2: N in a 4-membered ring (N_PL3 vs N_3 planarity) and O=P (O_CO2 vs O_2, same params via EQ fallback). 165 tests pass, 0 warnings; all demo molecules converge (caffeine 428 iters). WASM rebuilt + demo re-staged.
- **MMFF atom/bond typing aligned with RDKit 2025.09.3** (per PLAN.md "Align WebMM MMFF Atom/Bond Types with RDKit"):
  - Built a WebMM-vs-RDKit type-diff harness: `scripts/gen_atom_type_ref.py` (RDKit types, sanitized), `examples/dump_atom_types.rs` (WebMM types), `scripts/diff_atom_types.py`. Mismatches **108 → 15 (2.5%)**.
  - `src/mmff/mod.rs` (`MMFFForceField::new`): `type_ids` previously computed as `mmff_type_id(base_type(at))`, flattening every subtype (so `atom_types` held `H_COOH`/`H_ONC`/etc. but `type_ids` reported generic `H(5)`). Now H subtypes (H_OH/H_ONC/H_COOH/H_OAR/H_N3/H_NAM/H_NIM/HS) + CR4R/CE4R keep their real MMFF IDs; 5-ring types (C5A/C5B/N5A/N5B/NPYL/OFUR) still collapse to aromatic base (C_AR/N_AR/O_3) because WebMM's parameter tables resolve aromatic-ring params via the base type. Fixed **14 H-subtype** mismatches.
  - `src/molecule/graph.rs` (`determine_hybridization`): 2-coordinate S was forced Sp2 (→ thione S_2); single-bonded 2-coord S is bent sp3 → now S_3. Fixed **3 thiol-S** cases (+ corrects thioethers).
  - New atom types: `HS` (thiol H, MMFF 71), `CR4R` (sp3 C in 4-ring, MMFF 20), `CE4R` (sp2 C in 4-ring, MMFF 30) — added to enum, `mmff_type_id`, `base_type`, `atom_types.rs` properties, and detection (`determine_h_subtype` S branch; 4-ring detection via `find_rings` in `assign_atom_types`). Fixed **3 thiol-H** + **7 cyclobutane/azetidine-C** cases.
  - `src/prop_tests.rs` (`gradient_finite_difference`): the test failed at symmetric points where the analytical gradient is exactly 0 (`rel_err = |0−noise|/noise = 1`); switched to a combined relative+absolute tolerance.
  - **Deferred** (CODE_STATUS known limitations): amide N→N_AM(10) — correct per RDKit but N_AM params under-calibrated (~12 kcal/mol worse caffeine energy); kept N_PL3(40) pending calibration (3 cases). `test.sdf` 10 mismatches are a 2D no-explicit-H input artifact (no implicit-H inference). N-in-4-ring + O=P: 1 case each, negligible param impact.
  - Verified: 165 tests pass, 0 warnings; caffeine E=−56.961 (RDKit −57.074); cyclobutane (CR4R) + thiols (S_3/HS) converge cleanly.
- **Caffeine (amide) non-convergence fixed — estimation force-constant units bug** (per PLAN.md "Fix Caffeine (Amide) Non-Convergence — Estimation Force-Constant Units Bug"):
  - Symptom: demo Caffeine hit the 1000-iter cap non-converged; WebMM reported E=+538.58 at coords where RDKit reports −29.09, with the bond-stretch term alone = 574 kcal/mol.
  - Root cause: `src/mmff/estimation.rs::estimate_bond_params` computed `k_bond = 71.9662 * (2·bc1·bc2)/(bc1+bc2)`. The `71.9662` (= C_bn/2) is the kcal conversion already re-applied in `bond_energy` (c1=143.9325), so the unit conversion ran **twice**, inflating `k_bond` ~34–65×. This fallback was hit for common-but-unlisted pairs like `C_2-N_PL3` (amide → kb=323.85 vs real ~4.5) and `C_AR-N_AR` aromatic (kb=287.87 vs real ~6), making single amide bonds cost ~210 kcal/mol and forcing ~2000 optimization iterations.
  - Fix: removed the spurious `* 71.9662`; estimator now returns the harmonic-mean-of-`bond_class` in mdyn/Å (1–17 range), e.g. `C_2-N_PL3`→4.5, `C_AR-N_AR`→4.0, `C_3-C_3`→2.0. Corrected the `test_c3_c3_single_bond` assertion that codified the bug (`>100` → range `1.0..=5.0`).
  - Verified: caffeine at given coords E +538.58→**−18.82** (RDKit −29.09); caffeine optimization **converges in 276 iters** (was ~2026), final E=**−56.961** vs RDKit-minimized −57.074 (Δ=0.11 kcal/mol). All demo molecules converge <1000 iters (aspirin 261, paracetamol 361, sulfanilamide 254, glucose 640, nicotine 216). 165 tests pass, 0 warnings, WASM rebuilt + re-staged.
  - Not fixed (separate issue): caffeine purine-ring typing (WebMM 5-ring heteroaromatic C5A/C5B/N5B vs RDKit amide C_2/N_PL3) — an aromaticity-perception gap causing the remaining ~10 kcal/mol initial-energy gap. Does not affect convergence or final-energy accuracy materially.
- **GitHub Pages deployment fixed** (per PLAN.md "Deploy WebMM Demo to GitHub Pages"):
  - Root cause: `pkg/` is gitignored (`pkg/.gitignore` = `*`), so the manually-edited `pkg/index.html` was never committed — a fresh CI checkout had no landing page.
  - `site/index.html`: rewrote as the single committed, deployment-ready demo page — imports `./webmm.js` (same-dir, correct when served from Pages root), primary "Optimize (MMFF94s)" button calls `optimize_from_sdf` (full pipeline: ETKDG for 2D, MMFF refine for 3D), "Generate 3D (ETKDG v3)" calls `generate_initial_coordinates_wasm`, 14 RDKit-generated 3D SDF fixtures (caffeine/aspirin/paracetamol/nicotine/sulfanilamide/pyridine-N-oxide/naphthalene/glucose/DMSO/cyclohexane + water/ethanol/aniline/benzene), energy+iterations+atoms+time readout, GitHub + RDKit links.
  - `.github/workflows/pages.yml`: added a "Stage landing page" step (`cp site/index.html pkg/index.html`) after `wasm-pack build`; split into `build` + `deploy` jobs (official `upload-pages-artifact`/`deploy-pages` pattern); trigger now fires on push to `main` (path-filtered: `src/**`, `site/**`, `Cargo.*`, workflow) + existing `v*` tags + `workflow_dispatch`.
  - `README.md`: updated local-dev instructions to the build → `cp site/index.html pkg/index.html` → serve `pkg/` flow that mirrors CI.
  - Verified: `wasm-pack build --target web --release` OK, workflow YAML validates, all referenced WASM exports present in `pkg/webmm.d.ts`, local HTTP smoke test returns 200 for `/`, `/index.html`, `/webmm.js`, `/webmm_bg.wasm`. No Rust source changes this task.
- **Dead code & build-warning cleanup** (per PLAN.md "Code Cleanup: Dead Code, Build Warning, Outdated Docs"):
  - `Cargo.toml`: Removed `[[bin]] name="benchmark"` target — it duplicated the auto-discovered `benches/benchmark.rs` bench target, producing a Cargo "found to be present in multiple build targets" warning
  - `src/mmff/torsion.rs`: Removed dead `calculate_dihedral` function (never called; torsion gradient computes cos(φ)/sin(φ) inline). Updated adjacent comment that referenced the removed function.
  - `src/lib.rs` (`diag_aniline` test module): Removed dead `make_aniline` and `rdkit_opt_coords` helpers and the now-unused `use crate::molecule::{Atom, Bond, BondStereo, BondType, Molecule}` import
  - `src/lib.rs` (`type_audit5` test module): Removed dead `all_type_ids` helper
  - `PROJECT_STATUS.md`: Replaced severely outdated content (claimed 5/7 tests passing, ETKDG/MMFF "not implemented") with a brief summary table + pointer to `CODE_STATUS.md`
  - **All 165 tests pass**, 0 warnings (was 3 dead-code warnings + 1 Cargo duplicate-target warning), 0 errors
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
- **Amide N typing**: RESOLVED 2026-07-24 — non-aromatic carbonyl-bonded N now typed N_AM(10) (was N_PL3); charges match RDKit exactly.
- **5-ring N classification**: RESOLVED 2026-07-24 — pyrrole-type N (tricoordinate) now NPYL(39) regardless of H vs substituent; was mis-typing N-alkyl pyrrole Ns as N5B with a −1.03 charge bug.
- **No implicit-H inference**: WebMM requires explicit-hydrogen 3D SDF inputs. 2D / no-H inputs (e.g. `test_molecules/test.sdf`) read 2-coordinate carbons as sp2 where RDKit infers sp3. Verified: adding explicit H to test.sdf makes WebMM types identical to RDKit. The type-diff harness skips no-H molecules as documented input artifacts.
- **Caffeine energy reference**: the true (sanitized) RDKit caffeine minimized energy is **−123.5** kcal/mol; WebMM achieves −115.5 (Δ≈8). The earlier-cited −57.07 was from non-sanitized RDKit (wrong aromaticity) and was an invalid reference.
- **Atom typing vs RDKit: 0.0% mismatch** on all 577 atoms of the 43 valid (explicit-H) test molecules. The only excluded file is `test.sdf` (no explicit H — input artifact).

## Fixed Issues (Post-RDKit Review)
- 2026-04-24 — Fixed stretch-bend conversion factor: added missing `0.5 * C_SB` where `C_SB = 143.9324 * π/180 ≈ 2.512` (src/mmff/stretch_bend.rs). Was producing energies ~5× too small.
- 2026-04-24 — Fixed OOP angle calculation: rewrote `calculate_oop_angle` to use proper Wilson angle (RDKit-style) with normalized vectors from central atom (src/mmff/oop.rs).
- 2026-04-24 — Fixed MMFF94s variant ignored: `get_torsion_params` now accepts and uses `mmff_variant` parameter, with MMFF94s-specific overrides for amide and peptide torsions (src/mmff/torsion.rs).
- 2026-04-27 — Completed M6-M12 gaps: SP3D/SP3D2 atom types, atropisomer chirality, fragment embedding, coordMap constraints, timeout support, et_version usage, configurable triangle smoothing epsilon
- 2026-04-27 — Fixed acetic acid carboxyl planarity: `build_planarity_constraints` now adds improper torsions and exocyclic torsions for non-aromatic sp2 atoms (carboxyl, carbonyl, amide groups). Carboxyl group now planar (dihedrals ~0°/180°). Un-ignored `test_etkdg_acetic_acid_geometry`.
- All 135 tests pass after fixes, WASM build verified.
