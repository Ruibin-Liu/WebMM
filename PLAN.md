# PLAN.md

## Task: ETKDGv3 Complete Gap Analysis — WebMM vs RDKit

### Summary
Line-by-line comparison of `src/etkdg/mod.rs` (2016 lines) against RDKit's `BoundsMatrixBuilder.cpp`,
`Embedder.cpp`, and `TorsionPreferences.cpp`. **6 critical, 8 significant, 12 minor gaps found.**

---

## CRITICAL GAPS

### C1. Torsion Preferences: 15 patterns vs ~567 SMARTS patterns
- **RDKit**: ~267 V2 chain + ~300 macrocycle SMARTS patterns, matched via substructure search.
  Patterns differentiated by H-count (NX3H0 vs NX3H1), aromatic substitution (cH0 vs cH1),
  neighbor types, ring size (r5, r{9-}), CH2/CH1/CH0. First match wins per central bond.
- **Ours**: 15 hand-coded categories (`match_torsion_pattern` lines 1418-1584). No SMARTS, no
  H-count differentiation, no neighbor-type specificity.
- **Impact**: Wrong torsion preferences for substituted aromatics, amides, esters, heterocycles.

### C2. Bond Length: Covalent radii vs UFF formula
- **RDKit**: `r = r_UFF1 + r_UFF2 - 0.01 * ln(bondOrder) * (sqrt(r_UFF1) + sqrt(r_UFF2))`
  with per-type UFF radii (C_3=0.757, C_2=0.732, N_3=0.700, O_3=0.658, etc.)
- **Ours**: `single = cov_radius(e1) + cov_radius(e2)`, then `* 0.86/0.78/0.93` by bond order.
- **Impact**: Bond lengths differ by 0.01-0.05 Å, propagating through all distance bounds.

### C3. Missing First 4D Minimization
- **RDKit**: After MDS → 4D minimize (bounds w=1.0, chiral w=1.0, 4th-dim w=0.1, 400 iter/pass,
  repeat until converged) → energy/atom < 0.05 check → tetrahedral check → chiral check →
  minimize 4th dim (chiral w=0.2, 4th-dim w=1.0, 200 iter/pass) → 3D ETKDG minimize (300 iter).
- **Ours**: After MDS → minimize 4th dim only (simple gradient) → 3D minimize (400 iter) →
  checks → flatten aromatic → 3D minimize again (300 iter).
- **Missing**: First 4D minimization entirely. Energy-per-atom rejection. Chiral weight staging.
- **Impact**: Initial coords far from satisfying bounds/chirality; more failed attempts.

### C4. 1-4 Bounds: Missing Many Topology Cases
- **RDKit** has 7 distinct handlers: all-in-same-ring (size 3-5: no bounds; size 6-8: cis/sp2;
  size >=9: macrocycle), two-in-same-ring, two-in-different-rings, share-ring-bond, chain.
- **Ours**: Ring paths (generic) + chain paths (generic). Missing: two-in-same-ring,
  two-in-different-rings, share-ring-bond, cumulated doubles, S-S 90° torsion,
  stereo-tagged bonds, macrocycle amide +0.1 offset, fused ring short-circuit.
- **Impact**: Wrong distance bounds for polycyclic, bridged, fused ring, and stereo molecules.

### C5. Missing Basin Threshold for Distance Bounds
- **RDKit**: Only includes pairs where `upper - lower <= 5.0` in force field.
- **Ours**: All pairs included regardless of bounds width.
- **Impact**: Loosely-constrained pairs produce noisy gradients, poor convergence.

### C6. Force Field Structure Mismatch
- **RDKit**: Separate terms with distinct forces — torsion (Fourier), improper (UFF inversion,
  scaling=10.0), 1-2 (force=100), 1-3 (force=100), long-range (force=10), optional CPCI.
- **Ours**: Distance bounds (1/range weight), 1-2/1-3 (force=100, **BUT always satisfied** due
  to computing bounds from current distance), chiral, planarity, torsion.
- **Bug on lines 1682-1696**: `lo = d - 0.01; hi = d + 0.01` where d is current distance,
  so `lo <= d <= hi` always true → 1-2 and 1-3 constraints do NOTHING.
- **Missing**: Long-range distance constraints (force=10.0 from bounds matrix).

---

## SIGNIFICANT GAPS

| # | Gap | RDKit | Ours |
|---|-----|-------|------|
| S1 | Conjugated 5-ring squish | extraSquish=0.2 for Z>10 in conj 5-ring | Missing |
| S2 | Fused small ring volume relaxation | volScale=0.25 for atoms in >1 ring <5 | Missing |
| S3 | Non-SP3 hybridization angles | SP3D=105°, SP3D2=135°, degree>4 special | Only SP1/SP2/SP3 |
| S4 | Amide/ester 1-5 neighbor detection | `_checkAmideEster15` | Missing |
| S5 | Bridged ring torsion exclusion | Exclude fused/bridged ring bonds | Only >3 ring bonds |
| S6 | Planarity check method | Improper energy > n*0.7 | PCA deviation > 0.1Å |
| S7 | randNegEig handling | Random coords for negative eigenvalues | Matches (random) |
| S8 | Embedding attempts | 10*N_atoms | Always 10 |

---

## MINOR GAPS

| # | Gap | Notes |
|---|-----|-------|
| M1 | VDW clash check | We have it, RDKit doesn't (extra, not wrong) |
| M2 | Bond length validation | We have it, RDKit doesn't |
| M3 | Explicit aromatic flattening | We have it, RDKit uses improper energy instead |
| M4 | sameSide tolerance | We use 0.30 everywhere; RDKit uses 0.10 for chiral |
| M5 | Random seed default | We use 42; RDKit uses -1 (random) |
| M6 | SP3D/SP3D2 hybridization | Not supported in our code |
| M7 | Atropisomer chirality | Not supported |
| M8 | Fragment embedding | Not supported (embeds all together) |
| M9 | CoordMap (constrained atoms) | Not supported |
| M10 | Timeout | Not supported |
| M11 | et_version field | Defined but unused (no v1 fallback) |
| M12 | Triangle smoothing convergence | May have minor epsilon differences |

---

## CONSTANT COMPARISON

| Constant | RDKit | Ours | Match? |
|----------|-------|------|--------|
| DIST12_DELTA | 0.01 | 0.01 | YES |
| DIST13_TOL | 0.04 | 0.04 | YES |
| GEN_DIST_TOL | 0.06 | 0.06 | YES |
| DIST15_TOL | 0.08 | 0.08 | YES |
| VDW_SCALE_15 | 0.7 | 0.7 | YES |
| H_BOND_LENGTH | 1.8 | 1.8 | YES |
| MAX_UPPER | 1000.0 | 1000.0 | YES |
| MIN_TETRAHEDRAL_CHIRAL_VOL | 0.50 | 0.50 | YES |
| TETRAHEDRAL_CENTERINVOLUME_TOL | 0.30 | 0.30 | YES |
| MAX_MINIMIZED_E_PER_ATOM | 0.05 | 0.05 | YES |
| minMacrocycleRingSize | 9 | 9 | YES |
| First min iterations | 400/pass, repeat | MISSING | NO |
| 4th dim min iterations | 200/pass, repeat | 200 total | DIFFERENT |
| ETKDG 3D iterations | 300 single pass | 400+300 | DIFFERENT |
| Basin threshold | 5.0 | MISSING | NO |
| KNOWN_DIST_FORCE_CONSTANT | 100.0 | 100.0 | YES |
| KNOWN_DIST_TOL | 0.01 | 0.01 | YES |
| Improper force scaling | 10.0 (uniform) | Element-dependent (400/300/800/200) | DIFFERENT |
| Planarity tolerance | 0.7 energy/atom | 0.1 Å deviation | DIFFERENT |
| Long-range force | 10.0 | MISSING | NO |
| Linear double bond tol | 1e-3 | 1e-3 | YES |

---

## PRIORITIZED FIX PLAN

### Phase 1: Fix Critical Pipeline Issues ✅
1. ✅ Add first 4D minimization (bounds w=1.0, chiral w=1.0, 4th-dim w=0.1, 400 iter/pass)
2. ✅ Add energy-per-atom check (0.05) after first minimization
3. ✅ Fix 1-2/1-3 constraint bug (use TARGET distance, not current distance)
4. ✅ Add basin threshold (5.0) to distance bounds force field
5. ✅ Add long-range distance constraints (force=10.0 from bounds matrix)
6. ✅ Scale attempts to molecule size (10 * n_atoms)

### Phase 2: Fix Distance Bounds ✅
7. ~~Replace covalent radius bond lengths with UFF formula or table~~ (NOT DONE)
7. ✅ Replaced with UFF bond length formula: Pauling correction + O'Keefe/Breese electronegativity correction
8. ✅ Add conjugated 5-ring squish (0.2)
9. ✅ Add cumulated double bond detection for 1-4 (force cis)
10. ✅ Add S-S 90-degree torsion for 1-4
11. ~~Add stereo-tagged bond handling for 1-4~~ (NOT DONE — requires bond stereo data in parser)
12. ✅ Add two-in-same-ring (SP2→trans) and two-in-different-rings 1-4 cases
13. ✅ Add macrocycle amide +0.1 offset
14. ✅ Improve amide/ester 1-4 detection with H-count (H→trans, non-H→cis)

### Phase 3: Fix Torsion Preferences ✅
15. ✅ Expanded from 15 to ~50+ hand-coded patterns with atom-environment matching
16. SMARTS matching engine — NOT DONE (using hand-coded patterns instead)
17. ✅ Macrocycle torsion patterns (separate handling for macrocycle, 5-ring, 3-4 ring, 6-8 ring)
18. ✅ Bridged ring exclusion
19. ✅ Per-bond-done tracking

### Phase 4: Fix Secondary Issues ✅
20. ✅ Added fused small ring volume relaxation (volScale=0.25)
21. ✅ Fixed sameSide tolerance (0.10 for chiral, 0.30 for tetrahedral)
22. ✅ Planarity check now uses improper energy method (energy/atom < 0.7, matching RDKit)
23. ✅ Degree > 4 handling via SP3D/SP3D2 hybridization
24. ✅ SP3D (105°) and SP3D2 (90°) support added

### Phase 5: Edge Cases
25. Add fragment embedding (NOT DONE)
26. Add amide 1-5 neighbor detection (NOT DONE)
27. Add fused ring short-circuit in 1-4 (NOT DONE)
28. Add atropisomer chirality (NOT DONE)

## Verification
After each phase: `cargo test` and `cargo clippy` must pass cleanly.
