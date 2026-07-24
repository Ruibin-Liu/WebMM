# Comprehensive Review: WebMM ETKDGv3 & MMFF94s vs. RDKit

## Executive Summary

Line-by-line comparison with RDKit source code reveals significant discrepancies
that have now been **mostly resolved**. All critical, high, and medium issues fixed.
One remaining low-priority issue (E9) is intentionally kept.

---

## 1. ETKDGv3 Issues

### ~~CRITICAL Issues~~ (ALL FIXED)

#### ~~E1: UFF Electronegativity Values Wrong~~ (FIXED)

#### ~~E2: No Hybridization-Specific UFF Radii~~ (FIXED)

#### ~~E3: Chiral Center R/S Detection Incomplete~~ (FIXED — stereo parity parsing + negative bounds)

#### ~~E4: Planarity Check Wrong~~ (FIXED — divides by improper count, checks before minimization)

#### ~~E5: Triangle Smoothing Bug~~ (FIXED — uses `upper[k][j]`)

### ~~HIGH Issues~~ (ALL FIXED)

#### ~~E6: ET Version Logic Inverted~~ (FIXED — `et_version >= 3` not `>= 2`)

#### ~~E7: Torsion Pattern Ordering Wrong~~ (FIXED — chain patterns first for all bonds)

#### ~~E8: Bridged Ring Exclusion Incomplete~~ (FIXED — excludes ALL bonds of non-macrocycle rings)

#### E9: Force Constants Differ (etkdg/mod.rs) — INTENTIONAL

Uses `LONG_RANGE_FORCE=10.0` and `K_12=100.0` vs RDKit's weight=1.0 for bounds.
These were empirically calibrated. RDKit's ForceField contributions have internal
energy functions where weight=1.0 is a multiplier, not the total force constant.
Changing to 1.0 broke ETKDG conformer quality (reverted after testing).

### ~~MEDIUM Issues~~ (ALL FIXED)

#### ~~E10: Timeout Unit Mismatch~~ (FIXED — now uses seconds, matching RDKit)

#### ~~E11: RNG Different~~ (FIXED — replaced Xoshiro256** with MT19937 Mersenne Twister, matching std::mt19937)

### ~~LOW Issues~~ (ALL FIXED)

#### ~~E12: is_larger_sp2_atom Missing Ring Check~~ (FIXED)

---

## 2. MMFF94s Issues

### ~~CRITICAL Issues~~ (ALL FIXED)

#### ~~M1: Stretch-Bend Energy~~ (FIXED — removed spurious PI/180 factor)

#### ~~M2: No Linear Angle Formula~~ (FIXED — added `E = 143.9325 * ka * (1 + cos(theta))` for theta0 >= 179)

#### ~~M3: VDW Combining Rules~~ (FIXED — Gaussian R* correction + Slater-Kirkwood epsilon)

### ~~HIGH Issues~~ (ALL FIXED)

#### ~~M4: VDW 1-4 Scaling Disagreement~~ (FIXED)

#### ~~M5: Torsion Parameter Coverage~~ (FIXED — 926 entries with 5-stage equivalence protocol)

#### ~~M6: OOP Params~~ (FIXED — 117 entries with 4-type equivalence protocol)

#### ~~M7: BCI Charge Table~~ (FIXED — 498 entries from full RDKit parameter set)

### LOW Issues

#### ~~M8: MDYNE_A_TO_KCAL_MOL~~ (FIXED — 143.9325 throughout)

---

## 3. RDKit Consistency Status

### All Issues Resolved Except One

- **E1-E8, E10-E12**: Fixed (all ETKDGv3 critical/high/medium/low issues)
- **E9**: Intentional — force constants empirically calibrated (10.0/100.0 vs RDKit's 1.0 weight)
- **M1-M8**: Fixed (all MMFF94s issues)
- **M5-M7**: Fixed (full parameter tables from RDKit)

### Exact Matches with RDKit

- All 17 ETKDGv3 distance bounds constants
- 1-2, 1-3, 1-4, 1-5 distance bounds geometric formulas
- UFF electronegativity values (all 12 elements)
- Hybridization-specific UFF radii
- VDW H-bond detection and topological distance scaling
- Bond stretch energy formula (anharmonic cubic)
- Angle bend energy formula (non-linear case + linear case)
- Stretch-bend coupling energy formula
- Torsion energy formula (3-term Fourier)
- OOP Wilson angle calculation
- Buffered 14-7 VDW energy formula with Gaussian R* + Slater-Kirkwood
- Electrostatics formula: `332.0716 * qi*qj / (r+0.05)`
- Tetrahedral/chiral volume check thresholds with R/S stereochemistry
- 4D embedding pipeline step ordering
- Weight staging for first/second minimization
- Triangle smoothing using `upper[k][j]`
- ET version logic (v3 uses softer barriers)
- Torsion pattern ordering (chain first, ring fallback)
- Bridged ring exclusion (all bonds of non-macrocycle rings)
- MT19937 Mersenne Twister RNG (matches std::mt19937)
- Timeout in seconds (matches RDKit)
- MDYNE_A_TO_KCAL_MOL = 143.9325

### Full Parameter Tables

- Full MMFF94 torsion parameters: 926 entries with 5-stage equivalence protocol
- Full MMFF94s torsion parameters: separate 926-entry table (42 entries differ)
- Full MMFF94 OOP parameters: 117 entries with 4-type lookup + 4-stage equivalence
- Full MMFF94s OOP parameters: separate 117-entry table (11 entries differ)
- Full BCI charge parameters: 498 entries covering all bond types and atom types
- Equivalence level table: 99 atom types × 4 levels from RDKit's defaultMMFFDef
- Torsion type classification: single(0), aromatic(1), double(2), ring(4), generic(5)

### Architectural Differences (acceptable)

- Gradient: Numerical for torsion/OOP/VDW vs analytical (equivalent results)
- No threading/multi-conformer support
- No RMSD pruning

---

## 4. Fix Status

All items resolved except E9 (intentional):

| Issue | Status |
|-------|--------|
| E1 Electronegativity | ✅ Fixed |
| E2 Hybridization radii | ✅ Fixed |
| E3 Chiral R/S | ✅ Fixed |
| E4 Planarity check | ✅ Fixed |
| E5 Triangle smoothing | ✅ Fixed |
| E6 ET version logic | ✅ Fixed |
| E7 Torsion pattern ordering | ✅ Fixed |
| E8 Bridged ring exclusion | ✅ Fixed |
| E9 Force constants | ⚠️ Intentional (empirically calibrated) |
| E10 Timeout units | ✅ Fixed (seconds) |
| E11 RNG | ✅ Fixed (MT19937) |
| E12 is_larger_sp2_atom | ✅ Fixed |
| M1 Stretch-bend formula | ✅ Fixed |
| M2 Linear angle formula | ✅ Fixed |
| M3 VDW combining rules | ✅ Fixed |
| M4 VDW 1-4 scaling | ✅ Fixed |
| M5 Torsion parameters | ✅ Fixed (926 entries) |
| M6 OOP parameters | ✅ Fixed (117 entries) |
| M7 BCI charges | ✅ Fixed (498 entries) |
| M8 MDYNE constant | ✅ Fixed (143.9325) |
