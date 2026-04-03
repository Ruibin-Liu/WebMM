# Plan: Fix MMFF94 Angle Energy and Parameter Accuracy

## Problem
Angle energy contributes nearly zero to total energy due to wrong constant.
Bond and angle parameters differ significantly from RDKit's MMFF94 reference values,
causing incorrect optimized geometries and energy breakdowns.

## Root Cause Analysis

### Angle Energy Constant (Biggest Bug)
The angle energy formula uses constant `0.000043945` which is off by factor of ~1.6M.
RDKit's MMFF94 implementation uses `71.9662 = C_bn/2 = 143.9324/2`.
Verified by numerical second derivative of RDKit's angle-only energy for water (3 atoms):
- d2E/dtheta2 = 94.707, Z_IJK = 0.658
- C_eff = 94.707 / (2 * 0.658) = 71.966 (matches exactly)

### Bond Parameters
Several bond kb/r0 values differ from RDKit's MMFF94 (extracted via `GetMMFFBondStretchParams`):
- O_3-H: ours 5.5/0.960 vs RDKit 7.880/0.969 (kb 30% low)
- C_3-H: ours 4.5/1.113 vs RDKit 4.766/1.093 (r0 too high)
- C_3-C_3: ours 4.7/1.526 vs RDKit 4.258/1.508
- C_2=O_2: ours 10.5/1.22 vs RDKit 12.95/1.222 (kb 19% low)
- C_AR-H: ours 4.3/1.080 vs RDKit 5.306/1.084 (kb 19% low)
- C_AR-C_AR: ours 6.0/1.39 vs RDKit 5.573/1.374
- C_3-O_3: ours 5.5/1.43 vs RDKit 5.047/1.418

### Angle Parameters
Z_IJK and theta0 values differ from RDKit:
- H-O_3-H: ours Z=0.80/t0=109.47 vs RDKit Z=0.658/t0=103.978
- H-C_3-H: ours Z=0.77/t0=109.47 vs RDKit Z=0.516/t0=108.836
- H-C_2-H: ours Z=0.45/t0=120.0 vs RDKit Z=0.594/t0=116.699
- C_AR-C_AR-C_AR: ours Z=1.05/t0=120.0 vs RDKit Z=0.669/t0=119.977

## Fix

### Step 1: Fix angle energy and gradient constants
- `src/mmff/angle.rs`: Change `0.000043945` to `71.9662` in both `angle_energy()` and `angle_gradient()`

### Step 2: Update bond parameters to match RDKit MMFF94
Update `src/mmff/bond.rs` hardcoded entries based on RDKit reference data.

### Step 3: Update angle parameters to match RDKit MMFF94
Update `src/mmff/angle.rs` hardcoded entries based on RDKit reference data.

### Step 4: Run tests and fix failures
Fix any tests that depend on the old (wrong) parameter values.

### Step 5: Rebuild WASM and verify in browser

## Changes
- `src/mmff/angle.rs`: Fix constants and update angle parameter table
- `src/mmff/bond.rs`: Update bond parameter table
- `src/lib.rs`: Fix tests that depend on old parameters
- `CODE_STATUS.md`: Update with changes

## Reference Data
RDKit MMFF94 parameters extracted via Python API for 7 molecules:
Water, Methane, Formaldehyde, Ethane, Ethanol, Acetic Acid, Benzene.
