# Plan: val_set_new4 — Final 6 Exotic Atom Types

## Goal
Achieve 191/191 molecules matching RDKit 2025.09.3 to <0.01 kcal/mol by fixing
the remaining 6 molecules in `scripts/val_set_new4/` that use exotic atom types
(48, 54, 56, 67, 80, 81, 82).

## Status: COMPLETE ✅

All 6/6 molecules now match. Total: **191/191**.

## Completed Steps

### 1. Aromaticity Perception (DONE)
- `src/molecule/graph.rs`: Charged N in 5-ring contributes 1 pi electron
  (was 2). Enables imidazolium/N-oxide 5-ring aromaticity detection.

### 2. Atom Type Detection (DONE)
- `src/mmff/mod.rs`:
  - 5-ring N: if any ring N charged → N_5POS (81) for all tricoordinate N
  - 5-ring C: adjacent to N_5POS → C5A_M (78); between two → C_IM (80)
  - S with S=N + S=O → S_O2 (18)
  - O on S with 2+ double bonds → O_CO2 (32)
  - Non-aromatic sp2 N+ with N=O → N_5OX (67)
  - H on charged N with C=N double bond → HNRP (36)

### 3. Charge Distribution (DONE)
- `src/mmff/charges.rs`:
  - Type 56 (N_GD): shares +1 from CGD+ C equally among NCN+/N_GD N's
  - Type 81 (N_5POS): uses SDF formal charge

### 4. Parameter Tables (DONE)
- `src/mmff/bond.rs`: ~30 new entries for types 54, 56, 67, 71, 78, 80, 81, 82
- `src/mmff/mmff_tables.rs`: ~30 new angle entries + 1 STBN entry
- Angle entries at correct angle_type=1,2 for N_5OX (67 is sbmb=true)
- Bond entries for both H (type 5) and H_NAM (type 28) variants

### 5. Validation (DONE)
- 191/191 molecules match <0.01 kcal/mol
- 191 tests pass, 0 warnings, 0 regressions
