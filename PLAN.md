# Plan: Fix BCI Charge Calculation to Match RDKit MMFF94

## Problem
BCI charge calculation produces inaccurate charges for acetic acid and other molecules,
causing electrostatic energy to be ~60 kcal/mol off from RDKit's MMFF94 reference.

## Root Cause Analysis

Our current BCI implementation is fundamentally different from RDKit's:

### Our Implementation (WRONG)
1. Simple lookup table of ~37 BCI values by atom type pair
2. Ad-hoc fallback: `(fbci_j - fbci_i).abs() * bond_order * 0.3 / (degree_i + degree_j)`
3. Charge direction determined by comparing fbci values
4. Post-hoc neutralization to match formal charge

### RDKit's Implementation (CORRECT - Halgren 1996 Eq. 15)
1. Large BCI table of 498 entries indexed by (bondType, iAtomType, jAtomType)
2. PBCI (partial bond charge increment) per atom type for fallback: `BCI_ij = pbci_i - pbci_j`
3. Directional BCI lookup with sign convention: if i > j, swap and negate sign
4. Full Eq. 15: `q_i = (1 - M_i * v_i) * q0_i + v_i * SUM(q0_j) + SUM(BCI_ij)`
   - For neutral molecules with v_i = 0: `q_i = q0_i + SUM(BCI_ij)`
   - q0_i = MMFF formal charge (usually 0 for neutral atoms)
   - M_i = coordination number, v_i = formal charge adjustment factor
5. No separate neutralization step needed - the formula inherently produces correct total charge

### Key Differences for Acetic Acid
For acetic acid (C_3-C_2(=O_2)-O_3-H_COOH):
- Our O_3-H BCI: 0.40 (same for all H-O types)
- RDKit O_3-H_COOH BCI: 0.50 (type 6,24 entry: +0.5000)
- Our C_2=O_2 BCI: 0.42
- RDKit C_2-O_2 BCI: -0.5700 (much larger magnitude, different sign convention)
- Our C_2-O_3 BCI: 0.35
- RDKit C_2-O_3 BCI: -0.1500

## Fix

### Step 1: Replace BCI table in `src/mmff/charges.rs`
- Replace the 37-entry `get_bci()` function with a complete BCI table from RDKit's MMFFCHG.PAR
- Use MMFF numeric atom types as keys (1-99)
- Include bond type (0=single, 1=double, 4=aromatic) in lookup
- Implement directional lookup with sign convention (swap + negate if i > j)

### Step 2: Add PBCI table for fallback
- Add PBCI values from RDKit's MMFFPBCI.PAR
- When no explicit BCI exists, use: `BCI_ij = pbci_i - pbci_j`
- Also store `fcadj` (formal charge adjustment) and `crd` (coordination number) per type

### Step 3: Implement Eq. 15 charge calculation
- For neutral molecules (v=0): `q_i = SUM(BCI_ij for all bonds j of atom i)`
- Add MMFF formal charge handling for charged species
- Remove old neutralization step

### Step 4: Add atom type number mapping
- Map our MMFFAtomType enum to RDKit numeric types (1-99)
- Use this mapping for BCI/PBCI lookup

### Step 5: Run tests and fix failures

### Step 6: Rebuild WASM

## Changes
- `src/mmff/charges.rs`: Replace with RDKit-compatible BCI implementation
- `src/mmff/atom_types.rs`: Add numeric type ID to AtomTypeProperties
- `CODE_STATUS.md`: Update with changes

## Reference Data
RDKit MMFF94 BCI/PBCI parameters from Code/ForceField/MMFF/Params.cpp
