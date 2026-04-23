# Plan: Fill All Remaining Gaps — MMFF Atom Type Expansion + ETKDG Refinement

## Problem
Our MMFF atom type enum is missing several types that RDKit uses for 5-membered heteroaromatic rings (pyrrole, furan, thiophene, imidazole) and some simple molecules (water, alcohols). This causes incorrect force field parameters and ETKDG geometry.

## Missing Types (from RDKit/Towhee MMFF94)

| Num | Name | Where Used |
|-----|------|-----------|
| 21 | `HOR` / `H_ONC` | H on alcohol O (methanol, ethanol) — currently mis-assigned as `H_OH` (31) |
| 39 | `NPYL` | Pyrrole N, imidazole N with H |
| 59 | `OFUR` | Furan O |
| 63 | `C5A` | 5-ring C alpha to heteroatom |
| 64 | `C5B` | 5-ring C beta to heteroatom |
| 65 | `N5A` | 5-ring N alpha position (rare) |
| 66 | `N5B` | 5-ring N beta position (imidazole pyridine-like N) |
| 70 | `OH2` | Water O — currently mis-assigned as `O_3` (6) |

## Implementation Steps

### Step 1: Add enum variants
- Add `NPYL`, `OFUR`, `C5A`, `C5B`, `N5A`, `N5B` to `MMFFAtomType` enum
- Update `mmff_type_id()` mapping in `charges.rs`
- Ensure existing `N_AR=38`, `O_R=70`, `S_AR=44` are correct

### Step 2: Fix atom type assignment logic
- **Water**: O with 2 H neighbors → `OH2` (70), not `O_3` (6)
- **Alcohols**: H on O bonded to C → `HOR` (21), not `H_OH` (31)
- **5-membered heteroaromatics**: 
  - Detect 5-membered rings with heteroatoms (N, O, S)
  - C bonded to heteroatom → `C5A` (63)
  - C not bonded to heteroatom → `C5B` (64)
  - N with H in 5-ring → `NPYL` (39)
  - N without H, double-bonded in 5-ring → `N5B` (66)
  - O in 5-ring → `OFUR` (59)
  - S in 5-ring → `S_AR` (44) — already correct

### Step 3: Update property table
- Add property entries for new types in `atom_types.rs`
- Use values from MMFF literature / RDKit reference

### Step 4: Update parameter tables
- Ensure bond, angle, torsion, VDW, OOP parameter lookups include new types
- Add fallback/estimation rules for new type combinations

### Step 5: Validation
- Compare our types vs RDKit for: pyrrole, furan, thiophene, imidazole, water, methanol, ethanol, phenol
- All tests pass
- WASM build verified
