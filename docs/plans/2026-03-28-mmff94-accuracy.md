# MMFF94 Accuracy: H Subtypes, Charges, and Stretch-Bend Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Match RDKit MMFF94 energies within ~5 kcal/mol for small molecules by fixing H atom subtypes, BCI charges, and adding the stretch-bend coupling term.

**Architecture:** Add context-dependent H subtypes to the atom type enum (H_OH, H_ONC, H_COOH, H_N3, H_NAM, H_OAR). Assign based on 2-hop neighbor context. Calibrate fbci/BCI values to produce RDKit-matching charges. Add stretch-bend energy term with numerical gradients.

**Tech Stack:** Rust, RDKit Python (reference data generation only)

---

### Task 1: Add H Subtypes to Atom Type Enum

**Files:**
- Modify: `src/mmff/mod.rs:34-88` (MMFFAtomType enum)

**Step 1: Add new H variant types**

Add after `H` in the enum:
```rust
H_OH,    // H in water (H-O-H, both neighbors are H)
H_ONC,   // H bonded to O_3 where O bonded to C_3 (ethanol)
H_COOH,  // H bonded to O_3 where O bonded to C_2 (carboxylic acid -OH)
H_OAR,   // H bonded to O_3 where O bonded to C_AR (phenol)
H_N3,    // H bonded to N_3 (ammonia)
H_NAM,   // H bonded to N_AM or N_AR (aniline, acetamide)
```

**Step 2: Run `cargo check` to find all match arm errors**

Run: `cargo check 2>&1 | grep "non-exhaustive\|H_OH\|H_ONC"`
Expected: Compilation errors in charges.rs, bond.rs, angle.rs, atom_types.rs, mod.rs where match on MMFFAtomType needs new arms.

**Step 3: Add placeholder match arms everywhere**

In every `match` on MMFFAtomType, add the new H subtypes alongside the existing `H` arm. Bond/angle params can fall through to `H` initially.

**Step 4: Verify compilation**

Run: `cargo check`
Expected: Compiles with warnings only.

---

### Task 2: Implement Context-Based H Type Assignment

**Files:**
- Modify: `src/mmff/mod.rs:137-269` (assign_atom_types method)

**Step 1: Replace simple `H` assignment with context-dependent logic**

The hydrogen match `(1, _, _, _) => MMFFAtomType::H` must become context-aware:

```rust
// Hydrogen - context-dependent subtyping
(1, _, _, _) => {
    let neighbor = mol.adjacency[idx]
        .first()
        .expect("H must have a neighbor");
    let neighbor_atom = &mol.atoms[*neighbor];
    let neighbor_type = self.determine_h_subtype(idx, mol);
    neighbor_type
}
```

**Step 2: Add `determine_h_subtype` helper method**

Logic based on RDKit observations:
- H bonded to O: check what else O is bonded to
  - O bonded only to H → H_OH (water)
  - O bonded to C_AR → H_OAR (phenol)
  - O bonded to C_2 (has C=O bond) → H_COOH (carboxylic acid)
  - O bonded to C_3 → H_ONC (ethanol)
- H bonded to N:
  - N is N_AR or N_AM → H_NAM (aniline, acetamide)
  - N is N_3 → H_N3 (ammonia)
  - N is N_2 → H_NAM (amine-like)
- H bonded to C: H (standard)

**Step 3: Run tests**

Run: `cargo test`
Expected: All 96 tests pass (H subtypes should produce same atom type for C-bonded H).

---

### Task 3: Add Properties for H Subtypes

**Files:**
- Modify: `src/mmff/atom_types.rs` (AtomTypeProperties table)

**Step 1: Add property entries for each H subtype**

Based on RDKit reference charges and standard MMFF94 values:
- H: existing (fbci=0.0, bond_class=1)
- H_OH: fbci adjusted to produce q=+0.43 (water H)
- H_ONC: fbci adjusted to produce q=+0.40 (ethanol H)
- H_COOH: fbci adjusted to produce q=+0.50 (COOH H)
- H_OAR: fbci adjusted to produce q=+0.45 (phenol H)
- H_N3: fbci adjusted to produce q=+0.36 (ammonia H)
- H_NAM: fbci adjusted to produce q=+0.37 (amide/aniline H)

The fbci values control charge initialization. The BCI table entries then determine charge transfer. We need to calibrate fbci + BCI to match RDKit charges.

**Step 2: Add bond/angle params for H subtypes in bond.rs and angle.rs**

Most bond and angle parameters for H subtypes are the same as H (bond strength doesn't depend on the rest of the molecule). Copy H entries for each subtype.

**Step 3: Add BCI table entries for H subtypes in charges.rs**

Add BCI entries for each H_x - O_3, H_x - N_x bond pair.

**Step 4: Calibrate charges against RDKit**

Use Python to compute what fbci/BCI values produce RDKit charges for water, ethanol, acetic acid, ammonia, phenol, aniline.

Target charges from RDKit:
- Water: O=-0.86, H=+0.43
- Ethanol: O=-0.68, H(OH)=+0.40
- Acetic Acid: O(COOH)=-0.65, H(OH)=+0.50
- Ammonia: N=-1.08, H=+0.36
- Phenol: O=-0.53, H(OH)=+0.45
- Aniline: N=-0.90, H=+0.40

**Step 5: Add O_3 subtype for water oxygen**

RDKit uses type 70 for water oxygen (O bonded only to H). Add O_H2O subtype or handle in charge calculation.

---

### Task 4: Add Stretch-Bend Coupling Term

**Files:**
- Modify: `src/mmff/mod.rs` (calculate_energy_and_gradient, calculate_energy_breakdown)
- Create: `src/mmff/stretch_bend.rs`

**Step 1: Implement stretch-bend energy formula**

MMFF94 stretch-bend energy for angle I-J-K:
```
E_sb = kba_IJK * (r_IJ - r0_IJ) * (theta_IJK - theta0_IJK)
     + kba_KJI * (r_KJ - r0_KJ) * (theta_IJK - theta0_IJK)
```

Where kba_IJK and kba_KJI are stretch-bend force constants (typically 0.1-0.3 mdyn·Å⁻¹·rad⁻¹).

**Step 2: Use numerical gradients for stretch-bend**

Since analytical gradients are complex (coupling bond stretch and angle bend), use numerical finite differences:
```
grad_i[dim] = (E(x_i[dim]+eps) - E(x_i[dim]-eps)) / (2*eps)
```

This is safe because stretch-bend energy is always small (< 1 kcal/mol).

**Step 3: Add stretch-bend params from RDKit**

Extract kba values from RDKit for common angle types (already available from the reference data extraction).

**Step 4: Integrate into MMFFForceField**

Add stretch-bend loop in `calculate_energy_and_gradient` and `calculate_energy_breakdown`.

---

### Task 5: Validate Against RDKit

**Files:**
- Modify: `src/lib.rs` (RDKit comparison tests)

**Step 1: Generate RDKit reference data for 12 molecules**

Use Python to extract: atom types, charges, bond/angle/sb params, energy breakdown, gradients for Water, Methane, Ammonia, Formaldehyde, Ethane, Ethanol, Acetic Acid, Benzene, Phenol, Aniline, Dimethyl Ether, Acetamide.

**Step 2: Add tight comparison tests**

For each molecule at optimized geometry:
- Total energy within 2 kcal/mol of RDKit
- Per-term breakdown within 1 kcal/mol
- Charges within 0.05 of RDKit
- Gradient components within 5 kcal/mol/Å

**Step 3: Add distorted-geometry comparison tests**

For Water, Methane, Formaldehyde, Ethane, Acetic Acid:
- Distort coordinates by 0.2-0.5 Å
- Compare energy and forces at distorted geometry
- Tolerances: energy within 5 kcal/mol, forces within 50 kcal/mol/Å

**Step 4: Run all tests**

Run: `cargo test`
Expected: All tests pass.

---

### Task 6: Rebuild WASM and Final Verification

**Files:**
- Build: `pkg/` (wasm-pack build)

**Step 1: Build WASM**

Run: `wasm-pack build --target web --out-dir pkg`

**Step 2: Run clippy**

Run: `cargo clippy`
Expected: No new warnings.

**Step 3: Update CODE_STATUS.md and PLAN.md**

Document all changes made.
