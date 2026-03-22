# Full MMFF94 Compliance Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Implement full MMFF94 spec compliance: ring detection, proper atom typing, bond charge increments, parameter estimation rules, comprehensive parameter tables, and ETKDG improvements.

**Architecture:** Implement MMFF94 parameter estimation rules (Halgren eqs 10-16) as the primary parameter lookup mechanism, with comprehensive hardcoded tables for common drug-like molecule combinations as overrides. Add SSSR ring detection and aromaticity to enable context-sensitive atom typing. Implement bond charge increment method for non-zero electrostatics.

**Tech Stack:** Rust, serde_json

---

## Task 1: Ring Detection (SSSR)

**Files:**
- Modify: `src/molecule/graph.rs:111-116`

Implement Smallest Set of Smallest Rings (SSSR) using the BFS-based algorithm:

1. Build the adjacency list from `mol.adjacency`
2. For each atom, perform BFS to find all shortest paths between atom pairs
3. Detect rings using the path-finding approach: for each edge, find if there's a path between the two endpoints that doesn't use that edge
4. Deduplicate to get the smallest set

**Algorithm:** Use the Horton/SSSR approach:
- For each vertex v, perform BFS to build a shortest-path tree
- For each edge (u, w) not in the tree, the path tree[v]->u + (u,w) + tree[v]->w forms a ring
- Collect all rings, sort by size, remove non-primitive (composable) rings

Add test: cyclohexane (6-ring), benzene (6-ring), furan (5-ring), naphthalene (fused 6,6).

## Task 2: Aromaticity Detection (Huckel Rule)

**Files:**
- Modify: `src/molecule/graph.rs:97-109`

Replace placeholder aromaticity with Huckel 4n+2 rule:
1. Find all rings (using Task 1's `find_rings`)
2. For each ring, count pi electrons:
   - Each aromatic bond contributes 1 pi electron
   - Each sp2 atom with lone pair in ring contributes 2 (pyrrole N, furan O)
3. Ring is aromatic if pi_electrons == 4n+2 for some integer n
4. Aromaticity propagates to fused ring systems (sharing 2+ atoms)

## Task 3: Context-Sensitive Atom Typing

**Files:**
- Modify: `src/mmff/mod.rs:130-189`

Expand `assign_atom_types` to use ring membership, aromaticity, neighbor context, and formal charge:

1. Use `find_rings` to determine ring membership for each atom
2. Use proper `is_aromatic` from Task 2
3. Check `atom.charge` for formal charge (parsed from V2000 charge field)
4. Check neighbor element context:
   - N_PL3: nitrogen bonded to C=O (amide) with sp3 geometry
   - O_CO2: oxygen double-bonded to carbon
   - C_CAT/C_AN: formal charge +1/-1 on carbon
   - O_R: ether oxygen (sp3 O with 2 C neighbors)
5. Extend coverage to ~50 common types

## Task 4: MMFF94 Atom Type Property Table

**Files:**
- Create: `src/mmff/atom_types.rs`
- Modify: `src/mmff/mod.rs`

Create a comprehensive table of MMFF94 atom type properties per Halgren Table II:

```rust
pub struct AtomTypeProperties {
    pub bond_class: i32,      // b_c
    pub angle_class: i32,     // a_c
    pub fbci: f64,            // formal bond charge increment
    pub crd: f64,             // covalent radius
    pub vdw_radius: f64,      // R*
    pub vdw_epsilon: f64,     // epsilon
    pub vdw_alpha: f64,       // alpha (N/12)
    pub oop_k: f64,           // out-of-plane force constant
    pub sbl: f64,             // bond stretch-bend constant
    pub zeta: f64,            // electrostatic scaling
}

pub fn get_atom_type_props(atom_type: MMFFAtomType) -> AtomTypeProperties {
    match atom_type {
        MMFFAtomType::H => AtomTypeProperties { bond_class: 1, angle_class: 1, fbci: 0.0, crd: 0.35, vdw_radius: 1.20, vdw_epsilon: 0.0157, vdw_alpha: 0.0833, oop_k: 0.0, sbl: 0.0, zeta: 0.5 },
        MMFFAtomType::C_3 => AtomTypeProperties { bond_class: 2, angle_class: 2, fbci: 0.0, crd: 0.77, vdw_radius: 1.70, vdw_epsilon: 0.0100, vdw_alpha: 0.1667, oop_k: 0.0, sbl: 0.0, zeta: 0.5 },
        MMFFAtomType::C_2 => AtomTypeProperties { bond_class: 5, angle_class: 6, fbci: 0.0, crd: 0.67, vdw_radius: 1.68, vdw_epsilon: 0.0080, vdw_alpha: 0.1667, oop_k: 0.008, sbl: 0.0, zeta: 0.5 },
        // ... all types
    }
}
```

Values from Halgren 1996 Table II. Cover all ~34 types in our enum.

## Task 5: Bond Charge Increment Method

**Files:**
- Create: `src/mmff/charges.rs`
- Modify: `src/mmff/mod.rs:191-194`

Implement the BCI charge method (Halgren 1996 Section II.D):

1. **FBCI table**: Each atom type has a base formal bond charge (q_FBCI) from Task 4's `fbci` field
2. **BCI table**: For each pair of bonded atom types, there's a charge increment delta_ij
3. **Accumulation**: `q_i = q_FBCI(type_i) + sum over bonded j of BCI(type_i, type_j, bond_order)`
4. **Neutralization**: If total molecular charge != expected, distribute the residual proportionally

The BCI table needs ~200-300 entries for common type pairs. Use the Halgren Table VI values:
- C_3-C_3: 0.00
- C_3-H: 0.00  
- C_3-N_3: 0.08
- C_3-O_3: 0.13
- C_2=O_2: 0.10
- N_3-H: 0.00
- etc.

For missing combinations, use the estimation rule: `BCI(i,j) = (fbci_i - fbci_j) / (num_bonds_i + num_bonds_j)`.

## Task 6: Parameter Estimation Rules

**Files:**
- Create: `src/mmff/estimation.rs`

Implement Halgren equations 10-16 for computing bond and angle parameters from class values:

**Bond estimation (Eqs 10-13):**
```
kb_ij = (2 * G * bc_i * bc_j) / (bc_i + bc_j)  where G = 143.9324 * 0.5 = 71.9662
r0_ij = r_i + r_j - 0.01 * |crd_i - crd_j|^2  (Eq 11 for single bonds)
```

For double bonds: apply correction factor ~0.94 to r0
For triple bonds: apply correction factor ~0.90 to r0

**Angle estimation (Eqs 14-16):**
```
ktheta_ijk = C * sqrt(ac_j * ac_i * ac_k) / 2  where C varies by central atom type
theta0_ijk = 360.0 / (n_neighbors_central) + correction based on pi bonds
```

This gives reasonable parameters for ANY type combination without exhaustive tables.

## Task 7: Expanded Bond Parameters

**Files:**
- Modify: `src/mmff/bond.rs`
- Modify: `data/mmff94_sample_parameters.json`

Add comprehensive bond parameters from Halgren Table III for drug-like molecules (~100 entries):

**Critical additions:**
- All H-X bonds (H-C_3, H-C_2, H-C_AR, H-N_3, H-N_2, H-N_AR, H-O_3, H-O_2, H-S_3)
- Aryl-alkyl bonds (C_AR-C_3, C_AR-C_2, C_AR-N_3, C_AR-O_3, C_AR-S_3)
- Aromatic heterocycle bonds (C_AR-N_AR, C_AR-O_R, C_AR-S_AR)
- Carbonyl bonds (C_2=O_2, C_AR=O_2, C_2=O_R)
- N-containing bonds (N_3-C_3, N_2-C_2, N_1-C_1, N_3-N_3, N_AR-C_AR)
- Halogen bonds (C_3-F, C_3-Cl, C_3-Br, C_3-I, C_AR-F, C_AR-Cl)
- S/P bonds (C_3-S_3, C_2=S_2, C_3-P_3, P_3-O_3, P_4=O_2)

Use parameter estimation (Task 6) as fallback for missing combinations.

## Task 8: Expanded Angle Parameters

**Files:**
- Modify: `src/mmff/angle.rs`
- Modify: `data/mmff94_sample_parameters.json`

Add comprehensive angle parameters (~150 entries):

**Critical additions:**
- H-X-H and H-X-Y angles (H-C_3-H, H-C_2-H, H-N_3-H, H-O_3-H, H-C_3-N_3, etc.)
- Aryl angles (C_AR-C_AR-C_AR, C_AR-C_AR-C_3, C_3-C_AR-C_3, C_AR-C_AR-N_AR, etc.)
- Carbonyl angles (C_2-C_3-O_2, C_3-C_2=O_2, C_AR-C_2=O_2, O_2-C_2-O_2)
- Amide angles (N_3-C_2=O_2, C_3-N_3-C_2, C_3-N_3-C_3, H-N_3-C_2)
- Ether/ester angles (C_3-O_3-C_3, C_3-O_2=C_2, C_3-O_3-H)
- Heteroaromatic angles (C_AR-N_AR-C_AR, C_AR-O_R-C_AR, etc.)

Use angle estimation (Task 6) as fallback.

## Task 9: Expanded Torsion Parameters

**Files:**
- Modify: `src/mmff/torsion.rs`
- Modify: `data/mmff94_sample_parameters.json`

Add comprehensive torsion parameters (~100 entries):

**Critical additions:**
- X-C_3-C_3-Y (ethane-like, V1=0, V2=0.2, V3=0.15 for sp3-sp3)
- X-C_2=C_2-Y (ethylene-like, V1=0, V2=10-15, V3=0 for conjugated)
- X-C_AR-C_AR-Y (aryl-aryl, V1=0, V2=0-2, V3=0-1)
- X-C_3-C_2=Y (partial double bond)
- C_3-N_3-C_3-C_3 (amine, V1=0, V2=0.2, V3=0.0)
- C_3-C_2=O_2 (amide, V1=0, V2=10, V3=0.5)
- X-C_3-O_3-C_3 (ether, V1=0, V2=0.5, V3=0.0)
- X-C_AR-C_AR-X (diaryl, V1=0, V2=0-5, V3=0)
- H-C_3-C_3-H (ethane, V1=0, V2=0.2, V3=0.15)

Use default V1=0, V2=0, V3=0 for unknown torsions (flat energy landscape).

## Task 10: Per-Type VDW and OOP Parameters

**Files:**
- Modify: `src/mmff/vdw.rs`
- Modify: `src/mmff/oop.rs`

Replace per-element VDW parameters with per-type values from atom type properties (Task 4):
- C_3: R*=1.70, eps=0.010, alpha=0.1667
- C_2: R*=1.68, eps=0.008, alpha=0.1667
- C_AR: R*=1.6925, eps=0.005, alpha=0.1667
- N_3: R*=1.55, eps=0.020, alpha=0.1667
- N_AR: R*=1.55, eps=0.020, alpha=0.1667
- O_3: R*=1.52, eps=0.016, alpha=0.1667
- O_2: R*=1.50, eps=0.015, alpha=0.1667
- etc.

For OOP: per-type k_oop values:
- C_AR: 0.04
- C_2: 0.008
- N_AR: 0.02
- N_2: 0.02
- O_R: 0.02
- All others: 0.0

## Task 11: ETKDG 1-3/1-4 Distance Bounds

**Files:**
- Modify: `src/etkdg/mod.rs`

Add angle-derived (1-3) and torsion-derived (1-4) distance bounds:

1. **1-3 bounds**: For each angle (i,j,k), compute:
   ```
   r_ik^2 = r_ij^2 + r_jk^2 - 2*r_ij*r_jk*cos(theta0)
   r_ik_lower = sqrt(r_ik^2) - tolerance
   r_ik_upper = sqrt(r_ik^2) + tolerance
   ```
   Where r_ij, r_jk come from bond parameters and theta0 from angle parameters.

2. **1-4 bounds**: For each torsion (i,j,k,l), use the 1-3 bounds and the torsion's preferred angle to compute:
   ```
   For staggered (phi=180): r_il_upper = r_ij + r_jk + r_kl
   For eclipsed (phi=0): r_il_lower = sqrt(r_ij^2 + r_jk^2 + r_kl^2 - 2*r_ij*r_jk - 2*r_jk*r_kl + 2*r_ij*r_kl)
   ```

3. **Ring closure**: After finding rings, tighten upper bounds for ring atoms:
   ```
   For an n-membered ring: sum of n bond lengths = upper bound
   ```

## Task 12: ETKDG Proper 4D-to-3D Projection

**Files:**
- Modify: `src/etkdg/mod.rs:253-262`

Replace trivial "drop 4th coordinate" with eigenvector projection:
1. Compute the 4D Gram matrix G where G[i][j] = |r_i - r_j|^2
2. Convert to distance geometry metric matrix B: B[i][j] = -0.5 * (G[i][j] - G[i][0] - G[0][j])
3. Find the 3 largest eigenvalues/eigenvectors of B
4. Project coordinates onto these 3 eigenvectors

## Task 13: ETKDG FF-Based Refinement

**Files:**
- Modify: `src/etkdg/mod.rs:265-343`

Replace the simplified DGFF with the actual MMFF94 force field + L-BFGS optimizer:
1. After 4D embedding and projection, run the existing `optimize()` function
2. Use the MMFF94 force field (already implemented) instead of the custom DGFF
3. This gives much better initial coordinates since it uses the real energy function

## Task 14: ETKDG Multi-Conformer

**Files:**
- Modify: `src/etkdg/mod.rs`

Generate `max_attempts` conformers and return the best:
1. Loop N times: generate 4D coords -> project -> MMFF optimize
2. Score each by final MMFF energy
3. Return the lowest-energy conformer's coordinates
4. Add diversity check: skip conformers too similar to previously accepted ones

## Task 15: Integration Tests

**Files:**
- Modify: `src/lib.rs`

Add end-to-end tests for the new features:
1. Ring detection with cyclohexane, benzene, naphthalene
2. Atom typing correctness for aromatic N, amide N, carbonyl O
3. BCI charges: verify amine O has negative charge, amine N has partial negative
4. Full optimization of a real drug-like molecule (e.g., acetaminophen/paracetamol)
5. ETKDG quality: verify bond angles are reasonable for a benzene ring
6. Multi-conformer: verify different conformers are generated

## Task 16: Update Documentation

**Files:**
- Modify: `README.md`
- Modify: `CODE_STATUS.md`

Update README to reflect full MMFF94 compliance. Remove all TODO/partial markers.

## Out of Scope
- MMFF94s charge equilibration (complex, optional)
- CSD torsion knowledge terms (requires external statistical database)
- Stretch-bend coupling term
- Performance optimization
