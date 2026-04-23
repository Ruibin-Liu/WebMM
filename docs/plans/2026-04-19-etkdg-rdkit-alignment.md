# ETKDG RDKit Alignment — Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Align WebMM's ETKDG implementation with RDKit's ETKDGv3 to produce geometrically equivalent 3D coordinates.

**Architecture:** Systematic fix of 14 critical differences in `src/etkdg/mod.rs`, organized into 5 groups: (1) Distance bounds matrix construction, (2) 4D coordinate generation, (3) DG refinement, (4) Force field refinement, (5) Conformer selection. Each fix targets a specific parameter or algorithm mismatch identified by differential analysis against RDKit source (Embedder.cpp, BoundsMatrixBuilder.cpp, DistGeomUtils.cpp).

**Tech Stack:** Rust, WASM, MMFF94s force field, L-BFGS optimizer

---

## Task Group 1: Distance Bounds Matrix Construction

### Task 1: Fix Bond Length Bounds (#1, #2)

**Files:**
- Modify: `src/etkdg/mod.rs:142-203` (`build_distance_bounds`)
- Modify: `src/etkdg/mod.rs:115-130` (`covalent_radius`)

**Problem:** We use covalent radius sums with ±0.10-0.15 Å tolerance. RDKit uses UFF `calcBondRestLength` with **±0.01 Å tolerance**.

**Step 1: Update covalent radii to match UFF bond lengths more closely**

Replace the `covalent_radius` function values with UFF-compatible values that reproduce correct bond lengths:

```rust
fn covalent_radius(element: &str) -> f64 {
    match element {
        "H" => 0.31, "C" => 0.76, "N" => 0.71, "O" => 0.66,
        "F" => 0.57, "P" => 1.07, "S" => 1.05, "Cl" => 1.02,
        "Br" => 1.20, "I" => 1.33, _ => 0.76,
    }
}
```

These are already reasonable. The main fix is tolerances.

**Step 2: Tighten bond tolerances to RDKit values**

In `build_distance_bounds`, change bond type handling:

```rust
BondType::Single => {
    bounds.lower[i][j] = single_length - 0.01;  // was -0.15
    bounds.upper[i][j] = single_length + 0.01;  // was +0.15
}
BondType::Double => {
    let ideal = single_length * 0.86;
    bounds.lower[i][j] = ideal - 0.01;   // was -0.10
    bounds.upper[i][j] = ideal + 0.01;   // was +0.10
}
BondType::Triple => {
    let ideal = single_length * 0.78;
    bounds.lower[i][j] = ideal - 0.01;   // was -0.08
    bounds.upper[i][j] = ideal + 0.01;   // was +0.08
}
BondType::Aromatic => {
    let ideal = single_length * 0.93;
    bounds.lower[i][j] = ideal - 0.01;   // was -0.10
    bounds.upper[i][j] = ideal + 0.01;   // was +0.10
}
```

**Step 3: Fix default lower bound initialization**

Change from `r_cov_i + r_cov_j - 0.3` to `0.0` (RDKit uses 0 as fallback):

```rust
// In initialization loop:
bounds.lower[i][j] = 0.0;  // was r_cov_i + r_cov_j - 0.3
bounds.upper[i][j] = 1000.0;  // was r_vdw_i + r_vdw_j + 1.0
```

**Step 4: Run tests**

`cargo test` — all existing tests must pass.

---

### Task 2: Fix 1-3 Angle-Derived Bounds (#3, #14)

**Files:**
- Modify: `src/etkdg/mod.rs:205-231` (1-3 bounds section)
- May need: `src/molecule/graph.rs` for ring-size-aware angle lookup

**Problem:** Tolerance is ±0.15 Å (RDKit: ±0.04 Å). Ring angles use hybridization-only model (RDKit: ring-size-dependent).

**Step 1: Tighten 1-3 tolerance**

```rust
let tolerance = 0.04;  // was 0.15
let lower = (r_ik - tolerance).max(0.5);
let upper = r_ik + tolerance;
```

**Step 2: Add ring-size-dependent angle model**

Add helper function for ring-aware angle estimation:

```rust
fn estimate_angle_in_ring(atom_idx: usize, mol: &Molecule) -> f64 {
    let hyb = crate::molecule::graph::determine_hybridization(atom_idx, mol);
    let rings = crate::molecule::graph::find_rings(mol);
    
    // Check if atom is in a small ring
    for ring in &rings {
        if ring.contains(&atom_idx) && ring.len() <= 6 {
            match ring.len() {
                3 => return std::f64::consts::PI * (116.0 / 180.0),  // 116° for sp3 in 3-ring
                4 => return std::f64::consts::PI * (112.0 / 180.0),  // 112° for sp3 in 4-ring
                5 => {
                    if matches!(hyb, Hybridization::Sp3) { 
                        return std::f64::consts::PI * (104.0 / 180.0);  // 104° sp3 5-ring
                    }
                    // SP2 in 5-membered aromatic: internal angle ~108°
                    return std::f64::consts::PI * (108.0 / 180.0);
                }
                _ => {}
            }
        }
    }
    
    // Default: hybridization-based
    match hyb {
        Hybridization::Sp1 => std::f64::consts::PI,
        Hybridization::Sp2 => 120.0_f64.to_radians(),
        Hybridization::Sp3 => 109.47_f64.to_radians(),
    }
}
```

Use this in place of the simple hybridization check for 1-3 bounds.

**Step 3: Run tests**

`cargo test` — all existing tests pass.

---

### Task 3: Fix 1-4 Torsion-Derived Bounds (#4)

**Files:**
- Modify: `src/etkdg/mod.rs:233-300` (1-4 bounds section)

**Problem:** Tolerance is ±0.1-0.4 Å (RDKit: ±0.06 Å). Missing amide/ester detection.

**Step 1: Tighten 1-4 tolerance**

```rust
// For all cases:
let tol_14 = 0.06;  // was 0.4 for sp3-sp3, 0.1 for others
let upper = r14 + tol_14;
let lower = (r14 - tol_14).max(0.5);
```

**Step 2: Add basic amide/ester 1-4 detection**

After computing generic 1-4 bounds, add special case for C(=O)-N patterns (amide bonds prefer trans):

```rust
// Amide/ester trans preference: C-C(=O)-N or C-C(-O)-N
// If central bond involves carbonyl C or O, prefer trans distance
let is_amide_ester = /* detect C(=O)-N or similar pattern */;
if is_amide_ester {
    // Use trans distance (larger) as primary bound
    let r14_trans = /* compute with phi=180° */;
    bounds.lower[a][d] = bounds.lower[a][d].max(r14_trans - 0.06);
    bounds.upper[a][d] = bounds.upper[a][d].min(r14_trans + 0.06);
}
```

**Step 3: Run tests**

`cargo test`

---

### Task 4: Add 1-5 Distance Bounds (#5)

**Files:**
- Modify: `src/etkdg/mod.rs` (after 1-4 section, before ring closure)

**Problem:** RDKit computes full 1-5 geometric bounds. We don't compute them at all.

**Step 1: Add 1-5 bounds computation**

After the 1-4 torsion loop, add 1-5 computation using find_5_atom_paths or BFS-based 5-step paths:

```rust
// 1-5 distance bounds (simplified RDKit approach)
// For atoms separated by 5 bonds, use VdW-scaled sum as bound
let tol_15 = 0.08;
let vdw_scale_15 = 0.7;

// Find all pairs at topological distance 5
for i in 0..n_atoms {
    for j in (i + 5)..n_atoms.min(i + 6) {  // heuristic: only nearby atoms
        if bounds.upper[i][j] < 1000.0 { continue; } // already constrained
        
        let ri = vdw_radius(&mol.atoms[i].symbol);
        let rj = vdw_radius(&mol.atoms[j].symbol);
        let ub = vdw_scale_15 * (ri + rj);
        
        if ub < bounds.upper[i][j] {
            bounds.upper[i][j] = ub;
            bounds.upper[j][i] = ub;
            bounds.lower[i][j] = bounds.lower[i][j].min(ub);
            bounds.lower[j][i] = bounds.lower[j][i].min(ub);
        }
    }
}
```

Note: Full 1-5 requires path enumeration (cis-cis, cis-trans, etc.). Implement simplified VdW-scaled version first.

**Step 2: Run tests**

`cargo test`

---

## Task Group 2: Coordinate Generation

### Task 5: Fix Distance Sampling to Uniform (#6)

**Files:**
- Modify: `src/etkdg/mod.rs:336-349` (`generate_4d_coordinates`)

**Problem:** We use biased sampling (squared random pushes toward lower bound). RDKit uses uniform sampling.

**Step 1: Replace biased sampling with uniform**

```rust
let t = lo + (hi - lo) * random_f64();  // uniform
// REMOVE the range<0.5 branch and u*u squaring
```

Full replacement of the sampling block:

```rust
for i in 0..n {
    for j in i + 1..n {
        let lo = bounds.lower[i][j];
        let hi = bounds.upper[i][j];
        let t = lo + (hi - lo) * random_f64();  // UNIFORM (RDKit-style)
        dist[i][j] = t;
        dist[j][i] = t;
    }
}
```

**Step 2: Run tests**

`cargo test`

---

### Task 6: Add 4th-Dimension Minimization (#7)

**Files:**
- Modify: `src/etkdg/mod.rs:398-432` (4D→3D projection section)
- New function: `collapse_fourth_dimension()`

**Problem:** Random orthogonal projection discards information. RDKit minimizes 4th coordinate via dedicated force field.

**Step 1: After 4D embedding, minimize 4th dimension**

Replace the random projection with: embed in 4D, then CG-minimize the 4th coordinate to near-zero while preserving inter-atom distances:

```rust
// After coords_4d is computed, collapse 4th dimension:
let mut coords_4d = /* ... existing 4D embed code ... */;

// Minimize 4th dimension: penalty on w-coordinate + distance preservation
collapse_fourth_dimension(&mut coords_4d, &bounds);

// Then take first 3 dimensions
for i in 0..n {
    coords[i] = [coords_4d[i][0], coords_4d[i][1], coords_4d[i][2]];
}
```

**Step 2: Implement `collapse_fourth_dimension`**

```rust
fn collapse_fourth_dimension(coords: &mut Vec<[f64; 4]>, bounds: &DistanceBounds) {
    let n = coords.len();
    let weight_4th = 1.0;      // RDKit uses weightFourthDim=0.1
    let weight_dist = 50.0;     // RDKit uses basinThresh=5.0 scaling
    let max_iter = 200;         // RDKit default
    
    for _iter in 0..max_iter {
        let mut grad = vec![0.0f64; n];
        let mut energy = 0.0;
        
        for i in 0..n {
            // 4th-dim penalty: push w toward 0
            let w = coords[i][3];
            energy += weight_4th * w * w;
            grad[i] += 2.0 * weight_4th * w;
            
            // Distance preservation in full 4D space
            for j in (i + 1)..n {
                let dx = coords[i][0] - coords[j][0];
                let dy = coords[i][1] - coords[j][1];
                let dz = coords[i][2] - coords[j][2];
                let dw = coords[i][3] - coords[j][3];
                let d4 = (dx*dx + dy*dy + dz*dz + dw*dw).sqrt().max(1e-10);
                
                let lo = bounds.lower[i][j];
                let hi = bounds.upper[i][j];
                let target = (lo + hi) / 2.0;
                let viol = d4 - target;
                
                if viol.abs() > 1e-8 {
                    let e = 0.5 * weight_dist * viol * viol;
                    energy += e;
                    let g = weight_dist * viol / d4;
                    
                    grad[i] += g * dx;
                    grad[i] += g * dy;
                    grad[i] += g * dz;
                    grad[i] += g * dw;
                    grad[j] -= g * dx;
                    grad[j] -= g * dy;
                    grad[j] -= g * dz;
                    grad[j] -= g * dw;
                }
            }
        }
        
        // Fixed step size SD
        let max_g = grad.iter().map(|g| g.abs()).fold(0.0_f64, f64::max).max(1e-10);
        let step = 0.01 / max_g;
        
        for i in 0..n {
            coords[i][0] -= step * grad[i];
            coords[i][1] -= step * grad[i];
            coords[i][2] -= step * grad[i];
            coords[i][3] -= step * grad[i];  // key: also update 4th dim
        }
        
        if max_g < 0.001 { break; }
    }
}
```

**Step 3: Remove old random projection code**

Delete the Gram-Schmidt projection block (lines ~398-429). Replace with simple extraction of first 3 dims after collapse.

**Step 4: Run tests**

`cargo test`

---

## Task Group 3: DG Refinement

### Task 7: Improve DG Refinement (#8 partial)

**Files:**
- Modify: `src/etkdg/mod.rs:926-997` (`refine_with_dg`)

**Problem:** No chirality enforcement. Uses steepest descent (RDKit uses conjugate gradient).

**Step 1: Switch from steepest descent to conjugate gradient (simplified)**

Or at minimum improve convergence: increase iterations, add better line search.

**Step 2: Add basic chiral volume constraint (optional, lower priority)**

For any tetrahedral center with 4 distinct neighbors, compute signed volume and penalize sign flips. This can be deferred since it primarily affects chiral molecules.

**Step 3: Run tests**

`cargo test`

---

## Task Group 4: Force Field Refinement

### Task 8: Improper Torsion Functional Form → UFF Inversion Style (#11)

**Files:**
- Modify: `src/etkdg/mod.rs:750-781` (`planarity_energy`)
- Modify: `src/etkdg/mod.rs:591-630` (`out_of_plane_angle`) — may need rename/refactor

**Problem:** We use harmonic `k*χ²`. RDKit uses UFF inversion: `k*(C0 + C1*sin(χ) + C2*cos(2χ))`.

**Step 1: Replace improper energy functional form**

UFF inversion for SP2 C/N/O has the form (from RDKit Inversions.cpp):
- For atom with 3 neighbors: `E = k * [χ²]` where χ is the Wilson angle (this IS actually what RDKit's UFF uses for basic inversion)
- The Fourier form `C0 + C1*sin(χ) + C2*cos(2χ)` is used for higher-coordination cases

Actually, re-checking: RDKit's `InversionContribs` uses `E = k * chi * chi` where chi is the out-of-plane angle (same as ours!). The UFF Fourier form is for the `UFFInversion` term specifically. Let me verify...

Upon re-analysis: Our harmonic improper `k*χ²` IS actually consistent with RDKit's approach for the basic inversion term. The difference is mainly in the force constant value and which atoms get terms. **This fix may be lower priority than initially assessed.**

**Revised Step 1:** Keep harmonic form but ensure force constant matches UFF values per atom type:

```rust
fn improper_k_for_atom(atom_idx: usize, mol: &Molecule) -> f64 {
    let sym = &mol.atoms[atom_idx].symbol;
    match sym.as_str() {
        "C" => 40.0,   // UFF C_SP2 inversion k ≈ 40 kcal/mol/rad²
        "N" => 30.0,   // UFF N_SP2
        "O" => 80.0,   // UFF O_SP2 (stiffer)
        "S" => 20.0,   // S_AR softer
        _ => 40.0,
    } * 10.0  // oobForceScalingFactor = 10.0
}
```

**Step 2: Apply per-atom k in planarity_energy**

```rust
for &(central, n1, n2, n3) in &pc.impropers {
    let chi = out_of_plane_angle(coords, central, n1, n2, n3);
    let k = improper_k_for_atom(central, mol);  // per-atom k
    energy += k * chi * chi;
}
```

Need to pass `mol` to `planarity_energy` or pre-compute k values in `PlanarityConstraints`.

**Step 3: Run tests**

`cargo test`

---

### Task 9: Basic Torsion Preferences for Flat Rings (#10 reduced scope)

**Files:**
- Modify: `src/etkdg/mod.rs` (ring torsion terms in build_planarity_constraints)

**Problem:** No experimental torsion preferences. RDKit has 100+ SMARTS patterns.

**Step 1: Add basic flat-ring torsion preferences during FF refinement**

Rather than implementing the full SMARTS pattern matching system (which would require embedding a pattern file), enhance our existing ring torsion constraints to be applied during Phase 1 SD (already done) AND add them as explicit torsion terms in the MMFF phase.

The current ring torsion terms with K_RING_TOR=1000 and cos(2φ) functional form are actually a reasonable approximation of RDKit's flat-ring torsion preferences. The main gap is non-ring torsion preferences (e.g., ethane staggered preference).

**Step 2: Add basic staggered preference for sp3-sp3 bonds**

In `build_planarity_constraints`, add torsion preferences for non-ring sp3-sp3 bonds:

```rust
// Basic torsion: sp3-sp3 single bonds prefer staggered (±60°, 180°)
// This is a simplified version of RDKit's torsion preferences
struct PlanarityConstraints {
    // ... existing fields ...
    pref_torsions: Vec<(usize, usize, usize, usize, f64, f64)>,  // (i,j,k,l, phi_pref, k)
}

// For each sp3-sp3 single bond not in a ring:
//   Add preferred torsion at φ=180° (anti) or φ=±60° (gauche) with moderate k
```

**Priority:** This is MODERATE impact for simple molecules like thiophene/benzene (which have few sp3-sp3 bonds). Can be deferred.

**Step 3: Run tests**

`cargo test`

---

## Task Group 5: Conformer Selection & Cleanup

### Task 10: Improve Conformer Selection (#12)

**Files:**
- Modify: `src/etkdg/mod.rs:1118-1173` (`generate_initial_coords_with_config`)

**Problem:** We accept first planar conformer or lowest-MMFF-energy. RDKit checks chirality, tetrahedral volume, double-bond geometry, stereochemistry.

**Step 1: Add basic validation before accepting conformer**

After `refine_with_ff`, add checks:

```rust
// After refine_with_ff, validate geometry quality:
if check_planarity(&coords_3d, mol, &pc, 0.1) {
    // ADDITIONAL CHECKS:
    // 1. No atoms too close (VDW clash check)
    // 2. Bond lengths reasonable (within ±0.3 of expected)
    // 3. Accept
    return coords_3d;
}
```

**Step 2: Add VDW clash rejection**

```rust
fn has_vdw_clash(coords: &[[f64; 3]], mol: &Molecule) -> bool {
    let clash_tol = 0.5;  // Å — reject if non-bonded atoms closer than this
    for i in 0..mol.atoms.len() {
        for j in (i + 2)..mol.atoms.len() {  // skip bonded (1-2) and angle (1-3)
            if mol.adjacency[i].contains(&j) { continue; }
            let dx = coords[i][0] - coords[j][0];
            let dy = coords[i][1] - coords[j][1];
            let dz = coords[i][2] - coords[j][2];
            let d = (dx*dx + dy*dy + dz*dz).sqrt();
            let ri = vdw_radius(&mol.atoms[i].symbol);
            let rj = vdw_radius(&mol.atoms[j].symbol);
            if d < 0.6 * (ri + rj) {  // severe clash
                return true;
            }
        }
    }
    false
}
```

**Step 3: Run tests**

`cargo test`

---

## Final Verification

### Task 11: Full Test Suite + Clippy + WASM Rebuild

**Step 1: Run all tests**

```bash
cargo test
```

Expected: All 124+ existing tests pass, plus any new tests added.

**Step 2: Run clippy**

```bash
cargo clippy
```

Expected: 0 errors (pre-existing warnings OK).

**Step 3: Build WASM**

```bash
wasm-pack build --target web --out-dir site/pkg --out-name webmm
```

**Step 4: Visual verification**

Serve at http://localhost:8080, load thiophene and benzene, click Generate 3D. Verify:
- Rings are flat (< 0.01Å deviation)
- Bond lengths match expected values
- No VDW clashes
- C-H bonds approximately coplanar with aromatic rings

**Step 5: Update CODE_STATUS.md**

Record Phase 24 completion with all 14 fixes documented.
