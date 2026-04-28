# M6-M12 Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Implement remaining M6-M12 gaps: SP3D/SP3D2 MMFF typing, atropisomer chirality, fragment embedding, coordMap constraints, timeout, et_version usage, and triangle smoothing epsilon.

**Architecture:** Extend existing ETKDG and MMFF modules with minimal, focused additions. Each feature is self-contained.

**Tech Stack:** Rust, WASM-compatible (no_std for core math, web-time for timeout)

---

## Task 1: M12 — Triangle Smoothing Convergence Epsilon

**Files:**
- Modify: `src/etkdg/mod.rs:186-216` (`smooth_triangle_inequality`)

**Step 1: Make epsilon configurable**

```rust
impl DistanceBounds {
    pub fn smooth_triangle_inequality(&mut self, epsilon: f64) {
        let mut changed = true;
        let n = self.n_atoms;
        while changed {
            changed = false;
            for k in 0..n {
                for i in 0..n {
                    for j in (i + 1)..n {
                        let new_lower = (self.lower[i][k] - self.lower[k][j]).abs();
                        if new_lower > self.lower[i][j] + epsilon {
                            self.lower[i][j] = new_lower;
                            self.lower[j][i] = new_lower;
                            changed = true;
                        }
                        let new_upper = self.upper[i][k] + self.upper[k][j];
                        if new_upper < self.upper[i][j] - epsilon {
                            self.upper[i][j] = new_upper;
                            self.upper[j][i] = new_upper;
                            changed = true;
                        }
                        if self.lower[i][j] > self.upper[i][j] {
                            self.lower[i][j] = self.upper[i][j];
                            self.lower[j][i] = self.upper[i][j];
                            changed = true;
                        }
                    }
                }
            }
        }
    }
    
    // Convenience method with default epsilon (matching RDKit)
    pub fn smooth_triangle_inequality_default(&mut self) {
        self.smooth_triangle_inequality(1e-6);
    }
}
```

**Step 2: Update call site**

```rust
// In generate_initial_coords_with_config
bounds.smooth_triangle_inequality(1e-6);
```

**Step 3: Run tests**

```bash
cargo test --lib
```

Expected: All tests pass

---

## Task 2: M11 — et_version Field Usage

**Files:**
- Modify: `src/etkdg/mod.rs:2978-3090` (`build_torsion_preferences`)
- Modify: `src/etkdg/mod.rs:3569` (call site)

**Step 1: Add et_version parameter**

```rust
fn build_torsion_preferences(mol: &Molecule, et_version: u32) -> Vec<TorsionPreference> {
    // ... existing logic ...
    
    // ETversion-specific parameter adjustments
    let sp3_sp3_v3 = if et_version >= 2 { 5.0 } else { 7.0 };
    let sp3_sp3_ring_v3 = if et_version >= 2 { 1.5 } else { 2.0 };
    
    // In sp3-sp3 general pattern:
    if sp3(h2) && sp3(h3) {
        if is_ch2(mol, a2) && is_ch2(mol, a3) {
            return Some(([1, 0, 1, 1, 1, 1], [0.0, 0.0, sp3_sp3_v3 - 3.0, 0.0, 0.0, 0.0]));
        }
        let a2_in_ring = rings.iter().any(|r| r.contains(&a2));
        let a3_in_ring = rings.iter().any(|r| r.contains(&a3));
        if a2_in_ring || a3_in_ring {
            return Some(([1, 0, 1, 1, 1, 1], [0.0, 0.0, sp3_sp3_ring_v3, 0.0, 0.0, 0.0]));
        }
        return Some(([1, 0, 1, 1, 1, 1], [0.0, 0.0, sp3_sp3_v3, 0.0, 0.0, 0.0]));
    }
    
    // ... rest of patterns ...
}
```

**Step 2: Update call site**

```rust
let torsion_prefs = build_torsion_preferences(mol, config.et_version);
```

**Step 3: Run tests**

```bash
cargo test --lib
```

---

## Task 3: M10 — Timeout Support

**Files:**
- Modify: `src/etkdg/mod.rs:96-126` (`ETKDGConfig`)
- Modify: `src/etkdg/mod.rs:3575-3602` (`generate_initial_coords_with_config`)
- Modify: `Cargo.toml` (add web-time)

**Step 1: Add web-time dependency**

In `Cargo.toml`:
```toml
[dependencies]
web-time = "1.1"
```

**Step 2: Add timeout field to config**

```rust
pub struct ETKDGConfig {
    pub max_attempts: usize,
    pub convergence_threshold: f64,
    pub max_iterations: usize,
    pub vdw_scale: f64,
    pub random_seed: i64,
    pub force_trans_amides: bool,
    pub use_macrocycle_14config: bool,
    pub use_small_ring_torsions: bool,
    pub use_macrocycle_torsions: bool,
    pub et_version: u32,
    pub timeout_ms: u64, // NEW: 0 = no timeout
}

impl Default for ETKDGConfig {
    fn default() -> Self {
        Self {
            max_attempts: 0,
            convergence_threshold: 1e-6,
            max_iterations: 200,
            vdw_scale: 0.8,
            random_seed: -1,
            force_trans_amides: true,
            use_macrocycle_14config: true,
            use_small_ring_torsions: false,
            use_macrocycle_torsions: true,
            et_version: 2,
            timeout_ms: 0, // Default: no timeout
        }
    }
}
```

**Step 3: Add timeout checking in embedding loop**

```rust
use web_time::{Instant, Duration};

pub fn generate_initial_coords_with_config(mol: &Molecule, config: &ETKDGConfig) -> Vec<[f64; 3]> {
    // ... existing setup ...
    
    let start_time = Instant::now();
    let timeout = if config.timeout_ms > 0 {
        Some(Duration::from_millis(config.timeout_ms))
    } else {
        None
    };
    
    for _attempt in 0..max_attempts {
        // Check timeout
        if let Some(t) = timeout {
            if start_time.elapsed() > t {
                eprintln!("ETKDG timeout after {} attempts", _attempt);
                break;
            }
        }
        
        // ... rest of embedding logic ...
    }
    
    // ... return best coords ...
}
```

**Step 4: Run tests**

```bash
cargo test --lib
```

---

## Task 4: M6 — SP3D/SP3D2 in MMFF Atom Typing

**Files:**
- Modify: `src/mmff/mod.rs:52-122` (`MMFFAtomType` enum)
- Modify: `src/mmff/mod.rs:193-418` (`assign_atom_types`)
- Modify: `src/mmff/angle.rs` (add angle params)
- Test: `src/lib.rs`

**Step 1: Add new atom types**

```rust
pub enum MMFFAtomType {
    // ... existing types ...
    
    // SP3D (trigonal bipyramidal)
    P_3D,   // Phosphorus in trigonal bipyramidal geometry (e.g., PF5)
    S_3D,   // Sulfur in trigonal bipyramidal geometry
    
    // SP3D2 (octahedral)
    S_3D2,  // Sulfur in octahedral geometry (e.g., SF6)
}
```

**Step 2: Add to base_type function**

```rust
pub fn base_type(t: MMFFAtomType) -> MMFFAtomType {
    match t {
        // ... existing mappings ...
        MMFFAtomType::P_3D => MMFFAtomType::P_3,
        MMFFAtomType::S_3D => MMFFAtomType::S_3,
        MMFFAtomType::S_3D2 => MMFFAtomType::S_3,
        other => other,
    }
}
```

**Step 3: Add detection in assign_atom_types**

```rust
(15, Hybridization::Sp3D, _, 5) => MMFFAtomType::P_3D,
(16, Hybridization::Sp3D, _, 4..=5) => MMFFAtomType::S_3D,
(16, Hybridization::Sp3D2, _, 6) => MMFFAtomType::S_3D2,
```

**Step 4: Add angle parameters**

In `src/mmff/angle.rs`:
```rust
// SP3D angles (trigonal bipyramidal)
(MMFFAtomType::F, MMFFAtomType::P_3D, MMFFAtomType::F) => Some(AngleParams { 
    k_theta: 0.80, 
    theta0: 180.0 
}),
(MMFFAtomType::F, MMFFAtomType::P_3D, MMFFAtomType::C_3) => Some(AngleParams { 
    k_theta: 0.80, 
    theta0: 90.0 
}),

// SP3D2 angles (octahedral)
(MMFFAtomType::F, MMFFAtomType::S_3D2, MMFFAtomType::F) => Some(AngleParams { 
    k_theta: 0.80, 
    theta0: 90.0 
}),
```

**Step 5: Add test**

```rust
#[test]
fn test_sp3d_typing() {
    let sdf = r#"PF5
     RDKit          3D
  6  5  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 P   0  0  0  0  0  0  0  0  0
    1.5000    0.0000    0.0000 F   0  0  0  0  0  0  0  0  0
   -1.5000    0.0000    0.0000 F   0  0  0  0  0  0  0  0  0
    0.0000    1.5000    0.0000 F   0  0  0  0  0  0  0  0  0
    0.0000   -1.5000    0.0000 F   0  0  0  0  0  0  0  0  0
    0.0000    0.0000    1.8000 F   0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  1  3  1  0  0  0  0
  1  4  1  0  0  0  0
  1  5  1  0  0  0  0
  1  6  1  0  0  0  0
M  END"#;
    let mol = parse_sdf(sdf).expect("parse failed");
    let ff = MMFFForceField::new(&mol, MMFFVariant::MMFF94s);
    assert_eq!(ff.atom_types[0], MMFFAtomType::P_3D);
    assert!(ff.angles.len() > 0);
}
```

**Step 6: Run tests**

```bash
cargo test test_sp3d_typing -- --nocapture
```

---

## Task 5: M9 — CoordMap (Constrained Atoms)

**Files:**
- Modify: `src/etkdg/mod.rs:96-126` (`ETKDGConfig`)
- Modify: `src/etkdg/mod.rs:3108-3239` (`etkdg_energy`, `etkdg_gradient`)
- Modify: `src/etkdg/mod.rs:3441-3509` (`minimize_etkdg`)

**Step 1: Add coord_map to config**

```rust
use std::collections::HashMap;

pub struct ETKDGConfig {
    // ... existing fields ...
    pub coord_map: HashMap<usize, [f64; 3]>,
}

impl Default for ETKDGConfig {
    fn default() -> Self {
        Self {
            // ... existing defaults ...
            coord_map: HashMap::new(),
        }
    }
}
```

**Step 2: Fix coord_map atoms after embedding**

```rust
// After generating initial coords_4d
for (&atom_idx, &fixed_pos) in &config.coord_map {
    if atom_idx < coords_4d.len() {
        coords_4d[atom_idx] = [fixed_pos[0], fixed_pos[1], fixed_pos[2], 0.0];
    }
}
```

**Step 3: Skip gradient updates for fixed atoms in minimization**

```rust
fn minimize_etkdg(
    coords: &mut [[f64; 3]],
    bounds: &DistanceBounds,
    chiral_centers: &[ChiralCenter],
    tetrahedral: &[ChiralCenter],
    pc: &PlanarityConstraints,
    torsion_prefs: &[TorsionPreference],
    bonds_12: &[(usize, usize)],
    angles_13: &[(usize, usize, usize)],
    max_iter: usize,
    force_tol: f64,
    coord_map: &HashMap<usize, [f64; 3]>, // NEW parameter
) -> f64 {
    // ... existing logic ...
    for _ in 0..max_iter {
        let grad = etkdg_gradient(...);
        let max_g = grad.iter().map(|g| ...).fold(0.0f64, f64::max);
        if max_g < force_tol { break; }
        let step = 0.1 / max_g.max(1e-10);
        for i in 0..n {
            if !coord_map.contains_key(&i) {
                coords[i][0] -= step * grad[i][0];
                coords[i][1] -= step * grad[i][1];
                coords[i][2] -= step * grad[i][2];
            }
        }
        // ... energy calculation ...
    }
}
```

**Step 4: Add WASM API**

In `src/lib.rs`:
```rust
#[wasm_bindgen]
pub fn generate_initial_coordinates_wasm_with_constraints(
    sdf_content: &str,
    atom_indices: Vec<usize>,
    positions: Vec<f64>,
) -> Result<ETKDGResult, JsValue> {
    let mut config = ETKDGConfig::default();
    for (i, &idx) in atom_indices.iter().enumerate() {
        config.coord_map.insert(idx, [
            positions[i * 3],
            positions[i * 3 + 1],
            positions[i * 3 + 2],
        ]);
    }
    
    let mol = match crate::molecule::parser::parse_sdf(sdf_content.trim()) {
        Ok(m) => m,
        Err(e) => return Ok(ETKDGResult {
            coordinates: Vec::new(),
            n_atoms: 0,
            success: false,
            error: format!("Parse error: {}", e),
        }),
    };
    
    let coords = crate::etkdg::generate_initial_coords_with_config(&mol, &config
    );
    // ... rest same as existing ...
}
```

---

## Task 6: M8 — Fragment Embedding

**Files:**
- Modify: `src/etkdg/mod.rs` (add helper functions)
- Modify: `src/etkdg/mod.rs:3569-3714` (`generate_initial_coords_with_config`)

**Step 1: Add connected component detection**

```rust
fn find_connected_components(mol: &Molecule) -> Vec<Vec<usize>> {
    let n = mol.atoms.len();
    let mut visited = vec![false; n];
    let mut components = Vec::new();
    
    for start in 0..n {
        if visited[start] { continue; }
        let mut comp = Vec::new();
        let mut stack = vec![start];
        visited[start] = true;
        
        while let Some(atom) = stack.pop() {
            comp.push(atom);
            for &nbr in &mol.adjacency[atom] {
                if !visited[nbr] {
                    visited[nbr] = true;
                    stack.push(nbr);
                }
            }
        }
        components.push(comp);
    }
    components
}
```

**Step 2: Add fragment extraction**

```rust
fn extract_submol(mol: &Molecule, atoms: &[usize]) -> Molecule {
    let mut atom_map = std::collections::HashMap::new();
    let mut new_atoms = Vec::new();
    for (new_idx, &old_idx) in atoms.iter().enumerate() {
        atom_map.insert(old_idx, new_idx);
        new_atoms.push(mol.atoms[old_idx].clone());
    }
    
    let mut new_bonds = Vec::new();
    for bond in &mol.bonds {
        if let (Some(&a1), Some(&a2)) = (atom_map.get(&bond.atom1), atom_map.get(&bond.atom2)) {
            let mut new_bond = bond.clone();
            new_bond.atom1 = a1;
            new_bond.atom2 = a2;
            new_bonds.push(new_bond);
        }
    }
    
    Molecule {
        atoms: new_atoms,
        bonds: new_bonds,
        adjacency: crate::molecule::graph::build_adjacency(new_atoms.len(), &new_bonds),
    }
}
```

**Step 3: Add fragment spreading**

```rust
fn spread_fragments(coords: &mut [[f64; 3]], components: &[Vec<usize>]) {
    let mut offset_x = 0.0;
    for comp in components {
        let xs: Vec<f64> = comp.iter().map(|&i| coords[i][0]).collect();
        let min_x = xs.iter().fold(f64::INFINITY, |a, &b| a.min(b));
        let max_x = xs.iter().fold(f64::NEG_INFINITY, |a, &b| a.max(b));
        let width = max_x - min_x;
        
        for &atom_idx in comp {
            coords[atom_idx][0] += offset_x - min_x;
        }
        offset_x += width + 5.0; // 5Å gap
    }
}
```

**Step 4: Modify generate_initial_coords_with_config**

```rust
pub fn generate_initial_coords_with_config(mol: &Molecule, config: &ETKDGConfig) -> Vec<[f64; 3]> {
    if mol.atoms.is_empty() {
        return Vec::new();
    }
    
    let components = find_connected_components(mol);
    if components.len() > 1 {
        // Multi-fragment molecule: embed each separately
        let mut all_coords = vec![[0.0; 3]; mol.atoms.len()];
        
        for comp in &components {
            let submol = extract_submol(mol, comp);
            let subcoords = generate_initial_coords_with_config(&submol, config
            );
            for (i, &atom_idx) in comp.iter().enumerate() {
                if i < subcoords.len() {
                    all_coords[atom_idx] = subcoords[i];
                }
            }
        }
        
        spread_fragments(&mut all_coords, &components);
        return all_coords;
    }
    
    // ... existing single-fragment logic ...
}
```

---

## Task 7: M7 — Atropisomer Chirality

**Files:**
- Modify: `src/etkdg/mod.rs:2067-2110` (extend `find_stereo_double_bonds`)
- Modify: `src/etkdg/mod.rs:3569` (call site)

**Step 1: Add Atropisomer struct**

```rust
#[derive(Debug, Clone)]
struct Atropisomer {
    bond: (usize, usize),
    substituents: (usize, usize, usize, usize),
    sign: f64,
    vol_lower: f64,
    vol_upper: f64,
}
```

**Step 2: Add detection in find_stereo_bonds**

```rust
fn find_stereo_bonds(mol: &Molecule) -> (Vec<(usize, usize, usize)>, Vec<StereoDoubleBond>, Vec<Atropisomer>) {
    let mut double_bond_ends = Vec::new();
    let mut stereo_db = Vec::new();
    let mut atropisomers = Vec::new();
    
    for bond in &mol.bonds {
        // ... existing double bond logic ...
        
        // Atropisomer detection
        if bond.stereo == BondStereo::AtropCW || bond.stereo == BondStereo::AtropCCW {
            let a1 = bond.atom1;
            let a2 = bond.atom2;
            let nbrs1: Vec<usize> = mol.adjacency[a1].iter()
                .filter(|&&n| n != a2).copied().collect();
            let nbrs2: Vec<usize> = mol.adjacency[a2].iter()
                .filter(|&&n| n != a1).copied().collect();
            
            if nbrs1.len() >= 2 && nbrs2.len() >= 2 {
                let (vol_lower, vol_upper) = match bond.stereo {
                    BondStereo::AtropCW => (-100.0, -1.0),
                    BondStereo::AtropCCW => (1.0, 100.0),
                    _ => (0.0, 0.0),
                };
                atropisomers.push(Atropisomer {
                    bond: (a1, a2),
                    substituents: (nbrs1[0], nbrs1[1], nbrs2[0], nbrs2[1]),
                    sign: if bond.stereo == BondStereo::AtropCW { -1.0 } else { 1.0 },
                    vol_lower,
                    vol_upper,
                });
            }
        }
    }
    
    (double_bond_ends, stereo_db, atropisomers)
}
```

**Step 3: Add atropisomer validation**

```rust
fn check_atropisomers(coords: &[[f64; 3]], atropisomers: &[Atropisomer]) -> bool {
    for atrop in atropisomers {
        let (a, b, c, d) = atrop.substituents;
        let vol = chiral_volume(coords, a, b, c, d);
        if vol < atrop.vol_lower || vol > atrop.vol_upper {
            return false;
        }
    }
    true
}
```

**Step 4: Update validation pipeline**

```rust
let (double_bond_ends, stereo_dbs, atropisomers) = find_stereo_bonds(mol);

// In validation:
let atrop_ok = atropisomers.is_empty() || check_atropisomers(&coords_3d, &atropisomers);

if planar && db_geom_ok && chiral_ok && db_stereo_ok && atrop_ok && no_clash && bonds_ok {
    // ... accept conformer ...
}
```

---

## Task 8: Numerical Gradients Analysis

**Current Status:**
- **Bond**: Analytical ✓
- **Angle**: Analytical ✓  
- **Stretch-bend**: Analytical ✓
- **Torsion**: Numerical (finite difference, 12 evals per term)
- **OOP**: Numerical (finite difference, 12 evals per term)
- **VDW**: Numerical (finite difference, 6 evals per pair)
- **Electrostatics**: Analytical ✓

**Why Numerical?**
1. Torsion: Complex chain rule through dihedral angle → cos(φ) → energy. Derivable but error-prone.
2. OOP: Wilson angle involves normalized vectors + cross products. Analytical derivation is messy.
3. VDW: Buffered 14-7 has complex dE/dr involving t^7 * b where both t and b depend on r/r*.

**Performance Impact:**
- Torsion: ~5-10% of total gradient time
- OOP: ~1-2% (few terms)
- VDW: ~10-20% (O(N²) pairs, each needs 6 evaluations)

**Decision: KEEP numerical for now.** They are:
- Correct (verified by tests against finite difference)
- Only bottleneck for very large molecules (>100 atoms)
- Much simpler to maintain than analytical derivations
- Can be optimized later if profiling shows they're the bottleneck

**Future:** Implement analytical torsion gradient (standard formula from molecular mechanics textbooks) and VDW gradient (derivative of buffered 14-7) if needed.

---

## Task 9: WASM Compatibility Checklist

**Dependencies:**
- [ ] Add `web-time = "1.1"` to Cargo.toml
- [ ] Replace `std::time::SystemTime` with `web_time::SystemTime` in ETKDGConfig::default() for seed generation
- [ ] Ensure no `std::time::Instant` usage (already using custom Rng)

**API Surface:**
- [ ] `generate_initial_coordinates_wasm` — existing, works
- [ ] `optimize_from_sdf` — existing, works
- [ ] `generate_initial_coordinates_wasm_with_constraints` — NEW (M9)

**Build Verification:**
```bash
wasm-pack build --target web --out-dir site/pkg --out-name webmm
```

**Expected:** Clean build with 0 errors.

---

## Task 10: Update CODE_STATUS.md

After all tasks complete, update:

```markdown
## Phase 26: M6-M12 Gap Fixes (Completed)
- **M6**: Added SP3D/SP3D2 MMFF atom types (P_3D, S_3D, S_3D2)
- **M7**: Added atropisomer chirality detection with volume checks
- **M8**: Added fragment embedding for multi-component molecules
- **M9**: Added coordMap support for constrained atom positions
- **M10**: Added configurable timeout with web-time (WASM-compatible)
- **M11**: Added et_version usage in torsion parameter selection
- **M12**: Added configurable epsilon for triangle smoothing (default 1e-6)
- All tests pass, 0 clippy errors, WASM build verified
```

---

## Execution Order

1. **M12** (simplest — configurable epsilon)
2. **M11** (et_version usage)
3. **M10** (timeout + web-time dependency)
4. **M6** (SP3D/SP3D2 typing)
5. **M9** (coordMap)
6. **M8** (fragment embedding)
7. **M7** (atropisomer)
8. **WASM verification**
9. **Update CODE_STATUS.md**

Commit after each task.
