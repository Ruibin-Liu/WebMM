//! ETKDG v3 3D coordinate embedding
//!
//! This implementation follows RDKit's ETKDGv3 algorithm:
//! 1. Build distance bounds matrix from topology
//! 2. Triangle-smooth bounds
//! 3. For each attempt:
//!    a. Sample random distance matrix from bounds
//!    b. Compute initial coordinates via metric matrix embedding (4D)
//!    c. First minimization: distance bounds + chirality in 4D
//!    d. Check tetrahedral centers
//!    e. Check chiral centers
//!    f. Minimize 4th dimension (collapse w)
//!    g. Minimize with experimental torsions + improper torsions + distance constraints in 3D
//!    h. Check planarity
//!    i. Double bond geometry checks
//!    j. Final chiral checks
//!    k. Double bond stereo checks

// Index-based loops, multi-argument kernels, and large table types are
// pervasive in this distance-geometry numerical code; the indexed form is
// clearer for matrix/coordinate access, so these lints are allowed here.
#![allow(clippy::needless_range_loop, clippy::too_many_arguments, clippy::type_complexity)]

use crate::molecule::{BondStereo, BondType, Hybridization, Molecule};
use std::collections::HashSet;
use web_time::Instant;

// ============================================================================
// Constants (matching RDKit)
// ============================================================================

const DIST12_DELTA: f64 = 0.01;
const DIST13_TOL: f64 = 0.04;
const GEN_DIST_TOL: f64 = 0.06;
const DIST15_TOL: f64 = 0.08;
const VDW_SCALE_15: f64 = 0.7;
const H_BOND_LENGTH: f64 = 1.8;
const MAX_UPPER: f64 = 1000.0;

const MIN_TETRAHEDRAL_CHIRAL_VOL: f64 = 0.50;
const TETRAHEDRAL_CENTERINVOLUME_TOL: f64 = 0.30;
const CHIRAL_CENTERINVOLUME_TOL: f64 = 0.10;
const MAX_MINIMIZED_E_PER_ATOM: f64 = 0.05;
const BASIN_THRESH: f64 = 5.0;
const MIN_MACROCYCLE_SIZE: usize = 9;
const EXTRA_SQUISH: f64 = 0.2;
const LONG_RANGE_FORCE: f64 = 10.0;
const FIRST_MIN_WEIGHT_CHIRAL: f64 = 1.0;
const FIRST_MIN_WEIGHT_FOURTH: f64 = 0.1;
const FOURTH_MIN_WEIGHT_CHIRAL: f64 = 0.2;
const FOURTH_MIN_WEIGHT_FOURTH: f64 = 1.0;
const PLANARITY_ENERGY_TOL: f64 = 0.7;

// M2 bundle: master RDKit-faithful toggle. Read throughout the shared embed path
// (bounds, minimizer, 3D force field, workarounds). `generate_initial_coords_default`
// sets it false (shipped path byte-identical); `generate_initial_coords_rdkit` sets
// it true. Faithfulness is the goal — r-regression during the migration is
// expected and accepted, not a gate.
static EXP_RDKIT_ALL: std::sync::atomic::AtomicBool =
    std::sync::atomic::AtomicBool::new(false);
fn rdkit_all() -> bool {
    EXP_RDKIT_ALL.load(std::sync::atomic::Ordering::Relaxed)
}

// ============================================================================
// Seeded Random Number Generator
// ============================================================================

#[derive(Debug, Clone, Copy)]
struct Rng {
    s: [u64; 4],
}

impl Rng {
    fn new(seed: u64) -> Self {
        let mut z = seed.wrapping_add(0x9e3779b97f4a7c15);
        let mut s = [0u64; 4];
        for i in 0..4 {
            z = z.wrapping_add(0x9e3779b97f4a7c15);
            let mut x = z;
            x = (x ^ (x >> 30)).wrapping_mul(0xbf58476d1ce4e5b9);
            x = (x ^ (x >> 27)).wrapping_mul(0x94d049bb133111eb);
            s[i] = x ^ (x >> 31);
        }
        Self { s }
    }

    fn next_u64(&mut self) -> u64 {
        let result = self.rotl(self.s[1].wrapping_mul(5), 7).wrapping_mul(9);
        let t = self.s[1] << 17;
        self.s[2] ^= self.s[0];
        self.s[3] ^= self.s[1];
        self.s[1] ^= self.s[2];
        self.s[0] ^= self.s[3];
        self.s[2] ^= t;
        self.s[3] = self.rotl(self.s[3], 45);
        result
    }

    fn rotl(&self, x: u64, k: i32) -> u64 {
        x.rotate_left(k as u32)
    }

    fn random_f64(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 * (1.0 / ((1u64 << 53) as f64))
    }
}

// ============================================================================
// ETKDG Configuration
// ============================================================================

#[derive(Debug, Clone)]
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
    pub triangle_smoothing_epsilon: f64,
    pub timeout_ms: u64,
    pub coord_map: std::collections::HashMap<usize, [f64; 3]>,
    /// M0+: when true, use the RDKit-faithful embedding path (bounds + 4D L-BFGS +
    /// 3D force field) instead of the tuned default. Currently a passthrough.
    pub rdkit_faithful: bool,
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
            triangle_smoothing_epsilon: 1e-6,
            timeout_ms: 0,
            coord_map: std::collections::HashMap::new(),
            rdkit_faithful: false,
        }
    }
}

// ============================================================================
// Distance Bounds Matrix
// ============================================================================

#[derive(Debug, Clone)]
pub struct DistanceBounds {
    pub lower: Vec<Vec<f64>>,
    pub upper: Vec<Vec<f64>>,
    pub n_atoms: usize,
}

impl DistanceBounds {
    pub fn new(n_atoms: usize) -> Self {
        let mut lower = vec![vec![0.0; n_atoms]; n_atoms];
        let mut upper = vec![vec![MAX_UPPER; n_atoms]; n_atoms];
        for i in 0..n_atoms {
            lower[i][i] = 0.0;
            upper[i][i] = 0.0;
        }
        Self {
            lower,
            upper,
            n_atoms,
        }
    }

    fn set_lower(&mut self, i: usize, j: usize, val: f64) {
        let (a, b) = if i < j { (i, j) } else { (j, i) };
        if val > DIST12_DELTA && (self.lower[a][b] < DIST12_DELTA || val < self.lower[a][b]) {
            self.lower[a][b] = val;
            self.lower[b][a] = val;
        }
    }

    fn set_upper(&mut self, i: usize, j: usize, val: f64) {
        let (a, b) = if i < j { (i, j) } else { (j, i) };
        if val < MAX_UPPER && (self.upper[a][b] >= MAX_UPPER || val > self.upper[a][b]) {
            self.upper[a][b] = val;
            self.upper[b][a] = val;
        }
    }

    fn check_and_set(&mut self, i: usize, j: usize, lb: f64, ub: f64) {
        assert!(ub > lb, "upper bound not greater than lower bound");
        if rdkit_all() {
            // RDKit _checkAndSetBounds: conservative/widening (union of ranges).
            // set_lower/set_upper already widen monotonically (take smaller lower /
            // larger upper, with the DIST12_DELTA / MAX_UPPER guards RDKit applies).
            self.set_lower(i, j, lb);
            self.set_upper(i, j, ub);
            return;
        }
        let clb = self.lower[i.min(j)][i.max(j)];
        let cub = self.upper[i.min(j)][i.max(j)];
        let nlb = clb.max(lb);
        let nub = cub.min(ub);
        if nub <= nlb {
            let nlb2 = clb.min(lb);
            let nub2 = cub.max(ub);
            self.set_lower(i, j, nlb2);
            self.set_upper(i, j, nub2);
        } else {
            self.set_lower(i, j, nlb);
            self.set_upper(i, j, nub);
        }
    }

    pub fn smooth_triangle_inequality(&mut self, epsilon: f64) {
        let n = self.n_atoms;
        if rdkit_all() {
            // RDKit triangleSmoothBounds: a single Floyd-Warshall sweep over k
            // (upper = shortest path, converges in one sweep) with the correct
            // lower bound d_ij >= max(L_ik - U_kj, L_kj - U_ik, 0).
            for k in 0..n {
                for i in 0..n {
                    if i == k {
                        continue;
                    }
                    for j in (i + 1)..n {
                        if j == k {
                            continue;
                        }
                        let new_lower = (self.lower[i][k] - self.upper[k][j])
                            .max(self.lower[k][j] - self.upper[i][k]);
                        if new_lower > self.lower[i][j] + epsilon {
                            self.lower[i][j] = new_lower;
                            self.lower[j][i] = new_lower;
                        }
                        let new_upper = self.upper[i][k] + self.upper[k][j];
                        if new_upper < self.upper[i][j] - epsilon {
                            self.upper[i][j] = new_upper;
                            self.upper[j][i] = new_upper;
                        }
                        if self.lower[i][j] > self.upper[i][j] {
                            self.lower[i][j] = self.upper[i][j];
                            self.lower[j][i] = self.upper[i][j];
                        }
                    }
                }
            }
            return;
        }
        let mut changed = true;
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

    pub fn get_lower(&self, i: usize, j: usize) -> f64 {
        self.lower[i.min(j)][i.max(j)]
    }

    pub fn get_upper(&self, i: usize, j: usize) -> f64 {
        self.upper[i.min(j)][i.max(j)]
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
enum Path14Type {
    Cis,
    Trans,
    Other,
}

#[derive(Debug, Clone)]
struct Path14Config {
    bid1: usize,
    bid2: usize,
    bid3: usize,
    ptype: Path14Type,
}

#[derive(Debug, Clone)]
struct ComputedData {
    bond_lengths: Vec<f64>,
    bond_angles: Vec<Vec<f64>>,
    bond_adj: Vec<Vec<i32>>,
    paths14: Vec<Path14Config>,
    cis_paths: std::collections::HashSet<u64>,
    trans_paths: std::collections::HashSet<u64>,
    visited12: Vec<bool>,
    visited13: Vec<bool>,
    visited14: Vec<bool>,
    set15: Vec<bool>,
}

impl ComputedData {
    fn new(n_atoms: usize, n_bonds: usize) -> Self {
        Self {
            bond_lengths: vec![0.0; n_bonds],
            bond_angles: vec![vec![-1.0; n_bonds]; n_bonds],
            bond_adj: vec![vec![-1i32; n_bonds]; n_bonds],
            paths14: Vec::new(),
            cis_paths: std::collections::HashSet::new(),
            trans_paths: std::collections::HashSet::new(),
            visited12: vec![false; n_atoms * n_atoms],
            visited13: vec![false; n_atoms * n_atoms],
            visited14: vec![false; n_atoms * n_atoms],
            set15: vec![false; n_atoms * n_atoms],
        }
    }

    fn visited_bound(&self, pid: usize, max_dist: u8) -> bool {
        (max_dist >= 1 && self.visited12[pid])
            || (max_dist >= 2 && self.visited13[pid])
            || (max_dist >= 3 && self.visited14[pid])
    }
}

// ============================================================================
// Radii
// ============================================================================

fn vdw_radius(element: &str) -> f64 {
    match element {
        "H" => 1.20,
        "C" => 1.70,
        "N" => 1.55,
        "O" => 1.52,
        "F" => 1.47,
        "P" => 1.80,
        "S" => 1.80,
        "Cl" => 1.75,
        "Br" => 1.85,
        "I" => 1.98,
        "Si" => 2.10,
        "B" => 1.80,
        _ => 1.70,
    }
}

fn uff_radius(element: &str, hyb: Hybridization, is_aromatic: bool) -> f64 {
    match element {
        "B" => match hyb {
            Hybridization::Sp2 => 0.828,
            _ => 0.838,
        },
        "C" => {
            if is_aromatic {
                0.729
            } else {
                match hyb {
                    Hybridization::Sp2 => 0.732,
                    Hybridization::Sp1 => 0.706,
                    _ => 0.757,
                }
            }
        }
        "N" => {
            if is_aromatic {
                0.699
            } else {
                match hyb {
                    Hybridization::Sp2 => 0.685,
                    Hybridization::Sp1 => 0.656,
                    _ => 0.700,
                }
            }
        }
        "O" => {
            if is_aromatic {
                0.680
            } else {
                match hyb {
                    Hybridization::Sp2 => 0.634,
                    Hybridization::Sp1 => 0.639,
                    _ => 0.658,
                }
            }
        }
        "S" => {
            if is_aromatic {
                1.077
            } else {
                match hyb {
                    Hybridization::Sp2 => 0.854,
                    _ => 1.064,
                }
            }
        }
        "H" => 0.354,
        "F" => 0.668,
        "Si" => 1.117,
        "P" => 1.101,
        "Cl" => 1.044,
        "Br" => 1.192,
        "I" => 1.382,
        _ => 0.757,
    }
}

fn uff_electronegativity(element: &str) -> f64 {
    match element {
        "H" => 4.528,
        "B" => 5.11,
        "C" => 5.343,
        "N" => 6.899,
        "O" => 8.741,
        "F" => 10.874,
        "Si" => 4.168,
        "P" => 5.463,
        "S" => 6.928,
        "Cl" => 8.564,
        "Br" => 7.79,
        "I" => 6.822,
        _ => 5.343,
    }
}

const UFF_LAMBDA: f64 = 0.1332;

#[allow(dead_code)]
fn covalent_radius(element: &str) -> f64 {
    match element {
        "H" => 0.31,
        "C" => 0.76,
        "N" => 0.71,
        "O" => 0.66,
        "F" => 0.57,
        "P" => 1.07,
        "S" => 1.05,
        "Cl" => 1.02,
        "Br" => 1.20,
        "I" => 1.33,
        "Si" => 1.11,
        "B" => 0.84,
        _ => 0.76,
    }
}

fn compute_bond_length(
    element1: &str,
    hyb1: Hybridization,
    is_aromatic1: bool,
    element2: &str,
    hyb2: Hybridization,
    is_aromatic2: bool,
    bond_type: BondType,
) -> f64 {
    // Specific overrides for bonds where the UFF formula is significantly off
    // (avoids systematic S=O / C=S errors in the distance bounds).
    if bond_type == BondType::Double {
        let pair = (element1, element2);
        match pair {
            ("S", "O") | ("O", "S") => return 1.44, // MMFF sulfonyl S=O
            ("C", "S") | ("S", "C") => return 1.56, // MMFF thiocarbonyl C=S
            _ => {}
        }
    }
    let ri = uff_radius(element1, hyb1, is_aromatic1);
    let rj = uff_radius(element2, hyb2, is_aromatic2);
    let bo: f64 = if is_aromatic1 && is_aromatic2 {
        1.5
    } else {
        match bond_type {
            BondType::Single => 1.0,
            BondType::Double => 2.0,
            BondType::Triple => 3.0,
            BondType::Aromatic => 1.5,
        }
    };
    let r_bo = -UFF_LAMBDA * (ri + rj) * bo.ln();
    let xi = uff_electronegativity(element1);
    let xj = uff_electronegativity(element2);
    let sqrt_xi = xi.sqrt();
    let sqrt_xj = xj.sqrt();
    let r_en = ri * rj * (sqrt_xi - sqrt_xj) * (sqrt_xi - sqrt_xj) / (xi * ri + xj * rj);
    ri + rj + r_bo - r_en
}

fn element_to_atomic_num(element: &str) -> u8 {
    match element {
        "H" => 1,
        "C" => 6,
        "N" => 7,
        "O" => 8,
        "F" => 9,
        "P" => 15,
        "S" => 16,
        "Cl" => 17,
        "Br" => 35,
        "I" => 53,
        "Si" => 14,
        "B" => 5,
        _ => 6,
    }
}

// ============================================================================
// Geometry helpers
// ============================================================================

fn set_ring_angle(hyb: Hybridization, ring_size: usize) -> f64 {
    if (matches!(hyb, Hybridization::Sp2) && ring_size <= 8) || ring_size == 3 || ring_size == 4 {
        std::f64::consts::PI * (1.0 - 2.0 / ring_size as f64)
    } else if matches!(hyb, Hybridization::Sp3) {
        if ring_size == 5 {
            104.0_f64.to_radians()
        } else {
            109.5_f64.to_radians()
        }
    } else if matches!(hyb, Hybridization::Sp3D) {
        105.0_f64.to_radians()
    } else if matches!(hyb, Hybridization::Sp3D2) {
        // RDKit _setRingAngle: SP3D2 -> 90 deg (octahedral cis ligands).
        90.0_f64.to_radians()
    } else {
        120.0_f64.to_radians()
    }
}

fn compute_13_dist(bl1: f64, bl2: f64, angle: f64) -> f64 {
    (bl1 * bl1 + bl2 * bl2 - 2.0 * bl1 * bl2 * angle.cos()).sqrt()
}

fn compute_14_dist_cis(bl1: f64, bl2: f64, bl3: f64, ba12: f64, ba23: f64) -> f64 {
    let r = bl1 * bl1 + bl2 * bl2 + bl3 * bl3
        - 2.0 * bl1 * bl2 * ba12.cos()
        - 2.0 * bl2 * bl3 * ba23.cos()
        + 2.0 * bl1 * bl3 * (ba12.cos() * ba23.cos() - ba12.sin() * ba23.sin());
    if r > 0.0 {
        r.sqrt()
    } else {
        0.5
    }
}

fn compute_14_dist_trans(bl1: f64, bl2: f64, bl3: f64, ba12: f64, ba23: f64) -> f64 {
    let r = bl1 * bl1 + bl2 * bl2 + bl3 * bl3
        - 2.0 * bl1 * bl2 * ba12.cos()
        - 2.0 * bl2 * bl3 * ba23.cos()
        + 2.0 * bl1 * bl3 * (ba12.cos() * ba23.cos() + ba12.sin() * ba23.sin());
    if r > 0.0 {
        r.sqrt()
    } else {
        0.5
    }
}

fn compute_14_dist_at_dihedral(
    bl1: f64,
    bl2: f64,
    bl3: f64,
    ba12: f64,
    ba23: f64,
    dihedral: f64,
) -> f64 {
    let r = bl1 * bl1 + bl2 * bl2 + bl3 * bl3
        - 2.0 * bl1 * bl2 * ba12.cos()
        - 2.0 * bl2 * bl3 * ba23.cos()
        + 2.0 * bl1 * bl3 * (ba12.cos() * ba23.cos() + ba12.sin() * ba23.sin() * dihedral.cos());
    if r > 0.0 {
        r.sqrt()
    } else {
        0.5
    }
}

fn compute_15_dist_cis_cis(
    d1: f64,
    d2: f64,
    d3: f64,
    d4: f64,
    ang12: f64,
    ang23: f64,
    ang34: f64,
) -> f64 {
    let dx14 = d2 - d3 * ang23.cos() - d1 * ang12.cos();
    let dy14 = d3 * ang23.sin() - d1 * ang12.sin();
    let d14 = (dx14 * dx14 + dy14 * dy14).sqrt();
    let mut cval = (d3 - d2 * ang23.cos() + d1 * (ang12 + ang23).cos()) / d14;
    cval = cval.clamp(-1.0, 1.0);
    let ang143 = cval.acos();
    let ang145 = ang34 - ang143;
    compute_13_dist(d14, d4, ang145)
}

fn compute_15_dist_cis_trans(
    d1: f64,
    d2: f64,
    d3: f64,
    d4: f64,
    ang12: f64,
    ang23: f64,
    ang34: f64,
) -> f64 {
    let dx14 = d2 - d3 * ang23.cos() - d1 * ang12.cos();
    let dy14 = d3 * ang23.sin() - d1 * ang12.sin();
    let d14 = (dx14 * dx14 + dy14 * dy14).sqrt();
    let mut cval = (d3 - d2 * ang23.cos() + d1 * (ang12 + ang23).cos()) / d14;
    cval = cval.clamp(-1.0, 1.0);
    let ang143 = cval.acos();
    let ang145 = ang34 + ang143;
    compute_13_dist(d14, d4, ang145)
}

fn compute_15_dist_trans_trans(
    d1: f64,
    d2: f64,
    d3: f64,
    d4: f64,
    ang12: f64,
    ang23: f64,
    ang34: f64,
) -> f64 {
    let dx14 = d2 - d3 * ang23.cos() - d1 * ang12.cos();
    let dy14 = d3 * ang23.sin() + d1 * ang12.sin();
    let d14 = (dx14 * dx14 + dy14 * dy14).sqrt();
    let mut cval = (d3 - d2 * ang23.cos() + d1 * (ang12 - ang23).cos()) / d14;
    cval = cval.clamp(-1.0, 1.0);
    let ang143 = cval.acos();
    let ang145 = ang34 + ang143;
    compute_13_dist(d14, d4, ang145)
}

fn compute_15_dist_trans_cis(
    d1: f64,
    d2: f64,
    d3: f64,
    d4: f64,
    ang12: f64,
    ang23: f64,
    ang34: f64,
) -> f64 {
    let dx14 = d2 - d3 * ang23.cos() - d1 * ang12.cos();
    let dy14 = d3 * ang23.sin() + d1 * ang12.sin();
    let d14 = (dx14 * dx14 + dy14 * dy14).sqrt();
    let mut cval = (d3 - d2 * ang23.cos() + d1 * (ang12 - ang23).cos()) / d14;
    cval = cval.clamp(-1.0, 1.0);
    let ang143 = cval.acos();
    let ang145 = ang34 - ang143;
    compute_13_dist(d14, d4, ang145)
}

fn find_bond_index(mol: &Molecule, i: usize, j: usize) -> Option<usize> {
    mol.bonds
        .iter()
        .enumerate()
        .find(|(_, b)| (b.atom1 == i && b.atom2 == j) || (b.atom1 == j && b.atom2 == i))
        .map(|(idx, _)| idx)
}

fn is_larger_sp2_atom(atom_idx: usize, mol: &Molecule, rings: &[Vec<usize>]) -> bool {
    let atom = &mol.atoms[atom_idx];
    let hyb = crate::molecule::graph::determine_hybridization(atom_idx, mol);
    let atomic_num = element_to_atomic_num(&atom.symbol);
    // RDKit isLargerSP2Atom: atomicNum > 13 && SP2 && numAtomRings(idx) > 0.
    let in_ring = rings.iter().any(|r| r.contains(&atom_idx));
    atomic_num > 13 && matches!(hyb, Hybridization::Sp2) && in_ring
}

fn is_amide_bond(mol: &Molecule, a: usize, b: usize) -> bool {
    let sym_a = &mol.atoms[a].symbol;
    let sym_b = &mol.atoms[b].symbol;
    if !((sym_a == "C" && sym_b == "N") || (sym_a == "N" && sym_b == "C")) {
        return false;
    }
    let (c_atom, _) = if sym_a == "C" { (a, b) } else { (b, a) };
    mol.adjacency[c_atom].iter().any(|&nbr| {
        mol.bonds.iter().any(|bond| {
            (bond.atom1 == c_atom && bond.atom2 == nbr && bond.bond_type == BondType::Double)
                || (bond.atom2 == c_atom && bond.atom1 == nbr && bond.bond_type == BondType::Double)
        })
    })
}

fn is_ester_bond(mol: &Molecule, a: usize, b: usize) -> bool {
    let sym_a = &mol.atoms[a].symbol;
    let sym_b = &mol.atoms[b].symbol;
    if !((sym_a == "C" && sym_b == "O") || (sym_a == "O" && sym_b == "C")) {
        return false;
    }
    let (c_atom, o_atom) = if sym_a == "C" { (a, b) } else { (b, a) };
    // Carbon must have a double bond to another oxygen (C=O)
    let has_carbonyl = mol.adjacency[c_atom].iter().any(|&nbr| {
        nbr != o_atom
            && mol.bonds.iter().any(|bond| {
                (bond.atom1 == c_atom && bond.atom2 == nbr && bond.bond_type == BondType::Double)
                    || (bond.atom2 == c_atom
                        && bond.atom1 == nbr
                        && bond.bond_type == BondType::Double)
            })
    });
    if !has_carbonyl {
        return false;
    }
    // Single bond oxygen must have no other carbonyl neighbors (distinguish from acid anhydride)
    mol.adjacency[o_atom].iter().all(|&nbr| {
        nbr == c_atom
            || !mol.bonds.iter().any(|bond| {
                (bond.atom1 == o_atom && bond.atom2 == nbr && bond.bond_type == BondType::Double)
                    || (bond.atom2 == o_atom
                        && bond.atom1 == nbr
                        && bond.bond_type == BondType::Double)
            })
    })
}

fn is_double_bond(mol: &Molecule, i: usize, j: usize) -> bool {
    mol.bonds.iter().any(|b| {
        ((b.atom1 == i && b.atom2 == j) || (b.atom1 == j && b.atom2 == i))
            && b.bond_type == BondType::Double
    })
}

fn find_non_double_neighbor(mol: &Molecule, atom: usize, exclude: usize) -> Option<usize> {
    mol.adjacency[atom]
        .iter()
        .copied()
        .find(|&n| n != exclude && !is_double_bond(mol, atom, n))
}

fn atom_degree(atom_idx: usize, mol: &Molecule) -> usize {
    mol.adjacency[atom_idx].len()
}

fn atom_num(atom_idx: usize, mol: &Molecule) -> u8 {
    element_to_atomic_num(&mol.atoms[atom_idx].symbol)
}

#[allow(dead_code)]
fn is_carbonyl(mol: &Molecule, atom_idx: usize) -> bool {
    if atom_num(atom_idx, mol) != 6 {
        return false;
    }
    mol.adjacency[atom_idx].iter().any(|&nbr| {
        (atom_num(nbr, mol) == 8 || atom_num(nbr, mol) == 7) && is_double_bond(mol, atom_idx, nbr)
    })
}

fn total_h_count(atom_idx: usize, mol: &Molecule) -> u8 {
    mol.adjacency[atom_idx]
        .iter()
        .filter(|&&nbr| mol.atoms[nbr].symbol == "H")
        .count() as u8
}

#[allow(dead_code)]
fn bonds_share_ring(bid1: usize, bid2: usize, bond_rings: &[Vec<usize>]) -> bool {
    let b1 = if bid1 < bid2 { bid1 } else { bid2 };
    let b2 = if bid1 < bid2 { bid2 } else { bid1 };
    bond_rings
        .iter()
        .any(|ring| ring.contains(&b1) && ring.contains(&b2))
}

#[allow(dead_code)]
fn is_bond_in_ring(bid: usize, bond_rings: &[Vec<usize>]) -> bool {
    bond_rings.iter().any(|ring| ring.contains(&bid))
}

fn check_amide_ester_14(
    mol: &Molecule,
    aid1: usize,
    aid2: usize,
    aid3: usize,
    aid4: usize,
    bid1: usize,
    bid3: usize,
) -> bool {
    let a2_num = atom_num(aid2, mol);
    let a3_num = atom_num(aid3, mol);
    let a4_num = atom_num(aid4, mol);
    let b1_type = mol.bonds[bid1].bond_type;
    let b3_type = mol.bonds[bid3].bond_type;
    if a3_num == 6
        && b3_type == BondType::Double
        && (a4_num == 8 || a4_num == 7)
        && b1_type == BondType::Single
        && (a2_num == 8 || (a2_num == 7 && total_h_count(aid2, mol) == 1))
    {
        return true;
    }
    false
}

// ============================================================================
// Bounds Matrix Builder
// ============================================================================

fn is_conjugated_5ring_bond(mol: &Molecule, bond_idx: usize, rings: &[Vec<usize>]) -> bool {
    let bond = &mol.bonds[bond_idx];
    let a1 = bond.atom1;
    let a2 = bond.atom2;
    let z1 = element_to_atomic_num(&mol.atoms[a1].symbol);
    let z2 = element_to_atomic_num(&mol.atoms[a2].symbol);
    if (z1 <= 10 && z2 <= 10) || bond.bond_type == BondType::Single {
        return false;
    }
    rings
        .iter()
        .any(|r| r.len() == 5 && r.contains(&a1) && r.contains(&a2))
}

fn build_distance_bounds(mol: &Molecule, config: &ETKDGConfig) -> DistanceBounds {
    let n_atoms = mol.atoms.len();
    let mut bounds = DistanceBounds::new(n_atoms);
    let n_bonds = mol.bonds.len();
    let mut accum = ComputedData::new(n_atoms, n_bonds);

    let rings = crate::molecule::graph::find_rings(mol);

    // 1-2 bounds
    for (bond_idx, bond) in mol.bonds.iter().enumerate() {
        let i = bond.atom1;
        let j = bond.atom2;
        let hyb_i = crate::molecule::graph::determine_hybridization(i, mol);
        let hyb_j = crate::molecule::graph::determine_hybridization(j, mol);
        let aro_i = crate::molecule::graph::is_aromatic(i, mol);
        let aro_j = crate::molecule::graph::is_aromatic(j, mol);
        let bl = compute_bond_length(
            &mol.atoms[i].symbol,
            hyb_i,
            aro_i,
            &mol.atoms[j].symbol,
            hyb_j,
            aro_j,
            bond.bond_type,
        );
        accum.bond_lengths[bond_idx] = bl;
        let squish = if is_conjugated_5ring_bond(mol, bond_idx, &rings) {
            EXTRA_SQUISH
        } else {
            0.0
        };
        bounds.check_and_set(i, j, bl - squish - DIST12_DELTA, bl + squish + DIST12_DELTA);
        accum.visited12[i.min(j) * n_atoms + i.max(j)] = true;
    }

    // 1-3 bounds
    let mut visited_centers = vec![0usize; n_atoms];
    let mut angle_taken = vec![0.0f64; n_atoms];
    let mut done_paths = vec![false; n_bonds * n_bonds];

    let mut sorted_rings = rings.clone();
    sorted_rings.sort_by_key(|r| r.len());

    for ring in &sorted_rings {
        let rsize = ring.len();
        if rsize < 3 {
            continue;
        }
        let mut aid1 = ring[rsize - 1];
        for i in 0..rsize {
            let aid2 = ring[i];
            let aid3 = ring[(i + 1) % rsize];
            if let (Some(bid1), Some(bid2)) = (
                find_bond_index(mol, aid1, aid2),
                find_bond_index(mol, aid2, aid3),
            ) {
                let id1 = bid1 * n_bonds + bid2;
                let id2 = bid2 * n_bonds + bid1;
                let pid = aid1.min(aid3) * n_atoms + aid1.max(aid3);
                if !done_paths[id1] && !done_paths[id2] {
                    let hyb = crate::molecule::graph::determine_hybridization(aid2, mol);
                    let angle = set_ring_angle(hyb, rsize);
                    let dl =
                        compute_13_dist(accum.bond_lengths[bid1], accum.bond_lengths[bid2], angle);
                    let dist_tol = if is_larger_sp2_atom(aid1, mol, &rings)
                        || is_larger_sp2_atom(aid2, mol, &rings)
                        || is_larger_sp2_atom(aid3, mol, &rings)
                    {
                        DIST13_TOL * 2.0
                    } else {
                        DIST13_TOL
                    };
                    if !accum.visited12[pid] {
                        bounds.check_and_set(aid1, aid3, dl - dist_tol, dl + dist_tol);
                        accum.visited13[pid] = true;
                    }
                    accum.bond_angles[bid1][bid2] = angle;
                    accum.bond_angles[bid2][bid1] = angle;
                    accum.bond_adj[bid1][bid2] = aid2 as i32;
                    accum.bond_adj[bid2][bid1] = aid2 as i32;
                    visited_centers[aid2] += 1;
                    angle_taken[aid2] += angle;
                    done_paths[id1] = true;
                    done_paths[id2] = true;
                }
            }
            aid1 = aid2;
        }
    }

    for aid2 in 0..n_atoms {
        let deg = mol.adjacency[aid2].len();
        let n13 = deg * (deg - 1) / 2;
        if n13 == visited_centers[aid2] {
            continue;
        }
        let hyb = crate::molecule::graph::determine_hybridization(aid2, mol);
        let neighbors: Vec<usize> = mol.adjacency[aid2].to_vec();
        for i1 in 0..neighbors.len() {
            for i2 in 0..i1 {
                let aid1 = neighbors[i2];
                let aid3 = neighbors[i1];
                if let (Some(bid1), Some(bid2)) = (
                    find_bond_index(mol, aid1, aid2),
                    find_bond_index(mol, aid2, aid3),
                ) {
                    if accum.bond_angles[bid1][bid2] >= 0.0 {
                        continue;
                    }
                    let angle = {
                        // Helper to detect double bonds at central atom
                        let is_double = |a: usize| -> bool {
                            mol.bonds.iter().any(|b| {
                                (b.atom1 == aid2 && b.atom2 == a || b.atom2 == aid2 && b.atom1 == a)
                                    && (b.bond_type == BondType::Double
                                        || b.bond_type == BondType::Aromatic)
                            })
                        };

                        if visited_centers[aid2] >= 1 {
                            if matches!(hyb, Hybridization::Sp2) {
                                // For sp2 with mixed single/double bonds, distribute
                                // remaining angle intelligently
                                let remaining = 2.0 * std::f64::consts::PI - angle_taken[aid2];
                                let remaining_pairs = n13 - visited_centers[aid2];

                                if is_double(aid1) || is_double(aid3) {
                                    // This pair involves the double bond
                                    // Give it a larger share if the other remaining pair
                                    // is the single-single pair
                                    let other_is_single_single = if is_double(aid1) {
                                        !is_double(aid3) && remaining_pairs > 1
                                    } else {
                                        !is_double(aid1) && remaining_pairs > 1
                                    };
                                    if other_is_single_single {
                                        // Double-single angles get ~123° each
                                        123.0_f64.to_radians()
                                    } else {
                                        remaining / remaining_pairs as f64
                                    }
                                } else {
                                    // Single-single pair gets the smaller angle
                                    // (~114°, or whatever is left)
                                    if remaining_pairs == 1 {
                                        remaining
                                    } else {
                                        114.0_f64.to_radians()
                                    }
                                }
                            } else if matches!(hyb, Hybridization::Sp3) {
                                let mut a = 109.5_f64.to_radians();
                                for ring in &rings {
                                    if ring.contains(&aid2) && ring.len() == 3 {
                                        a = 116.0_f64.to_radians();
                                    } else if ring.contains(&aid2) && ring.len() == 4 {
                                        a = 112.0_f64.to_radians();
                                    }
                                }
                                a
                            } else if matches!(hyb, Hybridization::Sp3D) {
                                105.0_f64.to_radians()
                            } else if matches!(hyb, Hybridization::Sp3D2) {
                                90.0_f64.to_radians()
                            } else {
                                120.0_f64.to_radians()
                            }
                        } else {
                            match hyb {
                                Hybridization::Sp1 => std::f64::consts::PI,
                                Hybridization::Sp2 => {
                                    let nbrs = &mol.adjacency[aid2];
                                    let double_nbrs: Vec<usize> =
                                        nbrs.iter().filter(|&&n| is_double(n)).copied().collect();
                                    let single_nbrs: Vec<usize> =
                                        nbrs.iter().filter(|&&n| !is_double(n)).copied().collect();

                                    if double_nbrs.len() == 1 && single_nbrs.len() == 2 {
                                        let is_single_pair = !is_double(aid1) && !is_double(aid3);
                                        if is_single_pair {
                                            114.0_f64.to_radians()
                                        } else {
                                            123.0_f64.to_radians()
                                        }
                                    } else {
                                        120.0_f64.to_radians()
                                    }
                                }
                                Hybridization::Sp3 => 109.5_f64.to_radians(),
                                Hybridization::Sp3D => 105.0_f64.to_radians(),
                                Hybridization::Sp3D2 => 90.0_f64.to_radians(),
                            }
                        }
                    };
                    let pid = aid1.min(aid3) * n_atoms + aid1.max(aid3);
                    if !accum.visited12[pid] {
                        let dl = compute_13_dist(
                            accum.bond_lengths[bid1],
                            accum.bond_lengths[bid2],
                            angle,
                        );
                        let dist_tol = if is_larger_sp2_atom(aid1, mol, &rings)
                            || is_larger_sp2_atom(aid2, mol, &rings)
                            || is_larger_sp2_atom(aid3, mol, &rings)
                        {
                            DIST13_TOL * 2.0
                        } else {
                            DIST13_TOL
                        };
                        bounds.check_and_set(aid1, aid3, dl - dist_tol, dl + dist_tol);
                        accum.visited13[pid] = true;
                    }
                    accum.bond_angles[bid1][bid2] = angle;
                    accum.bond_angles[bid2][bid1] = angle;
                    accum.bond_adj[bid1][bid2] = aid2 as i32;
                    accum.bond_adj[bid2][bid1] = aid2 as i32;
                    angle_taken[aid2] += angle;
                    visited_centers[aid2] += 1;
                }
            }
        }
    }

    // Topological distance matrix
    let mut dist_mat = vec![vec![0usize; n_atoms]; n_atoms];
    for i in 0..n_atoms {
        for j in 0..n_atoms {
            dist_mat[i][j] = if i == j { 0 } else { usize::MAX };
        }
    }
    for bond in &mol.bonds {
        dist_mat[bond.atom1][bond.atom2] = 1;
        dist_mat[bond.atom2][bond.atom1] = 1;
    }
    for k in 0..n_atoms {
        for i in 0..n_atoms {
            for j in 0..n_atoms {
                if dist_mat[i][k] != usize::MAX && dist_mat[k][j] != usize::MAX {
                    let d = dist_mat[i][k] + dist_mat[k][j];
                    if d < dist_mat[i][j] {
                        dist_mat[i][j] = d;
                    }
                }
            }
        }
    }

    // 1-4 bounds
    let mut done_14_paths = HashSet::new();
    for ring in &sorted_rings {
        let rsize = ring.len();
        if rsize < 4 {
            continue;
        }
        let mut bid1 = find_bond_index(mol, ring[rsize - 1], ring[0]).unwrap();
        for i in 0..rsize {
            let bid2 = find_bond_index(mol, ring[i], ring[(i + 1) % rsize]).unwrap();
            let bid3 = find_bond_index(mol, ring[(i + 1) % rsize], ring[(i + 2) % rsize]).unwrap();
            let id1 = (bid1, bid2, bid3);
            let id2 = (bid3, bid2, bid1);
            done_14_paths.insert(id1);
            done_14_paths.insert(id2);
            let atm2 = accum.bond_adj[bid1][bid2];
            let atm3 = accum.bond_adj[bid2][bid3];
            if atm2 < 0 || atm3 < 0 {
                bid1 = bid2;
                continue;
            }
            let aid1 = mol.bonds[bid1].atom1 + mol.bonds[bid1].atom2 - atm2 as usize;
            let aid4 = mol.bonds[bid3].atom1 + mol.bonds[bid3].atom2 - atm3 as usize;
            let pid = aid1.min(aid4) * n_atoms + aid1.max(aid4);
            if accum.visited_bound(pid, 2) {
                bid1 = bid2;
                continue;
            }
            if dist_mat[aid1][aid4] < 3 {
                bid1 = bid2;
                continue;
            }
            let bl1 = accum.bond_lengths[bid1];
            let bl2 = accum.bond_lengths[bid2];
            let bl3 = accum.bond_lengths[bid3];
            let ba12 = accum.bond_angles[bid1][bid2];
            let ba23 = accum.bond_angles[bid2][bid3];
            if ba12 < 0.0 || ba23 < 0.0 {
                bid1 = bid2;
                continue;
            }
            let hyb2 = crate::molecule::graph::determine_hybridization(atm2 as usize, mol);
            let hyb3 = crate::molecule::graph::determine_hybridization(atm3 as usize, mol);

            let (prefer_cis, prefer_trans) = if config.use_macrocycle_14config && rsize >= 9 {
                // Macrocycle: check for amide bond with +0.1 offset
                if is_amide_bond(mol, aid1, aid4)
                    || is_ester_bond(mol, aid1, aid4)
                    || is_amide_bond(mol, mol.bonds[bid1].atom1, mol.bonds[bid1].atom2)
                    || is_amide_bond(mol, mol.bonds[bid3].atom1, mol.bonds[bid3].atom2)
                {
                    (false, true) // macrocycle amide -> trans with offset
                } else {
                    (false, false) // macrocycles don't assume cis
                }
            } else {
                let is_cis = rsize <= 8
                    && matches!(hyb2, Hybridization::Sp2)
                    && matches!(hyb3, Hybridization::Sp2);
                (is_cis, false)
            };

            let (dl, du, ptype) = if prefer_cis {
                let d = compute_14_dist_cis(bl1, bl2, bl3, ba12, ba23);
                (d - GEN_DIST_TOL, d + GEN_DIST_TOL, Path14Type::Cis)
            } else if prefer_trans {
                let mut d = compute_14_dist_trans(bl1, bl2, bl3, ba12, ba23);
                // Macrocycle amide offset: +0.1 to account for ring strain
                if config.use_macrocycle_14config && rsize >= 9 {
                    d += 0.1;
                }
                (d - GEN_DIST_TOL, d + GEN_DIST_TOL, Path14Type::Trans)
            } else {
                let dc = compute_14_dist_cis(bl1, bl2, bl3, ba12, ba23);
                let dt = compute_14_dist_trans(bl1, bl2, bl3, ba12, ba23);
                let dl = dc.min(dt);
                let du = dc.max(dt);
                if (du - dl).abs() < DIST12_DELTA {
                    (dl - GEN_DIST_TOL, du + GEN_DIST_TOL, Path14Type::Other)
                } else {
                    (dl, du, Path14Type::Other)
                }
            };
            bounds.check_and_set(aid1, aid4, dl, du);
            accum.visited14[pid] = true;
            accum.paths14.push(Path14Config {
                bid1,
                bid2,
                bid3,
                ptype,
            });
            let path_id = bid1 as u64 * n_bonds as u64 * n_bonds as u64
                + bid2 as u64 * n_bonds as u64
                + bid3 as u64;
            match ptype {
                Path14Type::Cis => {
                    accum.cis_paths.insert(path_id);
                }
                Path14Type::Trans => {
                    accum.trans_paths.insert(path_id);
                }
                _ => {}
            }
            bid1 = bid2;
        }
    }

    // Helper: check if a bond index is in any ring
    let bond_in_ring: Vec<bool> = (0..n_bonds)
        .map(|bid| {
            sorted_rings.iter().any(|ring| {
                (0..ring.len())
                    .any(|i| find_bond_index(mol, ring[i], ring[(i + 1) % ring.len()]) == Some(bid))
            })
        })
        .collect();

    // Helper: check if two bonds share any ring
    let bonds_share_a_ring = |b1: usize, b2: usize| -> bool {
        sorted_rings.iter().any(|ring| {
            let ring_bonds: Vec<usize> = (0..ring.len())
                .filter_map(|i| find_bond_index(mol, ring[i], ring[(i + 1) % ring.len()]))
                .collect();
            ring_bonds.contains(&b1) && ring_bonds.contains(&b2)
        })
    };

    for (bid2, bond2) in mol.bonds.iter().enumerate() {
        let aid2 = bond2.atom1;
        let aid3 = bond2.atom2;
        for &nbr2 in &mol.adjacency[aid2] {
            if nbr2 == aid3 {
                continue;
            }
            let Some(bid1) = find_bond_index(mol, nbr2, aid2) else {
                continue;
            };
            for &nbr3 in &mol.adjacency[aid3] {
                if nbr3 == aid2 || nbr3 == nbr2 {
                    continue;
                }
                let Some(bid3) = find_bond_index(mol, aid3, nbr3) else {
                    continue;
                };
                let id1 = (bid1, bid2, bid3);
                let id2 = (bid3, bid2, bid1);
                if done_14_paths.contains(&id1) || done_14_paths.contains(&id2) {
                    continue;
                }
                done_14_paths.insert(id1);
                done_14_paths.insert(id2);
                let aid1 = nbr2;
                let aid4 = nbr3;
                let pid = aid1.min(aid4) * n_atoms + aid1.max(aid4);
                if accum.visited_bound(pid, 2) {
                    continue;
                }
                let bl1 = accum.bond_lengths[bid1];
                let bl2 = accum.bond_lengths[bid2];
                let bl3 = accum.bond_lengths[bid3];
                let ba12 = accum.bond_angles[bid1][bid2];
                let ba23 = accum.bond_angles[bid2][bid3];
                if ba12 < 0.0 || ba23 < 0.0 {
                    continue;
                }

                let b1_b2_share = bonds_share_a_ring(bid1, bid2);
                let b2_b3_share = bonds_share_a_ring(bid2, bid3);
                let b1_in_ring = bond_in_ring[bid1];
                let b2_in_ring = bond_in_ring[bid2];
                let b3_in_ring = bond_in_ring[bid3];

                let (dl, du, ptype) = if b1_b2_share || b2_b3_share {
                    // Two bonds in the same ring (but not all three in same ring)
                    let hyb2 = crate::molecule::graph::determine_hybridization(aid2, mol);
                    let hyb3 = crate::molecule::graph::determine_hybridization(aid3, mol);
                    if matches!(hyb2, Hybridization::Sp2) && matches!(hyb3, Hybridization::Sp2) {
                        let dt = compute_14_dist_trans(bl1, bl2, bl3, ba12, ba23);
                        (dt - GEN_DIST_TOL, dt + GEN_DIST_TOL, Path14Type::Trans)
                    } else {
                        let dc = compute_14_dist_cis(bl1, bl2, bl3, ba12, ba23);
                        let dt = compute_14_dist_trans(bl1, bl2, bl3, ba12, ba23);
                        let dl = dc.min(dt);
                        let du = dc.max(dt);
                        if (du - dl).abs() < DIST12_DELTA {
                            (dl - GEN_DIST_TOL, du + GEN_DIST_TOL, Path14Type::Other)
                        } else {
                            (dl, du, Path14Type::Other)
                        }
                    }
                } else if b2_in_ring && (b1_in_ring || b3_in_ring) {
                    // Two bonds in different rings
                    let dc = compute_14_dist_cis(bl1, bl2, bl3, ba12, ba23);
                    let dt = compute_14_dist_trans(bl1, bl2, bl3, ba12, ba23);
                    let dl = dc.min(dt);
                    let du = dc.max(dt);
                    if (du - dl).abs() < DIST12_DELTA {
                        (dl - GEN_DIST_TOL, du + GEN_DIST_TOL, Path14Type::Other)
                    } else {
                        (dl, du, Path14Type::Other)
                    }
                } else {
                    // Chain path (no ring bonds or only middle bond in ring)
                    match bond2.bond_type {
                        BondType::Double => {
                            if mol.bonds[bid1].bond_type == BondType::Double
                                || mol.bonds[bid3].bond_type == BondType::Double
                            {
                                let d = compute_14_dist_cis(bl1, bl2, bl3, ba12, ba23);
                                (d - GEN_DIST_TOL, d + GEN_DIST_TOL, Path14Type::Cis)
                            } else {
                                let dc = compute_14_dist_cis(bl1, bl2, bl3, ba12, ba23);
                                let dt = compute_14_dist_trans(bl1, bl2, bl3, ba12, ba23);
                                let dl = dc.min(dt);
                                let du = dc.max(dt);
                                if (du - dl).abs() < DIST12_DELTA {
                                    (dl - GEN_DIST_TOL, du + GEN_DIST_TOL, Path14Type::Other)
                                } else {
                                    (dl, du, Path14Type::Other)
                                }
                            }
                        }
                        BondType::Single => {
                            let a2_num = atom_num(aid2, mol);
                            let a3_num = atom_num(aid3, mol);
                            if a2_num == 16
                                && a3_num == 16
                                && atom_degree(aid2, mol) == 2
                                && atom_degree(aid3, mol) == 2
                            {
                                let d = compute_14_dist_at_dihedral(
                                    bl1,
                                    bl2,
                                    bl3,
                                    ba12,
                                    ba23,
                                    std::f64::consts::FRAC_PI_2,
                                );
                                (d - GEN_DIST_TOL, d + GEN_DIST_TOL, Path14Type::Other)
                            } else if check_amide_ester_14(mol, aid1, aid2, aid3, aid4, bid1, bid3)
                                || check_amide_ester_14(mol, aid4, aid3, aid2, aid1, bid3, bid1)
                            {
                                if config.force_trans_amides {
                                    let a1_num = atom_num(aid1, mol);
                                    let a4_num = atom_num(aid4, mol);
                                    if (a1_num == 1
                                        && a2_num == 7
                                        && atom_degree(aid2, mol) == 3
                                        && total_h_count(aid2, mol) == 1)
                                        || (a4_num == 1
                                            && a3_num == 7
                                            && atom_degree(aid3, mol) == 3
                                            && total_h_count(aid3, mol) == 1)
                                    {
                                        let dt = compute_14_dist_trans(bl1, bl2, bl3, ba12, ba23);
                                        (dt - GEN_DIST_TOL, dt + GEN_DIST_TOL, Path14Type::Trans)
                                    } else {
                                        let dc = compute_14_dist_cis(bl1, bl2, bl3, ba12, ba23);
                                        (dc - GEN_DIST_TOL, dc + GEN_DIST_TOL, Path14Type::Cis)
                                    }
                                } else {
                                    let dc = compute_14_dist_cis(bl1, bl2, bl3, ba12, ba23);
                                    let dt = compute_14_dist_trans(bl1, bl2, bl3, ba12, ba23);
                                    let dl = dc.min(dt);
                                    let du = dc.max(dt);
                                    if (du - dl).abs() < DIST12_DELTA {
                                        (dl - GEN_DIST_TOL, du + GEN_DIST_TOL, Path14Type::Other)
                                    } else {
                                        (dl, du, Path14Type::Other)
                                    }
                                }
                            } else {
                                let dc = compute_14_dist_cis(bl1, bl2, bl3, ba12, ba23);
                                let dt = compute_14_dist_trans(bl1, bl2, bl3, ba12, ba23);
                                let mut dl = dc.min(dt);
                                let mut du = dc.max(dt);
                                if (du - dl).abs() < DIST12_DELTA {
                                    dl -= GEN_DIST_TOL;
                                    du += GEN_DIST_TOL;
                                }
                                (dl, du, Path14Type::Other)
                            }
                        }
                        _ => {
                            let dc = compute_14_dist_cis(bl1, bl2, bl3, ba12, ba23);
                            let dt = compute_14_dist_trans(bl1, bl2, bl3, ba12, ba23);
                            let mut dl = dc.min(dt);
                            let mut du = dc.max(dt);
                            if (du - dl).abs() < DIST12_DELTA {
                                dl -= GEN_DIST_TOL;
                                du += GEN_DIST_TOL;
                            }
                            (dl, du, Path14Type::Other)
                        }
                    }
                };
                bounds.check_and_set(aid1, aid4, dl, du);
                accum.visited14[pid] = true;
                accum.paths14.push(Path14Config {
                    bid1,
                    bid2,
                    bid3,
                    ptype,
                });
            }
        }
    }

    // 1-5 bounds
    let paths14_copy: Vec<Path14Config> = accum.paths14.clone();
    for path14 in paths14_copy {
        set_15_bounds_helper(
            mol,
            &mut accum,
            &mut bounds,
            &path14,
            &dist_mat,
            n_atoms,
            n_bonds,
            false,
        );
        let rev = Path14Config {
            bid1: path14.bid3,
            bid2: path14.bid2,
            bid3: path14.bid1,
            ptype: path14.ptype,
        };
        set_15_bounds_helper(
            mol,
            &mut accum,
            &mut bounds,
            &rev,
            &dist_mat,
            n_atoms,
            n_bonds,
            true,
        );
    }

    // VDW lower bounds with H-bond detection
    let mut h_bond_donor_h = vec![false; n_atoms];
    let mut h_bond_acceptor = vec![false; n_atoms];
    for i in 0..n_atoms {
        let sym = &mol.atoms[i].symbol;
        if *sym == "H" {
            for &nbr in &mol.adjacency[i] {
                let nbr_sym = &mol.atoms[nbr].symbol;
                if *nbr_sym == "N" || *nbr_sym == "O" {
                    h_bond_donor_h[i] = true;
                    break;
                }
            }
        }
        if *sym == "N" || *sym == "O" {
            h_bond_acceptor[i] = true;
        }
    }

    for i in 1..n_atoms {
        for j in 0..i {
            if bounds.lower[j][i] < DIST12_DELTA {
                let vw1 = vdw_radius(&mol.atoms[i].symbol);
                let vw2 = vdw_radius(&mol.atoms[j].symbol);
                let d = dist_mat[i][j];
                let lb = if (h_bond_donor_h[i] && h_bond_acceptor[j])
                    || (h_bond_acceptor[i] && h_bond_donor_h[j])
                {
                    H_BOND_LENGTH
                } else if d == 4 {
                    VDW_SCALE_15 * (vw1 + vw2)
                } else if d == 5 {
                    (VDW_SCALE_15 + 0.5 * (1.0 - VDW_SCALE_15)) * (vw1 + vw2)
                } else {
                    vw1 + vw2
                };
                bounds.set_lower(i, j, lb);
            }
        }
    }

    bounds
}

fn set_15_bounds_helper(
    mol: &Molecule,
    accum: &mut ComputedData,
    bounds: &mut DistanceBounds,
    path14: &Path14Config,
    dist_mat: &[Vec<usize>],
    n_atoms: usize,
    n_bonds: usize,
    reversed: bool,
) {
    let bid1 = path14.bid1;
    let bid2 = path14.bid2;
    let bid3 = path14.bid3;
    let aid2 = accum.bond_adj[bid1][bid2] as usize;
    let aid3 = accum.bond_adj[bid2][bid3] as usize;
    let aid1 = if mol.bonds[bid1].atom1 == aid2 {
        mol.bonds[bid1].atom2
    } else {
        mol.bonds[bid1].atom1
    };
    let aid4 = if mol.bonds[bid3].atom1 == aid3 {
        mol.bonds[bid3].atom2
    } else {
        mol.bonds[bid3].atom1
    };

    let d1 = accum.bond_lengths[bid1];
    let d2 = accum.bond_lengths[bid2];
    let d3 = accum.bond_lengths[bid3];
    let ang12 = accum.bond_angles[bid1][bid2];
    let ang23 = accum.bond_angles[bid2][bid3];
    if ang12 < 0.0 || ang23 < 0.0 {
        return;
    }

    for i in 0..mol.bonds.len() {
        if accum.bond_adj[bid3][i] == aid4 as i32 {
            let aid5 = if mol.bonds[i].atom1 == aid4 {
                mol.bonds[i].atom2
            } else {
                mol.bonds[i].atom1
            };
            if aid1 == aid5 {
                continue;
            }
            let pid = aid1.min(aid5) * n_atoms + aid1.max(aid5);
            if accum.visited_bound(pid, 3) {
                continue;
            }
            if dist_mat[aid1.max(aid5)][aid1.min(aid5)] < 4 {
                continue;
            }
            if accum.set15[pid] {
                continue;
            }

            let d4 = accum.bond_lengths[i];
            let ang34 = accum.bond_angles[bid3][i];
            if ang34 < 0.0 {
                continue;
            }

            let (dl, du) = if reversed {
                match path14.ptype {
                    Path14Type::Cis => {
                        let dl = compute_15_dist_cis_cis(d4, d3, d2, d1, ang34, ang23, ang12)
                            - DIST15_TOL;
                        let du = compute_15_dist_cis_trans(d4, d3, d2, d1, ang34, ang23, ang12)
                            + DIST15_TOL;
                        (dl, du)
                    }
                    Path14Type::Trans => {
                        let dl = compute_15_dist_trans_cis(d4, d3, d2, d1, ang34, ang23, ang12)
                            - DIST15_TOL;
                        let du = compute_15_dist_trans_trans(d4, d3, d2, d1, ang34, ang23, ang12)
                            + DIST15_TOL;
                        (dl, du)
                    }
                    Path14Type::Other => {
                        let path_id = bid2 as u64 * n_bonds as u64 * n_bonds as u64
                            + bid3 as u64 * n_bonds as u64
                            + i as u64;
                        if accum.cis_paths.contains(&path_id) {
                            let dl = compute_15_dist_cis_cis(d4, d3, d2, d1, ang34, ang23, ang12)
                                - DIST15_TOL;
                            let du = compute_15_dist_cis_trans(d4, d3, d2, d1, ang34, ang23, ang12)
                                + DIST15_TOL;
                            (dl, du)
                        } else if accum.trans_paths.contains(&path_id) {
                            let dl = compute_15_dist_trans_cis(d4, d3, d2, d1, ang34, ang23, ang12)
                                - DIST15_TOL;
                            let du =
                                compute_15_dist_trans_trans(d4, d3, d2, d1, ang34, ang23, ang12)
                                    + DIST15_TOL;
                            (dl, du)
                        } else {
                            let vw1 = vdw_radius(&mol.atoms[aid1].symbol);
                            let vw5 = vdw_radius(&mol.atoms[aid5].symbol);
                            (VDW_SCALE_15 * (vw1 + vw5), MAX_UPPER)
                        }
                    }
                }
            } else {
                match path14.ptype {
                    Path14Type::Cis => {
                        let dl = compute_15_dist_cis_cis(d1, d2, d3, d4, ang12, ang23, ang34)
                            - DIST15_TOL;
                        let du = compute_15_dist_cis_trans(d1, d2, d3, d4, ang12, ang23, ang34)
                            + DIST15_TOL;
                        (dl, du)
                    }
                    Path14Type::Trans => {
                        let dl = compute_15_dist_trans_cis(d1, d2, d3, d4, ang12, ang23, ang34)
                            - DIST15_TOL;
                        let du = compute_15_dist_trans_trans(d1, d2, d3, d4, ang12, ang23, ang34)
                            + DIST15_TOL;
                        (dl, du)
                    }
                    Path14Type::Other => {
                        let path_id = bid2 as u64 * n_bonds as u64 * n_bonds as u64
                            + bid3 as u64 * n_bonds as u64
                            + i as u64;
                        if accum.cis_paths.contains(&path_id) {
                            let dl = compute_15_dist_cis_cis(d1, d2, d3, d4, ang12, ang23, ang34)
                                - DIST15_TOL;
                            let du = compute_15_dist_cis_trans(d1, d2, d3, d4, ang12, ang23, ang34)
                                + DIST15_TOL;
                            (dl, du)
                        } else if accum.trans_paths.contains(&path_id) {
                            let dl = compute_15_dist_trans_cis(d1, d2, d3, d4, ang12, ang23, ang34)
                                - DIST15_TOL;
                            let du =
                                compute_15_dist_trans_trans(d1, d2, d3, d4, ang12, ang23, ang34)
                                    + DIST15_TOL;
                            (dl, du)
                        } else {
                            let vw1 = vdw_radius(&mol.atoms[aid1].symbol);
                            let vw5 = vdw_radius(&mol.atoms[aid5].symbol);
                            (VDW_SCALE_15 * (vw1 + vw5), MAX_UPPER)
                        }
                    }
                }
            };
            bounds.check_and_set(aid1, aid5, dl, du);
            accum.set15[pid] = true;
        }
    }
}

// ============================================================================
// Metric Matrix Embedding
// ============================================================================
// Metric Matrix Embedding
// ============================================================================

fn generate_initial_coords_from_bounds(bounds: &DistanceBounds, rng: &mut Rng) -> Vec<[f64; 4]> {
    let n = bounds.n_atoms;
    if n < 2 {
        return vec![[0.0; 4]; n];
    }
    let mut dist = vec![vec![0.0f64; n]; n];
    for i in 0..n {
        for j in (i + 1)..n {
            let lo = bounds.lower[i][j];
            let hi = bounds.upper[i][j];
            let t = lo + (hi - lo) * rng.random_f64();
            dist[i][j] = t;
            dist[j][i] = t;
        }
    }
    let n_f = n as f64;
    let mut d2 = vec![vec![0.0f64; n]; n];
    let mut row_sums = vec![0.0f64; n];
    let mut total_sum = 0.0f64;
    for i in 0..n {
        for j in 0..n {
            d2[i][j] = dist[i][j] * dist[i][j];
            row_sums[i] += d2[i][j];
            total_sum += d2[i][j];
        }
    }
    let row_means: Vec<f64> = row_sums.iter().map(|s| s / n_f).collect();
    let grand_mean = total_sum / (n_f * n_f);
    let mut b = vec![vec![0.0f64; n]; n];
    for i in 0..n {
        for j in 0..n {
            b[i][j] = -0.5 * (d2[i][j] - row_means[i] - row_means[j] + grand_mean);
        }
    }
    let (eigenvalues, eigenvectors) = jacobi_eigen(&b, 100);
    let mut indices: Vec<usize> = (0..n).collect();
    indices.sort_by(|&a, &b| eigenvalues[b].partial_cmp(&eigenvalues[a]).unwrap());
    let dim = 4.min(n);
    let mut coords_4d = vec![[0.0f64; 4]; n];
    for i in 0..n {
        for d in 0..dim {
            if d < indices.len() {
                let idx = indices[d];
                let ev = eigenvalues[idx];
                coords_4d[i][d] = if ev > 0.0 {
                    eigenvectors[idx][i] * ev.sqrt()
                } else {
                    1.0 - 2.0 * rng.random_f64()
                };
            }
        }
    }
    coords_4d
}

fn jacobi_eigen(matrix: &[Vec<f64>], max_sweeps: usize) -> (Vec<f64>, Vec<Vec<f64>>) {
    let n = matrix.len();
    let mut a = matrix.to_vec();
    let mut v = vec![vec![0.0f64; n]; n];
    for i in 0..n {
        v[i][i] = 1.0;
    }
    let mut d = vec![0.0f64; n];
    for i in 0..n {
        d[i] = a[i][i];
    }
    let mut b = d.clone();
    let mut z = vec![0.0f64; n];
    for _sweep in 0..max_sweeps {
        let mut sum = 0.0;
        for i in 0..n - 1 {
            for j in i + 1..n {
                sum += a[i][j].abs();
            }
        }
        if sum < 1e-12 {
            break;
        }
        let threshold = if _sweep < 3 {
            0.2 * sum / (n * n) as f64
        } else {
            0.0
        };
        for p in 0..n - 1 {
            for q in p + 1..n {
                let apq = a[p][q].abs();
                let g = 100.0 * apq;
                if _sweep > 3 && d[p].abs() + g == d[p].abs() && d[q].abs() + g == d[q].abs() {
                    a[p][q] = 0.0;
                    continue;
                }
                if apq <= threshold {
                    continue;
                }
                let h = d[q] - d[p];
                let t = if h.abs() + g == h.abs() {
                    a[p][q] / h
                } else {
                    let theta = 0.5 * h / a[p][q];
                    let mut tt = 1.0 / (theta.abs() + (1.0 + theta * theta).sqrt());
                    if theta < 0.0 {
                        tt = -tt;
                    }
                    tt
                };
                let c = 1.0 / (1.0 + t * t).sqrt();
                let s = t * c;
                let tau = s / (1.0 + c);
                let h = t * a[p][q];
                z[p] -= h;
                z[q] += h;
                d[p] -= h;
                d[q] += h;
                a[p][q] = 0.0;
                for r in 0..p {
                    let g = a[r][p];
                    let h = a[r][q];
                    a[r][p] = g - s * (h + g * tau);
                    a[r][q] = h + s * (g - h * tau);
                }
                for r in p + 1..q {
                    let g = a[p][r];
                    let h = a[r][q];
                    a[p][r] = g - s * (h + g * tau);
                    a[r][q] = h + s * (g - h * tau);
                }
                for r in q + 1..n {
                    let g = a[p][r];
                    let h = a[q][r];
                    a[p][r] = g - s * (h + g * tau);
                    a[q][r] = h + s * (g - h * tau);
                }
                for r in 0..n {
                    let g = v[r][p];
                    let h = v[r][q];
                    v[r][p] = g - s * (h + g * tau);
                    v[r][q] = h + s * (g - h * tau);
                }
            }
        }
        for i in 0..n {
            b[i] += z[i];
            d[i] = b[i];
            z[i] = 0.0;
        }
    }
    let mut evecs = vec![vec![0.0f64; n]; n];
    for i in 0..n {
        for r in 0..n {
            evecs[i][r] = v[r][i];
        }
    }
    (d, evecs)
}

fn minimize_4d_first(
    coords_4d: &mut [[f64; 4]],
    bounds: &DistanceBounds,
    chiral_centers: &[ChiralCenter],
    max_iter: usize,
) -> f64 {
    let n = coords_4d.len();
    let mut best_energy = f64::INFINITY;
    let mut best_coords = coords_4d.to_vec();
    for _ in 0..max_iter {
        let mut grad = vec![[0.0f64; 4]; n];
        let mut energy = 0.0;
        for i in 0..n {
            for j in (i + 1)..n {
                let lo = bounds.lower[i][j];
                let hi = bounds.upper[i][j];
                if hi >= MAX_UPPER || (hi - lo) > BASIN_THRESH {
                    continue;
                }
                let dx = coords_4d[i][0] - coords_4d[j][0];
                let dy = coords_4d[i][1] - coords_4d[j][1];
                let dz = coords_4d[i][2] - coords_4d[j][2];
                let dw = coords_4d[i][3] - coords_4d[j][3];
                let d4 = (dx * dx + dy * dy + dz * dz + dw * dw).sqrt().max(1e-10);
                let lo_viol = (lo - d4).max(0.0);
                let hi_viol = (d4 - hi).max(0.0);
                energy += lo_viol * lo_viol + hi_viol * hi_viol;
                let f = -2.0 * (lo_viol - hi_viol) / d4;
                for dim in 0..4 {
                    let diff = coords_4d[i][dim] - coords_4d[j][dim];
                    grad[i][dim] += f * diff;
                    grad[j][dim] -= f * diff;
                }
            }
        }
        energy += chiral_4d_penalty(coords_4d, chiral_centers, FIRST_MIN_WEIGHT_CHIRAL);
        chiral_4d_gradient(
            coords_4d,
            chiral_centers,
            FIRST_MIN_WEIGHT_CHIRAL,
            &mut grad,
        );
        for i in 0..n {
            energy += FIRST_MIN_WEIGHT_FOURTH * coords_4d[i][3] * coords_4d[i][3];
            grad[i][3] += 2.0 * FIRST_MIN_WEIGHT_FOURTH * coords_4d[i][3];
        }
        if energy < best_energy {
            best_energy = energy;
            best_coords = coords_4d.to_vec();
        }
        let max_g = grad
            .iter()
            .map(|g| (g[0] * g[0] + g[1] * g[1] + g[2] * g[2] + g[3] * g[3]).sqrt())
            .fold(0.0f64, f64::max);
        if max_g < 1e-3 {
            break;
        }
        let step = 0.1 / max_g.max(1e-10);
        for i in 0..n {
            for dim in 0..4 {
                coords_4d[i][dim] -= step * grad[i][dim];
            }
        }
    }
    coords_4d.copy_from_slice(&best_coords);
    best_energy
}

fn minimize_4d_collapse(
    coords_4d: &mut [[f64; 4]],
    bounds: &DistanceBounds,
    chiral_centers: &[ChiralCenter],
    max_iter: usize,
) {
    let n = coords_4d.len();
    for _ in 0..max_iter {
        let mut grad = vec![[0.0f64; 4]; n];
        for i in 0..n {
            for j in (i + 1)..n {
                let lo = bounds.lower[i][j];
                let hi = bounds.upper[i][j];
                if hi >= MAX_UPPER || (hi - lo) > BASIN_THRESH {
                    continue;
                }
                let dx = coords_4d[i][0] - coords_4d[j][0];
                let dy = coords_4d[i][1] - coords_4d[j][1];
                let dz = coords_4d[i][2] - coords_4d[j][2];
                let dw = coords_4d[i][3] - coords_4d[j][3];
                let d4 = (dx * dx + dy * dy + dz * dz + dw * dw).sqrt().max(1e-10);
                let lo_viol = (lo - d4).max(0.0);
                let hi_viol = (d4 - hi).max(0.0);
                let f = -2.0 * (lo_viol - hi_viol) / d4;
                for dim in 0..4 {
                    let diff = coords_4d[i][dim] - coords_4d[j][dim];
                    grad[i][dim] += f * diff;
                    grad[j][dim] -= f * diff;
                }
            }
        }
        chiral_4d_gradient(
            coords_4d,
            chiral_centers,
            FOURTH_MIN_WEIGHT_CHIRAL,
            &mut grad,
        );
        for i in 0..n {
            grad[i][3] += 2.0 * FOURTH_MIN_WEIGHT_FOURTH * coords_4d[i][3];
        }
        let max_g = grad
            .iter()
            .map(|g| (g[0] * g[0] + g[1] * g[1] + g[2] * g[2] + g[3] * g[3]).sqrt())
            .fold(0.0f64, f64::max);
        if max_g < 1e-3 {
            break;
        }
        let step = 0.1 / max_g.max(1e-10);
        for i in 0..n {
            for dim in 0..4 {
                coords_4d[i][dim] -= step * grad[i][dim];
            }
        }
    }
}

fn chiral_4d_penalty(coords_4d: &[[f64; 4]], chiral_centers: &[ChiralCenter], weight: f64) -> f64 {
    let mut energy = 0.0;
    for cc in chiral_centers {
        let c = cc.central;
        let n1 = cc.neighbors[0];
        let n2 = cc.neighbors[1];
        let n3 = cc.neighbors[2];
        let v1 = [
            coords_4d[n1][0] - coords_4d[c][0],
            coords_4d[n1][1] - coords_4d[c][1],
            coords_4d[n1][2] - coords_4d[c][2],
        ];
        let v2 = [
            coords_4d[n2][0] - coords_4d[c][0],
            coords_4d[n2][1] - coords_4d[c][1],
            coords_4d[n2][2] - coords_4d[c][2],
        ];
        let v3 = [
            coords_4d[n3][0] - coords_4d[c][0],
            coords_4d[n3][1] - coords_4d[c][1],
            coords_4d[n3][2] - coords_4d[c][2],
        ];
        let vol = v1[0] * (v2[1] * v3[2] - v2[2] * v3[1])
            + v1[1] * (v2[2] * v3[0] - v2[0] * v3[2])
            + v1[2] * (v2[0] * v3[1] - v2[1] * v3[0]);
        let penalty = if vol < cc.vol_lower {
            (vol - cc.vol_lower) * (vol - cc.vol_lower)
        } else if vol > cc.vol_upper {
            (vol - cc.vol_upper) * (vol - cc.vol_upper)
        } else {
            0.0
        };
        energy += weight * penalty;
    }
    energy
}

fn chiral_4d_gradient(
    coords_4d: &[[f64; 4]],
    chiral_centers: &[ChiralCenter],
    weight: f64,
    grad: &mut [[f64; 4]],
) {
    for cc in chiral_centers {
        let c = cc.central;
        let n1 = cc.neighbors[0];
        let n2 = cc.neighbors[1];
        let n3 = cc.neighbors[2];
        let v1 = [
            coords_4d[n1][0] - coords_4d[c][0],
            coords_4d[n1][1] - coords_4d[c][1],
            coords_4d[n1][2] - coords_4d[c][2],
        ];
        let v2 = [
            coords_4d[n2][0] - coords_4d[c][0],
            coords_4d[n2][1] - coords_4d[c][1],
            coords_4d[n2][2] - coords_4d[c][2],
        ];
        let v3 = [
            coords_4d[n3][0] - coords_4d[c][0],
            coords_4d[n3][1] - coords_4d[c][1],
            coords_4d[n3][2] - coords_4d[c][2],
        ];
        let vol = v1[0] * (v2[1] * v3[2] - v2[2] * v3[1])
            + v1[1] * (v2[2] * v3[0] - v2[0] * v3[2])
            + v1[2] * (v2[0] * v3[1] - v2[1] * v3[0]);
        let mut dvol = 0.0;
        if vol < cc.vol_lower {
            dvol = 2.0 * weight * (vol - cc.vol_lower);
        } else if vol > cc.vol_upper {
            dvol = 2.0 * weight * (vol - cc.vol_upper);
        }
        if dvol != 0.0 {
            let g_n1 = cross_product(v2, v3);
            let g_n2 = cross_product(v3, v1);
            let g_n3 = cross_product(v1, v2);
            let g_c = [
                -(g_n1[0] + g_n2[0] + g_n3[0]),
                -(g_n1[1] + g_n2[1] + g_n3[1]),
                -(g_n1[2] + g_n2[2] + g_n3[2]),
            ];
            for dim in 0..3 {
                grad[c][dim] += dvol * g_c[dim];
                grad[n1][dim] += dvol * g_n1[dim];
                grad[n2][dim] += dvol * g_n2[dim];
                grad[n3][dim] += dvol * g_n3[dim];
            }
        }
    }
}

// ============================================================================
// Chirality
// ============================================================================

#[derive(Debug, Clone)]
struct ChiralCenter {
    central: usize,
    neighbors: [usize; 4],
    vol_lower: f64,
    vol_upper: f64,
    vol_scale: f64,
}

fn find_chiral_centers(mol: &Molecule) -> (Vec<ChiralCenter>, Vec<ChiralCenter>) {
    let rings = crate::molecule::graph::find_rings(mol);
    let mut chiral = Vec::new();
    let mut tetrahedral = Vec::new();
    for atom_idx in 0..mol.atoms.len() {
        if mol.atoms[atom_idx].symbol == "H" {
            continue;
        }
        let hyb = crate::molecule::graph::determine_hybridization(atom_idx, mol);
        if !matches!(hyb, Hybridization::Sp3) {
            continue;
        }
        let neighbors: Vec<usize> = mol.adjacency[atom_idx]
            .iter()
            .copied()
            .filter(|&n| mol.atoms[n].symbol != "H")
            .collect();
        if neighbors.len() < 3 {
            continue;
        }
        let small_ring_count = rings
            .iter()
            .filter(|r| r.len() < 5 && r.contains(&atom_idx))
            .count();
        let vol_scale = if small_ring_count > 1 { 0.25 } else { 1.0 };
        let mut nbrs = neighbors.clone();
        let (vol_lower, vol_upper) = if nbrs.len() == 4 {
            (5.0, 100.0)
        } else {
            nbrs.push(atom_idx);
            (2.0, 100.0)
        };
        let cc = ChiralCenter {
            central: atom_idx,
            neighbors: [nbrs[0], nbrs[1], nbrs[2], nbrs[3]],
            vol_lower,
            vol_upper,
            vol_scale,
        };
        tetrahedral.push(cc.clone());
        // Basic topological chirality: 4 distinct non-H neighbors by element symbol
        if nbrs.len() == 4 {
            let mut distinct = true;
            for i in 1..4 {
                for j in 0..i {
                    if mol.atoms[nbrs[i]].symbol == mol.atoms[nbrs[j]].symbol {
                        distinct = false;
                        break;
                    }
                }
                if !distinct {
                    break;
                }
            }
            if distinct {
                chiral.push(cc);
            }
        }
    }
    (chiral, tetrahedral)
}

pub(crate) fn chiral_volume(coords: &[[f64; 3]], c: usize, n1: usize, n2: usize, n3: usize) -> f64 {
    let v1 = [
        coords[n1][0] - coords[c][0],
        coords[n1][1] - coords[c][1],
        coords[n1][2] - coords[c][2],
    ];
    let v2 = [
        coords[n2][0] - coords[c][0],
        coords[n2][1] - coords[c][1],
        coords[n2][2] - coords[c][2],
    ];
    let v3 = [
        coords[n3][0] - coords[c][0],
        coords[n3][1] - coords[c][1],
        coords[n3][2] - coords[c][2],
    ];
    let cross = [
        v2[1] * v3[2] - v2[2] * v3[1],
        v2[2] * v3[0] - v2[0] * v3[2],
        v2[0] * v3[1] - v2[1] * v3[0],
    ];
    v1[0] * cross[0] + v1[1] * cross[1] + v1[2] * cross[2]
}

fn normalize_vec(v: [f64; 3]) -> [f64; 3] {
    let n = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    if n < 1e-10 {
        return v;
    }
    [v[0] / n, v[1] / n, v[2] / n]
}

trait Vec3Ops {
    fn dot(self, other: [f64; 3]) -> f64;
}
impl Vec3Ops for [f64; 3] {
    fn dot(self, other: [f64; 3]) -> f64 {
        self[0] * other[0] + self[1] * other[1] + self[2] * other[2]
    }
}

fn cross_product(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn volume_test(coords: &[[f64; 3]], cc: &ChiralCenter) -> bool {
    let p0 = [
        coords[cc.central][0],
        coords[cc.central][1],
        coords[cc.central][2],
    ];
    let p1 = [
        coords[cc.neighbors[0]][0],
        coords[cc.neighbors[0]][1],
        coords[cc.neighbors[0]][2],
    ];
    let p2 = [
        coords[cc.neighbors[1]][0],
        coords[cc.neighbors[1]][1],
        coords[cc.neighbors[1]][2],
    ];
    let p3 = [
        coords[cc.neighbors[2]][0],
        coords[cc.neighbors[2]][1],
        coords[cc.neighbors[2]][2],
    ];
    let p4 = [
        coords[cc.neighbors[3]][0],
        coords[cc.neighbors[3]][1],
        coords[cc.neighbors[3]][2],
    ];
    let v1 = normalize_vec([p0[0] - p1[0], p0[1] - p1[1], p0[2] - p1[2]]);
    let v2 = normalize_vec([p0[0] - p2[0], p0[1] - p2[1], p0[2] - p2[2]]);
    let v3 = normalize_vec([p0[0] - p3[0], p0[1] - p3[1], p0[2] - p3[2]]);
    let v4 = normalize_vec([p0[0] - p4[0], p0[1] - p4[1], p0[2] - p4[2]]);
    let vol_thresh = MIN_TETRAHEDRAL_CHIRAL_VOL * cc.vol_scale;
    let cross = cross_product(v1, v2);
    if cross.dot(v3).abs() < vol_thresh {
        return false;
    }
    let cross = cross_product(v1, v2);
    if cross.dot(v4).abs() < vol_thresh {
        return false;
    }
    let cross = cross_product(v1, v3);
    if cross.dot(v4).abs() < vol_thresh {
        return false;
    }
    let cross = cross_product(v2, v3);
    if cross.dot(v4).abs() < vol_thresh {
        return false;
    }
    true
}

#[allow(dead_code)]
fn check_tetrahedral(coords: &[[f64; 3]], tetrahedral: &[ChiralCenter]) -> bool {
    for tc in tetrahedral {
        if !volume_test(coords, tc) || !center_in_volume(coords, tc) {
            return false;
        }
    }
    true
}

fn check_chiral_centers(coords: &[[f64; 3]], chiral_centers: &[ChiralCenter]) -> bool {
    for cc in chiral_centers {
        let vol = chiral_volume(
            coords,
            cc.central,
            cc.neighbors[0],
            cc.neighbors[1],
            cc.neighbors[2],
        );
        if cc.vol_lower > 0.0
            && vol < cc.vol_lower
            && (vol / cc.vol_lower < 0.8 || have_opposite_sign(vol, cc.vol_lower))
        {
            return false;
        }
        if cc.vol_upper < 0.0
            && vol > cc.vol_upper
            && (vol / cc.vol_upper < 0.8 || have_opposite_sign(vol, cc.vol_upper))
        {
            return false;
        }
    }
    true
}

fn have_opposite_sign(a: f64, b: f64) -> bool {
    (a < 0.0) != (b < 0.0)
}

#[allow(dead_code)]
fn center_in_volume(coords: &[[f64; 3]], cc: &ChiralCenter) -> bool {
    center_in_volume_tol(coords, cc, TETRAHEDRAL_CENTERINVOLUME_TOL)
}

fn center_in_volume_tol(coords: &[[f64; 3]], cc: &ChiralCenter, tol: f64) -> bool {
    if cc.central == cc.neighbors[3] {
        return true;
    }
    same_side_tol(
        coords[cc.neighbors[0]],
        coords[cc.neighbors[1]],
        coords[cc.neighbors[2]],
        coords[cc.neighbors[3]],
        coords[cc.central],
        tol,
    ) && same_side_tol(
        coords[cc.neighbors[1]],
        coords[cc.neighbors[2]],
        coords[cc.neighbors[3]],
        coords[cc.neighbors[0]],
        coords[cc.central],
        tol,
    ) && same_side_tol(
        coords[cc.neighbors[2]],
        coords[cc.neighbors[3]],
        coords[cc.neighbors[0]],
        coords[cc.neighbors[1]],
        coords[cc.central],
        tol,
    ) && same_side_tol(
        coords[cc.neighbors[3]],
        coords[cc.neighbors[0]],
        coords[cc.neighbors[1]],
        coords[cc.neighbors[2]],
        coords[cc.central],
        tol,
    )
}

#[allow(dead_code)]
fn same_side(v1: [f64; 3], v2: [f64; 3], v3: [f64; 3], v4: [f64; 3], p0: [f64; 3]) -> bool {
    same_side_tol(v1, v2, v3, v4, p0, TETRAHEDRAL_CENTERINVOLUME_TOL)
}

fn same_side_tol(
    v1: [f64; 3],
    v2: [f64; 3],
    v3: [f64; 3],
    v4: [f64; 3],
    p0: [f64; 3],
    tol: f64,
) -> bool {
    let normal = cross_product(
        [v2[0] - v1[0], v2[1] - v1[1], v2[2] - v1[2]],
        [v3[0] - v1[0], v3[1] - v1[1], v3[2] - v1[2]],
    );
    let d1 =
        normal[0] * (v4[0] - v1[0]) + normal[1] * (v4[1] - v1[1]) + normal[2] * (v4[2] - v1[2]);
    let d2 =
        normal[0] * (p0[0] - v1[0]) + normal[1] * (p0[1] - v1[1]) + normal[2] * (p0[2] - v1[2]);
    if d1.abs() < tol || d2.abs() < tol {
        return false;
    }
    !((d1 < 0.0) ^ (d2 < 0.0))
}

// ============================================================================
// Double Bond Stereo
// ============================================================================

#[derive(Debug, Clone)]
struct StereoDoubleBond {
    atoms: [usize; 4],
    sign: i8,
}

#[derive(Debug, Clone)]
struct Atropisomer {
    _bond: (usize, usize),
    substituents: (usize, usize, usize, usize),
    _sign: f64,
    vol_lower: f64,
    vol_upper: f64,
}

fn find_stereo_bonds(
    mol: &Molecule,
) -> (
    Vec<(usize, usize, usize)>,
    Vec<StereoDoubleBond>,
    Vec<Atropisomer>,
) {
    let mut double_bond_ends = Vec::new();
    let mut stereo_db = Vec::new();
    let mut atropisomers = Vec::new();
    for bond in &mol.bonds {
        if bond.bond_type == BondType::Double {
            let a1 = bond.atom1;
            let a2 = bond.atom2;
            for &nbr in &mol.adjacency[a1] {
                if nbr != a2 && !is_double_bond(mol, a1, nbr) {
                    double_bond_ends.push((nbr, a1, a2));
                }
            }
            for &nbr in &mol.adjacency[a2] {
                if nbr != a1 && !is_double_bond(mol, a2, nbr) {
                    double_bond_ends.push((nbr, a2, a1));
                }
            }
            if bond.stereo != BondStereo::None {
                if let (Some(n1), Some(n2)) = (
                    find_non_double_neighbor(mol, a1, a2),
                    find_non_double_neighbor(mol, a2, a1),
                ) {
                    let sign = match bond.stereo {
                        BondStereo::Trans => 1,
                        BondStereo::Cis => -1,
                        _ => 1,
                    };
                    stereo_db.push(StereoDoubleBond {
                        atoms: [n1, a1, a2, n2],
                        sign,
                    });
                }
            }
        }

        // Atropisomer detection
        if bond.stereo == BondStereo::AtropCW || bond.stereo == BondStereo::AtropCCW {
            let a1 = bond.atom1;
            let a2 = bond.atom2;
            let nbrs1: Vec<usize> = mol.adjacency[a1]
                .iter()
                .filter(|&&n| n != a2)
                .copied()
                .collect();
            let nbrs2: Vec<usize> = mol.adjacency[a2]
                .iter()
                .filter(|&&n| n != a1)
                .copied()
                .collect();

            if nbrs1.len() >= 2 && nbrs2.len() >= 2 {
                let (vol_lower, vol_upper) = match bond.stereo {
                    BondStereo::AtropCW => (-100.0, -1.0),
                    BondStereo::AtropCCW => (1.0, 100.0),
                    _ => (0.0, 0.0),
                };
                atropisomers.push(Atropisomer {
                    _bond: (a1, a2),
                    substituents: (nbrs1[0], nbrs1[1], nbrs2[0], nbrs2[1]),
                    _sign: if bond.stereo == BondStereo::AtropCW {
                        -1.0
                    } else {
                        1.0
                    },
                    vol_lower,
                    vol_upper,
                });
            }
        }
    }
    (double_bond_ends, stereo_db, atropisomers)
}

fn dihedral_angle4(p0: [f64; 3], p1: [f64; 3], p2: [f64; 3], p3: [f64; 3]) -> f64 {
    let b1 = [p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]];
    let b2 = [p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2]];
    let b3 = [p3[0] - p2[0], p3[1] - p2[1], p3[2] - p2[2]];
    let n1 = cross_product(b1, b2);
    let n2 = cross_product(b2, b3);
    let b2_norm = (b2[0] * b2[0] + b2[1] * b2[1] + b2[2] * b2[2]).sqrt();
    if b2_norm < 1e-10 {
        return 0.0;
    }
    let m1 = cross_product(n1, [b2[0] / b2_norm, b2[1] / b2_norm, b2[2] / b2_norm]);
    let x = n1.dot(n2);
    let y = m1.dot(n2);
    y.atan2(x)
}

fn dihedral_angle(coords: &[[f64; 3]], i: usize, j: usize, k: usize, l: usize) -> f64 {
    dihedral_angle4(coords[i], coords[j], coords[k], coords[l])
}

fn check_double_bond_stereo(coords: &[[f64; 3]], stereo_dbs: &[StereoDoubleBond]) -> bool {
    for sdb in stereo_dbs {
        let dihedral = dihedral_angle(
            coords,
            sdb.atoms[0],
            sdb.atoms[1],
            sdb.atoms[2],
            sdb.atoms[3],
        );
        if (dihedral - std::f64::consts::FRAC_PI_2) * (sdb.sign as f64) < 0.0 {
            return false;
        }
    }
    true
}

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

fn double_bond_geometry_checks(
    coords: &[[f64; 3]],
    double_bond_ends: &[(usize, usize, usize)],
) -> bool {
    const LINEAR_TOL: f64 = 1e-3;
    for &(a0, a1, a2) in double_bond_ends {
        let v1 = normalize_vec([
            coords[a1][0] - coords[a0][0],
            coords[a1][1] - coords[a0][1],
            coords[a1][2] - coords[a0][2],
        ]);
        let v2 = normalize_vec([
            coords[a1][0] - coords[a2][0],
            coords[a1][1] - coords[a2][1],
            coords[a1][2] - coords[a2][2],
        ]);
        if v1.dot(v2) + 1.0 < LINEAR_TOL {
            return false;
        }
    }
    true
}

// ============================================================================
// Planarity
// ============================================================================

#[derive(Debug, Clone)]
struct PlanarityConstraints {
    impropers: Vec<(usize, usize, usize, usize, f64)>,
    ring_torsions: Vec<(usize, usize, usize, usize)>,
    exocyclic_torsions: Vec<(usize, usize, usize, usize)>,
    aromatic_atoms: HashSet<usize>,
}

fn build_planarity_constraints(mol: &Molecule) -> PlanarityConstraints {
    use crate::molecule::graph::{get_aromatic_atoms, get_neighbors};
    let aromatic_atoms = get_aromatic_atoms(mol);
    let rings = crate::molecule::graph::find_rings(mol);
    let mut impropers = Vec::new();
    let mut ring_torsions = Vec::new();
    let mut exocyclic_torsions = Vec::new();

    for ring in &rings {
        if !(4..=6).contains(&ring.len()) {
            continue;
        }
        let all_aromatic = ring.iter().all(|a| aromatic_atoms.contains(a));
        if all_aromatic {
            let rsize = ring.len();
            for start in 0..rsize {
                ring_torsions.push((
                    ring[start],
                    ring[(start + 1) % rsize],
                    ring[(start + 2) % rsize],
                    ring[(start + 3) % rsize],
                ));
            }
        }
    }

    // Iterate aromatic atoms in a DETERMINISTIC (sorted) order. aromatic_atoms is
    // a HashSet (Rust RandomState -> randomized per process); iterating it raw
    // made the planarity-constraint Vec order — and thus the embedding —
    // non-deterministic run-to-run (every aromatic molecule varied).
    let mut arom_atoms: Vec<usize> = aromatic_atoms.iter().copied().collect();
    arom_atoms.sort_unstable();
    for &atom_idx in &arom_atoms {
        let neighbors = get_neighbors(atom_idx, mol);
        if neighbors.len() >= 3 {
            let k_improper = improper_k_for_atom(atom_idx, mol);
            for i in 0..neighbors.len() {
                for j in (i + 1)..neighbors.len() {
                    for k in (j + 1)..neighbors.len() {
                        impropers.push((
                            atom_idx,
                            neighbors[i],
                            neighbors[j],
                            neighbors[k],
                            k_improper,
                        ));
                    }
                }
            }
        }
        let ring_neighbors: Vec<usize> = neighbors
            .iter()
            .filter(|&&n| aromatic_atoms.contains(&n))
            .copied()
            .collect();
        let exo_neighbors: Vec<usize> = neighbors
            .iter()
            .filter(|&&n| !aromatic_atoms.contains(&n))
            .copied()
            .collect();
        if ring_neighbors.len() >= 2 {
            for &exo in &exo_neighbors {
                exocyclic_torsions.push((exo, atom_idx, ring_neighbors[0], ring_neighbors[1]));
            }
        }
    }

    // Non-aromatic sp2 atoms: carboxyl, carbonyl, amide, etc.
    for atom_idx in 0..mol.atoms.len() {
        if aromatic_atoms.contains(&atom_idx) {
            continue; // already handled above
        }
        let hyb = crate::molecule::graph::determine_hybridization(atom_idx, mol);
        if !matches!(hyb, Hybridization::Sp2) {
            continue;
        }
        let neighbors = get_neighbors(atom_idx, mol);
        if neighbors.len() >= 3 {
            let k_improper = improper_k_for_atom(atom_idx, mol);
            for i in 0..neighbors.len() {
                for j in (i + 1)..neighbors.len() {
                    for k in (j + 1)..neighbors.len() {
                        impropers.push((
                            atom_idx,
                            neighbors[i],
                            neighbors[j],
                            neighbors[k],
                            k_improper,
                        ));
                    }
                }
            }
        }
        // Exocyclic torsions for sp2 atoms with double bonds
        // Keep substituents planar with the double bond
        let double_bond_neighbors: Vec<usize> = neighbors
            .iter()
            .filter(|&&n| {
                mol.bonds.iter().any(|b| {
                    (b.atom1 == atom_idx && b.atom2 == n || b.atom2 == atom_idx && b.atom1 == n)
                        && (b.bond_type == BondType::Double || b.bond_type == BondType::Aromatic)
                })
            })
            .copied()
            .collect();
        if !double_bond_neighbors.is_empty() && neighbors.len() >= 3 {
            let other_neighbors: Vec<usize> = neighbors
                .iter()
                .filter(|&&n| !double_bond_neighbors.contains(&n))
                .copied()
                .collect();
            if other_neighbors.len() >= 2 {
                for &db_nbr in &double_bond_neighbors {
                    exocyclic_torsions.push((
                        other_neighbors[0],
                        atom_idx,
                        db_nbr,
                        other_neighbors[1],
                    ));
                }
            }
        }
    }

    // Hydroxyl/amino H on an sp2 center (carboxyl OH, phenol/enol OH, amide NH,
    // sulfonamide NH): pin the H into the sp2 plane as an IMPROPER
    // (central=H, refs = sp2-C, its =O neighbor, and the heteroatom X), so
    // O=C-O-H / C=O-N-H stays planar. Done as an improper (not an exocyclic
    // torsion) because impropers use a central-difference gradient (reliable);
    // the dihedral-gradient path has a known sign bug (see dbg_dihedral_grad_fd).
    const K_IMPROPER_OH: f64 = 3.0;
    for atom_idx in 0..mol.atoms.len() {
        let sym = mol.atoms[atom_idx].symbol.as_str();
        if sym != "O" && sym != "N" && sym != "S" {
            continue;
        }
        let neighbors = get_neighbors(atom_idx, mol);
        if neighbors.len() != 2 {
            continue;
        }
        let h_idx = neighbors
            .iter()
            .find(|&&n| mol.atoms[n].symbol == "H")
            .copied();
        let c_idx = neighbors.iter().find(|&&n| {
            mol.atoms[n].symbol == "C"
                && matches!(
                    crate::molecule::graph::determine_hybridization(n, mol),
                    Hybridization::Sp2
                )
        });
        let (Some(h), Some(&c)) = (h_idx, c_idx) else { continue; };
        let ref_atom = mol.bonds.iter().find_map(|b| {
            let other = if b.atom1 == c {
                b.atom2
            } else if b.atom2 == c {
                b.atom1
            } else {
                return None;
            };
            if other == atom_idx {
                return None;
            }
            if b.bond_type == BondType::Double {
                Some(other)
            } else {
                None
            }
        });
        if let Some(r) = ref_atom {
            // improper (central=h, n1=c, n2=r(=O), n3=atom_idx(X)): out-of-plane of H
            impropers.push((h, c, r, atom_idx, K_IMPROPER_OH));
        }
    }

    PlanarityConstraints {
        impropers,
        ring_torsions,
        exocyclic_torsions,
        aromatic_atoms,
    }
}

fn improper_k_for_atom(atom_idx: usize, mol: &Molecule) -> f64 {
    let base_k = match mol.atoms[atom_idx].symbol.as_str() {
        "C" => 40.0,
        "N" => 30.0,
        "O" => 80.0,
        "S" => 20.0,
        _ => 40.0,
    };
    base_k * 10.0
}

fn out_of_plane_angle(coords: &[[f64; 3]], central: usize, n1: usize, n2: usize, n3: usize) -> f64 {
    let v1 = [
        coords[n1][0] - coords[central][0],
        coords[n1][1] - coords[central][1],
        coords[n1][2] - coords[central][2],
    ];
    let v2 = [
        coords[n2][0] - coords[central][0],
        coords[n2][1] - coords[central][1],
        coords[n2][2] - coords[central][2],
    ];
    let v3 = [
        coords[n3][0] - coords[central][0],
        coords[n3][1] - coords[central][1],
        coords[n3][2] - coords[central][2],
    ];
    let normal = cross_product(
        [v1[0] - v2[0], v1[1] - v2[1], v1[2] - v2[2]],
        [v1[0] - v3[0], v1[1] - v3[1], v1[2] - v3[2]],
    );
    let normal_norm =
        (normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]).sqrt();
    if normal_norm < 1e-10 {
        return 0.0;
    }
    let v1_norm = (v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2]).sqrt();
    if v1_norm < 1e-10 {
        return 0.0;
    }
    let dot = (v1[0] * normal[0] + v1[1] * normal[1] + v1[2] * normal[2]) / normal_norm;
    dot.abs().asin()
}

fn planarity_energy(coords: &[[f64; 3]], pc: &PlanarityConstraints) -> f64 {
    const K_RING_TOR: f64 = 10.0;
    const K_EXO_TOR: f64 = 2.0;
    let mut energy = 0.0;
    for &(central, n1, n2, n3, k_imp) in &pc.impropers {
        let chi = out_of_plane_angle(coords, central, n1, n2, n3);
        energy += k_imp * chi * chi;
    }
    for &(i, j, k, l) in &pc.ring_torsions {
        let phi = dihedral_angle(coords, i, j, k, l);
        energy += K_RING_TOR * (1.0 - (2.0 * phi).cos());
    }
    for &(i, j, k, l) in &pc.exocyclic_torsions {
        let phi = dihedral_angle(coords, i, j, k, l);
        let phi_wrapped = if phi > std::f64::consts::FRAC_PI_2 {
            phi - std::f64::consts::PI
        } else if phi < -std::f64::consts::FRAC_PI_2 {
            phi + std::f64::consts::PI
        } else {
            phi
        };
        energy += K_EXO_TOR * phi_wrapped * phi_wrapped;
    }
    energy
}

fn check_planarity(
    coords: &[[f64; 3]],
    mol: &Molecule,
    pc: &PlanarityConstraints,
    _threshold: f64,
) -> bool {
    if pc.aromatic_atoms.is_empty() {
        return true;
    }
    let n_atoms = mol.atoms.len() as f64;
    let e_per_atom = planarity_energy(coords, pc) / n_atoms;
    e_per_atom < PLANARITY_ENERGY_TOL
}

/// Check for severe van-der-Waals clashes between non-bonded atoms.
/// Rejects conformers where any non-bonded pair is closer than 60% of the VDW sum.
fn has_vdw_clash(coords: &[[f64; 3]], mol: &Molecule) -> bool {
    let n = coords.len();
    for i in 0..n {
        for j in (i + 1)..n {
            // Skip directly bonded atoms (1-2)
            if mol.adjacency[i].contains(&j) {
                continue;
            }
            let dx = coords[i][0] - coords[j][0];
            let dy = coords[i][1] - coords[j][1];
            let dz = coords[i][2] - coords[j][2];
            let d = (dx * dx + dy * dy + dz * dz).sqrt();
            let ri = vdw_radius(&mol.atoms[i].symbol);
            let rj = vdw_radius(&mol.atoms[j].symbol);
            if d < 0.6 * (ri + rj) {
                return true;
            }
        }
    }
    false
}

const BOND_LENGTH_TOLERANCE: f64 = 0.30;

fn bond_lengths_reasonable(coords: &[[f64; 3]], mol: &Molecule) -> bool {
    for bond in &mol.bonds {
        let i = bond.atom1;
        let j = bond.atom2;
        let dx = coords[i][0] - coords[j][0];
        let dy = coords[i][1] - coords[j][1];
        let dz = coords[i][2] - coords[j][2];
        let actual = (dx * dx + dy * dy + dz * dz).sqrt();
        let hyb_i = crate::molecule::graph::determine_hybridization(i, mol);
        let hyb_j = crate::molecule::graph::determine_hybridization(j, mol);
        let aro_i = crate::molecule::graph::is_aromatic(i, mol);
        let aro_j = crate::molecule::graph::is_aromatic(j, mol);
        let expected = compute_bond_length(
            &mol.atoms[i].symbol,
            hyb_i,
            aro_i,
            &mol.atoms[j].symbol,
            hyb_j,
            aro_j,
            bond.bond_type,
        );
        if (actual - expected).abs() > BOND_LENGTH_TOLERANCE {
            return false;
        }
    }
    true
}

fn flatten_aromatic_rings(coords: &mut [[f64; 3]], mol: &Molecule, pc: &PlanarityConstraints) {
    let mut visited = HashSet::new();
    let mut components: Vec<Vec<usize>> = Vec::new();
    // Deterministic iteration order (see build_planarity_constraints).
    let mut arom_starts: Vec<usize> = pc.aromatic_atoms.iter().copied().collect();
    arom_starts.sort_unstable();
    for &start in &arom_starts {
        if visited.contains(&start) {
            continue;
        }
        let mut comp = Vec::new();
        let mut queue = std::collections::VecDeque::new();
        queue.push_back(start);
        visited.insert(start);
        comp.push(start);
        while let Some(idx) = queue.pop_front() {
            for &nbr in &mol.adjacency[idx] {
                if !pc.aromatic_atoms.contains(&nbr) || visited.contains(&nbr) {
                    continue;
                }
                visited.insert(nbr);
                queue.push_back(nbr);
                comp.push(nbr);
            }
        }
        if comp.len() >= 3 {
            components.push(comp);
        }
    }
    for component in &components {
        let cx = component.iter().map(|&i| coords[i][0]).sum::<f64>() / component.len() as f64;
        let cy = component.iter().map(|&i| coords[i][1]).sum::<f64>() / component.len() as f64;
        let cz = component.iter().map(|&i| coords[i][2]).sum::<f64>() / component.len() as f64;
        let mut cov = [[0.0f64; 3]; 3];
        for &idx in component {
            let dx = coords[idx][0] - cx;
            let dy = coords[idx][1] - cy;
            let dz = coords[idx][2] - cz;
            cov[0][0] += dx * dx;
            cov[0][1] += dx * dy;
            cov[0][2] += dx * dz;
            cov[1][1] += dy * dy;
            cov[1][2] += dy * dz;
            cov[2][2] += dz * dz;
        }
        let n = component.len() as f64;
        for row in cov.iter_mut() {
            for val in row.iter_mut() {
                *val /= n;
            }
        }
        cov[1][0] = cov[0][1];
        cov[2][0] = cov[0][2];
        cov[2][1] = cov[1][2];
        let normal = eigenvector_smallest_eigenvalue_3x3(&cov);
        let mut atoms_to_flatten: Vec<usize> = component.clone();
        // First pass: add direct neighbors of aromatic atoms
        for &ring_atom in component {
            for &nbr in &mol.adjacency[ring_atom] {
                if !pc.aromatic_atoms.contains(&nbr) && !atoms_to_flatten.contains(&nbr) {
                    atoms_to_flatten.push(nbr);
                }
            }
        }
        // Second pass: add H atoms bonded to flattened atoms that should be planar
        // (e.g., H atoms on aniline N, which is conjugated to the aromatic ring)
        let flattened_non_arom: Vec<usize> = atoms_to_flatten
            .iter()
            .filter(|&&a| !pc.aromatic_atoms.contains(&a))
            .copied()
            .collect();
        for &flattened_atom in &flattened_non_arom {
            // Include H neighbors of conjugated heteroatoms (N, O) adjacent to aromatic rings
            let sym = &mol.atoms[flattened_atom].symbol;
            if *sym == "N" || *sym == "O" {
                for &nbr in &mol.adjacency[flattened_atom] {
                    if mol.atoms[nbr].symbol == "H" && !atoms_to_flatten.contains(&nbr) {
                        atoms_to_flatten.push(nbr);
                    }
                }
            }
        }
        for &idx in &atoms_to_flatten {
            let dist = (coords[idx][0] - cx) * normal[0]
                + (coords[idx][1] - cy) * normal[1]
                + (coords[idx][2] - cz) * normal[2];
            coords[idx][0] -= dist * normal[0];
            coords[idx][1] -= dist * normal[1];
            coords[idx][2] -= dist * normal[2];
        }
    }
}

pub fn eigenvector_smallest_eigenvalue_3x3(m: &[[f64; 3]; 3]) -> [f64; 3] {
    let c1 = m[0][0] + m[1][1] + m[2][2];
    let c2 = (m[0][0] * m[1][1] - m[0][1].powi(2))
        + (m[1][1] * m[2][2] - m[1][2].powi(2))
        + (m[0][0] * m[2][2] - m[0][2].powi(2));
    let c3 = m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][2] * m[2][0] - m[1][0] * m[2][2])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
    let mut lambda = m[0][0].min(m[1][1]).min(m[2][2]);
    for _ in 0..20 {
        let f = lambda.powi(3) - c1 * lambda.powi(2) + c2 * lambda - c3;
        let fp = 3.0 * lambda.powi(2) - 2.0 * c1 * lambda + c2;
        if fp.abs() < 1e-12 {
            break;
        }
        lambda -= f / fp;
    }
    let mut mat = [[0.0f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            mat[i][j] = m[i][j];
        }
        mat[i][i] -= lambda;
    }
    let mut v = [
        mat[0][1] * mat[1][2] - mat[0][2] * mat[1][1],
        mat[0][2] * mat[1][0] - mat[0][0] * mat[1][2],
        mat[0][0] * mat[1][1] - mat[0][1] * mat[1][0],
    ];
    let norm = (v[0].powi(2) + v[1].powi(2) + v[2].powi(2)).sqrt();
    if norm > 1e-10 {
        v[0] /= norm;
        v[1] /= norm;
        v[2] /= norm;
    } else {
        v = [0.0, 0.0, 1.0];
    }
    v
}

// ============================================================================
// Torsion Preferences
// ============================================================================

#[derive(Debug, Clone)]
struct TorsionPreference {
    i: usize,
    j: usize,
    k: usize,
    l: usize,
    signs: [i32; 6],
    v: [f64; 6],
}

fn sym(mol: &Molecule, a: usize) -> &str {
    &mol.atoms[a].symbol
}
fn get_hyb(mol: &Molecule, a: usize) -> Hybridization {
    crate::molecule::graph::determine_hybridization(a, mol)
}
fn is_arom(mol: &Molecule, a: usize) -> bool {
    crate::molecule::graph::is_aromatic(a, mol)
}
fn h_count(mol: &Molecule, a: usize) -> usize {
    mol.adjacency[a]
        .iter()
        .filter(|&&n| mol.atoms[n].symbol == "H")
        .count()
}
fn has_db_to_o(mol: &Molecule, c: usize) -> bool {
    mol.adjacency[c]
        .iter()
        .any(|&nbr| mol.atoms[nbr].symbol == "O" && is_double_bond(mol, c, nbr))
}
fn has_db_to_n(mol: &Molecule, c: usize) -> bool {
    mol.adjacency[c]
        .iter()
        .any(|&nbr| mol.atoms[nbr].symbol == "N" && is_double_bond(mol, c, nbr))
}
fn is_carbonyl_c(mol: &Molecule, c: usize) -> bool {
    sym(mol, c) == "C" && has_db_to_o(mol, c)
}
fn is_sulfonyl_s(mol: &Molecule, s: usize) -> bool {
    if sym(mol, s) != "S" {
        return false;
    }
    let o_dbls = mol.adjacency[s]
        .iter()
        .filter(|&&nbr| mol.atoms[nbr].symbol == "O" && is_double_bond(mol, s, nbr))
        .count();
    o_dbls >= 2
}
fn is_nitro_n(mol: &Molecule, n: usize) -> bool {
    if sym(mol, n) != "N" {
        return false;
    }
    let o_dbls = mol.adjacency[n]
        .iter()
        .filter(|&&nbr| mol.atoms[nbr].symbol == "O" && is_double_bond(mol, n, nbr))
        .count();
    o_dbls >= 1
        && mol.adjacency[n]
            .iter()
            .any(|&nbr| mol.atoms[nbr].symbol == "O" && !is_double_bond(mol, n, nbr))
}
fn bond_in_ring(rings: &[Vec<usize>], a: usize, b: usize) -> bool {
    rings.iter().any(|r| r.contains(&a) && r.contains(&b))
}
fn ring_size_of_bond(rings: &[Vec<usize>], a: usize, b: usize) -> Option<usize> {
    rings
        .iter()
        .find(|r| r.contains(&a) && r.contains(&b))
        .map(|r| r.len())
}

/// Match a torsion pattern for atoms a1-a2-a3-a4.
/// Returns Fourier coefficients (signs, V) if a pattern matches.
/// Patterns are checked in order of specificity — first match wins.
fn is_nh1(mol: &Molecule, a: usize) -> bool {
    sym(mol, a) == "N" && h_count(mol, a) == 1
}
fn is_nh0(mol: &Molecule, a: usize) -> bool {
    sym(mol, a) == "N" && h_count(mol, a) == 0
}
fn is_ch2(mol: &Molecule, a: usize) -> bool {
    sym(mol, a) == "C" && matches!(get_hyb(mol, a), Hybridization::Sp3) && h_count(mol, a) == 2
}
#[allow(dead_code)]
fn is_ch1(mol: &Molecule, a: usize) -> bool {
    sym(mol, a) == "C" && matches!(get_hyb(mol, a), Hybridization::Sp3) && h_count(mol, a) == 1
}
fn is_ch0_sp3(mol: &Molecule, a: usize) -> bool {
    sym(mol, a) == "C" && matches!(get_hyb(mol, a), Hybridization::Sp3) && h_count(mol, a) == 0
}
fn is_sx2(mol: &Molecule, a: usize) -> bool {
    sym(mol, a) == "S" && mol.adjacency[a].len() == 2
}
fn is_sx3(mol: &Molecule, a: usize) -> bool {
    sym(mol, a) == "S" && mol.adjacency[a].len() == 3
}
fn is_sx4(mol: &Molecule, a: usize) -> bool {
    sym(mol, a) == "S" && mol.adjacency[a].len() == 4
}
fn has_cao_nbr(mol: &Molecule, a: usize) -> bool {
    mol.adjacency[a]
        .iter()
        .any(|&nbr| (nbr != a) && is_carbonyl_c(mol, nbr))
}
#[allow(dead_code)]
fn is_ester_link(mol: &Molecule, a: usize) -> bool {
    sym(mol, a) == "O" && mol.adjacency[a].iter().any(|&nbr| is_carbonyl_c(mol, nbr))
}

fn match_torsion_pattern(
    mol: &Molecule,
    a1: usize,
    a2: usize,
    a3: usize,
    a4: usize,
    rings: &[Vec<usize>],
    et_version: u32,
) -> Option<([i32; 6], [f64; 6])> {
    let s1 = sym(mol, a1);
    let s2 = sym(mol, a2);
    let s3 = sym(mol, a3);
    let s4 = sym(mol, a4);
    let h1 = get_hyb(mol, a1);
    let h2 = get_hyb(mol, a2);
    let h3 = get_hyb(mol, a3);
    let h4 = get_hyb(mol, a4);
    let ar1 = is_arom(mol, a1);
    let ar2 = is_arom(mol, a2);
    let ar3 = is_arom(mol, a3);
    let ar4 = is_arom(mol, a4);
    let hc1 = h_count(mol, a1);
    let hc2 = h_count(mol, a2);
    let hc3 = h_count(mol, a3);
    let hc4 = h_count(mol, a4);
    let b23_ring = bond_in_ring(rings, a2, a3);

    let sp2 = |h: Hybridization| matches!(h, Hybridization::Sp2);
    let sp3 = |h: Hybridization| matches!(h, Hybridization::Sp3);

    // ETversion-specific parameter adjustments
    let sp3_sp3_v3 = if et_version >= 2 { 5.0 } else { 7.0 };
    let sp3_sp3_ring_v3 = if et_version >= 2 { 1.5 } else { 2.0 };

    if b23_ring {
        if let Some(rsize) = ring_size_of_bond(rings, a2, a3) {
            if rsize >= MIN_MACROCYCLE_SIZE {
                // Macrocycle ring torsions
                if is_amide_bond(mol, a2, a3) {
                    return Some(([1, -1, 1, 1, 1, 1], [0.0, 8.0, 0.0, 0.0, 0.0, 0.0]));
                }
                if sp2(h1) && sp2(h2) && sp2(h3) && sp2(h4) {
                    return Some(([1, -1, 1, 1, 1, 1], [0.0, 100.0, 0.0, 0.0, 0.0, 0.0]));
                }
                if sp3(h2) && sp3(h3) {
                    return Some(([1, 0, 1, 1, 1, 1], [0.0, 0.0, 7.0, 0.0, 0.0, 0.0]));
                }
                if (sp2(h2) && sp3(h3)) || (sp3(h2) && sp2(h3)) {
                    return Some(([1, 0, 1, 1, 1, 1], [0.0, 0.0, 5.0, 0.0, 0.0, 0.0]));
                }
                return None;
            }
            if (3..=8).contains(&rsize) {
                if sp2(h1) && sp2(h2) && sp2(h3) && sp2(h4) {
                    return Some(([1, -1, 1, 1, 1, 1], [0.0, 10.0, 0.0, 0.0, 0.0, 0.0]));
                }
                if rsize == 3 || rsize == 4 {
                    return Some(([-1, 1, 1, 1, 1, 1], [30.0, 0.0, 0.0, 0.0, 0.0, 0.0]));
                }
                if rsize == 5 {
                    if sp3(h2) && sp3(h3) {
                        return Some(([1, 0, 1, 1, 1, 1], [0.0, 0.0, 30.0, 0.0, 0.0, 0.0]));
                    }
                    if sp2(h2) && sp3(h3) {
                        return Some(([1, 0, 1, 1, 1, -1], [0.0, 0.0, 0.0, 0.0, 0.0, 15.0]));
                    }
                    if sp3(h2) && sp2(h3) {
                        return Some(([1, 0, 1, 1, 1, -1], [0.0, 0.0, 0.0, 0.0, 0.0, 15.0]));
                    }
                    return Some(([1, 0, 1, 1, 1, 1], [0.0, 0.0, 5.0, 0.0, 0.0, 0.0]));
                }
                // 6-8 ring
                if sp3(h2) && sp3(h3) {
                    return Some(([-1, 1, 1, 1, 1, 1], [20.0, 0.0, 0.0, 0.0, 0.0, 0.0]));
                }
                return Some(([1, 0, 1, 1, 1, 1], [0.0, 0.0, 5.0, 0.0, 0.0, 0.0]));
            }
        }
        return None;
    }

    // ---- Chain (non-ring) patterns below ----

    // 1. Conjugated C=C=C (cumulated double) → cis
    if is_double_bond(mol, a2, a3) {
        let a2_has_other_db = mol.adjacency[a2]
            .iter()
            .any(|&nbr| nbr != a3 && is_double_bond(mol, a2, nbr));
        let a3_has_other_db = mol.adjacency[a3]
            .iter()
            .any(|&nbr| nbr != a2 && is_double_bond(mol, a3, nbr));
        if a2_has_other_db || a3_has_other_db {
            return Some(([1, -1, 1, 1, 1, 1], [0.0, 100.0, 0.0, 0.0, 0.0, 0.0]));
        }
    }

    // 2. Conjugated C=C-C=C diene: central bond SINGLE between two sp2 carbons,
    //    each carrying its own C=C -> trans + planar. (Old comment mislabeled the
    //    cumulated case below as "conjugated diene".) The torsion-snap basin-hop
    //    in generate_initial_coords_with_config now crosses the ~90° barrier that
    //    previously left this pref unable to take effect.
    let a2_has_cc = mol.adjacency[a2]
        .iter()
        .any(|&nb| nb != a3 && mol.atoms[nb].symbol == "C" && is_double_bond(mol, a2, nb));
    let a3_has_cc = mol.adjacency[a3]
        .iter()
        .any(|&nb| nb != a2 && mol.atoms[nb].symbol == "C" && is_double_bond(mol, a3, nb));
    if !is_double_bond(mol, a2, a3)
        && s2 == "C"
        && s3 == "C"
        && sp2(h2)
        && sp2(h3)
        && !ar2
        && !ar3
        && a2_has_cc
        && a3_has_cc
    {
        // k=1 favors trans (180°); k=2 enforces planar (0/180°).
        return Some(([1, -1, 1, 1, 1, 1], [2.0, 4.0, 0.0, 0.0, 0.0, 0.0]));
    }

    // 2b. Cumulated C=C=C (central bond double, both sp2) -> cis/planar
    if is_double_bond(mol, a2, a3) && sp2(h2) && sp2(h3) {
        return Some(([1, -1, 1, 1, 1, 1], [0.0, 100.0, 0.0, 0.0, 0.0, 0.0]));
    }

    // 3. C=C-CH2-O (vinyl ether) → gauche
    if is_double_bond(mol, a2, a3) && sp2(h2) && sp3(h3) && s3 == "C" && hc3 >= 1 {
        if s4 == "O" {
            return Some(([1, 0, 1, -1, 1, 1], [0.0, 0.0, 0.0, 20.0, 0.0, 0.0]));
        }
        if s4 == "C" {
            return Some(([1, 0, 1, -1, 1, 1], [0.0, 0.0, 0.0, 1.5, 0.0, 0.0]));
        }
        return Some(([1, 0, 1, -1, 1, 1], [0.0, 0.0, 0.0, 4.0, 0.0, 0.0]));
    }

    // 4. C=C-CH2-C (vinyl to sp3 chain)
    if is_double_bond(mol, a2, a3) && sp2(h2) && sp3(h3) {
        return Some(([1, 0, 1, -1, 1, 1], [0.0, 0.0, 0.0, 4.0, 0.0, 0.0]));
    }

    // 5. Ester: O=C-O-* (specific by end atom)
    if is_carbonyl_c(mol, a2) && s3 == "O" {
        if s4 == "C" && is_ch0_sp3(mol, a4) {
            return Some(([1, -1, 1, 1, 1, 1], [0.0, 78.2, 0.0, 0.0, 0.0, 0.0]));
        }
        if s4 == "C" {
            return Some(([1, -1, 1, 1, 1, 1], [0.0, 100.0, 0.0, 0.0, 0.0, 0.0]));
        }
        return Some(([1, -1, 1, 1, 1, 1], [0.0, 100.0, 0.0, 0.0, 0.0, 0.0]));
    }
    if s2 == "O" && is_carbonyl_c(mol, a3) {
        return Some(([1, -1, 1, 1, 1, 1], [0.0, 100.0, 0.0, 0.0, 0.0, 0.0]));
    }

    // 6. Amide: C(=O)-NH1 → strong trans
    if is_carbonyl_c(mol, a2) && is_nh1(mol, a3) {
        if has_cao_nbr(mol, a3) {
            return Some(([1, -1, 1, 1, 1, 1], [0.0, 100.0, 0.0, 0.0, 0.0, 0.0]));
        }
        return Some(([1, -1, 1, 1, 1, 1], [0.0, 100.0, 0.0, 0.0, 0.0, 0.0]));
    }
    if is_nh1(mol, a2) && is_carbonyl_c(mol, a3) {
        return Some(([1, -1, 1, 1, 1, 1], [0.0, 100.0, 0.0, 0.0, 0.0, 0.0]));
    }

    // 7. Amide: C(=O)-N(C=O) → trans
    if is_carbonyl_c(mol, a2) && is_nh0(mol, a3) && has_cao_nbr(mol, a3) {
        return Some(([1, -1, 1, 1, 1, 1], [0.0, 10.0, 0.0, 0.0, 0.0, 0.0]));
    }

    // 8. Amide: C(=O)-NH0-Ar → planar
    if is_carbonyl_c(mol, a2) && is_nh0(mol, a3) && ar4 {
        return Some(([1, -1, 1, 1, 1, 1], [0.0, 8.0, 0.0, 0.0, 0.0, 0.0]));
    }

    // 9. Amide side chain: C(=O)-N-CH2-*
    if is_carbonyl_c(mol, a2) && s3 == "N" && is_ch2(mol, a3) {
        if is_nh0(mol, a3) {
            return Some(([1, -1, 1, 1, 1, 1], [0.0, 15.0, 0.0, 0.0, 0.0, 0.0]));
        }
        return Some(([1, 0, 1, 1, 1, 1], [0.0, 0.0, 2.0, 0.0, 0.0, 0.0]));
    }
    if s2 == "N" && is_carbonyl_c(mol, a3) && is_ch2(mol, a2) {
        if is_nh0(mol, a2) {
            return Some(([1, -1, 1, 1, 1, 1], [0.0, 15.0, 0.0, 0.0, 0.0, 0.0]));
        }
        return Some(([1, 0, 1, 1, 1, 1], [0.0, 0.0, 2.0, 0.0, 0.0, 0.0]));
    }

    // 10. Sulfonamide: S(=O)2-N
    if is_sulfonyl_s(mol, a2) && s3 == "N" {
        if is_nh1(mol, a3) && !is_ch2(mol, a3) {
            return Some(([1, -1, 1, 1, 1, 1], [0.0, 16.0, 5.0, 7.0, 0.0, 0.0]));
        }
        return Some(([1, -1, 1, 1, 1, 1], [0.0, 16.0, 5.0, 7.0, 0.0, 0.0]));
    }
    if s2 == "N" && is_sulfonyl_s(mol, a3) {
        return Some(([1, -1, 1, 1, 1, 1], [0.0, 16.0, 5.0, 7.0, 0.0, 0.0]));
    }

    // 11. Sulfonamide: S(=O)2-N-CH2-*
    if is_sulfonyl_s(mol, a2) && s3 == "N" && is_ch2(mol, a3) {
        if is_nh0(mol, a3) {
            return Some(([1, 0, 1, 1, 1, 1], [0.0, 0.0, 1.9, 0.0, 0.0, 0.0]));
        }
        return Some(([1, -1, 1, 1, 1, 1], [17.9, -13.3, 9.2, 4.7, -2.3, -0.9]));
    }
    if s2 == "N" && is_sulfonyl_s(mol, a3) && is_ch2(mol, a2) {
        if is_nh0(mol, a2) {
            return Some(([1, 0, 1, 1, 1, 1], [0.0, 0.0, 1.9, 0.0, 0.0, 0.0]));
        }
        return Some(([1, -1, 1, 1, 1, 1], [17.9, -13.3, 9.2, 4.7, -2.3, -0.9]));
    }

    // 12. Nitro: Ar-NO2
    if ar2 && is_nitro_n(mol, a3) {
        return Some(([1, -1, 1, 1, 1, 1], [0.0, 11.8, 0.0, 0.0, 0.0, 0.0]));
    }
    if is_nitro_n(mol, a2) && ar3 {
        return Some(([1, -1, 1, 1, 1, 1], [0.0, 11.8, 0.0, 0.0, 0.0, 0.0]));
    }

    // 13. Biaryl: Ar-Ar (non-ring)
    if ar2 && ar3 {
        // Differentiated by aromatic substitution patterns
        if !ar1 && !ar4 {
            return Some(([1, 0, 1, 0, 1, -1], [0.1, 0.9, 0.1, -0.4, 0.0, 0.3]));
        }
        return Some(([1, -1, 1, 1, 1, 1], [0.0, 30.0, 0.0, 0.0, 0.0, 0.0]));
    }

    // 14. Aromatic-O-C (ether to aromatic)
    if ar2 && s3 == "O" {
        if s4 == "S" {
            return Some(([1, -1, 1, 1, 1, 1], [0.0, 3.2, 0.0, 0.0, 0.0, 0.0]));
        }
        if ar4 {
            return Some(([1, 0, -1, 1, 1, 1], [0.0, 0.0, 3.0, 0.0, 0.0, 0.0]));
        }
        if is_carbonyl_c(mol, a4) {
            return Some(([1, -1, 1, 1, 1, 1], [0.0, 0.8, 0.0, 0.0, 0.0, 0.0]));
        }
        return Some(([1, 0, 1, 1, 1, 1], [0.0, 0.0, 4.0, 0.0, 0.0, 0.0]));
    }
    if s2 == "O" && ar3 {
        if s1 == "S" {
            return Some(([1, -1, 1, 1, 1, 1], [0.0, 3.2, 0.0, 0.0, 0.0, 0.0]));
        }
        if ar1 {
            return Some(([1, 0, -1, 1, 1, 1], [0.0, 0.0, 3.0, 0.0, 0.0, 0.0]));
        }
        return Some(([1, 0, 1, 1, 1, 1], [0.0, 0.0, 4.0, 0.0, 0.0, 0.0]));
    }

    // 15. Aromatic-S-C (thioether to aromatic)
    if ar2 && s3 == "S" && !is_sx3(mol, a3) && !is_sx4(mol, a3) {
        return Some(([1, -1, 1, 1, 1, 1], [0.0, 39.4, 0.0, 0.0, 0.0, 0.0]));
    }

    // 16. S-S (disulfide)
    if is_sx2(mol, a2) && is_sx2(mol, a3) {
        return Some(([1, -1, 1, 1, 1, 1], [0.0, 12.9, 0.0, 0.0, 0.0, 0.0]));
    }

    // 17. C-SX3 (sulfoxide/sulfone)
    if sp2(h2) && is_sx3(mol, a3) {
        return Some(([1, -1, 1, 1, 1, 1], [0.0, 16.0, 5.0, 7.0, 0.0, 0.0]));
    }
    if sp3(h2) && is_sx3(mol, a3) {
        return Some(([1, 0, 1, 1, 1, 1], [0.0, 0.0, 10.4, 0.0, 0.0, 0.0]));
    }
    if is_sx3(mol, a2) && sp3(h3) {
        return Some(([1, 0, 1, 1, 1, 1], [0.0, 0.0, 10.4, 0.0, 0.0, 0.0]));
    }
    if sp2(h2) && is_sx4(mol, a3) {
        return Some(([1, -1, 1, 1, 1, 1], [0.0, 3.8, 0.0, 0.0, 0.0, 0.0]));
    }
    if sp3(h2) && is_sx4(mol, a3) {
        return Some(([1, 0, 1, 1, 1, 1], [0.0, 0.0, 5.0, 0.0, 0.0, 0.0]));
    }

    // 18. Aromatic-C(=O) (carbonyl next to aromatic)
    if ar2 && is_carbonyl_c(mol, a3) {
        return Some(([1, -1, 1, 1, 1, 1], [0.0, 8.0, 0.0, 0.0, 0.0, 0.0]));
    }
    if is_carbonyl_c(mol, a2) && ar3 {
        return Some(([1, -1, 1, 1, 1, 1], [0.0, 8.0, 0.0, 0.0, 0.0, 0.0]));
    }

    // 19. Aromatic-C(=N) (C=N next to aromatic)
    if ar2 && has_db_to_n(mol, a3) && s3 == "C" {
        return Some(([1, -1, 1, 1, 1, 1], [0.0, 15.0, 0.0, 0.0, 0.0, 0.0]));
    }

    // 20. C(sp3)-C(=O) (carbonyl-adjacent sp3)
    if sp3(h2) && is_carbonyl_c(mol, a3) {
        if is_ch2(mol, a2) {
            return Some(([1, 0, 1, 0, 1, 0], [0.0, -0.7, 0.0, -0.5, 0.0, -0.6]));
        }
        return Some(([1, 0, -1, 1, 1, 1], [0.0, 0.0, 1.0, 0.0, 0.0, 0.0]));
    }
    if is_carbonyl_c(mol, a2) && sp3(h3) {
        if is_ch2(mol, a3) {
            return Some(([1, 0, 1, 0, 1, 0], [0.0, -0.7, 0.0, -0.5, 0.0, -0.6]));
        }
        return Some(([1, 0, -1, 1, 1, 1], [0.0, 0.0, 1.0, 0.0, 0.0, 0.0]));
    }

    // 21. Alkyl ether: *-CH2-O-C
    if is_ch2(mol, a2) && s3 == "O" && sp3(h4) {
        return Some(([1, 0, 1, 1, 1, 1], [0.0, 0.0, 2.5, 0.0, 0.0, 0.0]));
    }
    if s2 == "O" && is_ch2(mol, a3) && sp3(h1) {
        return Some(([1, 0, 1, 1, 1, 1], [0.0, 0.0, 2.5, 0.0, 0.0, 0.0]));
    }
    // General ether
    if sp3(h2) && s3 == "O" {
        return Some(([1, -1, 1, 1, 1, 1], [0.0, 8.0, 0.0, 0.0, 0.0, 0.0]));
    }
    if s2 == "O" && sp3(h3) {
        return Some(([1, -1, 1, 1, 1, 1], [0.0, 8.0, 0.0, 0.0, 0.0, 0.0]));
    }

    // 22. Alkyl amine: *-CH2-N-C
    if sp3(h2) && s3 == "N" {
        if is_ch2(mol, a2) && is_ch2(mol, a3) {
            return Some(([1, 1, 1, -1, 1, 1], [4.0, 3.1, 3.9, -0.8, 0.0, 0.7]));
        }
        return Some(([1, 0, 1, 1, 1, 1], [0.0, 0.0, 1.0, 0.0, 0.0, 0.0]));
    }
    if s2 == "N" && sp3(h3) {
        if is_ch2(mol, a2) && is_ch2(mol, a3) {
            return Some(([1, 1, 1, -1, 1, 1], [4.0, 3.1, 3.9, -0.8, 0.0, 0.7]));
        }
        return Some(([1, 0, 1, 1, 1, 1], [0.0, 0.0, 1.0, 0.0, 0.0, 0.0]));
    }

    // 23. sp3-sp3 general C-C single bond
    if sp3(h2) && sp3(h3) {
        // *-CH2-CH2-* (ethylene bridge)
        if is_ch2(mol, a2) && is_ch2(mol, a3) {
            return Some((
                [1, 0, 1, 1, 1, 1],
                [0.0, 0.0, sp3_sp3_v3 - 3.0, 0.0, 0.0, 0.0],
            ));
        }
        // Ring sp3
        let a2_in_ring = rings.iter().any(|r| r.contains(&a2));
        let a3_in_ring = rings.iter().any(|r| r.contains(&a3));
        if a2_in_ring || a3_in_ring {
            return Some((
                [1, 0, 1, 1, 1, 1],
                [0.0, 0.0, sp3_sp3_ring_v3, 0.0, 0.0, 0.0],
            ));
        }
        return Some(([1, 0, 1, 1, 1, 1], [0.0, 0.0, sp3_sp3_v3, 0.0, 0.0, 0.0]));
    }

    // 24. Aromatic-C(sp3) (aryl-alkyl)
    if ar2 && sp3(h3) {
        return Some(([1, -1, 1, 1, 1, 1], [0.0, 3.0, 0.0, 0.0, 0.0, 0.0]));
    }
    if sp3(h2) && ar3 {
        return Some(([1, -1, 1, 1, 1, 1], [0.0, 3.0, 0.0, 0.0, 0.0, 0.0]));
    }

    // 25. sp2-sp3 general
    if sp2(h2) && sp3(h3) {
        return Some(([1, 0, 1, 1, 1, 1], [0.0, 0.0, 1.9, 0.0, 0.0, 0.0]));
    }
    if sp3(h2) && sp2(h3) {
        return Some(([1, 0, 1, 1, 1, 1], [0.0, 0.0, 1.9, 0.0, 0.0, 0.0]));
    }

    // 26. Catch-all double bond
    if is_double_bond(mol, a2, a3) {
        return Some(([1, -1, 1, 1, 1, 1], [0.0, 100.0, 0.0, 0.0, 0.0, 0.0]));
    }

    None
}

fn build_torsion_preferences(mol: &Molecule, et_version: u32) -> Vec<TorsionPreference> {
    let rings = crate::molecule::graph::find_rings(mol);
    let mut prefs = Vec::new();

    // ETversion-specific parameter adjustments
    let sp3_sp3_v3 = if et_version >= 2 { 5.0 } else { 7.0 };
    let sp3_sp3_ring_v3 = if et_version >= 2 { 1.5 } else { 2.0 };
    let mut done_bonds = vec![false; mol.bonds.len()];

    // Bridged ring exclusion: bonds shared by >1 ring where at least one is small (<9)
    // should be excluded from torsion preferences
    let mut excluded_bonds = vec![false; mol.bonds.len()];
    for bi in 0..rings.len() {
        for bj in (bi + 1)..rings.len() {
            if rings[bi].len() >= MIN_MACROCYCLE_SIZE && rings[bj].len() >= MIN_MACROCYCLE_SIZE {
                continue;
            }
            let shared_bonds: Vec<usize> = mol
                .bonds
                .iter()
                .enumerate()
                .filter(|(_, b)| {
                    rings[bi].contains(&b.atom1)
                        && rings[bi].contains(&b.atom2)
                        && rings[bj].contains(&b.atom1)
                        && rings[bj].contains(&b.atom2)
                })
                .map(|(idx, _)| idx)
                .collect();
            if shared_bonds.len() > 1 {
                for &bid in &shared_bonds {
                    excluded_bonds[bid] = true;
                }
            }
        }
    }

    for (bond_idx, bond) in mol.bonds.iter().enumerate() {
        let a2 = bond.atom1;
        let a3 = bond.atom2;

        let num_ring_bonds = rings
            .iter()
            // count only SMALL rings toward the fusion limit: a bond shared by
            // many small (rigid) rings is highly-fused and should be skipped, but
            // macrocycle rings (>=9) are flexible and must still receive torsion
            // prefs (otherwise macrocycle amide C-N bonds get skipped and twist).
            .filter(|r| r.contains(&a2) && r.contains(&a3) && r.len() < MIN_MACROCYCLE_SIZE)
            .count();
        if num_ring_bonds > 3 || excluded_bonds[bond_idx] {
            continue;
        }

        for &a1 in &mol.adjacency[a2] {
            if a1 == a3 {
                continue;
            }
            for &a4 in &mol.adjacency[a3] {
                if a4 == a2 || a4 == a1 {
                    continue;
                }
                if let Some((signs, v)) =
                    match_torsion_pattern(mol, a1, a2, a3, a4, &rings, et_version)
                {
                    if !done_bonds[bond_idx] {
                        prefs.push(TorsionPreference {
                            i: a1,
                            j: a2,
                            k: a3,
                            l: a4,
                            signs,
                            v,
                        });
                        done_bonds[bond_idx] = true;
                        break;
                    }
                }
            }
            if done_bonds[bond_idx] {
                break;
            }
        }
    }

    // Basic knowledge: flat aromatic / sp2 ring torsions (RDKit logic)
    for ring in &rings {
        let rsize = ring.len();
        if !(4..=6).contains(&rsize) {
            continue;
        }
        for i in 0..rsize {
            let a1 = ring[i];
            let a2 = ring[(i + 1) % rsize];
            let a3 = ring[(i + 2) % rsize];
            let a4 = ring[(i + 3) % rsize];
            let h1 = get_hyb(mol, a1);
            let h2 = get_hyb(mol, a2);
            let h3 = get_hyb(mol, a3);
            let h4 = get_hyb(mol, a4);
            if matches!(h1, Hybridization::Sp2)
                && matches!(h2, Hybridization::Sp2)
                && matches!(h3, Hybridization::Sp2)
                && matches!(h4, Hybridization::Sp2)
            {
                if let Some(bid) = find_bond_index(mol, a2, a3) {
                    if !done_bonds[bid] {
                        prefs.push(TorsionPreference {
                            i: a1,
                            j: a2,
                            k: a3,
                            l: a4,
                            signs: [1, -1, 1, 1, 1, 1],
                            v: [0.0, 100.0, 0.0, 0.0, 0.0, 0.0],
                        });
                        done_bonds[bid] = true;
                    }
                }
            }
        }
    }

    prefs
}

fn torsion_pref_energy(coords: &[[f64; 3]], prefs: &[TorsionPreference]) -> f64 {
    let mut energy = 0.0;
    for p in prefs {
        let phi = dihedral_angle(coords, p.i, p.j, p.k, p.l);
        for k in 0..6 {
            let m = (k + 1) as f64;
            energy += p.v[k] * (1.0 + p.signs[k] as f64 * (m * phi).cos());
        }
    }
    energy
}

// ============================================================================
// ETKDG Force Field
// ============================================================================

fn etkdg_energy(
    coords: &[[f64; 3]],
    bounds: &DistanceBounds,
    chiral_centers: &[ChiralCenter],
    tetrahedral: &[ChiralCenter],
    pc: &PlanarityConstraints,
    torsion_prefs: &[TorsionPreference],
    bonds_12: &[(usize, usize)],
    angles_13: &[(usize, usize, usize)],
) -> f64 {
    let mut energy = 0.0;
    for i in 0..bounds.n_atoms {
        for j in (i + 1)..bounds.n_atoms {
            let lo = bounds.lower[i][j];
            let hi = bounds.upper[i][j];
            if hi >= MAX_UPPER || (hi - lo) > BASIN_THRESH {
                continue;
            }
            let dx = coords[i][0] - coords[j][0];
            let dy = coords[i][1] - coords[j][1];
            let dz = coords[i][2] - coords[j][2];
            let d = (dx * dx + dy * dy + dz * dz).sqrt().max(1e-10);
            let lo_viol = (lo - d).max(0.0);
            let hi_viol = (d - hi).max(0.0);
            energy += LONG_RANGE_FORCE * (lo_viol * lo_viol + hi_viol * hi_viol);
        }
    }
    // 1-2 constraints: lock bond lengths to bounds-matrix target
    const K_12: f64 = 100.0;
    for &(i, j) in bonds_12 {
        let lo = bounds.lower[i.min(j)][i.max(j)];
        let hi = bounds.upper[i.min(j)][i.max(j)];
        if hi >= MAX_UPPER {
            continue;
        }
        let dx = coords[i][0] - coords[j][0];
        let dy = coords[i][1] - coords[j][1];
        let dz = coords[i][2] - coords[j][2];
        let d = (dx * dx + dy * dy + dz * dz).sqrt().max(1e-10);
        let lo_viol = (lo - d).max(0.0);
        let hi_viol = (d - hi).max(0.0);
        energy += K_12 * (lo_viol * lo_viol + hi_viol * hi_viol);
    }
    // 1-3 constraints: lock angle-derived distances
    for &(i, _, k) in angles_13 {
        let lo = bounds.lower[i.min(k)][i.max(k)];
        let hi = bounds.upper[i.min(k)][i.max(k)];
        if hi >= MAX_UPPER {
            continue;
        }
        let dx = coords[i][0] - coords[k][0];
        let dy = coords[i][1] - coords[k][1];
        let dz = coords[i][2] - coords[k][2];
        let d = (dx * dx + dy * dy + dz * dz).sqrt().max(1e-10);
        let lo_viol = (lo - d).max(0.0);
        let hi_viol = (d - hi).max(0.0);
        energy += K_12 * (lo_viol * lo_viol + hi_viol * hi_viol);
    }
    for cc in chiral_centers {
        let vol = chiral_volume(
            coords,
            cc.central,
            cc.neighbors[0],
            cc.neighbors[1],
            cc.neighbors[2],
        );
        if vol < cc.vol_lower {
            energy += (vol - cc.vol_lower) * (vol - cc.vol_lower);
        } else if vol > cc.vol_upper {
            energy += (vol - cc.vol_upper) * (vol - cc.vol_upper);
        }
    }
    for tc in tetrahedral {
        if !volume_test(coords, tc) {
            energy += 10.0;
        }
    }
    energy += planarity_energy(coords, pc);
    energy += torsion_pref_energy(coords, torsion_prefs);
    energy
}

fn dihedral_gradient_contrib(
    coords: &[[f64; 3]],
    i: usize,
    j: usize,
    k: usize,
    l: usize,
    dedphi: f64,
) -> ([f64; 3], [f64; 3], [f64; 3], [f64; 3]) {
    // Central-difference gradient of dihedral_angle. The previous closed-form had
    // a sign/convention mismatch vs dihedral_angle (see dbg_dihedral_grad_fd);
    // numerical is correct and cheap (dihedral is O(1), only the 4 atoms involved).
    let mut q = [coords[i], coords[j], coords[k], coords[l]];
    let phi_of = |p: &[[f64; 3]; 4]| dihedral_angle4(p[0], p[1], p[2], p[3]);
    let eps = 1e-6;
    let two_pi = 2.0 * std::f64::consts::PI;
    let mut g = [[0.0f64; 3]; 4];
    for a in 0..4 {
        for d in 0..3 {
            let orig = q[a][d];
            q[a][d] = orig + eps;
            let pp = phi_of(&q);
            q[a][d] = orig - eps;
            let pm = phi_of(&q);
            q[a][d] = orig;
            // unwrap atan2 branch cut near +/-pi
            let mut diff = pp - pm;
            if diff > std::f64::consts::PI {
                diff -= two_pi;
            } else if diff < -std::f64::consts::PI {
                diff += two_pi;
            }
            g[a][d] = dedphi * diff / (2.0 * eps);
        }
    }
    (g[0], g[1], g[2], g[3])
}

fn etkdg_gradient(
    coords: &[[f64; 3]],
    bounds: &DistanceBounds,
    chiral_centers: &[ChiralCenter],
    _tetrahedral: &[ChiralCenter],
    pc: &PlanarityConstraints,
    torsion_prefs: &[TorsionPreference],
    _bonds_12: &[(usize, usize)],
    _angles_13: &[(usize, usize, usize)],
) -> Vec<[f64; 3]> {
    let n = coords.len();
    let mut grad = vec![[0.0f64; 3]; n];

    // Distance bounds gradient (with basin threshold + long-range force)
    for i in 0..bounds.n_atoms {
        for j in (i + 1)..bounds.n_atoms {
            let lo = bounds.lower[i][j];
            let hi = bounds.upper[i][j];
            if hi >= MAX_UPPER || (hi - lo) > BASIN_THRESH {
                continue;
            }
            let dx = coords[i][0] - coords[j][0];
            let dy = coords[i][1] - coords[j][1];
            let dz = coords[i][2] - coords[j][2];
            let d = (dx * dx + dy * dy + dz * dz).sqrt().max(1e-10);
            let mut dedd = 0.0;
            if d < lo {
                dedd += -2.0 * LONG_RANGE_FORCE * (lo - d);
            }
            if d > hi {
                dedd += 2.0 * LONG_RANGE_FORCE * (d - hi);
            }
            let gx = dedd * dx / d;
            let gy = dedd * dy / d;
            let gz = dedd * dz / d;
            grad[i][0] += gx;
            grad[i][1] += gy;
            grad[i][2] += gz;
            grad[j][0] -= gx;
            grad[j][1] -= gy;
            grad[j][2] -= gz;
        }
    }

    // Chiral volume gradient
    for cc in chiral_centers {
        let c = cc.central;
        let n1 = cc.neighbors[0];
        let n2 = cc.neighbors[1];
        let n3 = cc.neighbors[2];
        let v1 = [
            coords[n1][0] - coords[c][0],
            coords[n1][1] - coords[c][1],
            coords[n1][2] - coords[c][2],
        ];
        let v2 = [
            coords[n2][0] - coords[c][0],
            coords[n2][1] - coords[c][1],
            coords[n2][2] - coords[c][2],
        ];
        let v3 = [
            coords[n3][0] - coords[c][0],
            coords[n3][1] - coords[c][1],
            coords[n3][2] - coords[c][2],
        ];
        let vol = v1[0] * (v2[1] * v3[2] - v2[2] * v3[1])
            + v1[1] * (v2[2] * v3[0] - v2[0] * v3[2])
            + v1[2] * (v2[0] * v3[1] - v2[1] * v3[0]);
        let mut dvol = 0.0;
        if vol < cc.vol_lower {
            dvol = 2.0 * (vol - cc.vol_lower);
        } else if vol > cc.vol_upper {
            dvol = 2.0 * (vol - cc.vol_upper);
        }
        if dvol != 0.0 {
            let g_n1 = cross_product(v2, v3);
            let g_n2 = cross_product(v3, v1);
            let g_n3 = cross_product(v1, v2);
            let g_c = [
                -(g_n1[0] + g_n2[0] + g_n3[0]),
                -(g_n1[1] + g_n2[1] + g_n3[1]),
                -(g_n1[2] + g_n2[2] + g_n3[2]),
            ];
            for dim in 0..3 {
                grad[c][dim] += dvol * g_c[dim];
                grad[n1][dim] += dvol * g_n1[dim];
                grad[n2][dim] += dvol * g_n2[dim];
                grad[n3][dim] += dvol * g_n3[dim];
            }
        }
    }

    // Planarity gradient: ring torsions
    const K_RING_TOR: f64 = 10.0;
    for &(i, j, k, l) in &pc.ring_torsions {
        let phi = dihedral_angle(coords, i, j, k, l);
        let dedphi = 2.0 * K_RING_TOR * (2.0 * phi).sin();
        let (g1, g2, g3, g4) = dihedral_gradient_contrib(coords, i, j, k, l, dedphi);
        for dim in 0..3 {
            grad[i][dim] += g1[dim];
            grad[j][dim] += g2[dim];
            grad[k][dim] += g3[dim];
            grad[l][dim] += g4[dim];
        }
    }

    // Planarity gradient: exocyclic torsions
    const K_EXO_TOR: f64 = 2.0;
    for &(i, j, k, l) in &pc.exocyclic_torsions {
        let phi = dihedral_angle(coords, i, j, k, l);
        let phi_wrapped = if phi > std::f64::consts::FRAC_PI_2 {
            phi - std::f64::consts::PI
        } else if phi < -std::f64::consts::FRAC_PI_2 {
            phi + std::f64::consts::PI
        } else {
            phi
        };
        let dedphi = 2.0 * K_EXO_TOR * phi_wrapped;
        let (g1, g2, g3, g4) = dihedral_gradient_contrib(coords, i, j, k, l, dedphi);
        for dim in 0..3 {
            grad[i][dim] += g1[dim];
            grad[j][dim] += g2[dim];
            grad[k][dim] += g3[dim];
            grad[l][dim] += g4[dim];
        }
    }

    // Planarity gradient: impropers (small-angle approx)
    for &(central, n1, n2, n3, k_imp) in &pc.impropers {
        let v1 = [
            coords[n1][0] - coords[central][0],
            coords[n1][1] - coords[central][1],
            coords[n1][2] - coords[central][2],
        ];
        let v2 = [
            coords[n2][0] - coords[central][0],
            coords[n2][1] - coords[central][1],
            coords[n2][2] - coords[central][2],
        ];
        let v3 = [
            coords[n3][0] - coords[central][0],
            coords[n3][1] - coords[central][1],
            coords[n3][2] - coords[central][2],
        ];
        let a = [v1[0] - v2[0], v1[1] - v2[1], v1[2] - v2[2]];
        let b = [v1[0] - v3[0], v1[1] - v3[1], v1[2] - v3[2]];
        let normal = cross_product(a, b);
        let n_norm = (normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]).sqrt();
        let v1_norm = (v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2]).sqrt();
        if n_norm < 1e-10 || v1_norm < 1e-10 {
            continue;
        }
        let sin_oop =
            (normal[0] * v1[0] + normal[1] * v1[1] + normal[2] * v1[2]) / (n_norm * v1_norm);
        let chi = sin_oop.clamp(-1.0, 1.0).asin();
        let dchi = 2.0 * k_imp * chi / (1.0 - sin_oop * sin_oop).max(1e-10).sqrt();
        // Gradient of sin_oop w.r.t. v1, v2, v3 (central atom affects all three)
        // For simplicity and robustness, use central-difference for improper gradients
        // (typically very few impropers, so cost is negligible)
        let eps = 1e-7;
        for atom in [central, n1, n2, n3] {
            for dim in 0..3 {
                let mut cp = coords.to_vec();
                cp[atom][dim] += eps;
                let ep = out_of_plane_angle(&cp, central, n1, n2, n3);
                let em = out_of_plane_angle(coords, central, n1, n2, n3);
                grad[atom][dim] += k_imp * (ep * ep - em * em) / eps;
            }
        }
    }

    // Torsion preferences gradient
    for p in torsion_prefs {
        let phi = dihedral_angle(coords, p.i, p.j, p.k, p.l);
        let mut dedphi = 0.0;
        for k in 0..6 {
            let m = (k + 1) as f64;
            dedphi += -p.v[k] * p.signs[k] as f64 * m * (m * phi).sin();
        }
        let (g1, g2, g3, g4) = dihedral_gradient_contrib(coords, p.i, p.j, p.k, p.l, dedphi);
        for dim in 0..3 {
            grad[p.i][dim] += g1[dim];
            grad[p.j][dim] += g2[dim];
            grad[p.k][dim] += g3[dim];
            grad[p.l][dim] += g4[dim];
        }
    }

    grad
}

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
    coord_map: &std::collections::HashMap<usize, [f64; 3]>,
) -> f64 {
    let n = coords.len();
    if n == 0 {
        return 0.0;
    }
    // L-BFGS (two-loop recursion + backtracking line search) over etkdg_energy /
    // etkdg_gradient. The prior fixed-step normalized descent (`step = 0.1/max_g`)
    // diverged as the gradient shrank (step -> inf near a minimum) and never
    // converged — stalling in bad local minima even for tiny molecules (butadiene
    // E~282 vs min ~6; 2000 iters didn't help). L-BFGS uses curvature history to
    // take well-scaled steps and actually reach the force tolerance.
    let dim = 3 * n;
    let mut x = vec![0.0f64; dim];
    for i in 0..n {
        x[3 * i] = coords[i][0];
        x[3 * i + 1] = coords[i][1];
        x[3 * i + 2] = coords[i][2];
    }
    let fixed: Vec<bool> = (0..n).map(|i| coord_map.contains_key(&i)).collect();
    let energy_at = |xx: &[f64]| -> f64 {
        let c: Vec<[f64; 3]> = (0..n).map(|i| [xx[3 * i], xx[3 * i + 1], xx[3 * i + 2]]).collect();
        etkdg_energy(&c, bounds, chiral_centers, tetrahedral, pc, torsion_prefs, bonds_12, angles_13)
    };
    let gradient_at = |xx: &[f64]| -> Vec<f64> {
        let c: Vec<[f64; 3]> = (0..n).map(|i| [xx[3 * i], xx[3 * i + 1], xx[3 * i + 2]]).collect();
        let g3 = etkdg_gradient(&c, bounds, chiral_centers, tetrahedral, pc, torsion_prefs, bonds_12, angles_13);
        let mut gx = vec![0.0f64; dim];
        for i in 0..n {
            if fixed[i] {
                continue;
            }
            gx[3 * i] = g3[i][0];
            gx[3 * i + 1] = g3[i][1];
            gx[3 * i + 2] = g3[i][2];
        }
        gx
    };

    let mut f = energy_at(&x);
    let mut g = gradient_at(&x);
    const M: usize = 8;
    let mut s_hist: Vec<Vec<f64>> = Vec::with_capacity(M);
    let mut y_hist: Vec<Vec<f64>> = Vec::with_capacity(M);
    let mut rho_hist: Vec<f64> = Vec::with_capacity(M);
    let mut iters_done = 0usize;
    let mut converged = false;
    let mut stall = 0usize;
    for _iter in 0..max_iter {
        iters_done += 1;
        // convergence: max per-atom force (excluding fixed atoms)
        let mut max_g = 0.0f64;
        for i in 0..n {
            if fixed[i] {
                continue;
            }
            let gm = (g[3 * i] * g[3 * i] + g[3 * i + 1] * g[3 * i + 1] + g[3 * i + 2] * g[3 * i + 2]).sqrt();
            if gm > max_g {
                max_g = gm;
            }
        }
        if max_g < force_tol {
            converged = true;
            break;
        }

        // L-BFGS two-loop recursion -> q = H_k·g ; direction = -q
        let k = s_hist.len();
        let mut q = g.clone();
        let mut alpha = vec![0.0f64; k];
        for j in (0..k).rev() {
            let sjq: f64 = s_hist[j].iter().zip(q.iter()).map(|(a, b)| a * b).sum();
            alpha[j] = rho_hist[j] * sjq;
            let yj = &y_hist[j];
            for d in 0..dim {
                q[d] -= alpha[j] * yj[d];
            }
        }
        let gamma = if k > 0 {
            let sy: f64 = s_hist[k - 1].iter().zip(y_hist[k - 1].iter()).map(|(a, b)| a * b).sum();
            let yy: f64 = y_hist[k - 1].iter().map(|v| v * v).sum();
            if yy > 1e-20 { sy / yy } else { 1.0 }
        } else {
            1.0
        };
        for d in 0..dim {
            q[d] *= gamma;
        }
        for j in 0..k {
            let yjq: f64 = y_hist[j].iter().zip(q.iter()).map(|(a, b)| a * b).sum();
            let beta = rho_hist[j] * yjq;
            let sj = &s_hist[j];
            for d in 0..dim {
                q[d] += sj[d] * (alpha[j] - beta);
            }
        }
        let mut dir: Vec<f64> = q.iter().map(|qi| -qi).collect();
        let mut dg: f64 = dir.iter().zip(g.iter()).map(|(a, b)| a * b).sum();
        // If not a descent direction (rare, bad curvature), reset and use steepest descent.
        if dg >= 0.0 {
            s_hist.clear();
            y_hist.clear();
            rho_hist.clear();
            dir = g.iter().map(|gi| -gi).collect();
            dg = dir.iter().zip(g.iter()).map(|(a, b)| a * b).sum();
        }

        // Backtracking Armijo line search: f(x + step·dir) <= f + c·step·dg
        let c_armijo = 1e-4;
        let mut step = if k == 0 { 1.0 / max_g.max(1e-10) } else { 1.0 };
        let mut accepted = false;
        for _ls in 0..25 {
            let x_new: Vec<f64> = (0..dim).map(|d| x[d] + step * dir[d]).collect();
            let f_new = energy_at(&x_new);
            if f_new.is_finite() && f_new <= f + c_armijo * step * dg {
                let g_new = gradient_at(&x_new);
                let s: Vec<f64> = (0..dim).map(|d| x_new[d] - x[d]).collect();
                let y: Vec<f64> = (0..dim).map(|d| g_new[d] - g[d]).collect();
                let sy: f64 = s.iter().zip(y.iter()).map(|(a, b)| a * b).sum();
                if sy > 1e-20 {
                    if s_hist.len() >= M {
                        s_hist.remove(0);
                        y_hist.remove(0);
                        rho_hist.remove(0);
                    }
                    s_hist.push(s);
                    y_hist.push(y);
                    rho_hist.push(1.0 / sy);
                }
                x = x_new;
                g = g_new;
                f = f_new;
                accepted = true;
                break;
            }
            step *= 0.5;
        }
        if !accepted {
            // line search failed: reset L-BFGS history and keep going with fresh
            // steepest-descent steps; only give up after several consecutive
            // failures. (Previously broke on the first failure, leaving strained-
            // ring molecules under-relaxed after ~2 iters.)
            stall += 1;
            s_hist.clear();
            y_hist.clear();
            rho_hist.clear();
            if stall >= 5 {
                break;
            }
        } else {
            stall = 0;
        }
    }

    eprintln!("MIN n={n} iters={iters_done}/{max_iter} converged={converged} E={f:.1}");
    for i in 0..n {
        coords[i][0] = x[3 * i];
        coords[i][1] = x[3 * i + 1];
        coords[i][2] = x[3 * i + 2];
    }
    f
}

fn bounds_fulfilled(coords: &[[f64; 3]], bounds: &DistanceBounds, atom_indices: &[usize]) -> bool {
    for i in 0..atom_indices.len() {
        for j in (i + 1)..atom_indices.len() {
            let a1 = atom_indices[i];
            let a2 = atom_indices[j];
            let dx = coords[a1][0] - coords[a2][0];
            let dy = coords[a1][1] - coords[a2][1];
            let dz = coords[a1][2] - coords[a2][2];
            let d = (dx * dx + dy * dy + dz * dz).sqrt();
            let lb = bounds.get_lower(a1, a2);
            let ub = bounds.get_upper(a1, a2);
            if ((d < lb) && (d - lb).abs() > 0.1 * ub) || ((d > ub) && (d - ub).abs() > 0.1 * ub) {
                return false;
            }
        }
    }
    true
}

fn collect_bonds(mol: &Molecule) -> Vec<(usize, usize)> {
    mol.bonds.iter().map(|b| (b.atom1, b.atom2)).collect()
}

fn collect_angles(mol: &Molecule) -> Vec<(usize, usize, usize)> {
    let mut angles = Vec::new();
    let n_bonds = mol.bonds.len();
    for i in 0..n_bonds {
        for j in (i + 1)..n_bonds {
            let b1 = &mol.bonds[i];
            let b2 = &mol.bonds[j];
            let (a11, a12, a21, a22) = (b1.atom1, b1.atom2, b2.atom1, b2.atom2);
            let mut found = None;
            if a12 == a21 {
                found = Some((a11, a12, a22));
            } else if a12 == a22 {
                found = Some((a11, a12, a21));
            } else if a11 == a21 {
                found = Some((a12, a11, a22));
            } else if a11 == a22 {
                found = Some((a12, a11, a21));
            }
            if let Some(angle) = found {
                angles.push(angle);
            }
        }
    }
    angles
}

// ============================================================================
// Fragment Embedding Helpers
// ============================================================================

fn find_connected_components(mol: &Molecule) -> Vec<Vec<usize>> {
    let n = mol.atoms.len();
    let mut visited = vec![false; n];
    let mut components = Vec::new();

    for start in 0..n {
        if visited[start] {
            continue;
        }
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

    let n_atoms = new_atoms.len();
    let adjacency = crate::molecule::graph::build_adjacency_list_from_bonds(&n_atoms, &new_bonds);
    Molecule {
        atoms: new_atoms,
        bonds: new_bonds,
        name: mol.name.clone(),
        adjacency,
    }
}

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

// ============================================================================
// Main Public API
// ============================================================================

/// Fix aniline-like NH2 geometries after ETKDG minimization.
/// RDKit ETKDG produces slightly pyramidal aniline N (H-N-H ~117.5°,
/// H atoms ±0.84 Å out of the ring plane) rather than fully planar.
#[allow(dead_code)]
fn fix_aniline_nh2_geometry(
    coords: &mut [[f64; 3]],
    mol: &Molecule,
    pc: &PlanarityConstraints,
) {
    let aromatic_atoms = &pc.aromatic_atoms;
    for atom_idx in 0..mol.atoms.len() {
        // Only consider non-aromatic sp2 nitrogen atoms
        if aromatic_atoms.contains(&atom_idx) {
            continue;
        }
        if mol.atoms[atom_idx].symbol != "N" {
            continue;
        }
        let hyb = crate::molecule::graph::determine_hybridization(atom_idx, mol);
        if !matches!(hyb, Hybridization::Sp2) {
            continue;
        }
        // Must be directly bonded to at least one aromatic atom
        let has_aromatic_neighbor = mol.adjacency[atom_idx]
            .iter()
            .any(|n| aromatic_atoms.contains(n));
        if !has_aromatic_neighbor {
            continue;
        }
        // Must have exactly 2 H neighbors and 1 heavy neighbor
        let h_neighbors: Vec<usize> = mol.adjacency[atom_idx]
            .iter()
            .filter(|&&n| mol.atoms[n].symbol == "H")
            .copied()
            .collect();
        let heavy_neighbors: Vec<usize> = mol.adjacency[atom_idx]
            .iter()
            .filter(|&&n| mol.atoms[n].symbol != "H")
            .copied()
            .collect();
        if h_neighbors.len() != 2 || heavy_neighbors.len() != 1 {
            continue;
        }
        let heavy = heavy_neighbors[0];
        let h1 = h_neighbors[0];
        let h2 = h_neighbors[1];

        // Compute ring-plane normal at the N position.
        // Use the aromatic atoms that are neighbors of the heavy atom
        // (the ring carbon attached to N) to define the local plane.
        let ring_nbrs: Vec<usize> = mol.adjacency[heavy]
            .iter()
            .filter(|&&n| aromatic_atoms.contains(&n) && n != atom_idx)
            .copied()
            .collect();
        if ring_nbrs.len() < 2 {
            continue;
        }
        let a = coords[ring_nbrs[0]];
        let b = coords[heavy];
        let c = coords[ring_nbrs[1]];
        let u = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
        let v = [c[0] - b[0], c[1] - b[1], c[2] - b[2]];
        let mut normal = [
            u[1] * v[2] - u[2] * v[1],
            u[2] * v[0] - u[0] * v[2],
            u[0] * v[1] - u[1] * v[0],
        ];
        let n_norm = (normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]).sqrt();
        if n_norm < 1e-10 {
            continue;
        }
        for i in 0..3 {
            normal[i] /= n_norm;
        }

        // Ensure the normal points toward the side where the N atom currently sits
        let n_to_plane = (coords[atom_idx][0] - b[0]) * normal[0]
            + (coords[atom_idx][1] - b[1]) * normal[1]
            + (coords[atom_idx][2] - b[2]) * normal[2];
        if n_to_plane < 0.0 {
            for i in 0..3 {
                normal[i] = -normal[i];
            }
        }

        // Direction from N to heavy atom (in-plane)
        let mut to_heavy = [
            coords[heavy][0] - coords[atom_idx][0],
            coords[heavy][1] - coords[atom_idx][1],
            coords[heavy][2] - coords[atom_idx][2],
        ];
        let th_len = (to_heavy[0] * to_heavy[0] + to_heavy[1] * to_heavy[1] + to_heavy[2] * to_heavy[2]).sqrt();
        if th_len < 1e-10 {
            continue;
        }
        for i in 0..3 {
            to_heavy[i] /= th_len;
        }

        // Perpendicular direction in the ring plane
        let mut perp = [
            normal[1] * to_heavy[2] - normal[2] * to_heavy[1],
            normal[2] * to_heavy[0] - normal[0] * to_heavy[2],
            normal[0] * to_heavy[1] - normal[1] * to_heavy[0],
        ];
        let p_len = (perp[0] * perp[0] + perp[1] * perp[1] + perp[2] * perp[2]).sqrt();
        if p_len < 1e-10 {
            continue;
        }
        for i in 0..3 {
            perp[i] /= p_len;
        }

        // RDKit ETKDGv3 aniline NH2: N-H avg 1.049, H-N-H 118.17°
        // With PYRAMIDAL_TILT = 90°, formula simplifies to exact H-N-H = target
        const N_H_BOND: f64 = 1.049;
        const H_N_H_ANGLE: f64 = 118.0_f64.to_radians();
        const PYRAMIDAL_TILT: f64 = 90.0_f64.to_radians();

        let half_angle = H_N_H_ANGLE / 2.0;
        let cos_h = half_angle.cos();
        let sin_h = half_angle.sin();
        let cos_t = PYRAMIDAL_TILT.cos();
        let sin_t = PYRAMIDAL_TILT.sin();

        // H1 (above plane)
        coords[h1][0] = coords[atom_idx][0]
            + N_H_BOND * (cos_h * to_heavy[0] + sin_h * (cos_t * perp[0] + sin_t * normal[0]));
        coords[h1][1] = coords[atom_idx][1]
            + N_H_BOND * (cos_h * to_heavy[1] + sin_h * (cos_t * perp[1] + sin_t * normal[1]));
        coords[h1][2] = coords[atom_idx][2]
            + N_H_BOND * (cos_h * to_heavy[2] + sin_h * (cos_t * perp[2] + sin_t * normal[2]));

        // H2 (below plane)
        coords[h2][0] = coords[atom_idx][0]
            + N_H_BOND * (cos_h * to_heavy[0] + sin_h * (cos_t * perp[0] - sin_t * normal[0]));
        coords[h2][1] = coords[atom_idx][1]
            + N_H_BOND * (cos_h * to_heavy[1] + sin_h * (cos_t * perp[1] - sin_t * normal[1]));
        coords[h2][2] = coords[atom_idx][2]
            + N_H_BOND * (cos_h * to_heavy[2] + sin_h * (cos_t * perp[2] - sin_t * normal[2]));
    }
}

// ============================================================================

/// Rotate the k-side fragment of the j-k bond around the j->k axis by `angle`
/// (Rodrigues). The fragment = atoms reachable from k without crossing j.
/// Used by the torsion-snap basin-hop to cross torsional barriers.
fn rotate_fragment_around_bond(coords: &mut [[f64; 3]], mol: &Molecule, j: usize, k: usize, angle: f64) {
    let n = mol.atoms.len();
    let mut visited = vec![false; n];
    let mut frag: Vec<usize> = Vec::new();
    let mut stack = vec![k];
    visited[j] = true; // barrier — don't cross j
    visited[k] = true;
    while let Some(a) = stack.pop() {
        frag.push(a);
        for &nb in &mol.adjacency[a] {
            if !visited[nb] {
                visited[nb] = true;
                stack.push(nb);
            }
        }
    }
    let p1 = coords[j];
    let axis = [coords[k][0] - p1[0], coords[k][1] - p1[1], coords[k][2] - p1[2]];
    let axlen = (axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2]).sqrt();
    if axlen < 1e-10 {
        return;
    }
    let u = [axis[0] / axlen, axis[1] / axlen, axis[2] / axlen];
    let ca = angle.cos();
    let sa = angle.sin();
    let omca = 1.0 - ca;
    for &a in &frag {
        let v = [coords[a][0] - p1[0], coords[a][1] - p1[1], coords[a][2] - p1[2]];
        let nxv = [
            u[1] * v[2] - u[2] * v[1],
            u[2] * v[0] - u[0] * v[2],
            u[0] * v[1] - u[1] * v[0],
        ];
        let ndv = u[0] * v[0] + u[1] * v[1] + u[2] * v[2];
        coords[a][0] = p1[0] + v[0] * ca + nxv[0] * sa + u[0] * ndv * omca;
        coords[a][1] = p1[1] + v[1] * ca + nxv[1] * sa + u[1] * ndv * omca;
        coords[a][2] = p1[2] + v[2] * ca + nxv[2] * sa + u[2] * ndv * omca;
    }
}

/// 3-sphere intersection: given centers p1,p2,p3 and radii r1,r2,r3, return
/// the two intersection points (or None if degenerate).
fn trilaterate(p1: [f64; 3], p2: [f64; 3], p3: [f64; 3], r1: f64, r2: f64, r3: f64) -> Option<([f64; 3], [f64; 3])> {
    let d = ((p2[0] - p1[0]).powi(2) + (p2[1] - p1[1]).powi(2) + (p2[2] - p1[2]).powi(2)).sqrt();
    if d < 1e-10 {
        return None;
    }
    let ex = [(p2[0] - p1[0]) / d, (p2[1] - p1[1]) / d, (p2[2] - p1[2]) / d];
    let p3p1 = [p3[0] - p1[0], p3[1] - p1[1], p3[2] - p1[2]];
    let i_val = ex[0] * p3p1[0] + ex[1] * p3p1[1] + ex[2] * p3p1[2];
    let p3_perp = [p3p1[0] - i_val * ex[0], p3p1[1] - i_val * ex[1], p3p1[2] - i_val * ex[2]];
    let ppl = (p3_perp[0] * p3_perp[0] + p3_perp[1] * p3_perp[1] + p3_perp[2] * p3_perp[2]).sqrt();
    if ppl < 1e-10 {
        return None;
    }
    let ey = [p3_perp[0] / ppl, p3_perp[1] / ppl, p3_perp[2] / ppl];
    let ez = cross_product(ex, ey);
    let x = (d * d + r1 * r1 - r2 * r2) / (2.0 * d);
    let j_val = ey[0] * p3p1[0] + ey[1] * p3p1[1] + ey[2] * p3p1[2];
    if j_val.abs() < 1e-10 {
        return None;
    }
    let y = (i_val * i_val + j_val * j_val + r1 * r1 - r3 * r3 - 2.0 * i_val * x) / (2.0 * j_val);
    let z2 = r1 * r1 - x * x - y * y;
    if z2 < 0.0 {
        return None;
    }
    let z = z2.sqrt();
    let base = [p1[0] + x * ex[0] + y * ey[0], p1[1] + x * ex[1] + y * ey[1], p1[2] + x * ex[2] + y * ey[2]];
    Some((
        [base[0] + z * ez[0], base[1] + z * ez[1], base[2] + z * ez[2]],
        [base[0] - z * ez[0], base[1] - z * ez[1], base[2] - z * ez[2]],
    ))
}

/// Reposition misplaced H atoms via 3-sphere trilateration from the ETKDG
/// distance bounds (1-2 bond + 1-3 angle constraints). For each H with ≥2 heavy
/// neighbors on its bonded atom, analytically compute the correct position from
/// the bounds — gives the minimizer a good starting point it can't reach from
/// the 4D projection.
/// Snap badly-violated 1-2 bonds (off > 0.15 Å) by translating the smaller
/// fragment along the bond axis to the target distance. Helps the minimizer
/// escape ring-closure failures (e.g. cyclopentene C4-C0 at 1.80 instead of 1.54).
/// Guarded by accept-only-if-ETKDG-energy-drops at the call site.
fn snap_bond_lengths(coords: &mut [[f64; 3]], mol: &Molecule, bounds: &DistanceBounds) {
    for b in &mol.bonds {
        let (i, j) = (b.atom1, b.atom2);
        let (lo, hi) = (i.min(j), i.max(j));
        let target = (bounds.lower[lo][hi] + bounds.upper[lo][hi]) / 2.0;
        let cur = atom_dist(coords, i, j);
        if (cur - target).abs() < 0.12 || cur < 1e-10 {
            continue;
        }
        // BFS fragment on each side; snap the smaller one
        let frag_j = bfs_fragment(mol, j, i);
        let frag_i = bfs_fragment(mol, i, j);
        let (frag, anchor, mover) = if frag_j.len() <= frag_i.len() {
            (frag_j, i, j)
        } else {
            (frag_i, j, i)
        };
        let axis = [
            coords[mover][0] - coords[anchor][0],
            coords[mover][1] - coords[anchor][1],
            coords[mover][2] - coords[anchor][2],
        ];
        let unit = [axis[0] / cur, axis[1] / cur, axis[2] / cur];
        let shift = target - cur; // negative = shorten (move mover-side toward anchor)
        for &a in &frag {
            coords[a][0] += shift * unit[0];
            coords[a][1] += shift * unit[1];
            coords[a][2] += shift * unit[2];
        }
    }
}

fn bfs_fragment(mol: &Molecule, start: usize, exclude: usize) -> Vec<usize> {
    let mut visited = vec![false; mol.atoms.len()];
    visited[exclude] = true;
    visited[start] = true;
    let mut frag = vec![start];
    let mut stack = vec![start];
    while let Some(a) = stack.pop() {
        for &nb in &mol.adjacency[a] {
            if !visited[nb] {
                visited[nb] = true;
                frag.push(nb);
                stack.push(nb);
            }
        }
    }
    frag
}

fn atom_dist(coords: &[[f64; 3]], a: usize, b: usize) -> f64 {
    let d = [coords[a][0] - coords[b][0], coords[a][1] - coords[b][1], coords[a][2] - coords[b][2]];
    (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt()
}

fn dist_points(a: [f64; 3], b: [f64; 3]) -> f64 {
    let d = [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
    (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt()
}

fn trilaterate_hydrogens(coords: &mut [[f64; 3]], mol: &Molecule, bounds: &DistanceBounds) {
    let mut processed_a = vec![false; mol.atoms.len()];
    for h_idx in 0..mol.atoms.len() {
        if mol.atoms[h_idx].symbol != "H" || processed_a[h_idx] {
            continue;
        }
        let a_idx = match mol.adjacency[h_idx].iter().find(|&&n| mol.atoms[n].symbol != "H") {
            Some(&n) => n,
            None => continue,
        };
        if processed_a[a_idx] {
            continue;
        }
        processed_a[a_idx] = true;
        let r0_ha = {
            let (lo, hi) = (h_idx.min(a_idx), h_idx.max(a_idx));
            (bounds.lower[lo][hi] + bounds.upper[lo][hi]) / 2.0
        };
        let h_atoms: Vec<usize> = mol.adjacency[a_idx]
            .iter()
            .filter(|&&n| mol.atoms[n].symbol == "H")
            .copied()
            .collect();
        if !h_atoms.iter().any(|&h| (atom_dist(coords, a_idx, h) - r0_ha).abs() > 0.05) {
            continue;
        }
        let heavy_nbrs: Vec<usize> = mol.adjacency[a_idx]
            .iter()
            .filter(|&&n| n != h_idx && mol.atoms[n].symbol != "H")
            .copied()
            .collect();
        if heavy_nbrs.len() < 2 {
            continue;
        }
        let n1 = heavy_nbrs[0];
        let n2 = heavy_nbrs[1];
        let mid_13 = |n: usize| -> f64 {
            let (lo, hi) = (h_idx.min(n), h_idx.max(n));
            let upper = bounds.upper[lo][hi];
            if !(0.1..MAX_UPPER).contains(&upper) {
                let r_an = atom_dist(coords, a_idx, n);
                let cos_a = -1.0_f64 / 3.0;
                (r0_ha * r0_ha + r_an * r_an - 2.0 * r0_ha * r_an * cos_a).sqrt()
            } else {
                (bounds.lower[lo][hi] + upper) / 2.0
            }
        };
        let r0_hn1 = mid_13(n1);
        let r0_hn2 = mid_13(n2);
        if let Some((sol1, sol2)) =
            trilaterate(coords[a_idx], coords[n1], coords[n2], r0_ha, r0_hn1, r0_hn2)
        {
            match h_atoms.len() {
                1 => {
                    let d1 = dist_points(sol1, coords[h_atoms[0]]);
                    let d2 = dist_points(sol2, coords[h_atoms[0]]);
                    coords[h_atoms[0]] = if d1 < d2 { sol1 } else { sol2 };
                }
                _ => {
                    let d1_0 = dist_points(sol1, coords[h_atoms[0]]);
                    let d1_1 = dist_points(sol1, coords[h_atoms[1]]);
                    if d1_0 <= d1_1 {
                        coords[h_atoms[0]] = sol1;
                        coords[h_atoms[1]] = sol2;
                    } else {
                        coords[h_atoms[0]] = sol2;
                        coords[h_atoms[1]] = sol1;
                    }
                }
            }
        }
    }
}

pub fn generate_initial_coords(mol: &Molecule) -> Vec<[f64; 3]> {
    let config = ETKDGConfig::default();
    generate_initial_coords_with_config(mol, &config)
}

pub fn generate_initial_coords_with_config(mol: &Molecule, config: &ETKDGConfig) -> Vec<[f64; 3]> {
    // M0 scaffolding: dispatch on the rdkit_faithful flag. The default path is
    // unchanged; the _rdkit path is currently a passthrough (the real RDKit-faithful
    // bounds/minimizer/force-field land in M1+, behind this same flag).
    if config.rdkit_faithful {
        generate_initial_coords_rdkit(mol, config)
    } else {
        generate_initial_coords_default(mol, config)
    }
}

fn generate_initial_coords_default(mol: &Molecule, config: &ETKDGConfig) -> Vec<[f64; 3]> {
    EXP_RDKIT_ALL.store(false, std::sync::atomic::Ordering::Relaxed);
    embed_impl(mol, config)
}

fn generate_initial_coords_rdkit(mol: &Molecule, config: &ETKDGConfig) -> Vec<[f64; 3]> {
    EXP_RDKIT_ALL.store(true, std::sync::atomic::Ordering::Relaxed);
    embed_impl(mol, config)
}

fn embed_impl(mol: &Molecule, config: &ETKDGConfig) -> Vec<[f64; 3]> {
    if mol.atoms.is_empty() {
        return Vec::new();
    }

    // Handle multi-fragment molecules
    let components = find_connected_components(mol);
    if components.len() > 1 {
        let mut all_coords = vec![[0.0; 3]; mol.atoms.len()];

        for comp in &components {
            let submol = extract_submol(mol, comp);
            let subcoords = generate_initial_coords_with_config(&submol, config);
            for (i, &atom_idx) in comp.iter().enumerate() {
                if i < subcoords.len() {
                    all_coords[atom_idx] = subcoords[i];
                }
            }
        }

        spread_fragments(&mut all_coords, &components);
        return all_coords;
    }

    let n_atoms = mol.atoms.len();

    let mut rng = if config.random_seed >= 0 {
        Rng::new(config.random_seed as u64)
    } else {
        let seed = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap_or_default()
            .as_nanos() as u64;
        Rng::new(seed)
    };

    let mut bounds = build_distance_bounds(mol, config);
    bounds.smooth_triangle_inequality(config.triangle_smoothing_epsilon);

    let pc = build_planarity_constraints(mol);
    let (chiral_centers, tetrahedral) = find_chiral_centers(mol);
    let (double_bond_ends, stereo_dbs, atropisomers) = find_stereo_bonds(mol);
    let torsion_prefs = build_torsion_preferences(mol, config.et_version);
    let bonds_12 = collect_bonds(mol);
    let angles_13 = collect_angles(mol);

    let mut best_coords = Vec::new();
    let mut best_energy = f64::INFINITY;
    let max_attempts = if config.max_attempts > 0 {
        config.max_attempts
    } else {
        (10 * n_atoms).max(10)
    };

    let start_time = Instant::now();
    let timeout_duration = if config.timeout_ms > 0 {
        Some(std::time::Duration::from_millis(config.timeout_ms))
    } else {
        None
    };

    for _attempt in 0..max_attempts {
        // Check timeout
        if let Some(duration) = timeout_duration {
            if start_time.elapsed() > duration {
                eprintln!("ETKDG timeout after {} attempts", _attempt);
                break;
            }
        }
        let mut coords_4d = generate_initial_coords_from_bounds(&bounds, &mut rng);

        // Apply coord_map constraints
        for (&atom_idx, &fixed_pos) in &config.coord_map {
            if atom_idx < coords_4d.len() {
                coords_4d[atom_idx] = [fixed_pos[0], fixed_pos[1], fixed_pos[2], 0.0];
            }
        }

        // Step 1: First 4D minimization — distance bounds + chirality + light 4th-dim penalty
        let first_e = minimize_4d_first(&mut coords_4d, &bounds, &chiral_centers, 400);
        let e_per_atom_4d = first_e / n_atoms as f64;
        if e_per_atom_4d >= MAX_MINIMIZED_E_PER_ATOM {
            if first_e < best_energy {
                best_energy = first_e;
                best_coords = coords_4d.iter().map(|c| [c[0], c[1], c[2]]).collect();
            }
            continue;
        }

        // Steps 2-3 (tetrahedral + chiral pre-checks) used to be hard gates that
        // `continue`d on failure. For macrocycles the 4D->3D projection frequently
        // flattens a stereocenter (e.g. RDKit #9143's [C@@H]), failing
        // check_tetrahedral on every attempt and aborting before the 3D refinement
        // (minimize_etkdg) — whose chiral-volume term is designed to fix exactly
        // that — could run. Pre-checks are now advisory: proceed to refinement and
        // let the post-refinement acceptance gate below (check_chiral_centers +
        // center_in_volume_tol + planarity + clashes + energy) judge the result.

        // Step 4: Minimize 4th dimension — strong 4th-dim weight, relaxed chirality
        minimize_4d_collapse(&mut coords_4d, &bounds, &chiral_centers, 200);

        // Step 5: Flatten aromatic rings before torsion minimization
        let mut coords_3d: Vec<[f64; 3]> = coords_4d.iter().map(|c| [c[0], c[1], c[2]]).collect();
        if !pc.aromatic_atoms.is_empty() {
            flatten_aromatic_rings(&mut coords_3d, mol, &pc);
        }

        // Re-apply coord_map after projection to 3D
        for (&atom_idx, &fixed_pos) in &config.coord_map {
            if atom_idx < coords_3d.len() {
                coords_3d[atom_idx] = fixed_pos;
            }
        }

        let mut energy = minimize_etkdg(
            &mut coords_3d,
            &bounds,
            &chiral_centers,
            &tetrahedral,
            &pc,
            &torsion_prefs,
            &bonds_12,
            &angles_13,
            300,
            1e-3,
            &config.coord_map,
        );

        // Torsion barrier-crossing: L-BFGS gets stuck at torsional barriers from
        // the 4D-projection start (e.g. butadiene ~100° instead of 180°). Snap the
        // most-twisted torsion prefs toward a planar/trans minimum (rotate the
        // k-side fragment around the central bond) and re-minimize; accept only if
        // total ETKDG energy drops. Bounded to 3 snaps/attempt.
        let mut snaps = 0usize;
        while snaps < 3 {
            let mut best_snap: Option<(usize, f64)> = None;
            for (pi, p) in torsion_prefs.iter().enumerate() {
                let phi = dihedral_angle(&coords_3d, p.i, p.j, p.k, p.l);
                let target = if phi.abs() > std::f64::consts::FRAC_PI_2 {
                    std::f64::consts::PI
                } else {
                    0.0
                };
                let delta = target - phi;
                if delta.abs() > std::f64::consts::FRAC_PI_4
                    && best_snap.is_none_or(|(_, d)| delta.abs() > d.abs())
                {
                    best_snap = Some((pi, delta));
                }
            }
            let Some((pi, delta)) = best_snap else { break };
            let p = &torsion_prefs[pi];
            let mut trial = coords_3d.clone();
            rotate_fragment_around_bond(&mut trial, mol, p.j, p.k, delta);
            let te = minimize_etkdg(
                &mut trial,
                &bounds,
                &chiral_centers,
                &tetrahedral,
                &pc,
                &torsion_prefs,
                &bonds_12,
                &angles_13,
                300,
                1e-3,
                &config.coord_map,
            );
            if te < energy {
                energy = te;
                coords_3d = trial;
                snaps += 1;
            } else {
                break;
            }
        }

        // Bond-snap: fix badly-violated 1-2 bonds (ring-closure failures) by
        // translating the smaller fragment to the target distance. Guarded.
        {
            let e_before = etkdg_energy(
                &coords_3d, &bounds, &chiral_centers, &tetrahedral, &pc,
                &torsion_prefs, &bonds_12, &angles_13,
            );
            let mut trial = coords_3d.clone();
            snap_bond_lengths(&mut trial, mol, &bounds);
            let e_after = etkdg_energy(
                &trial, &bounds, &chiral_centers, &tetrahedral, &pc,
                &torsion_prefs, &bonds_12, &angles_13,
            );
            if e_after < e_before {
                coords_3d = trial;
            }
        }

        // H trilateration: analytically place misplaced H atoms via 3-sphere
        // intersection from their ETKDG 1-2/1-3 distance bounds. Gives the
        // H-only relaxation below a good starting point.
        trilaterate_hydrogens(&mut coords_3d, mol, &bounds);

        // H-only relaxation: fix all non-H atoms, re-minimize to let H atoms
        // settle into their ideal positions (satisfying 1-2/1-3 bounds) without
        // the heavy-atom landscape trapping them. Cheap (few DOF) and can't
        // regress (heavy atoms are fixed).
        {
            let mut h_coord_map = std::collections::HashMap::new();
            for i in 0..mol.atoms.len() {
                if mol.atoms[i].symbol != "H" {
                    h_coord_map.insert(i, coords_3d[i]);
                }
            }
            for (&k, &v) in &config.coord_map {
                h_coord_map.insert(k, v);
            }
            let he = minimize_etkdg(
                &mut coords_3d,
                &bounds,
                &chiral_centers,
                &tetrahedral,
                &pc,
                &torsion_prefs,
                &bonds_12,
                &angles_13,
                300,
                1e-3,
                &h_coord_map,
            );
            if he < energy {
                energy = he;
            }
        }

        // Post-process aniline-like NH2 groups: RDKit ETKDG gives slightly
        // pyramidal geometry (H-N-H ~117.5°, H atoms ±0.84 Å out of ring plane)
        // instead of fully planar.  Re-place the H atoms explicitly.
        // fix_aniline_nh2_geometry was disabled: it placed NH2 H's at wrong C-N-H
        // angles (59° instead of ~112°) — the H-trilateration + H-only relaxation
        // now place them correctly, and the ad-hoc fix was raising aniline's
        // energy from 12 to 127 kcal/mol.

        let planar = check_planarity(&coords_3d, mol, &pc, 0.1);
        let db_geom_ok = double_bond_geometry_checks(&coords_3d, &double_bond_ends);

        let mut chiral_ok = true;
        if !chiral_centers.is_empty() {
            chiral_ok = check_chiral_centers(&coords_3d, &chiral_centers);
            if chiral_ok {
                let mut atoms = std::collections::HashSet::new();
                for cc in &chiral_centers {
                    if cc.central != cc.neighbors[3] {
                        atoms.insert(cc.central);
                        atoms.insert(cc.neighbors[0]);
                        atoms.insert(cc.neighbors[1]);
                        atoms.insert(cc.neighbors[2]);
                        atoms.insert(cc.neighbors[3]);
                    }
                }
                let atoms_vec: Vec<usize> = atoms.into_iter().collect();
                if !atoms_vec.is_empty() && !bounds_fulfilled(&coords_3d, &bounds, &atoms_vec) {
                    chiral_ok = false;
                }
            }
            if chiral_ok {
                for cc in &chiral_centers {
                    if !center_in_volume_tol(&coords_3d, cc, CHIRAL_CENTERINVOLUME_TOL) {
                        chiral_ok = false;
                        break;
                    }
                }
            }
        }

        let db_stereo_ok =
            stereo_dbs.is_empty() || check_double_bond_stereo(&coords_3d, &stereo_dbs);
        let atrop_ok = atropisomers.is_empty() || check_atropisomers(&coords_3d, &atropisomers);
        let no_clash = !has_vdw_clash(&coords_3d, mol);
        let bonds_ok = bond_lengths_reasonable(&coords_3d, mol);

        if planar && db_geom_ok && chiral_ok && db_stereo_ok && atrop_ok && no_clash && bonds_ok {
            let e_per_atom = energy / coords_3d.len() as f64;
            if e_per_atom < MAX_MINIMIZED_E_PER_ATOM {
                return coords_3d;
            }
        }

        if energy < best_energy {
            best_energy = energy;
            best_coords = coords_3d;
        }
    }

    if !best_coords.is_empty() {
        best_coords
    } else {
        let coords_4d = generate_initial_coords_from_bounds(&bounds, &mut rng);
        coords_4d.iter().map(|c| [c[0], c[1], c[2]]).collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn dbg_dihedral_grad_fd() {
        // 4 atoms forming a non-trivial (~40deg) dihedral; check analytical
        // dihedral_gradient_contrib(dedphi=1) == numerical d(phi)/dx.
        let coords: Vec<[f64; 3]> = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.5, 1.0, 0.0],
            [2.0, 1.4, 0.5],
        ];
        let phi0 = dihedral_angle(&coords, 0, 1, 2, 3);
        let (g1, g2, g3, g4) = dihedral_gradient_contrib(&coords, 0, 1, 2, 3, 1.0);
        let grads = [g1, g2, g3, g4];
        let eps = 1e-6;
        let mut max_err = 0.0f64;
        let mut sign_mismatch = 0usize;
        for atom in 0..4 {
            for dim in 0..3 {
                let mut cp = coords.clone();
                cp[atom][dim] += eps;
                let num = (dihedral_angle(&cp, 0, 1, 2, 3) - phi0) / eps;
                let ana = grads[atom][dim];
                eprintln!("  atom{atom} dim{dim}: ana={:+.4} num={:+.4}", ana, num);
                max_err = max_err.max((ana - num).abs());
                if ana.abs() > 1e-4 && num.abs() > 1e-4 && (ana.signum() != num.signum()) {
                    sign_mismatch += 1;
                }
            }
        }
        eprintln!("DBG dihedral: phi0={:.3}rad max|ana-num|={:.2e} sign_mismatches={sign_mismatch}", phi0, max_err);
        assert!(max_err < 1e-3 && sign_mismatch == 0,
            "dihedral_gradient_contrib vs finite-difference mismatch: max_err={max_err:.2e} sign_mismatch={sign_mismatch}");
    }

    #[test]
    fn test_aniline_hh_bounds() {
        let sdf = "Aniline\n     RDKit          3D\n\n 14 14  0  0  0  0  0  0  0  0999 V2000\n   -1.8551    0.3019   -0.2147 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.9433    1.3121   -0.5108 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.4265    1.0872   -0.3490 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.9000   -0.1487    0.0976 C   0  0  0  0  0  0  0  0  0  0  0  0\n    2.2537   -0.3576    0.2878 N   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.0248   -1.1486    0.4072 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.3958   -0.9291    0.2472 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -2.9206    0.4752   -0.3382 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.2957    2.2773   -0.8642 H   0  0  0  0  0  0  0  0  0  0  0  0\n    1.1231    1.8892   -0.5767 H   0  0  0  0  0  0  0  0  0  0  0  0\n    2.5964   -1.2716    0.5480 H   0  0  0  0  0  0  0  0  0  0  0  0\n    2.9224    0.3435    0.0017 H   0  0  0  0  0  0  0  0  0  0  0  0\n    0.3154   -2.1119    0.7767 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -2.1023   -1.7188    0.4874 H   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  2  0\n  2  3  1  0\n  3  4  2  0\n  4  5  1  0\n  4  6  1  0\n  6  7  2  0\n  7  1  1  0\n  1  8  1  0\n  2  9  1  0\n  3 10  1  0\n  5 11  1  0\n  5 12  1  0\n  6 13  1  0\n  7 14  1  0\nM  END";
        let mol = crate::molecule::parser::parse_sdf(sdf).expect("parse");
        let config = ETKDGConfig::default();
        let bounds = build_distance_bounds(&mol, &config);
        
        println!("\nAniline bounds:");
        println!("H10-H11: [{:.4}, {:.4}]", bounds.get_lower(10, 11), bounds.get_upper(10, 11));
        println!("C3-H10:  [{:.4}, {:.4}]", bounds.get_lower(3, 10), bounds.get_upper(3, 10));
        println!("C3-H11:  [{:.4}, {:.4}]", bounds.get_lower(3, 11), bounds.get_upper(3, 11));
        println!("N4-H10:  [{:.4}, {:.4}]", bounds.get_lower(4, 10), bounds.get_upper(4, 10));
        println!("N4-H11:  [{:.4}, {:.4}]", bounds.get_lower(4, 11), bounds.get_upper(4, 11));
        
        // With sp2 N, 120° angle, N-H ~1.04:
        // H-H should be ~2*1.04*sin(60°) = 1.80
        assert!(bounds.get_lower(10, 11) > 1.6, "H-H lower bound too small: {:.4}", bounds.get_lower(10, 11));
    }

    #[test]
    fn test_aniline_geometry() {
        let sdf = "Aniline\n     RDKit          3D\n\n 14 14  0  0  0  0  0  0  0  0999 V2000\n   -1.8551    0.3019   -0.2147 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.9433    1.3121   -0.5108 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.4265    1.0872   -0.3490 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.9000   -0.1487    0.0976 C   0  0  0  0  0  0  0  0  0  0  0  0\n    2.2537   -0.3576    0.2878 N   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.0248   -1.1486    0.4072 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.3958   -0.9291    0.2472 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -2.9206    0.4752   -0.3382 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.2957    2.2773   -0.8642 H   0  0  0  0  0  0  0  0  0  0  0  0\n    1.1231    1.8892   -0.5767 H   0  0  0  0  0  0  0  0  0  0  0  0\n    2.5964   -1.2716    0.5480 H   0  0  0  0  0  0  0  0  0  0  0  0\n    2.9224    0.3435    0.0017 H   0  0  0  0  0  0  0  0  0  0  0  0\n    0.3154   -2.1119    0.7767 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -2.1023   -1.7188    0.4874 H   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  2  0\n  2  3  1  0\n  3  4  2  0\n  4  5  1  0\n  4  6  1  0\n  6  7  2  0\n  7  1  1  0\n  1  8  1  0\n  2  9  1  0\n  3 10  1  0\n  5 11  1  0\n  5 12  1  0\n  6 13  1  0\n  7 14  1  0\nM  END";
        let mol = crate::molecule::parser::parse_sdf(sdf).expect("parse");
        let config = ETKDGConfig::default();
        let coords = generate_initial_coords_with_config(&mol, &config);
        
        // Compute H-N-H angle (atoms 10, 4, 11 in 0-based)
        let n = coords[4];
        let h1 = coords[10];
        let h2 = coords[11];
        let v1 = [h1[0]-n[0], h1[1]-n[1], h1[2]-n[2]];
        let v2 = [h2[0]-n[0], h2[1]-n[1], h2[2]-n[2]];
        let d1 = (v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2]).sqrt();
        let d2 = (v2[0]*v2[0] + v2[1]*v2[1] + v2[2]*v2[2]).sqrt();
        let dot = v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
        let angle = (dot / (d1 * d2)).acos().to_degrees();
        
        // H-H distance
        let hh_dx = h1[0] - h2[0];
        let hh_dy = h1[1] - h2[1];
        let hh_dz = h1[2] - h2[2];
        let hh_dist = (hh_dx*hh_dx + hh_dy*hh_dy + hh_dz*hh_dz).sqrt();
        
        // N-C bond distance (atom 4 to 3)
        let nc_dx = coords[4][0] - coords[3][0];
        let nc_dy = coords[4][1] - coords[3][1];
        let nc_dz = coords[4][2] - coords[3][2];
        let nc_dist = (nc_dx*nc_dx + nc_dy*nc_dy + nc_dz*nc_dz).sqrt();
        
        println!("\nAniline geometry:");
        println!("H-N-H angle: {:.2}°", angle);
        println!("H-H distance: {:.4} Å", hh_dist);
        println!("N-C distance: {:.4} Å", nc_dist);
        println!("N-H1 distance: {:.4} Å", d1);
        println!("N-H2 distance: {:.4} Å", d2);
        
        // Check planarity: distance of N and H's from ring plane
        let ring_atoms = [0, 1, 2, 3, 5, 6];
        let mut center = [0.0f64; 3];
        for &a in &ring_atoms {
            center[0] += coords[a][0];
            center[1] += coords[a][1];
            center[2] += coords[a][2];
        }
        for i in 0..3 { center[i] /= ring_atoms.len() as f64; }
        
        let mut normal = [0.0f64; 3];
        for i in 0..ring_atoms.len() {
            let a = ring_atoms[i];
            let b = ring_atoms[(i+1) % ring_atoms.len()];
            let u = [coords[b][0]-coords[a][0], coords[b][1]-coords[a][1], coords[b][2]-coords[a][2]];
            let v = [center[0]-coords[a][0], center[1]-coords[a][1], center[2]-coords[a][2]];
            normal[0] += u[1]*v[2] - u[2]*v[1];
            normal[1] += u[2]*v[0] - u[0]*v[2];
            normal[2] += u[0]*v[1] - u[1]*v[0];
        }
        let n_norm = (normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2]).sqrt();
        for i in 0..3 { normal[i] /= n_norm; }
        
        let n_dist = (coords[4][0]-center[0])*normal[0] + (coords[4][1]-center[1])*normal[1] + (coords[4][2]-center[2])*normal[2];
        let h1_dist = (coords[10][0]-center[0])*normal[0] + (coords[10][1]-center[1])*normal[1] + (coords[10][2]-center[2])*normal[2];
        let h2_dist = (coords[11][0]-center[0])*normal[0] + (coords[11][1]-center[1])*normal[1] + (coords[11][2]-center[2])*normal[2];
        
        println!("N from plane: {:.4} Å", n_dist);
        println!("H1 from plane: {:.4} Å", h1_dist);
        println!("H2 from plane: {:.4} Å", h2_dist);
        
        // RDKit-like values: H-N-H ~117.5°, N ~-0.03Å, H's ~±0.84Å
        assert!(angle > 110.0 && angle < 125.0, 
                "H-N-H angle out of range: {:.2}° (expected 110-125°)", angle);
    }
}
