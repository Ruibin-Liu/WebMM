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

use crate::molecule::{BondStereo, BondType, Hybridization, Molecule};
use std::collections::HashSet;

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
const MAX_MINIMIZED_E_PER_ATOM: f64 = 0.05;

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
        (x << k) | (x >> (64 - k))
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
}

impl Default for ETKDGConfig {
    fn default() -> Self {
        Self {
            max_attempts: 10,
            convergence_threshold: 1e-6,
            max_iterations: 200,
            vdw_scale: 0.8,
            random_seed: 42,
            force_trans_amides: true,
            use_macrocycle_14config: true,
            use_small_ring_torsions: false,
            use_macrocycle_torsions: true,
            et_version: 2,
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
        Self { lower, upper, n_atoms }
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

    pub fn smooth_triangle_inequality(&mut self) {
        let mut changed = true;
        let epsilon = 1e-6;
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

    pub fn get_lower(&self, i: usize, j: usize) -> f64 {
        self.lower[i.min(j)][i.max(j)]
    }

    pub fn get_upper(&self, i: usize, j: usize) -> f64 {
        self.upper[i.min(j)][i.max(j)]
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
enum Path14Type { Cis, Trans, Other }

#[derive(Debug, Clone)]
struct Path14Config { bid1: usize, bid2: usize, bid3: usize, ptype: Path14Type }

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
        (max_dist >= 1 && self.visited12[pid]) ||
        (max_dist >= 2 && self.visited13[pid]) ||
        (max_dist >= 3 && self.visited14[pid])
    }

}

// ============================================================================
// Radii
// ============================================================================

fn vdw_radius(element: &str) -> f64 {
    match element {
        "H" => 1.20, "C" => 1.70, "N" => 1.55, "O" => 1.52,
        "F" => 1.47, "P" => 1.80, "S" => 1.80, "Cl" => 1.75,
        "Br" => 1.85, "I" => 1.98, "Si" => 2.10, "B" => 1.80,
        _ => 1.70,
    }
}

fn covalent_radius(element: &str) -> f64 {
    match element {
        "H" => 0.31, "C" => 0.76, "N" => 0.71, "O" => 0.66,
        "F" => 0.57, "P" => 1.07, "S" => 1.05, "Cl" => 1.02,
        "Br" => 1.20, "I" => 1.33, "Si" => 1.11, "B" => 0.84,
        _ => 0.76,
    }
}

fn compute_bond_length(element1: &str, element2: &str, bond_type: BondType) -> f64 {
    let r1 = covalent_radius(element1);
    let r2 = covalent_radius(element2);
    let single = r1 + r2;
    match bond_type {
        BondType::Single => single,
        BondType::Double => single * 0.86,
        BondType::Triple => single * 0.78,
        BondType::Aromatic => single * 0.93,
    }
}

fn element_to_atomic_num(element: &str) -> u8 {
    match element {
        "H" => 1, "C" => 6, "N" => 7, "O" => 8, "F" => 9,
        "P" => 15, "S" => 16, "Cl" => 17, "Br" => 35, "I" => 53,
        "Si" => 14, "B" => 5, _ => 6,
    }
}

// ============================================================================
// Geometry helpers
// ============================================================================

fn set_ring_angle(hyb: Hybridization, ring_size: usize) -> f64 {
    if (matches!(hyb, Hybridization::Sp2) && ring_size <= 8) || ring_size == 3 || ring_size == 4 {
        std::f64::consts::PI * (1.0 - 2.0 / ring_size as f64)
    } else if matches!(hyb, Hybridization::Sp3) {
        if ring_size == 5 { 104.0_f64.to_radians() } else { 109.5_f64.to_radians() }
    } else {
        120.0_f64.to_radians()
    }
}

fn compute_13_dist(bl1: f64, bl2: f64, angle: f64) -> f64 {
    (bl1 * bl1 + bl2 * bl2 - 2.0 * bl1 * bl2 * angle.cos()).sqrt()
}

fn compute_14_dist_cis(bl1: f64, bl2: f64, bl3: f64, ba12: f64, ba23: f64) -> f64 {
    let r = bl1*bl1 + bl2*bl2 + bl3*bl3
        - 2.0*bl1*bl2*ba12.cos() - 2.0*bl2*bl3*ba23.cos()
        + 2.0*bl1*bl3*(ba12.cos()*ba23.cos() - ba12.sin()*ba23.sin());
    if r > 0.0 { r.sqrt() } else { 0.5 }
}

fn compute_14_dist_trans(bl1: f64, bl2: f64, bl3: f64, ba12: f64, ba23: f64) -> f64 {
    let r = bl1*bl1 + bl2*bl2 + bl3*bl3
        - 2.0*bl1*bl2*ba12.cos() - 2.0*bl2*bl3*ba23.cos()
        + 2.0*bl1*bl3*(ba12.cos()*ba23.cos() + ba12.sin()*ba23.sin());
    if r > 0.0 { r.sqrt() } else { 0.5 }
}

fn compute_15_dist_cis_cis(d1: f64, d2: f64, d3: f64, d4: f64, ang12: f64, ang23: f64, ang34: f64) -> f64 {
    let dx14 = d2 - d3 * ang23.cos() - d1 * ang12.cos();
    let dy14 = d3 * ang23.sin() - d1 * ang12.sin();
    let d14 = (dx14 * dx14 + dy14 * dy14).sqrt();
    let mut cval = (d3 - d2 * ang23.cos() + d1 * (ang12 + ang23).cos()) / d14;
    cval = cval.clamp(-1.0, 1.0);
    let ang143 = cval.acos();
    let ang145 = ang34 - ang143;
    compute_13_dist(d14, d4, ang145)
}

fn compute_15_dist_cis_trans(d1: f64, d2: f64, d3: f64, d4: f64, ang12: f64, ang23: f64, ang34: f64) -> f64 {
    let dx14 = d2 - d3 * ang23.cos() - d1 * ang12.cos();
    let dy14 = d3 * ang23.sin() - d1 * ang12.sin();
    let d14 = (dx14 * dx14 + dy14 * dy14).sqrt();
    let mut cval = (d3 - d2 * ang23.cos() + d1 * (ang12 + ang23).cos()) / d14;
    cval = cval.clamp(-1.0, 1.0);
    let ang143 = cval.acos();
    let ang145 = ang34 + ang143;
    compute_13_dist(d14, d4, ang145)
}

fn compute_15_dist_trans_trans(d1: f64, d2: f64, d3: f64, d4: f64, ang12: f64, ang23: f64, ang34: f64) -> f64 {
    let dx14 = d2 - d3 * ang23.cos() - d1 * ang12.cos();
    let dy14 = d3 * ang23.sin() + d1 * ang12.sin();
    let d14 = (dx14 * dx14 + dy14 * dy14).sqrt();
    let mut cval = (d3 - d2 * ang23.cos() + d1 * (ang12 - ang23).cos()) / d14;
    cval = cval.clamp(-1.0, 1.0);
    let ang143 = cval.acos();
    let ang145 = ang34 + ang143;
    compute_13_dist(d14, d4, ang145)
}

fn compute_15_dist_trans_cis(d1: f64, d2: f64, d3: f64, d4: f64, ang12: f64, ang23: f64, ang34: f64) -> f64 {
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
    mol.bonds.iter().enumerate().find(|(_, b)| (b.atom1 == i && b.atom2 == j) || (b.atom1 == j && b.atom2 == i)).map(|(idx, _)| idx)
}

fn is_larger_sp2_atom(atom_idx: usize, mol: &Molecule) -> bool {
    let atom = &mol.atoms[atom_idx];
    let hyb = crate::molecule::graph::determine_hybridization(atom_idx, mol);
    let atomic_num = element_to_atomic_num(&atom.symbol);
    atomic_num > 13 && matches!(hyb, Hybridization::Sp2)
}

fn is_amide_bond(mol: &Molecule, a: usize, b: usize) -> bool {
    let sym_a = &mol.atoms[a].symbol;
    let sym_b = &mol.atoms[b].symbol;
    if !((sym_a == "C" && sym_b == "N") || (sym_a == "N" && sym_b == "C")) { return false; }
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
    if !((sym_a == "C" && sym_b == "O") || (sym_a == "O" && sym_b == "C")) { return false; }
    let (c_atom, o_atom) = if sym_a == "C" { (a, b) } else { (b, a) };
    // Carbon must have a double bond to another oxygen (C=O)
    let has_carbonyl = mol.adjacency[c_atom].iter().any(|&nbr| {
        nbr != o_atom && mol.bonds.iter().any(|bond| {
            (bond.atom1 == c_atom && bond.atom2 == nbr && bond.bond_type == BondType::Double)
                || (bond.atom2 == c_atom && bond.atom1 == nbr && bond.bond_type == BondType::Double)
        })
    });
    if !has_carbonyl { return false; }
    // Single bond oxygen must have no other carbonyl neighbors (distinguish from acid anhydride)
    mol.adjacency[o_atom].iter().all(|&nbr| {
        nbr == c_atom || !mol.bonds.iter().any(|bond| {
            (bond.atom1 == o_atom && bond.atom2 == nbr && bond.bond_type == BondType::Double)
                || (bond.atom2 == o_atom && bond.atom1 == nbr && bond.bond_type == BondType::Double)
        })
    })
}

fn is_double_bond(mol: &Molecule, i: usize, j: usize) -> bool {
    mol.bonds.iter().any(|b| ((b.atom1 == i && b.atom2 == j) || (b.atom1 == j && b.atom2 == i)) && b.bond_type == BondType::Double)
}

fn find_non_double_neighbor(mol: &Molecule, atom: usize, exclude: usize) -> Option<usize> {
    mol.adjacency[atom].iter().copied().find(|&n| n != exclude && !is_double_bond(mol, atom, n))
}

// ============================================================================
// Bounds Matrix Builder
// ============================================================================

fn build_distance_bounds(mol: &Molecule, config: &ETKDGConfig) -> DistanceBounds {
    let n_atoms = mol.atoms.len();
    let mut bounds = DistanceBounds::new(n_atoms);
    let n_bonds = mol.bonds.len();
    let mut accum = ComputedData::new(n_atoms, n_bonds);

    // 1-2 bounds
    for (bond_idx, bond) in mol.bonds.iter().enumerate() {
        let i = bond.atom1;
        let j = bond.atom2;
        let bl = compute_bond_length(&mol.atoms[i].symbol, &mol.atoms[j].symbol, bond.bond_type);
        accum.bond_lengths[bond_idx] = bl;
        bounds.check_and_set(i, j, bl - DIST12_DELTA, bl + DIST12_DELTA);
        accum.visited12[i.min(j) * n_atoms + i.max(j)] = true;
    }

    // 1-3 bounds
    let rings = crate::molecule::graph::find_rings(mol);
    let mut visited_centers = vec![0usize; n_atoms];
    let mut angle_taken = vec![0.0f64; n_atoms];
    let mut done_paths = vec![false; n_bonds * n_bonds];

    let mut sorted_rings = rings.clone();
    sorted_rings.sort_by_key(|r| r.len());

    for ring in &sorted_rings {
        let rsize = ring.len();
        if rsize < 3 { continue; }
        let mut aid1 = ring[rsize - 1];
        for i in 0..rsize {
            let aid2 = ring[i];
            let aid3 = ring[(i + 1) % rsize];
            if let (Some(bid1), Some(bid2)) = (find_bond_index(mol, aid1, aid2), find_bond_index(mol, aid2, aid3)) {
                let id1 = bid1 * n_bonds + bid2;
                let id2 = bid2 * n_bonds + bid1;
                let pid = aid1.min(aid3) * n_atoms + aid1.max(aid3);
                if !done_paths[id1] && !done_paths[id2] {
                    let hyb = crate::molecule::graph::determine_hybridization(aid2, mol);
                    let angle = set_ring_angle(hyb, rsize);
                    let dl = compute_13_dist(accum.bond_lengths[bid1], accum.bond_lengths[bid2], angle);
                    let dist_tol = if is_larger_sp2_atom(aid1, mol) || is_larger_sp2_atom(aid2, mol) || is_larger_sp2_atom(aid3, mol) {
                        DIST13_TOL * 2.0
                    } else { DIST13_TOL };
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
        if n13 == visited_centers[aid2] { continue; }
        let hyb = crate::molecule::graph::determine_hybridization(aid2, mol);
        let neighbors: Vec<usize> = mol.adjacency[aid2].iter().copied().collect();
        for i1 in 0..neighbors.len() {
            for i2 in 0..i1 {
                let aid1 = neighbors[i2];
                let aid3 = neighbors[i1];
                if let (Some(bid1), Some(bid2)) = (find_bond_index(mol, aid1, aid2), find_bond_index(mol, aid2, aid3)) {
                    if accum.bond_angles[bid1][bid2] >= 0.0 { continue; }
                    let angle = if visited_centers[aid2] >= 1 {
                        if matches!(hyb, Hybridization::Sp2) {
                            (2.0 * std::f64::consts::PI - angle_taken[aid2]) / (n13 - visited_centers[aid2]) as f64
                        } else if matches!(hyb, Hybridization::Sp3) {
                            let mut a = 109.5_f64.to_radians();
                            for ring in &rings {
                                if ring.contains(&aid2) && ring.len() == 3 { a = 116.0_f64.to_radians(); }
                                else if ring.contains(&aid2) && ring.len() == 4 { a = 112.0_f64.to_radians(); }
                            }
                            a
                        } else { 120.0_f64.to_radians() }
                    } else {
                        match hyb {
                            Hybridization::Sp1 => std::f64::consts::PI,
                            Hybridization::Sp2 => 120.0_f64.to_radians(),
                            Hybridization::Sp3 => 109.5_f64.to_radians(),
                        }
                    };
                    let pid = aid1.min(aid3) * n_atoms + aid1.max(aid3);
                    if !accum.visited12[pid] {
                        let dl = compute_13_dist(accum.bond_lengths[bid1], accum.bond_lengths[bid2], angle);
                        let dist_tol = if is_larger_sp2_atom(aid1, mol) || is_larger_sp2_atom(aid2, mol) || is_larger_sp2_atom(aid3, mol) {
                            DIST13_TOL * 2.0
                        } else { DIST13_TOL };
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
    for i in 0..n_atoms { for j in 0..n_atoms { dist_mat[i][j] = if i == j { 0 } else { usize::MAX }; } }
    for bond in &mol.bonds { dist_mat[bond.atom1][bond.atom2] = 1; dist_mat[bond.atom2][bond.atom1] = 1; }
    for k in 0..n_atoms { for i in 0..n_atoms { for j in 0..n_atoms {
        if dist_mat[i][k] != usize::MAX && dist_mat[k][j] != usize::MAX {
            let d = dist_mat[i][k] + dist_mat[k][j];
            if d < dist_mat[i][j] { dist_mat[i][j] = d; }
        }
    }}}

    // 1-4 bounds
    let mut done_14_paths = HashSet::new();
    for ring in &sorted_rings {
        let rsize = ring.len();
        if rsize < 4 { continue; }
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
            if atm2 < 0 || atm3 < 0 { bid1 = bid2; continue; }
            let aid1 = mol.bonds[bid1].atom1 + mol.bonds[bid1].atom2 - atm2 as usize;
            let aid4 = mol.bonds[bid3].atom1 + mol.bonds[bid3].atom2 - atm3 as usize;
            let pid = aid1.min(aid4) * n_atoms + aid1.max(aid4);
            if accum.visited_bound(pid, 2) { bid1 = bid2; continue; }
            if dist_mat[aid1][aid4] < 3 { bid1 = bid2; continue; }
            let bl1 = accum.bond_lengths[bid1];
            let bl2 = accum.bond_lengths[bid2];
            let bl3 = accum.bond_lengths[bid3];
            let ba12 = accum.bond_angles[bid1][bid2];
            let ba23 = accum.bond_angles[bid2][bid3];
            if ba12 < 0.0 || ba23 < 0.0 { bid1 = bid2; continue; }
            let hyb2 = crate::molecule::graph::determine_hybridization(atm2 as usize, mol);
            let hyb3 = crate::molecule::graph::determine_hybridization(atm3 as usize, mol);
            
            let (prefer_cis, prefer_trans) = if config.use_macrocycle_14config && rsize >= 9 {
                (false, false) // macrocycles don't assume cis
            } else {
                let is_cis = rsize <= 8 && matches!(hyb2, Hybridization::Sp2) && matches!(hyb3, Hybridization::Sp2);
                (is_cis, false)
            };
            
            let (dl, du, ptype) = if prefer_cis {
                let d = compute_14_dist_cis(bl1, bl2, bl3, ba12, ba23);
                (d - GEN_DIST_TOL, d + GEN_DIST_TOL, Path14Type::Cis)
            } else if prefer_trans {
                let d = compute_14_dist_trans(bl1, bl2, bl3, ba12, ba23);
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
            accum.paths14.push(Path14Config { bid1, bid2, bid3, ptype });
            let path_id = bid1 as u64 * n_bonds as u64 * n_bonds as u64 + bid2 as u64 * n_bonds as u64 + bid3 as u64;
            match ptype {
                Path14Type::Cis => { accum.cis_paths.insert(path_id); }
                Path14Type::Trans => { accum.trans_paths.insert(path_id); }
                _ => {}
            }
            bid1 = bid2;
        }
    }

    for (bid2, bond2) in mol.bonds.iter().enumerate() {
        let aid2 = bond2.atom1;
        let aid3 = bond2.atom2;
        for &nbr2 in &mol.adjacency[aid2] {
            if nbr2 == aid3 { continue; }
            let Some(bid1) = find_bond_index(mol, nbr2, aid2) else { continue; };
            for &nbr3 in &mol.adjacency[aid3] {
                if nbr3 == aid2 || nbr3 == nbr2 { continue; }
                let Some(bid3) = find_bond_index(mol, aid3, nbr3) else { continue; };
                let id1 = (bid1, bid2, bid3);
                let id2 = (bid3, bid2, bid1);
                if done_14_paths.contains(&id1) || done_14_paths.contains(&id2) { continue; }
                done_14_paths.insert(id1);
                done_14_paths.insert(id2);
                let aid1 = nbr2;
                let aid4 = nbr3;
                let pid = aid1.min(aid4) * n_atoms + aid1.max(aid4);
                if accum.visited_bound(pid, 2) { continue; }
                let bl1 = accum.bond_lengths[bid1];
                let bl2 = accum.bond_lengths[bid2];
                let bl3 = accum.bond_lengths[bid3];
                let ba12 = accum.bond_angles[bid1][bid2];
                let ba23 = accum.bond_angles[bid2][bid3];
                if ba12 < 0.0 || ba23 < 0.0 { continue; }
                
                let (dl, du, ptype) = if config.force_trans_amides && is_amide_bond(mol, aid2, aid3) {
                    let dt = compute_14_dist_trans(bl1, bl2, bl3, ba12, ba23);
                    (dt - GEN_DIST_TOL, dt + GEN_DIST_TOL, Path14Type::Trans)
                } else if config.force_trans_amides && is_ester_bond(mol, aid2, aid3) {
                    let dt = compute_14_dist_trans(bl1, bl2, bl3, ba12, ba23);
                    (dt - GEN_DIST_TOL, dt + GEN_DIST_TOL, Path14Type::Trans)
                } else {
                    let dc = compute_14_dist_cis(bl1, bl2, bl3, ba12, ba23);
                    let dt = compute_14_dist_trans(bl1, bl2, bl3, ba12, ba23);
                    let mut dl = dc.min(dt);
                    let mut du = dc.max(dt);
                    if (du - dl).abs() < DIST12_DELTA { dl -= GEN_DIST_TOL; du += GEN_DIST_TOL; }
                    (dl, du, Path14Type::Other)
                };
                bounds.check_and_set(aid1, aid4, dl, du);
                accum.visited14[pid] = true;
                accum.paths14.push(Path14Config { bid1, bid2, bid3, ptype });
            }
        }
    }

    // 1-5 bounds
    let paths14_copy: Vec<Path14Config> = accum.paths14.clone();
    for path14 in paths14_copy {
        set_15_bounds_helper(mol, &mut accum, &mut bounds, &path14, &dist_mat, n_atoms, n_bonds, false);
        let rev = Path14Config { bid1: path14.bid3, bid2: path14.bid2, bid3: path14.bid1, ptype: path14.ptype };
        set_15_bounds_helper(mol, &mut accum, &mut bounds, &rev, &dist_mat, n_atoms, n_bonds, true);
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
                let lb = if (h_bond_donor_h[i] && h_bond_acceptor[j]) || (h_bond_acceptor[i] && h_bond_donor_h[j]) {
                    H_BOND_LENGTH
                } else if d == 4 { VDW_SCALE_15 * (vw1 + vw2) }
                else if d == 5 { (VDW_SCALE_15 + 0.5 * (1.0 - VDW_SCALE_15)) * (vw1 + vw2) }
                else { vw1 + vw2 };
                bounds.set_lower(i, j, lb);
            }
        }
    }

    bounds
}

fn set_15_bounds_helper(mol: &Molecule, accum: &mut ComputedData, bounds: &mut DistanceBounds, path14: &Path14Config, dist_mat: &[Vec<usize>], n_atoms: usize, n_bonds: usize, reversed: bool) {
    let bid1 = path14.bid1;
    let bid2 = path14.bid2;
    let bid3 = path14.bid3;
    let aid2 = accum.bond_adj[bid1][bid2] as usize;
    let aid3 = accum.bond_adj[bid2][bid3] as usize;
    let aid1 = if mol.bonds[bid1].atom1 == aid2 { mol.bonds[bid1].atom2 } else { mol.bonds[bid1].atom1 };
    let aid4 = if mol.bonds[bid3].atom1 == aid3 { mol.bonds[bid3].atom2 } else { mol.bonds[bid3].atom1 };
    
    let d1 = accum.bond_lengths[bid1];
    let d2 = accum.bond_lengths[bid2];
    let d3 = accum.bond_lengths[bid3];
    let ang12 = accum.bond_angles[bid1][bid2];
    let ang23 = accum.bond_angles[bid2][bid3];
    if ang12 < 0.0 || ang23 < 0.0 { return; }
    
    for i in 0..mol.bonds.len() {
        if accum.bond_adj[bid3][i] == aid4 as i32 {
            let aid5 = if mol.bonds[i].atom1 == aid4 { mol.bonds[i].atom2 } else { mol.bonds[i].atom1 };
            if aid1 == aid5 { continue; }
            let pid = aid1.min(aid5) * n_atoms + aid1.max(aid5);
            if accum.visited_bound(pid, 3) { continue; }
            if dist_mat[aid1.max(aid5)][aid1.min(aid5)] < 4 { continue; }
            if accum.set15[pid] { continue; }
            
            let d4 = accum.bond_lengths[i];
            let ang34 = accum.bond_angles[bid3][i];
            if ang34 < 0.0 { continue; }
            
            let (dl, du) = if reversed {
                match path14.ptype {
                    Path14Type::Cis => {
                        let dl = compute_15_dist_cis_cis(d4, d3, d2, d1, ang34, ang23, ang12) - DIST15_TOL;
                        let du = compute_15_dist_cis_trans(d4, d3, d2, d1, ang34, ang23, ang12) + DIST15_TOL;
                        (dl, du)
                    }
                    Path14Type::Trans => {
                        let dl = compute_15_dist_trans_cis(d4, d3, d2, d1, ang34, ang23, ang12) - DIST15_TOL;
                        let du = compute_15_dist_trans_trans(d4, d3, d2, d1, ang34, ang23, ang12) + DIST15_TOL;
                        (dl, du)
                    }
                    Path14Type::Other => {
                        let path_id = bid2 as u64 * n_bonds as u64 * n_bonds as u64 + bid3 as u64 * n_bonds as u64 + i as u64;
                        if accum.cis_paths.contains(&path_id) {
                            let dl = compute_15_dist_cis_cis(d4, d3, d2, d1, ang34, ang23, ang12) - DIST15_TOL;
                            let du = compute_15_dist_cis_trans(d4, d3, d2, d1, ang34, ang23, ang12) + DIST15_TOL;
                            (dl, du)
                        } else if accum.trans_paths.contains(&path_id) {
                            let dl = compute_15_dist_trans_cis(d4, d3, d2, d1, ang34, ang23, ang12) - DIST15_TOL;
                            let du = compute_15_dist_trans_trans(d4, d3, d2, d1, ang34, ang23, ang12) + DIST15_TOL;
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
                        let dl = compute_15_dist_cis_cis(d1, d2, d3, d4, ang12, ang23, ang34) - DIST15_TOL;
                        let du = compute_15_dist_cis_trans(d1, d2, d3, d4, ang12, ang23, ang34) + DIST15_TOL;
                        (dl, du)
                    }
                    Path14Type::Trans => {
                        let dl = compute_15_dist_trans_cis(d1, d2, d3, d4, ang12, ang23, ang34) - DIST15_TOL;
                        let du = compute_15_dist_trans_trans(d1, d2, d3, d4, ang12, ang23, ang34) + DIST15_TOL;
                        (dl, du)
                    }
                    Path14Type::Other => {
                        let path_id = bid2 as u64 * n_bonds as u64 * n_bonds as u64 + bid3 as u64 * n_bonds as u64 + i as u64;
                        if accum.cis_paths.contains(&path_id) {
                            let dl = compute_15_dist_cis_cis(d1, d2, d3, d4, ang12, ang23, ang34) - DIST15_TOL;
                            let du = compute_15_dist_cis_trans(d1, d2, d3, d4, ang12, ang23, ang34) + DIST15_TOL;
                            (dl, du)
                        } else if accum.trans_paths.contains(&path_id) {
                            let dl = compute_15_dist_trans_cis(d1, d2, d3, d4, ang12, ang23, ang34) - DIST15_TOL;
                            let du = compute_15_dist_trans_trans(d1, d2, d3, d4, ang12, ang23, ang34) + DIST15_TOL;
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
    if n < 2 { return vec![[0.0; 4]; n]; }
    let mut dist = vec![vec![0.0f64; n]; n];
    for i in 0..n { for j in (i + 1)..n {
        let lo = bounds.lower[i][j]; let hi = bounds.upper[i][j];
        let t = lo + (hi - lo) * rng.random_f64();
        dist[i][j] = t; dist[j][i] = t;
    }}
    let n_f = n as f64;
    let mut d2 = vec![vec![0.0f64; n]; n];
    let mut row_sums = vec![0.0f64; n];
    let mut total_sum = 0.0f64;
    for i in 0..n { for j in 0..n {
        d2[i][j] = dist[i][j] * dist[i][j];
        row_sums[i] += d2[i][j]; total_sum += d2[i][j];
    }}
    let row_means: Vec<f64> = row_sums.iter().map(|s| s / n_f).collect();
    let grand_mean = total_sum / (n_f * n_f);
    let mut b = vec![vec![0.0f64; n]; n];
    for i in 0..n { for j in 0..n {
        b[i][j] = -0.5 * (d2[i][j] - row_means[i] - row_means[j] + grand_mean);
    }}
    let (eigenvalues, eigenvectors) = jacobi_eigen(&b, 100);
    let mut indices: Vec<usize> = (0..n).collect();
    indices.sort_by(|&a, &b| eigenvalues[b].partial_cmp(&eigenvalues[a]).unwrap());
    let dim = 4.min(n);
    let mut coords_4d = vec![[0.0f64; 4]; n];
    for i in 0..n { for d in 0..dim {
        if d < indices.len() {
            let idx = indices[d]; let ev = eigenvalues[idx];
            coords_4d[i][d] = if ev > 0.0 { eigenvectors[idx][i] * ev.sqrt() } else { 1.0 - 2.0 * rng.random_f64() };
        }
    }}
    coords_4d
}

fn jacobi_eigen(matrix: &[Vec<f64>], max_sweeps: usize) -> (Vec<f64>, Vec<Vec<f64>>) {
    let n = matrix.len();
    let mut a = matrix.to_vec();
    let mut v = vec![vec![0.0f64; n]; n];
    for i in 0..n { v[i][i] = 1.0; }
    let mut d = vec![0.0f64; n];
    for i in 0..n { d[i] = a[i][i]; }
    let mut b = d.clone();
    let mut z = vec![0.0f64; n];
    for _sweep in 0..max_sweeps {
        let mut sum = 0.0;
        for i in 0..n - 1 { for j in i + 1..n { sum += a[i][j].abs(); } }
        if sum < 1e-12 { break; }
        let threshold = if _sweep < 3 { 0.2 * sum / (n * n) as f64 } else { 0.0 };
        for p in 0..n - 1 { for q in p + 1..n {
            let apq = a[p][q].abs(); let g = 100.0 * apq;
            if _sweep > 3 && d[p].abs() + g == d[p].abs() && d[q].abs() + g == d[q].abs() { a[p][q] = 0.0; continue; }
            if apq <= threshold { continue; }
            let h = d[q] - d[p];
            let t = if h.abs() + g == h.abs() { a[p][q] / h } else {
                let theta = 0.5 * h / a[p][q];
                let mut tt = 1.0 / (theta.abs() + (1.0 + theta * theta).sqrt());
                if theta < 0.0 { tt = -tt; }
                tt
            };
            let c = 1.0 / (1.0 + t * t).sqrt(); let s = t * c; let tau = s / (1.0 + c); let h = t * a[p][q];
            z[p] -= h; z[q] += h; d[p] -= h; d[q] += h; a[p][q] = 0.0;
            for r in 0..p { let g = a[r][p]; let h = a[r][q]; a[r][p] = g - s * (h + g * tau); a[r][q] = h + s * (g - h * tau); }
            for r in p + 1..q { let g = a[p][r]; let h = a[r][q]; a[p][r] = g - s * (h + g * tau); a[r][q] = h + s * (g - h * tau); }
            for r in q + 1..n { let g = a[p][r]; let h = a[q][r]; a[p][r] = g - s * (h + g * tau); a[q][r] = h + s * (g - h * tau); }
            for r in 0..n { let g = v[r][p]; let h = v[r][q]; v[r][p] = g - s * (h + g * tau); v[r][q] = h + s * (g - h * tau); }
        }}
        for i in 0..n { b[i] += z[i]; d[i] = b[i]; z[i] = 0.0; }
    }
    let mut evecs = vec![vec![0.0f64; n]; n];
    for i in 0..n { for r in 0..n { evecs[i][r] = v[r][i]; } }
    (d, evecs)
}

fn minimize_fourth_dimension(coords_4d: &mut [[f64; 4]], bounds: &DistanceBounds, max_iter: usize) {
    let n = coords_4d.len();
    const K_W: f64 = 1.0;
    for _ in 0..max_iter {
        let mut grad = vec![[0.0f64; 4]; n];
        for i in 0..n {
            for j in (i + 1)..n {
                let dx = coords_4d[i][0] - coords_4d[j][0];
                let dy = coords_4d[i][1] - coords_4d[j][1];
                let dz = coords_4d[i][2] - coords_4d[j][2];
                let dw = coords_4d[i][3] - coords_4d[j][3];
                let d4 = (dx*dx + dy*dy + dz*dz + dw*dw).sqrt().max(1e-10);
                let lo = bounds.lower[i][j];
                let hi = bounds.upper[i][j];
                if hi >= MAX_UPPER { continue; }
                let range = (hi - lo).max(0.01);
                let w = 1.0 / range;
                let lo_viol = (lo - d4).max(0.0);
                let hi_viol = (d4 - hi).max(0.0);
                let f = -2.0 * w * (lo_viol - hi_viol) / d4;
                for dim in 0..4 {
                    let diff = coords_4d[i][dim] - coords_4d[j][dim];
                    grad[i][dim] += f * diff;
                    grad[j][dim] -= f * diff;
                }
            }
        }
        for i in 0..n {
            grad[i][3] += 2.0 * K_W * coords_4d[i][3];
        }
        let max_g = grad.iter().map(|g| (g[0]*g[0] + g[1]*g[1] + g[2]*g[2] + g[3]*g[3]).sqrt()).fold(0.0f64, f64::max);
        if max_g < 1e-6 { break; }
        let step = 0.01 / max_g.max(1e-10);
        for i in 0..n {
            for dim in 0..4 {
                coords_4d[i][dim] -= step * grad[i][dim];
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
}

fn find_chiral_centers(mol: &Molecule) -> (Vec<ChiralCenter>, Vec<ChiralCenter>) {
    let mut chiral = Vec::new();
    let mut tetrahedral = Vec::new();
    for atom_idx in 0..mol.atoms.len() {
        if mol.atoms[atom_idx].symbol == "H" { continue; }
        let hyb = crate::molecule::graph::determine_hybridization(atom_idx, mol);
        if !matches!(hyb, Hybridization::Sp3) { continue; }
        let neighbors: Vec<usize> = mol.adjacency[atom_idx].iter().copied().filter(|&n| mol.atoms[n].symbol != "H").collect();
        if neighbors.len() < 3 { continue; }
        let mut nbrs = neighbors.clone();
        let (vol_lower, vol_upper) = if nbrs.len() == 4 { (5.0, 100.0) } else { nbrs.push(atom_idx); (2.0, 100.0) };
        let cc = ChiralCenter { central: atom_idx, neighbors: [nbrs[0], nbrs[1], nbrs[2], nbrs[3]], vol_lower, vol_upper };
        tetrahedral.push(cc.clone());
        // Basic topological chirality: 4 distinct non-H neighbors by element symbol
        if nbrs.len() == 4 {
            let mut distinct = true;
            for i in 1..4 {
                for j in 0..i {
                    if mol.atoms[nbrs[i]].symbol == mol.atoms[nbrs[j]].symbol {
                        distinct = false; break;
                    }
                }
                if !distinct { break; }
            }
            if distinct {
                chiral.push(cc);
            }
        }
    }
    (chiral, tetrahedral)
}

fn chiral_volume(coords: &[[f64; 3]], c: usize, n1: usize, n2: usize, n3: usize) -> f64 {
    let v1 = [coords[n1][0]-coords[c][0], coords[n1][1]-coords[c][1], coords[n1][2]-coords[c][2]];
    let v2 = [coords[n2][0]-coords[c][0], coords[n2][1]-coords[c][1], coords[n2][2]-coords[c][2]];
    let v3 = [coords[n3][0]-coords[c][0], coords[n3][1]-coords[c][1], coords[n3][2]-coords[c][2]];
    let cross = [v2[1]*v3[2]-v2[2]*v3[1], v2[2]*v3[0]-v2[0]*v3[2], v2[0]*v3[1]-v2[1]*v3[0]];
    v1[0]*cross[0] + v1[1]*cross[1] + v1[2]*cross[2]
}

fn normalize_vec(v: [f64; 3]) -> [f64; 3] {
    let n = (v[0]*v[0] + v[1]*v[1] + v[2]*v[2]).sqrt();
    if n < 1e-10 { return v; }
    [v[0]/n, v[1]/n, v[2]/n]
}

trait Vec3Ops { fn dot(self, other: [f64; 3]) -> f64; }
impl Vec3Ops for [f64; 3] { fn dot(self, other: [f64; 3]) -> f64 { self[0]*other[0] + self[1]*other[1] + self[2]*other[2] } }

fn cross_product(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]]
}

fn volume_test(coords: &[[f64; 3]], cc: &ChiralCenter) -> bool {
    let p0 = [coords[cc.central][0], coords[cc.central][1], coords[cc.central][2]];
    let p1 = [coords[cc.neighbors[0]][0], coords[cc.neighbors[0]][1], coords[cc.neighbors[0]][2]];
    let p2 = [coords[cc.neighbors[1]][0], coords[cc.neighbors[1]][1], coords[cc.neighbors[1]][2]];
    let p3 = [coords[cc.neighbors[2]][0], coords[cc.neighbors[2]][1], coords[cc.neighbors[2]][2]];
    let p4 = [coords[cc.neighbors[3]][0], coords[cc.neighbors[3]][1], coords[cc.neighbors[3]][2]];
    let v1 = normalize_vec([p0[0]-p1[0], p0[1]-p1[1], p0[2]-p1[2]]);
    let v2 = normalize_vec([p0[0]-p2[0], p0[1]-p2[1], p0[2]-p2[2]]);
    let v3 = normalize_vec([p0[0]-p3[0], p0[1]-p3[1], p0[2]-p3[2]]);
    let v4 = normalize_vec([p0[0]-p4[0], p0[1]-p4[1], p0[2]-p4[2]]);
    let cross = cross_product(v1, v2);
    if cross.dot(v3).abs() < MIN_TETRAHEDRAL_CHIRAL_VOL { return false; }
    let cross = cross_product(v1, v2);
    if cross.dot(v4).abs() < MIN_TETRAHEDRAL_CHIRAL_VOL { return false; }
    let cross = cross_product(v1, v3);
    if cross.dot(v4).abs() < MIN_TETRAHEDRAL_CHIRAL_VOL { return false; }
    let cross = cross_product(v2, v3);
    if cross.dot(v4).abs() < MIN_TETRAHEDRAL_CHIRAL_VOL { return false; }
    true
}

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
        let vol = chiral_volume(coords, cc.central, cc.neighbors[0], cc.neighbors[1], cc.neighbors[2]);
        if cc.vol_lower > 0.0 && vol < cc.vol_lower && (vol / cc.vol_lower < 0.8 || have_opposite_sign(vol, cc.vol_lower)) { return false; }
        if cc.vol_upper < 0.0 && vol > cc.vol_upper && (vol / cc.vol_upper < 0.8 || have_opposite_sign(vol, cc.vol_upper)) { return false; }
    }
    true
}

fn have_opposite_sign(a: f64, b: f64) -> bool { (a < 0.0) != (b < 0.0) }

fn center_in_volume(coords: &[[f64; 3]], cc: &ChiralCenter) -> bool {
    if cc.central == cc.neighbors[3] { return true; }
    same_side(coords[cc.neighbors[0]], coords[cc.neighbors[1]], coords[cc.neighbors[2]], coords[cc.neighbors[3]], coords[cc.central])
        && same_side(coords[cc.neighbors[1]], coords[cc.neighbors[2]], coords[cc.neighbors[3]], coords[cc.neighbors[0]], coords[cc.central])
        && same_side(coords[cc.neighbors[2]], coords[cc.neighbors[3]], coords[cc.neighbors[0]], coords[cc.neighbors[1]], coords[cc.central])
        && same_side(coords[cc.neighbors[3]], coords[cc.neighbors[0]], coords[cc.neighbors[1]], coords[cc.neighbors[2]], coords[cc.central])
}

fn same_side(v1: [f64; 3], v2: [f64; 3], v3: [f64; 3], v4: [f64; 3], p0: [f64; 3]) -> bool {
    let normal = cross_product([v2[0]-v1[0], v2[1]-v1[1], v2[2]-v1[2]], [v3[0]-v1[0], v3[1]-v1[1], v3[2]-v1[2]]);
    let d1 = normal[0]*(v4[0]-v1[0]) + normal[1]*(v4[1]-v1[1]) + normal[2]*(v4[2]-v1[2]);
    let d2 = normal[0]*(p0[0]-v1[0]) + normal[1]*(p0[1]-v1[1]) + normal[2]*(p0[2]-v1[2]);
    if d1.abs() < TETRAHEDRAL_CENTERINVOLUME_TOL || d2.abs() < TETRAHEDRAL_CENTERINVOLUME_TOL { return false; }
    !((d1 < 0.0) ^ (d2 < 0.0))
}

// ============================================================================
// Double Bond Stereo
// ============================================================================

#[derive(Debug, Clone)]
struct StereoDoubleBond { atoms: [usize; 4], sign: i8 }

fn find_stereo_double_bonds(mol: &Molecule) -> (Vec<(usize, usize, usize)>, Vec<StereoDoubleBond>) {
    let mut double_bond_ends = Vec::new();
    let mut stereo_db = Vec::new();
    for bond in &mol.bonds {
        if bond.bond_type != BondType::Double { continue; }
        let a1 = bond.atom1; let a2 = bond.atom2;
        for &nbr in &mol.adjacency[a1] { if nbr != a2 && !is_double_bond(mol, a1, nbr) { double_bond_ends.push((nbr, a1, a2)); } }
        for &nbr in &mol.adjacency[a2] { if nbr != a1 && !is_double_bond(mol, a2, nbr) { double_bond_ends.push((nbr, a2, a1)); } }
        if bond.stereo != BondStereo::None {
            if let (Some(n1), Some(n2)) = (find_non_double_neighbor(mol, a1, a2), find_non_double_neighbor(mol, a2, a1)) {
                let sign = match bond.stereo { BondStereo::Trans => 1, BondStereo::Cis => -1, _ => 1 };
                stereo_db.push(StereoDoubleBond { atoms: [n1, a1, a2, n2], sign });
            }
        }
    }
    (double_bond_ends, stereo_db)
}

fn dihedral_angle(coords: &[[f64; 3]], i: usize, j: usize, k: usize, l: usize) -> f64 {
    let b1 = [coords[j][0]-coords[i][0], coords[j][1]-coords[i][1], coords[j][2]-coords[i][2]];
    let b2 = [coords[k][0]-coords[j][0], coords[k][1]-coords[j][1], coords[k][2]-coords[j][2]];
    let b3 = [coords[l][0]-coords[k][0], coords[l][1]-coords[k][1], coords[l][2]-coords[k][2]];
    let n1 = cross_product(b1, b2);
    let n2 = cross_product(b2, b3);
    let b2_norm = (b2[0]*b2[0] + b2[1]*b2[1] + b2[2]*b2[2]).sqrt();
    if b2_norm < 1e-10 { return 0.0; }
    let m1 = cross_product(n1, [b2[0]/b2_norm, b2[1]/b2_norm, b2[2]/b2_norm]);
    let x = n1.dot(n2); let y = m1.dot(n2);
    y.atan2(x)
}

fn check_double_bond_stereo(coords: &[[f64; 3]], stereo_dbs: &[StereoDoubleBond]) -> bool {
    for sdb in stereo_dbs {
        let dihedral = dihedral_angle(coords, sdb.atoms[0], sdb.atoms[1], sdb.atoms[2], sdb.atoms[3]);
        if (dihedral - std::f64::consts::FRAC_PI_2) * (sdb.sign as f64) < 0.0 { return false; }
    }
    true
}

fn double_bond_geometry_checks(coords: &[[f64; 3]], double_bond_ends: &[(usize, usize, usize)]) -> bool {
    const LINEAR_TOL: f64 = 1e-3;
    for &(a0, a1, a2) in double_bond_ends {
        let v1 = normalize_vec([coords[a1][0]-coords[a0][0], coords[a1][1]-coords[a0][1], coords[a1][2]-coords[a0][2]]);
        let v2 = normalize_vec([coords[a1][0]-coords[a2][0], coords[a1][1]-coords[a2][1], coords[a1][2]-coords[a2][2]]);
        if v1.dot(v2) + 1.0 < LINEAR_TOL { return false; }
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
        if !(4..=6).contains(&ring.len()) { continue; }
        let all_aromatic = ring.iter().all(|a| aromatic_atoms.contains(a));
        if all_aromatic {
            let rsize = ring.len();
            for start in 0..rsize {
                ring_torsions.push((ring[start], ring[(start+1)%rsize], ring[(start+2)%rsize], ring[(start+3)%rsize]));
            }
        }
    }

    for &atom_idx in aromatic_atoms.iter() {
        let neighbors = get_neighbors(atom_idx, mol);
        if neighbors.len() >= 3 {
            let k_improper = improper_k_for_atom(atom_idx, mol);
            for i in 0..neighbors.len() {
                for j in (i+1)..neighbors.len() {
                    for k in (j+1)..neighbors.len() {
                        impropers.push((atom_idx, neighbors[i], neighbors[j], neighbors[k], k_improper));
                    }
                }
            }
        }
        let ring_neighbors: Vec<usize> = neighbors.iter().filter(|&&n| aromatic_atoms.contains(&n)).copied().collect();
        let exo_neighbors: Vec<usize> = neighbors.iter().filter(|&&n| !aromatic_atoms.contains(&n)).copied().collect();
        if ring_neighbors.len() >= 2 {
            for &exo in &exo_neighbors {
                exocyclic_torsions.push((exo, atom_idx, ring_neighbors[0], ring_neighbors[1]));
            }
        }
    }

    PlanarityConstraints { impropers, ring_torsions, exocyclic_torsions, aromatic_atoms }
}

fn improper_k_for_atom(atom_idx: usize, mol: &Molecule) -> f64 {
    let base_k = match mol.atoms[atom_idx].symbol.as_str() {
        "C" => 40.0, "N" => 30.0, "O" => 80.0, "S" => 20.0, _ => 40.0,
    };
    base_k * 10.0
}

fn out_of_plane_angle(coords: &[[f64; 3]], central: usize, n1: usize, n2: usize, n3: usize) -> f64 {
    let v1 = [coords[n1][0]-coords[central][0], coords[n1][1]-coords[central][1], coords[n1][2]-coords[central][2]];
    let v2 = [coords[n2][0]-coords[central][0], coords[n2][1]-coords[central][1], coords[n2][2]-coords[central][2]];
    let v3 = [coords[n3][0]-coords[central][0], coords[n3][1]-coords[central][1], coords[n3][2]-coords[central][2]];
    let normal = cross_product([v1[0]-v2[0], v1[1]-v2[1], v1[2]-v2[2]], [v1[0]-v3[0], v1[1]-v3[1], v1[2]-v3[2]]);
    let normal_norm = (normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2]).sqrt();
    if normal_norm < 1e-10 { return 0.0; }
    let v1_norm = (v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2]).sqrt();
    if v1_norm < 1e-10 { return 0.0; }
    let dot = (v1[0]*normal[0] + v1[1]*normal[1] + v1[2]*normal[2]) / normal_norm;
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
        let phi_wrapped = if phi > std::f64::consts::FRAC_PI_2 { phi - std::f64::consts::PI }
            else if phi < -std::f64::consts::FRAC_PI_2 { phi + std::f64::consts::PI } else { phi };
        energy += K_EXO_TOR * phi_wrapped * phi_wrapped;
    }
    energy
}

fn check_planarity(coords: &[[f64; 3]], mol: &Molecule, pc: &PlanarityConstraints, threshold: f64) -> bool {
    let mut visited = HashSet::new();
    let mut components: Vec<Vec<usize>> = Vec::new();
    for &start in &pc.aromatic_atoms {
        if visited.contains(&start) { continue; }
        let mut comp = Vec::new();
        let mut queue = std::collections::VecDeque::new();
        queue.push_back(start);
        visited.insert(start);
        comp.push(start);
        while let Some(idx) = queue.pop_front() {
            for &nbr in &mol.adjacency[idx] {
                if !pc.aromatic_atoms.contains(&nbr) || visited.contains(&nbr) { continue; }
                visited.insert(nbr);
                queue.push_back(nbr);
                comp.push(nbr);
            }
        }
        if comp.len() >= 3 { components.push(comp); }
    }
    for component in &components {
        let cx = component.iter().map(|&i| coords[i][0]).sum::<f64>() / component.len() as f64;
        let cy = component.iter().map(|&i| coords[i][1]).sum::<f64>() / component.len() as f64;
        let cz = component.iter().map(|&i| coords[i][2]).sum::<f64>() / component.len() as f64;
        let mut cov = [[0.0f64; 3]; 3];
        for &idx in component {
            let dx = coords[idx][0] - cx; let dy = coords[idx][1] - cy; let dz = coords[idx][2] - cz;
            cov[0][0] += dx * dx; cov[0][1] += dx * dy; cov[0][2] += dx * dz;
            cov[1][1] += dy * dy; cov[1][2] += dy * dz; cov[2][2] += dz * dz;
        }
        let n = component.len() as f64;
        for row in cov.iter_mut() { for val in row.iter_mut() { *val /= n; } }
        cov[1][0] = cov[0][1]; cov[2][0] = cov[0][2]; cov[2][1] = cov[1][2];
        let normal = eigenvector_smallest_eigenvalue_3x3(&cov);
        for &idx in component {
            let dist = ((coords[idx][0] - cx) * normal[0] + (coords[idx][1] - cy) * normal[1] + (coords[idx][2] - cz) * normal[2]).abs();
            if dist > threshold { return false; }
        }
    }
    true
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
        let expected = compute_bond_length(&mol.atoms[i].symbol, &mol.atoms[j].symbol, bond.bond_type);
        if (actual - expected).abs() > BOND_LENGTH_TOLERANCE {
            return false;
        }
    }
    true
}

fn flatten_aromatic_rings(coords: &mut [[f64; 3]], mol: &Molecule, pc: &PlanarityConstraints) {
    let mut visited = HashSet::new();
    let mut components: Vec<Vec<usize>> = Vec::new();
    for &start in &pc.aromatic_atoms {
        if visited.contains(&start) { continue; }
        let mut comp = Vec::new();
        let mut queue = std::collections::VecDeque::new();
        queue.push_back(start);
        visited.insert(start);
        comp.push(start);
        while let Some(idx) = queue.pop_front() {
            for &nbr in &mol.adjacency[idx] {
                if !pc.aromatic_atoms.contains(&nbr) || visited.contains(&nbr) { continue; }
                visited.insert(nbr);
                queue.push_back(nbr);
                comp.push(nbr);
            }
        }
        if comp.len() >= 3 { components.push(comp); }
    }
    for component in &components {
        let cx = component.iter().map(|&i| coords[i][0]).sum::<f64>() / component.len() as f64;
        let cy = component.iter().map(|&i| coords[i][1]).sum::<f64>() / component.len() as f64;
        let cz = component.iter().map(|&i| coords[i][2]).sum::<f64>() / component.len() as f64;
        let mut cov = [[0.0f64; 3]; 3];
        for &idx in component {
            let dx = coords[idx][0] - cx; let dy = coords[idx][1] - cy; let dz = coords[idx][2] - cz;
            cov[0][0] += dx * dx; cov[0][1] += dx * dy; cov[0][2] += dx * dz;
            cov[1][1] += dy * dy; cov[1][2] += dy * dz; cov[2][2] += dz * dz;
        }
        let n = component.len() as f64;
        for row in cov.iter_mut() { for val in row.iter_mut() { *val /= n; } }
        cov[1][0] = cov[0][1]; cov[2][0] = cov[0][2]; cov[2][1] = cov[1][2];
        let normal = eigenvector_smallest_eigenvalue_3x3(&cov);
        let mut atoms_to_flatten: Vec<usize> = component.clone();
        // First pass: add direct neighbors of aromatic atoms
        for &ring_atom in component {
            for &nbr in &mol.adjacency[ring_atom] {
                if !pc.aromatic_atoms.contains(&nbr) && !atoms_to_flatten.contains(&nbr) { atoms_to_flatten.push(nbr); }
            }
        }
        // Second pass: add H atoms bonded to flattened atoms that should be planar
        // (e.g., H atoms on aniline N, which is conjugated to the aromatic ring)
        let flattened_non_arom: Vec<usize> = atoms_to_flatten.iter().filter(|&&a| !pc.aromatic_atoms.contains(&a)).copied().collect();
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
            let dist = (coords[idx][0] - cx) * normal[0] + (coords[idx][1] - cy) * normal[1] + (coords[idx][2] - cz) * normal[2];
            coords[idx][0] -= dist * normal[0];
            coords[idx][1] -= dist * normal[1];
            coords[idx][2] -= dist * normal[2];
        }
    }
}

pub fn eigenvector_smallest_eigenvalue_3x3(m: &[[f64; 3]; 3]) -> [f64; 3] {
    let c1 = m[0][0] + m[1][1] + m[2][2];
    let c2 = (m[0][0]*m[1][1] - m[0][1].powi(2)) + (m[1][1]*m[2][2] - m[1][2].powi(2)) + (m[0][0]*m[2][2] - m[0][2].powi(2));
    let c3 = m[0][0]*(m[1][1]*m[2][2] - m[1][2]*m[2][1]) - m[0][1]*(m[1][2]*m[2][0] - m[1][0]*m[2][2]) + m[0][2]*(m[1][0]*m[2][1] - m[1][1]*m[2][0]);
    let mut lambda = m[0][0].min(m[1][1]).min(m[2][2]);
    for _ in 0..20 {
        let f = lambda.powi(3) - c1*lambda.powi(2) + c2*lambda - c3;
        let fp = 3.0*lambda.powi(2) - 2.0*c1*lambda + c2;
        if fp.abs() < 1e-12 { break; }
        lambda -= f / fp;
    }
    let mut mat = [[0.0f64; 3]; 3];
    for i in 0..3 { for j in 0..3 { mat[i][j] = m[i][j]; } mat[i][i] -= lambda; }
    let mut v = [mat[0][1]*mat[1][2] - mat[0][2]*mat[1][1], mat[0][2]*mat[1][0] - mat[0][0]*mat[1][2], mat[0][0]*mat[1][1] - mat[0][1]*mat[1][0]];
    let norm = (v[0].powi(2) + v[1].powi(2) + v[2].powi(2)).sqrt();
    if norm > 1e-10 { v[0] /= norm; v[1] /= norm; v[2] /= norm; } else { v = [0.0, 0.0, 1.0]; }
    v
}

// ============================================================================
// Torsion Preferences
// ============================================================================

#[derive(Debug, Clone)]
struct TorsionPreference { i: usize, j: usize, k: usize, l: usize, signs: [i32; 6], v: [f64; 6] }

fn sym(mol: &Molecule, a: usize) -> &str { &mol.atoms[a].symbol }
fn get_hyb(mol: &Molecule, a: usize) -> Hybridization { crate::molecule::graph::determine_hybridization(a, mol) }
fn is_arom(mol: &Molecule, a: usize) -> bool { crate::molecule::graph::is_aromatic(a, mol) }
fn h_count(mol: &Molecule, a: usize) -> usize { mol.adjacency[a].iter().filter(|&&n| mol.atoms[n].symbol == "H").count() }
fn has_db_to_o(mol: &Molecule, c: usize) -> bool {
    mol.adjacency[c].iter().any(|&nbr| {
        mol.atoms[nbr].symbol == "O" && is_double_bond(mol, c, nbr)
    })
}
fn has_db_to_n(mol: &Molecule, c: usize) -> bool {
    mol.adjacency[c].iter().any(|&nbr| {
        mol.atoms[nbr].symbol == "N" && is_double_bond(mol, c, nbr)
    })
}
fn is_carbonyl_c(mol: &Molecule, c: usize) -> bool { sym(mol, c) == "C" && has_db_to_o(mol, c) }
fn is_sulfonyl_s(mol: &Molecule, s: usize) -> bool {
    if sym(mol, s) != "S" { return false; }
    let o_dbls = mol.adjacency[s].iter().filter(|&&nbr| {
        mol.atoms[nbr].symbol == "O" && is_double_bond(mol, s, nbr)
    }).count();
    o_dbls >= 2
}
fn is_nitro_n(mol: &Molecule, n: usize) -> bool {
    if sym(mol, n) != "N" { return false; }
    let o_dbls = mol.adjacency[n].iter().filter(|&&nbr| {
        mol.atoms[nbr].symbol == "O" && is_double_bond(mol, n, nbr)
    }).count();
    o_dbls >= 1 && mol.adjacency[n].iter().any(|&nbr| mol.atoms[nbr].symbol == "O" && !is_double_bond(mol, n, nbr))
}
fn bond_in_ring(rings: &[Vec<usize>], a: usize, b: usize) -> bool {
    rings.iter().any(|r| r.contains(&a) && r.contains(&b))
}
fn ring_size_of_bond(rings: &[Vec<usize>], a: usize, b: usize) -> Option<usize> {
    rings.iter().find(|r| r.contains(&a) && r.contains(&b)).map(|r| r.len())
}

/// Match a torsion pattern for atoms a1-a2-a3-a4.
/// Returns Fourier coefficients (signs, V) if a pattern matches.
/// Patterns are checked in order of specificity — first match wins.
fn match_torsion_pattern(mol: &Molecule, a1: usize, a2: usize, a3: usize, a4: usize, rings: &[Vec<usize>])
    -> Option<([i32; 6], [f64; 6])>
{
    let s1 = sym(mol, a1); let s2 = sym(mol, a2); let s3 = sym(mol, a3); let s4 = sym(mol, a4);
    let h1 = get_hyb(mol, a1); let h2 = get_hyb(mol, a2); let h3 = get_hyb(mol, a3); let h4 = get_hyb(mol, a4);
    let ar1 = is_arom(mol, a1); let ar2 = is_arom(mol, a2); let ar3 = is_arom(mol, a3); let ar4 = is_arom(mol, a4);
    let hc1 = h_count(mol, a1); let hc2 = h_count(mol, a2); let hc3 = h_count(mol, a3); let hc4 = h_count(mol, a4);
    let b23_ring = bond_in_ring(rings, a2, a3);

    // ================================================================
    // 1. AMIDE / PEPTIDE: C(=O)-N  → trans
    // ================================================================
    if is_amide_bond(mol, a2, a3) {
        return Some(([1, -1, 1, 1, 1, 1], [0.0, 10.0, 0.0, 0.0, 0.0, 0.0]));
    }

    // ================================================================
    // 2. ESTER / CARBOXYLIC: O=C-O-*
    // ================================================================
    if is_carbonyl_c(mol, a2) && s3 == "O" && !b23_ring {
        // O=C-O-C (general ester)
        return Some(([1, -1, 1, 1, 1, 1], [0.0, 39.4, 0.0, 0.0, 0.0, 0.0]));
    }
    if s2 == "O" && is_carbonyl_c(mol, a3) && !b23_ring {
        // C-O-C=O (reverse ester)
        return Some(([1, -1, 1, 1, 1, 1], [0.0, 39.4, 0.0, 0.0, 0.0, 0.0]));
    }

    // ================================================================
    // 3. SULFONAMIDE: S(=O)2-N
    // ================================================================
    if is_sulfonyl_s(mol, a2) && s3 == "N" && !b23_ring {
        return Some(([1, -1, 1, 1, 1, 1], [0.0, 16.0, 5.0, 7.0, 0.0, 0.0]));
    }
    if s2 == "N" && is_sulfonyl_s(mol, a3) && !b23_ring {
        return Some(([1, -1, 1, 1, 1, 1], [0.0, 16.0, 5.0, 7.0, 0.0, 0.0]));
    }

    // ================================================================
    // 4. NITRO on aromatic: Ar-NO2
    // ================================================================
    if ar2 && is_nitro_n(mol, a3) && !b23_ring {
        return Some(([1, -1, 1, 1, 1, 1], [0.0, 11.8, 0.0, 0.0, 0.0, 0.0]));
    }
    if is_nitro_n(mol, a2) && ar3 && !b23_ring {
        return Some(([1, -1, 1, 1, 1, 1], [0.0, 11.8, 0.0, 0.0, 0.0, 0.0]));
    }

    // ================================================================
    // 5. AROMATIC-AROMATIC (biaryl) — planar
    // ================================================================
    if ar2 && ar3 && !b23_ring {
        return Some(([1, -1, 1, 1, 1, 1], [0.0, 30.0, 0.0, 0.0, 0.0, 0.0]));
    }

    // ================================================================
    // 6. sp2-sp2 double bond — planar
    // ================================================================
    if is_double_bond(mol, a2, a3) && matches!(h1, Hybridization::Sp2) && matches!(h4, Hybridization::Sp2) && !b23_ring {
        return Some(([1, -1, 1, 1, 1, 1], [0.0, 100.0, 0.0, 0.0, 0.0, 0.0]));
    }

    // ================================================================
    // 7. C=C-C(sp3) and C=C-C(sp2) — specific preferences
    // ================================================================
    if is_double_bond(mol, a2, a3) && matches!(h2, Hybridization::Sp2) && matches!(h3, Hybridization::Sp2) && !b23_ring {
        if matches!(h1, Hybridization::Sp3) && matches!(h4, Hybridization::Sp3) {
            return Some(([1, 0.0_f64 as i32, 1, 1, 1, 1], [0.0, 0.0, 5.0, 0.0, 0.0, 0.0]));
        }
        if matches!(h1, Hybridization::Sp3) || matches!(h4, Hybridization::Sp3) {
            return Some(([1, 0.0_f64 as i32, 1, 1, 1, 1], [0.0, 0.0, 4.0, 0.0, 0.0, 0.0]));
        }
    }

    // ================================================================
    // 8. Aromatic-C(=O) — planar
    // ================================================================
    if (ar2 && is_carbonyl_c(mol, a3)) || (is_carbonyl_c(mol, a2) && ar3) && !b23_ring {
        return Some(([1, -1, 1, 1, 1, 1], [0.0, 8.0, 0.0, 0.0, 0.0, 0.0]));
    }

    // ================================================================
    // 9. Aromatic-C(=N) — planar
    // ================================================================
    if (ar2 && has_db_to_n(mol, a3) && s3 == "C") || (has_db_to_n(mol, a2) && s2 == "C" && ar3) && !b23_ring {
        return Some(([1, -1, 1, 1, 1, 1], [0.0, 15.0, 0.0, 0.0, 0.0, 0.0]));
    }

    // ================================================================
    // 10. Aromatic-O-C / Aromatic-S-C — preferences
    // ================================================================
    if ar2 && s3 == "O" && s4 == "C" && !b23_ring {
        return Some(([1, 0.0_f64 as i32, 1, 1, 1, 1], [0.0, 0.0, 4.0, 0.0, 0.0, 0.0]));
    }
    if s2 == "O" && ar3 && s1 == "C" && !b23_ring {
        return Some(([1, 0.0_f64 as i32, 1, 1, 1, 1], [0.0, 0.0, 4.0, 0.0, 0.0, 0.0]));
    }
    if ar2 && s3 == "S" && s4 == "C" && !b23_ring {
        return Some(([1, -1, 1, 1, 1, 1], [0.0, 39.4, 0.0, 0.0, 0.0, 0.0]));
    }

    // ================================================================
    // 11. C(sp3)-C(=O) — staggered/gauche preferences
    // ================================================================
    if matches!(h2, Hybridization::Sp3) && is_carbonyl_c(mol, a3) && !b23_ring {
        return Some(([1, -1, 1, 1, 1, 1], [0.0, 1.0, 0.0, 0.0, 0.0, 0.0]));
    }
    if is_carbonyl_c(mol, a2) && matches!(h3, Hybridization::Sp3) && !b23_ring {
        return Some(([1, -1, 1, 1, 1, 1], [0.0, 1.0, 0.0, 0.0, 0.0, 0.0]));
    }

    // ================================================================
    // 12. sp3-sp3 single bonds — staggered
    // ================================================================
    if matches!(h2, Hybridization::Sp3) && matches!(h3, Hybridization::Sp3) && !b23_ring {
        // O-C(sp3)-C(sp3)-*
        if s2 == "O" || s3 == "O" {
            return Some(([1, 0.0_f64 as i32, 1, 1, 1, 1], [0.0, 0.0, 2.5, 0.0, 0.0, 0.0]));
        }
        // N-C(sp3)-C(sp3)-*
        if s2 == "N" || s3 == "N" {
            return Some(([1, 0.0_f64 as i32, 1, 1, 1, 1], [0.0, 0.0, 1.0, 0.0, 0.0, 0.0]));
        }
        // General C-C single bond
        return Some(([1, 0.0_f64 as i32, 1, 1, 1, 1], [0.0, 0.0, 4.0, 0.0, 0.0, 0.0]));
    }

    // ================================================================
    // 13. Aromatic-C(sp3) — perpendicular / specific
    // ================================================================
    if (ar2 && matches!(h3, Hybridization::Sp3)) || (matches!(h2, Hybridization::Sp3) && ar3) && !b23_ring {
        return Some(([1, -1, 1, 1, 1, 1], [0.0, 3.0, 0.0, 0.0, 0.0, 0.0]));
    }

    // ================================================================
    // 14. C(sp2)-C(sp3) — gauche preferences
    // ================================================================
    if (matches!(h2, Hybridization::Sp2) && matches!(h3, Hybridization::Sp3)) ||
       (matches!(h2, Hybridization::Sp3) && matches!(h3, Hybridization::Sp2)) && !b23_ring {
        return Some(([1, 0.0_f64 as i32, 1, 1, 1, 1], [0.0, 0.0, 1.9, 0.0, 0.0, 0.0]));
    }

    // ================================================================
    // 15. Small ring torsions (3-8 membered rings)
    // ================================================================
    if b23_ring {
        if let Some(rsize) = ring_size_of_bond(rings, a2, a3) {
            if rsize >= 3 && rsize <= 8 {
                // All-sp2 ring sequences get strong planar preference
                if matches!(h1, Hybridization::Sp2) && matches!(h2, Hybridization::Sp2)
                    && matches!(h3, Hybridization::Sp2) && matches!(h4, Hybridization::Sp2) {
                    return Some(([1, -1, 1, 1, 1, 1], [0.0, 10.0, 0.0, 0.0, 0.0, 0.0]));
                }
                // sp3-sp3 in ring
                if matches!(h2, Hybridization::Sp3) && matches!(h3, Hybridization::Sp3) {
                    return Some(([1, 0.0_f64 as i32, 1, 1, 1, 1], [0.0, 0.0, 5.0, 0.0, 0.0, 0.0]));
                }
                // sp2-sp3 in ring
                if (matches!(h2, Hybridization::Sp2) && matches!(h3, Hybridization::Sp3)) ||
                   (matches!(h2, Hybridization::Sp3) && matches!(h3, Hybridization::Sp2)) {
                    return Some(([1, 0.0_f64 as i32, 1, 1, 1, 1], [0.0, 0.0, 5.0, 0.0, 0.0, 0.0]));
                }
            }
        }
    }

    None
}

fn build_torsion_preferences(mol: &Molecule) -> Vec<TorsionPreference> {
    let rings = crate::molecule::graph::find_rings(mol);
    let mut prefs = Vec::new();
    let mut done_bonds = vec![false; mol.bonds.len()];

    // Iterate all bonds as potential central bonds (like RDKit)
    for (bond_idx, bond) in mol.bonds.iter().enumerate() {
        let a2 = bond.atom1;
        let a3 = bond.atom2;

        // Skip if bond is in >3 rings (RDKit logic)
        let num_ring_bonds = rings.iter().filter(|r| r.contains(&a2) && r.contains(&a3)).count();
        if num_ring_bonds > 3 { continue; }

        // Try all combinations of a1 from a2-neighbors and a4 from a3-neighbors
        for &a1 in &mol.adjacency[a2] {
            if a1 == a3 { continue; }
            for &a4 in &mol.adjacency[a3] {
                if a4 == a2 || a4 == a1 { continue; }

                if let Some((signs, v)) = match_torsion_pattern(mol, a1, a2, a3, a4, &rings) {
                    if !done_bonds[bond_idx] {
                        prefs.push(TorsionPreference { i: a1, j: a2, k: a3, l: a4, signs, v });
                        done_bonds[bond_idx] = true;
                        break;
                    }
                }
            }
            if done_bonds[bond_idx] { break; }
        }
    }

    // Basic knowledge: flat aromatic / sp2 ring torsions (RDKit logic)
    for ring in &rings {
        let rsize = ring.len();
        if rsize < 4 || rsize > 6 { continue; }
        for i in 0..rsize {
            let a1 = ring[i];
            let a2 = ring[(i + 1) % rsize];
            let a3 = ring[(i + 2) % rsize];
            let a4 = ring[(i + 3) % rsize];
            let h1 = get_hyb(mol, a1);
            let h2 = get_hyb(mol, a2);
            let h3 = get_hyb(mol, a3);
            let h4 = get_hyb(mol, a4);
            if matches!(h1, Hybridization::Sp2) && matches!(h2, Hybridization::Sp2)
                && matches!(h3, Hybridization::Sp2) && matches!(h4, Hybridization::Sp2)
            {
                if let Some(bid) = find_bond_index(mol, a2, a3) {
                    if !done_bonds[bid] {
                        prefs.push(TorsionPreference {
                            i: a1, j: a2, k: a3, l: a4,
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

fn etkdg_energy(coords: &[[f64; 3]], bounds: &DistanceBounds, chiral_centers: &[ChiralCenter], tetrahedral: &[ChiralCenter], pc: &PlanarityConstraints, torsion_prefs: &[TorsionPreference], bonds_12: &[(usize, usize)], angles_13: &[(usize, usize, usize)]) -> f64 {
    let mut energy = 0.0;
    for i in 0..bounds.n_atoms { for j in (i + 1)..bounds.n_atoms {
        let dx = coords[i][0] - coords[j][0]; let dy = coords[i][1] - coords[j][1]; let dz = coords[i][2] - coords[j][2];
        let d = (dx*dx + dy*dy + dz*dz).sqrt().max(1e-10);
        let lo = bounds.lower[i][j]; let hi = bounds.upper[i][j];
        if hi >= MAX_UPPER { continue; }
        let range = (hi - lo).max(0.01); let w = 1.0 / range;
        let lo_viol = (lo - d).max(0.0); let hi_viol = (d - hi).max(0.0);
        energy += w * (lo_viol * lo_viol + hi_viol * hi_viol);
    }}
    // 1-2 constraints: lock bond lengths
    const K_12: f64 = 100.0;
    const TOL_12: f64 = 0.01;
    for &(i, j) in bonds_12 {
        let dx = coords[i][0] - coords[j][0]; let dy = coords[i][1] - coords[j][1]; let dz = coords[i][2] - coords[j][2];
        let d = (dx*dx + dy*dy + dz*dz).sqrt();
        let lo = d - TOL_12; let hi = d + TOL_12;
        let lo_viol = (lo - d).max(0.0); let hi_viol = (d - hi).max(0.0);
        energy += K_12 * (lo_viol * lo_viol + hi_viol * hi_viol);
    }
    // 1-3 constraints: lock angle-derived distances
    for &(i, _, k) in angles_13 {
        let dx = coords[i][0] - coords[k][0]; let dy = coords[i][1] - coords[k][1]; let dz = coords[i][2] - coords[k][2];
        let d = (dx*dx + dy*dy + dz*dz).sqrt();
        let lo = d - TOL_12; let hi = d + TOL_12;
        let lo_viol = (lo - d).max(0.0); let hi_viol = (d - hi).max(0.0);
        energy += K_12 * (lo_viol * lo_viol + hi_viol * hi_viol);
    }
    for cc in chiral_centers {
        let vol = chiral_volume(coords, cc.central, cc.neighbors[0], cc.neighbors[1], cc.neighbors[2]);
        if vol < cc.vol_lower { energy += (vol - cc.vol_lower) * (vol - cc.vol_lower); }
        else if vol > cc.vol_upper { energy += (vol - cc.vol_upper) * (vol - cc.vol_upper); }
    }
    for tc in tetrahedral { if !volume_test(coords, tc) { energy += 10.0; } }
    energy += planarity_energy(coords, pc);
    energy += torsion_pref_energy(coords, torsion_prefs);
    energy
}

fn dihedral_gradient_contrib(coords: &[[f64; 3]], i: usize, j: usize, k: usize, l: usize, dedphi: f64) -> ([f64; 3], [f64; 3], [f64; 3], [f64; 3]) {
    let b1 = [coords[j][0]-coords[i][0], coords[j][1]-coords[i][1], coords[j][2]-coords[i][2]];
    let b2 = [coords[k][0]-coords[j][0], coords[k][1]-coords[j][1], coords[k][2]-coords[j][2]];
    let b3 = [coords[l][0]-coords[k][0], coords[l][1]-coords[k][1], coords[l][2]-coords[k][2]];
    let n1 = cross_product(b1, b2);
    let n2 = cross_product(b2, b3);
    let n1_sq = n1[0]*n1[0] + n1[1]*n1[1] + n1[2]*n1[2];
    let n2_sq = n2[0]*n2[0] + n2[1]*n2[1] + n2[2]*n2[2];
    let b2_norm = (b2[0]*b2[0] + b2[1]*b2[1] + b2[2]*b2[2]).sqrt();
    if n1_sq < 1e-20 || n2_sq < 1e-20 || b2_norm < 1e-10 { return ([0.0; 3], [0.0; 3], [0.0; 3], [0.0; 3]); }
    let b2n = [b2[0]/b2_norm, b2[1]/b2_norm, b2[2]/b2_norm];
    let m1 = cross_product(n1, b2n);
    let m2 = cross_product(n2, b2n);
    let g1 = [-m1[0]/n1_sq * dedphi, -m1[1]/n1_sq * dedphi, -m1[2]/n1_sq * dedphi];
    let g4 = [m2[0]/n2_sq * dedphi, m2[1]/n2_sq * dedphi, m2[2]/n2_sq * dedphi];
    let b1_dot_b2 = b1[0]*b2[0] + b1[1]*b2[1] + b1[2]*b2[2];
    let b2_dot_b3 = b2[0]*b3[0] + b2[1]*b3[1] + b2[2]*b3[2];
    let b2_sq = b2_norm * b2_norm;
    let g2 = [
        -g1[0] + (b1_dot_b2 / b2_sq) * g1[0] - (b2_dot_b3 / b2_sq) * g4[0],
        -g1[1] + (b1_dot_b2 / b2_sq) * g1[1] - (b2_dot_b3 / b2_sq) * g4[1],
        -g1[2] + (b1_dot_b2 / b2_sq) * g1[2] - (b2_dot_b3 / b2_sq) * g4[2],
    ];
    let g3 = [
        -g4[0] - (b1_dot_b2 / b2_sq) * g1[0] + (b2_dot_b3 / b2_sq) * g4[0],
        -g4[1] - (b1_dot_b2 / b2_sq) * g1[1] + (b2_dot_b3 / b2_sq) * g4[1],
        -g4[2] - (b1_dot_b2 / b2_sq) * g1[2] + (b2_dot_b3 / b2_sq) * g4[2],
    ];
    (g1, g2, g3, g4)
}

fn etkdg_gradient(coords: &[[f64; 3]], bounds: &DistanceBounds, chiral_centers: &[ChiralCenter], _tetrahedral: &[ChiralCenter], pc: &PlanarityConstraints, torsion_prefs: &[TorsionPreference], _bonds_12: &[(usize, usize)], _angles_13: &[(usize, usize, usize)]) -> Vec<[f64; 3]> {
    let n = coords.len();
    let mut grad = vec![[0.0f64; 3]; n];

    // Distance bounds gradient
    for i in 0..bounds.n_atoms {
        for j in (i + 1)..bounds.n_atoms {
            let dx = coords[i][0] - coords[j][0];
            let dy = coords[i][1] - coords[j][1];
            let dz = coords[i][2] - coords[j][2];
            let d = (dx*dx + dy*dy + dz*dz).sqrt().max(1e-10);
            let lo = bounds.lower[i][j];
            let hi = bounds.upper[i][j];
            if hi >= MAX_UPPER { continue; }
            let range = (hi - lo).max(0.01);
            let w = 1.0 / range;
            let mut dedd = 0.0;
            if d < lo { dedd += -2.0 * w * (lo - d); }
            if d > hi { dedd += 2.0 * w * (d - hi); }
            let gx = dedd * dx / d;
            let gy = dedd * dy / d;
            let gz = dedd * dz / d;
            grad[i][0] += gx; grad[i][1] += gy; grad[i][2] += gz;
            grad[j][0] -= gx; grad[j][1] -= gy; grad[j][2] -= gz;
        }
    }

    // Chiral volume gradient
    for cc in chiral_centers {
        let c = cc.central;
        let n1 = cc.neighbors[0];
        let n2 = cc.neighbors[1];
        let n3 = cc.neighbors[2];
        let v1 = [coords[n1][0]-coords[c][0], coords[n1][1]-coords[c][1], coords[n1][2]-coords[c][2]];
        let v2 = [coords[n2][0]-coords[c][0], coords[n2][1]-coords[c][1], coords[n2][2]-coords[c][2]];
        let v3 = [coords[n3][0]-coords[c][0], coords[n3][1]-coords[c][1], coords[n3][2]-coords[c][2]];
        let vol = v1[0]*(v2[1]*v3[2]-v2[2]*v3[1]) + v1[1]*(v2[2]*v3[0]-v2[0]*v3[2]) + v1[2]*(v2[0]*v3[1]-v2[1]*v3[0]);
        let mut dvol = 0.0;
        if vol < cc.vol_lower { dvol = 2.0 * (vol - cc.vol_lower); }
        else if vol > cc.vol_upper { dvol = 2.0 * (vol - cc.vol_upper); }
        if dvol != 0.0 {
            let g_n1 = cross_product(v2, v3);
            let g_n2 = cross_product(v3, v1);
            let g_n3 = cross_product(v1, v2);
            let g_c = [-(g_n1[0]+g_n2[0]+g_n3[0]), -(g_n1[1]+g_n2[1]+g_n3[1]), -(g_n1[2]+g_n2[2]+g_n3[2])];
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
            grad[i][dim] += g1[dim]; grad[j][dim] += g2[dim];
            grad[k][dim] += g3[dim]; grad[l][dim] += g4[dim];
        }
    }

    // Planarity gradient: exocyclic torsions
    const K_EXO_TOR: f64 = 2.0;
    for &(i, j, k, l) in &pc.exocyclic_torsions {
        let phi = dihedral_angle(coords, i, j, k, l);
        let phi_wrapped = if phi > std::f64::consts::FRAC_PI_2 { phi - std::f64::consts::PI }
            else if phi < -std::f64::consts::FRAC_PI_2 { phi + std::f64::consts::PI } else { phi };
        let dedphi = 2.0 * K_EXO_TOR * phi_wrapped;
        let (g1, g2, g3, g4) = dihedral_gradient_contrib(coords, i, j, k, l, dedphi);
        for dim in 0..3 {
            grad[i][dim] += g1[dim]; grad[j][dim] += g2[dim];
            grad[k][dim] += g3[dim]; grad[l][dim] += g4[dim];
        }
    }

    // Planarity gradient: impropers (small-angle approx)
    for &(central, n1, n2, n3, k_imp) in &pc.impropers {
        let v1 = [coords[n1][0]-coords[central][0], coords[n1][1]-coords[central][1], coords[n1][2]-coords[central][2]];
        let v2 = [coords[n2][0]-coords[central][0], coords[n2][1]-coords[central][1], coords[n2][2]-coords[central][2]];
        let v3 = [coords[n3][0]-coords[central][0], coords[n3][1]-coords[central][1], coords[n3][2]-coords[central][2]];
        let a = [v1[0]-v2[0], v1[1]-v2[1], v1[2]-v2[2]];
        let b = [v1[0]-v3[0], v1[1]-v3[1], v1[2]-v3[2]];
        let normal = cross_product(a, b);
        let n_norm = (normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2]).sqrt();
        let v1_norm = (v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2]).sqrt();
        if n_norm < 1e-10 || v1_norm < 1e-10 { continue; }
        let sin_oop = (normal[0]*v1[0] + normal[1]*v1[1] + normal[2]*v1[2]) / (n_norm * v1_norm);
        let chi = sin_oop.clamp(-1.0, 1.0).asin();
        let dchi = 2.0 * k_imp * chi / (1.0 - sin_oop*sin_oop).max(1e-10).sqrt();
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
                grad[atom][dim] += k_imp * (ep*ep - em*em) / eps;
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
            grad[p.i][dim] += g1[dim]; grad[p.j][dim] += g2[dim];
            grad[p.k][dim] += g3[dim]; grad[p.l][dim] += g4[dim];
        }
    }

    grad
}

fn minimize_etkdg(coords: &mut [[f64; 3]], bounds: &DistanceBounds, chiral_centers: &[ChiralCenter], tetrahedral: &[ChiralCenter], pc: &PlanarityConstraints, torsion_prefs: &[TorsionPreference], bonds_12: &[(usize, usize)], angles_13: &[(usize, usize, usize)], max_iter: usize, force_tol: f64) -> f64 {
    let n = coords.len();
    if n == 0 { return 0.0; }
    let mut best_coords = coords.to_vec();
    let mut best_energy = etkdg_energy(coords, bounds, chiral_centers, tetrahedral, pc, torsion_prefs, bonds_12, angles_13);
    for _ in 0..max_iter {
        let grad = etkdg_gradient(coords, bounds, chiral_centers, tetrahedral, pc, torsion_prefs, bonds_12, angles_13);
        let max_g = grad.iter().map(|g| (g[0]*g[0] + g[1]*g[1] + g[2]*g[2]).sqrt()).fold(0.0f64, f64::max);
        if max_g < force_tol { break; }
        let step = 0.1 / max_g.max(1e-10);
        for i in 0..n { coords[i][0] -= step * grad[i][0]; coords[i][1] -= step * grad[i][1]; coords[i][2] -= step * grad[i][2]; }
        let energy = etkdg_energy(coords, bounds, chiral_centers, tetrahedral, pc, torsion_prefs, bonds_12, angles_13);
        if energy < best_energy { best_energy = energy; best_coords = coords.to_vec(); }
    }
    coords.copy_from_slice(&best_coords);
    best_energy
}

fn bounds_fulfilled(coords: &[[f64; 3]], bounds: &DistanceBounds, atom_indices: &[usize]) -> bool {
    for i in 0..atom_indices.len() { for j in (i + 1)..atom_indices.len() {
        let a1 = atom_indices[i]; let a2 = atom_indices[j];
        let dx = coords[a1][0] - coords[a2][0]; let dy = coords[a1][1] - coords[a2][1]; let dz = coords[a1][2] - coords[a2][2];
        let d = (dx*dx + dy*dy + dz*dz).sqrt();
        let lb = bounds.get_lower(a1, a2); let ub = bounds.get_upper(a1, a2);
        if ((d < lb) && (d - lb).abs() > 0.1 * ub) || ((d > ub) && (d - ub).abs() > 0.1 * ub) { return false; }
    }}
    true
}

fn collect_bonds(mol: &Molecule) -> Vec<(usize, usize)> {
    mol.bonds.iter().map(|b| (b.atom1, b.atom2)).collect()
}

fn collect_angles(mol: &Molecule) -> Vec<(usize, usize, usize)> {
    let mut angles = Vec::new();
    let n_bonds = mol.bonds.len();
    for i in 0..n_bonds {
        for j in (i+1)..n_bonds {
            let b1 = &mol.bonds[i];
            let b2 = &mol.bonds[j];
            let (a11, a12, a21, a22) = (b1.atom1, b1.atom2, b2.atom1, b2.atom2);
            let mut found = None;
            if a12 == a21 { found = Some((a11, a12, a22)); }
            else if a12 == a22 { found = Some((a11, a12, a21)); }
            else if a11 == a21 { found = Some((a12, a11, a22)); }
            else if a11 == a22 { found = Some((a12, a11, a21)); }
            if let Some(angle) = found {
                angles.push(angle);
            }
        }
    }
    angles
}

// ============================================================================
// Main Public API
// ============================================================================

pub fn generate_initial_coords(mol: &Molecule) -> Vec<[f64; 3]> {
    let config = ETKDGConfig::default();
    generate_initial_coords_with_config(mol, &config)
}

pub fn generate_initial_coords_with_config(mol: &Molecule, config: &ETKDGConfig) -> Vec<[f64; 3]> {
    if mol.atoms.is_empty() { return Vec::new(); }

    let mut rng = if config.random_seed >= 0 {
        Rng::new(config.random_seed as u64)
    } else {
        let seed = std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).unwrap_or_default().as_nanos() as u64;
        Rng::new(seed)
    };

    let mut bounds = build_distance_bounds(mol, config);
    bounds.smooth_triangle_inequality();

    let pc = build_planarity_constraints(mol);
    let (chiral_centers, tetrahedral) = find_chiral_centers(mol);
    let (double_bond_ends, stereo_dbs) = find_stereo_double_bonds(mol);
    let torsion_prefs = build_torsion_preferences(mol);
    let bonds_12 = collect_bonds(mol);
    let angles_13 = collect_angles(mol);

    let mut best_coords = Vec::new();
    let mut best_energy = f64::INFINITY;
    let max_attempts = config.max_attempts.max(10);

    for _attempt in 0..max_attempts {
        let mut coords_4d = generate_initial_coords_from_bounds(&bounds, &mut rng);
        minimize_fourth_dimension(&mut coords_4d, &bounds, 200);
        let mut coords_3d: Vec<[f64; 3]> = coords_4d.iter().map(|c| [c[0], c[1], c[2]]).collect();

        let energy = minimize_etkdg(&mut coords_3d, &bounds, &chiral_centers, &tetrahedral, &pc, &torsion_prefs, &bonds_12, &angles_13, 400, 1e-3);

        if !check_tetrahedral(&coords_3d, &tetrahedral) {
            if energy < best_energy { best_energy = energy; best_coords = coords_3d; }
            continue;
        }

        if !chiral_centers.is_empty() && !check_chiral_centers(&coords_3d, &chiral_centers) {
            if energy < best_energy { best_energy = energy; best_coords = coords_3d; }
            continue;
        }

        if !pc.aromatic_atoms.is_empty() {
            flatten_aromatic_rings(&mut coords_3d, mol, &pc);
        }

        let energy = minimize_etkdg(&mut coords_3d, &bounds, &chiral_centers, &tetrahedral, &pc, &torsion_prefs, &bonds_12, &angles_13, 300, 1e-3);

        let planar = check_planarity(&coords_3d, mol, &pc, 0.1);
        let db_geom_ok = double_bond_geometry_checks(&coords_3d, &double_bond_ends);

        let mut chiral_ok = true;
        if !chiral_centers.is_empty() {
            chiral_ok = check_chiral_centers(&coords_3d, &chiral_centers);
            if chiral_ok {
                let mut atoms = std::collections::HashSet::new();
                for cc in &chiral_centers {
                    if cc.central != cc.neighbors[3] {
                        atoms.insert(cc.central); atoms.insert(cc.neighbors[0]);
                        atoms.insert(cc.neighbors[1]); atoms.insert(cc.neighbors[2]); atoms.insert(cc.neighbors[3]);
                    }
                }
                let atoms_vec: Vec<usize> = atoms.into_iter().collect();
                if !atoms_vec.is_empty() && !bounds_fulfilled(&coords_3d, &bounds, &atoms_vec) { chiral_ok = false; }
            }
            if chiral_ok {
                for cc in &chiral_centers { if !center_in_volume(&coords_3d, cc) { chiral_ok = false; break; } }
            }
        }

        let db_stereo_ok = stereo_dbs.is_empty() || check_double_bond_stereo(&coords_3d, &stereo_dbs);
        let no_clash = !has_vdw_clash(&coords_3d, mol);
        let bonds_ok = bond_lengths_reasonable(&coords_3d, mol);

        if planar && db_geom_ok && chiral_ok && db_stereo_ok && no_clash && bonds_ok {
            let e_per_atom = energy / coords_3d.len() as f64;
            if e_per_atom < MAX_MINIMIZED_E_PER_ATOM { return coords_3d; }
        }

        if energy < best_energy { best_energy = energy; best_coords = coords_3d; }
    }

    if !best_coords.is_empty() {
        best_coords
    } else {
        let coords_4d = generate_initial_coords_from_bounds(&bounds, &mut rng);
        coords_4d.iter().map(|c| [c[0], c[1], c[2]]).collect()
    }
}
