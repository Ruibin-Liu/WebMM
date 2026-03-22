//! ETKDG v3 3D coordinate embedding

use crate::molecule::{BondType, Molecule};

fn random_f64() -> f64 {
    let mut buf = [0u8; 8];
    let _ = getrandom::getrandom(&mut buf);
    let val = u64::from_ne_bytes(buf);
    val as f64 / u64::MAX as f64
}

/// ETKDG v3 configuration
#[derive(Debug, Clone)]
pub struct ETKDGConfig {
    /// Number of random conformations to generate
    pub max_attempts: usize,
    /// Convergence threshold for DGFF
    pub convergence_threshold: f64,
    /// Maximum iterations for DGFF minimization
    pub max_iterations: usize,
    /// Scale factor for van der Waals radii
    pub vdw_scale: f64,
}

impl Default for ETKDGConfig {
    fn default() -> Self {
        Self {
            max_attempts: 10,
            convergence_threshold: 1e-6,
            max_iterations: 200,
            vdw_scale: 0.8,
        }
    }
}

/// Distance bounds matrix
#[derive(Debug, Clone)]
pub struct DistanceBounds {
    /// Lower bounds for each atom pair
    pub lower: Vec<Vec<f64>>,
    /// Upper bounds for each atom pair
    pub upper: Vec<Vec<f64>>,
    /// Number of atoms
    pub n_atoms: usize,
}

impl DistanceBounds {
    pub fn new(n_atoms: usize) -> Self {
        let lower = vec![vec![0.0; n_atoms]; n_atoms];
        let upper = vec![vec![f64::INFINITY; n_atoms]; n_atoms];
        Self {
            lower,
            upper,
            n_atoms,
        }
    }

    /// Apply triangle inequality smoothing
    pub fn smooth_triangle_inequality(&mut self) {
        let mut changed = true;
        let epsilon = 1e-6;

        while changed {
            changed = false;
            for k in 0..self.n_atoms {
                for i in 0..self.n_atoms {
                    for j in i + 1..self.n_atoms {
                        // Lower bound: |L(i,k) - L(k,j)| <= L(i,j)
                        let new_lower = (self.lower[i][k] - self.lower[k][j]).abs();
                        if new_lower > self.lower[i][j] + epsilon {
                            self.lower[i][j] = new_lower;
                            self.lower[j][i] = new_lower;
                            changed = true;
                        }

                        // Upper bound: U(i,j) <= U(i,k) + U(k,j)
                        let new_upper = self.upper[i][k] + self.upper[k][j];
                        if new_upper < self.upper[i][j] - epsilon {
                            self.upper[i][j] = new_upper;
                            self.upper[j][i] = new_upper;
                            changed = true;
                        }

                        // Consistency check: L(i,j) <= U(i,j)
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
}

/// Van der Waals radii for common elements (Angstroms)
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
        _ => 1.70,
    }
}

/// Covalent radii for common elements (Angstroms)
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
        _ => 0.76,
    }
}

/// Build distance bounds matrix from molecular topology
fn build_distance_bounds(mol: &Molecule, config: &ETKDGConfig) -> DistanceBounds {
    let n_atoms = mol.atoms.len();
    let mut bounds = DistanceBounds::new(n_atoms);

    // Initialize with reasonable bounds
    for i in 0..n_atoms {
        for j in 0..n_atoms {
            if i == j {
                bounds.lower[i][j] = 0.0;
                bounds.upper[i][j] = 0.0;
            } else {
                let r_i = vdw_radius(&mol.atoms[i].symbol);
                let r_j = vdw_radius(&mol.atoms[j].symbol);

                // Lower bound: sum of covalent radii
                let r_cov_i = covalent_radius(&mol.atoms[i].symbol);
                let r_cov_j = covalent_radius(&mol.atoms[j].symbol);
                bounds.lower[i][j] = r_cov_i + r_cov_j - 0.2;

                // Upper bound: sum of van der Waals radii
                bounds.upper[i][j] = (r_i + r_j) * config.vdw_scale;
            }
        }
    }

    // Apply bond constraints
    for bond in &mol.bonds {
        let i = bond.atom1;
        let j = bond.atom2;

        let r_cov_i = covalent_radius(&mol.atoms[i].symbol);
        let r_cov_j = covalent_radius(&mol.atoms[j].symbol);
        let ideal_length = r_cov_i + r_cov_j;

        match bond.bond_type {
            BondType::Single => {
                bounds.lower[i][j] = ideal_length - 0.1;
                bounds.upper[i][j] = ideal_length + 0.1;
            }
            BondType::Double => {
                bounds.lower[i][j] = ideal_length - 0.08;
                bounds.upper[i][j] = ideal_length + 0.08;
            }
            BondType::Triple => {
                bounds.lower[i][j] = ideal_length - 0.06;
                bounds.upper[i][j] = ideal_length + 0.06;
            }
            BondType::Aromatic => {
                bounds.lower[i][j] = ideal_length - 0.09;
                bounds.upper[i][j] = ideal_length + 0.09;
            }
        }

        bounds.lower[j][i] = bounds.lower[i][j];
        bounds.upper[j][i] = bounds.upper[i][j];
    }

    // 1-3 bounds (angle-derived): law of cosines
    let angle_info = crate::molecule::graph::find_angles(mol);
    for angle in &angle_info {
        let (a_min, a_max) = (angle.atom1.min(angle.atom2), angle.atom1.max(angle.atom2));
        let (b_min, b_max) = (angle.atom2.min(angle.atom3), angle.atom2.max(angle.atom3));
        let r_ij = bounds.lower[a_min][a_max];
        let r_jk = bounds.lower[b_min][b_max];

        // Default theta0 = 109.47 degrees (tetrahedral)
        let theta0_rad = 109.47_f64.to_radians();
        let r_ik_sq = r_ij * r_ij + r_jk * r_jk - 2.0 * r_ij * r_jk * theta0_rad.cos();
        if r_ik_sq > 0.0 {
            let r_ik = r_ik_sq.sqrt();
            let tolerance = 0.2;
            let lower = (r_ik - tolerance).max(0.5);
            let upper = r_ik + tolerance;
            let (a, c) = (angle.atom1.min(angle.atom3), angle.atom1.max(angle.atom3));
            bounds.lower[a][c] = bounds.lower[a][c].max(lower);
            bounds.upper[a][c] = bounds.upper[a][c].min(upper);
            bounds.lower[c][a] = bounds.lower[a][c];
            bounds.upper[c][a] = bounds.upper[a][c];
        }
    }

    // 1-4 bounds (torsion-derived)
    let torsion_info = crate::molecule::graph::find_torsions(mol);
    for torsion in &torsion_info {
        let r_ij = (bounds.lower[torsion.atom1.min(torsion.atom2)]
            [torsion.atom1.max(torsion.atom2)]
            + bounds.upper[torsion.atom1.min(torsion.atom2)][torsion.atom1.max(torsion.atom2)])
            / 2.0;
        let r_jk = (bounds.lower[torsion.atom2.min(torsion.atom3)]
            [torsion.atom2.max(torsion.atom3)]
            + bounds.upper[torsion.atom2.min(torsion.atom3)][torsion.atom2.max(torsion.atom3)])
            / 2.0;
        let r_kl = (bounds.lower[torsion.atom3.min(torsion.atom4)]
            [torsion.atom3.max(torsion.atom4)]
            + bounds.upper[torsion.atom3.min(torsion.atom4)][torsion.atom3.max(torsion.atom4)])
            / 2.0;

        let upper = r_ij + r_jk + r_kl + 0.1;
        let inner = r_ij * r_ij + r_jk * r_jk + r_kl * r_kl - 2.0 * r_ij * r_jk - 2.0 * r_jk * r_kl
            + 2.0 * r_ij * r_kl;
        let lower = if inner > 0.0 { inner.sqrt() - 0.1 } else { 0.5 };

        let (a, d) = (
            torsion.atom1.min(torsion.atom4),
            torsion.atom1.max(torsion.atom4),
        );
        bounds.lower[a][d] = bounds.lower[a][d].max(lower);
        bounds.upper[a][d] = bounds.upper[a][d].min(upper);
        bounds.lower[d][a] = bounds.lower[a][d];
        bounds.upper[d][a] = bounds.upper[a][d];
    }

    // Ring closure: tighten upper bounds for ring bonds
    let rings = crate::molecule::graph::find_rings(mol);
    for ring in &rings {
        let total_upper: f64 = (0..ring.len())
            .map(|i| {
                let a = ring[i];
                let b = ring[(i + 1) % ring.len()];
                bounds.upper[a.min(b)][a.max(b)]
            })
            .sum();
        let per_atom = total_upper / ring.len() as f64;
        for w in 0..ring.len() {
            let a = ring[w];
            let b = ring[(w + 1) % ring.len()];
            let (a2, b2) = (a.min(b), a.max(b));
            bounds.upper[a2][b2] = bounds.upper[a2][b2].min(per_atom + 0.3);
            bounds.upper[b2][a2] = bounds.upper[a2][b2];
        }
    }

    bounds
}

/// Generate initial 4D coordinates using metric matrix approach
fn generate_4d_coordinates(bounds: &DistanceBounds) -> Vec<[f64; 4]> {
    let n_atoms = bounds.n_atoms;

    // Start with small random coordinates
    let mut coords = Vec::with_capacity(n_atoms);
    for _ in 0..n_atoms {
        coords.push([
            (random_f64() - 0.5) * 0.1,
            (random_f64() - 0.5) * 0.1,
            (random_f64() - 0.5) * 0.1,
            (random_f64() - 0.5) * 0.1,
        ]);
    }

    // Simple embedding using stochastic proximity embedding
    for iteration in 0..100 {
        let learning_rate = 0.1 * (1.0 - iteration as f64 / 100.0);

        for i in 0..n_atoms {
            for j in i + 1..n_atoms {
                let dx = coords[i][0] - coords[j][0];
                let dy = coords[i][1] - coords[j][1];
                let dz = coords[i][2] - coords[j][2];
                let dw = coords[i][3] - coords[j][3];

                let current_dist = (dx * dx + dy * dy + dz * dz + dw * dw).sqrt();

                let target_dist = if current_dist < bounds.lower[i][j] {
                    bounds.lower[i][j]
                } else if current_dist > bounds.upper[i][j] {
                    bounds.upper[i][j].min(10.0)
                } else {
                    continue;
                };

                if current_dist >= 1e-8 {
                    let scale = learning_rate * (target_dist - current_dist) / current_dist;
                    let dx_scale = dx * scale;
                    let dy_scale = dy * scale;
                    let dz_scale = dz * scale;
                    let dw_scale = dw * scale;

                    coords[i][0] += dx_scale;
                    coords[i][1] += dy_scale;
                    coords[i][2] += dz_scale;
                    coords[i][3] += dw_scale;

                    coords[j][0] -= dx_scale;
                    coords[j][1] -= dy_scale;
                    coords[j][2] -= dz_scale;
                    coords[j][3] -= dw_scale;
                }
            }
        }
    }

    coords
}

/// Project 4D coordinates to 3D using eigenvector projection
fn project_to_3d(coords_4d: &[[f64; 4]]) -> Vec<[f64; 3]> {
    let n = coords_4d.len();
    if n < 3 {
        return coords_4d.iter().map(|c| [c[0], c[1], c[2]]).collect();
    }

    let centroid: [f64; 4] = {
        let mut c = [0.0; 4];
        for coord in coords_4d {
            for d in 0..4 {
                c[d] += coord[d];
            }
        }
        for val in c.iter_mut() {
            *val /= n as f64;
        }
        c
    };

    let centered: Vec<[f64; 4]> = coords_4d
        .iter()
        .map(|c| {
            [
                c[0] - centroid[0],
                c[1] - centroid[1],
                c[2] - centroid[2],
                c[3] - centroid[3],
            ]
        })
        .collect();

    let mut gram = vec![vec![0.0f64; n]; n];
    for i in 0..n {
        for j in 0..n {
            let dx = centered[i][0] - centered[j][0];
            let dy = centered[i][1] - centered[j][1];
            let dz = centered[i][2] - centered[j][2];
            let dw = centered[i][3] - centered[j][3];
            gram[i][j] = dx * dx + dy * dy + dz * dz + dw * dw;
        }
    }

    let mut b = vec![vec![0.0f64; n]; n];
    for i in 0..n {
        for j in 0..n {
            b[i][j] = -0.5 * (gram[i][j] - gram[i][0] - gram[0][j] + gram[0][0]);
        }
    }

    let mut eigenvectors = vec![vec![0.0; n]; 3];
    let mut eigenvalues = [0.0f64; 3];

    for k in 0..3 {
        let mut v = vec![0.0f64; n];
        for vi in v.iter_mut() {
            *vi = random_f64() - 0.5;
        }
        let v_norm: f64 = v.iter().map(|x| x * x).sum::<f64>().sqrt();
        for vi in v.iter_mut() {
            *vi /= v_norm;
        }

        for _ in 0..200 {
            let mut new_v = vec![0.0f64; n];
            for i in 0..n {
                for j in 0..n {
                    new_v[i] += b[i][j] * v[j];
                }
            }

            for prev_eigvec in eigenvectors.iter().take(k) {
                let dot: f64 = new_v
                    .iter()
                    .zip(prev_eigvec.iter())
                    .map(|(a, b)| a * b)
                    .sum();
                for (new_vi, prev_ei) in new_v.iter_mut().zip(prev_eigvec.iter()) {
                    *new_vi -= dot * prev_ei;
                }
            }

            let norm: f64 = new_v.iter().map(|x| x * x).sum::<f64>().sqrt();
            if norm < 1e-10 {
                break;
            }
            for vi in new_v.iter_mut() {
                *vi /= norm;
            }

            v = new_v;
        }

        let mut lambda = 0.0;
        for i in 0..n {
            for j in 0..n {
                lambda += v[i] * b[i][j] * v[j];
            }
        }

        eigenvalues[k] = lambda;
        eigenvectors[k] = v;
    }

    let mut coords_3d = vec![[0.0; 3]; n];
    for i in 0..n {
        for d in 0..3 {
            let ev = eigenvalues[d];
            coords_3d[i][d] = if ev > 0.0 {
                eigenvectors[d][i] * ev.sqrt()
            } else {
                eigenvectors[d][i] * 0.001
            };
        }
    }

    coords_3d
}

/// Refine coordinates using actual MMFF94 force field + L-BFGS
fn refine_with_ff(mol: &Molecule, coords: &mut [[f64; 3]], config: &ETKDGConfig) {
    let variant = crate::MMFFVariant::MMFF94s;
    let ff = crate::mmff::MMFFForceField::new(mol, variant);

    let conv = crate::ConvergenceOptions {
        max_force: 0.05,
        rms_force: 0.005,
        energy_change: 1e-6,
        max_iterations: config.max_iterations,
    };

    let result = crate::optimizer::optimize(&ff, coords, &conv);

    for (i, coord) in result.optimized_coords.iter().enumerate() {
        coords[i] = *coord;
    }
}

/// Generate initial 3D coordinates using ETKDG v3
pub fn generate_initial_coords(mol: &Molecule) -> Vec<[f64; 3]> {
    let config = ETKDGConfig::default();
    generate_initial_coords_with_config(mol, &config)
}

/// Generate initial 3D coordinates using ETKDG v3 with custom configuration
pub fn generate_initial_coords_with_config(mol: &Molecule, config: &ETKDGConfig) -> Vec<[f64; 3]> {
    if mol.atoms.is_empty() {
        return Vec::new();
    }

    let mut bounds = build_distance_bounds(mol, config);
    bounds.smooth_triangle_inequality();

    let mut best_coords = Vec::new();
    let mut best_energy = f64::INFINITY;

    for _attempt in 0..config.max_attempts {
        let coords_4d = generate_4d_coordinates(&bounds);
        let mut coords_3d = project_to_3d(&coords_4d);
        refine_with_ff(mol, &mut coords_3d, config);

        let variant = crate::MMFFVariant::MMFF94s;
        let ff = crate::mmff::MMFFForceField::new(mol, variant);
        let energy = ff.calculate_energy(&coords_3d);

        if energy < best_energy {
            best_energy = energy;
            best_coords = coords_3d;
        }
    }

    best_coords
}
