//! ETKDG v3 3D coordinate embedding

use crate::molecule::{Bond, BondType, Molecule};

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
            max_attempts: 50,
            convergence_threshold: 1e-6,
            max_iterations: 1000,
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
        _ => 1.70, // Default for unknown elements
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
        _ => 0.76, // Default for unknown elements
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
                bounds.lower[i][j] = r_cov_i + r_cov_j - 0.2; // Allow slight compression

                // Upper bound: sum of van der Waals radii
                bounds.upper[i][j] = (r_i + r_j) * config.vdw_scale;
            }
        }
    }

    // Apply bond constraints
    for bond in &mol.bonds {
        let i = bond.atom1 as usize;
        let j = bond.atom2 as usize;

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

    bounds
}

/// Generate initial 4D coordinates using metric matrix approach
fn generate_4d_coordinates(bounds: &DistanceBounds) -> Vec<[f64; 4]> {
    let n_atoms = bounds.n_atoms;

    // Start with small random coordinates
    let mut coords = Vec::with_capacity(n_atoms);
    for i in 0..n_atoms {
        coords.push([
            (js_sys::Math::random() - 0.5) * 0.1,
            (js_sys::Math::random() - 0.5) * 0.1,
            (js_sys::Math::random() - 0.5) * 0.1,
            (js_sys::Math::random() - 0.5) * 0.1,
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
                    bounds.upper[i][j].min(10.0) // Cap to reasonable value
                } else {
                    continue; // Within bounds
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

/// Project 4D coordinates to 3D
fn project_to_3d(coords_4d: &[[f64; 4]]) -> Vec<[f64; 3]> {
    let n_atoms = coords_4d.len();
    let mut coords_3d = Vec::with_capacity(n_atoms);

    for coord in coords_4d {
        coords_3d.push([coord[0], coord[1], coord[2]]);
    }

    coords_3d
}

/// Distance Geometry Force Field (DGFF) minimization
fn dgff_minimization(
    mol: &Molecule,
    coords: &mut [[f64; 3]],
    bounds: &DistanceBounds,
    config: &ETKDGConfig,
) {
    let n_atoms = coords.len();
    let k_bond = 100.0; // Bond force constant (reduced for stability)
    let k_vdw = 10.0; // Van der Waals force constant (reduced for stability)

    for iteration in 0..config.max_iterations {
        let mut total_gradient = 0.0;
        let mut forces = vec![[0.0; 3]; n_atoms];

        // Calculate forces from distance constraints
        for i in 0..n_atoms {
            for j in i + 1..n_atoms {
                let dx = coords[j][0] - coords[i][0];
                let dy = coords[j][1] - coords[i][1];
                let dz = coords[j][2] - coords[i][2];
                let dist_sq = dx * dx + dy * dy + dz * dz;
                let dist = dist_sq.sqrt();

                if dist < 1e-8 {
                    continue;
                }

                let mut force = 0.0;

                // Van der Waals repulsion
                if dist < bounds.upper[i][j] {
                    force += k_vdw * (bounds.upper[i][j] - dist) / dist;
                }

                // Bond constraints (stronger force)
                for bond in &mol.bonds {
                    if (bond.atom1 as usize == i && bond.atom2 as usize == j)
                        || (bond.atom1 as usize == j && bond.atom2 as usize == i)
                    {
                        let ideal_dist = (bounds.lower[i][j] + bounds.upper[i][j]) / 2.0;
                        force += k_bond * (ideal_dist - dist) / dist;
                        break;
                    }
                }

                let fx = force * dx;
                let fy = force * dy;
                let fz = force * dz;

                forces[i][0] -= fx;
                forces[i][1] -= fy;
                forces[i][2] -= fz;

                forces[j][0] += fx;
                forces[j][1] += fy;
                forces[j][2] += fz;

                total_gradient += fx * fx + fy * fy + fz * fz;
            }
        }

        // Update coordinates (simple gradient descent with damping)
        let step_size = 0.01;
        for i in 0..n_atoms {
            // Limit force magnitude to prevent instability
            let fx = forces[i][0].max(-1.0).min(1.0);
            let fy = forces[i][1].max(-1.0).min(1.0);
            let fz = forces[i][2].max(-1.0).min(1.0);

            coords[i][0] += step_size * fx;
            coords[i][1] += step_size * fy;
            coords[i][2] += step_size * fz;
        }

        // Check convergence
        if total_gradient.sqrt() < config.convergence_threshold {
            break;
        }
    }
}

/// Generate initial 3D coordinates using ETKDG v3
pub fn generate_initial_coords(mol: &Molecule) -> Vec<[f64; 3]> {
    let config = ETKDGConfig::default();
    generate_initial_coords_with_config(mol, &config)
}

/// Generate initial 3D coordinates using ETKDG v3 with custom configuration
pub fn generate_initial_coords_with_config(mol: &Molecule, config: &ETKDGConfig) -> Vec<[f64; 3]> {
    if mol.atoms.len() < 1 {
        return Vec::new();
    }

    // Step 1: Build distance bounds matrix
    let mut bounds = build_distance_bounds(mol, config);

    // Step 2: Triangle inequality smoothing
    bounds.smooth_triangle_inequality();

    // Step 3: Generate 4D coordinates
    let coords_4d = generate_4d_coordinates(&bounds);

    // Step 4: Project to 3D
    let mut coords_3d = project_to_3d(&coords_4d);

    // Step 5: DGFF minimization
    dgff_minimization(mol, &mut coords_3d, &bounds, config);

    coords_3d
}
