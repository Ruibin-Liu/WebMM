//! ETKDG v3 3D coordinate embedding

use crate::molecule::{BondType, Hybridization, Molecule};

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

/// Estimate the bond angle at the central atom using hybridization
fn estimate_angle(a: usize, b: usize, c: usize, mol: &Molecule, bounds: &DistanceBounds) -> f64 {
    let hyb = crate::molecule::graph::determine_hybridization(b, mol);
    match hyb {
        Hybridization::Sp1 => std::f64::consts::PI,
        Hybridization::Sp2 => 120.0_f64.to_radians(),
        Hybridization::Sp3 => 109.47_f64.to_radians(),
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

                let r_cov_i = covalent_radius(&mol.atoms[i].symbol);
                let r_cov_j = covalent_radius(&mol.atoms[j].symbol);
                bounds.lower[i][j] = r_cov_i + r_cov_j - 0.3;

                // Upper bound: generous VdW sum to accommodate non-bonded contacts
                bounds.upper[i][j] = r_i + r_j + 1.0;
            }
        }
    }

    // Apply bond constraints
    for bond in &mol.bonds {
        let i = bond.atom1;
        let j = bond.atom2;

        let r_cov_i = covalent_radius(&mol.atoms[i].symbol);
        let r_cov_j = covalent_radius(&mol.atoms[j].symbol);
        let single_length = r_cov_i + r_cov_j;

        match bond.bond_type {
            BondType::Single => {
                bounds.lower[i][j] = single_length - 0.15;
                bounds.upper[i][j] = single_length + 0.15;
            }
            BondType::Double => {
                // Double bonds are ~86% of single bond length
                let ideal = single_length * 0.86;
                bounds.lower[i][j] = ideal - 0.10;
                bounds.upper[i][j] = ideal + 0.10;
            }
            BondType::Triple => {
                // Triple bonds are ~78% of single bond length
                let ideal = single_length * 0.78;
                bounds.lower[i][j] = ideal - 0.08;
                bounds.upper[i][j] = ideal + 0.08;
            }
            BondType::Aromatic => {
                // Aromatic bonds are intermediate
                let ideal = single_length * 0.93;
                bounds.lower[i][j] = ideal - 0.10;
                bounds.upper[i][j] = ideal + 0.10;
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
        let r_ij = (bounds.lower[a_min][a_max] + bounds.upper[a_min][a_max]) / 2.0;
        let r_jk = (bounds.lower[b_min][b_max] + bounds.upper[b_min][b_max]) / 2.0;

        let hyb = crate::molecule::graph::determine_hybridization(angle.atom2, mol);
        let theta0_rad: f64 = match hyb {
            Hybridization::Sp1 => 180.0_f64.to_radians(),
            Hybridization::Sp2 => 120.0_f64.to_radians(),
            Hybridization::Sp3 => 109.47_f64.to_radians(),
        };
        let r_ik_sq = r_ij * r_ij + r_jk * r_jk - 2.0 * r_ij * r_jk * theta0_rad.cos();
        if r_ik_sq > 0.0 {
            let r_ik = r_ik_sq.sqrt();
            let tolerance = 0.15;
            let lower = (r_ik - tolerance).max(0.5);
            let upper = r_ik + tolerance;
            let (a, c) = (angle.atom1.min(angle.atom3), angle.atom1.max(angle.atom3));
            bounds.lower[a][c] = bounds.lower[a][c].max(lower);
            bounds.upper[a][c] = bounds.upper[a][c].min(upper);
            bounds.lower[c][a] = bounds.lower[a][c];
            bounds.upper[c][a] = bounds.upper[a][c];
        }
    }

    // 1-4 bounds (torsion-derived) with torsion angle preferences
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

        // Determine preferred torsion angle based on hybridization of central atoms
        let hyb_j = crate::molecule::graph::determine_hybridization(torsion.atom2, mol);
        let hyb_k = crate::molecule::graph::determine_hybridization(torsion.atom3, mol);

        let (upper, lower) = match (hyb_j, hyb_k) {
            // sp3-sp3: staggered (~60°), use cosine law with theta=60°
            (Hybridization::Sp3, Hybridization::Sp3) => {
                let angle_ij =
                    estimate_angle(torsion.atom1, torsion.atom2, torsion.atom3, mol, &bounds);
                let angle_kl =
                    estimate_angle(torsion.atom2, torsion.atom3, torsion.atom4, mol, &bounds);
                let phi = 60.0_f64.to_radians();
                let r14_sq = r_ij * r_ij + r_jk * r_jk + r_kl * r_kl
                    - 2.0 * r_ij * r_jk * angle_ij.cos()
                    - 2.0 * r_jk * r_kl * angle_kl.cos()
                    + 2.0
                        * r_ij
                        * r_kl
                        * (angle_ij.cos() * angle_kl.cos()
                            - angle_ij.sin() * angle_kl.sin() * phi.cos());
                let r14 = if r14_sq > 0.0 { r14_sq.sqrt() } else { 0.5 };
                (r14 + 0.4, (r14 - 0.4).max(0.5))
            }
            // sp2-sp2: planar (0° or 180°)
            (Hybridization::Sp2, Hybridization::Sp2) => {
                let upper = r_ij + r_jk + r_kl + 0.1;
                let inner =
                    r_ij * r_ij + r_jk * r_jk + r_kl * r_kl - 2.0 * r_ij * r_jk - 2.0 * r_jk * r_kl
                        + 2.0 * r_ij * r_kl;
                let lower = if inner > 0.0 { inner.sqrt() - 0.1 } else { 0.5 };
                (upper, lower)
            }
            // Default: wide bounds
            _ => {
                let upper = r_ij + r_jk + r_kl + 0.1;
                let inner =
                    r_ij * r_ij + r_jk * r_jk + r_kl * r_kl - 2.0 * r_ij * r_jk - 2.0 * r_jk * r_kl
                        + 2.0 * r_ij * r_kl;
                let lower = if inner > 0.0 { inner.sqrt() - 0.1 } else { 0.5 };
                (upper, lower)
            }
        };

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

/// Generate coordinates using classical distance geometry (metric matrix + eigenvector method)
/// Uses 4D embedding then random projection to 3D (standard ETKDG approach)
fn generate_4d_coordinates(bounds: &DistanceBounds) -> Vec<[f64; 3]> {
    let n = bounds.n_atoms;
    if n < 3 {
        return vec![[0.0; 3]; n];
    }

    // Build a target distance matrix with biased sampling
    let mut dist = vec![vec![0.0f64; n]; n];
    for i in 0..n {
        for j in i + 1..n {
            let lo = bounds.lower[i][j];
            let hi = bounds.upper[i][j];
            let range = hi - lo;
            let t = if range < 0.5 {
                lo + range * (0.3 + 0.4 * random_f64())
            } else {
                let u = random_f64();
                lo + range * u * u
            };
            dist[i][j] = t;
            dist[j][i] = t;
        }
    }

    // Build squared distance matrix and double-center for Gram matrix
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

    // Double-centering: B = -0.5 * (D2 - row_means - col_means + grand_mean)
    let mut b = vec![vec![0.0f64; n]; n];
    for i in 0..n {
        for j in 0..n {
            b[i][j] = -0.5 * (d2[i][j] - row_means[i] - row_means[j] + grand_mean);
        }
    }

    // Jacobi eigenvalue decomposition
    let (eigenvalues, eigenvectors) = jacobi_eigen(&b, 100);

    // Sort by eigenvalue (largest first)
    let mut indices: Vec<usize> = (0..n).collect();
    indices.sort_by(|&a, &b| eigenvalues[b].partial_cmp(&eigenvalues[a]).unwrap());

    // Embed in 4D using top 4 eigenvalues
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
                    eigenvectors[idx][i] * 0.001
                };
            }
        }
    }

    // Random projection from 4D to 3D using a random rotation matrix
    // Generate a random 3x4 projection matrix with orthogonal rows
    let mut proj = [[0.0f64; 4]; 3];
    for r in 0..3 {
        for c in 0..4 {
            proj[r][c] = random_f64() * 2.0 - 1.0;
        }
    }

    // Gram-Schmidt orthogonalization of projection rows
    for r in 0..3 {
        for prev in 0..r {
            let dot: f64 = (0..4).map(|c| proj[r][c] * proj[prev][c]).sum();
            for c in 0..4 {
                proj[r][c] -= dot * proj[prev][c];
            }
        }
        let norm: f64 = (0..4).map(|c| proj[r][c] * proj[r][c]).sum::<f64>().sqrt();
        if norm > 1e-10 {
            for c in 0..4 {
                proj[r][c] /= norm;
            }
        }
    }

    // Project 4D coords to 3D
    let mut coords = vec![[0.0; 3]; n];
    for i in 0..n {
        for r in 0..3 {
            coords[i][r] = (0..4).map(|c| proj[r][c] * coords_4d[i][c]).sum();
        }
    }

    coords
}

/// Jacobi eigenvalue decomposition for symmetric matrices
/// Returns (eigenvalues, eigenvectors) where eigenvectors[i] is the i-th eigenvector
fn jacobi_eigen(matrix: &[Vec<f64>], max_sweeps: usize) -> (Vec<f64>, Vec<Vec<f64>>) {
    let n = matrix.len();
    let mut a = matrix.to_vec(); // work on a copy
    let n = a.len();
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

    // v[r][i] contains the i-th eigenvector's r-th component
    // We want eigenvectors[i] = column i of v
    let mut evecs = vec![vec![0.0f64; n]; n];
    for i in 0..n {
        for r in 0..n {
            evecs[i][r] = v[r][i];
        }
    }

    (d, evecs)
}

/// Scale coordinates so that average bonded distance matches expected bond lengths.
/// This corrects for eigenvector embeddings that are too large or too small.
fn scale_to_bonds(bounds: &DistanceBounds, coords: &mut [[f64; 3]]) {
    let n = bounds.n_atoms;
    let mut expected_sum = 0.0f64;
    let mut actual_sum = 0.0f64;
    let mut count = 0usize;

    for i in 0..n {
        for j in (i + 1)..n {
            let lo = bounds.lower[i][j];
            let hi = bounds.upper[i][j];
            if hi - lo < 0.5 {
                // Tight bounds = bonded pair
                let target = (lo + hi) / 2.0;
                let dx = coords[i][0] - coords[j][0];
                let dy = coords[i][1] - coords[j][1];
                let dz = coords[i][2] - coords[j][2];
                let d = (dx * dx + dy * dy + dz * dz).sqrt();
                if d > 1e-10 {
                    expected_sum += target;
                    actual_sum += d;
                    count += 1;
                }
            }
        }
    }

    if count > 0 && actual_sum > 1e-10 {
        let scale = expected_sum / actual_sum;
        // Only apply if scale is reasonable (not NaN/inf, not wildly off)
        if scale.is_finite() && scale > 0.001 && scale < 1000.0 {
            for c in coords.iter_mut() {
                c[0] *= scale;
                c[1] *= scale;
                c[2] *= scale;
            }
        }
    }
}

/// Refine coordinates by minimizing distance bound violations (steepest descent)
/// Applied before force field refinement to establish correct local geometry
fn refine_with_dg(bounds: &DistanceBounds, coords: &mut [[f64; 3]], max_iterations: usize) {
    let n = bounds.n_atoms;
    if n < 2 {
        return;
    }

    // First scale to match expected bond lengths
    scale_to_bonds(bounds, coords);

    let max_atom_step = 0.02; // Max displacement per atom per iteration (Å)

    for _iter in 0..max_iterations {
        let mut grad = vec![[0.0f64; 3]; n];
        let mut error = 0.0f64;

        for i in 0..n {
            for j in (i + 1)..n {
                let dx = coords[i][0] - coords[j][0];
                let dy = coords[i][1] - coords[j][1];
                let dz = coords[i][2] - coords[j][2];
                let d = (dx * dx + dy * dy + dz * dz).sqrt().max(1e-10);

                let lo = bounds.lower[i][j];
                let hi = bounds.upper[i][j];

                let lo_viol = (lo - d).max(0.0);
                let hi_viol = (d - hi).max(0.0);

                if lo_viol < 1e-8 && hi_viol < 1e-8 {
                    continue;
                }

                let range = (hi - lo).max(0.01);
                let w = 1.0 / range;

                error += w * (lo_viol * lo_viol + hi_viol * hi_viol);

                let coeff = 2.0 * w * (hi_viol - lo_viol);

                let fx = coeff * dx / d;
                let fy = coeff * dy / d;
                let fz = coeff * dz / d;

                grad[i][0] += fx;
                grad[i][1] += fy;
                grad[i][2] += fz;
                grad[j][0] -= fx;
                grad[j][1] -= fy;
                grad[j][2] -= fz;
            }
        }

        if error < 1e-10 {
            break;
        }

        // Per-atom step: cap each atom's displacement at max_atom_step
        let max_g = grad
            .iter()
            .map(|g| (g[0] * g[0] + g[1] * g[1] + g[2] * g[2]).sqrt())
            .fold(0.0f64, f64::max)
            .max(1e-10);

        let step = max_atom_step / max_g;

        for i in 0..n {
            coords[i][0] -= step * grad[i][0];
            coords[i][1] -= step * grad[i][1];
            coords[i][2] -= step * grad[i][2];
        }
    }
}

/// Refine coordinates using MMFF94 force field.
/// Two-phase approach: steepest descent for robust initial convergence,
/// then L-BFGS for efficient final optimization.
fn refine_with_ff(mol: &Molecule, coords: &mut [[f64; 3]], config: &ETKDGConfig) {
    let variant = crate::MMFFVariant::MMFF94s;
    let ff = crate::mmff::MMFFForceField::new(mol, variant);
    let n = mol.atoms.len();
    let max_total_iters = config.max_iterations.max(1000);

    // Phase 1: Steepest descent with Armijo line search
    // Robust for bad starting geometries where L-BFGS struggles
    let sd_iters = max_total_iters / 2;

    for _iter in 0..sd_iters {
        let (energy, grad) = ff.calculate_energy_and_gradient(coords);

        // Check convergence by max force
        let max_g = grad
            .iter()
            .map(|g| (g[0] * g[0] + g[1] * g[1] + g[2] * g[2]).sqrt())
            .fold(0.0f64, f64::max);

        if max_g < 0.5 {
            break;
        }

        // Descent direction: d = -g
        // Compute g^T * d = -|g|^2
        let g_dot_d: f64 = -grad
            .iter()
            .flat_map(|g| g.iter())
            .map(|gi| gi * gi)
            .sum::<f64>();

        // Armijo backtracking: start with step that gives max 0.1 Å displacement
        let max_component = grad
            .iter()
            .flat_map(|g| g.iter())
            .map(|gi| gi.abs())
            .fold(0.0f64, f64::max)
            .max(1e-10);
        let mut alpha = 0.1 / max_component;

        let mut accepted = false;
        for _backtrack in 0..30 {
            let mut new_coords = coords.to_vec();
            for i in 0..n {
                new_coords[i][0] -= alpha * grad[i][0];
                new_coords[i][1] -= alpha * grad[i][1];
                new_coords[i][2] -= alpha * grad[i][2];
            }
            let new_energy = ff.calculate_energy(&new_coords);

            // Armijo condition: f(x + alpha*d) <= f(x) + c1 * alpha * g^T * d
            if new_energy <= energy + 1e-4 * alpha * g_dot_d {
                for i in 0..n {
                    coords[i] = new_coords[i];
                }
                accepted = true;
                break;
            }
            alpha *= 0.5;
        }

        if !accepted {
            break; // Can't make progress
        }
    }

    // Phase 2: L-BFGS for fine-tuning
    let conv = crate::ConvergenceOptions {
        max_force: 0.01,
        rms_force: 0.001,
        energy_change: 1e-6,
        max_iterations: max_total_iters / 2,
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
        let mut coords_3d = generate_4d_coordinates(&bounds);
        refine_with_dg(&bounds, &mut coords_3d, 500);
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
