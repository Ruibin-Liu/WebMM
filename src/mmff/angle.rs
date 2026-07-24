//! Angle bending term for MMFF94

use super::mmff_tables;
use super::MMFFAtomType;

#[derive(Debug, Clone, Copy)]
pub struct AngleParams {
    pub k_theta: f64,
    pub theta0: f64,
}

pub fn get_angle_params(
    type1: MMFFAtomType,
    type2: MMFFAtomType,
    type3: MMFFAtomType,
) -> Option<AngleParams> {
    match (type1, type2, type3) {
        (MMFFAtomType::F, MMFFAtomType::P_3D, MMFFAtomType::F) => {
            return Some(AngleParams {
                k_theta: 0.80,
                theta0: 180.0,
            });
        }
        (MMFFAtomType::F, MMFFAtomType::P_3D, MMFFAtomType::C_3)
        | (MMFFAtomType::C_3, MMFFAtomType::P_3D, MMFFAtomType::F) => {
            return Some(AngleParams {
                k_theta: 0.80,
                theta0: 90.0,
            });
        }
        (MMFFAtomType::F, MMFFAtomType::S_3D2, MMFFAtomType::F) => {
            return Some(AngleParams {
                k_theta: 0.80,
                theta0: 90.0,
            });
        }
        _ => {}
    }

    let bt1 = super::base_type(type1);
    let bt2 = super::base_type(type2);
    let bt3 = super::base_type(type3);

    let r0_ij = super::bond::get_bond_params(bt1, bt2, super::super::molecule::BondType::Single)
        .map(|p| p.r0)
        .unwrap_or(1.5);
    let r0_jk = super::bond::get_bond_params(bt2, bt3, super::super::molecule::BondType::Single)
        .map(|p| p.r0)
        .unwrap_or(1.5);

    get_angle_params_from_table(type1, type2, type3, 0, r0_ij, r0_jk, 0)
}

#[allow(clippy::too_many_arguments)]
pub fn get_angle_params_with_bond_info(
    type1: MMFFAtomType,
    type2: MMFFAtomType,
    type3: MMFFAtomType,
    bond_type_ij: u8,
    bond_type_jk: u8,
    ring_size: u8,
    r0_ij: f64,
    r0_jk: f64,
) -> Option<AngleParams> {
    match (type1, type2, type3) {
        (MMFFAtomType::F, MMFFAtomType::P_3D, MMFFAtomType::F) => {
            return Some(AngleParams {
                k_theta: 0.80,
                theta0: 180.0,
            });
        }
        (MMFFAtomType::F, MMFFAtomType::P_3D, MMFFAtomType::C_3)
        | (MMFFAtomType::C_3, MMFFAtomType::P_3D, MMFFAtomType::F) => {
            return Some(AngleParams {
                k_theta: 0.80,
                theta0: 90.0,
            });
        }
        (MMFFAtomType::F, MMFFAtomType::S_3D2, MMFFAtomType::F) => {
            return Some(AngleParams {
                k_theta: 0.80,
                theta0: 90.0,
            });
        }
        _ => {}
    }

    let angle_type = mmff_tables::compute_angle_type(bond_type_ij, bond_type_jk, ring_size);
    get_angle_params_from_table(type1, type2, type3, angle_type, r0_ij, r0_jk, ring_size)
}

fn get_angle_params_from_table(
    type1: MMFFAtomType,
    type2: MMFFAtomType,
    type3: MMFFAtomType,
    angle_type: u8,
    r0_ij: f64,
    r0_jk: f64,
    ring_size: u8,
) -> Option<AngleParams> {
    let ti = super::params::mmff_type_id(type1);
    let tj = super::params::mmff_type_id(type2);
    let tk = super::params::mmff_type_id(type3);

    if let Some((ka, t0)) = mmff_tables::lookup_angle_params(ti, tj, tk, angle_type) {
        if ka.abs() > 1e-10 {
            return Some(AngleParams {
                k_theta: ka,
                theta0: t0,
            });
        }
    }

    if let Some((ka, theta0)) =
        mmff_tables::empirical_angle_params(ti, tj, tk, r0_ij, r0_jk, ring_size)
    {
        return Some(AngleParams {
            k_theta: ka,
            theta0,
        });
    }

    None
}

/// Calculate angle bending energy
///
/// MMFF94 anharmonic angle bend (RDKit-compatible):
///   E = 0.5 * c2 * ka * dtheta^2 * (1 + cb * dtheta)
/// where cb = -0.006981317, c2 = 143.9325 * (pi/180)^2 = 0.04385
/// dtheta is in degrees
///
/// For linear angles (theta0 >= 179.0), uses:
///   E = 143.9325 * ka * (1 + cos(theta))
pub fn angle_energy(
    coords: &[[f64; 3]],
    i: usize,
    j: usize,
    k: usize,
    params: &AngleParams,
) -> f64 {
    let theta_rad = calculate_angle(coords, i, j, k);

    if params.theta0 >= 179.0 {
        let cos_theta = theta_rad.cos();
        143.9325 * params.k_theta * (1.0 + cos_theta)
    } else {
        let theta_deg = theta_rad.to_degrees();
        let dtheta = theta_deg - params.theta0;

        let c2 = 143.9325 * (std::f64::consts::PI / 180.0).powi(2);
        let cb = -0.006981317;

        0.5 * c2 * params.k_theta * dtheta * dtheta * (1.0 + cb * dtheta)
    }
}

/// Calculate angle between three atoms
fn calculate_angle(coords: &[[f64; 3]], i: usize, j: usize, k: usize) -> f64 {
    let v1 = [
        coords[i][0] - coords[j][0],
        coords[i][1] - coords[j][1],
        coords[i][2] - coords[j][2],
    ];
    let v2 = [
        coords[k][0] - coords[j][0],
        coords[k][1] - coords[j][1],
        coords[k][2] - coords[j][2],
    ];

    let v1_norm = (v1[0].powi(2) + v1[1].powi(2) + v1[2].powi(2)).sqrt();
    let v2_norm = (v2[0].powi(2) + v2[1].powi(2) + v2[2].powi(2)).sqrt();

    if v1_norm == 0.0 || v2_norm == 0.0 {
        return std::f64::consts::PI;
    }

    let dot = v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
    let cos_theta = dot / (v1_norm * v2_norm);

    cos_theta.clamp(-1.0, 1.0).acos()
}

pub fn angle_gradient(
    coords: &[[f64; 3]],
    atom1: usize,
    atom2: usize,
    atom3: usize,
    params: &AngleParams,
) -> ([f64; 3], [f64; 3], [f64; 3]) {
    let r1 = [
        coords[atom1][0] - coords[atom2][0],
        coords[atom1][1] - coords[atom2][1],
        coords[atom1][2] - coords[atom2][2],
    ];

    let r2 = [
        coords[atom3][0] - coords[atom2][0],
        coords[atom3][1] - coords[atom2][1],
        coords[atom3][2] - coords[atom2][2],
    ];

    let r1_norm = (r1[0].powi(2) + r1[1].powi(2) + r1[2].powi(2)).sqrt();
    let r2_norm = (r2[0].powi(2) + r2[1].powi(2) + r2[2].powi(2)).sqrt();

    if r1_norm < 1e-10 || r2_norm < 1e-10 {
        return ([0.0; 3], [0.0; 3], [0.0; 3]);
    }

    let dot = r1[0] * r2[0] + r1[1] * r2[1] + r1[2] * r2[2];
    let cos_theta = dot / (r1_norm * r2_norm);
    let theta = cos_theta.clamp(-1.0, 1.0).acos();

    let sin_theta = theta.sin();

    if sin_theta.abs() < 1e-10 {
        return ([0.0; 3], [0.0; 3], [0.0; 3]);
    }

    let d_e_dtheta = if params.theta0 >= 179.0 {
        -143.9325 * params.k_theta * sin_theta
    } else {
        let theta_deg = theta.to_degrees();
        let angle_term = theta_deg - params.theta0;
        let c2 = 143.9325 * (std::f64::consts::PI / 180.0).powi(2);
        let cb = -0.006981317;
        let rad2deg = 180.0 / std::f64::consts::PI;
        rad2deg * c2 * params.k_theta * angle_term * (1.0 + 1.5 * cb * angle_term)
    };

    // dtheta/dr1 = -1/(|r1|*sin(theta)) * (r2/|r2| - cos(theta)*r1/|r1|)
    // dtheta/dr3 = -1/(|r2|*sin(theta)) * (r1/|r1| - cos(theta)*r2/|r2|)
    // dtheta/dr2 = -(dtheta/dr1 + dtheta/dr3)
    // dE/dr = dE/dtheta * dtheta/dr

    let inv_r1_sq = 1.0 / (r1_norm * r1_norm);
    let inv_r2_sq = 1.0 / (r2_norm * r2_norm);
    let inv_r1_r2 = 1.0 / (r1_norm * r2_norm);

    // Gradient for atom1
    let grad1 = [
        -d_e_dtheta * (r2[0] * inv_r1_r2 - cos_theta * r1[0] * inv_r1_sq) / sin_theta,
        -d_e_dtheta * (r2[1] * inv_r1_r2 - cos_theta * r1[1] * inv_r1_sq) / sin_theta,
        -d_e_dtheta * (r2[2] * inv_r1_r2 - cos_theta * r1[2] * inv_r1_sq) / sin_theta,
    ];

    // Gradient for atom3
    let grad3 = [
        -d_e_dtheta * (r1[0] * inv_r1_r2 - cos_theta * r2[0] * inv_r2_sq) / sin_theta,
        -d_e_dtheta * (r1[1] * inv_r1_r2 - cos_theta * r2[1] * inv_r2_sq) / sin_theta,
        -d_e_dtheta * (r1[2] * inv_r1_r2 - cos_theta * r2[2] * inv_r2_sq) / sin_theta,
    ];

    // Gradient for atom2 (center of angle)
    let grad2 = [
        -(grad1[0] + grad3[0]),
        -(grad1[1] + grad3[1]),
        -(grad1[2] + grad3[2]),
    ];

    (grad1, grad2, grad3)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_angle_energy() {
        let coords = vec![
            [1.526, 0.0, 0.0],   // H
            [0.0, 0.0, 0.0],     // C
            [1.526, 0.934, 0.0], // H
        ];

        let params = AngleParams {
            k_theta: 0.8,
            theta0: 109.47,
        };
        let energy = angle_energy(&coords, 0, 1, 2, &params);

        assert!(energy.is_finite());
    }

    #[test]
    fn test_angle_energy_zero_at_equilibrium() {
        let theta0_deg = 108.836;
        let half = (theta0_deg / 2.0_f64).to_radians();
        let r = 1.093;
        let coords = vec![
            [r * half.cos(), r * half.sin(), 0.0],
            [0.0, 0.0, 0.0],
            [r * half.cos(), -r * half.sin(), 0.0],
        ];

        let params = AngleParams {
            k_theta: 0.516,
            theta0: theta0_deg,
        };
        let energy = angle_energy(&coords, 0, 1, 2, &params);

        assert!(
            energy.abs() < 0.01,
            "Angle energy at equilibrium should be ~0, got {}",
            energy
        );
    }

    #[test]
    fn test_angle_energy_straight_line() {
        // Collinear atoms: angle = 180 degrees, should be higher energy
        let coords = vec![[-2.0, 0.0, 0.0], [0.0, 0.0, 0.0], [2.0, 0.0, 0.0]];

        let params = AngleParams {
            k_theta: 0.8,
            theta0: 109.47,
        };
        let energy = angle_energy(&coords, 0, 1, 2, &params);

        assert!(energy.is_finite());
        assert!(energy > 0.0, "Deviated angle should have positive energy");
    }

    #[test]
    fn test_angle_gradient_numerical() {
        // Use a non-equilibrium angle to test the gradient
        let theta0_deg = 108.836;
        let actual_deg = 115.0; // Perturbed from equilibrium
        let half = (actual_deg / 2.0_f64).to_radians();
        let r = 1.093;
        let coords = vec![
            [r * half.cos(), r * half.sin(), 0.0],
            [0.0, 0.0, 0.0],
            [r * half.cos(), -r * half.sin(), 0.0],
        ];

        let params = AngleParams {
            k_theta: 0.516,
            theta0: theta0_deg,
        };

        let (g1, g2, g3) = angle_gradient(&coords, 0, 1, 2, &params);

        let eps = 1e-7;
        let e0 = angle_energy(&coords, 0, 1, 2, &params);

        for (atom_idx, grad) in [(0, g1), (1, g2), (2, g3)] {
            for dim in 0..3 {
                let mut coords_p = coords.clone();
                coords_p[atom_idx][dim] += eps;
                let e_plus = angle_energy(&coords_p, 0, 1, 2, &params);
                let num_grad = (e_plus - e0) / eps;
                let abs_diff = (grad[dim] - num_grad).abs();
                let max_mag = grad[dim].abs().max(num_grad.abs());
                let rel_err = if max_mag > 1e-4 {
                    abs_diff / max_mag
                } else {
                    abs_diff
                };
                assert!(
                    rel_err < 1e-3,
                    "Analytical grad[{}] = {} but numerical = {} for atom {} (rel_err={})",
                    dim,
                    grad[dim],
                    num_grad,
                    atom_idx,
                    rel_err
                );
            }
        }
    }

    #[test]
    fn test_angle_gradient_coincident_atoms() {
        // Two atoms at the same position should not panic
        let coords = vec![[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];

        let params = AngleParams {
            k_theta: 0.8,
            theta0: 109.47,
        };

        let (g1, g2, g3) = angle_gradient(&coords, 0, 1, 2, &params);
        for grad in [g1, g2, g3] {
            for &v in &grad {
                assert!(v.is_finite(), "Gradient should be finite, got NaN/Inf");
            }
        }
    }

    #[test]
    fn test_angle_gradient_linear_atoms() {
        // Collinear atoms (sin_theta ~ 0) should not panic
        let coords = vec![[-1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];

        let params = AngleParams {
            k_theta: 0.8,
            theta0: 109.47,
        };

        let (g1, g2, g3) = angle_gradient(&coords, 0, 1, 2, &params);
        for grad in [g1, g2, g3] {
            for &v in &grad {
                assert!(v.is_finite(), "Gradient should be finite for linear case");
            }
        }
    }
}
