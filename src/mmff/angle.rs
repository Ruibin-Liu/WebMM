//! Angle bending term for MMFF94

use super::MMFFAtomType;

/// Angle bending parameters
#[derive(Debug, Clone, Copy)]
pub struct AngleParams {
    pub k_theta: f64, // mdyn·rad⁻²
    pub theta0: f64,  // degrees
}

/// Get angle parameters for atom types
pub fn get_angle_params(
    type1: MMFFAtomType,
    type2: MMFFAtomType,
    type3: MMFFAtomType,
) -> Option<AngleParams> {
    // TODO: Load from parameter tables
    // For now, return simple default values
    match (type1, type2, type3) {
        // C-C-C angles
        (MMFFAtomType::C_3, MMFFAtomType::C_3, MMFFAtomType::C_3) => Some(AngleParams {
            k_theta: 1.05,
            theta0: 109.47,
        }),
        (MMFFAtomType::C_AR, MMFFAtomType::C_AR, MMFFAtomType::C_AR) => Some(AngleParams {
            k_theta: 1.05,
            theta0: 120.0,
        }),
        (MMFFAtomType::C_2, MMFFAtomType::C_3, MMFFAtomType::C_3) => Some(AngleParams {
            k_theta: 1.15,
            theta0: 121.0,
        }),

        // H-C-X angles (with symmetric matching for H at either end)
        (MMFFAtomType::C_3, MMFFAtomType::C_3, MMFFAtomType::H)
        | (MMFFAtomType::H, MMFFAtomType::C_3, MMFFAtomType::C_3) => Some(AngleParams {
            k_theta: 0.80,
            theta0: 109.47,
        }),
        (MMFFAtomType::H, MMFFAtomType::C_3, MMFFAtomType::H) => Some(AngleParams {
            k_theta: 0.77,
            theta0: 109.47,
        }),
        (MMFFAtomType::C_3, MMFFAtomType::N_3, MMFFAtomType::H)
        | (MMFFAtomType::H, MMFFAtomType::N_3, MMFFAtomType::C_3) => Some(AngleParams {
            k_theta: 0.60,
            theta0: 109.47,
        }),
        (MMFFAtomType::C_3, MMFFAtomType::O_3, MMFFAtomType::H)
        | (MMFFAtomType::H, MMFFAtomType::O_3, MMFFAtomType::C_3) => Some(AngleParams {
            k_theta: 0.75,
            theta0: 109.47,
        }),

        // H-C_2-X angles
        (MMFFAtomType::C_3, MMFFAtomType::C_2, MMFFAtomType::H)
        | (MMFFAtomType::H, MMFFAtomType::C_2, MMFFAtomType::C_3) => Some(AngleParams {
            k_theta: 0.50,
            theta0: 120.0,
        }),
        (MMFFAtomType::H, MMFFAtomType::C_2, MMFFAtomType::H) => Some(AngleParams {
            k_theta: 0.45,
            theta0: 120.0,
        }),

        // C_2-C_3-C_3
        (MMFFAtomType::C_3, MMFFAtomType::C_3, MMFFAtomType::C_2) => Some(AngleParams {
            k_theta: 1.15,
            theta0: 121.0,
        }),

        _ => None,
    }
}

/// Calculate angle bending energy
pub fn angle_energy(
    coords: &[[f64; 3]],
    i: usize,
    j: usize, // Central atom
    k: usize,
    params: &AngleParams,
) -> f64 {
    let theta_rad = calculate_angle(coords, i, j, k);
    let theta_deg = theta_rad.to_degrees();
    let dtheta = theta_deg - params.theta0;

    // Convert to radians
    let dtheta_rad = dtheta.to_radians();

    // E_angle = 0.000043945 * k_theta * (theta - theta0)^2  (kcal/mol)
    0.000043945 * params.k_theta * dtheta_rad * dtheta_rad
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

    // dE/dtheta = 2 * 0.000043945 * k_theta * (theta - theta0)
    let prefactor = 2.0 * 0.000043945 * params.k_theta * (theta - params.theta0.to_radians());

    // Gradient for atom1
    let grad1 = [
        prefactor
            * (r2[0] / (r1_norm * r1_norm * sin_theta) - cos_theta * r1[0] / (r1_norm * r1_norm)),
        prefactor
            * (r2[1] / (r1_norm * r1_norm * sin_theta) - cos_theta * r1[1] / (r1_norm * r1_norm)),
        prefactor
            * (r2[2] / (r1_norm * r1_norm * sin_theta) - cos_theta * r1[2] / (r1_norm * r1_norm)),
    ];

    // Gradient for atom2 (center of angle)
    let grad2 = [
        -grad1[0]
            - prefactor
                * (r1[0] / (r2_norm * r2_norm * sin_theta)
                    - cos_theta * r2[0] / (r2_norm * r2_norm)),
        -grad1[1]
            - prefactor
                * (r1[1] / (r2_norm * r2_norm * sin_theta)
                    - cos_theta * r2[1] / (r2_norm * r2_norm)),
        -grad1[2]
            - prefactor
                * (r1[2] / (r2_norm * r2_norm * sin_theta)
                    - cos_theta * r2[2] / (r2_norm * r2_norm)),
    ];

    // Gradient for atom3
    let grad3 = [
        prefactor
            * (r1[0] / (r2_norm * r2_norm * sin_theta) - cos_theta * r2[0] / (r2_norm * r2_norm)),
        prefactor
            * (r1[1] / (r2_norm * r2_norm * sin_theta) - cos_theta * r2[1] / (r2_norm * r2_norm)),
        prefactor
            * (r1[2] / (r2_norm * r2_norm * sin_theta) - cos_theta * r2[2] / (r2_norm * r2_norm)),
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
}
