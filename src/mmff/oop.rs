//! Out-of-plane bending term for MMFF94

use super::super::mmff::MMFFAtomType;
use super::super::mmff::MMFFVariant;

/// Out-of-plane parameters
#[derive(Debug, Clone, Copy)]
pub struct OOPParams {
    pub k_oop: f64, // kcal/mol/rad²
}

/// Get OOP parameters for central atom type
pub fn get_oop_params(central_type: MMFFAtomType, mmff_variant: MMFFVariant) -> OOPParams {
    // TODO: Load from parameter tables
    // For now, return simple default values
    match central_type {
        MMFFAtomType::C_3 | MMFFAtomType::C_2 | MMFFAtomType::C_AR => OOPParams { k_oop: 0.04 },
        MMFFAtomType::N_3 | MMFFAtomType::N_2 | MMFFAtomType::N_AR => OOPParams { k_oop: 0.04 },
        MMFFAtomType::O_2 | MMFFAtomType::O_3 => OOPParams { k_oop: 0.04 },
        _ => OOPParams { k_oop: 0.04 },
    }
}

/// Calculate out-of-plane bending energy
pub fn oop_energy(
    coords: &[[f64; 3]],
    central: usize,
    i: usize,
    j: usize,
    k: usize,
    params: &OOPParams,
) -> f64 {
    let chi = calculate_oop_angle(coords, central, i, j, k);

    // E_oop = k_oop * sum[n=1 to 4] (C_n * (1 - cos(n*chi)))
    let cos_chi = chi.cos();
    let cos_2chi = (2.0 * chi).cos();
    let cos_3chi = (3.0 * chi).cos();
    let cos_4chi = (4.0 * chi).cos();

    params.k_oop * (1.0 - cos_chi + 1.0 - cos_2chi + 1.0 - cos_3chi + 1.0 - cos_4chi)
}

/// Calculate out-of-plane angle (deviation from plane)
fn calculate_oop_angle(coords: &[[f64; 3]], central: usize, i: usize, j: usize, k: usize) -> f64 {
    // Vectors from central to bonded atoms
    let v1 = [
        coords[i][0] - coords[central][0],
        coords[i][1] - coords[central][1],
        coords[i][2] - coords[central][2],
    ];
    let v2 = [
        coords[j][0] - coords[central][0],
        coords[j][1] - coords[central][1],
        coords[j][2] - coords[central][2],
    ];
    let v3 = [
        coords[k][0] - coords[central][0],
        coords[k][1] - coords[central][1],
        coords[k][2] - coords[central][2],
    ];

    // Plane normal = (v1 - v2) cross (v1 - v3)
    let diff12 = [v1[0] - v2[0], v1[1] - v2[1], v1[2] - v2[2]];
    let diff13 = [v1[0] - v3[0], v1[1] - v3[1], v1[2] - v3[2]];

    let normal = [
        diff12[1] * diff13[2] - diff12[2] * diff13[1],
        diff12[2] * diff13[0] - diff12[0] * diff13[2],
        diff12[0] * diff13[1] - diff12[1] * diff13[0],
    ];

    let normal_norm = (normal[0].powi(2) + normal[1].powi(2) + normal[2].powi(2)).sqrt();

    if normal_norm == 0.0 {
        return 0.0;
    }

    // Normalize normal
    let normal = [
        normal[0] / normal_norm,
        normal[1] / normal_norm,
        normal[2] / normal_norm,
    ];

    // OOP angle is the angle between v1 and the plane
    // Use average of angles with v2 and v3
    let v1_norm = (v1[0].powi(2) + v1[1].powi(2) + v1[2].powi(2)).sqrt();
    if v1_norm == 0.0 {
        return 0.0;
    }

    let v1_unit = [v1[0] / v1_norm, v1[1] / v1_norm, v1[2] / v1_norm];

    // Project v1 onto normal (distance from plane)
    let dot = v1_unit[0] * normal[0] + v1_unit[1] * normal[1] + v1_unit[2] * normal[2];
    let dist_from_plane = dot.abs();

    // OOP angle = arcsin(distance)
    dist_from_plane.asin()
}

pub fn oop_gradient(
    coords: &[[f64; 3]],
    central: usize,
    atom1: usize,
    atom2: usize,
    atom3: usize,
    params: &OOPParams,
) -> ([f64; 3], [f64; 3], [f64; 3], [f64; 3]) {
    // Simplified OOP gradient
    let chi = oop_energy(coords, central, atom1, atom2, atom3, params);
    // chi is already the angle in radians
    let sin_chi = chi.sin();

    // Force magnitude
    let force = 2.0 * params.k_oop * sin_chi;

    // Normal to plane
    let v1 = [
        coords[atom1][0] - coords[central][0],
        coords[atom1][1] - coords[central][1],
        coords[atom1][2] - coords[central][2],
    ];

    let v2 = [
        coords[atom2][0] - coords[central][0],
        coords[atom2][1] - coords[central][1],
        coords[atom2][2] - coords[central][2],
    ];

    let v3 = [
        coords[atom3][0] - coords[central][0],
        coords[atom3][1] - coords[central][1],
        coords[atom3][2] - coords[central][2],
    ];

    // Cross product to get normal
    let normal = [
        v1[1] * v2[2] - v1[2] * v2[1],
        v1[2] * v2[0] - v1[0] * v2[2],
        v1[0] * v2[1] - v1[1] * v2[0],
    ];

    let normal_norm = (normal[0].powi(2) + normal[1].powi(2) + normal[2].powi(2)).sqrt();
    let normal_unit = [
        normal[0] / normal_norm,
        normal[1] / normal_norm,
        normal[2] / normal_norm,
    ];

    // Gradient for central atom
    let grad_central = [
        -force * normal_unit[0],
        -force * normal_unit[1],
        -force * normal_unit[2],
    ];

    // Gradient for out-of-plane atoms (simplified)
    let grad1 = [force * 0.3, force * 0.3, force * 0.3];
    let grad2 = [force * 0.3, force * 0.3, force * 0.3];
    let grad3 = [force * 0.3, force * 0.3, force * 0.3];

    (grad_central, grad1, grad2, grad3)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mmff::MMFFVariant;

    #[test]
    fn test_oop_energy() {
        let coords = vec![
            [1.526, 0.0, 0.0],    // H
            [0.0, 0.0, 0.0],      // C
            [1.526, 0.934, 0.0],  // H
            [0.0, 0.934, -1.526], // H (out of plane)
        ];

        let params = OOPParams { k_oop: 0.04 };
        let energy = oop_energy(&coords, 1, 0, 2, 3, &params);

        assert!(energy.is_finite());
    }
}
