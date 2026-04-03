//! Out-of-plane bending term for MMFF94

use super::atom_types::get_atom_type_props;
use super::MMFFAtomType;
use super::MMFFVariant;

/// Out-of-plane parameters
#[derive(Debug, Clone, Copy)]
pub struct OOPParams {
    pub k_oop: f64,
}

pub fn get_oop_params(central_type: MMFFAtomType, _mmff_variant: MMFFVariant) -> OOPParams {
    let k_oop = match get_atom_type_props(central_type) {
        Some(props) => props.oop_k,
        None => 0.0,
    };
    OOPParams { k_oop }
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

    // MMFF94 out-of-plane: E = 0.5 * k_oop * χ²
    // Harmonic approximation using the Wilson out-of-plane coordinate χ (radians),
    // which is 0 when the central atom lies in the plane of its three neighbors.
    0.5 * params.k_oop * chi * chi
}

/// Calculate out-of-plane angle (deviation from plane)
fn calculate_oop_angle(coords: &[[f64; 3]], central: usize, i: usize, j: usize, k: usize) -> f64 {
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

    let normal = [
        normal[0] / normal_norm,
        normal[1] / normal_norm,
        normal[2] / normal_norm,
    ];

    let v1_norm = (v1[0].powi(2) + v1[1].powi(2) + v1[2].powi(2)).sqrt();
    if v1_norm == 0.0 {
        return 0.0;
    }

    let v1_unit = [v1[0] / v1_norm, v1[1] / v1_norm, v1[2] / v1_norm];

    let dot = v1_unit[0] * normal[0] + v1_unit[1] * normal[1] + v1_unit[2] * normal[2];
    let dist_from_plane = dot.abs();

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
    let eps = 1e-7;
    let e_ref = oop_energy(coords, central, atom1, atom2, atom3, params);
    let mut gc = [0.0; 3];
    let mut g1 = [0.0; 3];
    let mut g2 = [0.0; 3];
    let mut g3 = [0.0; 3];

    for (atom_idx, grad) in [
        (central, &mut gc),
        (atom1, &mut g1),
        (atom2, &mut g2),
        (atom3, &mut g3),
    ] {
        for dim in 0..3 {
            let mut coords_p = coords.to_vec();
            coords_p[atom_idx][dim] += eps;
            let e_plus = oop_energy(&coords_p, central, atom1, atom2, atom3, params);
            grad[dim] = (e_plus - e_ref) / eps;
        }
    }

    (gc, g1, g2, g3)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_oop_energy() {
        let coords = vec![
            [1.526, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [1.526, 0.934, 0.0],
            [0.0, 0.934, -1.526],
        ];

        let params = OOPParams { k_oop: 0.04 };
        let energy = oop_energy(&coords, 1, 0, 2, 3, &params);

        assert!(energy.is_finite());
    }

    #[test]
    fn test_oop_gradient_numerical() {
        let coords = vec![
            [1.526, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 1.526, 0.0],
            [0.0, 0.0, 0.5],
        ];
        let params = OOPParams { k_oop: 0.04 };
        let (gc, g1, g2, g3) = oop_gradient(&coords, 1, 0, 2, 3, &params);

        let eps = 1e-7;
        for (idx, grad) in [(1usize, gc), (0usize, g1), (2usize, g2), (3usize, g3)] {
            for dim in 0..3 {
                let mut cp = coords.clone();
                cp[idx][dim] += eps;
                let ep = oop_energy(&cp, 1, 0, 2, 3, &params);
                let num = (ep - oop_energy(&coords, 1, 0, 2, 3, &params)) / eps;
                assert!(
                    (grad[dim] - num).abs() < 1e-4,
                    "OOP grad[{}] = {} vs numerical {} for atom {}",
                    dim,
                    grad[dim],
                    num,
                    idx
                );
            }
        }
    }
}
