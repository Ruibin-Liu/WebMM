//! van der Waals (buffered 14-7 potential) term for MMFF94

use super::atom_types::get_atom_type_props;
use super::MMFFAtomType;

/// van der Waals parameters
#[derive(Debug, Clone, Copy)]
pub struct VDWParams {
    pub r0: f64,
    pub epsilon: f64,
}

pub fn get_vdw_params(atom_type: MMFFAtomType) -> VDWParams {
    let atom_type = super::base_type(atom_type);
    match get_atom_type_props(atom_type) {
        Some(props) => VDWParams {
            r0: props.vdw_r,
            epsilon: props.vdw_eps,
        },
        None => VDWParams {
            r0: 1.7,
            epsilon: 0.07,
        },
    }
}

/// Calculate VDW energy and gradient
///
/// MMFF94 Buffered 14-7 potential (RDKit-compatible):
///   a = 1.07 * R* / (R + 0.07 * R*)
///   b = 1.12 * R*^7 / (R^7 + 0.12 * R*^7) - 2.0
///   E = ε * a^7 * b
///
/// Minimum at R = R* gives E = -ε.
/// For 1-4 interactions, the energy is scaled by 0.75.
pub fn vdw_energy_and_gradient(
    coords: &[[f64; 3]],
    i: usize,
    j: usize,
    params_i: &VDWParams,
    params_j: &VDWParams,
    is_14: bool,
) -> (f64, [f64; 3], [f64; 3]) {
    let r_vec = [
        coords[j][0] - coords[i][0],
        coords[j][1] - coords[i][1],
        coords[j][2] - coords[i][2],
    ];
    let r_sq = r_vec[0] * r_vec[0] + r_vec[1] * r_vec[1] + r_vec[2] * r_vec[2];
    let r = r_sq.sqrt();

    if r < 1e-10 {
        return (0.0, [0.0; 3], [0.0; 3]);
    }

    // R*_ij = 0.5 * (R*_ii + R*_jj)
    let r_star = 0.5 * (params_i.r0 + params_j.r0);
    let epsilon = (params_i.epsilon * params_j.epsilon).sqrt();

    // RDKit buffered 14-7 with damping:
    // a = 1.07 * R* / (R + 0.07 * R*)
    // b = 1.12 * R*^7 / (R^7 + 0.12 * R*^7) - 2.0
    // E = ε * a^7 * b
    let vdw1 = 1.07;
    let vdw1m1 = 0.07;
    let vdw2 = 1.12;

    let q = r / r_star;
    let q7 = q.powi(7);

    let t = vdw1 / (q + vdw1m1);
    let t7 = t.powi(7);

    let b = vdw2 / (q7 + 0.12) - 2.0;
    let scale = if is_14 { 0.75 } else { 1.0 };
    let energy = scale * epsilon * t7 * b;

    // Gradient via numerical differentiation
    let eps = 1e-7;
    let mut grad_i = [0.0; 3];
    let mut grad_j = [0.0; 3];
    for dim in 0..3 {
        let mut cp: Vec<[f64; 3]> = coords.to_vec();
        cp[i][dim] += eps;
        let (ep, _, _) = vdw_energy_only(&cp, i, j, r_star, epsilon, scale);
        let (e0, _, _) = vdw_energy_only(coords, i, j, r_star, epsilon, scale);
        let num = (ep - e0) / eps;
        grad_i[dim] = num;
        grad_j[dim] = -num;
    }

    (energy, grad_i, grad_j)
}

/// VDW energy only (for numerical gradient)
fn vdw_energy_only(
    coords: &[[f64; 3]],
    i: usize,
    j: usize,
    r_star: f64,
    epsilon: f64,
    scale: f64,
) -> (f64, [f64; 3], [f64; 3]) {
    let dx = coords[j][0] - coords[i][0];
    let dy = coords[j][1] - coords[i][1];
    let dz = coords[j][2] - coords[i][2];
    let r = (dx * dx + dy * dy + dz * dz).sqrt();

    if r < 1e-10 {
        return (0.0, [0.0; 3], [0.0; 3]);
    }

    let vdw1 = 1.07;
    let vdw1m1 = 0.07;
    let vdw2 = 1.12;

    let q = r / r_star;
    let q7 = q.powi(7);

    let t = vdw1 / (q + vdw1m1);
    let t7 = t.powi(7);

    let b = vdw2 / (q7 + 0.12) - 2.0;
    let energy = scale * epsilon * t7 * b;

    (energy, [0.0; 3], [0.0; 3])
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_params(r0: f64, epsilon: f64) -> VDWParams {
        VDWParams { r0, epsilon }
    }

    #[test]
    fn test_vdw_energy_has_minimum() {
        let pi = make_params(1.7, 0.07);
        let pj = make_params(1.7, 0.07);
        let coords_close = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let coords_eq = [[0.0, 0.0, 0.0], [1.7, 0.0, 0.0]];
        let coords_far = [[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]];

        let (e_close, _, _) = vdw_energy_and_gradient(&coords_close, 0, 1, &pi, &pj, false);
        let (e_eq, _, _) = vdw_energy_and_gradient(&coords_eq, 0, 1, &pi, &pj, false);
        let (e_far, _, _) = vdw_energy_and_gradient(&coords_far, 0, 1, &pi, &pj, false);

        assert!(
            e_close > 0.0,
            "VDW should be repulsive at short range: {}",
            e_close
        );
        assert!(
            (e_eq - (-0.07)).abs() < 0.01,
            "VDW should have attractive well near r0: got {}, expected ~-0.07",
            e_eq
        );
        assert!(
            e_far > e_eq,
            "VDW should rise toward zero at long range: e_far={} > e_eq={}",
            e_far,
            e_eq
        );
    }

    #[test]
    fn test_vdw_14_scaling() {
        let pi = make_params(1.7, 0.07);
        let pj = make_params(1.7, 0.07);
        let coords = [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]];

        let (e_full, _, _) = vdw_energy_and_gradient(&coords, 0, 1, &pi, &pj, false);
        let (e_14, _, _) = vdw_energy_and_gradient(&coords, 0, 1, &pi, &pj, true);

        assert!(
            (e_14 - 0.75 * e_full).abs() < 1e-10,
            "1-4 VDW should be 0.75 * full: got {}, expected {}",
            e_14,
            0.75 * e_full
        );
    }

    #[test]
    fn test_vdw_gradient_numerical() {
        let pi = make_params(1.7, 0.07);
        let pj = make_params(1.55, 0.17);
        let coords = [[0.0, 0.0, 0.0], [1.8, 0.0, 0.0]];

        let (_, gi, gj) = vdw_energy_and_gradient(&coords, 0, 1, &pi, &pj, false);

        let eps = 1e-7;
        for (atom_idx, grad) in [(0usize, gi), (1usize, gj)] {
            for dim in 0..3 {
                let mut cp2: Vec<[f64; 3]> = coords.to_vec();
                cp2[atom_idx][dim] += eps;
                let (ep, _, _) = vdw_energy_and_gradient(&cp2, 0, 1, &pi, &pj, false);
                let (e0, _, _) = vdw_energy_and_gradient(&coords, 0, 1, &pi, &pj, false);
                let num = (ep - e0) / eps;
                assert!(
                    (grad[dim] - num).abs() < 1e-4,
                    "VDW grad[{}] = {} vs numerical {} for atom {}",
                    dim,
                    grad[dim],
                    num,
                    atom_idx
                );
            }
        }
    }
}
