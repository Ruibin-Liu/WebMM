//! van der Waals (buffered 14-7 potential) term for MMFF94

use super::atom_types::get_atom_type_props;
use super::MMFFAtomType;

/// van der Waals parameters
#[derive(Debug, Clone, Copy)]
pub struct VDWParams {
    pub r0: f64,
    pub epsilon: f64,
    pub alpha: f64,
    pub beta: f64,
}

pub fn get_vdw_params(atom_type: MMFFAtomType) -> VDWParams {
    match get_atom_type_props(atom_type) {
        Some(props) => VDWParams {
            r0: props.vdw_r,
            epsilon: props.vdw_eps,
            alpha: props.vdw_alpha,
            beta: 2.25,
        },
        None => VDWParams {
            r0: 1.7,
            epsilon: 0.07,
            alpha: 0.083,
            beta: 2.25,
        },
    }
}

/// Calculate VDW energy and gradient
pub fn vdw_energy_and_gradient(
    coords: &[[f64; 3]],
    i: usize,
    j: usize,
    params_i: &VDWParams,
    params_j: &VDWParams,
) -> (f64, [f64; 3], [f64; 3]) {
    let r_vec = [
        coords[j][0] - coords[i][0],
        coords[j][1] - coords[i][1],
        coords[j][2] - coords[i][2],
    ];
    let r = (r_vec[0].powi(2) + r_vec[1].powi(2) + r_vec[2].powi(2)).sqrt();

    if r < 1e-10 {
        return (0.0, [0.0; 3], [0.0; 3]);
    }

    let r0 = 0.5 * (params_i.r0 + params_j.r0);
    let epsilon = (params_i.epsilon * params_j.epsilon).sqrt();

    // Buffered 14-7: E = epsilon * (ratio14 - ratio7 * (1 + ga)) / (1 + ga)
    // ratio14 = repulsive (dominates at short range), ratio7 = attractive
    // gamma = 0.07, ga = gamma * ratio7
    let gamma = 0.07;
    let r_ratio = r0 / r;
    let ratio7 = r_ratio.powi(7);
    let ratio14 = ratio7 * ratio7;
    let ga = gamma * ratio7;
    let one_plus_ga = 1.0 + ga;

    let energy = epsilon * (ratio14 - ratio7 * one_plus_ga) / one_plus_ga;

    // Gradient via numerical differentiation
    let eps = 1e-7;
    let r_p = r + eps;
    let r_ratio_p = r0 / r_p;
    let ratio7_p = r_ratio_p.powi(7);
    let ratio14_p = ratio7_p * ratio7_p;
    let ga_p = gamma * ratio7_p;
    let one_plus_ga_p = 1.0 + ga_p;
    let e_p = epsilon * (ratio14_p - ratio7_p * one_plus_ga_p) / one_plus_ga_p;

    let d_e_dr = (e_p - energy) / eps;

    let grad_i = [
        -d_e_dr * r_vec[0] / r,
        -d_e_dr * r_vec[1] / r,
        -d_e_dr * r_vec[2] / r,
    ];
    let grad_j = [
        d_e_dr * r_vec[0] / r,
        d_e_dr * r_vec[1] / r,
        d_e_dr * r_vec[2] / r,
    ];

    (energy, grad_i, grad_j)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_params(r0: f64, epsilon: f64) -> VDWParams {
        VDWParams {
            r0,
            epsilon,
            alpha: 0.08333,
            beta: 2.25,
        }
    }

    #[test]
    fn test_vdw_energy_has_minimum() {
        let pi = make_params(1.7, 0.07);
        let pj = make_params(1.7, 0.07);
        let coords_close = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let coords_eq = [[0.0, 0.0, 0.0], [1.7, 0.0, 0.0]];
        let coords_far = [[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]];

        let (e_close, _, _) = vdw_energy_and_gradient(&coords_close, 0, 1, &pi, &pj);
        let (e_eq, _, _) = vdw_energy_and_gradient(&coords_eq, 0, 1, &pi, &pj);
        let (e_far, _, _) = vdw_energy_and_gradient(&coords_far, 0, 1, &pi, &pj);

        assert!(
            e_close > 0.0,
            "VDW should be repulsive at short range: {}",
            e_close
        );
        assert!(
            e_eq < 0.0,
            "VDW should have attractive well near r0: {}",
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
    fn test_vdw_gradient_numerical() {
        let pi = make_params(1.7, 0.07);
        let pj = make_params(1.55, 0.17);
        let coords = [[0.0, 0.0, 0.0], [1.8, 0.0, 0.0]];

        let (_, gi, gj) = vdw_energy_and_gradient(&coords, 0, 1, &pi, &pj);

        let eps = 1e-7;
        for (atom_idx, grad) in [(0usize, gi), (1usize, gj)] {
            for dim in 0..3 {
                let mut cp2: Vec<[f64; 3]> = coords.to_vec();
                cp2[atom_idx][dim] += eps;
                let (ep, _, _) = vdw_energy_and_gradient(&cp2, 0, 1, &pi, &pj);
                let (e0, _, _) = vdw_energy_and_gradient(&coords, 0, 1, &pi, &pj);
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
