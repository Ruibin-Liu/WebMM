//! van der Waals (14-7 potential) term for MMFF94

use super::MMFFAtomType;

/// van der Waals parameters
#[derive(Debug, Clone, Copy)]
pub struct VDWParams {
    pub r0: f64,      // Å (vdW radius)
    pub epsilon: f64, // kcal/mol
    pub alpha: f64,   // 1/12 to 1/13 (exponent for r0/r term)
    pub beta: f64,    // Buffered potential constant (typically 2.25)
}

/// Get VDW parameters for atom type
pub fn get_vdw_params(atom_type: MMFFAtomType) -> VDWParams {
    // TODO: Load from parameter tables
    // For now, return simple default values
    match atom_type {
        MMFFAtomType::C_3 | MMFFAtomType::C_2 | MMFFAtomType::C_1 | MMFFAtomType::C_AR => {
            VDWParams {
                r0: 1.7,
                epsilon: 0.07,
                alpha: 0.08333,
                beta: 2.25,
            }
        }
        MMFFAtomType::N_3 | MMFFAtomType::N_2 | MMFFAtomType::N_1 | MMFFAtomType::N_AR => {
            VDWParams {
                r0: 1.55,
                epsilon: 0.17,
                alpha: 0.08333,
                beta: 2.25,
            }
        }
        MMFFAtomType::O_3 | MMFFAtomType::O_2 => VDWParams {
            r0: 1.52,
            epsilon: 0.095,
            alpha: 0.08333,
            beta: 2.25,
        },
        MMFFAtomType::F => VDWParams {
            r0: 1.47,
            epsilon: 0.11,
            alpha: 0.08333,
            beta: 2.25,
        },
        MMFFAtomType::Cl => VDWParams {
            r0: 1.75,
            epsilon: 0.266,
            alpha: 0.08333,
            beta: 2.25,
        },
        MMFFAtomType::Br => VDWParams {
            r0: 1.85,
            epsilon: 0.401,
            alpha: 0.08333,
            beta: 2.25,
        },
        MMFFAtomType::I => VDWParams {
            r0: 1.98,
            epsilon: 0.50,
            alpha: 0.08333,
            beta: 2.25,
        },
        MMFFAtomType::S_3 | MMFFAtomType::S_2 | MMFFAtomType::S_AR => VDWParams {
            r0: 1.8,
            epsilon: 0.30,
            alpha: 0.08333,
            beta: 2.25,
        },
        MMFFAtomType::P_3 | MMFFAtomType::P_4 => VDWParams {
            r0: 1.8,
            epsilon: 0.20,
            alpha: 0.08333,
            beta: 2.25,
        },
        _ => {
            // Default to carbon
            VDWParams {
                r0: 1.7,
                epsilon: 0.07,
                alpha: 0.08333,
                beta: 2.25,
            }
        }
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

    // Combined parameters
    let r0 = 0.5 * (params_i.r0 + params_j.r0);
    let epsilon = (params_i.epsilon * params_j.epsilon).sqrt();
    let _alpha = 0.5 * (params_i.alpha + params_j.alpha);
    let beta = 0.5 * (params_i.beta + params_j.beta);

    let r_ratio = r0 / r;
    let ratio7 = r_ratio.powi(7);
    let ratio12 = r_ratio.powi(12);
    let beta_term = beta * ratio7;

    // Energy: epsilon * [(r0/r)^7 * (1 + beta) - (1 + beta*(r0/r)^7)]
    let energy = epsilon * (ratio7 * (1.0 + beta_term) - (1.0 + beta_term * ratio7));

    // Gradient: dE/dx = (dE/dr) * (dr/dx)
    // dE/dr = epsilon * (r0^7 / r^8) * [(1+beta)*(7-12*beta*(r0/r)^7) + 12*beta*(r0/r)^12]
    let term1 = 7.0 * ratio7 * (1.0 + beta_term);
    let term2 = 12.0 * ratio12 * (1.0 + beta_term);
    let dE_dr = epsilon * (term1 - term2) / r;

    if r < 1e-10 {
        return (energy, [0.0; 3], [0.0; 3]);
    }

    let factor = dE_dr / r;
    let grad_i = [
        -factor * r_vec[0] / r,
        -factor * r_vec[1] / r,
        -factor * r_vec[2] / r,
    ];
    let grad_j = [
        factor * r_vec[0] / r,
        factor * r_vec[1] / r,
        factor * r_vec[2] / r,
    ];

    (energy, grad_i, grad_j)
}
