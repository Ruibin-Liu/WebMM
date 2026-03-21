//! Bond stretching term for MMFF94

use super::MMFFAtomType;
use crate::molecule::BondType;

/// Bond stretching parameters
#[derive(Debug, Clone, Copy)]
pub struct BondParams {
    pub k_bond: f64, // mdyn/Å²
    pub r0: f64,     // Å
}

/// Get bond parameters for atom types
pub fn get_bond_params(
    type1: MMFFAtomType,
    type2: MMFFAtomType,
    bond_type: BondType,
) -> Option<BondParams> {
    // TODO: Load from parameter tables
    // For now, return simple default values
    match (type1, type2, bond_type) {
        // C-C bonds
        (MMFFAtomType::C_3, MMFFAtomType::C_3, BondType::Single) => Some(BondParams {
            k_bond: 4.7,
            r0: 1.526,
        }),
        (MMFFAtomType::C_2, MMFFAtomType::C_3, BondType::Single) => Some(BondParams {
            k_bond: 5.5,
            r0: 1.501,
        }),
        (MMFFAtomType::C_2, MMFFAtomType::C_2, BondType::Double) => Some(BondParams {
            k_bond: 6.7,
            r0: 1.339,
        }),
        (MMFFAtomType::C_AR, MMFFAtomType::C_AR, BondType::Aromatic) => Some(BondParams {
            k_bond: 6.0,
            r0: 1.39,
        }),

        // C-N bonds
        (MMFFAtomType::C_3, MMFFAtomType::N_3, BondType::Single) => Some(BondParams {
            k_bond: 4.7,
            r0: 1.45,
        }),
        (MMFFAtomType::C_2, MMFFAtomType::N_2, BondType::Double) => Some(BondParams {
            k_bond: 7.0,
            r0: 1.27,
        }),

        // C-O bonds
        (MMFFAtomType::C_3, MMFFAtomType::O_3, BondType::Single) => Some(BondParams {
            k_bond: 5.5,
            r0: 1.43,
        }),
        (MMFFAtomType::C_2, MMFFAtomType::O_2, BondType::Double) => Some(BondParams {
            k_bond: 10.5,
            r0: 1.22,
        }),
        (MMFFAtomType::C_2, MMFFAtomType::O_R, BondType::Double) => Some(BondParams {
            k_bond: 10.5,
            r0: 1.23,
        }),

        // C-S bonds
        (MMFFAtomType::C_3, MMFFAtomType::S_3, BondType::Single) => Some(BondParams {
            k_bond: 2.7,
            r0: 1.82,
        }),

        // H-X bonds (symmetric matching)
        (MMFFAtomType::H, MMFFAtomType::C_3, BondType::Single)
        | (MMFFAtomType::C_3, MMFFAtomType::H, BondType::Single) => Some(BondParams {
            k_bond: 4.5,
            r0: 1.113,
        }),
        (MMFFAtomType::H, MMFFAtomType::N_3, BondType::Single)
        | (MMFFAtomType::N_3, MMFFAtomType::H, BondType::Single) => Some(BondParams {
            k_bond: 5.5,
            r0: 1.012,
        }),
        (MMFFAtomType::H, MMFFAtomType::O_3, BondType::Single)
        | (MMFFAtomType::O_3, MMFFAtomType::H, BondType::Single) => Some(BondParams {
            k_bond: 5.5,
            r0: 0.960,
        }),
        (MMFFAtomType::H, MMFFAtomType::S_3, BondType::Single)
        | (MMFFAtomType::S_3, MMFFAtomType::H, BondType::Single) => Some(BondParams {
            k_bond: 4.0,
            r0: 1.336,
        }),

        _ => {
            let t1_name = format!("{:?}", type1);
            let t2_name = format!("{:?}", type2);
            let bt_name = format!("{:?}", bond_type);
            crate::utils::get_bond_params_from_json(&t1_name, &t2_name, &bt_name)
        }
    }
}

/// Calculate bond stretching energy
pub fn bond_energy(coords: &[[f64; 3]], i: usize, j: usize, params: &BondParams) -> f64 {
    let r_vec = [
        coords[j][0] - coords[i][0],
        coords[j][1] - coords[i][1],
        coords[j][2] - coords[i][2],
    ];
    let r = (r_vec[0].powi(2) + r_vec[1].powi(2) + r_vec[2].powi(2)).sqrt();
    let dr = r - params.r0;

    // E_bond = 143.9324 * k_bond * (r - r0)^2  (kcal/mol)
    143.9324 * params.k_bond * dr * dr
}

/// Calculate bond stretching gradient (forces on atoms i and j)
pub fn bond_gradient(
    coords: &[[f64; 3]],
    i: usize,
    j: usize,
    params: &BondParams,
) -> ([f64; 3], [f64; 3]) {
    let r_vec = [
        coords[j][0] - coords[i][0],
        coords[j][1] - coords[i][1],
        coords[j][2] - coords[i][2],
    ];
    let r = (r_vec[0].powi(2) + r_vec[1].powi(2) + r_vec[2].powi(2)).sqrt();
    let dr = r - params.r0;

    if r < 1e-10 {
        return ([0.0; 3], [0.0; 3]);
    }

    // dE/dr = 2 * 143.9324 * k_bond * (r - r0)
    let d_e_dr = 2.0 * 143.9324 * params.k_bond * dr;

    // grad_i = -dE/dr * r_vec / r, grad_j = +dE/dr * r_vec / r
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

    (grad_i, grad_j)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bond_energy() {
        let coords = vec![[0.0, 0.0, 0.0], [1.526, 0.0, 0.0]];
        let params = BondParams {
            k_bond: 4.7,
            r0: 1.526,
        };
        let energy = bond_energy(&coords, 0, 1, &params);
        assert!(energy.is_finite());
        assert!(
            (energy - 0.0).abs() < 1e-10,
            "Energy should be zero at equilibrium"
        );
    }

    #[test]
    fn test_bond_gradient_direction() {
        let coords = vec![[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let params = BondParams {
            k_bond: 4.7,
            r0: 1.526,
        };
        let (gi, gj) = bond_gradient(&coords, 0, 1, &params);

        // Stretched bond: gradient on atom i is negative x (descent moves i toward j)
        assert!(
            gi[0] < 0.0,
            "Stretched bond: grad_i[0] should be negative (descent pulls i toward j)"
        );
        assert!(
            gj[0] > 0.0,
            "Stretched bond: grad_j[0] should be positive (descent pulls j toward i)"
        );
        assert!(
            (gi[1].abs() < 1e-10) && (gi[2].abs() < 1e-10),
            "No y/z force for bond along x"
        );
    }

    #[test]
    fn test_bond_gradient_equilibrium() {
        let coords = vec![[0.0, 0.0, 0.0], [1.526, 0.0, 0.0]];
        let params = BondParams {
            k_bond: 4.7,
            r0: 1.526,
        };
        let (gi, gj) = bond_gradient(&coords, 0, 1, &params);

        for d in 0..3 {
            assert!(
                gi[d].abs() < 1e-10,
                "Gradient should be zero at equilibrium, gi[{}]={}",
                d,
                gi[d]
            );
            assert!(
                gj[d].abs() < 1e-10,
                "Gradient should be zero at equilibrium, gj[{}]={}",
                d,
                gj[d]
            );
        }
    }

    #[test]
    fn test_bond_gradient_numerical() {
        let coords = vec![[0.0, 0.0, 0.0], [1.8, 0.0, 0.0]];
        let params = BondParams {
            k_bond: 4.7,
            r0: 1.526,
        };
        let (gi, gj) = bond_gradient(&coords, 0, 1, &params);

        let eps = 1e-7;
        for (atom_idx, grad) in [(0usize, gi), (1usize, gj)] {
            for dim in 0..3 {
                let mut coords_p = coords.clone();
                coords_p[atom_idx][dim] += eps;
                let e_plus = bond_energy(&coords_p, 0, 1, &params);
                let e_ref = bond_energy(&coords, 0, 1, &params);
                let num_grad = (e_plus - e_ref) / eps;
                assert!(
                    (grad[dim] - num_grad).abs() < 1e-4,
                    "Analytical grad[{}] = {} but numerical = {} for atom {}",
                    dim,
                    grad[dim],
                    num_grad,
                    atom_idx
                );
            }
        }
    }
}
