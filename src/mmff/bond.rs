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
    let type1 = super::base_type(type1);
    let type2 = super::base_type(type2);
    match (type1, type2, bond_type) {
        // C-C bonds
        (MMFFAtomType::C_3, MMFFAtomType::C_3, BondType::Single) => Some(BondParams {
            k_bond: 4.258,
            r0: 1.508,
        }),
        (MMFFAtomType::C_2, MMFFAtomType::C_3, BondType::Single)
        | (MMFFAtomType::C_3, MMFFAtomType::C_2, BondType::Single) => Some(BondParams {
            k_bond: 4.19,
            r0: 1.492,
        }),
        (MMFFAtomType::C_2, MMFFAtomType::C_2, BondType::Double) => Some(BondParams {
            k_bond: 9.505,
            r0: 1.333,
        }),
        (MMFFAtomType::C_1, MMFFAtomType::C_1, BondType::Triple) => Some(BondParams {
            k_bond: 15.206,
            r0: 1.2,
        }),
        (MMFFAtomType::C_AR, MMFFAtomType::C_AR, BondType::Aromatic) => Some(BondParams {
            k_bond: 5.573,
            r0: 1.374,
        }),
        (MMFFAtomType::C_3, MMFFAtomType::C_AR, BondType::Single)
        | (MMFFAtomType::C_AR, MMFFAtomType::C_3, BondType::Single) => Some(BondParams {
            k_bond: 4.957,
            r0: 1.486,
        }),
        (MMFFAtomType::C_2, MMFFAtomType::C_AR, BondType::Single)
        | (MMFFAtomType::C_AR, MMFFAtomType::C_2, BondType::Single) => Some(BondParams {
            k_bond: 5.0,
            r0: 1.484,
        }),
        (MMFFAtomType::C_AR, MMFFAtomType::C_AR, BondType::Single) => Some(BondParams {
            k_bond: 5.0,
            r0: 1.484,
        }),

        // C-N bonds
        (MMFFAtomType::C_3, MMFFAtomType::N_3, BondType::Single)
        | (MMFFAtomType::N_3, MMFFAtomType::C_3, BondType::Single) => Some(BondParams {
            k_bond: 4.7,
            r0: 1.45,
        }),
        (MMFFAtomType::C_3, MMFFAtomType::N_2, BondType::Single)
        | (MMFFAtomType::N_2, MMFFAtomType::C_3, BondType::Single) => Some(BondParams {
            k_bond: 5.0,
            r0: 1.42,
        }),
        (MMFFAtomType::C_3, MMFFAtomType::N_AR, BondType::Single)
        | (MMFFAtomType::N_AR, MMFFAtomType::C_3, BondType::Single) => Some(BondParams {
            k_bond: 4.5,
            r0: 1.42,
        }),
        (MMFFAtomType::C_3, MMFFAtomType::N_PL3, BondType::Single)
        | (MMFFAtomType::N_PL3, MMFFAtomType::C_3, BondType::Single) => Some(BondParams {
            k_bond: 4.7,
            r0: 1.45,
        }),
        (MMFFAtomType::C_3, MMFFAtomType::N_AM, BondType::Single)
        | (MMFFAtomType::N_AM, MMFFAtomType::C_3, BondType::Single) => Some(BondParams {
            k_bond: 4.7,
            r0: 1.42,
        }),
        (MMFFAtomType::C_2, MMFFAtomType::N_2, BondType::Double)
        | (MMFFAtomType::N_2, MMFFAtomType::C_2, BondType::Double) => Some(BondParams {
            k_bond: 7.0,
            r0: 1.27,
        }),
        (MMFFAtomType::C_2, MMFFAtomType::N_AR, BondType::Double)
        | (MMFFAtomType::N_AR, MMFFAtomType::C_2, BondType::Double) => Some(BondParams {
            k_bond: 7.0,
            r0: 1.28,
        }),
        (MMFFAtomType::C_AR, MMFFAtomType::N_AR, BondType::Aromatic) => Some(BondParams {
            k_bond: 6.0,
            r0: 1.34,
        }),
        (MMFFAtomType::C_1, MMFFAtomType::N_1, BondType::Triple)
        | (MMFFAtomType::N_1, MMFFAtomType::C_1, BondType::Triple) => Some(BondParams {
            k_bond: 10.0,
            r0: 1.15,
        }),
        (MMFFAtomType::C_AR, MMFFAtomType::N_AM, BondType::Aromatic)
        | (MMFFAtomType::N_AM, MMFFAtomType::C_AR, BondType::Aromatic) => Some(BondParams {
            k_bond: 5.5,
            r0: 1.37,
        }),

        // C-O bonds
        (MMFFAtomType::C_3, MMFFAtomType::O_3, BondType::Single)
        | (MMFFAtomType::O_3, MMFFAtomType::C_3, BondType::Single) => Some(BondParams {
            k_bond: 5.047,
            r0: 1.418,
        }),
        (MMFFAtomType::C_3, MMFFAtomType::O_2, BondType::Single)
        | (MMFFAtomType::O_2, MMFFAtomType::C_3, BondType::Single) => Some(BondParams {
            k_bond: 5.5,
            r0: 1.40,
        }),
        (MMFFAtomType::C_3, MMFFAtomType::O_R, BondType::Single)
        | (MMFFAtomType::O_R, MMFFAtomType::C_3, BondType::Single) => Some(BondParams {
            k_bond: 5.0,
            r0: 1.43,
        }),
        (MMFFAtomType::C_3, MMFFAtomType::O_CO2, BondType::Single)
        | (MMFFAtomType::O_CO2, MMFFAtomType::C_3, BondType::Single) => Some(BondParams {
            k_bond: 5.5,
            r0: 1.40,
        }),
        (MMFFAtomType::C_2, MMFFAtomType::O_2, BondType::Double)
        | (MMFFAtomType::O_2, MMFFAtomType::C_2, BondType::Double) => Some(BondParams {
            k_bond: 12.95,
            r0: 1.222,
        }),
        (MMFFAtomType::C_2, MMFFAtomType::O_R, BondType::Double)
        | (MMFFAtomType::O_R, MMFFAtomType::C_2, BondType::Double) => Some(BondParams {
            k_bond: 10.5,
            r0: 1.23,
        }),
        (MMFFAtomType::C_AR, MMFFAtomType::O_CO2, BondType::Double)
        | (MMFFAtomType::O_CO2, MMFFAtomType::C_AR, BondType::Double) => Some(BondParams {
            k_bond: 10.0,
            r0: 1.23,
        }),
        (MMFFAtomType::C_2, MMFFAtomType::O_CO2, BondType::Double)
        | (MMFFAtomType::O_CO2, MMFFAtomType::C_2, BondType::Double) => Some(BondParams {
            k_bond: 12.95,
            r0: 1.222,
        }),
        (MMFFAtomType::C_2, MMFFAtomType::O_3, BondType::Single)
        | (MMFFAtomType::O_3, MMFFAtomType::C_2, BondType::Single) => Some(BondParams {
            k_bond: 5.801,
            r0: 1.355,
        }),
        (MMFFAtomType::C_2, MMFFAtomType::O_R, BondType::Single)
        | (MMFFAtomType::O_R, MMFFAtomType::C_2, BondType::Single) => Some(BondParams {
            k_bond: 5.801,
            r0: 1.355,
        }),
        (MMFFAtomType::C_AR, MMFFAtomType::O_R, BondType::Aromatic)
        | (MMFFAtomType::O_R, MMFFAtomType::C_AR, BondType::Aromatic) => Some(BondParams {
            k_bond: 5.0,
            r0: 1.37,
        }),
        (MMFFAtomType::C_AR, MMFFAtomType::O_3, BondType::Single)
        | (MMFFAtomType::O_3, MMFFAtomType::C_AR, BondType::Single) => Some(BondParams {
            k_bond: 5.614,
            r0: 1.376,
        }),

        // C-S bonds
        (MMFFAtomType::C_3, MMFFAtomType::S_3, BondType::Single)
        | (MMFFAtomType::S_3, MMFFAtomType::C_3, BondType::Single) => Some(BondParams {
            k_bond: 2.7,
            r0: 1.82,
        }),
        (MMFFAtomType::C_2, MMFFAtomType::S_2, BondType::Double)
        | (MMFFAtomType::S_2, MMFFAtomType::C_2, BondType::Double) => Some(BondParams {
            k_bond: 5.5,
            r0: 1.61,
        }),
        (MMFFAtomType::C_AR, MMFFAtomType::S_AR, BondType::Aromatic)
        | (MMFFAtomType::S_AR, MMFFAtomType::C_AR, BondType::Aromatic) => Some(BondParams {
            k_bond: 4.0,
            r0: 1.71,
        }),

        // C-H bonds (symmetric)
        (MMFFAtomType::H, MMFFAtomType::C_3, BondType::Single)
        | (MMFFAtomType::C_3, MMFFAtomType::H, BondType::Single) => Some(BondParams {
            k_bond: 4.766,
            r0: 1.093,
        }),
        (MMFFAtomType::H, MMFFAtomType::C_2, BondType::Single)
        | (MMFFAtomType::C_2, MMFFAtomType::H, BondType::Single) => Some(BondParams {
            k_bond: 4.65,
            r0: 1.101,
        }),
        (MMFFAtomType::H, MMFFAtomType::C_1, BondType::Single)
        | (MMFFAtomType::C_1, MMFFAtomType::H, BondType::Single) => Some(BondParams {
            k_bond: 5.726,
            r0: 1.065,
        }),
        (MMFFAtomType::H, MMFFAtomType::C_AR, BondType::Single)
        | (MMFFAtomType::C_AR, MMFFAtomType::H, BondType::Single) => Some(BondParams {
            k_bond: 5.306,
            r0: 1.084,
        }),
        (MMFFAtomType::H, MMFFAtomType::C_CAT, BondType::Single)
        | (MMFFAtomType::C_CAT, MMFFAtomType::H, BondType::Single) => Some(BondParams {
            k_bond: 4.766,
            r0: 1.093,
        }),
        (MMFFAtomType::H, MMFFAtomType::C_AN, BondType::Single)
        | (MMFFAtomType::C_AN, MMFFAtomType::H, BondType::Single) => Some(BondParams {
            k_bond: 4.766,
            r0: 1.093,
        }),

        // N-H bonds (symmetric)
        (MMFFAtomType::H, MMFFAtomType::N_3, BondType::Single)
        | (MMFFAtomType::N_3, MMFFAtomType::H, BondType::Single) => Some(BondParams {
            k_bond: 5.5,
            r0: 1.012,
        }),
        (MMFFAtomType::H, MMFFAtomType::N_2, BondType::Single)
        | (MMFFAtomType::N_2, MMFFAtomType::H, BondType::Single) => Some(BondParams {
            k_bond: 5.5,
            r0: 1.000,
        }),
        (MMFFAtomType::H, MMFFAtomType::N_AR, BondType::Single)
        | (MMFFAtomType::N_AR, MMFFAtomType::H, BondType::Single) => Some(BondParams {
            k_bond: 5.0,
            r0: 1.010,
        }),
        (MMFFAtomType::H, MMFFAtomType::N_PL3, BondType::Single)
        | (MMFFAtomType::N_PL3, MMFFAtomType::H, BondType::Single) => Some(BondParams {
            k_bond: 5.5,
            r0: 1.012,
        }),
        (MMFFAtomType::H, MMFFAtomType::N_AM, BondType::Single)
        | (MMFFAtomType::N_AM, MMFFAtomType::H, BondType::Single) => Some(BondParams {
            k_bond: 5.0,
            r0: 1.010,
        }),
        (MMFFAtomType::H, MMFFAtomType::N_4, BondType::Single)
        | (MMFFAtomType::N_4, MMFFAtomType::H, BondType::Single) => Some(BondParams {
            k_bond: 5.5,
            r0: 1.012,
        }),

        // O-H bonds (symmetric)
        (MMFFAtomType::H, MMFFAtomType::O_3, BondType::Single)
        | (MMFFAtomType::O_3, MMFFAtomType::H, BondType::Single) => Some(BondParams {
            k_bond: 7.794,
            r0: 0.972,
        }),
        (MMFFAtomType::H, MMFFAtomType::O_2, BondType::Single)
        | (MMFFAtomType::O_2, MMFFAtomType::H, BondType::Single) => Some(BondParams {
            k_bond: 7.794,
            r0: 0.972,
        }),
        (MMFFAtomType::H, MMFFAtomType::O_R, BondType::Single)
        | (MMFFAtomType::O_R, MMFFAtomType::H, BondType::Single) => Some(BondParams {
            k_bond: 7.794,
            r0: 0.972,
        }),

        // S-H bonds (symmetric)
        (MMFFAtomType::H, MMFFAtomType::S_3, BondType::Single)
        | (MMFFAtomType::S_3, MMFFAtomType::H, BondType::Single) => Some(BondParams {
            k_bond: 4.0,
            r0: 1.336,
        }),
        (MMFFAtomType::H, MMFFAtomType::S_2, BondType::Single)
        | (MMFFAtomType::S_2, MMFFAtomType::H, BondType::Single) => Some(BondParams {
            k_bond: 4.0,
            r0: 1.336,
        }),

        // Halogen bonds (symmetric)
        (MMFFAtomType::C_3, MMFFAtomType::F, BondType::Single)
        | (MMFFAtomType::F, MMFFAtomType::C_3, BondType::Single) => Some(BondParams {
            k_bond: 6.0,
            r0: 1.38,
        }),
        (MMFFAtomType::C_3, MMFFAtomType::Cl, BondType::Single)
        | (MMFFAtomType::Cl, MMFFAtomType::C_3, BondType::Single) => Some(BondParams {
            k_bond: 2.893,
            r0: 1.805,
        }),
        (MMFFAtomType::C_3, MMFFAtomType::Br, BondType::Single)
        | (MMFFAtomType::Br, MMFFAtomType::C_3, BondType::Single) => Some(BondParams {
            k_bond: 3.0,
            r0: 1.94,
        }),
        (MMFFAtomType::C_3, MMFFAtomType::I, BondType::Single)
        | (MMFFAtomType::I, MMFFAtomType::C_3, BondType::Single) => Some(BondParams {
            k_bond: 2.5,
            r0: 2.14,
        }),
        (MMFFAtomType::C_AR, MMFFAtomType::F, BondType::Single)
        | (MMFFAtomType::F, MMFFAtomType::C_AR, BondType::Single) => Some(BondParams {
            k_bond: 5.5,
            r0: 1.33,
        }),
        (MMFFAtomType::C_AR, MMFFAtomType::Cl, BondType::Single)
        | (MMFFAtomType::Cl, MMFFAtomType::C_AR, BondType::Single) => Some(BondParams {
            k_bond: 3.5,
            r0: 1.72,
        }),
        (MMFFAtomType::C_AR, MMFFAtomType::Br, BondType::Single)
        | (MMFFAtomType::Br, MMFFAtomType::C_AR, BondType::Single) => Some(BondParams {
            k_bond: 3.0,
            r0: 1.87,
        }),
        (MMFFAtomType::C_2, MMFFAtomType::F, BondType::Single)
        | (MMFFAtomType::F, MMFFAtomType::C_2, BondType::Single) => Some(BondParams {
            k_bond: 6.0,
            r0: 1.30,
        }),
        (MMFFAtomType::C_2, MMFFAtomType::Cl, BondType::Single)
        | (MMFFAtomType::Cl, MMFFAtomType::C_2, BondType::Single) => Some(BondParams {
            k_bond: 3.5,
            r0: 1.72,
        }),

        // N-N bonds
        (MMFFAtomType::N_3, MMFFAtomType::N_3, BondType::Single) => Some(BondParams {
            k_bond: 3.5,
            r0: 1.45,
        }),
        (MMFFAtomType::N_2, MMFFAtomType::N_2, BondType::Double) => Some(BondParams {
            k_bond: 6.0,
            r0: 1.25,
        }),
        (MMFFAtomType::N_3, MMFFAtomType::N_AR, BondType::Single)
        | (MMFFAtomType::N_AR, MMFFAtomType::N_3, BondType::Single) => Some(BondParams {
            k_bond: 4.0,
            r0: 1.40,
        }),
        (MMFFAtomType::N_AR, MMFFAtomType::N_AR, BondType::Aromatic) => Some(BondParams {
            k_bond: 5.5,
            r0: 1.33,
        }),
        (MMFFAtomType::N_3, MMFFAtomType::C_2, BondType::Single)
        | (MMFFAtomType::C_2, MMFFAtomType::N_3, BondType::Single) => Some(BondParams {
            k_bond: 5.0,
            r0: 1.42,
        }),

        // O-O bonds
        (MMFFAtomType::O_3, MMFFAtomType::O_3, BondType::Single) => Some(BondParams {
            k_bond: 4.0,
            r0: 1.48,
        }),
        (MMFFAtomType::O_3, MMFFAtomType::O_2, BondType::Single)
        | (MMFFAtomType::O_2, MMFFAtomType::O_3, BondType::Single) => Some(BondParams {
            k_bond: 4.5,
            r0: 1.45,
        }),

        // P bonds
        (MMFFAtomType::C_3, MMFFAtomType::P_3, BondType::Single)
        | (MMFFAtomType::P_3, MMFFAtomType::C_3, BondType::Single) => Some(BondParams {
            k_bond: 3.0,
            r0: 1.87,
        }),
        (MMFFAtomType::P_3, MMFFAtomType::O_3, BondType::Single)
        | (MMFFAtomType::O_3, MMFFAtomType::P_3, BondType::Single) => Some(BondParams {
            k_bond: 4.0,
            r0: 1.60,
        }),
        (MMFFAtomType::P_4, MMFFAtomType::O_2, BondType::Double)
        | (MMFFAtomType::O_2, MMFFAtomType::P_4, BondType::Double) => Some(BondParams {
            k_bond: 8.5,
            r0: 1.48,
        }),

        _ => {
            if let Some((kb, r0)) = super::estimation::estimate_bond_params(type1, type2, bond_type)
            {
                Some(BondParams { k_bond: kb, r0 })
            } else {
                let t1_name = format!("{:?}", type1);
                let t2_name = format!("{:?}", type2);
                let bt_name = format!("{:?}", bond_type);
                crate::utils::get_bond_params_from_json(&t1_name, &t2_name, &bt_name)
            }
        }
    }
}

/// Calculate bond stretching energy
///
/// MMFF94 anharmonic bond stretch (RDKit-compatible):
///   E = 0.5 * c1 * kb * dr² * (1 + cs * dr + c3 * cs² * dr²)
/// where cs = -2.0, c3 = 7/12, c1 = 143.9324
pub fn bond_energy(coords: &[[f64; 3]], i: usize, j: usize, params: &BondParams) -> f64 {
    let r_vec = [
        coords[j][0] - coords[i][0],
        coords[j][1] - coords[i][1],
        coords[j][2] - coords[i][2],
    ];
    let r = (r_vec[0].powi(2) + r_vec[1].powi(2) + r_vec[2].powi(2)).sqrt();
    let dr = r - params.r0;

    // RDKit anharmonic bond stretch
    let c1 = 143.9324;
    let cs = -2.0;
    let c3 = 7.0 / 12.0;
    let dr2 = dr * dr;

    c1 * params.k_bond * dr2 * (1.0 + cs * dr + c3 * cs * cs * dr2) / 2.0
}

/// Calculate bond stretching gradient (forces on atoms i and j)
///
/// Uses numerical differentiation for the anharmonic term.
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

    // dE/dr for anharmonic bond:
    // E = 0.5 * c1 * kb * dr² * (1 + cs * dr + c3 * cs² * dr²)
    // dE/dr = c1 * kb * dr * (1 + 1.5 * cs * dr + 2.0 * c3 * cs² * dr²)
    let c1 = 143.9324;
    let cs = -2.0;
    let c3 = 7.0 / 12.0;

    let d_e_dr = c1 * params.k_bond * dr * (1.0 + 1.5 * cs * dr + 2.0 * c3 * cs * cs * dr * dr);

    // grad_i = -dE/dr * r_vec / r, grad_j = +dE/dr * r_vec / r
    // (gradient points in direction of increasing energy, so descent moves opposite)
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
