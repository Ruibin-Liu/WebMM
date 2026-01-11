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
        (MMFFAtomType::C_3, MMFFAtomType::C_3, BondType::Single) => Some(BondParams {
            k_bond: 4.7,
            r0: 1.526,
        }),
        (MMFFAtomType::C_AR, MMFFAtomType::C_AR, BondType::Aromatic) => Some(BondParams {
            k_bond: 6.0,
            r0: 1.39,
        }),
        (MMFFAtomType::C_3, MMFFAtomType::N_3, BondType::Single) => Some(BondParams {
            k_bond: 4.7,
            r0: 1.45,
        }),
        (MMFFAtomType::C_3, MMFFAtomType::O_3, BondType::Single) => Some(BondParams {
            k_bond: 5.5,
            r0: 1.43,
        }),
        (MMFFAtomType::C_3, MMFFAtomType::S_3, BondType::Single) => Some(BondParams {
            k_bond: 2.7,
            r0: 1.82,
        }),
        _ => None,
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

    if r == 0.0 {
        return ([0.0; 3], [0.0; 3]);
    }

    // Gradient: dE/dx = (dE/dr) * (dr/dx) = (dE/dr) * (r_vec_x / r)
    // dE/dr = 2 * 143.9324 * k_bond * (r - r0)
    let dE_dr = 2.0 * 143.9324 * params.k_bond * dr;
    let factor = dE_dr / r;

    // Force on atom i is negative gradient, on atom j is positive
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

    (grad_i, grad_j)
}
