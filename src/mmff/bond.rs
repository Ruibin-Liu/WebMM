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
    if let Some(p) = lookup_bond_params_exact(type1, type2, bond_type) {
        return Some(p);
    }
    let base1 = super::base_type(type1);
    let base2 = super::base_type(type2);
    if base1 != type1 || base2 != type2 {
        if let Some(p) = lookup_bond_params_exact(base1, base2, bond_type) {
            return Some(p);
        }
    }
    if let Some((kb, r0)) = super::estimation::estimate_bond_params(base1, base2, bond_type) {
        return Some(BondParams { k_bond: kb, r0 });
    }
    let t1_name = format!("{:?}", base1);
    let t2_name = format!("{:?}", base2);
    let bt_name = format!("{:?}", bond_type);
    crate::utils::get_bond_params_from_json(&t1_name, &t2_name, &bt_name)
}

fn lookup_bond_params_exact(
    type1: MMFFAtomType,
    type2: MMFFAtomType,
    bond_type: BondType,
) -> Option<BondParams> {
    match (type1, type2, bond_type) {
        // C-C bonds
        (MMFFAtomType::C_3, MMFFAtomType::C_3, BondType::Single) => Some(BondParams {
            k_bond: 4.258,
            r0: 1.508,
        }),
        (MMFFAtomType::C_2, MMFFAtomType::C_2, BondType::Single) => Some(BondParams {
            k_bond: 4.418,
            r0: 1.489,
        }),
        (MMFFAtomType::C_VIN, MMFFAtomType::C_VIN, BondType::Single) => Some(BondParams {
            k_bond: 5.310,
            r0: 1.430,
        }),
        (MMFFAtomType::C_VIN, MMFFAtomType::N_2, BondType::Single)
        | (MMFFAtomType::N_2, MMFFAtomType::C_VIN, BondType::Single) => Some(BondParams {
            k_bond: 6.385,
            r0: 1.360,
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
        (MMFFAtomType::C_3, MMFFAtomType::C5B, BondType::Single)
        | (MMFFAtomType::C5B, MMFFAtomType::C_3, BondType::Single) => Some(BondParams {
            k_bond: 4.518,
            r0: 1.469,
        }),
        (MMFFAtomType::C_2, MMFFAtomType::C_AR, BondType::Single)
        | (MMFFAtomType::C_AR, MMFFAtomType::C_2, BondType::Single) => Some(BondParams {
            k_bond: 4.488,
            r0: 1.457,
        }),
        (MMFFAtomType::C_AR, MMFFAtomType::C_AR, BondType::Single) => Some(BondParams {
            k_bond: 5.0,
            r0: 1.484,
        }),

        // Vinyl (alkene, MMFF type 2) bonds — RDKit-extracted values
        (MMFFAtomType::C_VIN, MMFFAtomType::C_VIN, BondType::Double) => Some(BondParams {
            k_bond: 9.505,
            r0: 1.333,
        }),
        (MMFFAtomType::C_VIN, MMFFAtomType::C_3, BondType::Single)
        | (MMFFAtomType::C_3, MMFFAtomType::C_VIN, BondType::Single) => Some(BondParams {
            k_bond: 4.539,
            r0: 1.482,
        }),
        (MMFFAtomType::C_VIN, MMFFAtomType::C_2, BondType::Single)
        | (MMFFAtomType::C_2, MMFFAtomType::C_VIN, BondType::Single) => Some(BondParams {
            k_bond: 4.565,
            r0: 1.468,
        }),
        (MMFFAtomType::C_VIN, MMFFAtomType::C_AR, BondType::Single)
        | (MMFFAtomType::C_AR, MMFFAtomType::C_VIN, BondType::Single) => Some(BondParams {
            k_bond: 5.007,
            r0: 1.449,
        }),
        (MMFFAtomType::C_VIN, MMFFAtomType::H, BondType::Single)
        | (MMFFAtomType::H, MMFFAtomType::C_VIN, BondType::Single) => Some(BondParams {
            k_bond: 5.17,
            r0: 1.083,
        }),
        (MMFFAtomType::C_VIN, MMFFAtomType::O_R, BondType::Single)
        | (MMFFAtomType::O_R, MMFFAtomType::C_VIN, BondType::Single) => Some(BondParams {
            k_bond: 5.52,
            r0: 1.373,
        }),

        // Oxidized sulfur bonds — RDKit-extracted values
        (MMFFAtomType::S_OX, MMFFAtomType::C_3, BondType::Single)
        | (MMFFAtomType::C_3, MMFFAtomType::S_OX, BondType::Single) => Some(BondParams {
            k_bond: 2.841,
            r0: 1.813,
        }),
        (MMFFAtomType::S_OX, MMFFAtomType::C_AR, BondType::Single)
        | (MMFFAtomType::C_AR, MMFFAtomType::S_OX, BondType::Single) => Some(BondParams {
            k_bond: 3.098,
            r0: 1.787,
        }),
        (MMFFAtomType::S_OX, MMFFAtomType::O_2, BondType::Double)
        | (MMFFAtomType::O_2, MMFFAtomType::S_OX, BondType::Double) => Some(BondParams {
            k_bond: 8.77,
            r0: 1.5,
        }),
        (MMFFAtomType::S_O2, MMFFAtomType::C_3, BondType::Single)
        | (MMFFAtomType::C_3, MMFFAtomType::S_O2, BondType::Single) => Some(BondParams {
            k_bond: 3.258,
            r0: 1.772,
        }),
        (MMFFAtomType::S_O2, MMFFAtomType::C_AR, BondType::Single)
        | (MMFFAtomType::C_AR, MMFFAtomType::S_O2, BondType::Single) => Some(BondParams {
            k_bond: 3.281,
            r0: 1.77,
        }),
        (MMFFAtomType::S_O2, MMFFAtomType::O_CO2, BondType::Double)
        | (MMFFAtomType::O_CO2, MMFFAtomType::S_O2, BondType::Double) => Some(BondParams {
            k_bond: 10.748,
            r0: 1.45,
        }),
        (MMFFAtomType::S_O2, MMFFAtomType::N_SO2, BondType::Single)
        | (MMFFAtomType::N_SO2, MMFFAtomType::S_O2, BondType::Single) => Some(BondParams {
            k_bond: 3.301,
            r0: 1.71,
        }),
        (MMFFAtomType::S_O2, MMFFAtomType::O_3, BondType::Single)
        | (MMFFAtomType::O_3, MMFFAtomType::S_O2, BondType::Single) => Some(BondParams {
            k_bond: 5.326,
            r0: 1.630,
        }),
        (MMFFAtomType::N_SO2, MMFFAtomType::C_3, BondType::Single)
        | (MMFFAtomType::C_3, MMFFAtomType::N_SO2, BondType::Single) => Some(BondParams {
            k_bond: 3.971,
            r0: 1.472,
        }),

        // Nitro group bonds — RDKit-extracted values
        // (N_NO2-O_CO2 params are identical for Single and Double in RDKit)
        (MMFFAtomType::N_NO2, MMFFAtomType::C_3, BondType::Single)
        | (MMFFAtomType::C_3, MMFFAtomType::N_NO2, BondType::Single) => Some(BondParams {
            k_bond: 3.844,
            r0: 1.48,
        }),
        (MMFFAtomType::N_NO2, MMFFAtomType::C_AR, BondType::Single)
        | (MMFFAtomType::C_AR, MMFFAtomType::N_NO2, BondType::Single) => Some(BondParams {
            k_bond: 4.705,
            r0: 1.431,
        }),
        (MMFFAtomType::N_NO2, MMFFAtomType::O_CO2, BondType::Double)
        | (MMFFAtomType::O_CO2, MMFFAtomType::N_NO2, BondType::Double)
        | (MMFFAtomType::N_NO2, MMFFAtomType::O_CO2, BondType::Single)
        | (MMFFAtomType::O_CO2, MMFFAtomType::N_NO2, BondType::Single) => Some(BondParams {
            k_bond: 9.42,
            r0: 1.233,
        }),

        // Carboxylate bonds — RDKit-extracted values (MMFF94s)
        // (1,41) C_3-CO2M, (32,41) O_CO2-CO2M
        (MMFFAtomType::C_3, MMFFAtomType::C_CO2, BondType::Single)
        | (MMFFAtomType::C_CO2, MMFFAtomType::C_3, BondType::Single) => Some(BondParams {
            k_bond: 3.830,
            r0: 1.510,
        }),
        (MMFFAtomType::C_AR, MMFFAtomType::C_CO2, BondType::Single)
        | (MMFFAtomType::C_CO2, MMFFAtomType::C_AR, BondType::Single) => Some(BondParams {
            k_bond: 4.537,
            r0: 1.468,
        }),
        (MMFFAtomType::O_CO2, MMFFAtomType::C_CO2, BondType::Double)
        | (MMFFAtomType::C_CO2, MMFFAtomType::O_CO2, BondType::Double)
        | (MMFFAtomType::O_CO2, MMFFAtomType::C_CO2, BondType::Single)
        | (MMFFAtomType::C_CO2, MMFFAtomType::O_CO2, BondType::Single) => Some(BondParams {
            k_bond: 9.756,
            r0: 1.261,
        }),

        // Pyridine N-oxide bonds — RDKit-extracted values (MMFF94s)
        // (37,69) C_AR-NPOX, (32,69) O_CO2-NPOX
        (MMFFAtomType::C_AR, MMFFAtomType::N_POX, BondType::Aromatic)
        | (MMFFAtomType::N_POX, MMFFAtomType::C_AR, BondType::Aromatic) => Some(BondParams {
            k_bond: 5.396,
            r0: 1.352,
        }),
        (MMFFAtomType::O_CO2, MMFFAtomType::N_POX, BondType::Single)
        | (MMFFAtomType::N_POX, MMFFAtomType::O_CO2, BondType::Single) => Some(BondParams {
            k_bond: 6.098,
            r0: 1.261,
        }),

        // Silicon bonds — RDKit-extracted values (MMFF94s)
        // (1,19) C_3-Si, (6,19) O_3-Si
        (MMFFAtomType::C_3, MMFFAtomType::Si, BondType::Single)
        | (MMFFAtomType::Si, MMFFAtomType::C_3, BondType::Single) => Some(BondParams {
            k_bond: 2.866,
            r0: 1.830,
        }),
        (MMFFAtomType::O_3, MMFFAtomType::Si, BondType::Single)
        | (MMFFAtomType::Si, MMFFAtomType::O_3, BondType::Single) => Some(BondParams {
            k_bond: 4.661,
            r0: 1.660,
        }),

        // Imine N-H bond — RDKit-extracted value (MMFF94s)
        // (9,27) N_2-H_NIM
        (MMFFAtomType::N_2, MMFFAtomType::H_NIM, BondType::Single)
        | (MMFFAtomType::H_NIM, MMFFAtomType::N_2, BondType::Single) => Some(BondParams {
            k_bond: 6.230,
            r0: 1.026,
        }),

        // C-N bonds
        (MMFFAtomType::C_3, MMFFAtomType::N_3, BondType::Single)
        | (MMFFAtomType::N_3, MMFFAtomType::C_3, BondType::Single) => Some(BondParams {
            k_bond: 5.084,
            r0: 1.451,
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
        (MMFFAtomType::C_3, MMFFAtomType::NPYL, BondType::Single)
        | (MMFFAtomType::NPYL, MMFFAtomType::C_3, BondType::Single) => Some(BondParams {
            k_bond: 6.114,
            r0: 1.445,
        }),
        (MMFFAtomType::C_3, MMFFAtomType::N_PL3, BondType::Single)
        | (MMFFAtomType::N_PL3, MMFFAtomType::C_3, BondType::Single) => Some(BondParams {
            k_bond: 4.7,
            r0: 1.45,
        }),
        (MMFFAtomType::C_AR, MMFFAtomType::N_PL3, BondType::Single)
        | (MMFFAtomType::N_PL3, MMFFAtomType::C_AR, BondType::Single) => Some(BondParams {
            k_bond: 6.168,
            r0: 1.398,
        }),
        (MMFFAtomType::C_3, MMFFAtomType::N_AM, BondType::Single)
        | (MMFFAtomType::N_AM, MMFFAtomType::C_3, BondType::Single) => Some(BondParams {
            k_bond: 4.664,
            r0: 1.436,
        }),
        (MMFFAtomType::C_2, MMFFAtomType::N_2, BondType::Double)
        | (MMFFAtomType::N_2, MMFFAtomType::C_2, BondType::Double) => Some(BondParams {
            k_bond: 10.077,
            r0: 1.290,
        }),
        (MMFFAtomType::C_2, MMFFAtomType::C_1, BondType::Double)
        | (MMFFAtomType::C_1, MMFFAtomType::C_2, BondType::Double) => Some(BondParams {
            k_bond: 9.538,
            r0: 1.297,
        }),
        (MMFFAtomType::C_VIN, MMFFAtomType::C_1, BondType::Double)
        | (MMFFAtomType::C_1, MMFFAtomType::C_VIN, BondType::Double) => Some(BondParams {
            k_bond: 9.538,
            r0: 1.297,
        }),
        (MMFFAtomType::C_2, MMFFAtomType::N_2, BondType::Single)
        | (MMFFAtomType::N_2, MMFFAtomType::C_2, BondType::Single) => Some(BondParams {
            k_bond: 6.273,
            r0: 1.364,
        }),
        (MMFFAtomType::C_VIN, MMFFAtomType::N_AM, BondType::Single)
        | (MMFFAtomType::N_AM, MMFFAtomType::C_VIN, BondType::Single) => Some(BondParams {
            k_bond: 6.329,
            r0: 1.362,
        }),
        (MMFFAtomType::C_2, MMFFAtomType::N_AR, BondType::Double)
        | (MMFFAtomType::N_AR, MMFFAtomType::C_2, BondType::Double) => Some(BondParams {
            k_bond: 7.0,
            r0: 1.28,
        }),
        (MMFFAtomType::C_AR, MMFFAtomType::N_AR, BondType::Aromatic)
        | (MMFFAtomType::N_AR, MMFFAtomType::C_AR, BondType::Aromatic) => Some(BondParams {
            k_bond: 5.737,
            r0: 1.333,
        }),
        (MMFFAtomType::N_AR, MMFFAtomType::C5A, BondType::Aromatic)
        | (MMFFAtomType::C5A, MMFFAtomType::N_AR, BondType::Aromatic) => Some(BondParams {
            k_bond: 7.299,
            r0: 1.330,
        }),
        (MMFFAtomType::C_1, MMFFAtomType::N_1, BondType::Triple)
        | (MMFFAtomType::N_1, MMFFAtomType::C_1, BondType::Triple) => Some(BondParams {
            k_bond: 16.582,
            r0: 1.160,
        }),
        (MMFFAtomType::C_3, MMFFAtomType::C_1, BondType::Single)
        | (MMFFAtomType::C_1, MMFFAtomType::C_3, BondType::Single) => Some(BondParams {
            k_bond: 4.707,
            r0: 1.459,
        }),
        // 3-ring (cyclopropane/epoxide/aziridine) specific bond params (CR3R=22)
        (MMFFAtomType::CR3R, MMFFAtomType::CR3R, BondType::Single) => Some(BondParams {
            k_bond: 3.969,
            r0: 1.499,
        }),
        (MMFFAtomType::CR3R, MMFFAtomType::O_3, BondType::Single)
        | (MMFFAtomType::O_3, MMFFAtomType::CR3R, BondType::Single) => Some(BondParams {
            k_bond: 4.556,
            r0: 1.433,
        }),
        (MMFFAtomType::CR3R, MMFFAtomType::O_R, BondType::Single)
        | (MMFFAtomType::O_R, MMFFAtomType::CR3R, BondType::Single) => Some(BondParams {
            k_bond: 4.556,
            r0: 1.433,
        }),
        (MMFFAtomType::CR3R, MMFFAtomType::H, BondType::Single)
        | (MMFFAtomType::H, MMFFAtomType::CR3R, BondType::Single) => Some(BondParams {
            k_bond: 5.191,
            r0: 1.082,
        }),
        (MMFFAtomType::C_2, MMFFAtomType::CR3R, BondType::Single)
        | (MMFFAtomType::CR3R, MMFFAtomType::C_2, BondType::Single) => Some(BondParams {
            k_bond: 4.926,
            r0: 1.448,
        }),
        (MMFFAtomType::C_VIN, MMFFAtomType::CR3R, BondType::Single)
        | (MMFFAtomType::CR3R, MMFFAtomType::C_VIN, BondType::Single) => Some(BondParams {
            k_bond: 4.926,
            r0: 1.448,
        }),
        // 4-ring (cyclobutane) specific bond params (CR4R=20)
        (MMFFAtomType::CR4R, MMFFAtomType::CR4R, BondType::Single) => Some(BondParams {
            k_bond: 3.663,
            r0: 1.526,
        }),
        (MMFFAtomType::CR4R, MMFFAtomType::H, BondType::Single)
        | (MMFFAtomType::H, MMFFAtomType::CR4R, BondType::Single) => Some(BondParams {
            k_bond: 4.852,
            r0: 1.093,
        }),
        (MMFFAtomType::CR3R, MMFFAtomType::N_3, BondType::Single)
        | (MMFFAtomType::N_3, MMFFAtomType::CR3R, BondType::Single) => Some(BondParams {
            k_bond: 4.223,
            r0: 1.457,
        }),
        (MMFFAtomType::C_AR, MMFFAtomType::N_AM, BondType::Aromatic)
        | (MMFFAtomType::N_AM, MMFFAtomType::C_AR, BondType::Aromatic) => Some(BondParams {
            k_bond: 5.5,
            r0: 1.37,
        }),
        (MMFFAtomType::C_AR, MMFFAtomType::N_AM, BondType::Single)
        | (MMFFAtomType::N_AM, MMFFAtomType::C_AR, BondType::Single) => Some(BondParams {
            k_bond: 5.482,
            r0: 1.395,
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
            k_bond: 5.047,
            r0: 1.418,
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
            k_bond: 2.893,
            r0: 1.805,
        }),
        (MMFFAtomType::S_3, MMFFAtomType::S_3, BondType::Single) => Some(BondParams {
            k_bond: 2.531,
            r0: 2.050,
        }),
        (MMFFAtomType::C_2, MMFFAtomType::S_2, BondType::Double)
        | (MMFFAtomType::S_2, MMFFAtomType::C_2, BondType::Double) => Some(BondParams {
            k_bond: 4.735,
            r0: 1.665,
        }),
        (MMFFAtomType::C_AR, MMFFAtomType::S_AR, BondType::Aromatic)
        | (MMFFAtomType::S_AR, MMFFAtomType::C_AR, BondType::Aromatic) => Some(BondParams {
            k_bond: 4.0,
            r0: 1.71,
        }),
        // 5-ring heteroaromatic specific bond params (C5A=63, C5B=64) from RDKit verbose
        (MMFFAtomType::C5A, MMFFAtomType::S_AR, BondType::Aromatic)
        | (MMFFAtomType::S_AR, MMFFAtomType::C5A, BondType::Aromatic) => Some(BondParams {
            k_bond: 3.589,
            r0: 1.717,
        }),
        (MMFFAtomType::C5B, MMFFAtomType::C5B, BondType::Aromatic) => Some(BondParams {
            k_bond: 4.313,
            r0: 1.418,
        }),
        (MMFFAtomType::NPYL, MMFFAtomType::N5A, BondType::Aromatic)
        | (MMFFAtomType::N5A, MMFFAtomType::NPYL, BondType::Aromatic) => Some(BondParams {
            k_bond: 5.513,
            r0: 1.339,
        }),
        (MMFFAtomType::N5A, MMFFAtomType::C5B, BondType::Aromatic)
        | (MMFFAtomType::C5B, MMFFAtomType::N5A, BondType::Aromatic) => Some(BondParams {
            k_bond: 8.258,
            r0: 1.335,
        }),
        (MMFFAtomType::C5A, MMFFAtomType::OFUR, BondType::Aromatic)
        | (MMFFAtomType::OFUR, MMFFAtomType::C5A, BondType::Aromatic) => Some(BondParams {
            k_bond: 5.787,
            r0: 1.360,
        }),
        (MMFFAtomType::C5A, MMFFAtomType::H, BondType::Single)
        | (MMFFAtomType::H, MMFFAtomType::C5A, BondType::Single) => Some(BondParams {
            k_bond: 5.531,
            r0: 1.080,
        }),
        (MMFFAtomType::C5B, MMFFAtomType::H, BondType::Single)
        | (MMFFAtomType::H, MMFFAtomType::C5B, BondType::Single) => Some(BondParams {
            k_bond: 5.506,
            r0: 1.080,
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
            k_bond: 6.490,
            r0: 1.019,
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
            k_bond: 6.576,
            r0: 1.018,
        }),
        (MMFFAtomType::H, MMFFAtomType::N_AM, BondType::Single)
        | (MMFFAtomType::N_AM, MMFFAtomType::H, BondType::Single) => Some(BondParams {
            k_bond: 6.663,
            r0: 1.015,
        }),
        (MMFFAtomType::N_SO2, MMFFAtomType::H, BondType::Single)
        | (MMFFAtomType::H, MMFFAtomType::N_SO2, BondType::Single) => Some(BondParams {
            k_bond: 6.265,
            r0: 1.028,
        }),
        (MMFFAtomType::H, MMFFAtomType::N_4, BondType::Single)
        | (MMFFAtomType::N_4, MMFFAtomType::H, BondType::Single) => Some(BondParams {
            k_bond: 6.163,
            r0: 1.028,
        }),
        (MMFFAtomType::C_3, MMFFAtomType::N_4, BondType::Single)
        | (MMFFAtomType::N_4, MMFFAtomType::C_3, BondType::Single) => Some(BondParams {
            k_bond: 3.844,
            r0: 1.480,
        }),

        // O-H bonds (symmetric)
        (MMFFAtomType::H_OH, MMFFAtomType::OH2, BondType::Single)
        | (MMFFAtomType::OH2, MMFFAtomType::H_OH, BondType::Single) => Some(BondParams {
            k_bond: 7.88,
            r0: 0.969,
        }),
        (MMFFAtomType::H_ONC, MMFFAtomType::O_3, BondType::Single)
        | (MMFFAtomType::O_3, MMFFAtomType::H_ONC, BondType::Single) => Some(BondParams {
            k_bond: 7.794,
            r0: 0.972,
        }),
        (MMFFAtomType::H_OH, MMFFAtomType::O_3, BondType::Single)
        | (MMFFAtomType::O_3, MMFFAtomType::H_OH, BondType::Single) => Some(BondParams {
            k_bond: 7.88,
            r0: 0.969,
        }),
        (MMFFAtomType::H_COOH, MMFFAtomType::O_3, BondType::Single)
        | (MMFFAtomType::O_3, MMFFAtomType::H_COOH, BondType::Single) => Some(BondParams {
            k_bond: 7.403,
            r0: 0.981,
        }),
        (MMFFAtomType::H_COOH, MMFFAtomType::O_R, BondType::Single)
        | (MMFFAtomType::O_R, MMFFAtomType::H_COOH, BondType::Single) => Some(BondParams {
            k_bond: 7.403,
            r0: 0.981,
        }),
        (MMFFAtomType::HOS, MMFFAtomType::O_3, BondType::Single)
        | (MMFFAtomType::O_3, MMFFAtomType::HOS, BondType::Single) => Some(BondParams {
            k_bond: 7.143,
            r0: 0.986,
        }),
        (MMFFAtomType::HOS, MMFFAtomType::O_R, BondType::Single)
        | (MMFFAtomType::O_R, MMFFAtomType::HOS, BondType::Single) => Some(BondParams {
            k_bond: 7.143,
            r0: 0.986,
        }),
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
            k_bond: 4.014,
            r0: 1.341,
        }),
        (MMFFAtomType::H, MMFFAtomType::S_2, BondType::Single)
        | (MMFFAtomType::S_2, MMFFAtomType::H, BondType::Single) => Some(BondParams {
            k_bond: 4.0,
            r0: 1.336,
        }),

        // Halogen bonds (symmetric)
        (MMFFAtomType::C_3, MMFFAtomType::F, BondType::Single)
        | (MMFFAtomType::F, MMFFAtomType::C_3, BondType::Single) => Some(BondParams {
            k_bond: 6.011,
            r0: 1.360,
        }),
        (MMFFAtomType::C_3, MMFFAtomType::Cl, BondType::Single)
        | (MMFFAtomType::Cl, MMFFAtomType::C_3, BondType::Single) => Some(BondParams {
            k_bond: 2.974,
            r0: 1.773,
        }),
        (MMFFAtomType::C_VIN, MMFFAtomType::Cl, BondType::Single)
        | (MMFFAtomType::Cl, MMFFAtomType::C_VIN, BondType::Single) => Some(BondParams {
            k_bond: 3.390,
            r0: 1.720,
        }),
        (MMFFAtomType::C_2, MMFFAtomType::Cl, BondType::Single)
        | (MMFFAtomType::Cl, MMFFAtomType::C_2, BondType::Single) => Some(BondParams {
            k_bond: 3.449,
            r0: 1.715,
        }),
        (MMFFAtomType::C_3, MMFFAtomType::Br, BondType::Single)
        | (MMFFAtomType::Br, MMFFAtomType::C_3, BondType::Single) => Some(BondParams {
            k_bond: 3.0,
            r0: 1.94,
        }),
        (MMFFAtomType::C_3, MMFFAtomType::I, BondType::Single)
        | (MMFFAtomType::I, MMFFAtomType::C_3, BondType::Single) => Some(BondParams {
            k_bond: 1.706,
            r0: 2.090,
        }),
        (MMFFAtomType::C_AR, MMFFAtomType::I, BondType::Single)
        | (MMFFAtomType::I, MMFFAtomType::C_AR, BondType::Single) => Some(BondParams {
            k_bond: 1.781,
            r0: 2.075,
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
            k_bond: 3.031,
            r0: 1.891,
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
            k_bond: 5.002,
            r0: 1.246,
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
            k_bond: 2.790,
            r0: 1.830,
        }),
        (MMFFAtomType::P_3, MMFFAtomType::O_3, BondType::Single)
        | (MMFFAtomType::O_3, MMFFAtomType::P_3, BondType::Single) => Some(BondParams {
            k_bond: 4.0,
            r0: 1.60,
        }),
        (MMFFAtomType::P_4, MMFFAtomType::O_2, BondType::Double)
        | (MMFFAtomType::O_2, MMFFAtomType::P_4, BondType::Double) => Some(BondParams {
            k_bond: 9.020,
            r0: 1.496,
        }),
        // CS2 / S=C=S — S2CM (72) to C_1 (4) double (RDKit-extracted via symmetric scan)
        (MMFFAtomType::S2CM, MMFFAtomType::C_1, BondType::Double)
        | (MMFFAtomType::C_1, MMFFAtomType::S2CM, BondType::Double) => Some(BondParams {
            k_bond: 2.982,
            r0: 1.798,
        }),

        // 5-ring heteroaromatic bond params (RDKit verbose-extracted)
        (MMFFAtomType::C_2, MMFFAtomType::C5A, BondType::Aromatic)
        | (MMFFAtomType::C5A, MMFFAtomType::C_2, BondType::Aromatic)
        | (MMFFAtomType::C_2, MMFFAtomType::C5A, BondType::Single)
        | (MMFFAtomType::C5A, MMFFAtomType::C_2, BondType::Single) => Some(BondParams {
            k_bond: 5.468, r0: 1.423,
        }),
        (MMFFAtomType::C5A, MMFFAtomType::NPYL, BondType::Aromatic)
        | (MMFFAtomType::NPYL, MMFFAtomType::C5A, BondType::Aromatic)
        | (MMFFAtomType::C5A, MMFFAtomType::NPYL, BondType::Single)
        | (MMFFAtomType::NPYL, MMFFAtomType::C5A, BondType::Single) => Some(BondParams {
            k_bond: 6.301, r0: 1.364,
        }),
        (MMFFAtomType::C5A, MMFFAtomType::N5B, BondType::Aromatic)
        | (MMFFAtomType::N5B, MMFFAtomType::C5A, BondType::Aromatic)
        | (MMFFAtomType::C5A, MMFFAtomType::N5B, BondType::Single)
        | (MMFFAtomType::N5B, MMFFAtomType::C5A, BondType::Single) => Some(BondParams {
            k_bond: 8.326, r0: 1.313,
        }),
        (MMFFAtomType::C5B, MMFFAtomType::C5A, BondType::Aromatic)
        | (MMFFAtomType::C5A, MMFFAtomType::C5B, BondType::Aromatic)
        | (MMFFAtomType::C5B, MMFFAtomType::C5A, BondType::Single)
        | (MMFFAtomType::C5A, MMFFAtomType::C5B, BondType::Single) => Some(BondParams {
            k_bond: 7.118, r0: 1.377,
        }),
        (MMFFAtomType::N5B, MMFFAtomType::C5B, BondType::Aromatic)
        | (MMFFAtomType::C5B, MMFFAtomType::N5B, BondType::Aromatic)
        | (MMFFAtomType::N5B, MMFFAtomType::C5B, BondType::Single)
        | (MMFFAtomType::C5B, MMFFAtomType::N5B, BondType::Single) => Some(BondParams {
            k_bond: 4.456, r0: 1.369,
        }),
        (MMFFAtomType::N_AM, MMFFAtomType::C_2, BondType::Single)
        | (MMFFAtomType::C_2, MMFFAtomType::N_AM, BondType::Single) => Some(BondParams {
            k_bond: 5.829, r0: 1.369,
        }),
        (MMFFAtomType::C5B, MMFFAtomType::N_AM, BondType::Single)
        | (MMFFAtomType::N_AM, MMFFAtomType::C5B, BondType::Single) => Some(BondParams {
            k_bond: 5.952, r0: 1.376,
        }),
        (MMFFAtomType::NPYL, MMFFAtomType::H_N3, BondType::Single)
        | (MMFFAtomType::H_N3, MMFFAtomType::NPYL, BondType::Single) => Some(BondParams {
            k_bond: 7.112, r0: 1.012,
        }),
        // Purine 6-ring N_PL3-C bonds (RDKit verbose-estimated)
        (MMFFAtomType::N_PL3, MMFFAtomType::C_2, BondType::Single)
        | (MMFFAtomType::C_2, MMFFAtomType::N_PL3, BondType::Single) => Some(BondParams {
            k_bond: 6.110, r0: 1.370,
        }),
        (MMFFAtomType::N_PL3, MMFFAtomType::C_VIN, BondType::Single)
        | (MMFFAtomType::C_VIN, MMFFAtomType::N_PL3, BondType::Single) => Some(BondParams {
            k_bond: 6.110, r0: 1.370,
        }),
        // C_AR-C_1 single (aryl to nitrile C) — RDKit verbose-extracted
        (MMFFAtomType::C_AR, MMFFAtomType::C_1, BondType::Single)
        | (MMFFAtomType::C_1, MMFFAtomType::C_AR, BondType::Single) => Some(BondParams {
            k_bond: 5.445, r0: 1.424,
        }),
        (MMFFAtomType::C_3, MMFFAtomType::P_4, BondType::Single)
        | (MMFFAtomType::P_4, MMFFAtomType::C_3, BondType::Single) => Some(BondParams {
            k_bond: 2.980,
            r0: 1.810,
        }),
        (MMFFAtomType::P_4, MMFFAtomType::O_3, BondType::Single)
        | (MMFFAtomType::O_3, MMFFAtomType::P_4, BondType::Single) => Some(BondParams {
            k_bond: 5.243,
            r0: 1.630,
        }),
        (MMFFAtomType::P_4, MMFFAtomType::O_CO2, BondType::Double)
        | (MMFFAtomType::O_CO2, MMFFAtomType::P_4, BondType::Double) => Some(BondParams {
            k_bond: 8.296,
            r0: 1.510,
        }),

        _ => None,
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
    let c1 = 143.9325;
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
    let c1 = 143.9325;
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
