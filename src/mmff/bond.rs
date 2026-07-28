//! Bond stretching term for MMFF94

use super::MMFFAtomType;
use crate::molecule::BondType;

/// Bond stretching parameters
#[derive(Debug, Clone, Copy)]
pub struct BondParams {
    pub k_bond: f64, // mdyn/Å²
    pub r0: f64,     // Å
    pub cb: f64,     // cubic stretch constant (cs = -2*cb); 1.0 for most bonds
}

/// Per-bond-type cubic stretch constant (cb) overrides.
/// Derived from RDKit verbose energy inversion (bonds with |dr|>0.02).
/// cb=1.0 is the default; these override for specific type pairs where
/// RDKit's MMFFBOND table uses a different value.
const CB_OVERRIDES: &[(u8, u8, f64)] = &[
    // (type_min, type_max, cb) — populated from RDKit MMFFBOND table when available.
    // Derived-from-energy values were too noisy (biguanide regressed +0.43).
    // cb=1.0 is the default and gives RMSD=0.014; overrides need exact table values.
];

fn apply_cb_override(params: &mut BondParams, type1: MMFFAtomType, type2: MMFFAtomType) {
    use super::params::mmff_type_id;
    let ti = mmff_type_id(type1);
    let tk = mmff_type_id(type2);
    let (lo, hi) = if ti <= tk { (ti, tk) } else { (tk, ti) };
    for &(t_lo, t_hi, cb) in CB_OVERRIDES {
        if t_lo == lo && t_hi == hi {
            params.cb = cb;
            return;
        }
    }
}

/// Get bond parameters for atom types
pub fn get_bond_params(
    type1: MMFFAtomType,
    type2: MMFFAtomType,
    bond_type: BondType,
) -> Option<BondParams> {
    if let Some(mut p) = lookup_bond_params_exact(type1, type2, bond_type) {
        apply_cb_override(&mut p, type1, type2);
        return Some(p);
    }
    let base1 = super::base_type(type1);
    let base2 = super::base_type(type2);
    if base1 != type1 || base2 != type2 {
        if let Some(mut p) = lookup_bond_params_exact(base1, base2, bond_type) {
            apply_cb_override(&mut p, type1, type2);
            return Some(p);
        }
    }
    if let Some((kb, r0)) = super::estimation::estimate_bond_params(base1, base2, bond_type) {
        return Some(BondParams { k_bond: kb, r0, cb: 1.0 });
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
        cb: 1.0,
        }),
        (MMFFAtomType::C_2, MMFFAtomType::C_2, BondType::Single) => Some(BondParams {
            k_bond: 4.418,
            r0: 1.489,
        cb: 1.0,
        }),
        (MMFFAtomType::C_VIN, MMFFAtomType::C_VIN, BondType::Single) => Some(BondParams {
            k_bond: 5.310,
            r0: 1.430,
        cb: 1.0,
        }),
        (MMFFAtomType::C_VIN, MMFFAtomType::N_2, BondType::Single)
        | (MMFFAtomType::N_2, MMFFAtomType::C_VIN, BondType::Single) => Some(BondParams {
            k_bond: 6.385,
            r0: 1.360,
        cb: 1.0,
        }),
        (MMFFAtomType::C_2, MMFFAtomType::C_3, BondType::Single)
        | (MMFFAtomType::C_3, MMFFAtomType::C_2, BondType::Single) => Some(BondParams {
            k_bond: 4.19,
            r0: 1.492,
        cb: 1.0,
        }),
        (MMFFAtomType::C_2, MMFFAtomType::C_2, BondType::Double) => Some(BondParams {
            k_bond: 9.505,
            r0: 1.333,
        cb: 1.0,
        }),
        (MMFFAtomType::C_1, MMFFAtomType::C_1, BondType::Triple) => Some(BondParams {
            k_bond: 15.206,
            r0: 1.2,
        cb: 1.0,
        }),
        (MMFFAtomType::C_AR, MMFFAtomType::C_AR, BondType::Aromatic) => Some(BondParams {
            k_bond: 5.573,
            r0: 1.374,
        cb: 1.0,
        }),
        (MMFFAtomType::C_3, MMFFAtomType::C_AR, BondType::Single)
        | (MMFFAtomType::C_AR, MMFFAtomType::C_3, BondType::Single) => Some(BondParams {
            k_bond: 4.957,
            r0: 1.486,
        cb: 1.0,
        }),
        (MMFFAtomType::C_3, MMFFAtomType::C5B, BondType::Single)
        | (MMFFAtomType::C5B, MMFFAtomType::C_3, BondType::Single) => Some(BondParams {
            k_bond: 4.518,
            r0: 1.469,
        cb: 1.0,
        }),
        (MMFFAtomType::C_2, MMFFAtomType::C_AR, BondType::Single)
        | (MMFFAtomType::C_AR, MMFFAtomType::C_2, BondType::Single) => Some(BondParams {
            k_bond: 4.488,
            r0: 1.457,
        cb: 1.0,
        }),
        (MMFFAtomType::C_AR, MMFFAtomType::C_AR, BondType::Single) => Some(BondParams {
            k_bond: 5.0,
            r0: 1.484,
        cb: 1.0,
        }),

        // Vinyl (alkene, MMFF type 2) bonds — RDKit-extracted values
        (MMFFAtomType::C_VIN, MMFFAtomType::C_VIN, BondType::Double) => Some(BondParams {
            k_bond: 9.505,
            r0: 1.333,
        cb: 1.0,
        }),
        (MMFFAtomType::C_VIN, MMFFAtomType::C_3, BondType::Single)
        | (MMFFAtomType::C_3, MMFFAtomType::C_VIN, BondType::Single) => Some(BondParams {
            k_bond: 4.539,
            r0: 1.482,
        cb: 1.0,
        }),
        (MMFFAtomType::C_VIN, MMFFAtomType::C_2, BondType::Single)
        | (MMFFAtomType::C_2, MMFFAtomType::C_VIN, BondType::Single) => Some(BondParams {
            k_bond: 4.565,
            r0: 1.468,
        cb: 1.0,
        }),
        (MMFFAtomType::C_VIN, MMFFAtomType::C_AR, BondType::Single)
        | (MMFFAtomType::C_AR, MMFFAtomType::C_VIN, BondType::Single) => Some(BondParams {
            k_bond: 5.007,
            r0: 1.449,
        cb: 1.0,
        }),
        (MMFFAtomType::C_VIN, MMFFAtomType::H, BondType::Single)
        | (MMFFAtomType::H, MMFFAtomType::C_VIN, BondType::Single) => Some(BondParams {
            k_bond: 5.17,
            r0: 1.083,
        cb: 1.0,
        }),
        (MMFFAtomType::C_VIN, MMFFAtomType::O_R, BondType::Single)
        | (MMFFAtomType::O_R, MMFFAtomType::C_VIN, BondType::Single) => Some(BondParams {
            k_bond: 5.52,
            r0: 1.373,
        cb: 1.0,
        }),

        // Oxidized sulfur bonds — RDKit-extracted values
        (MMFFAtomType::S_OX, MMFFAtomType::C_3, BondType::Single)
        | (MMFFAtomType::C_3, MMFFAtomType::S_OX, BondType::Single) => Some(BondParams {
            k_bond: 2.841,
            r0: 1.813,
        cb: 1.0,
        }),
        (MMFFAtomType::S_OX, MMFFAtomType::C_AR, BondType::Single)
        | (MMFFAtomType::C_AR, MMFFAtomType::S_OX, BondType::Single) => Some(BondParams {
            k_bond: 3.098,
            r0: 1.787,
        cb: 1.0,
        }),
        (MMFFAtomType::S_OX, MMFFAtomType::O_2, BondType::Double)
        | (MMFFAtomType::O_2, MMFFAtomType::S_OX, BondType::Double) => Some(BondParams {
            k_bond: 8.77,
            r0: 1.5,
        cb: 1.0,
        }),
        (MMFFAtomType::S_O2, MMFFAtomType::C_3, BondType::Single)
        | (MMFFAtomType::C_3, MMFFAtomType::S_O2, BondType::Single) => Some(BondParams {
            k_bond: 3.258,
            r0: 1.772,
        cb: 1.0,
        }),
        (MMFFAtomType::S_O2, MMFFAtomType::C_AR, BondType::Single)
        | (MMFFAtomType::C_AR, MMFFAtomType::S_O2, BondType::Single) => Some(BondParams {
            k_bond: 3.281,
            r0: 1.77,
        cb: 1.0,
        }),
        (MMFFAtomType::S_O2, MMFFAtomType::O_CO2, BondType::Double)
        | (MMFFAtomType::O_CO2, MMFFAtomType::S_O2, BondType::Double) => Some(BondParams {
            k_bond: 10.748,
            r0: 1.45,
        cb: 1.0,
        }),
        // S-F (hypervalent S type 18) — from RDKit verbose (SF4)
        (MMFFAtomType::S_O2, MMFFAtomType::F, BondType::Single)
        | (MMFFAtomType::F, MMFFAtomType::S_O2, BondType::Single) => Some(BondParams {
            k_bond: 5.827,
            r0: 1.594,
        cb: 1.0,
        }),
        (MMFFAtomType::S_O2, MMFFAtomType::N_SO2, BondType::Single)
        | (MMFFAtomType::N_SO2, MMFFAtomType::S_O2, BondType::Single) => Some(BondParams {
            k_bond: 3.301,
            r0: 1.71,
        cb: 1.0,
        }),
        (MMFFAtomType::S_O2, MMFFAtomType::O_3, BondType::Single)
        | (MMFFAtomType::O_3, MMFFAtomType::S_O2, BondType::Single) => Some(BondParams {
            k_bond: 5.326,
            r0: 1.630,
        cb: 1.0,
        }),
        (MMFFAtomType::N_SO2, MMFFAtomType::C_3, BondType::Single)
        | (MMFFAtomType::C_3, MMFFAtomType::N_SO2, BondType::Single) => Some(BondParams {
            k_bond: 3.971,
            r0: 1.472,
        cb: 1.0,
        }),

        // Nitro group bonds — RDKit-extracted values
        // (N_NO2-O_CO2 params are identical for Single and Double in RDKit)
        (MMFFAtomType::N_NO2, MMFFAtomType::C_3, BondType::Single)
        | (MMFFAtomType::C_3, MMFFAtomType::N_NO2, BondType::Single) => Some(BondParams {
            k_bond: 3.844,
            r0: 1.48,
        cb: 1.0,
        }),
        (MMFFAtomType::N_NO2, MMFFAtomType::C_AR, BondType::Single)
        | (MMFFAtomType::C_AR, MMFFAtomType::N_NO2, BondType::Single) => Some(BondParams {
            k_bond: 4.705,
            r0: 1.431,
        cb: 1.0,
        }),
        (MMFFAtomType::N_NO2, MMFFAtomType::O_CO2, BondType::Double)
        | (MMFFAtomType::O_CO2, MMFFAtomType::N_NO2, BondType::Double)
        | (MMFFAtomType::N_NO2, MMFFAtomType::O_CO2, BondType::Single)
        | (MMFFAtomType::O_CO2, MMFFAtomType::N_NO2, BondType::Single) => Some(BondParams {
            k_bond: 9.42,
            r0: 1.233,
        cb: 1.0,
        }),
        // Nitroso N bonds — from RDKit verbose (nitrosomethane)
        (MMFFAtomType::N_NITROSO, MMFFAtomType::C_3, BondType::Single)
        | (MMFFAtomType::C_3, MMFFAtomType::N_NITROSO, BondType::Single) => Some(BondParams {
            k_bond: 3.813,
            r0: 1.482,
        cb: 1.0,
        }),
        (MMFFAtomType::N_NITROSO, MMFFAtomType::O_2, BondType::Double)
        | (MMFFAtomType::O_2, MMFFAtomType::N_NITROSO, BondType::Double) => Some(BondParams {
            k_bond: 9.329,
            r0: 1.235,
        cb: 1.0,
        }),

        // Carboxylate bonds — RDKit-extracted values (MMFF94s)
        // (1,41) C_3-CO2M, (32,41) O_CO2-CO2M
        (MMFFAtomType::C_3, MMFFAtomType::C_CO2, BondType::Single)
        | (MMFFAtomType::C_CO2, MMFFAtomType::C_3, BondType::Single) => Some(BondParams {
            k_bond: 3.830,
            r0: 1.510,
        cb: 1.0,
        }),
        (MMFFAtomType::C_AR, MMFFAtomType::C_CO2, BondType::Single)
        | (MMFFAtomType::C_CO2, MMFFAtomType::C_AR, BondType::Single) => Some(BondParams {
            k_bond: 4.537,
            r0: 1.468,
        cb: 1.0,
        }),
        (MMFFAtomType::O_CO2, MMFFAtomType::C_CO2, BondType::Double)
        | (MMFFAtomType::C_CO2, MMFFAtomType::O_CO2, BondType::Double)
        | (MMFFAtomType::O_CO2, MMFFAtomType::C_CO2, BondType::Single)
        | (MMFFAtomType::C_CO2, MMFFAtomType::O_CO2, BondType::Single) => Some(BondParams {
            k_bond: 9.756,
            r0: 1.261,
        cb: 1.0,
        }),

        // Pyridine N-oxide bonds — RDKit-extracted values (MMFF94s)
        // (37,69) C_AR-NPOX, (32,69) O_CO2-NPOX
        (MMFFAtomType::C_AR, MMFFAtomType::N_POX, BondType::Aromatic)
        | (MMFFAtomType::N_POX, MMFFAtomType::C_AR, BondType::Aromatic) => Some(BondParams {
            k_bond: 5.396,
            r0: 1.352,
        cb: 1.0,
        }),
        (MMFFAtomType::O_CO2, MMFFAtomType::N_POX, BondType::Single)
        | (MMFFAtomType::N_POX, MMFFAtomType::O_CO2, BondType::Single) => Some(BondParams {
            k_bond: 6.098,
            r0: 1.261,
        cb: 1.0,
        }),

        // Silicon bonds — RDKit-extracted values (MMFF94s)
        // (1,19) C_3-Si, (6,19) O_3-Si
        (MMFFAtomType::C_3, MMFFAtomType::Si, BondType::Single)
        | (MMFFAtomType::Si, MMFFAtomType::C_3, BondType::Single) => Some(BondParams {
            k_bond: 2.866,
            r0: 1.830,
        cb: 1.0,
        }),
        (MMFFAtomType::O_3, MMFFAtomType::Si, BondType::Single)
        | (MMFFAtomType::Si, MMFFAtomType::O_3, BondType::Single) => Some(BondParams {
            k_bond: 4.661,
            r0: 1.660,
        cb: 1.0,
        }),
        // Si-N bond (types 19,8) — from RDKit verbose
        (MMFFAtomType::Si, MMFFAtomType::N_3, BondType::Single)
        | (MMFFAtomType::N_3, MMFFAtomType::Si, BondType::Single) => Some(BondParams {
            k_bond: 4.254,
            r0: 1.700,
        cb: 1.0,
        }),

        // Imine N-H bond — RDKit-extracted value (MMFF94s)
        // (9,27) N_2-H_NIM
        (MMFFAtomType::N_2, MMFFAtomType::H_NIM, BondType::Single)
        | (MMFFAtomType::H_NIM, MMFFAtomType::N_2, BondType::Single) => Some(BondParams {
            k_bond: 6.230,
            r0: 1.026,
        cb: 1.0,
        }),

        // C-N bonds
        (MMFFAtomType::C_3, MMFFAtomType::N_3, BondType::Single)
        | (MMFFAtomType::N_3, MMFFAtomType::C_3, BondType::Single) => Some(BondParams {
            k_bond: 5.084,
            r0: 1.451,
        cb: 1.0,
        }),
        (MMFFAtomType::C_3, MMFFAtomType::N_2, BondType::Single)
        | (MMFFAtomType::N_2, MMFFAtomType::C_3, BondType::Single) => Some(BondParams {
            k_bond: 5.0,
            r0: 1.42,
        cb: 1.0,
        }),
        (MMFFAtomType::C_3, MMFFAtomType::N_AR, BondType::Single)
        | (MMFFAtomType::N_AR, MMFFAtomType::C_3, BondType::Single) => Some(BondParams {
            k_bond: 4.5,
            r0: 1.42,
        cb: 1.0,
        }),
        (MMFFAtomType::C_3, MMFFAtomType::NPYL, BondType::Single)
        | (MMFFAtomType::NPYL, MMFFAtomType::C_3, BondType::Single) => Some(BondParams {
            k_bond: 6.114,
            r0: 1.445,
        cb: 1.0,
        }),
        (MMFFAtomType::C_3, MMFFAtomType::N_PL3, BondType::Single)
        | (MMFFAtomType::N_PL3, MMFFAtomType::C_3, BondType::Single) => Some(BondParams {
            k_bond: 4.922,
            r0: 1.446,
        cb: 1.0,
        }),
        (MMFFAtomType::C_AR, MMFFAtomType::N_PL3, BondType::Single)
        | (MMFFAtomType::N_PL3, MMFFAtomType::C_AR, BondType::Single) => Some(BondParams {
            k_bond: 6.168,
            r0: 1.398,
        cb: 1.0,
        }),
        (MMFFAtomType::C_3, MMFFAtomType::N_AM, BondType::Single)
        | (MMFFAtomType::N_AM, MMFFAtomType::C_3, BondType::Single) => Some(BondParams {
            k_bond: 4.664,
            r0: 1.436,
        cb: 1.0,
        }),
        (MMFFAtomType::C_2, MMFFAtomType::N_2, BondType::Double)
        | (MMFFAtomType::N_2, MMFFAtomType::C_2, BondType::Double) => Some(BondParams {
            k_bond: 10.077,
            r0: 1.290,
        cb: 1.0,
        }),
        (MMFFAtomType::C_2, MMFFAtomType::C_1, BondType::Double)
        | (MMFFAtomType::C_1, MMFFAtomType::C_2, BondType::Double) => Some(BondParams {
            k_bond: 9.538,
            r0: 1.297,
        cb: 1.0,
        }),
        (MMFFAtomType::C_VIN, MMFFAtomType::C_1, BondType::Double)
        | (MMFFAtomType::C_1, MMFFAtomType::C_VIN, BondType::Double) => Some(BondParams {
            k_bond: 9.538,
            r0: 1.297,
        cb: 1.0,
        }),
        // Cumulated double bonds: C_1 (CSP, type 4) double bonds
        // C_1=O_2 (ketene C=C=O) — from RDKit verbose
        (MMFFAtomType::C_1, MMFFAtomType::O_2, BondType::Double)
        | (MMFFAtomType::O_2, MMFFAtomType::C_1, BondType::Double) => Some(BondParams {
            k_bond: 14.916,
            r0: 1.176,
        cb: 1.0,
        }),
        // C_1=N_2 (carbodiimide, isocyanate, isothiocyanate N=C=X) — from RDKit verbose
        (MMFFAtomType::C_1, MMFFAtomType::N_2, BondType::Double)
        | (MMFFAtomType::N_2, MMFFAtomType::C_1, BondType::Double) => Some(BondParams {
            k_bond: 15.589,
            r0: 1.172,
        cb: 1.0,
        }),
        // C_1=S_2 (isothiocyanate N=C=S) — from RDKit verbose
        (MMFFAtomType::C_1, MMFFAtomType::S_2, BondType::Double)
        | (MMFFAtomType::S_2, MMFFAtomType::C_1, BondType::Double) => Some(BondParams {
            k_bond: 2.982,
            r0: 1.798,
        cb: 1.0,
        }),
        (MMFFAtomType::C_2, MMFFAtomType::N_2, BondType::Single)
        | (MMFFAtomType::N_2, MMFFAtomType::C_2, BondType::Single) => Some(BondParams {
            k_bond: 6.273,
            r0: 1.364,
        cb: 1.0,
        }),
        (MMFFAtomType::C_VIN, MMFFAtomType::N_AM, BondType::Single)
        | (MMFFAtomType::N_AM, MMFFAtomType::C_VIN, BondType::Single) => Some(BondParams {
            k_bond: 6.329,
            r0: 1.362,
        cb: 1.0,
        }),
        (MMFFAtomType::C_2, MMFFAtomType::N_AR, BondType::Double)
        | (MMFFAtomType::N_AR, MMFFAtomType::C_2, BondType::Double) => Some(BondParams {
            k_bond: 7.0,
            r0: 1.28,
        cb: 1.0,
        }),
        (MMFFAtomType::C_AR, MMFFAtomType::N_AR, BondType::Aromatic)
        | (MMFFAtomType::N_AR, MMFFAtomType::C_AR, BondType::Aromatic) => Some(BondParams {
            k_bond: 5.737,
            r0: 1.333,
        cb: 1.0,
        }),
        (MMFFAtomType::N_AR, MMFFAtomType::C5A, BondType::Aromatic)
        | (MMFFAtomType::C5A, MMFFAtomType::N_AR, BondType::Aromatic) => Some(BondParams {
            k_bond: 7.299,
            r0: 1.330,
        cb: 1.0,
        }),
        (MMFFAtomType::C_AR, MMFFAtomType::C5A, BondType::Aromatic)
        | (MMFFAtomType::C5A, MMFFAtomType::C_AR, BondType::Aromatic) => Some(BondParams {
            k_bond: 6.095,
            r0: 1.372,
        cb: 1.0,
        }),
        (MMFFAtomType::C5B, MMFFAtomType::C_AR, BondType::Aromatic)
        | (MMFFAtomType::C_AR, MMFFAtomType::C5B, BondType::Aromatic) => Some(BondParams {
            k_bond: 6.161,
            r0: 1.379,
        cb: 1.0,
        }),
        (MMFFAtomType::C_1, MMFFAtomType::N_1, BondType::Triple)
        | (MMFFAtomType::N_1, MMFFAtomType::C_1, BondType::Triple) => Some(BondParams {
            k_bond: 16.582,
            r0: 1.160,
        cb: 1.0,
        }),
        // Isonitrile CID≡NID triple bond — from RDKit verbose
        (MMFFAtomType::CID, MMFFAtomType::NID, BondType::Triple)
        | (MMFFAtomType::NID, MMFFAtomType::CID, BondType::Triple) => Some(BondParams {
            k_bond: 15.749,
            r0: 1.170,
        cb: 1.0,
        }),
        // C_3-NID (methyl isocyanide) — from RDKit verbose
        (MMFFAtomType::C_3, MMFFAtomType::NID, BondType::Single)
        | (MMFFAtomType::NID, MMFFAtomType::C_3, BondType::Single) => Some(BondParams {
            k_bond: 4.845,
            r0: 1.424,
        cb: 1.0,
        }),
        // Cumulated N=N bonds (diazomethane/azide) — from RDKit verbose
        (MMFFAtomType::C_2, MMFFAtomType::N_2Z, BondType::Double)
        | (MMFFAtomType::N_2Z, MMFFAtomType::C_2, BondType::Double) => Some(BondParams {
            k_bond: 7.637,
            r0: 1.320,
        cb: 1.0,
        }),
        (MMFFAtomType::N_2Z, MMFFAtomType::N_1M, BondType::Double)
        | (MMFFAtomType::N_1M, MMFFAtomType::N_2Z, BondType::Double) => Some(BondParams {
            k_bond: 12.192,
            r0: 1.140,
        cb: 1.0,
        }),
        // N_2=N_2Z (azide central bond) — from RDKit verbose
        (MMFFAtomType::N_2, MMFFAtomType::N_2Z, BondType::Double)
        | (MMFFAtomType::N_2Z, MMFFAtomType::N_2, BondType::Double) => Some(BondParams {
            k_bond: 7.291,
            r0: 1.242,
        cb: 1.0,
        }),
        (MMFFAtomType::C_3, MMFFAtomType::C_1, BondType::Single)
        | (MMFFAtomType::C_1, MMFFAtomType::C_3, BondType::Single) => Some(BondParams {
            k_bond: 4.707,
            r0: 1.459,
        cb: 1.0,
        }),
        // 3-ring (cyclopropane/epoxide/aziridine) specific bond params (CR3R=22)
        (MMFFAtomType::CR3R, MMFFAtomType::CR3R, BondType::Single) => Some(BondParams {
            k_bond: 3.969,
            r0: 1.499,
        cb: 1.0,
        }),
        (MMFFAtomType::CR3R, MMFFAtomType::O_3, BondType::Single)
        | (MMFFAtomType::O_3, MMFFAtomType::CR3R, BondType::Single) => Some(BondParams {
            k_bond: 4.556,
            r0: 1.433,
        cb: 1.0,
        }),
        (MMFFAtomType::CR3R, MMFFAtomType::O_R, BondType::Single)
        | (MMFFAtomType::O_R, MMFFAtomType::CR3R, BondType::Single) => Some(BondParams {
            k_bond: 4.556,
            r0: 1.433,
        cb: 1.0,
        }),
        (MMFFAtomType::CR3R, MMFFAtomType::H, BondType::Single)
        | (MMFFAtomType::H, MMFFAtomType::CR3R, BondType::Single) => Some(BondParams {
            k_bond: 5.191,
            r0: 1.082,
        cb: 1.0,
        }),
        (MMFFAtomType::C_2, MMFFAtomType::CR3R, BondType::Single)
        | (MMFFAtomType::CR3R, MMFFAtomType::C_2, BondType::Single) => Some(BondParams {
            k_bond: 4.926,
            r0: 1.448,
        cb: 1.0,
        }),
        (MMFFAtomType::C_VIN, MMFFAtomType::CR3R, BondType::Single)
        | (MMFFAtomType::CR3R, MMFFAtomType::C_VIN, BondType::Single) => Some(BondParams {
            k_bond: 4.926,
            r0: 1.448,
        cb: 1.0,
        }),
        // 4-ring (cyclobutane) specific bond params (CR4R=20)
        (MMFFAtomType::CR4R, MMFFAtomType::CR4R, BondType::Single) => Some(BondParams {
            k_bond: 3.663,
            r0: 1.526,
        cb: 1.0,
        }),
        (MMFFAtomType::CR4R, MMFFAtomType::H, BondType::Single)
        | (MMFFAtomType::H, MMFFAtomType::CR4R, BondType::Single) => Some(BondParams {
            k_bond: 4.852,
            r0: 1.093,
        cb: 1.0,
        }),
        // CR4R-S (4-ring S) — from RDKit verbose (thietane)
        (MMFFAtomType::CR4R, MMFFAtomType::S_3, BondType::Single)
        | (MMFFAtomType::S_3, MMFFAtomType::CR4R, BondType::Single) => Some(BondParams {
            k_bond: 2.757,
            r0: 1.822,
        cb: 1.0,
        }),
        // CR4R-O_R (4-ring O) — from RDKit verbose (oxetane)
        (MMFFAtomType::CR4R, MMFFAtomType::O_R, BondType::Single)
        | (MMFFAtomType::O_R, MMFFAtomType::CR4R, BondType::Single) => Some(BondParams {
            k_bond: 5.623,
            r0: 1.433,
        cb: 1.0,
        }),
        (MMFFAtomType::CR3R, MMFFAtomType::N_3, BondType::Single)
        | (MMFFAtomType::N_3, MMFFAtomType::CR3R, BondType::Single) => Some(BondParams {
            k_bond: 4.223,
            r0: 1.457,
        cb: 1.0,
        }),
        (MMFFAtomType::C_AR, MMFFAtomType::N_AM, BondType::Aromatic)
        | (MMFFAtomType::N_AM, MMFFAtomType::C_AR, BondType::Aromatic) => Some(BondParams {
            k_bond: 5.5,
            r0: 1.37,
        cb: 1.0,
        }),
        (MMFFAtomType::C_AR, MMFFAtomType::N_AM, BondType::Single)
        | (MMFFAtomType::N_AM, MMFFAtomType::C_AR, BondType::Single) => Some(BondParams {
            k_bond: 5.482,
            r0: 1.395,
        cb: 1.0,
        }),

        // C-O bonds
        (MMFFAtomType::C_3, MMFFAtomType::O_3, BondType::Single)
        | (MMFFAtomType::O_3, MMFFAtomType::C_3, BondType::Single) => Some(BondParams {
            k_bond: 5.047,
            r0: 1.418,
        cb: 1.0,
        }),
        (MMFFAtomType::C_3, MMFFAtomType::O_2, BondType::Single)
        | (MMFFAtomType::O_2, MMFFAtomType::C_3, BondType::Single) => Some(BondParams {
            k_bond: 5.5,
            r0: 1.40,
        cb: 1.0,
        }),
        (MMFFAtomType::C_3, MMFFAtomType::O_R, BondType::Single)
        | (MMFFAtomType::O_R, MMFFAtomType::C_3, BondType::Single) => Some(BondParams {
            k_bond: 5.047,
            r0: 1.418,
        cb: 1.0,
        }),
        (MMFFAtomType::C_3, MMFFAtomType::O_CO2, BondType::Single)
        | (MMFFAtomType::O_CO2, MMFFAtomType::C_3, BondType::Single) => Some(BondParams {
            k_bond: 5.5,
            r0: 1.40,
        cb: 1.0,
        }),
        (MMFFAtomType::C_2, MMFFAtomType::O_2, BondType::Double)
        | (MMFFAtomType::O_2, MMFFAtomType::C_2, BondType::Double) => Some(BondParams {
            k_bond: 12.95,
            r0: 1.222,
        cb: 1.0,
        }),
        (MMFFAtomType::C_2, MMFFAtomType::O_R, BondType::Double)
        | (MMFFAtomType::O_R, MMFFAtomType::C_2, BondType::Double) => Some(BondParams {
            k_bond: 10.5,
            r0: 1.23,
        cb: 1.0,
        }),
        (MMFFAtomType::C_AR, MMFFAtomType::O_CO2, BondType::Double)
        | (MMFFAtomType::O_CO2, MMFFAtomType::C_AR, BondType::Double) => Some(BondParams {
            k_bond: 10.0,
            r0: 1.23,
        cb: 1.0,
        }),
        (MMFFAtomType::C_2, MMFFAtomType::O_CO2, BondType::Double)
        | (MMFFAtomType::O_CO2, MMFFAtomType::C_2, BondType::Double) => Some(BondParams {
            k_bond: 12.95,
            r0: 1.222,
        cb: 1.0,
        }),
        (MMFFAtomType::C_2, MMFFAtomType::O_3, BondType::Single)
        | (MMFFAtomType::O_3, MMFFAtomType::C_2, BondType::Single) => Some(BondParams {
            k_bond: 5.801,
            r0: 1.355,
        cb: 1.0,
        }),
        (MMFFAtomType::C_2, MMFFAtomType::O_R, BondType::Single)
        | (MMFFAtomType::O_R, MMFFAtomType::C_2, BondType::Single) => Some(BondParams {
            k_bond: 5.801,
            r0: 1.355,
        cb: 1.0,
        }),
        (MMFFAtomType::C_AR, MMFFAtomType::O_R, BondType::Aromatic)
        | (MMFFAtomType::O_R, MMFFAtomType::C_AR, BondType::Aromatic) => Some(BondParams {
            k_bond: 5.0,
            r0: 1.37,
        cb: 1.0,
        }),
        (MMFFAtomType::C_AR, MMFFAtomType::O_3, BondType::Single)
        | (MMFFAtomType::O_3, MMFFAtomType::C_AR, BondType::Single) => Some(BondParams {
            k_bond: 5.614,
            r0: 1.376,
        cb: 1.0,
        }),

        // C-S bonds
        (MMFFAtomType::C_3, MMFFAtomType::S_3, BondType::Single)
        | (MMFFAtomType::S_3, MMFFAtomType::C_3, BondType::Single) => Some(BondParams {
            k_bond: 2.893,
            r0: 1.805,
        cb: 1.0,
        }),
        (MMFFAtomType::S_3, MMFFAtomType::S_3, BondType::Single) => Some(BondParams {
            k_bond: 2.531,
            r0: 2.050,
        cb: 1.0,
        }),
        (MMFFAtomType::C_2, MMFFAtomType::S_2, BondType::Double)
        | (MMFFAtomType::S_2, MMFFAtomType::C_2, BondType::Double) => Some(BondParams {
            k_bond: 4.735,
            r0: 1.665,
        cb: 1.0,
        }),
        (MMFFAtomType::C_AR, MMFFAtomType::S_AR, BondType::Aromatic)
        | (MMFFAtomType::S_AR, MMFFAtomType::C_AR, BondType::Aromatic) => Some(BondParams {
            k_bond: 4.0,
            r0: 1.71,
        cb: 1.0,
        }),
        // 5-ring heteroaromatic specific bond params (C5A=63, C5B=64) from RDKit verbose
        (MMFFAtomType::C5A, MMFFAtomType::S_AR, BondType::Aromatic)
        | (MMFFAtomType::S_AR, MMFFAtomType::C5A, BondType::Aromatic) => Some(BondParams {
            k_bond: 3.589,
            r0: 1.717,
        cb: 1.0,
        }),
        (MMFFAtomType::C5B, MMFFAtomType::C5B, BondType::Aromatic) => Some(BondParams {
            k_bond: 4.313,
            r0: 1.418,
        cb: 1.0,
        }),
        (MMFFAtomType::NPYL, MMFFAtomType::N5A, BondType::Aromatic)
        | (MMFFAtomType::N5A, MMFFAtomType::NPYL, BondType::Aromatic) => Some(BondParams {
            k_bond: 5.513,
            r0: 1.339,
        cb: 1.0,
        }),
        (MMFFAtomType::N5A, MMFFAtomType::C5B, BondType::Aromatic)
        | (MMFFAtomType::C5B, MMFFAtomType::N5A, BondType::Aromatic) => Some(BondParams {
            k_bond: 8.258,
            r0: 1.335,
        cb: 1.0,
        }),
        (MMFFAtomType::C5A, MMFFAtomType::OFUR, BondType::Aromatic)
        | (MMFFAtomType::OFUR, MMFFAtomType::C5A, BondType::Aromatic) => Some(BondParams {
            k_bond: 5.787,
            r0: 1.360,
        cb: 1.0,
        }),
        // N5A-OFUR (oxazole N-O) — from RDKit verbose
        (MMFFAtomType::N5A, MMFFAtomType::OFUR, BondType::Single)
        | (MMFFAtomType::N5A, MMFFAtomType::OFUR, BondType::Aromatic)
        | (MMFFAtomType::OFUR, MMFFAtomType::N5A, BondType::Single)
        | (MMFFAtomType::OFUR, MMFFAtomType::N5A, BondType::Aromatic) => Some(BondParams {
            k_bond: 4.756,
            r0: 1.388,
        cb: 1.0,
        }),
        (MMFFAtomType::C5A, MMFFAtomType::H, BondType::Single)
        | (MMFFAtomType::H, MMFFAtomType::C5A, BondType::Single) => Some(BondParams {
            k_bond: 5.531,
            r0: 1.080,
        cb: 1.0,
        }),
        (MMFFAtomType::C5B, MMFFAtomType::H, BondType::Single)
        | (MMFFAtomType::H, MMFFAtomType::C5B, BondType::Single) => Some(BondParams {
            k_bond: 5.506,
            r0: 1.080,
        cb: 1.0,
        }),

        // C-H bonds (symmetric)
        (MMFFAtomType::H, MMFFAtomType::C_3, BondType::Single)
        | (MMFFAtomType::C_3, MMFFAtomType::H, BondType::Single) => Some(BondParams {
            k_bond: 4.766,
            r0: 1.093,
        cb: 1.0,
        }),
        (MMFFAtomType::H, MMFFAtomType::C_2, BondType::Single)
        | (MMFFAtomType::C_2, MMFFAtomType::H, BondType::Single) => Some(BondParams {
            k_bond: 4.65,
            r0: 1.101,
        cb: 1.0,
        }),
        (MMFFAtomType::H, MMFFAtomType::C_1, BondType::Single)
        | (MMFFAtomType::C_1, MMFFAtomType::H, BondType::Single) => Some(BondParams {
            k_bond: 5.726,
            r0: 1.065,
        cb: 1.0,
        }),
        (MMFFAtomType::H, MMFFAtomType::C_AR, BondType::Single)
        | (MMFFAtomType::C_AR, MMFFAtomType::H, BondType::Single) => Some(BondParams {
            k_bond: 5.306,
            r0: 1.084,
        cb: 1.0,
        }),
        (MMFFAtomType::H, MMFFAtomType::C_CAT, BondType::Single)
        | (MMFFAtomType::C_CAT, MMFFAtomType::H, BondType::Single) => Some(BondParams {
            k_bond: 4.766,
            r0: 1.093,
        cb: 1.0,
        }),
        (MMFFAtomType::H, MMFFAtomType::C_AN, BondType::Single)
        | (MMFFAtomType::C_AN, MMFFAtomType::H, BondType::Single) => Some(BondParams {
            k_bond: 4.766,
            r0: 1.093,
        cb: 1.0,
        }),
        // Guanidinium CGD+-NCN+ bond — from RDKit verbose (same for Single/Double)
        (MMFFAtomType::C_AN, MMFFAtomType::NCN_PLUS, BondType::Single)
        | (MMFFAtomType::NCN_PLUS, MMFFAtomType::C_AN, BondType::Single)
        | (MMFFAtomType::C_AN, MMFFAtomType::NCN_PLUS, BondType::Double)
        | (MMFFAtomType::NCN_PLUS, MMFFAtomType::C_AN, BondType::Double) => Some(BondParams {
            k_bond: 7.227,
            r0: 1.319,
        cb: 1.0,
        }),
        // NCN+-HNRP bond — from RDKit verbose
        (MMFFAtomType::NCN_PLUS, MMFFAtomType::HNRP, BondType::Single)
        | (MMFFAtomType::HNRP, MMFFAtomType::NCN_PLUS, BondType::Single) => Some(BondParams {
            k_bond: 6.744,
            r0: 1.014,
        cb: 1.0,
        }),

        // N-H bonds (symmetric)
        (MMFFAtomType::H, MMFFAtomType::N_3, BondType::Single)
        | (MMFFAtomType::N_3, MMFFAtomType::H, BondType::Single) => Some(BondParams {
            k_bond: 6.490,
            r0: 1.019,
        cb: 1.0,
        }),
        (MMFFAtomType::H, MMFFAtomType::N_2, BondType::Single)
        | (MMFFAtomType::N_2, MMFFAtomType::H, BondType::Single) => Some(BondParams {
            k_bond: 5.5,
            r0: 1.000,
        cb: 1.0,
        }),
        (MMFFAtomType::H, MMFFAtomType::N_AR, BondType::Single)
        | (MMFFAtomType::N_AR, MMFFAtomType::H, BondType::Single) => Some(BondParams {
            k_bond: 5.0,
            r0: 1.010,
        cb: 1.0,
        }),
        (MMFFAtomType::H, MMFFAtomType::N_PL3, BondType::Single)
        | (MMFFAtomType::N_PL3, MMFFAtomType::H, BondType::Single) => Some(BondParams {
            k_bond: 6.576,
            r0: 1.018,
        cb: 1.0,
        }),
        (MMFFAtomType::H, MMFFAtomType::N_AM, BondType::Single)
        | (MMFFAtomType::N_AM, MMFFAtomType::H, BondType::Single) => Some(BondParams {
            k_bond: 6.663,
            r0: 1.015,
        cb: 1.0,
        }),
        (MMFFAtomType::N_SO2, MMFFAtomType::H, BondType::Single)
        | (MMFFAtomType::H, MMFFAtomType::N_SO2, BondType::Single) => Some(BondParams {
            k_bond: 6.265,
            r0: 1.028,
        cb: 1.0,
        }),
        (MMFFAtomType::H, MMFFAtomType::N_4, BondType::Single)
        | (MMFFAtomType::N_4, MMFFAtomType::H, BondType::Single) => Some(BondParams {
            k_bond: 6.163,
            r0: 1.028,
        cb: 1.0,
        }),
        (MMFFAtomType::C_3, MMFFAtomType::N_4, BondType::Single)
        | (MMFFAtomType::N_4, MMFFAtomType::C_3, BondType::Single) => Some(BondParams {
            k_bond: 3.844,
            r0: 1.480,
        cb: 1.0,
        }),

        // O-H bonds (symmetric)
        (MMFFAtomType::H_OH, MMFFAtomType::OH2, BondType::Single)
        | (MMFFAtomType::OH2, MMFFAtomType::H_OH, BondType::Single) => Some(BondParams {
            k_bond: 7.88,
            r0: 0.969,
        cb: 1.0,
        }),
        (MMFFAtomType::H_ONC, MMFFAtomType::O_3, BondType::Single)
        | (MMFFAtomType::O_3, MMFFAtomType::H_ONC, BondType::Single) => Some(BondParams {
            k_bond: 7.794,
            r0: 0.972,
        cb: 1.0,
        }),
        (MMFFAtomType::H_OH, MMFFAtomType::O_3, BondType::Single)
        | (MMFFAtomType::O_3, MMFFAtomType::H_OH, BondType::Single) => Some(BondParams {
            k_bond: 7.88,
            r0: 0.969,
        cb: 1.0,
        }),
        (MMFFAtomType::H_COOH, MMFFAtomType::O_3, BondType::Single)
        | (MMFFAtomType::O_3, MMFFAtomType::H_COOH, BondType::Single) => Some(BondParams {
            k_bond: 7.403,
            r0: 0.981,
        cb: 1.0,
        }),
        (MMFFAtomType::H_COOH, MMFFAtomType::O_R, BondType::Single)
        | (MMFFAtomType::O_R, MMFFAtomType::H_COOH, BondType::Single) => Some(BondParams {
            k_bond: 7.403,
            r0: 0.981,
        cb: 1.0,
        }),
        (MMFFAtomType::H_OAR, MMFFAtomType::O_3, BondType::Single)
        | (MMFFAtomType::O_3, MMFFAtomType::H_OAR, BondType::Single) => Some(BondParams {
            k_bond: 7.839,
            r0: 0.973,
        cb: 1.0,
        }),
        (MMFFAtomType::H_OAR, MMFFAtomType::O_R, BondType::Single)
        | (MMFFAtomType::O_R, MMFFAtomType::H_OAR, BondType::Single) => Some(BondParams {
            k_bond: 7.839,
            r0: 0.973,
        cb: 1.0,
        }),
        (MMFFAtomType::HOS, MMFFAtomType::O_3, BondType::Single)
        | (MMFFAtomType::O_3, MMFFAtomType::HOS, BondType::Single) => Some(BondParams {
            k_bond: 7.143,
            r0: 0.986,
        cb: 1.0,
        }),
        (MMFFAtomType::HOS, MMFFAtomType::O_R, BondType::Single)
        | (MMFFAtomType::O_R, MMFFAtomType::HOS, BondType::Single) => Some(BondParams {
            k_bond: 7.143,
            r0: 0.986,
        cb: 1.0,
        }),
        (MMFFAtomType::H, MMFFAtomType::O_3, BondType::Single)
        | (MMFFAtomType::O_3, MMFFAtomType::H, BondType::Single) => Some(BondParams {
            k_bond: 7.794,
            r0: 0.972,
        cb: 1.0,
        }),
        (MMFFAtomType::H, MMFFAtomType::O_2, BondType::Single)
        | (MMFFAtomType::O_2, MMFFAtomType::H, BondType::Single) => Some(BondParams {
            k_bond: 7.794,
            r0: 0.972,
        cb: 1.0,
        }),
        (MMFFAtomType::H, MMFFAtomType::O_R, BondType::Single)
        | (MMFFAtomType::O_R, MMFFAtomType::H, BondType::Single) => Some(BondParams {
            k_bond: 7.794,
            r0: 0.972,
        cb: 1.0,
        }),

        // S-H bonds (symmetric)
        (MMFFAtomType::H, MMFFAtomType::S_3, BondType::Single)
        | (MMFFAtomType::S_3, MMFFAtomType::H, BondType::Single) => Some(BondParams {
            k_bond: 4.014,
            r0: 1.341,
        cb: 1.0,
        }),
        (MMFFAtomType::H, MMFFAtomType::S_2, BondType::Single)
        | (MMFFAtomType::S_2, MMFFAtomType::H, BondType::Single) => Some(BondParams {
            k_bond: 4.0,
            r0: 1.336,
        cb: 1.0,
        }),

        // Halogen bonds (symmetric)
        (MMFFAtomType::C_3, MMFFAtomType::F, BondType::Single)
        | (MMFFAtomType::F, MMFFAtomType::C_3, BondType::Single) => Some(BondParams {
            k_bond: 6.011,
            r0: 1.360,
        cb: 1.0,
        }),
        (MMFFAtomType::C_3, MMFFAtomType::Cl, BondType::Single)
        | (MMFFAtomType::Cl, MMFFAtomType::C_3, BondType::Single) => Some(BondParams {
            k_bond: 2.974,
            r0: 1.773,
        cb: 1.0,
        }),
        (MMFFAtomType::C_VIN, MMFFAtomType::Cl, BondType::Single)
        | (MMFFAtomType::Cl, MMFFAtomType::C_VIN, BondType::Single) => Some(BondParams {
            k_bond: 3.390,
            r0: 1.720,
        cb: 1.0,
        }),
        (MMFFAtomType::C_2, MMFFAtomType::Cl, BondType::Single)
        | (MMFFAtomType::Cl, MMFFAtomType::C_2, BondType::Single) => Some(BondParams {
            k_bond: 3.449,
            r0: 1.715,
        cb: 1.0,
        }),
        (MMFFAtomType::C_3, MMFFAtomType::Br, BondType::Single)
        | (MMFFAtomType::Br, MMFFAtomType::C_3, BondType::Single) => Some(BondParams {
            k_bond: 2.529,
            r0: 1.949,
        cb: 1.0,
        }),
        (MMFFAtomType::C_3, MMFFAtomType::I, BondType::Single)
        | (MMFFAtomType::I, MMFFAtomType::C_3, BondType::Single) => Some(BondParams {
            k_bond: 1.706,
            r0: 2.090,
        cb: 1.0,
        }),
        (MMFFAtomType::C_AR, MMFFAtomType::I, BondType::Single)
        | (MMFFAtomType::I, MMFFAtomType::C_AR, BondType::Single) => Some(BondParams {
            k_bond: 1.781,
            r0: 2.075,
        cb: 1.0,
        }),
        (MMFFAtomType::C_AR, MMFFAtomType::F, BondType::Single)
        | (MMFFAtomType::F, MMFFAtomType::C_AR, BondType::Single) => Some(BondParams {
            k_bond: 5.5,
            r0: 1.33,
        cb: 1.0,
        }),
        (MMFFAtomType::C_AR, MMFFAtomType::Cl, BondType::Single)
        | (MMFFAtomType::Cl, MMFFAtomType::C_AR, BondType::Single) => Some(BondParams {
            k_bond: 3.5,
            r0: 1.72,
        cb: 1.0,
        }),
        (MMFFAtomType::C_AR, MMFFAtomType::Br, BondType::Single)
        | (MMFFAtomType::Br, MMFFAtomType::C_AR, BondType::Single) => Some(BondParams {
            k_bond: 3.031,
            r0: 1.891,
        cb: 1.0,
        }),
        (MMFFAtomType::C_2, MMFFAtomType::F, BondType::Single)
        | (MMFFAtomType::F, MMFFAtomType::C_2, BondType::Single) => Some(BondParams {
            k_bond: 6.0,
            r0: 1.30,
        cb: 1.0,
        }),
        // N-N bonds
        (MMFFAtomType::N_3, MMFFAtomType::N_3, BondType::Single) => Some(BondParams {
            k_bond: 3.5,
            r0: 1.45,
        cb: 1.0,
        }),
        (MMFFAtomType::N_2, MMFFAtomType::N_2, BondType::Double) => Some(BondParams {
            k_bond: 6.0,
            r0: 1.25,
        cb: 1.0,
        }),
        (MMFFAtomType::N_3, MMFFAtomType::N_AR, BondType::Single)
        | (MMFFAtomType::N_AR, MMFFAtomType::N_3, BondType::Single) => Some(BondParams {
            k_bond: 4.0,
            r0: 1.40,
        cb: 1.0,
        }),
        (MMFFAtomType::N_AR, MMFFAtomType::N_AR, BondType::Aromatic) => Some(BondParams {
            k_bond: 5.002,
            r0: 1.246,
        cb: 1.0,
        }),
        (MMFFAtomType::N_3, MMFFAtomType::C_2, BondType::Single)
        | (MMFFAtomType::C_2, MMFFAtomType::N_3, BondType::Single) => Some(BondParams {
            k_bond: 5.0,
            r0: 1.42,
        cb: 1.0,
        }),

        // O-O bonds
        (MMFFAtomType::O_3, MMFFAtomType::O_3, BondType::Single) => Some(BondParams {
            k_bond: 4.088,
            r0: 1.449,
        cb: 1.0,
        }),
        (MMFFAtomType::O_3, MMFFAtomType::O_2, BondType::Single)
        | (MMFFAtomType::O_2, MMFFAtomType::O_3, BondType::Single) => Some(BondParams {
            k_bond: 4.5,
            r0: 1.45,
        cb: 1.0,
        }),

        // P bonds
        (MMFFAtomType::C_3, MMFFAtomType::P_3, BondType::Single)
        | (MMFFAtomType::P_3, MMFFAtomType::C_3, BondType::Single) => Some(BondParams {
            k_bond: 2.790,
            r0: 1.830,
        cb: 1.0,
        }),
        (MMFFAtomType::P_3, MMFFAtomType::O_3, BondType::Single)
        | (MMFFAtomType::O_3, MMFFAtomType::P_3, BondType::Single) => Some(BondParams {
            k_bond: 4.0,
            r0: 1.60,
        cb: 1.0,
        }),
        // P-H bond (types 26,71) — from RDKit verbose
        (MMFFAtomType::P_3, MMFFAtomType::HS, BondType::Single)
        | (MMFFAtomType::HS, MMFFAtomType::P_3, BondType::Single) => Some(BondParams {
            k_bond: 2.959,
            r0: 1.415,
        cb: 1.0,
        }),
        (MMFFAtomType::P_4, MMFFAtomType::O_2, BondType::Double)
        | (MMFFAtomType::O_2, MMFFAtomType::P_4, BondType::Double) => Some(BondParams {
            k_bond: 9.020,
            r0: 1.496,
        cb: 1.0,
        }),
        // CS2 / S=C=S — S2CM (72) to C_1 (4) double
        // Params from empirical rule: fitted to RDKit binary (verbose truncates to 3 dp)
        (MMFFAtomType::S2CM, MMFFAtomType::C_1, BondType::Double)
        | (MMFFAtomType::C_1, MMFFAtomType::S2CM, BondType::Double) => Some(BondParams {
            k_bond: 2.982333,
            r0: 1.798345,
        cb: 1.0,
        }),

        // 5-ring heteroaromatic bond params (RDKit verbose-extracted)
        (MMFFAtomType::C_2, MMFFAtomType::C5A, BondType::Aromatic)
        | (MMFFAtomType::C5A, MMFFAtomType::C_2, BondType::Aromatic)
        | (MMFFAtomType::C_2, MMFFAtomType::C5A, BondType::Single)
        | (MMFFAtomType::C5A, MMFFAtomType::C_2, BondType::Single) => Some(BondParams {
            k_bond: 5.468, r0: 1.423,
        cb: 1.0,
        }),
        (MMFFAtomType::C5A, MMFFAtomType::NPYL, BondType::Aromatic)
        | (MMFFAtomType::NPYL, MMFFAtomType::C5A, BondType::Aromatic)
        | (MMFFAtomType::C5A, MMFFAtomType::NPYL, BondType::Single)
        | (MMFFAtomType::NPYL, MMFFAtomType::C5A, BondType::Single) => Some(BondParams {
            k_bond: 6.301, r0: 1.364,
        cb: 1.0,
        }),
        (MMFFAtomType::C5A, MMFFAtomType::N5B, BondType::Aromatic)
        | (MMFFAtomType::N5B, MMFFAtomType::C5A, BondType::Aromatic)
        | (MMFFAtomType::C5A, MMFFAtomType::N5B, BondType::Single)
        | (MMFFAtomType::N5B, MMFFAtomType::C5A, BondType::Single) => Some(BondParams {
            k_bond: 8.326, r0: 1.313,
        cb: 1.0,
        }),
        (MMFFAtomType::C5B, MMFFAtomType::C5A, BondType::Aromatic)
        | (MMFFAtomType::C5A, MMFFAtomType::C5B, BondType::Aromatic)
        | (MMFFAtomType::C5B, MMFFAtomType::C5A, BondType::Single)
        | (MMFFAtomType::C5A, MMFFAtomType::C5B, BondType::Single) => Some(BondParams {
            k_bond: 7.118, r0: 1.377,
        cb: 1.0,
        }),
        (MMFFAtomType::N5B, MMFFAtomType::C5B, BondType::Aromatic)
        | (MMFFAtomType::C5B, MMFFAtomType::N5B, BondType::Aromatic)
        | (MMFFAtomType::N5B, MMFFAtomType::C5B, BondType::Single)
        | (MMFFAtomType::C5B, MMFFAtomType::N5B, BondType::Single) => Some(BondParams {
            k_bond: 4.456, r0: 1.369,
        cb: 1.0,
        }),
        (MMFFAtomType::N_AM, MMFFAtomType::C_2, BondType::Single)
        | (MMFFAtomType::C_2, MMFFAtomType::N_AM, BondType::Single) => Some(BondParams {
            k_bond: 5.829, r0: 1.369,
        cb: 1.0,
        }),
        (MMFFAtomType::C5B, MMFFAtomType::N_AM, BondType::Single)
        | (MMFFAtomType::N_AM, MMFFAtomType::C5B, BondType::Single) => Some(BondParams {
            k_bond: 5.952, r0: 1.376,
        cb: 1.0,
        }),
        (MMFFAtomType::NPYL, MMFFAtomType::H_N3, BondType::Single)
        | (MMFFAtomType::H_N3, MMFFAtomType::NPYL, BondType::Single) => Some(BondParams {
            k_bond: 7.112, r0: 1.012,
        cb: 1.0,
        }),
        // Purine 6-ring N_PL3-C bonds (RDKit verbose-estimated)
        (MMFFAtomType::N_PL3, MMFFAtomType::C_2, BondType::Single)
        | (MMFFAtomType::C_2, MMFFAtomType::N_PL3, BondType::Single) => Some(BondParams {
            k_bond: 6.110, r0: 1.370,
        cb: 1.0,
        }),
        (MMFFAtomType::N_PL3, MMFFAtomType::C_VIN, BondType::Single)
        | (MMFFAtomType::C_VIN, MMFFAtomType::N_PL3, BondType::Single) => Some(BondParams {
            k_bond: 6.110, r0: 1.370,
        cb: 1.0,
        }),
        // C_AR-C_1 single (aryl to nitrile C) — RDKit verbose-extracted
        (MMFFAtomType::C_AR, MMFFAtomType::C_1, BondType::Single)
        | (MMFFAtomType::C_1, MMFFAtomType::C_AR, BondType::Single) => Some(BondParams {
            k_bond: 5.445, r0: 1.424,
        cb: 1.0,
        }),
        (MMFFAtomType::C_3, MMFFAtomType::P_4, BondType::Single)
        | (MMFFAtomType::P_4, MMFFAtomType::C_3, BondType::Single) => Some(BondParams {
            k_bond: 2.980,
            r0: 1.810,
        cb: 1.0,
        }),
        (MMFFAtomType::P_4, MMFFAtomType::O_3, BondType::Single)
        | (MMFFAtomType::O_3, MMFFAtomType::P_4, BondType::Single) => Some(BondParams {
            k_bond: 5.243,
            r0: 1.630,
        cb: 1.0,
        }),
        (MMFFAtomType::P_4, MMFFAtomType::O_CO2, BondType::Double)
        | (MMFFAtomType::O_CO2, MMFFAtomType::P_4, BondType::Double) => Some(BondParams {
            k_bond: 8.296,
            r0: 1.510,
        cb: 1.0,
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
    let cs = -2.0 * params.cb;
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
    let cs = -2.0 * params.cb;
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
            cb: 1.0,
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
            cb: 1.0,
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
            cb: 1.0,
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
            cb: 1.0,
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
