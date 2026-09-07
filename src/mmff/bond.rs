//! Bond stretching term for MMFF94

use super::params::{get_mmff_bond_type, mmff_type_id};
use super::MMFFAtomType;
use crate::molecule::BondType;

/// Bond stretching parameters
#[derive(Debug, Clone, Copy)]
pub struct BondParams {
    pub k_bond: f64, // mdyn/Å²
    pub r0: f64,     // Å
    pub cb: f64,     // cubic stretch constant (cs = -2*cb); 1.0 for most bonds
}

/// MMFF94 uses a single universal cubic stretch constant `cs = -2.0` for every
/// bond (Halgren 1996 eq. 3); `cb = -cs/2 = 1.0` everywhere. Verified empirically
/// against RDKit 2025.09.3 by stretching isolated bonds (C-H, C-F, C-Cl, C-O,
/// C-N, C-C, C=C, C#N, P=O) — every one yields `cs = -2.0000`. So there are NO
/// per-bond `cb` values to override; this table stays empty by design.
const CB_OVERRIDES: &[(u8, u8, f64)] = &[];

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
        return Some(BondParams {
            k_bond: kb,
            r0,
            cb: 1.0,
        });
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
    // The complete MMFF94 release table takes priority over any hand-written
    // fallbacks below (which only exist as a safety net for exotic pairs).
    let class = get_mmff_bond_type(bond_type, mmff_type_id(type1), mmff_type_id(type2));
    if let Some((k_bond, r0)) = lookup_mmff94_bond(class, mmff_type_id(type1), mmff_type_id(type2))
    {
        return Some(BondParams {
            k_bond,
            r0,
            cb: 1.0,
        });
    }

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
        // C_2-CE4R (sp2 C in 4-ring) single bond
        (MMFFAtomType::C_2, MMFFAtomType::CE4R, BondType::Single)
        | (MMFFAtomType::CE4R, MMFFAtomType::C_2, BondType::Single) => Some(BondParams {
            k_bond: 4.481,
            r0: 1.471,
            cb: 1.0,
        }),
        (MMFFAtomType::C_AR, MMFFAtomType::C_AR, BondType::Single) => Some(BondParams {
            k_bond: 5.178,
            r0: 1.436,
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
            k_bond: 5.827362192285844,
            r0: 1.594267567692851,
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
        // O_R-N_NITROSO (nitrite O-N) — from RDKit verbose
        (MMFFAtomType::O_R, MMFFAtomType::N_NITROSO, BondType::Single)
        | (MMFFAtomType::N_NITROSO, MMFFAtomType::O_R, BondType::Single) => Some(BondParams {
            k_bond: 3.971,
            r0: 1.424,
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
        // N_RAD (nitrene N, type 62) bonds
        (MMFFAtomType::N_RAD, MMFFAtomType::C_3, BondType::Single)
        | (MMFFAtomType::C_3, MMFFAtomType::N_RAD, BondType::Single) => Some(BondParams {
            k_bond: 4.456,
            r0: 1.444,
            cb: 1.0,
        }),
        (MMFFAtomType::N_RAD, MMFFAtomType::H_N3, BondType::Single)
        | (MMFFAtomType::H_N3, MMFFAtomType::N_RAD, BondType::Single) => Some(BondParams {
            k_bond: 6.339,
            r0: 1.026,
            cb: 1.0,
        }),
        (MMFFAtomType::C_3, MMFFAtomType::N_2, BondType::Single)
        | (MMFFAtomType::N_2, MMFFAtomType::C_3, BondType::Single) => Some(BondParams {
            k_bond: 4.763,
            r0: 1.458,
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
            k_bond: 2.982332971118793,
            r0: 1.798344875466168,
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
        // Pyrrolide anion ring bonds (NPYL_M=76, C5A_M=78)
        (MMFFAtomType::NPYL_M, MMFFAtomType::C5A_M, BondType::Single)
        | (MMFFAtomType::C5A_M, MMFFAtomType::NPYL_M, BondType::Single)
        | (MMFFAtomType::NPYL_M, MMFFAtomType::C5A_M, BondType::Double)
        | (MMFFAtomType::C5A_M, MMFFAtomType::NPYL_M, BondType::Double)
        | (MMFFAtomType::NPYL_M, MMFFAtomType::C5A_M, BondType::Aromatic)
        | (MMFFAtomType::C5A_M, MMFFAtomType::NPYL_M, BondType::Aromatic) => Some(BondParams {
            k_bond: 6.824,
            r0: 1.345,
            cb: 1.0,
        }),
        // N_PYR (N-methylpyridinium N, type 58) bonds
        (MMFFAtomType::N_PYR, MMFFAtomType::C_3, BondType::Single)
        | (MMFFAtomType::C_3, MMFFAtomType::N_PYR, BondType::Single) => Some(BondParams {
            k_bond: 4.329,
            r0: 1.451,
            cb: 1.0,
        }),
        (MMFFAtomType::N_PYR, MMFFAtomType::C_AR, BondType::Aromatic)
        | (MMFFAtomType::C_AR, MMFFAtomType::N_PYR, BondType::Aromatic) => Some(BondParams {
            k_bond: 7.432,
            r0: 1.326,
            cb: 1.0,
        }),
        // P_ARM (aromatic phosphirene P, type 75) bonds
        (MMFFAtomType::P_ARM, MMFFAtomType::C_2, BondType::Double)
        | (MMFFAtomType::C_2, MMFFAtomType::P_ARM, BondType::Double) => Some(BondParams {
            k_bond: 4.191,
            r0: 1.710,
            cb: 1.0,
        }),
        (MMFFAtomType::P_ARM, MMFFAtomType::CE4R, BondType::Single)
        | (MMFFAtomType::CE4R, MMFFAtomType::P_ARM, BondType::Single) => Some(BondParams {
            k_bond: 2.761836,
            r0: 1.833069,
            cb: 1.0,
        }),
        // O_2P (furanium O+, type 51) bonds
        (MMFFAtomType::O_2P, MMFFAtomType::C_2, BondType::Double)
        | (MMFFAtomType::C_2, MMFFAtomType::O_2P, BondType::Double) => Some(BondParams {
            k_bond: 8.562,
            r0: 1.290,
            cb: 1.0,
        }),
        (MMFFAtomType::O_2P, MMFFAtomType::CE4R, BondType::Single)
        | (MMFFAtomType::CE4R, MMFFAtomType::O_2P, BondType::Single) => Some(BondParams {
            k_bond: 5.129116,
            r0: 1.405,
            cb: 1.0,
        }),
        // N_T3 (trimethylamine N-oxide N, type 68) bonds
        (MMFFAtomType::N_T3, MMFFAtomType::C_3, BondType::Single)
        | (MMFFAtomType::C_3, MMFFAtomType::N_T3, BondType::Single) => Some(BondParams {
            k_bond: 4.217,
            r0: 1.479,
            cb: 1.0,
        }),
        (MMFFAtomType::N_T3, MMFFAtomType::O_CO2, BondType::Single)
        | (MMFFAtomType::O_CO2, MMFFAtomType::N_T3, BondType::Single) => Some(BondParams {
            k_bond: 4.398,
            r0: 1.348,
            cb: 1.0,
        }),
        (MMFFAtomType::C5A_M, MMFFAtomType::C5A_M, BondType::Single)
        | (MMFFAtomType::C5A_M, MMFFAtomType::C5A_M, BondType::Double)
        | (MMFFAtomType::C5A_M, MMFFAtomType::C5A_M, BondType::Aromatic) => Some(BondParams {
            k_bond: 5.573,
            r0: 1.374,
            cb: 1.0,
        }),
        (MMFFAtomType::C5A_M, MMFFAtomType::H, BondType::Single)
        | (MMFFAtomType::H, MMFFAtomType::C5A_M, BondType::Single) => Some(BondParams {
            k_bond: 5.506,
            r0: 1.080,
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
        // C_1≡NID (nitrile/nitrile-oxide) — empirical rule from RDKit verbose
        (MMFFAtomType::C_1, MMFFAtomType::NID, BondType::Triple)
        | (MMFFAtomType::NID, MMFFAtomType::C_1, BondType::Triple) => Some(BondParams {
            k_bond: 4.148863731428771,
            r0: 1.4613059914169322,
            cb: 1.0,
        }),
        // NID-OXIDE (nitrile oxide N-O)
        (MMFFAtomType::NID, MMFFAtomType::OXIDE, BondType::Single)
        | (MMFFAtomType::OXIDE, MMFFAtomType::NID, BondType::Single) => Some(BondParams {
            k_bond: 3.971145145803896,
            r0: 1.4239219710253257,
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
        // CR3R-P_3 (3-ring P-C, e.g. phosphirane) — was falling back to C_3-P_3 (2.79/1.83);
        // RDKit's 3-ring-specific value. Found via full-precision param audit.
        (MMFFAtomType::CR3R, MMFFAtomType::P_3, BondType::Single)
        | (MMFFAtomType::P_3, MMFFAtomType::CR3R, BondType::Single) => Some(BondParams {
            k_bond: 2.7618,
            r0: 1.8331,
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
        // CR4R-N_3 (4-ring N) — from RDKit GetMMFFBondStretchParams
        (MMFFAtomType::CR4R, MMFFAtomType::N_3, BondType::Single)
        | (MMFFAtomType::N_3, MMFFAtomType::CR4R, BondType::Single) => Some(BondParams {
            k_bond: 5.107,
            r0: 1.456,
            cb: 1.0,
        }),
        // CE4R (sp2 C in 4-ring, type 30) bonds
        (MMFFAtomType::CE4R, MMFFAtomType::CE4R, BondType::Double)
        | (MMFFAtomType::CE4R, MMFFAtomType::CE4R, BondType::Aromatic) => Some(BondParams {
            k_bond: 9.579,
            r0: 1.343,
            cb: 1.0,
        }),
        (MMFFAtomType::CE4R, MMFFAtomType::CR4R, BondType::Single)
        | (MMFFAtomType::CR4R, MMFFAtomType::CE4R, BondType::Single) => Some(BondParams {
            k_bond: 3.977,
            r0: 1.507,
            cb: 1.0,
        }),
        (MMFFAtomType::CE4R, MMFFAtomType::H, BondType::Single)
        | (MMFFAtomType::H, MMFFAtomType::CE4R, BondType::Single) => Some(BondParams {
            k_bond: 5.176,
            r0: 1.086,
            cb: 1.0,
        }),
        // O_3P (oxonium O+, type 49) bonds
        (MMFFAtomType::O_3P, MMFFAtomType::C_3, BondType::Single)
        | (MMFFAtomType::C_3, MMFFAtomType::O_3P, BondType::Single) => Some(BondParams {
            k_bond: 5.129115902527102,
            r0: 1.405,
            cb: 1.0,
        }),
        (MMFFAtomType::O_3P, MMFFAtomType::H_OXP, BondType::Single)
        | (MMFFAtomType::H_OXP, MMFFAtomType::O_3P, BondType::Single) => Some(BondParams {
            k_bond: 6.812,
            r0: 0.991,
            cb: 1.0,
        }),
        // S_O3 (sulfite S, type 73) bonds
        (MMFFAtomType::S_O3, MMFFAtomType::C_3, BondType::Single)
        | (MMFFAtomType::C_3, MMFFAtomType::S_O3, BondType::Single) => Some(BondParams {
            k_bond: 2.608,
            r0: 1.839,
            cb: 1.0,
        }),
        (MMFFAtomType::S_O3, MMFFAtomType::O_CO2, BondType::Single)
        | (MMFFAtomType::O_CO2, MMFFAtomType::S_O3, BondType::Single)
        | (MMFFAtomType::S_O3, MMFFAtomType::O_CO2, BondType::Double)
        | (MMFFAtomType::O_CO2, MMFFAtomType::S_O3, BondType::Double) => Some(BondParams {
            k_bond: 8.427,
            r0: 1.510,
            cb: 1.0,
        }),
        // S_CSO (sulfene S, type 74) bonds
        (MMFFAtomType::S_CSO, MMFFAtomType::C_2, BondType::Double)
        | (MMFFAtomType::C_2, MMFFAtomType::S_CSO, BondType::Double) => Some(BondParams {
            k_bond: 5.204,
            r0: 1.639,
            cb: 1.0,
        }),
        (MMFFAtomType::S_CSO, MMFFAtomType::O_2, BondType::Double)
        | (MMFFAtomType::O_2, MMFFAtomType::S_CSO, BondType::Double) => Some(BondParams {
            k_bond: 9.129,
            r0: 1.490,
            cb: 1.0,
        }),
        // CL4 (perchlorate Cl, type 77) bonds
        (MMFFAtomType::CL4, MMFFAtomType::O_CO2, BondType::Single)
        | (MMFFAtomType::O_CO2, MMFFAtomType::CL4, BondType::Single) => Some(BondParams {
            k_bond: 10.648,
            r0: 1.450,
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
            k_bond: 5.633,
            r0: 1.076,
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
            k_bond: 3.378,
            r0: 1.721,
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
            k_bond: 6.570,
            r0: 1.340,
            cb: 1.0,
        }),
        // N-N bonds
        (MMFFAtomType::N_3, MMFFAtomType::N_3, BondType::Single) => Some(BondParams {
            k_bond: 3.5,
            r0: 1.45,
            cb: 1.0,
        }),
        (MMFFAtomType::N_2, MMFFAtomType::N_2, BondType::Double) => Some(BondParams {
            k_bond: 7.256,
            r0: 1.243,
            cb: 1.0,
        }),
        // N_2-C_AR Single (azo N to aromatic C)
        (MMFFAtomType::N_2, MMFFAtomType::C_AR, BondType::Single)
        | (MMFFAtomType::C_AR, MMFFAtomType::N_2, BondType::Single) => Some(BondParams {
            k_bond: 5.529,
            r0: 1.393,
            cb: 1.0,
        }),
        // C_2-S_3 Single (thioester/thioamide C-S)
        (MMFFAtomType::C_2, MMFFAtomType::S_3, BondType::Single)
        | (MMFFAtomType::S_3, MMFFAtomType::C_2, BondType::Single) => Some(BondParams {
            k_bond: 3.536,
            r0: 1.748,
            cb: 1.0,
        }),
        // N_PL3-H_NAM (amidine H on N_PL3); plain H variants are covered by the earlier N_PL3-H arm
        (MMFFAtomType::N_PL3, MMFFAtomType::H_NAM, BondType::Single)
        | (MMFFAtomType::H_NAM, MMFFAtomType::N_PL3, BondType::Single) => Some(BondParams {
            k_bond: 6.576,
            r0: 1.018,
            cb: 1.0,
        }),
        // N_3-N_AM (hydrazide N-N) and N_3-H (amine N-H)
        (MMFFAtomType::N_3, MMFFAtomType::N_AM, BondType::Single)
        | (MMFFAtomType::N_AM, MMFFAtomType::N_3, BondType::Single) => Some(BondParams {
            k_bond: 3.909,
            r0: 1.378,
            cb: 1.0,
        }),
        // S_3-C_AR (thioether to aromatic C)
        (MMFFAtomType::S_3, MMFFAtomType::C_AR, BondType::Single)
        | (MMFFAtomType::C_AR, MMFFAtomType::S_3, BondType::Single) => Some(BondParams {
            k_bond: 3.565,
            r0: 1.765,
            cb: 1.0,
        }),
        (MMFFAtomType::N_3, MMFFAtomType::H_N3, BondType::Single)
        | (MMFFAtomType::H_N3, MMFFAtomType::N_3, BondType::Single) => Some(BondParams {
            k_bond: 6.490,
            r0: 1.019,
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
            k_bond: 5.468,
            r0: 1.423,
            cb: 1.0,
        }),
        (MMFFAtomType::C5A, MMFFAtomType::NPYL, BondType::Aromatic)
        | (MMFFAtomType::NPYL, MMFFAtomType::C5A, BondType::Aromatic)
        | (MMFFAtomType::C5A, MMFFAtomType::NPYL, BondType::Single)
        | (MMFFAtomType::NPYL, MMFFAtomType::C5A, BondType::Single) => Some(BondParams {
            k_bond: 6.301,
            r0: 1.364,
            cb: 1.0,
        }),
        (MMFFAtomType::C5A, MMFFAtomType::N5B, BondType::Aromatic)
        | (MMFFAtomType::N5B, MMFFAtomType::C5A, BondType::Aromatic)
        | (MMFFAtomType::C5A, MMFFAtomType::N5B, BondType::Single)
        | (MMFFAtomType::N5B, MMFFAtomType::C5A, BondType::Single) => Some(BondParams {
            k_bond: 8.326,
            r0: 1.313,
            cb: 1.0,
        }),
        (MMFFAtomType::C5B, MMFFAtomType::C5A, BondType::Aromatic)
        | (MMFFAtomType::C5A, MMFFAtomType::C5B, BondType::Aromatic)
        | (MMFFAtomType::C5B, MMFFAtomType::C5A, BondType::Single)
        | (MMFFAtomType::C5A, MMFFAtomType::C5B, BondType::Single) => Some(BondParams {
            k_bond: 7.118,
            r0: 1.377,
            cb: 1.0,
        }),
        (MMFFAtomType::N5B, MMFFAtomType::C5B, BondType::Aromatic)
        | (MMFFAtomType::C5B, MMFFAtomType::N5B, BondType::Aromatic)
        | (MMFFAtomType::N5B, MMFFAtomType::C5B, BondType::Single)
        | (MMFFAtomType::C5B, MMFFAtomType::N5B, BondType::Single) => Some(BondParams {
            k_bond: 4.456,
            r0: 1.369,
            cb: 1.0,
        }),
        (MMFFAtomType::N_AM, MMFFAtomType::C_2, BondType::Single)
        | (MMFFAtomType::C_2, MMFFAtomType::N_AM, BondType::Single) => Some(BondParams {
            k_bond: 5.829,
            r0: 1.369,
            cb: 1.0,
        }),
        (MMFFAtomType::C5B, MMFFAtomType::N_AM, BondType::Single)
        | (MMFFAtomType::N_AM, MMFFAtomType::C5B, BondType::Single) => Some(BondParams {
            k_bond: 5.952,
            r0: 1.376,
            cb: 1.0,
        }),
        (MMFFAtomType::NPYL, MMFFAtomType::H_N3, BondType::Single)
        | (MMFFAtomType::H_N3, MMFFAtomType::NPYL, BondType::Single) => Some(BondParams {
            k_bond: 7.112,
            r0: 1.012,
            cb: 1.0,
        }),
        // Purine 6-ring N_PL3-C bonds (RDKit verbose-estimated)
        (MMFFAtomType::N_PL3, MMFFAtomType::C_2, BondType::Single)
        | (MMFFAtomType::C_2, MMFFAtomType::N_PL3, BondType::Single) => Some(BondParams {
            k_bond: 6.110,
            r0: 1.370,
            cb: 1.0,
        }),
        (MMFFAtomType::N_PL3, MMFFAtomType::C_VIN, BondType::Single)
        | (MMFFAtomType::C_VIN, MMFFAtomType::N_PL3, BondType::Single) => Some(BondParams {
            k_bond: 6.110,
            r0: 1.370,
            cb: 1.0,
        }),
        // C_AR-C_1 single (aryl to nitrile C) — RDKit verbose-extracted
        (MMFFAtomType::C_AR, MMFFAtomType::C_1, BondType::Single)
        | (MMFFAtomType::C_1, MMFFAtomType::C_AR, BondType::Single) => Some(BondParams {
            k_bond: 5.445,
            r0: 1.424,
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
        // P-C_AR (aryl phosphines / phosphine oxides) — was missing, fell back to
        // estimation and overstated bond energy by ~+9 kcal on Ph3P=O. RDKit values.
        (MMFFAtomType::P_4, MMFFAtomType::C_AR, BondType::Single)
        | (MMFFAtomType::C_AR, MMFFAtomType::P_4, BondType::Single) => Some(BondParams {
            k_bond: 3.586,
            r0: 1.755,
            cb: 1.0,
        }),
        (MMFFAtomType::P_3, MMFFAtomType::C_AR, BondType::Single)
        | (MMFFAtomType::C_AR, MMFFAtomType::P_3, BondType::Single) => Some(BondParams {
            k_bond: 3.207,
            r0: 1.788,
            cb: 1.0,
        }),

        // === val_set_new5 bond params ===
        // C_AR-O_2P (furanium ring O+=)
        // N_AM-H_NAM and N_AR-H_N3 bonds (indole/tryptophan NH); plain H variants are covered by the earlier N-H arms
        (MMFFAtomType::N_AM, MMFFAtomType::H_NAM, BondType::Single)
        | (MMFFAtomType::H_NAM, MMFFAtomType::N_AM, BondType::Single) => Some(BondParams {
            k_bond: 6.663,
            r0: 1.015,
            cb: 1.0,
        }),
        (MMFFAtomType::N_AR, MMFFAtomType::H_N3, BondType::Single)
        | (MMFFAtomType::H_N3, MMFFAtomType::N_AR, BondType::Single) => Some(BondParams {
            k_bond: 7.112,
            r0: 1.012,
            cb: 1.0,
        }),
        // N_AM-O_3 bond (hydroxamic acid N-O)
        (MMFFAtomType::N_AM, MMFFAtomType::O_3, BondType::Single)
        | (MMFFAtomType::O_3, MMFFAtomType::N_AM, BondType::Single) => Some(BondParams {
            k_bond: 5.982,
            r0: 1.410,
            cb: 1.0,
        }),
        // C_AR-O_2P (furanium ring O+=)
        (MMFFAtomType::C_AR, MMFFAtomType::O_2P, BondType::Aromatic)
        | (MMFFAtomType::O_2P, MMFFAtomType::C_AR, BondType::Aromatic) => Some(BondParams {
            k_bond: 5.129116,
            r0: 1.405,
            cb: 1.0,
        }),
        // O_2P (51) / H_OXP2 (52) bonds — oxenium
        (MMFFAtomType::O_2P, MMFFAtomType::H_OXP2, BondType::Single)
        | (MMFFAtomType::H_OXP2, MMFFAtomType::O_2P, BondType::Single) => Some(BondParams {
            k_bond: 7.100,
            r0: 0.987,
            cb: 1.0,
        }),
        // N5 (79) bonds — general 5-ring N
        (MMFFAtomType::C5A_M, MMFFAtomType::N5, BondType::Aromatic)
        | (MMFFAtomType::N5, MMFFAtomType::C5A_M, BondType::Aromatic) => Some(BondParams {
            k_bond: 8.890,
            r0: 1.287,
            cb: 1.0,
        }),
        (MMFFAtomType::N5, MMFFAtomType::N_5POS, BondType::Aromatic)
        | (MMFFAtomType::N_5POS, MMFFAtomType::N5, BondType::Aromatic) => Some(BondParams {
            k_bond: 4.305,
            r0: 1.356,
            cb: 1.0,
        }),
        (MMFFAtomType::C5B, MMFFAtomType::N5, BondType::Aromatic)
        | (MMFFAtomType::N5, MMFFAtomType::C5B, BondType::Aromatic) => Some(BondParams {
            k_bond: 4.148864,
            r0: 1.461306,
            cb: 1.0,
        }),
        (MMFFAtomType::N_5POS, MMFFAtomType::N_5POS, BondType::Aromatic) => Some(BondParams {
            k_bond: 2.763084,
            r0: 1.460,
            cb: 1.0,
        }),
        // HNRP-N_5POS bond (H on charged 5-ring N)
        (MMFFAtomType::HNRP, MMFFAtomType::N_5POS, BondType::Single)
        | (MMFFAtomType::N_5POS, MMFFAtomType::HNRP, BondType::Single)
        | (MMFFAtomType::H, MMFFAtomType::N_5POS, BondType::Single)
        | (MMFFAtomType::N_5POS, MMFFAtomType::H, BondType::Single) => Some(BondParams {
            k_bond: 6.980,
            r0: 1.016,
            cb: 1.0,
        }),

        // === val_set_new4 bond params ===
        // N_GD (56) bonds — guanidinium
        (MMFFAtomType::C_3, MMFFAtomType::N_GD, BondType::Single)
        | (MMFFAtomType::N_GD, MMFFAtomType::C_3, BondType::Single) => Some(BondParams {
            k_bond: 4.166,
            r0: 1.453,
            cb: 1.0,
        }),
        (MMFFAtomType::N_GD, MMFFAtomType::C_AN, BondType::Single)
        | (MMFFAtomType::C_AN, MMFFAtomType::N_GD, BondType::Single)
        | (MMFFAtomType::N_GD, MMFFAtomType::C_AN, BondType::Double)
        | (MMFFAtomType::C_AN, MMFFAtomType::N_GD, BondType::Double) => Some(BondParams {
            k_bond: 4.137,
            r0: 1.383,
            cb: 1.0,
        }),
        (MMFFAtomType::HNRP, MMFFAtomType::N_GD, BondType::Single)
        | (MMFFAtomType::N_GD, MMFFAtomType::HNRP, BondType::Single) => Some(BondParams {
            k_bond: 6.490,
            r0: 1.017,
            cb: 1.0,
        }),
        // N_5OX2 (82) bonds — isoxazole N-oxide 5-ring
        (MMFFAtomType::O_CO2, MMFFAtomType::N_5OX2, BondType::Single)
        | (MMFFAtomType::N_5OX2, MMFFAtomType::O_CO2, BondType::Single)
        | (MMFFAtomType::O_CO2, MMFFAtomType::N_5OX2, BondType::Double)
        | (MMFFAtomType::N_5OX2, MMFFAtomType::O_CO2, BondType::Double) => Some(BondParams {
            k_bond: 8.594,
            r0: 1.252,
            cb: 1.0,
        }),
        (MMFFAtomType::C5B, MMFFAtomType::N_5OX2, BondType::Aromatic)
        | (MMFFAtomType::N_5OX2, MMFFAtomType::C5B, BondType::Aromatic) => Some(BondParams {
            k_bond: 6.794,
            r0: 1.346,
            cb: 1.0,
        }),
        (MMFFAtomType::C5A, MMFFAtomType::N_5OX2, BondType::Aromatic)
        | (MMFFAtomType::N_5OX2, MMFFAtomType::C5A, BondType::Aromatic) => Some(BondParams {
            k_bond: 4.1489,
            r0: 1.4613,
            cb: 1.0,
        }),
        // N_5OX (67) bonds — pyridine N-oxide
        (MMFFAtomType::O_2, MMFFAtomType::N_5OX, BondType::Double)
        | (MMFFAtomType::N_5OX, MMFFAtomType::O_2, BondType::Double) => Some(BondParams {
            k_bond: 3.971145,
            r0: 1.423922,
            cb: 1.0,
        }),
        (MMFFAtomType::N_5OX, MMFFAtomType::C_VIN, BondType::Single)
        | (MMFFAtomType::C_VIN, MMFFAtomType::N_5OX, BondType::Single) => Some(BondParams {
            k_bond: 4.685,
            r0: 1.432,
            cb: 1.0,
        }),
        // N_IM (54) bonds — iminium
        (MMFFAtomType::C_2, MMFFAtomType::N_IM, BondType::Double)
        | (MMFFAtomType::N_IM, MMFFAtomType::C_2, BondType::Double) => Some(BondParams {
            k_bond: 10.333,
            r0: 1.280,
            cb: 1.0,
        }),
        (MMFFAtomType::HNRP, MMFFAtomType::N_IM, BondType::Single)
        | (MMFFAtomType::N_IM, MMFFAtomType::HNRP, BondType::Single) => Some(BondParams {
            k_bond: 6.529,
            r0: 1.022,
            cb: 1.0,
        }),
        // N_SO (48) bonds — sulfinylamine
        (MMFFAtomType::S_O2, MMFFAtomType::N_SO, BondType::Double)
        | (MMFFAtomType::N_SO, MMFFAtomType::S_O2, BondType::Double) => Some(BondParams {
            k_bond: 6.186,
            r0: 1.540,
            cb: 1.0,
        }),
        (MMFFAtomType::S_O2, MMFFAtomType::HS, BondType::Single)
        | (MMFFAtomType::HS, MMFFAtomType::S_O2, BondType::Single) => Some(BondParams {
            k_bond: 3.806451,
            r0: 1.353219,
            cb: 1.0,
        }),
        (MMFFAtomType::H_NAM, MMFFAtomType::N_SO, BondType::Single)
        | (MMFFAtomType::N_SO, MMFFAtomType::H_NAM, BondType::Single) => Some(BondParams {
            k_bond: 6.413,
            r0: 1.024,
            cb: 1.0,
        }),
        // N_5POS (81) / C5A_M (78) / C_IM (80) bonds — imidazolium
        (MMFFAtomType::C_3, MMFFAtomType::N_5POS, BondType::Single)
        | (MMFFAtomType::N_5POS, MMFFAtomType::C_3, BondType::Single) => Some(BondParams {
            k_bond: 4.512,
            r0: 1.441,
            cb: 1.0,
        }),
        (MMFFAtomType::H_NAM, MMFFAtomType::C5A_M, BondType::Single)
        | (MMFFAtomType::C5A_M, MMFFAtomType::H_NAM, BondType::Single) => Some(BondParams {
            k_bond: 5.506,
            r0: 1.080,
            cb: 1.0,
        }),
        (MMFFAtomType::H_NAM, MMFFAtomType::C_IM, BondType::Single)
        | (MMFFAtomType::C_IM, MMFFAtomType::H_NAM, BondType::Single)
        | (MMFFAtomType::H, MMFFAtomType::C_IM, BondType::Single)
        | (MMFFAtomType::C_IM, MMFFAtomType::H, BondType::Single) => Some(BondParams {
            k_bond: 5.633,
            r0: 1.076,
            cb: 1.0,
        }),
        (MMFFAtomType::C5A_M, MMFFAtomType::N_5POS, BondType::Aromatic)
        | (MMFFAtomType::N_5POS, MMFFAtomType::C5A_M, BondType::Aromatic) => Some(BondParams {
            k_bond: 5.046,
            r0: 1.381,
            cb: 1.0,
        }),
        (MMFFAtomType::C_IM, MMFFAtomType::N_5POS, BondType::Aromatic)
        | (MMFFAtomType::N_5POS, MMFFAtomType::C_IM, BondType::Aromatic) => Some(BondParams {
            k_bond: 8.237,
            r0: 1.335,
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

// Complete MMFF94 bond-stretch parameter table (493 rows), generated from
// the MMFF94 release tables embedded in RDKit's Code/ForceField/MMFF/Params.cpp
// (`defaultMMFFBond`). Row: (bond class, type A, type B, kb mdyn/Å, r0 Å).
// Types stored as (min, max) — the MMFF tables are symmetric.

#[allow(clippy::approx_constant)] // 6.283 is the release-table value
pub const MMFF94_BOND_TABLE: &[(u8, u8, u8, f64, f64)] = &[
    (0, 1, 1, 4.258, 1.508),
    (0, 1, 2, 4.539, 1.482),
    (0, 1, 3, 4.19, 1.492),
    (0, 1, 4, 4.707, 1.459),
    (0, 1, 5, 4.766, 1.093),
    (0, 1, 6, 5.047, 1.418),
    (0, 1, 8, 5.084, 1.451),
    (0, 1, 9, 4.763, 1.458),
    (0, 1, 10, 4.664, 1.436),
    (0, 1, 11, 6.011, 1.360),
    (0, 1, 12, 2.974, 1.773),
    (0, 1, 13, 2.529, 1.949),
    (0, 1, 14, 1.706, 2.090),
    (0, 1, 15, 2.893, 1.805),
    (0, 1, 17, 2.841, 1.813),
    (0, 1, 18, 3.258, 1.772),
    (0, 1, 19, 2.866, 1.830),
    (0, 1, 20, 4.65, 1.504),
    (0, 1, 22, 4.286, 1.482),
    (0, 1, 25, 2.98, 1.810),
    (0, 1, 26, 2.79, 1.830),
    (0, 1, 34, 3.844, 1.480),
    (0, 1, 35, 7.915, 1.307),
    (0, 1, 37, 4.957, 1.486),
    (0, 1, 39, 6.114, 1.445),
    (0, 1, 40, 4.922, 1.446),
    (0, 1, 41, 3.83, 1.510),
    (0, 1, 43, 3.971, 1.472),
    (0, 1, 45, 3.844, 1.480),
    (0, 1, 46, 3.813, 1.482),
    (0, 1, 54, 4.267, 1.461),
    (0, 1, 55, 4.646, 1.454),
    (0, 1, 56, 4.166, 1.453),
    (0, 1, 57, 4.669, 1.461),
    (0, 1, 58, 4.329, 1.451),
    (0, 1, 61, 4.845, 1.424),
    (0, 1, 62, 4.456, 1.444),
    (0, 1, 63, 4.481, 1.471),
    (0, 1, 64, 4.518, 1.469),
    (0, 1, 67, 4.188, 1.459),
    (0, 1, 68, 4.217, 1.479),
    (0, 1, 72, 2.956, 1.801),
    (0, 1, 73, 2.608, 1.839),
    (0, 1, 75, 2.547, 1.858),
    (0, 1, 78, 4.593, 1.465),
    (0, 1, 80, 4.373, 1.477),
    (0, 1, 81, 4.512, 1.441),
    (0, 2, 2, 9.505, 1.333),
    (0, 2, 4, 9.538, 1.297),
    (0, 2, 5, 5.17, 1.083),
    (0, 2, 6, 5.52, 1.373),
    (0, 2, 10, 6.329, 1.362),
    (0, 2, 11, 6.283, 1.350),
    (0, 2, 12, 3.39, 1.720),
    (0, 2, 13, 3.413, 1.854),
    (0, 2, 14, 2.062, 2.025),
    (0, 2, 15, 3.896, 1.720),
    (0, 2, 17, 3.247, 1.773),
    (0, 2, 18, 3.789, 1.728),
    (0, 2, 19, 3.052, 1.811),
    (0, 2, 20, 4.593, 1.465),
    (0, 2, 22, 4.926, 1.448),
    (0, 2, 25, 3.75, 1.742),
    (0, 2, 30, 8.166, 1.331),
    (0, 2, 34, 5.207, 1.407),
    (0, 2, 35, 10.343, 1.250),
    (0, 2, 40, 6.11, 1.370),
    (0, 2, 41, 3.746, 1.505),
    (0, 2, 43, 4.928, 1.420),
    (0, 2, 45, 4.725, 1.430),
    (0, 2, 46, 7.466, 1.325),
    (0, 2, 55, 6.164, 1.368),
    (0, 2, 56, 6.246, 1.365),
    (0, 2, 62, 7.105, 1.336),
    (0, 2, 72, 4.179, 1.700),
    (0, 3, 5, 4.65, 1.101),
    (0, 3, 6, 5.801, 1.355),
    (0, 3, 7, 12.95, 1.222),
    (0, 3, 9, 10.077, 1.290),
    (0, 3, 10, 5.829, 1.369),
    (0, 3, 11, 6.57, 1.340),
    (0, 3, 12, 3.449, 1.715),
    (0, 3, 15, 3.536, 1.748),
    (0, 3, 16, 4.735, 1.665),
    (0, 3, 17, 2.888, 1.808),
    (0, 3, 18, 3.394, 1.760),
    (0, 3, 20, 3.298, 1.530),
    (0, 3, 22, 4.593, 1.465),
    (0, 3, 25, 3.164, 1.792),
    (0, 3, 35, 11.012, 1.237),
    (0, 3, 40, 6.11, 1.370),
    (0, 3, 41, 4.286, 1.482),
    (0, 3, 43, 4.928, 1.420),
    (0, 3, 45, 4.531, 1.440),
    (0, 3, 48, 5.412, 1.398),
    (0, 3, 51, 8.562, 1.290),
    (0, 3, 53, 7.637, 1.320),
    (0, 3, 54, 10.333, 1.280),
    (0, 3, 55, 4.886, 1.422),
    (0, 3, 56, 4.907, 1.421),
    (0, 3, 62, 7.568, 1.322),
    (0, 3, 67, 8.217, 1.304),
    (0, 3, 74, 5.204, 1.639),
    (0, 3, 75, 4.191, 1.710),
    (0, 4, 4, 15.206, 1.200),
    (0, 4, 5, 5.726, 1.065),
    (0, 4, 6, 7.193, 1.328),
    (0, 4, 7, 14.916, 1.176),
    (0, 4, 9, 15.589, 1.172),
    (0, 4, 10, 6.824, 1.345),
    (0, 4, 15, 4.33, 1.690),
    (0, 4, 20, 5.178, 1.436),
    (0, 4, 22, 5.4, 1.426),
    (0, 4, 30, 10.227, 1.282),
    (0, 4, 42, 16.582, 1.160),
    (0, 4, 43, 6.947, 1.341),
    (0, 5, 19, 2.254, 1.485),
    (0, 5, 20, 4.852, 1.093),
    (0, 5, 22, 5.191, 1.082),
    (0, 5, 30, 5.176, 1.086),
    (0, 5, 37, 5.306, 1.084),
    (0, 5, 41, 3.256, 1.144),
    (0, 5, 57, 5.633, 1.076),
    (0, 5, 63, 5.531, 1.080),
    (0, 5, 64, 5.506, 1.080),
    (0, 5, 78, 5.506, 1.080),
    (0, 5, 80, 5.633, 1.076),
    (0, 6, 6, 4.088, 1.449),
    (0, 6, 8, 5.059, 1.450),
    (0, 6, 9, 4.491, 1.395),
    (0, 6, 10, 5.982, 1.410),
    (0, 6, 15, 4.757, 1.661),
    (0, 6, 17, 5.779, 1.608),
    (0, 6, 18, 5.326, 1.630),
    (0, 6, 19, 4.661, 1.660),
    (0, 6, 20, 5.623, 1.433),
    (0, 6, 21, 7.794, 0.972),
    (0, 6, 22, 4.556, 1.433),
    (0, 6, 24, 7.403, 0.981),
    (0, 6, 25, 5.243, 1.630),
    (0, 6, 26, 5.481, 1.618),
    (0, 6, 29, 7.839, 0.973),
    (0, 6, 30, 9.359, 1.271),
    (0, 6, 33, 7.143, 0.986),
    (0, 6, 37, 5.614, 1.376),
    (0, 6, 39, 4.629, 1.388),
    (0, 6, 40, 4.609, 1.389),
    (0, 6, 41, 6.754, 1.342),
    (0, 6, 43, 3.937, 1.426),
    (0, 6, 45, 4.321, 1.404),
    (0, 6, 54, 5.117, 1.365),
    (0, 6, 55, 4.772, 1.381),
    (0, 6, 57, 7.128, 1.330),
    (0, 6, 58, 4.792, 1.380),
    (0, 6, 63, 7.324, 1.324),
    (0, 6, 64, 6.664, 1.345),
    (0, 7, 17, 8.77, 1.500),
    (0, 7, 46, 9.329, 1.235),
    (0, 7, 74, 9.129, 1.490),
    (0, 8, 8, 3.264, 1.420),
    (0, 8, 9, 4.581, 1.342),
    (0, 8, 10, 3.909, 1.378),
    (0, 8, 12, 3.371, 1.761),
    (0, 8, 15, 4.06, 1.652),
    (0, 8, 17, 3.901, 1.663),
    (0, 8, 19, 4.254, 1.700),
    (0, 8, 20, 5.107, 1.456),
    (0, 8, 22, 4.223, 1.457),
    (0, 8, 23, 6.49, 1.019),
    (0, 8, 25, 4.629, 1.660),
    (0, 8, 26, 4.027, 1.699),
    (0, 8, 34, 3.775, 1.386),
    (0, 8, 39, 3.435, 1.408),
    (0, 8, 40, 3.71, 1.390),
    (0, 8, 43, 3.977, 1.374),
    (0, 8, 45, 4.267, 1.358),
    (0, 8, 46, 5.519, 1.301),
    (0, 8, 55, 4.229, 1.360),
    (0, 8, 56, 3.995, 1.373),
    (0, 9, 9, 7.256, 1.243),
    (0, 9, 10, 4.48, 1.347),
    (0, 9, 12, 3.635, 1.739),
    (0, 9, 15, 3.791, 1.671),
    (0, 9, 18, 4.465, 1.626),
    (0, 9, 19, 3.687, 1.741),
    (0, 9, 20, 4.401, 1.447),
    (0, 9, 25, 5.379, 1.619),
    (0, 9, 27, 6.23, 1.026),
    (0, 9, 34, 3.223, 1.423),
    (0, 9, 35, 5.095, 1.366),
    (0, 9, 40, 4.382, 1.352),
    (0, 9, 41, 5.65, 1.388),
    (0, 9, 45, 4.857, 1.329),
    (0, 9, 53, 7.291, 1.242),
    (0, 9, 54, 4.991, 1.323),
    (0, 9, 55, 3.825, 1.383),
    (0, 9, 56, 4.602, 1.341),
    (0, 9, 62, 4.749, 1.334),
    (0, 9, 67, 6.752, 1.258),
    (0, 10, 10, 3.977, 1.374),
    (0, 10, 13, 3.11, 1.878),
    (0, 10, 14, 1.967, 2.029),
    (0, 10, 15, 3.593, 1.686),
    (0, 10, 17, 3.93, 1.661),
    (0, 10, 20, 4.24, 1.456),
    (0, 10, 22, 4.97, 1.418),
    (0, 10, 25, 3.82, 1.714),
    (0, 10, 26, 3.651, 1.727),
    (0, 10, 28, 6.663, 1.015),
    (0, 10, 34, 3.96, 1.375),
    (0, 10, 35, 4.898, 1.375),
    (0, 10, 37, 5.482, 1.395),
    (0, 10, 39, 4.382, 1.352),
    (0, 10, 40, 3.841, 1.382),
    (0, 10, 41, 7.466, 1.325),
    (0, 10, 45, 3.524, 1.402),
    (0, 10, 63, 6.137, 1.369),
    (0, 10, 64, 5.952, 1.376),
    (0, 11, 20, 6.339, 1.348),
    (0, 11, 22, 5.296, 1.389),
    (0, 11, 25, 6.019, 1.583),
    (0, 11, 26, 6.204, 1.575),
    (0, 11, 37, 6.511, 1.342),
    (0, 11, 40, 4.187, 1.440),
    (0, 12, 15, 2.978, 2.031),
    (0, 12, 18, 2.808, 2.051),
    (0, 12, 19, 2.838, 2.050),
    (0, 12, 20, 2.859, 1.751),
    (0, 12, 22, 3.056, 1.750),
    (0, 12, 25, 3.063, 2.023),
    (0, 12, 26, 2.448, 2.100),
    (0, 12, 37, 3.378, 1.721),
    (0, 12, 40, 3.737, 1.731),
    (0, 12, 57, 3.714, 1.694),
    (0, 12, 63, 3.413, 1.718),
    (0, 12, 64, 3.649, 1.699),
    (0, 13, 20, 2.767, 1.920),
    (0, 13, 22, 2.928, 1.902),
    (0, 13, 37, 3.031, 1.891),
    (0, 13, 64, 3.031, 1.891),
    (0, 14, 20, 0.884, 2.332),
    (0, 14, 37, 1.781, 2.075),
    (0, 15, 15, 2.531, 2.050),
    (0, 15, 18, 2.214, 2.094),
    (0, 15, 19, 2.022, 2.146),
    (0, 15, 20, 2.757, 1.822),
    (0, 15, 22, 3.802, 1.727),
    (0, 15, 25, 2.319, 2.112),
    (0, 15, 26, 2.359, 2.106),
    (0, 15, 30, 3.75, 1.731),
    (0, 15, 37, 3.565, 1.765),
    (0, 15, 40, 3.859, 1.666),
    (0, 15, 43, 3.221, 1.717),
    (0, 15, 57, 3.993, 1.713),
    (0, 15, 63, 3.724, 1.733),
    (0, 15, 64, 3.548, 1.747),
    (0, 15, 71, 4.014, 1.341),
    (0, 17, 20, 2.397, 1.865),
    (0, 17, 22, 2.566, 1.844),
    (0, 17, 37, 3.098, 1.787),
    (0, 17, 43, 4.9, 1.601),
    (0, 18, 20, 3.172, 1.780),
    (0, 18, 22, 2.757, 1.822),
    (0, 18, 32, 10.748, 1.450),
    (0, 18, 37, 3.281, 1.770),
    (0, 18, 39, 3.504, 1.693),
    (0, 18, 43, 3.301, 1.710),
    (0, 18, 48, 6.186, 1.540),
    (0, 18, 55, 4.432, 1.628),
    (0, 18, 58, 2.568, 1.783),
    (0, 18, 62, 5.51, 1.570),
    (0, 18, 63, 3.524, 1.749),
    (0, 18, 64, 3.856, 1.723),
    (0, 18, 80, 4.15, 1.702),
    (0, 19, 20, 2.288, 1.900),
    (0, 19, 37, 3.072, 1.809),
    (0, 19, 40, 4.47, 1.686),
    (0, 19, 63, 3.219, 1.795),
    (0, 19, 75, 1.6, 2.226),
    (0, 20, 20, 3.663, 1.526),
    (0, 20, 22, 4.251, 1.484),
    (0, 20, 25, 2.718, 1.838),
    (0, 20, 26, 2.588, 1.853),
    (0, 20, 30, 3.977, 1.507),
    (0, 20, 34, 4.171, 1.460),
    (0, 20, 37, 3.74, 1.516),
    (0, 20, 40, 4.784, 1.427),
    (0, 20, 41, 4.286, 1.482),
    (0, 20, 43, 3.737, 1.487),
    (0, 20, 45, 3.844, 1.480),
    (0, 22, 22, 3.969, 1.499),
    (0, 22, 30, 3.785, 1.513),
    (0, 22, 34, 4.103, 1.464),
    (0, 22, 37, 4.481, 1.471),
    (0, 22, 40, 4.188, 1.459),
    (0, 22, 41, 5.071, 1.441),
    (0, 22, 43, 4.07, 1.466),
    (0, 22, 45, 4.311, 1.452),
    (0, 23, 39, 7.112, 1.012),
    (0, 23, 62, 6.339, 1.026),
    (0, 23, 67, 6.61, 1.019),
    (0, 23, 68, 5.899, 1.038),
    (0, 25, 25, 1.514, 2.253),
    (0, 25, 32, 8.296, 1.510),
    (0, 25, 37, 3.586, 1.755),
    (0, 25, 39, 4.37, 1.676),
    (0, 25, 40, 4.629, 1.660),
    (0, 25, 43, 3.237, 1.762),
    (0, 25, 57, 4.356, 1.699),
    (0, 25, 63, 3.711, 1.745),
    (0, 25, 71, 3.001, 1.411),
    (0, 25, 72, 3.744, 1.950),
    (0, 26, 26, 1.414, 2.279),
    (0, 26, 34, 3.395, 1.748),
    (0, 26, 37, 3.207, 1.788),
    (0, 26, 40, 4.87, 1.646),
    (0, 26, 71, 2.959, 1.415),
    (0, 28, 40, 6.576, 1.018),
    (0, 28, 43, 6.265, 1.028),
    (0, 28, 48, 6.413, 1.024),
    (0, 30, 30, 9.579, 1.343),
    (0, 30, 40, 8.447, 1.298),
    (0, 31, 70, 7.88, 0.969),
    (0, 32, 41, 9.756, 1.261),
    (0, 32, 45, 9.42, 1.233),
    (0, 32, 67, 7.926, 1.269),
    (0, 32, 68, 4.398, 1.348),
    (0, 32, 69, 6.098, 1.261),
    (0, 32, 73, 8.427, 1.510),
    (0, 32, 77, 10.648, 1.450),
    (0, 32, 82, 8.594, 1.252),
    (0, 34, 36, 6.163, 1.028),
    (0, 34, 37, 4.347, 1.450),
    (0, 34, 43, 4.401, 1.351),
    (0, 35, 37, 9.767, 1.262),
    (0, 35, 63, 12.76, 1.207),
    (0, 36, 54, 6.529, 1.022),
    (0, 36, 55, 6.744, 1.014),
    (0, 36, 56, 6.49, 1.017),
    (0, 36, 58, 6.61, 1.019),
    (0, 36, 81, 6.98, 1.016),
    (0, 37, 37, 5.573, 1.374),
    (0, 37, 38, 5.737, 1.333),
    (0, 37, 39, 5.978, 1.375),
    (0, 37, 40, 6.168, 1.398),
    (0, 37, 41, 4.537, 1.468),
    (0, 37, 43, 4.764, 1.428),
    (0, 37, 45, 4.705, 1.431),
    (0, 37, 46, 6.191, 1.367),
    (0, 37, 55, 6.615, 1.352),
    (0, 37, 56, 5.055, 1.414),
    (0, 37, 58, 7.432, 1.326),
    (0, 37, 61, 5.724, 1.385),
    (0, 37, 62, 7.137, 1.335),
    (0, 37, 63, 6.095, 1.372),
    (0, 37, 64, 6.161, 1.379),
    (0, 37, 69, 5.396, 1.352),
    (0, 37, 78, 6.719, 1.375),
    (0, 37, 81, 3.987, 1.471),
    (0, 38, 38, 5.002, 1.246),
    (0, 38, 63, 7.299, 1.330),
    (0, 38, 64, 6.978, 1.340),
    (0, 38, 69, 5.036, 1.321),
    (0, 38, 78, 6.218, 1.366),
    (0, 39, 40, 4.101, 1.367),
    (0, 39, 45, 3.524, 1.402),
    (0, 39, 63, 6.301, 1.364),
    (0, 39, 64, 6.357, 1.361),
    (0, 39, 65, 5.513, 1.339),
    (0, 39, 78, 6.137, 1.369),
    (0, 40, 40, 4.248, 1.359),
    (0, 40, 45, 4.305, 1.356),
    (0, 40, 46, 4.727, 1.335),
    (0, 40, 54, 6.817, 1.256),
    (0, 40, 63, 6.733, 1.348),
    (0, 40, 64, 6.644, 1.351),
    (0, 40, 78, 5.9, 1.378),
    (0, 41, 41, 5.029, 1.443),
    (0, 41, 55, 5.577, 1.391),
    (0, 41, 62, 7.137, 1.335),
    (0, 41, 72, 4.519, 1.678),
    (0, 41, 80, 5.222, 1.434),
    (0, 42, 61, 16.223, 1.087),
    (0, 43, 43, 4.211, 1.361),
    (0, 43, 45, 3.71, 1.390),
    (0, 43, 64, 5.389, 1.399),
    (0, 44, 63, 3.589, 1.717),
    (0, 44, 65, 3.374, 1.684),
    (0, 44, 78, 3.711, 1.734),
    (0, 44, 80, 3.91, 1.719),
    (0, 45, 63, 5.119, 1.411),
    (0, 45, 64, 5.076, 1.413),
    (0, 45, 78, 5.724, 1.385),
    (0, 47, 53, 12.192, 1.140),
    (0, 49, 50, 6.812, 0.991),
    (0, 51, 52, 7.1, 0.987),
    (0, 55, 57, 7.227, 1.319),
    (0, 55, 62, 3.977, 1.374),
    (0, 55, 64, 5.529, 1.393),
    (0, 55, 80, 7.5, 1.324),
    (0, 56, 57, 4.137, 1.383),
    (0, 56, 63, 5.9, 1.378),
    (0, 56, 80, 6.47, 1.357),
    (0, 58, 63, 6.794, 1.346),
    (0, 58, 64, 6.164, 1.368),
    (0, 59, 63, 5.787, 1.360),
    (0, 59, 65, 4.756, 1.388),
    (0, 59, 78, 6.127, 1.364),
    (0, 59, 80, 7.064, 1.332),
    (0, 59, 82, 3.855, 1.431),
    (0, 60, 61, 15.749, 1.170),
    (0, 62, 63, 6.947, 1.341),
    (0, 62, 64, 6.273, 1.364),
    (0, 63, 64, 7.118, 1.377),
    (0, 63, 66, 8.326, 1.313),
    (0, 63, 72, 4.503, 1.679),
    (0, 63, 78, 7.434, 1.352),
    (0, 63, 81, 7.778, 1.316),
    (0, 64, 64, 4.313, 1.418),
    (0, 64, 65, 8.258, 1.335),
    (0, 64, 66, 4.456, 1.369),
    (0, 64, 78, 5.492, 1.422),
    (0, 64, 81, 5.824, 1.381),
    (0, 64, 82, 6.794, 1.346),
    (0, 65, 66, 7.243, 1.323),
    (0, 65, 78, 8.447, 1.298),
    (0, 65, 81, 5.223, 1.313),
    (0, 65, 82, 5.622, 1.297),
    (0, 66, 66, 3.874, 1.368),
    (0, 66, 78, 6.385, 1.360),
    (0, 66, 81, 3.96, 1.375),
    (0, 67, 67, 6.085, 1.280),
    (0, 71, 75, 2.852, 1.423),
    (0, 72, 73, 2.628, 2.035),
    (0, 76, 76, 4.286, 1.357),
    (0, 76, 78, 6.824, 1.345),
    (0, 78, 78, 5.573, 1.374),
    (0, 78, 79, 8.89, 1.287),
    (0, 78, 81, 5.046, 1.381),
    (0, 79, 79, 6.408, 1.269),
    (0, 79, 81, 4.305, 1.356),
    (0, 80, 81, 8.237, 1.335),
    (1, 2, 2, 5.31, 1.430),
    (1, 2, 3, 4.565, 1.468),
    (1, 2, 4, 5.657, 1.415),
    (1, 2, 9, 6.385, 1.360),
    (1, 2, 37, 5.007, 1.449),
    (1, 2, 39, 6.164, 1.368),
    (1, 2, 63, 6.03, 1.400),
    (1, 2, 64, 5.754, 1.411),
    (1, 2, 67, 4.685, 1.432),
    (1, 2, 81, 6.357, 1.361),
    (1, 3, 3, 4.418, 1.489),
    (1, 3, 4, 5.135, 1.438),
    (1, 3, 9, 6.273, 1.364),
    (1, 3, 30, 4.481, 1.471),
    (1, 3, 37, 4.488, 1.457),
    (1, 3, 39, 5.978, 1.375),
    (1, 3, 54, 2.771, 1.563),
    (1, 3, 57, 5.492, 1.422),
    (1, 3, 58, 5.163, 1.409),
    (1, 3, 63, 5.468, 1.423),
    (1, 3, 64, 5.288, 1.431),
    (1, 3, 78, 5.705, 1.413),
    (1, 3, 80, 6.719, 1.375),
    (1, 4, 9, 7.041, 1.338),
    (1, 4, 37, 5.445, 1.424),
    (1, 4, 63, 5.633, 1.416),
    (1, 4, 64, 5.492, 1.422),
    (1, 9, 9, 3.808, 1.384),
    (1, 9, 37, 5.529, 1.393),
    (1, 9, 39, 4.685, 1.337),
    (1, 9, 57, 6.824, 1.345),
    (1, 9, 63, 6.824, 1.345),
    (1, 9, 64, 5.458, 1.396),
    (1, 9, 78, 6.644, 1.351),
    (1, 9, 81, 3.909, 1.378),
    (1, 30, 30, 5.355, 1.428),
    (1, 30, 67, 5.274, 1.404),
    (1, 37, 37, 5.178, 1.436),
    (1, 37, 39, 5.65, 1.388),
    (1, 37, 57, 5.092, 1.440),
    (1, 37, 58, 5.055, 1.414),
    (1, 37, 63, 5.178, 1.436),
    (1, 37, 64, 5.265, 1.432),
    (1, 37, 67, 4.725, 1.430),
    (1, 37, 81, 4.531, 1.440),
    (1, 39, 63, 6.137, 1.369),
    (1, 39, 64, 5.482, 1.395),
    (1, 57, 63, 5.4, 1.426),
    (1, 57, 64, 5.135, 1.438),
    (1, 63, 63, 5.729, 1.412),
    (1, 64, 64, 4.926, 1.448),
];

/// Look up the exact MMFF94 bond parameter row by (bond class, t1, t2).
pub fn lookup_mmff94_bond(bond_class: u8, t1: u8, t2: u8) -> Option<(f64, f64)> {
    let key = (bond_class, t1.min(t2), t1.max(t2));
    MMFF94_BOND_TABLE
        .binary_search_by_key(&key, |r| (r.0, r.1, r.2))
        .ok()
        .map(|i| (MMFF94_BOND_TABLE[i].3, MMFF94_BOND_TABLE[i].4))
}

#[cfg(test)]
mod tests {
    #[allow(unused_imports)]
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

    // Regression: P_4-C_AR / P_3-C_AR bond params must come from the RDKit table
    // (3.586/1.755 and 3.207/1.788), not the estimation fallback. Without these,
    // Ph3P=O diverged +9.4 kcal/mol (bond term overstated). RDKit-verified.
    #[test]
    fn test_p_car_bond_params() {
        let p4_car =
            get_bond_params(MMFFAtomType::P_4, MMFFAtomType::C_AR, BondType::Single).unwrap();
        assert!(
            (p4_car.k_bond - 3.586).abs() < 1e-3,
            "P_4-C_AR kb={}",
            p4_car.k_bond
        );
        assert!(
            (p4_car.r0 - 1.755).abs() < 1e-3,
            "P_4-C_AR r0={}",
            p4_car.r0
        );
        let p3_car =
            get_bond_params(MMFFAtomType::P_3, MMFFAtomType::C_AR, BondType::Single).unwrap();
        assert!(
            (p3_car.k_bond - 3.207).abs() < 1e-3,
            "P_3-C_AR kb={}",
            p3_car.k_bond
        );
        assert!(
            (p3_car.r0 - 1.788).abs() < 1e-3,
            "P_3-C_AR r0={}",
            p3_car.r0
        );
        // symmetric
        assert!(get_bond_params(MMFFAtomType::C_AR, MMFFAtomType::P_4, BondType::Single).is_some());
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
