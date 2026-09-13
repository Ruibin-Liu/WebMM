use super::MMFFAtomType;
use crate::molecule::BondType;

pub fn mmff_type_id(t: MMFFAtomType) -> u8 {
    match t {
        MMFFAtomType::H => 5,
        MMFFAtomType::H_OH => 31,
        MMFFAtomType::H_ONC => 21,
        MMFFAtomType::H_COOH => 24,
        MMFFAtomType::H_OAR => 29,
        MMFFAtomType::H_N3 => 23,
        MMFFAtomType::H_NAM => 28,
        MMFFAtomType::H_NIM => 27,
        MMFFAtomType::HS => 71,
        MMFFAtomType::CR4R => 20,
        MMFFAtomType::CE4R => 30,
        MMFFAtomType::CR3R => 22,
        MMFFAtomType::HNRP => 36,
        MMFFAtomType::S2CM => 72,
        MMFFAtomType::HOS => 33,
        MMFFAtomType::H_OXP => 50,
        MMFFAtomType::H_OXP2 => 52,
        MMFFAtomType::O_3P => 49,
        MMFFAtomType::O_2P => 51,
        MMFFAtomType::N_RAD => 62,
        MMFFAtomType::NPYL_M => 76,
        MMFFAtomType::N_PYR => 58,
        MMFFAtomType::N_T3 => 68,
        MMFFAtomType::N_POX2 => 69,
        MMFFAtomType::N_SO => 48,
        MMFFAtomType::N_IM => 54,
        MMFFAtomType::N_GD => 56,
        MMFFAtomType::N_5OX => 67,
        MMFFAtomType::N_5POS => 81,
        MMFFAtomType::N_5OX2 => 82,
        MMFFAtomType::CL4 => 77,
        MMFFAtomType::S_O3 => 73,
        MMFFAtomType::S_CSO => 74,
        MMFFAtomType::P_ARM => 75,
        MMFFAtomType::C5A_M => 78,
        MMFFAtomType::C_IM => 80,
        MMFFAtomType::C_3 => 1,
        MMFFAtomType::C_2 => 3,
        MMFFAtomType::C_VIN => 2,
        MMFFAtomType::C_CO2 => 41,
        MMFFAtomType::C_1 => 4,
        MMFFAtomType::C_AR => 37,
        MMFFAtomType::C5A => 63,
        MMFFAtomType::C5B => 64,
        MMFFAtomType::C_CAT => 56,
        MMFFAtomType::C_AN => 57,
        MMFFAtomType::CID => 60,
        MMFFAtomType::NID => 61,
        MMFFAtomType::NCN_PLUS => 55,
        MMFFAtomType::OXIDE => 35,
        MMFFAtomType::N_3 => 8,
        MMFFAtomType::N_2 => 9,
        MMFFAtomType::N_1 => 42,
        MMFFAtomType::N_AR => 38,
        MMFFAtomType::NPYL => 39,
        MMFFAtomType::N_PL3 => 40,
        MMFFAtomType::N_AM => 10,
        MMFFAtomType::N_4 => 34,
        MMFFAtomType::N_2Z => 53,
        MMFFAtomType::N_1M => 47,
        MMFFAtomType::N_SOM => 48,
        MMFFAtomType::N_NO2 => 45,
        MMFFAtomType::N_NITROSO => 46,
        MMFFAtomType::N_SO2 => 43,
        MMFFAtomType::N5A => 65,
        MMFFAtomType::N5B => 66,
        MMFFAtomType::N5 => 79,
        MMFFAtomType::N_POX => 69,
        MMFFAtomType::O_3 => 6,
        MMFFAtomType::O_2 => 7,
        MMFFAtomType::O_R => 6,
        MMFFAtomType::OH2 => 70,
        MMFFAtomType::OFUR => 59,
        MMFFAtomType::O_CO2 => 32,
        MMFFAtomType::O_3_Z => 35,
        MMFFAtomType::F => 11,
        MMFFAtomType::Cl => 12,
        MMFFAtomType::Br => 13,
        MMFFAtomType::I => 14,
        MMFFAtomType::F_M => 89,
        MMFFAtomType::CL_M => 90,
        MMFFAtomType::BR_M => 91,
        MMFFAtomType::S_3 => 15,
        MMFFAtomType::S_2 => 16,
        MMFFAtomType::S_AR => 44,
        MMFFAtomType::S_OX => 17,
        MMFFAtomType::S_O2 => 18,
        MMFFAtomType::P_3 => 26,
        MMFFAtomType::P_4 => 25,
        MMFFAtomType::Si => 19,
        MMFFAtomType::P_3D => 26,
        MMFFAtomType::S_3D => 15,
        MMFFAtomType::S_3D2 => 15,
        MMFFAtomType::Fe_P2 => 87,
        MMFFAtomType::Fe_P3 => 88,
        MMFFAtomType::Li => 92,
        MMFFAtomType::Na => 93,
        MMFFAtomType::K => 94,
        MMFFAtomType::Zn_P2 => 95,
        MMFFAtomType::Ca_P2 => 96,
        MMFFAtomType::Cu_P1 => 97,
        MMFFAtomType::Cu_P2 => 98,
        MMFFAtomType::Mg_P2 => 99,
    }
}

/// Equivalence levels for MMFF atom type lookup protocol.
/// Index by type_id. Returns [level1, level2, level3, level4].
/// Level 0 = wildcard, levels get progressively more specific.
pub fn get_eq_levels(type_id: u8) -> [u8; 4] {
    // The single source of truth is the RDKit MMFFDef table
    // (MMFF_DEF_EQ_LEVELS). The separate table previously kept here had its
    // columns shifted (it was missing the level-5 column), which silently
    // changed torsion wildcard lookups.
    let idx = type_id as usize;
    if idx < 100 {
        crate::mmff::mmff_tables::MMFF_DEF_EQ_LEVELS[idx]
    } else {
        [0, 0, 0, 0]
    }
}

#[rustfmt::skip]
const MMFFPROP_SBMB: [bool; 100] = [
    false, false, true, true, true, // 0-4
    false, false, false, false, true, // 5-9
    false, false, false, false, false, // 10-14
    false, false, false, false, false, // 15-19
    false, false, false, false, false, // 20-24
    false, false, false, false, false, // 25-29
    true, false, false, false, false, // 30-34
    false, false, true, false, true, // 35-39
    false, false, false, false, false, // 40-44
    false, false, false, false, false, // 45-49
    false, false, false, false, true, // 50-54
    false, false, true, true, false, // 55-59
    false, false, false, true, true, // 60-64
    false, false, true, false, false, // 65-69
    false, false, false, false, false, // 70-74
    true, false, false, true, false, // 75-79
    true, true, false, false, false, // 80-84
    false, false, false, false, false, // 85-89
    false, false, false, false, false, // 90-94
    false, false, false, false, false, // 95-99
];

#[rustfmt::skip]
const MMFFPROP_AROM: [bool; 100] = [
    false, false, false, false, false, // 0-4
    false, false, false, false, false, // 5-9
    false, false, false, false, false, // 10-14
    false, false, false, false, false, // 15-19
    false, false, false, false, false, // 20-24
    false, false, false, false, false, // 25-29
    false, false, false, false, false, // 30-34
    false, false, true, true, true, // 35-39
    false, false, false, false, true, // 40-44
    false, false, false, false, false, // 45-49
    false, false, false, false, false, // 50-54
    false, false, false, true, true, // 55-59
    false, false, false, true, true, // 60-64
    true, true, false, false, true, // 65-69
    false, false, false, false, false, // 70-74
    false, false, false, true, true, // 75-79
    false, true, true, false, false, // 80-84
    false, false, false, false, false, // 85-89
    false, false, false, false, false, // 90-94
    false, false, false, false, false, // 95-99
];

pub fn is_sbmb(type_id: u8) -> bool {
    let idx = type_id as usize;
    idx < 100 && MMFFPROP_SBMB[idx]
}

pub fn is_arom(type_id: u8) -> bool {
    let idx = type_id as usize;
    idx < 100 && MMFFPROP_AROM[idx]
}

/// RDKit's getMMFFBondType: returns 1 if bond is SINGLE and both atoms
/// have sbmb=1 or both atoms have arom=1, else 0.
pub fn get_mmff_bond_type(bond_type: BondType, type_id_a: u8, type_id_b: u8) -> u8 {
    if bond_type == BondType::Single {
        let a_sbmb = is_sbmb(type_id_a);
        let b_sbmb = is_sbmb(type_id_b);
        let a_arom = is_arom(type_id_a);
        let b_arom = is_arom(type_id_b);
        if (a_sbmb && b_sbmb) || (a_arom && b_arom) {
            return 1;
        }
    }
    0
}

/// RDKit's getMMFFTorsionType: determines (primary, secondary) torsion type.
/// bond_ij_type, bond_jk_type, bond_kl_type: from get_mmff_bond_type (0 or 1)
/// bond_jk_actual: the actual BondType of the central J-K bond
/// ring4: true if atoms i,j,k,l all belong to the same 4-membered ring
/// ring5: true if atoms i,j,k,l all belong to the same 5-membered ring
/// type_i..type_l: MMFF type IDs of the four atoms
#[allow(clippy::too_many_arguments)]
pub fn get_mmff_torsion_type(
    bond_ij_type: u8,
    bond_jk_type: u8,
    bond_kl_type: u8,
    bond_jk_actual: BondType,
    ring4: bool,
    ring5: bool,
    type_i: u8,
    type_j: u8,
    type_k: u8,
    type_l: u8,
) -> (u8, u8) {
    let mut torsion_type = bond_jk_type;
    let mut second_torsion_type = 0u8;

    if bond_jk_type == 0
        && bond_jk_actual == BondType::Single
        && (bond_ij_type == 1 || bond_kl_type == 1)
    {
        torsion_type = 2;
    }

    if ring4 {
        second_torsion_type = torsion_type;
        torsion_type = 4;
    } else if ring5 && (type_i == 1 || type_j == 1 || type_k == 1 || type_l == 1) {
        second_torsion_type = torsion_type;
        torsion_type = 5;
    }

    (torsion_type, second_torsion_type)
}

/// Element (atomic number) of an MMFF atom type — used by the RDKit
/// empirical bond rule, which keys on elements rather than MMFF types.
pub fn element_of(t: MMFFAtomType) -> i32 {
    use MMFFAtomType::*;
    match t {
        // hydrogen
        H | H_OH | H_ONC | H_COOH | H_OAR | H_N3 | H_NAM | H_NIM | HNRP | HS => 1,
        // carbon
        C_3 | C_2 | C_VIN | C_CO2 | C_1 | C_AR | C5A | C5B | C5A_M | C_IM | C_CAT | C_AN | CID
        | NID | CR4R | CE4R | CR3R => 6,
        // nitrogen
        N_3 | N_2 | N_1 | N_AR | NPYL | N_PL3 | N_AM | N_4 | N_2Z | N_1M | N_SOM | N_NO2
        | N_SO2 | N_NITROSO | N5A | N5B | N5 | N_POX | N_RAD | NPYL_M | N_PYR | N_T3 | N_POX2
        | N_SO | N_IM | N_GD | N_5OX | N_5POS | N_5OX2 | NCN_PLUS | OXIDE => 7,
        // oxygen
        O_3 | O_2 | O_R | OH2 | OFUR | O_CO2 | O_3_Z | O_3P | O_2P => 8,
        // fluorine
        F_M => 9,
        // silicon, phosphorus
        P_ARM | P_3 | P_4 | P_3D => 15,
        // sulfur
        S_3 | S_2 | S_AR | S_OX | S_O2 | S_O3 | S_CSO | S_3D | S_3D2 => 16,
        // chlorine, bromine (iodine has no dedicated type here)
        CL4 | CL_M => 17,
        BR_M => 35,
        S2CM => 16,
        HOS => 1,
        // fluorine (neutral F is its own variant in some paths)
        MMFFAtomType::F => 9,
        // hydrogen on oxidized O / anything else exotic (metals, dummies):
        // no element mapping — callers fall back to the class-based rule
        _ => -1,
    }
}
