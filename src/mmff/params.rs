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

#[rustfmt::skip]
const EQ_LEVELS: [[u8; 4]; 100] = [
    [ 0,  0,  0,  0], // 0 unused
    [ 1,  1,  1,  1], // 1 CR
    [ 2,  2,  2,  1], // 2 C=C
    [ 3,  3,  3,  1], // 3 C=O
    [ 4,  4,  4,  1], // 4 CSP
    [ 5,  5,  5,  5], // 5 HC
    [ 6,  6,  6,  6], // 6 OR
    [ 7,  7,  7,  6], // 7 O=C
    [ 8,  8,  8,  8], // 8 NR
    [ 9,  9,  9,  8], // 9 N=C
    [10, 10, 10,  8], // 10 NC=O
    [11, 11, 11, 11], // 11 F
    [12, 12, 12, 12], // 12 CL
    [13, 13, 13, 13], // 13 BR
    [14, 14, 14, 14], // 14 I
    [15, 15, 15, 15], // 15 S
    [16, 16, 16, 15], // 16 S=C
    [17, 17, 17, 15], // 17 S=O
    [18, 18, 18, 15], // 18 SO2
    [19, 19, 19, 19], // 19 SI
    [20, 20,  1,  1], // 20 CR4R
    [21, 21, 21,  5], // 21 HOR
    [22, 22, 22,  1], // 22 CR3R
    [23, 23, 23,  5], // 23 HNR
    [24, 24, 24,  5], // 24 HOCO
    [25, 25, 25, 25], // 25 PO4
    [26, 26, 26, 25], // 26 P
    [27, 27, 28,  5], // 27 HN=C
    [28, 28, 28,  5], // 28 HNCO
    [29, 29, 29,  5], // 29 HOCC
    [30, 30,  2,  1], // 30 CE4R
    [31, 31, 31, 31], // 31 HOH
    [32, 32,  7,  6], // 32 O2CM
    [33, 33, 21,  5], // 33 HOS
    [34, 34,  8,  8], // 34 NR+
    [35, 35,  6,  6], // 35 OM
    [36, 36, 36,  5], // 36 HNR+
    [37, 37,  2,  1], // 37 CB
    [38, 38,  9,  8], // 38 NPYD
    [39, 39, 10,  8], // 39 NPYL
    [40, 40, 10,  8], // 40 NC=C
    [41, 41,  3,  1], // 41 CO2M
    [42, 42, 42,  8], // 42 NSP
    [43, 43, 10,  8], // 43 NSO2
    [44, 44, 16, 15], // 44 STHI
    [45, 45, 10,  8], // 45 NO2
    [46, 46,  9,  8], // 46 N=O
    [47, 47, 42,  8], // 47 NAZT
    [48, 48,  9,  8], // 48 NSO
    [49, 49,  6,  6], // 49 O+
    [50, 50, 21,  5], // 50 HO+
    [51, 51,  7,  6], // 51 O=+
    [52, 52, 21,  5], // 52 HO=+
    [53, 53, 42,  8], // 53 =N=
    [54, 54,  9,  8], // 54 N+=C
    [55, 55, 10,  8], // 55 NCN+
    [56, 56, 10,  8], // 56 NGD+
    [57, 57,  2,  1], // 57 CGD+
    [58, 58, 10,  8], // 58 NPD+
    [59, 59,  6,  6], // 59 OFUR
    [60, 60,  4,  1], // 60 C%
    [61, 61, 42,  8], // 61 NR%
    [62, 62, 10,  8], // 62 NM
    [63, 63,  2,  1], // 63 C5A
    [64, 64,  2,  1], // 64 C5B
    [65, 65,  9,  8], // 65 N5A
    [66, 66,  9,  8], // 66 N5B
    [67, 67,  9,  8], // 67 N2OX
    [68, 68,  8,  8], // 68 N3OX
    [69, 69,  9,  8], // 69 NPOX
    [70, 70, 70, 70], // 70 OH2
    [71, 71,  5,  5], // 71 HS
    [72, 72, 16, 15], // 72 S2CM
    [73, 73, 18, 15], // 73 SO2M
    [74, 74, 17, 15], // 74 =S=O
    [75, 75, 26, 25], // 75 -P=C
    [76, 76,  9,  8], // 76 N5M
    [77, 77, 12, 12], // 77 CLO4
    [78, 78,  2,  1], // 78 C5
    [79, 79,  9,  8], // 79 N5
    [80, 80,  2,  1], // 80 CIM+
    [81, 81, 10,  8], // 81 NIM+
    [82, 82,  9,  8], // 82 N5AX
    [ 0,  0,  0,  0], // 83 gap
    [ 0,  0,  0,  0], // 84 gap
    [ 0,  0,  0,  0], // 85 gap
    [ 0,  0,  0,  0], // 86 gap
    [87, 87, 87, 87], // 87 FE+2
    [88, 88, 88, 88], // 88 FE+3
    [89, 89, 89, 89], // 89 F-
    [90, 90, 90, 90], // 90 CL-
    [91, 91, 91, 91], // 91 BR-
    [92, 92, 92, 92], // 92 LI+
    [93, 93, 93, 93], // 93 NA+
    [94, 94, 94, 94], // 94 K+
    [95, 95, 95, 95], // 95 ZN+2
    [96, 96, 96, 96], // 96 CA+2
    [97, 97, 97, 97], // 97 CU+1
    [98, 98, 98, 98], // 98 CU+2
    [99, 99, 99, 99], // 99 MG+2
];

/// Equivalence levels for MMFF atom type lookup protocol.
/// Index by type_id. Returns [level1, level2, level3, level4].
/// Level 0 = wildcard, levels get progressively more specific.
pub fn get_eq_levels(type_id: u8) -> [u8; 4] {
    let idx = type_id as usize;
    if idx < 100 {
        EQ_LEVELS[idx]
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
