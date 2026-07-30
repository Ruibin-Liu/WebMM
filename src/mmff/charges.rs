#![allow(clippy::approx_constant)]

use super::mmff_type_id;
use super::MMFFAtomType;
use crate::molecule::{BondType, Molecule};

struct PbciEntry {
    pbci: f64,
    fcadj: f64,
}

fn get_pbci(type_id: u8) -> PbciEntry {
    match type_id {
        1 => PbciEntry {
            pbci: 0.000,
            fcadj: 0.000,
        },
        2 => PbciEntry {
            pbci: -0.135,
            fcadj: 0.000,
        },
        3 => PbciEntry {
            pbci: -0.095,
            fcadj: 0.000,
        },
        4 => PbciEntry {
            pbci: -0.200,
            fcadj: 0.000,
        },
        5 => PbciEntry {
            pbci: -0.023,
            fcadj: 0.000,
        },
        6 => PbciEntry {
            pbci: -0.243,
            fcadj: 0.000,
        },
        7 => PbciEntry {
            pbci: -0.687,
            fcadj: 0.000,
        },
        8 => PbciEntry {
            pbci: -0.253,
            fcadj: 0.000,
        },
        9 => PbciEntry {
            pbci: -0.306,
            fcadj: 0.000,
        },
        10 => PbciEntry {
            pbci: -0.244,
            fcadj: 0.000,
        },
        11 => PbciEntry {
            pbci: -0.317,
            fcadj: 0.000,
        },
        12 => PbciEntry {
            pbci: -0.304,
            fcadj: 0.000,
        },
        13 => PbciEntry {
            pbci: -0.238,
            fcadj: 0.000,
        },
        14 => PbciEntry {
            pbci: -0.208,
            fcadj: 0.000,
        },
        15 => PbciEntry {
            pbci: -0.236,
            fcadj: 0.000,
        },
        16 => PbciEntry {
            pbci: -0.475,
            fcadj: 0.000,
        },
        17 => PbciEntry {
            pbci: -0.191,
            fcadj: 0.000,
        },
        18 => PbciEntry {
            pbci: -0.118,
            fcadj: 0.000,
        },
        19 => PbciEntry {
            pbci: 0.094,
            fcadj: 0.000,
        },
        20 => PbciEntry {
            pbci: -0.019,
            fcadj: 0.000,
        },
        21 => PbciEntry {
            pbci: 0.157,
            fcadj: 0.000,
        },
        22 => PbciEntry {
            pbci: -0.095,
            fcadj: 0.000,
        },
        23 => PbciEntry {
            pbci: 0.193,
            fcadj: 0.000,
        },
        24 => PbciEntry {
            pbci: 0.257,
            fcadj: 0.000,
        },
        25 => PbciEntry {
            pbci: 0.012,
            fcadj: 0.000,
        },
        26 => PbciEntry {
            pbci: -0.142,
            fcadj: 0.000,
        },
        27 => PbciEntry {
            pbci: 0.094,
            fcadj: 0.000,
        },
        28 => PbciEntry {
            pbci: 0.058,
            fcadj: 0.000,
        },
        29 => PbciEntry {
            pbci: 0.207,
            fcadj: 0.000,
        },
        30 => PbciEntry {
            pbci: -0.166,
            fcadj: 0.000,
        },
        31 => PbciEntry {
            pbci: 0.161,
            fcadj: 0.000,
        },
        32 => PbciEntry {
            pbci: -0.732,
            fcadj: 0.500,
        },
        33 => PbciEntry {
            pbci: 0.257,
            fcadj: 0.000,
        },
        34 => PbciEntry {
            pbci: -0.491,
            fcadj: 0.000,
        },
        35 => PbciEntry {
            pbci: -0.456,
            fcadj: 0.500,
        },
        36 => PbciEntry {
            pbci: -0.031,
            fcadj: 0.000,
        },
        37 => PbciEntry {
            pbci: -0.127,
            fcadj: 0.000,
        },
        38 => PbciEntry {
            pbci: -0.437,
            fcadj: 0.000,
        },
        39 => PbciEntry {
            pbci: -0.104,
            fcadj: 0.000,
        },
        40 => PbciEntry {
            pbci: -0.264,
            fcadj: 0.000,
        },
        41 => PbciEntry {
            pbci: 0.052,
            fcadj: 0.000,
        },
        42 => PbciEntry {
            pbci: -0.757,
            fcadj: 0.000,
        },
        43 => PbciEntry {
            pbci: -0.326,
            fcadj: 0.000,
        },
        44 => PbciEntry {
            pbci: -0.237,
            fcadj: 0.000,
        },
        45 => PbciEntry {
            pbci: -0.260,
            fcadj: 0.000,
        },
        46 => PbciEntry {
            pbci: -0.429,
            fcadj: 0.000,
        },
        47 => PbciEntry {
            pbci: -0.418,
            fcadj: 0.000,
        },
        48 => PbciEntry {
            pbci: -0.525,
            fcadj: 0.000,
        },
        49 => PbciEntry {
            pbci: -0.283,
            fcadj: 0.000,
        },
        50 => PbciEntry {
            pbci: 0.284,
            fcadj: 0.000,
        },
        51 => PbciEntry {
            pbci: -1.046,
            fcadj: 0.000,
        },
        52 => PbciEntry {
            pbci: -0.546,
            fcadj: 0.000,
        },
        53 => PbciEntry {
            pbci: -0.048,
            fcadj: 0.000,
        },
        54 => PbciEntry {
            pbci: -0.424,
            fcadj: 0.000,
        },
        55 => PbciEntry {
            pbci: -0.476,
            fcadj: 0.000,
        },
        56 => PbciEntry {
            pbci: -0.438,
            fcadj: 0.000,
        },
        57 => PbciEntry {
            pbci: -0.105,
            fcadj: 0.000,
        },
        58 => PbciEntry {
            pbci: -0.488,
            fcadj: 0.000,
        },
        59 => PbciEntry {
            pbci: -0.635,
            fcadj: 0.000,
        },
        60 => PbciEntry {
            pbci: -0.105,
            fcadj: 0.000,
        },
        61 => PbciEntry {
            pbci: -0.265,
            fcadj: 0.000,
        },
        62 => PbciEntry {
            pbci: -0.125,
            fcadj: 0.250,
        },
        63 => PbciEntry {
            pbci: -0.180,
            fcadj: 0.000,
        },
        64 => PbciEntry {
            pbci: -0.181,
            fcadj: 0.000,
        },
        65 => PbciEntry {
            pbci: -0.475,
            fcadj: 0.000,
        },
        66 => PbciEntry {
            pbci: -0.467,
            fcadj: 0.000,
        },
        67 => PbciEntry {
            pbci: -0.099,
            fcadj: 0.000,
        },
        68 => PbciEntry {
            pbci: -0.135,
            fcadj: 0.000,
        },
        69 => PbciEntry {
            pbci: -0.099,
            fcadj: 0.000,
        },
        70 => PbciEntry {
            pbci: -0.269,
            fcadj: 0.000,
        },
        71 => PbciEntry {
            pbci: -0.071,
            fcadj: 0.000,
        },
        72 => PbciEntry {
            pbci: -0.580,
            fcadj: 0.500,
        },
        73 => PbciEntry {
            pbci: -0.200,
            fcadj: 0.000,
        },
        74 => PbciEntry {
            pbci: -0.301,
            fcadj: 0.000,
        },
        75 => PbciEntry {
            pbci: -0.255,
            fcadj: 0.000,
        },
        76 => PbciEntry {
            pbci: -0.568,
            fcadj: 0.250,
        },
        77 => PbciEntry {
            pbci: -0.282,
            fcadj: 0.000,
        },
        78 => PbciEntry {
            pbci: -0.168,
            fcadj: 0.000,
        },
        79 => PbciEntry {
            pbci: -0.471,
            fcadj: 0.000,
        },
        80 => PbciEntry {
            pbci: -0.144,
            fcadj: 0.000,
        },
        81 => PbciEntry {
            pbci: -0.514,
            fcadj: 0.000,
        },
        82 => PbciEntry {
            pbci: -0.099,
            fcadj: 0.000,
        },
        87 => PbciEntry {
            pbci: 2.000,
            fcadj: 0.000,
        },
        88 => PbciEntry {
            pbci: 3.000,
            fcadj: 0.000,
        },
        89 => PbciEntry {
            pbci: -1.000,
            fcadj: 0.000,
        },
        90 => PbciEntry {
            pbci: -1.000,
            fcadj: 0.000,
        },
        91 => PbciEntry {
            pbci: -1.000,
            fcadj: 0.000,
        },
        92 => PbciEntry {
            pbci: 1.000,
            fcadj: 0.000,
        },
        93 => PbciEntry {
            pbci: 1.000,
            fcadj: 0.000,
        },
        94 => PbciEntry {
            pbci: 1.000,
            fcadj: 0.000,
        },
        95 => PbciEntry {
            pbci: 2.000,
            fcadj: 0.000,
        },
        96 => PbciEntry {
            pbci: 2.000,
            fcadj: 0.000,
        },
        97 => PbciEntry {
            pbci: 1.000,
            fcadj: 0.000,
        },
        98 => PbciEntry {
            pbci: 2.000,
            fcadj: 0.000,
        },
        99 => PbciEntry {
            pbci: 2.000,
            fcadj: 0.000,
        },
        _ => PbciEntry {
            pbci: 0.0,
            fcadj: 0.0,
        },
    }
}

/// MMFF bond type for BCI charge lookup (RDKit getMMFFBondType):
/// - bt=1 if SINGLE bond and both atoms are sbmb or both are aromatic
/// - bt=0 for everything else (including Double, Triple, Aromatic bonds)
fn mmff_bond_type(bt: BondType, type_i: u8, type_j: u8) -> u8 {
    match bt {
        BondType::Single => {
            let both_sbmb =
                crate::mmff::params::is_sbmb(type_i) && crate::mmff::params::is_sbmb(type_j);
            let both_arom = crate::mmff::params::is_arom(type_i) && crate::mmff::params::is_arom(type_j);
            if both_sbmb || both_arom {
                1
            } else {
                0
            }
        }
        _ => 0,
    }
}

fn lookup_bci(bt: u8, i_type: u8, j_type: u8) -> Option<f64> {
    let (can_i, can_j, sign) = if i_type <= j_type {
        (i_type, j_type, -1.0f64)
    } else {
        (j_type, i_type, 1.0f64)
    };

    let bci = lookup_bci_canonical(bt, can_i, can_j)?;
    Some(sign * bci)
}

#[rustfmt::skip]
const BCI_TABLE: &[(u8, u8, u8, f64)] = &[
    (0, 1, 1, 0.0000),
    (0, 1, 2, -0.1382),
    (0, 1, 3, -0.0610),
    (0, 1, 4, -0.2000),
    (0, 1, 5, 0.0000),
    (0, 1, 6, -0.2800),
    (0, 1, 8, -0.2700),
    (0, 1, 9, -0.2460),
    (0, 1, 10, -0.3001),
    (0, 1, 11, -0.3400),
    (0, 1, 12, -0.2900),
    (0, 1, 13, -0.2300),
    (0, 1, 14, -0.1900),
    (0, 1, 15, -0.2300),
    (0, 1, 17, -0.1935),
    (0, 1, 18, -0.1052),
    (0, 1, 19, 0.0805),
    (0, 1, 20, 0.0000),
    (0, 1, 22, -0.0950),
    (0, 1, 25, 0.0000),
    (0, 1, 26, -0.1669),
    (0, 1, 34, -0.5030),
    (0, 1, 35, -0.4274),
    (0, 1, 37, -0.1435),
    (0, 1, 39, -0.2556),
    (0, 1, 40, -0.3691),
    (0, 1, 41, 0.1060),
    (0, 1, 43, -0.3557),
    (0, 1, 45, -0.2402),
    (0, 1, 46, -0.3332),
    (0, 1, 54, -0.3461),
    (0, 1, 55, -0.4895),
    (0, 1, 56, -0.3276),
    (0, 1, 57, -0.1050),
    (0, 1, 58, -0.4880),
    (0, 1, 61, -0.2657),
    (0, 1, 62, -0.2000),
    (0, 1, 63, -0.1800),
    (0, 1, 64, -0.1810),
    (0, 1, 67, -0.0990),
    (0, 1, 68, -0.2560),
    (0, 1, 72, -0.5500),
    (0, 1, 73, -0.0877),
    (0, 1, 75, -0.2550),
    (0, 1, 78, -0.1680),
    (0, 1, 80, -0.1440),
    (0, 1, 81, -0.5140),
    (0, 2, 2, 0.0000),
    (0, 2, 4, -0.0650),
    (0, 2, 5, 0.1500),
    (0, 2, 6, -0.0767),
    (0, 2, 10, -0.1090),
    (0, 2, 11, -0.1495),
    (0, 2, 12, -0.1400),
    (0, 2, 13, -0.1100),
    (0, 2, 14, -0.0900),
    (0, 2, 15, -0.1010),
    (0, 2, 17, -0.0560),
    (0, 2, 18, 0.0170),
    (0, 2, 19, 0.2290),
    (0, 2, 20, 0.1160),
    (0, 2, 22, 0.0400),
    (0, 2, 25, 0.1470),
    (0, 2, 30, -0.0310),
    (0, 2, 34, -0.3560),
    (0, 2, 35, -0.3500),
    (0, 2, 40, -0.1000),
    (0, 2, 41, 0.2500),
    (0, 2, 43, -0.1910),
    (0, 2, 45, -0.2044),
    (0, 2, 46, -0.2940),
    (0, 2, 55, -0.3410),
    (0, 2, 56, -0.3030),
    (0, 2, 62, -0.0500),
    (0, 2, 72, -0.4500),
    (0, 3, 5, 0.0600),
    (0, 3, 6, -0.1500),
    (0, 3, 7, -0.5700),
    (0, 3, 9, -0.4500),
    (0, 3, 10, -0.0600),
    (0, 3, 11, -0.2220),
    (0, 3, 12, -0.2090),
    (0, 3, 15, -0.1410),
    (0, 3, 16, -0.3800),
    (0, 3, 17, -0.0960),
    (0, 3, 18, -0.0230),
    (0, 3, 20, 0.0530),
    (0, 3, 22, 0.0000),
    (0, 3, 25, 0.1070),
    (0, 3, 35, -0.3610),
    (0, 3, 40, -0.0500),
    (0, 3, 41, 0.1470),
    (0, 3, 43, -0.2363),
    (0, 3, 45, -0.1650),
    (0, 3, 48, -0.4300),
    (0, 3, 51, -0.9500),
    (0, 3, 53, -0.0134),
    (0, 3, 54, -0.4000),
    (0, 3, 55, -0.3810),
    (0, 3, 56, -0.3430),
    (0, 3, 62, -0.0300),
    (0, 3, 67, -0.0040),
    (0, 3, 74, -0.3190),
    (0, 3, 75, -0.2474),
    (0, 4, 5, 0.1770),
    (0, 4, 6, -0.0430),
    (0, 4, 7, -0.4870),
    (0, 4, 9, -0.3000),
    (0, 4, 10, -0.0440),
    (0, 4, 15, -0.0360),
    (0, 4, 20, 0.1810),
    (0, 4, 22, 0.1050),
    (0, 4, 30, 0.0340),
    (0, 4, 40, -0.0640),
    (0, 4, 42, -0.5571),
    (0, 4, 43, -0.1260),
    (0, 5, 19, 0.2000),
    (0, 5, 20, 0.0000),
    (0, 5, 22, -0.1000),
    (0, 5, 30, -0.1500),
    (0, 5, 37, -0.1500),
    (0, 5, 41, 0.2203),
    (0, 5, 57, -0.1500),
    (0, 5, 63, -0.1500),
    (0, 5, 64, -0.1500),
    (0, 5, 78, -0.1500),
    (0, 5, 80, -0.1500),
    (0, 6, 6, 0.0000),
    (0, 6, 8, -0.1000),
    (0, 6, 9, -0.0630),
    (0, 6, 10, 0.0355),
    (0, 6, 15, 0.0070),
    (0, 6, 17, 0.0520),
    (0, 6, 18, 0.1837),
    (0, 6, 19, 0.2974),
    (0, 6, 20, 0.2579),
    (0, 6, 21, 0.4000),
    (0, 6, 22, 0.1480),
    (0, 6, 24, 0.5000),
    (0, 6, 25, 0.2712),
    (0, 6, 26, 0.1010),
    (0, 6, 29, 0.4500),
    (0, 6, 30, 0.0770),
    (0, 6, 33, 0.5000),
    (0, 6, 37, 0.0825),
    (0, 6, 39, 0.1390),
    (0, 6, 40, -0.0210),
    (0, 6, 41, 0.2950),
    (0, 6, 43, -0.0830),
    (0, 6, 45, -0.0090),
    (0, 6, 54, -0.1810),
    (0, 6, 55, -0.2330),
    (0, 6, 57, 0.1380),
    (0, 6, 58, -0.2450),
    (0, 6, 63, 0.0630),
    (0, 6, 64, 0.0620),
    (0, 7, 17, 0.5000),
    (0, 7, 46, 0.1618),
    (0, 7, 74, 0.5000),
    (0, 8, 8, 0.0000),
    (0, 8, 9, -0.0530),
    (0, 8, 10, 0.0090),
    (0, 8, 12, -0.0510),
    (0, 8, 15, 0.0170),
    (0, 8, 17, 0.0620),
    (0, 8, 19, 0.3470),
    (0, 8, 20, 0.2096),
    (0, 8, 22, 0.1580),
    (0, 8, 23, 0.3600),
    (0, 8, 25, 0.2679),
    (0, 8, 26, 0.1110),
    (0, 8, 34, -0.2380),
    (0, 8, 39, 0.1490),
    (0, 8, 40, -0.0110),
    (0, 8, 43, -0.0730),
    (0, 8, 45, -0.0070),
    (0, 8, 46, -0.1760),
    (0, 8, 55, -0.2230),
    (0, 8, 56, -0.1850),
    (0, 9, 9, 0.0000),
    (0, 9, 10, 0.0620),
    (0, 9, 12, 0.0020),
    (0, 9, 15, 0.0700),
    (0, 9, 18, 0.1880),
    (0, 9, 19, 0.4000),
    (0, 9, 20, 0.2870),
    (0, 9, 25, 0.3180),
    (0, 9, 27, 0.4000),
    (0, 9, 34, -0.1850),
    (0, 9, 35, -0.1500),
    (0, 9, 40, 0.0420),
    (0, 9, 41, 0.3580),
    (0, 9, 45, 0.0460),
    (0, 9, 53, 0.3179),
    (0, 9, 54, -0.1180),
    (0, 9, 55, -0.1700),
    (0, 9, 56, -0.1320),
    (0, 9, 62, 0.1810),
    (0, 9, 67, 0.2070),
    (0, 10, 10, 0.0000),
    (0, 10, 13, 0.0060),
    (0, 10, 14, 0.0360),
    (0, 10, 15, 0.0080),
    (0, 10, 17, 0.0530),
    (0, 10, 20, 0.2250),
    (0, 10, 22, 0.1490),
    (0, 10, 25, 0.2560),
    (0, 10, 26, 0.1020),
    (0, 10, 28, 0.3700),
    (0, 10, 34, -0.2470),
    (0, 10, 35, -0.2120),
    (0, 10, 37, 0.1170),
    (0, 10, 39, 0.1400),
    (0, 10, 40, -0.0200),
    (0, 10, 41, 0.2960),
    (0, 10, 45, -0.0160),
    (0, 10, 63, 0.0640),
    (0, 10, 64, 0.0630),
    (0, 11, 20, 0.2980),
    (0, 11, 22, 0.2317),
    (0, 11, 25, 0.3290),
    (0, 11, 26, 0.1750),
    (0, 11, 37, 0.1900),
    (0, 11, 40, 0.0530),
    (0, 12, 15, 0.0680),
    (0, 12, 18, 0.1860),
    (0, 12, 19, 0.3701),
    (0, 12, 20, 0.2900),
    (0, 12, 22, 0.2273),
    (0, 12, 25, 0.3160),
    (0, 12, 26, 0.2112),
    (0, 12, 37, 0.1770),
    (0, 12, 40, 0.0400),
    (0, 12, 57, 0.1990),
    (0, 12, 63, 0.1240),
    (0, 12, 64, 0.1230),
    (0, 13, 20, 0.2190),
    (0, 13, 22, 0.1430),
    (0, 13, 37, 0.1110),
    (0, 13, 64, 0.0570),
    (0, 14, 20, 0.1890),
    (0, 14, 37, 0.0810),
    (0, 15, 15, 0.0000),
    (0, 15, 18, 0.1180),
    (0, 15, 19, 0.3300),
    (0, 15, 20, 0.2170),
    (0, 15, 22, 0.1410),
    (0, 15, 25, 0.2480),
    (0, 15, 26, 0.0940),
    (0, 15, 30, 0.0700),
    (0, 15, 37, 0.1015),
    (0, 15, 40, -0.0280),
    (0, 15, 43, -0.0900),
    (0, 15, 57, 0.1310),
    (0, 15, 63, 0.0560),
    (0, 15, 64, 0.0550),
    (0, 15, 71, 0.1800),
    (0, 16, 16, 0.0000),
    (0, 17, 17, 0.0000),
    (0, 17, 20, 0.1720),
    (0, 17, 22, 0.0960),
    (0, 17, 37, 0.0640),
    (0, 17, 43, -0.1350),
    (0, 18, 18, 0.0000),
    (0, 18, 20, 0.0990),
    (0, 18, 22, 0.0230),
    (0, 18, 32, -0.6500),
    (0, 18, 37, -0.0090),
    (0, 18, 39, 0.0140),
    (0, 18, 43, -0.1380),
    (0, 18, 48, -0.5895),
    (0, 18, 55, -0.3580),
    (0, 18, 58, -0.3700),
    (0, 18, 62, 0.2099),
    (0, 18, 63, -0.0620),
    (0, 18, 64, -0.0630),
    (0, 18, 80, -0.0260),
    (0, 19, 19, 0.0000),
    (0, 19, 20, -0.1130),
    (0, 19, 37, -0.2210),
    (0, 19, 40, -0.3580),
    (0, 19, 63, -0.2740),
    (0, 19, 75, -0.3490),
    (0, 20, 20, 0.0000),
    (0, 20, 22, -0.0760),
    (0, 20, 25, 0.0310),
    (0, 20, 26, -0.1230),
    (0, 20, 30, -0.1380),
    (0, 20, 34, -0.4720),
    (0, 20, 37, -0.1080),
    (0, 20, 40, -0.2450),
    (0, 20, 41, 0.0710),
    (0, 20, 43, -0.3070),
    (0, 20, 45, -0.2410),
    (0, 22, 22, 0.0000),
    (0, 22, 30, -0.0710),
    (0, 22, 34, -0.3960),
    (0, 22, 37, -0.0320),
    (0, 22, 40, -0.1690),
    (0, 22, 41, 0.1470),
    (0, 22, 43, -0.2310),
    (0, 22, 45, -0.1650),
    (0, 23, 39, -0.2700),
    (0, 23, 62, -0.4000),
    (0, 23, 67, -0.2920),
    (0, 23, 68, -0.3600),
    (0, 25, 25, 0.0000),
    (0, 25, 32, -0.7000),
    (0, 25, 37, -0.1390),
    (0, 25, 39, -0.1160),
    (0, 25, 40, -0.2760),
    (0, 25, 43, -0.3380),
    (0, 25, 57, -0.1170),
    (0, 25, 63, -0.1920),
    (0, 25, 71, -0.0362),
    (0, 25, 72, -0.6773),
    (0, 26, 26, 0.0000),
    (0, 26, 34, -0.3490),
    (0, 26, 37, 0.0150),
    (0, 26, 40, -0.1220),
    (0, 26, 71, 0.0960),
    (0, 28, 40, -0.4000),
    (0, 28, 43, -0.4200),
    (0, 28, 48, -0.4000),
    (0, 30, 30, 0.0000),
    (0, 30, 40, -0.0980),
    (0, 31, 70, -0.4300),
    (0, 32, 41, 0.6500),
    (0, 32, 45, 0.5200),
    (0, 32, 67, 0.6330),
    (0, 32, 68, 0.7500),
    (0, 32, 69, 0.7500),
    (0, 32, 73, 0.3500),
    (0, 32, 77, 0.4500),
    (0, 32, 82, 0.6330),
    (0, 34, 36, 0.4500),
    (0, 34, 37, 0.3640),
    (0, 34, 43, 0.1650),
    (0, 35, 37, 0.3290),
    (0, 35, 63, 0.2760),
    (0, 36, 54, -0.4000),
    (0, 36, 55, -0.4500),
    (0, 36, 56, -0.4500),
    (0, 36, 58, -0.4570),
    (0, 36, 81, -0.4500),
    (0, 37, 37, 0.0000),
    (0, 37, 38, -0.3100),
    (0, 37, 39, 0.0230),
    (0, 37, 40, -0.1000),
    (0, 37, 41, 0.1790),
    (0, 37, 43, -0.1990),
    (0, 37, 45, -0.1330),
    (0, 37, 46, -0.3020),
    (0, 37, 55, -0.3490),
    (0, 37, 56, -0.3110),
    (0, 37, 58, -0.3610),
    (0, 37, 61, -0.1380),
    (0, 37, 62, 0.0020),
    (0, 37, 63, 0.0000),
    (0, 37, 64, 0.0000),
    (0, 37, 69, -0.0895),
    (0, 37, 78, -0.0410),
    (0, 37, 81, -0.3870),
    (0, 38, 38, 0.0000),
    (0, 38, 63, 0.2570),
    (0, 38, 64, 0.2560),
    (0, 38, 69, 0.3380),
    (0, 38, 78, 0.2690),
    (0, 39, 40, -0.1600),
    (0, 39, 45, -0.1560),
    (0, 39, 63, -0.1516),
    (0, 39, 64, -0.0770),
    (0, 39, 65, -0.4180),
    (0, 39, 78, -0.0640),
    (0, 40, 40, 0.0000),
    (0, 40, 45, 0.0040),
    (0, 40, 46, -0.1650),
    (0, 40, 54, -0.1600),
    (0, 40, 63, 0.0840),
    (0, 40, 64, 0.0830),
    (0, 40, 78, 0.0960),
    (0, 41, 41, 0.0000),
    (0, 41, 55, -0.5280),
    (0, 41, 62, -0.1770),
    (0, 41, 72, -0.5000),
    (0, 41, 80, -0.1960),
    (0, 42, 61, 0.4920),
    (0, 43, 43, 0.0000),
    (0, 43, 45, 0.0660),
    (0, 43, 64, 0.1450),
    (0, 44, 63, 0.0400),
    (0, 44, 65, -0.2207),
    (0, 44, 78, 0.0690),
    (0, 44, 80, 0.0930),
    (0, 45, 63, 0.0800),
    (0, 45, 64, 0.0790),
    (0, 45, 78, 0.0920),
    (0, 47, 53, 0.3700),
    (0, 49, 50, 0.5673),
    (0, 51, 52, 0.5000),
    (0, 55, 57, 0.3544),
    (0, 55, 62, 0.3510),
    (0, 55, 64, 0.2950),
    (0, 55, 80, 0.3320),
    (0, 56, 57, 0.4000),
    (0, 56, 63, 0.2580),
    (0, 56, 80, 0.2700),
    (0, 58, 63, 0.3080),
    (0, 58, 64, 0.3070),
    (0, 59, 63, 0.1400),
    (0, 59, 65, -0.1209),
    (0, 59, 78, 0.1690),
    (0, 59, 80, 0.1930),
    (0, 59, 82, 0.2380),
    (0, 60, 61, 0.3700),
    (0, 62, 63, -0.0550),
    (0, 62, 64, -0.0560),
    (0, 63, 63, 0.0000),
    (0, 63, 64, 0.0000),
    (0, 63, 66, -0.3381),
    (0, 63, 72, -0.4000),
    (0, 63, 78, 0.0120),
    (0, 63, 81, -0.3340),
    (0, 64, 64, 0.0000),
    (0, 64, 65, -0.2888),
    (0, 64, 66, -0.2272),
    (0, 64, 78, 0.0130),
    (0, 64, 81, -0.3330),
    (0, 64, 82, 0.0820),
    (0, 65, 66, 0.0000),
    (0, 65, 78, 0.3070),
    (0, 65, 81, -0.0390),
    (0, 65, 82, 0.3760),
    (0, 66, 66, 0.0000),
    (0, 66, 78, 0.2990),
    (0, 66, 81, -0.0470),
    (0, 71, 75, -0.0958),
    (0, 72, 73, 0.4500),
    (0, 76, 76, 0.0000),
    (0, 76, 78, 0.4000),
    (0, 78, 78, 0.0000),
    (0, 78, 79, -0.3030),
    (0, 78, 81, -0.3500),
    (0, 79, 81, -0.0430),
    (0, 80, 81, -0.4000),
    (1, 2, 2, 0.0000),
    (1, 2, 3, -0.0144),
    (1, 2, 4, -0.0650),
    (1, 2, 9, -0.1710),
    (1, 2, 37, 0.0284),
    (1, 2, 39, 0.0310),
    (1, 2, 63, -0.0450),
    (1, 2, 64, -0.0460),
    (1, 2, 67, 0.0360),
    (1, 2, 81, -0.3790),
    (1, 3, 3, 0.0000),
    (1, 3, 4, -0.1050),
    (1, 3, 9, -0.2110),
    (1, 3, 30, -0.0710),
    (1, 3, 37, 0.0862),
    (1, 3, 39, -0.0090),
    (1, 3, 54, -0.3290),
    (1, 3, 57, -0.0100),
    (1, 3, 58, -0.3930),
    (1, 3, 63, -0.0850),
    (1, 3, 64, -0.0860),
    (1, 3, 78, -0.0730),
    (1, 3, 80, -0.0490),
    (1, 4, 9, -0.1060),
    (1, 4, 37, 0.0730),
    (1, 4, 63, 0.0200),
    (1, 4, 64, 0.0190),
    (1, 9, 37, 0.1790),
    (1, 9, 39, 0.2020),
    (1, 9, 57, 0.2010),
    (1, 9, 63, 0.1260),
    (1, 9, 64, 0.1250),
    (1, 9, 78, 0.1380),
    (1, 9, 81, -0.2080),
    (1, 30, 67, 0.0670),
    (1, 37, 37, 0.0000),
    (1, 37, 39, 0.0230),
    (1, 37, 57, 0.0220),
    (1, 37, 58, -0.3610),
    (1, 37, 63, -0.0530),
    (1, 37, 64, -0.0540),
    (1, 37, 67, 0.0280),
    (1, 37, 81, -0.3870),
    (1, 39, 39, 0.0000),
    (1, 39, 63, -0.0760),
    (1, 39, 64, -0.0770),
    (1, 57, 63, -0.0750),
    (1, 57, 64, -0.0760),
    (1, 63, 63, 0.0000),
    (1, 78, 78, 0.0000),
    (4, 36, 58, -0.4500),
    (4, 37, 58, -0.3500),
    (4, 57, 58, -0.4000),
];

fn lookup_bci_canonical(bt: u8, i: u8, j: u8) -> Option<f64> {
    BCI_TABLE
        .binary_search_by_key(&(bt, i, j), |&(b, ti, tj, _)| (b, ti, tj))
        .ok()
        .map(|idx| BCI_TABLE[idx].3)
}

/// Compute MMFF formal charges from molecular formal charges and atom types.
/// Follows RDKit's computeMMFFCharges formal-charge stage.
/// Most atom types simply carry their molecular formal charge; the exceptions
/// (types 32, 35, 62, 72, 76) share charge across equivalent siblings on the
/// same central atom.
fn compute_mmff_formal_charges(mol: &Molecule, type_ids: &[u8]) -> Vec<f64> {
    let n = mol.atoms.len();
    let mut charges: Vec<f64> = vec![0.0; n];

    // Helper: degree (neighbor count) of an atom
    let degree = |idx: usize| mol.adjacency[idx].len();
    // Helper: count terminal O/S atoms bonded to a central atom
    let count_term_os = |center: usize| {
        mol.adjacency[center]
            .iter()
            .filter(|&&s| {
                mol.adjacency[s].len() == 1
                    && (mol.atoms[s].atomic_number == 8 || mol.atoms[s].atomic_number == 16)
            })
            .count()
    };

    for i in 0..n {
        let ti = type_ids[i];
        charges[i] = match ti {
            // O2CM (32): carboxylate/sulfonate/phosphate/sulfate O — charge shared
            32 => {
                let mut fchg = 0.0;
                for &nbr in &mol.adjacency[i] {
                    let n_term_os = count_term_os(nbr);
                    if n_term_os > 0 {
                        let nbr_type = type_ids[nbr];
                        if mol.atoms[nbr].atomic_number == 6 {
                            fchg = if n_term_os == 1 {
                                -1.0
                            } else {
                                -((n_term_os - 1) as f64 / n_term_os as f64)
                            };
                        } else if nbr_type == 45 && n_term_os == 3 {
                            fchg = -1.0 / 3.0;
                        } else if nbr_type == 25 {
                            fchg = if n_term_os == 1 {
                                0.0
                            } else {
                                -((n_term_os - 1) as f64 / n_term_os as f64)
                            };
                        } else if nbr_type == 18 {
                            let n_sec_n = mol.adjacency[nbr]
                                .iter()
                                .filter(|&&s| {
                                    mol.atoms[s].atomic_number == 7 && mol.adjacency[s].len() == 2
                                })
                                .count();
                            let total = n_sec_n + n_term_os;
                            fchg = if total <= 2 {
                                0.0
                            } else {
                                -((total - 2) as f64 / n_term_os as f64)
                            };
                        } else if nbr_type == 77 {
                            fchg = -(1.0 / n_term_os as f64);
                        } else if nbr_type == 69 {
                            // N-oxide O — formal charge is 0 (neutralized)
                            fchg = 0.0;
                        }
                        break;
                    }
                }
                fchg
            }
            // OM/OM2 (35): oxide oxygen — formal -1
            35 => -1.0,
            // NM (62): anionic divalent nitrogen — formal -1
            62 => -1.0,
            // NPYL_M (76): pyrrolide anion N — formal -1
            76 => -1.0,
            // SM (72): anionic terminal sulfur
            72 => {
                let mut fchg = 0.0;
                for &nbr in &mol.adjacency[i] {
                    let n_term_os = count_term_os(nbr);
                    if mol.atoms[nbr].atomic_number == 6 && n_term_os > 0 {
                        fchg = if n_term_os == 1 {
                            -1.0
                        } else {
                            -((n_term_os - 1) as f64 / n_term_os as f64)
                        };
                        break;
                    }
                }
                fchg
            }
            // Cations
            34 | 49 | 51 | 54 | 58 | 92 | 93 | 94 | 97 => 1.0,
            87 | 95 | 96 | 98 | 99 => 2.0,
            88 => 3.0,
            // NPOX (69): pyridine N-oxide N — formal charge is 0 (neutralized)
            69 => 0.0,
            // Guanidinium CGD+ C (57): formal charge distributed to NCN+ N's
            57 => 0.0,
            // NCN+ N (55) and N_GD (56): shares +1 from CGD+ C equally among NCN+/N_GD N's
            55 | 56 => {
                let mut fc = 0.0;
                for &nbr in &mol.adjacency[i] {
                    if type_ids[nbr] == 57 {
                        let n_ncn = mol.adjacency[nbr]
                            .iter()
                            .filter(|&&m| type_ids[m] == 55 || type_ids[m] == 56)
                            .count()
                            .max(1);
                        fc += 1.0 / n_ncn as f64;
                    }
                }
                fc
            }
            // N_5POS (81): use SDF formal charge (imidazolium has charge delocalized)
            81 => mol.atoms[i].charge as f64,
            // Halide anions
            89..=91 => -1.0,
            // Everything else: formal charge is 0
            _ => 0.0,
        };
    }

    // suppress unused warning
    let _ = degree;
    charges
}

/// Compute MMFF partial charges using Halgren MMFF.V equation 15.
///
/// pChg = (1 - M*v) * q0 + v * sum(nbr_formal) + sum(bci)
///
/// where M=crd, v=fcadj, q0=MMFF formal charge (adjusted by nbr charges if v==0),
/// sum(nbr_formal) is the sum of neighbor formal charges, and sum(bci) is
/// the sum of bond charge increments.
pub fn calculate_bci_charges(mol: &Molecule, atom_types: &[MMFFAtomType]) -> Vec<f64> {
    let n = mol.atoms.len();
    let type_ids: Vec<u8> = atom_types.iter().map(|t| mmff_type_id(*t)).collect();
    let mmff_formal = compute_mmff_formal_charges(mol, &type_ids);
    let mut charges = vec![0.0; n];

    for i in 0..n {
        let ti = type_ids[i];

        // Get crd (M) from the MMFF property table (integer coordination number)
        let m_val = crate::mmff::mmff_tables::get_mmff_prop(ti)
            .map(|(_, crd, _, _, _, _, _, _)| crd as f64)
            .unwrap_or(0.0);
        let v = get_pbci(ti).fcadj;

        // q0 starts as the MMFF formal charge
        let mut q0 = mmff_formal[i];

        // If v == 0 (fcadj is 0), neighbors' negative formal charges influence q0
        if v.abs() < 1e-10 {
            for &nbr in &mol.adjacency[i] {
                let nbr_fc = mmff_formal[nbr];
                if nbr_fc < 0.0 {
                    let nbr_degree = mol.adjacency[nbr].len().max(1) as f64;
                    q0 += nbr_fc / (2.0 * nbr_degree);
                }
            }
        }

        // Special case: anionic divalent N (type 62) with positive neighbor
        if ti == 62 {
            for &nbr in &mol.adjacency[i] {
                let nbr_fc = mmff_formal[nbr];
                if nbr_fc > 0.0 {
                    q0 -= nbr_fc / 2.0;
                }
            }
        }

        // sum of neighbor formal charges
        let sum_formal: f64 = mol.adjacency[i].iter().map(|&nbr| mmff_formal[nbr]).sum();

        // sum of BCI increments from all bonds
        let mut sum_bci = 0.0;
        for &nbr in &mol.adjacency[i] {
            let tj = type_ids[nbr];
            // Find the bond between i and nbr
            if let Some(bond) = mol.bonds.iter().find(|b| {
                (b.atom1 == i && b.atom2 == nbr) || (b.atom1 == nbr && b.atom2 == i)
            }) {
                let bt = mmff_bond_type(bond.bond_type, ti, tj);
                let bci = match lookup_bci(bt, ti, tj) {
                    Some(v) => v,
                    None => {
                        let pi = get_pbci(ti);
                        let pj = get_pbci(tj);
                        pi.pbci - pj.pbci
                    }
                };
                sum_bci += bci;
            }
        }

        // Equation 15: pChg = (1 - M*v) * q0 + v * sumFormal + sumPartial
        charges[i] = (1.0 - m_val * v) * q0 + v * sum_formal + sum_bci;
    }

    charges
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::molecule::{Atom, Bond};

    fn make_molecule(atoms: Vec<Atom>, bonds: Vec<Bond>) -> Molecule {
        let n = atoms.len();
        let mut adjacency = vec![vec![]; n];
        for bond in &bonds {
            adjacency[bond.atom1].push(bond.atom2);
            adjacency[bond.atom2].push(bond.atom1);
        }
        Molecule {
            atoms,
            bonds,
            name: String::new(),
            adjacency,
        }
    }

    #[test]
    fn test_ethanol_partial_charges() {
        let atoms = vec![
            Atom {
                symbol: "C".into(),
                atomic_number: 6,
                mass: 12.0,
                charge: 0.0,
                position: [0.0; 3],
                index: 0,
                stereo_parity: 0,
            },
            Atom {
                symbol: "C".into(),
                atomic_number: 6,
                mass: 12.0,
                charge: 0.0,
                position: [0.0; 3],
                index: 1,
                stereo_parity: 0,
            },
            Atom {
                symbol: "O".into(),
                atomic_number: 8,
                mass: 16.0,
                charge: 0.0,
                position: [0.0; 3],
                index: 2,
                stereo_parity: 0,
            },
            Atom {
                symbol: "H".into(),
                atomic_number: 1,
                mass: 1.0,
                charge: 0.0,
                position: [0.0; 3],
                index: 3,
                stereo_parity: 0,
            },
            Atom {
                symbol: "H".into(),
                atomic_number: 1,
                mass: 1.0,
                charge: 0.0,
                position: [0.0; 3],
                index: 4,
                stereo_parity: 0,
            },
            Atom {
                symbol: "H".into(),
                atomic_number: 1,
                mass: 1.0,
                charge: 0.0,
                position: [0.0; 3],
                index: 5,
                stereo_parity: 0,
            },
            Atom {
                symbol: "H".into(),
                atomic_number: 1,
                mass: 1.0,
                charge: 0.0,
                position: [0.0; 3],
                index: 6,
                stereo_parity: 0,
            },
            Atom {
                symbol: "H".into(),
                atomic_number: 1,
                mass: 1.0,
                charge: 0.0,
                position: [0.0; 3],
                index: 7,
                stereo_parity: 0,
            },
            Atom {
                symbol: "H".into(),
                atomic_number: 1,
                mass: 1.0,
                charge: 0.0,
                position: [0.0; 3],
                index: 8,
                stereo_parity: 0,
            },
        ];
        let bonds = vec![
            Bond {
                atom1: 0,
                atom2: 1,
                bond_type: BondType::Single,
                ..Default::default()
            },
            Bond {
                atom1: 1,
                atom2: 2,
                bond_type: BondType::Single,
                ..Default::default()
            },
            Bond {
                atom1: 2,
                atom2: 3,
                bond_type: BondType::Single,
                ..Default::default()
            },
            Bond {
                atom1: 0,
                atom2: 4,
                bond_type: BondType::Single,
                ..Default::default()
            },
            Bond {
                atom1: 0,
                atom2: 5,
                bond_type: BondType::Single,
                ..Default::default()
            },
            Bond {
                atom1: 0,
                atom2: 6,
                bond_type: BondType::Single,
                ..Default::default()
            },
            Bond {
                atom1: 1,
                atom2: 7,
                bond_type: BondType::Single,
                ..Default::default()
            },
            Bond {
                atom1: 1,
                atom2: 8,
                bond_type: BondType::Single,
                ..Default::default()
            },
        ];
        let mol = make_molecule(atoms, bonds);
        let atom_types = vec![
            MMFFAtomType::C_3,
            MMFFAtomType::C_3,
            MMFFAtomType::O_3,
            MMFFAtomType::H_ONC,
            MMFFAtomType::H,
            MMFFAtomType::H,
            MMFFAtomType::H,
            MMFFAtomType::H,
            MMFFAtomType::H,
        ];
        let charges = calculate_bci_charges(&mol, &atom_types);

        assert!(
            charges[2] < 0.0,
            "oxygen should have negative partial charge, got {}",
            charges[2]
        );
        assert!(
            charges[3] > 0.0,
            "hydroxyl hydrogen should have positive partial charge, got {}",
            charges[3]
        );
    }

    #[test]
    fn test_water_charges_sum_to_zero() {
        let atoms = vec![
            Atom {
                symbol: "O".into(),
                atomic_number: 8,
                mass: 16.0,
                charge: 0.0,
                position: [0.0; 3],
                index: 0,
                stereo_parity: 0,
            },
            Atom {
                symbol: "H".into(),
                atomic_number: 1,
                mass: 1.0,
                charge: 0.0,
                position: [0.0; 3],
                index: 1,
                stereo_parity: 0,
            },
            Atom {
                symbol: "H".into(),
                atomic_number: 1,
                mass: 1.0,
                charge: 0.0,
                position: [0.0; 3],
                index: 2,
                stereo_parity: 0,
            },
        ];
        let bonds = vec![
            Bond {
                atom1: 0,
                atom2: 1,
                bond_type: BondType::Single,
                ..Default::default()
            },
            Bond {
                atom1: 0,
                atom2: 2,
                bond_type: BondType::Single,
                ..Default::default()
            },
        ];
        let mol = make_molecule(atoms, bonds);
        let atom_types = vec![MMFFAtomType::O_3, MMFFAtomType::H_OH, MMFFAtomType::H_OH];
        let charges = calculate_bci_charges(&mol, &atom_types);

        let total: f64 = charges.iter().sum();
        assert!(
            total.abs() < 1e-10,
            "water charges should sum to zero, got {}",
            total
        );
    }

    #[test]
    fn test_neutralization_with_formal_charge() {
        let atoms = vec![
            Atom {
                symbol: "O".into(),
                atomic_number: 8,
                mass: 16.0,
                charge: -1.0,
                position: [0.0; 3],
                index: 0,
                stereo_parity: 0,
            },
            Atom {
                symbol: "H".into(),
                atomic_number: 1,
                mass: 1.0,
                charge: 0.0,
                position: [0.0; 3],
                index: 1,
                stereo_parity: 0,
            },
        ];
        let bonds = vec![Bond {
            atom1: 0,
            atom2: 1,
            bond_type: BondType::Single,
            ..Default::default()
        }];
        let mol = make_molecule(atoms, bonds);
        let atom_types = vec![MMFFAtomType::O_3, MMFFAtomType::H];
        let charges = calculate_bci_charges(&mol, &atom_types);

        // With eq 15, charges are computed per-atom; verify they're finite
        assert!(
            charges.iter().all(|c| c.is_finite()),
            "hydroxide charges should be finite, got {:?}",
            charges
        );
    }

    #[test]
    fn test_ammonia_charges() {
        let atoms = vec![
            Atom {
                symbol: "N".into(),
                atomic_number: 7,
                mass: 14.0,
                charge: 0.0,
                position: [0.0; 3],
                index: 0,
                stereo_parity: 0,
            },
            Atom {
                symbol: "H".into(),
                atomic_number: 1,
                mass: 1.0,
                charge: 0.0,
                position: [1.0, 0.0, 0.0],
                index: 1,
                stereo_parity: 0,
            },
            Atom {
                symbol: "H".into(),
                atomic_number: 1,
                mass: 1.0,
                charge: 0.0,
                position: [-0.5, 0.866, 0.0],
                index: 2,
                stereo_parity: 0,
            },
            Atom {
                symbol: "H".into(),
                atomic_number: 1,
                mass: 1.0,
                charge: 0.0,
                position: [-0.5, -0.866, 0.0],
                index: 3,
                stereo_parity: 0,
            },
        ];
        let bonds = vec![
            Bond {
                atom1: 0,
                atom2: 1,
                bond_type: BondType::Single,
                ..Default::default()
            },
            Bond {
                atom1: 0,
                atom2: 2,
                bond_type: BondType::Single,
                ..Default::default()
            },
            Bond {
                atom1: 0,
                atom2: 3,
                bond_type: BondType::Single,
                ..Default::default()
            },
        ];
        let mol = make_molecule(atoms, bonds);
        let atom_types = vec![
            MMFFAtomType::N_3,
            MMFFAtomType::H_N3,
            MMFFAtomType::H_N3,
            MMFFAtomType::H_N3,
        ];
        let charges = calculate_bci_charges(&mol, &atom_types);

        let total: f64 = charges.iter().sum();
        assert!(
            total.abs() < 1e-10,
            "Ammonia charges should sum to zero, got {}",
            total
        );
    }

    #[test]
    fn test_single_atom_charge() {
        let atoms = vec![Atom {
            symbol: "C".into(),
            atomic_number: 6,
            mass: 12.0,
            charge: 0.0,
            position: [0.0; 3],
            index: 0,
            stereo_parity: 0,
        }];
        let bonds: Vec<Bond> = vec![];
        let mol = make_molecule(atoms, bonds);
        let atom_types = vec![MMFFAtomType::C_3];
        let charges = calculate_bci_charges(&mol, &atom_types);

        assert_eq!(charges.len(), 1);
        assert!(charges[0].is_finite());
    }

    #[test]
    fn test_acetic_acid_charges() {
        let atoms = vec![
            Atom {
                symbol: "C".into(),
                atomic_number: 6,
                mass: 12.0,
                charge: 0.0,
                position: [0.0; 3],
                index: 0,
                stereo_parity: 0,
            },
            Atom {
                symbol: "C".into(),
                atomic_number: 6,
                mass: 12.0,
                charge: 0.0,
                position: [0.0; 3],
                index: 1,
                stereo_parity: 0,
            },
            Atom {
                symbol: "O".into(),
                atomic_number: 8,
                mass: 16.0,
                charge: 0.0,
                position: [0.0; 3],
                index: 2,
                stereo_parity: 0,
            },
            Atom {
                symbol: "O".into(),
                atomic_number: 8,
                mass: 16.0,
                charge: 0.0,
                position: [0.0; 3],
                index: 3,
                stereo_parity: 0,
            },
            Atom {
                symbol: "H".into(),
                atomic_number: 1,
                mass: 1.0,
                charge: 0.0,
                position: [0.0; 3],
                index: 4,
                stereo_parity: 0,
            },
            Atom {
                symbol: "H".into(),
                atomic_number: 1,
                mass: 1.0,
                charge: 0.0,
                position: [0.0; 3],
                index: 5,
                stereo_parity: 0,
            },
            Atom {
                symbol: "H".into(),
                atomic_number: 1,
                mass: 1.0,
                charge: 0.0,
                position: [0.0; 3],
                index: 6,
                stereo_parity: 0,
            },
            Atom {
                symbol: "H".into(),
                atomic_number: 1,
                mass: 1.0,
                charge: 0.0,
                position: [0.0; 3],
                index: 7,
                stereo_parity: 0,
            },
        ];
        let bonds = vec![
            Bond {
                atom1: 0,
                atom2: 1,
                bond_type: BondType::Single,
                ..Default::default()
            },
            Bond {
                atom1: 1,
                atom2: 2,
                bond_type: BondType::Double,
                ..Default::default()
            },
            Bond {
                atom1: 1,
                atom2: 3,
                bond_type: BondType::Single,
                ..Default::default()
            },
            Bond {
                atom1: 0,
                atom2: 4,
                bond_type: BondType::Single,
                ..Default::default()
            },
            Bond {
                atom1: 0,
                atom2: 5,
                bond_type: BondType::Single,
                ..Default::default()
            },
            Bond {
                atom1: 0,
                atom2: 6,
                bond_type: BondType::Single,
                ..Default::default()
            },
            Bond {
                atom1: 3,
                atom2: 7,
                bond_type: BondType::Single,
                ..Default::default()
            },
        ];
        let mol = make_molecule(atoms, bonds);
        let atom_types = vec![
            MMFFAtomType::C_3,
            MMFFAtomType::C_2,
            MMFFAtomType::O_2,
            MMFFAtomType::O_3,
            MMFFAtomType::H,
            MMFFAtomType::H,
            MMFFAtomType::H,
            MMFFAtomType::H_COOH,
        ];
        let charges = calculate_bci_charges(&mol, &atom_types);

        let total: f64 = charges.iter().sum();
        assert!(
            total.abs() < 1e-10,
            "acetic acid charges should sum to zero, got {}",
            total
        );

        assert!(
            charges[2] < -0.3,
            "carbonyl O should have significant negative charge, got {}",
            charges[2]
        );
        assert!(
            charges[3] < -0.1,
            "hydroxyl O should have negative charge, got {}",
            charges[3]
        );
        assert!(
            charges[7] > 0.2,
            "acid H should have significant positive charge, got {}",
            charges[7]
        );
        assert!(
            charges[1] > 0.3,
            "carbonyl C should have positive charge, got {}",
            charges[1]
        );

        eprintln!("Acetic acid charges:");
        eprintln!("  C_3 (0): {:.4}", charges[0]);
        eprintln!("  C_2 (1): {:.4}", charges[1]);
        eprintln!("  O_2 (2): {:.4}", charges[2]);
        eprintln!("  O_3 (3): {:.4}", charges[3]);
        eprintln!("  H   (4): {:.4}", charges[4]);
        eprintln!("  H   (5): {:.4}", charges[5]);
        eprintln!("  H   (6): {:.4}", charges[6]);
        eprintln!("  H_COOH(7): {:.4}", charges[7]);
    }

    #[test]
    fn test_bci_lookup_sign_convention() {
        let bci_1_6 = lookup_bci(0, 1, 6).unwrap();
        let bci_6_1 = lookup_bci(0, 6, 1).unwrap();
        assert!(
            (bci_1_6 + bci_6_1).abs() < 1e-10,
            "BCI(1,6) + BCI(6,1) should be zero: {} + {} = {}",
            bci_1_6,
            bci_6_1,
            bci_1_6 + bci_6_1
        );

        assert!(
            bci_1_6 > 0.0,
            "BCI for C_3 bonded to O_3: C should gain positive charge, got {}",
            bci_1_6
        );
    }
}
