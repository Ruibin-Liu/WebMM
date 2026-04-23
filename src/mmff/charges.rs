use super::MMFFAtomType;
use crate::molecule::{BondType, Molecule};

fn mmff_type_id(t: MMFFAtomType) -> u8 {
    match t {
        MMFFAtomType::H => 5,
        MMFFAtomType::H_OH => 31,
        MMFFAtomType::H_ONC => 21,
        MMFFAtomType::H_COOH => 24,
        MMFFAtomType::H_OAR => 29,
        MMFFAtomType::H_N3 => 23,
        MMFFAtomType::H_NAM => 28,
        MMFFAtomType::C_3 => 1,
        MMFFAtomType::C_2 => 3,
        MMFFAtomType::C_1 => 4,
        MMFFAtomType::C_AR => 37,
        MMFFAtomType::C5A => 63,
        MMFFAtomType::C5B => 64,
        MMFFAtomType::C_CAT => 56,
        MMFFAtomType::C_AN => 57,
        MMFFAtomType::N_3 => 8,
        MMFFAtomType::N_2 => 9,
        MMFFAtomType::N_1 => 42,
        MMFFAtomType::N_AR => 38,
        MMFFAtomType::NPYL => 39,
        MMFFAtomType::N_PL3 => 40,
        MMFFAtomType::N_AM => 10,
        MMFFAtomType::N_4 => 34,
        MMFFAtomType::N_2Z => 53,
        MMFFAtomType::N_SOM => 48,
        MMFFAtomType::N5A => 65,
        MMFFAtomType::N5B => 66,
        MMFFAtomType::O_3 => 6,
        MMFFAtomType::O_2 => 7,
        MMFFAtomType::O_R => 6,  // alcohol/ether O is type 6, not 70
        MMFFAtomType::OH2 => 70,
        MMFFAtomType::OFUR => 59,
        MMFFAtomType::O_CO2 => 32,
        MMFFAtomType::O_3_Z => 35,
        MMFFAtomType::F => 11,
        MMFFAtomType::Cl => 12,
        MMFFAtomType::Br => 13,
        MMFFAtomType::I => 14,
        MMFFAtomType::S_3 => 15,
        MMFFAtomType::S_2 => 16,
        MMFFAtomType::S_AR => 44,
        MMFFAtomType::P_3 => 26,
        MMFFAtomType::P_4 => 25,
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

#[allow(dead_code)]
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

fn mmff_bond_type(bt: BondType) -> u8 {
    match bt {
        BondType::Single => 0,
        BondType::Double => 1,
        BondType::Triple => 2,
        BondType::Aromatic => 4,
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

fn lookup_bci_canonical(bt: u8, i: u8, j: u8) -> Option<f64> {
    let entries: BciEntries = BciEntries;
    entries.lookup(bt, i, j)
}

struct BciEntries;

impl BciEntries {
    fn lookup(&self, bt: u8, i: u8, j: u8) -> Option<f64> {
        match (bt, i, j) {
            (0, 1, 1) => Some(0.0000),
            (0, 1, 3) => Some(-0.0610),
            (0, 1, 5) => Some(0.0000),
            (0, 1, 6) => Some(-0.2800),
            (0, 1, 8) => Some(-0.2700),
            (0, 1, 11) => Some(-0.3400),
            (0, 1, 12) => Some(-0.2900),
            (0, 1, 13) => Some(-0.2300),
            (0, 1, 14) => Some(-0.1900),
            (0, 1, 15) => Some(-0.2300),
            (0, 1, 17) => Some(-0.1935),
            (0, 2, 5) => Some(0.1500),
            (0, 3, 5) => Some(0.0600),
            (0, 3, 6) => Some(-0.1500),
            (0, 3, 10) => Some(-0.0600),
            (0, 4, 5) => Some(0.1770),
            (0, 5, 37) => Some(-0.1500),
            (0, 6, 6) => Some(0.0000),
            (0, 6, 21) => Some(0.4000),
            (0, 6, 24) => Some(0.5000),
            (0, 6, 29) => Some(0.4500),
            (0, 6, 37) => Some(0.0825),
            (0, 8, 23) => Some(0.3600),
            (0, 10, 28) => Some(0.3700),
            (0, 15, 71) => Some(0.1800),
            (0, 40, 28) => Some(0.4000),
            (0, 70, 31) => Some(0.4300),
            (1, 2, 2) => Some(0.0000),
            (1, 3, 7) => Some(-0.5700),
            (1, 7, 17) => Some(0.5000),
            (2, 4, 4) => Some(0.0000),
            (4, 37, 37) => Some(0.0000),
            _ => None,
        }
    }
}

pub fn calculate_bci_charges(mol: &Molecule, atom_types: &[MMFFAtomType]) -> Vec<f64> {
    let n = mol.atoms.len();
    let type_ids: Vec<u8> = atom_types.iter().map(|t| mmff_type_id(*t)).collect();
    let mut charges: Vec<f64> = vec![0.0; n];

    for bond in &mol.bonds {
        let i = bond.atom1;
        let j = bond.atom2;
        let bt = mmff_bond_type(bond.bond_type);
        let ti = type_ids[i];
        let tj = type_ids[j];

        let bci = match lookup_bci(bt, ti, tj) {
            Some(v) => v,
            None => {
                let pi = get_pbci(ti);
                let pj = get_pbci(tj);
                pi.pbci - pj.pbci
            }
        };

        charges[i] += bci;
        charges[j] -= bci;
    }

    let total_formal: f64 = mol.atoms.iter().map(|a| a.charge).sum();
    let current_total: f64 = charges.iter().sum();
    let residual = total_formal - current_total;

    if n > 0 && residual.abs() > 1e-12 {
        let correction = residual / n as f64;
        for c in &mut charges {
            *c += correction;
        }
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
            },
            Atom {
                symbol: "C".into(),
                atomic_number: 6,
                mass: 12.0,
                charge: 0.0,
                position: [0.0; 3],
                index: 1,
            },
            Atom {
                symbol: "O".into(),
                atomic_number: 8,
                mass: 16.0,
                charge: 0.0,
                position: [0.0; 3],
                index: 2,
            },
            Atom {
                symbol: "H".into(),
                atomic_number: 1,
                mass: 1.0,
                charge: 0.0,
                position: [0.0; 3],
                index: 3,
            },
            Atom {
                symbol: "H".into(),
                atomic_number: 1,
                mass: 1.0,
                charge: 0.0,
                position: [0.0; 3],
                index: 4,
            },
            Atom {
                symbol: "H".into(),
                atomic_number: 1,
                mass: 1.0,
                charge: 0.0,
                position: [0.0; 3],
                index: 5,
            },
            Atom {
                symbol: "H".into(),
                atomic_number: 1,
                mass: 1.0,
                charge: 0.0,
                position: [0.0; 3],
                index: 6,
            },
            Atom {
                symbol: "H".into(),
                atomic_number: 1,
                mass: 1.0,
                charge: 0.0,
                position: [0.0; 3],
                index: 7,
            },
            Atom {
                symbol: "H".into(),
                atomic_number: 1,
                mass: 1.0,
                charge: 0.0,
                position: [0.0; 3],
                index: 8,
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
            },
            Atom {
                symbol: "H".into(),
                atomic_number: 1,
                mass: 1.0,
                charge: 0.0,
                position: [0.0; 3],
                index: 1,
            },
            Atom {
                symbol: "H".into(),
                atomic_number: 1,
                mass: 1.0,
                charge: 0.0,
                position: [0.0; 3],
                index: 2,
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
            },
            Atom {
                symbol: "H".into(),
                atomic_number: 1,
                mass: 1.0,
                charge: 0.0,
                position: [0.0; 3],
                index: 1,
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

        let total: f64 = charges.iter().sum();
        assert!(
            (total - (-1.0)).abs() < 1e-10,
            "hydroxide charges should sum to -1, got {}",
            total
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
            },
            Atom {
                symbol: "H".into(),
                atomic_number: 1,
                mass: 1.0,
                charge: 0.0,
                position: [1.0, 0.0, 0.0],
                index: 1,
            },
            Atom {
                symbol: "H".into(),
                atomic_number: 1,
                mass: 1.0,
                charge: 0.0,
                position: [-0.5, 0.866, 0.0],
                index: 2,
            },
            Atom {
                symbol: "H".into(),
                atomic_number: 1,
                mass: 1.0,
                charge: 0.0,
                position: [-0.5, -0.866, 0.0],
                index: 3,
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
            },
            Atom {
                symbol: "C".into(),
                atomic_number: 6,
                mass: 12.0,
                charge: 0.0,
                position: [0.0; 3],
                index: 1,
            },
            Atom {
                symbol: "O".into(),
                atomic_number: 8,
                mass: 16.0,
                charge: 0.0,
                position: [0.0; 3],
                index: 2,
            },
            Atom {
                symbol: "O".into(),
                atomic_number: 8,
                mass: 16.0,
                charge: 0.0,
                position: [0.0; 3],
                index: 3,
            },
            Atom {
                symbol: "H".into(),
                atomic_number: 1,
                mass: 1.0,
                charge: 0.0,
                position: [0.0; 3],
                index: 4,
            },
            Atom {
                symbol: "H".into(),
                atomic_number: 1,
                mass: 1.0,
                charge: 0.0,
                position: [0.0; 3],
                index: 5,
            },
            Atom {
                symbol: "H".into(),
                atomic_number: 1,
                mass: 1.0,
                charge: 0.0,
                position: [0.0; 3],
                index: 6,
            },
            Atom {
                symbol: "H".into(),
                atomic_number: 1,
                mass: 1.0,
                charge: 0.0,
                position: [0.0; 3],
                index: 7,
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
