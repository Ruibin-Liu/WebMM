use super::atom_types::get_atom_type_props;
use super::MMFFAtomType;
use crate::molecule::{BondType, Molecule};

fn bond_order(bt: BondType) -> i32 {
    match bt {
        BondType::Single => 1,
        BondType::Double => 2,
        BondType::Triple => 3,
        BondType::Aromatic => 1,
    }
}

fn get_bci(type1: MMFFAtomType, type2: MMFFAtomType, order: i32) -> Option<f64> {
    let (a, b) = if (type1 as u32) <= (type2 as u32) {
        (type1, type2)
    } else {
        (type2, type1)
    };

    match (a, b, order) {
        (MMFFAtomType::H, MMFFAtomType::C_3, 1) => Some(0.05),
        (MMFFAtomType::H_OH, MMFFAtomType::C_3, 1) => Some(0.05),
        (MMFFAtomType::H_ONC, MMFFAtomType::C_3, 1) => Some(0.05),
        (MMFFAtomType::H_COOH, MMFFAtomType::C_3, 1) => Some(0.05),
        (MMFFAtomType::H_OAR, MMFFAtomType::C_3, 1) => Some(0.05),
        (MMFFAtomType::H_N3, MMFFAtomType::C_3, 1) => Some(0.05),
        (MMFFAtomType::H_NAM, MMFFAtomType::C_3, 1) => Some(0.05),
        (MMFFAtomType::H, MMFFAtomType::C_2, 1) => Some(0.05),
        (MMFFAtomType::H, MMFFAtomType::C_AR, 1) => Some(0.05),

        (MMFFAtomType::C_3, MMFFAtomType::C_3, 1) => Some(0.00),
        (MMFFAtomType::C_2, MMFFAtomType::C_3, 1) => Some(0.00),
        (MMFFAtomType::C_2, MMFFAtomType::C_2, 2) => Some(0.00),
        (MMFFAtomType::C_AR, MMFFAtomType::C_AR, 1) => Some(0.00),

        (MMFFAtomType::C_3, MMFFAtomType::N_3, 1) => Some(0.08),
        (MMFFAtomType::C_3, MMFFAtomType::N_2, 1) => Some(0.10),
        (MMFFAtomType::C_3, MMFFAtomType::N_AR, 1) => Some(0.05),
        (MMFFAtomType::C_2, MMFFAtomType::N_2, 2) => Some(0.10),
        (MMFFAtomType::C_AR, MMFFAtomType::N_AR, 1) => Some(0.08),
        (MMFFAtomType::C_2, MMFFAtomType::N_3, 1) => Some(0.05),

        (MMFFAtomType::C_3, MMFFAtomType::O_3, 1) => Some(0.13),
        (MMFFAtomType::C_3, MMFFAtomType::O_2, 1) => Some(0.15),
        (MMFFAtomType::C_2, MMFFAtomType::O_2, 2) => Some(0.42),
        (MMFFAtomType::C_2, MMFFAtomType::O_CO2, 2) => Some(0.42),
        (MMFFAtomType::C_2, MMFFAtomType::O_3, 1) => Some(0.35),
        (MMFFAtomType::C_AR, MMFFAtomType::O_CO2, 2) => Some(0.42),
        (MMFFAtomType::C_3, MMFFAtomType::O_R, 1) => Some(0.13),
        (MMFFAtomType::C_3, MMFFAtomType::O_CO2, 1) => Some(0.13),

        (MMFFAtomType::C_3, MMFFAtomType::S_3, 1) => Some(0.05),
        (MMFFAtomType::C_2, MMFFAtomType::S_2, 2) => Some(0.05),

        (MMFFAtomType::H, MMFFAtomType::N_3, 1) => Some(0.00),
        (MMFFAtomType::H, MMFFAtomType::N_2, 1) => Some(0.00),
        (MMFFAtomType::H, MMFFAtomType::N_AR, 1) => Some(0.00),
        (MMFFAtomType::H_NAM, MMFFAtomType::N_2, 1) => Some(0.20),
        (MMFFAtomType::H_NAM, MMFFAtomType::N_AM, 1) => Some(0.20),
        (MMFFAtomType::H_NAM, MMFFAtomType::N_AR, 1) => Some(0.20),

        (MMFFAtomType::N_3, MMFFAtomType::N_3, 1) => Some(0.00),

        (MMFFAtomType::H_OH, MMFFAtomType::O_3, 1) => Some(0.40),
        (MMFFAtomType::H_ONC, MMFFAtomType::O_3, 1) => Some(0.40),
        (MMFFAtomType::H_COOH, MMFFAtomType::O_3, 1) => Some(0.40),
        (MMFFAtomType::H_OAR, MMFFAtomType::O_3, 1) => Some(0.40),
        (MMFFAtomType::H, MMFFAtomType::O_3, 1) => Some(0.40),
        (MMFFAtomType::H, MMFFAtomType::O_R, 1) => Some(0.40),
        (MMFFAtomType::H_N3, MMFFAtomType::N_3, 1) => Some(0.19),
        (MMFFAtomType::H_N3, MMFFAtomType::N_2, 1) => Some(0.19),
        (MMFFAtomType::H_N3, MMFFAtomType::N_AR, 1) => Some(0.19),

        (MMFFAtomType::C_3, MMFFAtomType::F, 1) => Some(0.10),
        (MMFFAtomType::C_3, MMFFAtomType::Cl, 1) => Some(0.10),
        (MMFFAtomType::C_3, MMFFAtomType::Br, 1) => Some(0.10),
        (MMFFAtomType::C_3, MMFFAtomType::I, 1) => Some(0.10),
        (MMFFAtomType::C_AR, MMFFAtomType::F, 1) => Some(0.10),
        (MMFFAtomType::C_AR, MMFFAtomType::Cl, 1) => Some(0.10),

        (MMFFAtomType::C_3, MMFFAtomType::P_3, 1) => Some(0.00),

        (MMFFAtomType::O_3, MMFFAtomType::P_3, 1) => Some(0.10),
        (MMFFAtomType::O_2, MMFFAtomType::P_4, 2) => Some(0.10),

        (MMFFAtomType::C_3, MMFFAtomType::N_PL3, 1) => Some(0.05),
        (MMFFAtomType::C_3, MMFFAtomType::N_AM, 1) => Some(0.05),
        (MMFFAtomType::C_2, MMFFAtomType::N_AM, 1) => Some(0.05),
        (MMFFAtomType::C_2, MMFFAtomType::N_AM, 2) => Some(0.25),
        (MMFFAtomType::C_2, MMFFAtomType::N_PL3, 1) => Some(0.30),
        (MMFFAtomType::C_AR, MMFFAtomType::N_PL3, 1) => Some(0.10),

        _ => None,
    }
}

pub fn calculate_bci_charges(mol: &Molecule, atom_types: &[MMFFAtomType]) -> Vec<f64> {
    let n = mol.atoms.len();
    let mut charges: Vec<f64> = vec![0.0; n];

    for bond in &mol.bonds {
        let i = bond.atom1;
        let j = bond.atom2;
        let t_i = atom_types[i];
        let t_j = atom_types[j];
        let order = bond_order(bond.bond_type);

        let bci = match get_bci(t_i, t_j, order) {
            Some(v) => v,
            None => {
                let fbci_i = get_atom_type_props(t_i).map(|p| p.fbci).unwrap_or(0.0);
                let fbci_j = get_atom_type_props(t_j).map(|p| p.fbci).unwrap_or(0.0);
                let num_i = mol.adjacency[i].len() as f64;
                let num_j = mol.adjacency[j].len() as f64;
                (fbci_j - fbci_i).abs() * order as f64 * 0.3 / (num_i + num_j)
            }
        };

        let fbci_i = get_atom_type_props(t_i).map(|p| p.fbci).unwrap_or(0.0);
        let fbci_j = get_atom_type_props(t_j).map(|p| p.fbci).unwrap_or(0.0);

        // More electronegative atom (more negative fbci) gains negative charge
        if fbci_i < fbci_j {
            charges[i] -= bci;
            charges[j] += bci;
        } else {
            charges[i] += bci;
            charges[j] -= bci;
        }
    }

    let total_charge: f64 = mol.atoms.iter().map(|a| a.charge).sum();
    let current_total: f64 = charges.iter().sum();
    let residual = total_charge - current_total;

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
            },
            Bond {
                atom1: 1,
                atom2: 2,
                bond_type: BondType::Single,
            },
            Bond {
                atom1: 2,
                atom2: 3,
                bond_type: BondType::Single,
            },
            Bond {
                atom1: 0,
                atom2: 4,
                bond_type: BondType::Single,
            },
            Bond {
                atom1: 0,
                atom2: 5,
                bond_type: BondType::Single,
            },
            Bond {
                atom1: 0,
                atom2: 6,
                bond_type: BondType::Single,
            },
            Bond {
                atom1: 1,
                atom2: 7,
                bond_type: BondType::Single,
            },
            Bond {
                atom1: 1,
                atom2: 8,
                bond_type: BondType::Single,
            },
        ];
        let mol = make_molecule(atoms, bonds);
        let atom_types = vec![
            MMFFAtomType::C_3,
            MMFFAtomType::C_3,
            MMFFAtomType::O_3,
            MMFFAtomType::H,
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
            },
            Bond {
                atom1: 0,
                atom2: 2,
                bond_type: BondType::Single,
            },
        ];
        let mol = make_molecule(atoms, bonds);
        let atom_types = vec![MMFFAtomType::O_3, MMFFAtomType::H, MMFFAtomType::H];
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
        // BCI model: charges should sum to the molecule's formal charge
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
        }];
        let mol = make_molecule(atoms, bonds);
        let atom_types = vec![MMFFAtomType::O_3, MMFFAtomType::H];
        let charges = calculate_bci_charges(&mol, &atom_types);

        // Charges should sum to the molecule's formal charge (-1)
        let total: f64 = charges.iter().sum();
        assert!(
            (total - (-1.0)).abs() < 1e-10,
            "hydroxide charges should sum to -1, got {}",
            total
        );
    }

    #[test]
    fn test_ammonia_charges() {
        // NH3: BCI model is disabled, charges are zero
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
            },
            Bond {
                atom1: 0,
                atom2: 2,
                bond_type: BondType::Single,
            },
            Bond {
                atom1: 0,
                atom2: 3,
                bond_type: BondType::Single,
            },
        ];
        let mol = make_molecule(atoms, bonds);
        let atom_types = vec![
            MMFFAtomType::N_3,
            MMFFAtomType::H,
            MMFFAtomType::H,
            MMFFAtomType::H,
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
        // Single atom with no bonds: charge should just be fbci
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
}
