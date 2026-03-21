//! Molecular graph analysis

use super::{BondType, Molecule};

/// Build adjacency list from bonds (standalone, takes counts and bonds directly)
pub fn build_adjacency_list_from_bonds(
    num_atoms: &usize,
    bonds: &[super::Bond],
) -> Vec<Vec<usize>> {
    let mut adj = vec![vec![]; *num_atoms];

    for bond in bonds {
        if bond.atom1 < *num_atoms && bond.atom2 < *num_atoms {
            adj[bond.atom1].push(bond.atom2);
            adj[bond.atom2].push(bond.atom1);
        }
    }

    adj
}

/// Build adjacency list from bonds (uses cached adjacency if available)
pub fn build_adjacency_list(mol: &Molecule) -> Vec<Vec<usize>> {
    mol.adjacency.clone()
}

/// Get bonded neighbors of an atom (returns cached reference)
pub fn get_neighbors(atom_idx: usize, mol: &Molecule) -> &[usize] {
    if atom_idx < mol.adjacency.len() {
        &mol.adjacency[atom_idx]
    } else {
        &[]
    }
}

/// Determine hybridization (sp3, sp2, sp1)
pub fn determine_hybridization(atom_idx: usize, mol: &Molecule) -> Hybridization {
    let neighbors = get_neighbors(atom_idx, mol);
    let num_bonds = neighbors.len();

    let atom = &mol.atoms[atom_idx];
    let symbol = &atom.symbol;

    let pi_bonds: u8 = mol
        .bonds
        .iter()
        .filter(|b| b.atom1 == atom_idx || b.atom2 == atom_idx)
        .map(|b| match b.bond_type {
            BondType::Single => 0,
            BondType::Double => 1,
            BondType::Triple => 2,
            BondType::Aromatic => 1,
        })
        .sum();

    if pi_bonds >= 2 {
        return Hybridization::Sp1;
    }
    if pi_bonds == 1 {
        return Hybridization::Sp2;
    }

    match symbol.as_str() {
        "C" => match num_bonds {
            1 => Hybridization::Sp1,
            2 => Hybridization::Sp2,
            3 => Hybridization::Sp3,
            4 => Hybridization::Sp3,
            _ => Hybridization::Sp3,
        },
        "N" => match num_bonds {
            1 => Hybridization::Sp1,
            2 => Hybridization::Sp2,
            3 => Hybridization::Sp3,
            _ => Hybridization::Sp3,
        },
        "O" => match num_bonds {
            1 => Hybridization::Sp2,
            2 => Hybridization::Sp3,
            3 => Hybridization::Sp3,
            _ => Hybridization::Sp3,
        },
        "P" => match num_bonds {
            2 => Hybridization::Sp2,
            3 => Hybridization::Sp3,
            _ => Hybridization::Sp3,
        },
        "S" => match num_bonds {
            2 => Hybridization::Sp2,
            3 => Hybridization::Sp3,
            _ => Hybridization::Sp3,
        },
        _ => Hybridization::Sp3,
    }
}

/// Check if atom is in an aromatic ring
pub fn is_aromatic(atom_idx: usize, mol: &Molecule) -> bool {
    // TODO: Implement aromaticity detection
    // For now, check if any bonded bond is aromatic
    for bond in &mol.bonds {
        if (bond.atom1 == atom_idx || bond.atom2 == atom_idx)
            && bond.bond_type == BondType::Aromatic
        {
            return true;
        }
    }
    false
}

/// Find rings in molecule
pub fn find_rings(_mol: &Molecule) -> Vec<Vec<usize>> {
    // TODO: Implement ring detection
    // Use depth-first search or union-find for ring detection
    vec![]
}

/// Find rotatable bonds (for torsion angles)
pub fn find_rotatable_bonds(mol: &Molecule) -> Vec<(usize, usize)> {
    let mut rotatable = Vec::new();

    for bond in &mol.bonds {
        // Bond is rotatable if both atoms have more than 1 heavy atom neighbor
        let neighbors1 = get_neighbors(bond.atom1, mol);
        let neighbors2 = get_neighbors(bond.atom2, mol);

        let heavy_neighbors1: usize = neighbors1
            .iter()
            .filter(|&n| mol.atoms[*n].atomic_number != 1)
            .count();

        let heavy_neighbors2: usize = neighbors2
            .iter()
            .filter(|&n| mol.atoms[*n].atomic_number != 1)
            .count();

        if heavy_neighbors1 > 1 && heavy_neighbors2 > 1 {
            rotatable.push((bond.atom1, bond.atom2));
        }
    }

    rotatable
}

/// Hybridization types
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Hybridization {
    Sp3,
    Sp2,
    Sp1,
}

/// Angle in molecule (for angle bending term)
pub struct Angle {
    pub atom1: usize,
    pub atom2: usize,
    pub atom3: usize,
}

/// Find all angles in molecule
pub fn find_angles(mol: &Molecule) -> Vec<Angle> {
    let mut angles = Vec::new();

    for bond in &mol.bonds {
        let neighbors1 = get_neighbors(bond.atom1, mol);
        let neighbors2 = get_neighbors(bond.atom2, mol);

        // Find angles centered on each atom of the bond
        for &n1 in neighbors1 {
            if n1 != bond.atom2 {
                angles.push(Angle {
                    atom1: n1,
                    atom2: bond.atom1,
                    atom3: bond.atom2,
                });
            }
        }

        for &n2 in neighbors2 {
            if n2 != bond.atom1 {
                angles.push(Angle {
                    atom1: bond.atom1,
                    atom2: bond.atom2,
                    atom3: n2,
                });
            }
        }
    }

    angles
}

/// Torsion in molecule (4 atoms in sequence)
pub struct Torsion {
    pub atom1: usize,
    pub atom2: usize,
    pub atom3: usize,
    pub atom4: usize,
}

/// Find all torsions in molecule
pub fn find_torsions(mol: &Molecule) -> Vec<Torsion> {
    let mut torsions = Vec::new();
    let rotatable_bonds = find_rotatable_bonds(mol);

    for (i, j) in rotatable_bonds {
        let neighbors_i = get_neighbors(i, mol);
        let neighbors_j = get_neighbors(j, mol);

        // Find atoms to form torsion i-k-j-l
        for &k in neighbors_i {
            if k != j {
                for &l in neighbors_j {
                    if l != i {
                        torsions.push(Torsion {
                            atom1: k,
                            atom2: i,
                            atom3: j,
                            atom4: l,
                        });
                    }
                }
            }
        }
    }

    torsions
}

/// Out-of-plane (central atom with 3 bonded atoms)
pub struct OutOfPlane {
    pub central: usize,
    pub atom1: usize,
    pub atom2: usize,
    pub atom3: usize,
}

/// Find all out-of-plane groups
pub fn find_out_of_planes(mol: &Molecule) -> Vec<OutOfPlane> {
    let mut oops = Vec::new();

    for atom_idx in 0..mol.atoms.len() {
        let neighbors = get_neighbors(atom_idx, mol);

        // Only atoms with 3+ neighbors can have out-of-plane bending
        if neighbors.len() >= 3 {
            // Only sp2 and aromatic atoms typically have significant OOP
            let hybrid = determine_hybridization(atom_idx, mol);
            if hybrid == Hybridization::Sp2 || is_aromatic(atom_idx, mol) {
                // Create OOP combinations (choose 3 out of N neighbors)
                for i in 0..neighbors.len() {
                    for j in (i + 1)..neighbors.len() {
                        for k in (j + 1)..neighbors.len() {
                            oops.push(OutOfPlane {
                                central: atom_idx,
                                atom1: neighbors[i],
                                atom2: neighbors[j],
                                atom3: neighbors[k],
                            });
                        }
                    }
                }
            }
        }
    }

    oops
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::molecule::{Atom, Bond};

    #[test]
    fn test_build_adjacency_list() {
        let atoms = vec![
            Atom {
                symbol: "C".to_string(),
                atomic_number: 6,
                mass: 12.011,
                charge: 0.0,
                position: [0.0; 3],
                index: 0,
            },
            Atom {
                symbol: "C".to_string(),
                atomic_number: 6,
                mass: 12.011,
                charge: 0.0,
                position: [0.0; 3],
                index: 1,
            },
            Atom {
                symbol: "H".to_string(),
                atomic_number: 1,
                mass: 1.008,
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
                atom1: 1,
                atom2: 2,
                bond_type: BondType::Single,
            },
        ];

        let mol = Molecule {
            atoms,
            bonds,
            name: "Ethane".to_string(),
            adjacency: vec![vec![1], vec![0, 2], vec![1]],
        };
        let adj = build_adjacency_list(&mol);

        assert_eq!(adj.len(), 3);
        assert_eq!(adj[0], vec![1]);
        assert_eq!(adj[1], vec![0, 2]);
        assert_eq!(adj[2], vec![1]);
    }
}
