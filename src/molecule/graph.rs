//! Molecular graph analysis

use std::collections::{HashSet, VecDeque};

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

/// Find all aromatic atoms in the molecule
pub fn get_aromatic_atoms(mol: &Molecule) -> HashSet<usize> {
    (0..mol.atoms.len())
        .filter(|&i| is_aromatic(i, mol))
        .collect()
}

/// Check if atom is in an aromatic ring (ring membership + Huckel rule)
pub fn is_aromatic(atom_idx: usize, mol: &Molecule) -> bool {
    let rings = find_rings(mol);
    let adj = &mol.adjacency;

    for ring in &rings {
        if !ring.contains(&atom_idx) {
            continue;
        }

        let mut pi_electrons = 0i32;
        let ring_set: HashSet<usize> = ring.iter().copied().collect();

        for &atom in ring {
            let ring_bonds: usize = adj[atom].iter().filter(|&&n| ring_set.contains(&n)).count();

            let mut double_bonds_in_ring = 0;
            for bond in &mol.bonds {
                if (bond.atom1 == atom || bond.atom2 == atom)
                    && ring_set.contains(&bond.atom1)
                    && ring_set.contains(&bond.atom2)
                {
                    match bond.bond_type {
                        BondType::Double | BondType::Aromatic | BondType::Triple => {
                            double_bonds_in_ring += 1;
                        }
                        BondType::Single => {}
                    }
                }
            }

            match mol.atoms[atom].atomic_number {
                6 => {
                    if double_bonds_in_ring >= 1 || ring_bonds == 3 {
                        pi_electrons += 1;
                    }
                }
                7 => {
                    if ring_bonds == 3 && double_bonds_in_ring == 0 {
                        pi_electrons += 2;
                    } else if double_bonds_in_ring >= 1 {
                        pi_electrons += 1;
                    }
                }
                8 => {
                    if ring.len() == 5 && ring_bonds == 2 {
                        pi_electrons += 2;
                    } else if double_bonds_in_ring >= 1 {
                        pi_electrons += 1;
                    }
                }
                16 => {
                    if ring.len() == 5 && ring_bonds == 2 {
                        pi_electrons += 2;
                    } else if double_bonds_in_ring >= 1 {
                        pi_electrons += 1;
                    }
                }
                _ => {
                    if double_bonds_in_ring >= 1 {
                        pi_electrons += 1;
                    }
                }
            }
        }

        if pi_electrons >= 6 && (pi_electrons - 2) % 4 == 0 {
            return true;
        }
    }
    false
}

/// Find smallest set of smallest rings (SSSR) using BFS
pub fn find_rings(mol: &Molecule) -> Vec<Vec<usize>> {
    let n = mol.atoms.len();
    if n < 3 {
        return vec![];
    }

    let adj = &mol.adjacency;
    let mut rings: Vec<Vec<usize>> = Vec::new();

    for start in 0..n {
        let mut parent: Vec<Option<usize>> = vec![None; n];
        let mut visited = vec![false; n];
        let mut queue = VecDeque::new();
        queue.push_back(start);
        visited[start] = true;

        while let Some(v) = queue.pop_front() {
            for &neighbor in &adj[v] {
                if !visited[neighbor] {
                    visited[neighbor] = true;
                    parent[neighbor] = Some(v);
                    queue.push_back(neighbor);
                }
            }
        }

        for bond in &mol.bonds {
            let (a, b) = (bond.atom1, bond.atom2);
            if parent[a] == Some(b) || parent[b] == Some(a) {
                continue;
            }
            if !visited[a] || !visited[b] {
                continue;
            }

            let mut path_a = Vec::new();
            let mut v = a;
            while v != start {
                path_a.push(v);
                v = parent[v].unwrap();
            }
            let mut path_b = Vec::new();
            v = b;
            while v != start {
                path_b.push(v);
                v = parent[v].unwrap();
            }

            path_b.reverse();
            let mut ring = path_a;
            ring.push(start);
            ring.extend(path_b);

            if ring.len() >= 3 {
                if let Some(&min_idx) = ring.iter().min() {
                    let pos = ring.iter().position(|&x| x == min_idx).unwrap();
                    ring.rotate_left(pos);
                }

                let mut is_valid = true;
                for w in 0..ring.len() {
                    let next_w = (w + 1) % ring.len();
                    if !adj[ring[w]].contains(&ring[next_w]) {
                        is_valid = false;
                        break;
                    }
                }
                if is_valid {
                    rings.push(ring);
                }
            }
        }
    }

    rings.sort_by_key(|r| r.len());

    let mut seen: HashSet<Vec<usize>> = HashSet::new();
    let mut unique: Vec<Vec<usize>> = Vec::new();
    for mut ring in rings {
        if let Some(&min_idx) = ring.iter().min() {
            let pos = ring.iter().position(|&x| x == min_idx).unwrap();
            ring.rotate_left(pos);
        }
        let rev: Vec<usize> = ring.iter().rev().copied().collect();
        let canonical = if ring < rev { ring.clone() } else { rev };
        if seen.insert(canonical.clone()) {
            unique.push(canonical);
        }
    }

    let mut primitive: Vec<Vec<usize>> = Vec::new();
    for ring in &unique {
        let ring_set: HashSet<usize> = ring.iter().copied().collect();
        let mut is_prim = true;
        for existing in &primitive {
            if existing.len() < ring.len() && existing.iter().all(|a| ring_set.contains(a)) {
                is_prim = false;
                break;
            }
        }
        if is_prim {
            primitive.push(ring.clone());
        }
    }

    primitive
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

    fn make_atom(symbol: &str, atomic_number: u8, index: usize) -> Atom {
        Atom {
            symbol: symbol.to_string(),
            atomic_number,
            mass: 0.0,
            charge: 0.0,
            position: [0.0; 3],
            index,
        }
    }

    #[test]
    fn test_find_rings_cyclohexane() {
        let atoms: Vec<Atom> = (0..6).map(|i| make_atom("C", 6, i)).collect();
        let mut adjacency = vec![vec![]; 6];
        let bonds: Vec<Bond> = (0..6)
            .map(|i| {
                let j = (i + 1) % 6;
                adjacency[i].push(j);
                adjacency[j].push(i);
                Bond {
                    atom1: i,
                    atom2: j,
                    bond_type: BondType::Single,
                }
            })
            .collect();

        let mol = Molecule {
            atoms,
            bonds,
            name: "Cyclohexane".to_string(),
            adjacency,
        };

        let rings = find_rings(&mol);
        assert_eq!(rings.len(), 1);
        assert_eq!(rings[0].len(), 6);
    }

    #[test]
    fn test_find_rings_benzene() {
        let atoms: Vec<Atom> = (0..6).map(|i| make_atom("C", 6, i)).collect();
        let mut adjacency = vec![vec![]; 6];
        let bonds: Vec<Bond> = (0..6)
            .map(|i| {
                let j = (i + 1) % 6;
                adjacency[i].push(j);
                adjacency[j].push(i);
                Bond {
                    atom1: i,
                    atom2: j,
                    bond_type: BondType::Aromatic,
                }
            })
            .collect();

        let mol = Molecule {
            atoms,
            bonds,
            name: "Benzene".to_string(),
            adjacency,
        };

        let rings = find_rings(&mol);
        assert_eq!(rings.len(), 1);
        assert_eq!(rings[0].len(), 6);
    }

    #[test]
    fn test_find_rings_naphthalene() {
        let atoms: Vec<Atom> = (0..10).map(|i| make_atom("C", 6, i)).collect();
        let mut adjacency = vec![vec![]; 10];
        let mut bonds = Vec::new();

        // Ring 1: 0-1-2-3-4-5-0
        for (i, j) in [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0)] {
            adjacency[i].push(j);
            adjacency[j].push(i);
            bonds.push(Bond {
                atom1: i,
                atom2: j,
                bond_type: BondType::Aromatic,
            });
        }
        // Ring 2: 5-6-7-8-9-4-5
        for (i, j) in [(5, 6), (6, 7), (7, 8), (8, 9), (9, 4)] {
            adjacency[i].push(j);
            adjacency[j].push(i);
            bonds.push(Bond {
                atom1: i,
                atom2: j,
                bond_type: BondType::Aromatic,
            });
        }

        let mol = Molecule {
            atoms,
            bonds,
            name: "Naphthalene".to_string(),
            adjacency,
        };

        let rings = find_rings(&mol);
        assert_eq!(rings.len(), 2);
        assert!(rings.iter().all(|r| r.len() == 6));
    }

    #[test]
    fn test_find_rings_water() {
        let atoms = vec![
            make_atom("O", 8, 0),
            make_atom("H", 1, 1),
            make_atom("H", 1, 2),
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
        let mol = Molecule {
            atoms,
            bonds,
            name: "Water".to_string(),
            adjacency: vec![vec![1, 2], vec![0], vec![0]],
        };

        let rings = find_rings(&mol);
        assert!(rings.is_empty());
    }

    #[test]
    fn test_is_aromatic_benzene() {
        let atoms: Vec<Atom> = (0..6).map(|i| make_atom("C", 6, i)).collect();
        let mut adjacency = vec![vec![]; 6];
        let bonds: Vec<Bond> = (0..6)
            .map(|i| {
                let j = (i + 1) % 6;
                adjacency[i].push(j);
                adjacency[j].push(i);
                Bond {
                    atom1: i,
                    atom2: j,
                    bond_type: BondType::Aromatic,
                }
            })
            .collect();

        let mol = Molecule {
            atoms,
            bonds,
            name: "Benzene".to_string(),
            adjacency,
        };

        let aromatic = get_aromatic_atoms(&mol);
        assert_eq!(aromatic.len(), 6);
    }

    #[test]
    fn test_is_aromatic_cyclohexane() {
        let atoms: Vec<Atom> = (0..6).map(|i| make_atom("C", 6, i)).collect();
        let mut adjacency = vec![vec![]; 6];
        let bonds: Vec<Bond> = (0..6)
            .map(|i| {
                let j = (i + 1) % 6;
                adjacency[i].push(j);
                adjacency[j].push(i);
                Bond {
                    atom1: i,
                    atom2: j,
                    bond_type: BondType::Single,
                }
            })
            .collect();

        let mol = Molecule {
            atoms,
            bonds,
            name: "Cyclohexane".to_string(),
            adjacency,
        };

        let aromatic = get_aromatic_atoms(&mol);
        assert!(aromatic.is_empty());
    }

    #[test]
    fn test_is_aromatic_pyridine() {
        // Pyridine: C5H5N — 6-membered ring with one N replacing C
        // Atoms: 0-4 are C, 5 is N
        let mut atoms: Vec<Atom> = (0..5).map(|i| make_atom("C", 6, i)).collect();
        atoms.push(make_atom("N", 7, 5));
        let mut adjacency = vec![vec![]; 6];
        let mut bonds = Vec::new();

        // Ring: 0-1-2-3-4-5-0
        for (i, j) in [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0)] {
            adjacency[i].push(j);
            adjacency[j].push(i);
            bonds.push(Bond {
                atom1: i,
                atom2: j,
                bond_type: BondType::Aromatic,
            });
        }

        let mol = Molecule {
            atoms,
            bonds,
            name: "Pyridine".to_string(),
            adjacency,
        };

        let aromatic = get_aromatic_atoms(&mol);
        assert_eq!(aromatic.len(), 6);
    }

    #[test]
    fn test_is_aromatic_thiophene() {
        // Thiophene: 5-membered ring with S
        // Atoms: 0 is S, 1-4 are C
        let mut atoms: Vec<Atom> = vec![make_atom("S", 16, 0)];
        for i in 1..=4 {
            atoms.push(make_atom("C", 6, i));
        }
        let mut adjacency = vec![vec![]; 5];
        let mut bonds = Vec::new();

        // Ring: 0(S)-1(C)-2(C)-3(C)-4(C)-0
        for (i, j) in [(0, 1), (1, 2), (2, 3), (3, 4), (4, 0)] {
            adjacency[i].push(j);
            adjacency[j].push(i);
            bonds.push(Bond {
                atom1: i,
                atom2: j,
                bond_type: BondType::Aromatic,
            });
        }
        // Alternating double bonds: 1=2, 3=4
        bonds[1].bond_type = BondType::Double;
        bonds[3].bond_type = BondType::Double;

        let mol = Molecule {
            atoms,
            bonds,
            name: "Thiophene".to_string(),
            adjacency,
        };

        let aromatic = get_aromatic_atoms(&mol);
        assert_eq!(aromatic.len(), 5);
    }
}
