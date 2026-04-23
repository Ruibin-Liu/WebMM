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

/// Convert bonds in aromatic rings to BondType::Aromatic
/// This matches RDKit behavior where Kekulé structures are perceived as aromatic
pub fn perceive_aromatic_bonds(mol: &mut Molecule) {
    let aromatic_atoms = get_aromatic_atoms(mol);
    let rings = find_rings(mol);

    for ring in &rings {
        let ring_set: HashSet<usize> = ring.iter().copied().collect();
        let all_aromatic = ring.iter().all(|a| aromatic_atoms.contains(a));
        if !all_aromatic {
            continue;
        }

        for bond in mol.bonds.iter_mut() {
            if ring_set.contains(&bond.atom1) && ring_set.contains(&bond.atom2) {
                bond.bond_type = BondType::Aromatic;
            }
        }
    }
}

/// Check if atom has a multiple bond (double, triple, or aromatic) within the given ring
fn has_multiple_bond_in_ring(atom_idx: usize, ring_set: &HashSet<usize>, mol: &Molecule) -> bool {
    mol.bonds.iter().any(|bond| {
        (bond.atom1 == atom_idx || bond.atom2 == atom_idx)
            && ring_set.contains(&bond.atom1)
            && ring_set.contains(&bond.atom2)
            && matches!(bond.bond_type, BondType::Double | BondType::Triple | BondType::Aromatic)
    })
}

/// Check if a ring contains a heteroatom that can donate lone pairs (N, O, S)
fn ring_has_heteroatom(ring: &[usize], mol: &Molecule) -> bool {
    ring.iter()
        .any(|&a| matches!(mol.atoms[a].atomic_number, 7 | 8 | 16))
}

/// Estimate total neighbor count including implicit hydrogens
fn estimate_total_neighbors(atom_idx: usize, mol: &Molecule) -> usize {
    let explicit_neighbors = mol.adjacency[atom_idx].len();
    let explicit_valence: u32 = mol
        .bonds
        .iter()
        .filter(|b| b.atom1 == atom_idx || b.atom2 == atom_idx)
        .map(|b| match b.bond_type {
            BondType::Single => 1,
            BondType::Double => 2,
            BondType::Triple => 3,
            BondType::Aromatic => 1,
        })
        .sum();

    let typical_valence = match mol.atoms[atom_idx].atomic_number {
        6 => 4,
        7 => 3,
        8 => 2,
        16 => 2,
        _ => return explicit_neighbors,
    };

    let implicit_h = (typical_valence as i32 - explicit_valence as i32).max(0) as usize;
    explicit_neighbors + implicit_h
}

/// Determine if an atom is a candidate for aromaticity in the given ring.
/// Matches RDKit behavior where only atoms that can contribute electrons
/// (or have empty p-orbitals) are considered.
fn is_aromatic_candidate(atom_idx: usize, ring: &[usize], mol: &Molecule) -> bool {
    let atom = &mol.atoms[atom_idx];
    let ring_set: HashSet<usize> = ring.iter().copied().collect();
    let ring_bonds = mol.adjacency[atom_idx]
        .iter()
        .filter(|&&n| ring_set.contains(&n))
        .count();
    let has_multiple = has_multiple_bond_in_ring(atom_idx, &ring_set, mol);

    match atom.atomic_number {
        6 => {
            if has_multiple {
                true
            } else if ring.len() == 5 && ring_has_heteroatom(ring, mol) {
                // In 5-membered heteroaromatics (furan, thiophene, imidazole),
                // a C with only single bonds may still be sp2 if the ring
                // has enough multiple bonds to suggest conjugation.
                let multiple_bonds_in_ring = mol
                    .bonds
                    .iter()
                    .filter(|b| ring_set.contains(&b.atom1) && ring_set.contains(&b.atom2))
                    .filter(|b| {
                        matches!(
                            b.bond_type,
                            BondType::Double | BondType::Triple | BondType::Aromatic
                        )
                    })
                    .count();
                multiple_bonds_in_ring >= 2
            } else {
                // Saturated C (e.g. cyclohexane, cyclohexene sp3 C's)
                estimate_total_neighbors(atom_idx, mol) <= 3
            }
        }
        7 => {
            ring.len() == 5 && ring_bonds == 2
                || has_multiple
                || (ring_bonds == 3 && !has_multiple)
        }
        8 | 16 => {
            ring.len() == 5 && ring_bonds == 2 || has_multiple
        }
        _ => has_multiple,
    }
}

/// Count pi electrons contributed by an atom to the aromatic ring.
/// Assumes the atom is already confirmed as an aromatic candidate.
/// Works correctly even after bonds have been upgraded to Aromatic.
fn count_pi_electrons(atom_idx: usize, ring: &[usize], mol: &Molecule) -> i32 {
    let atom = &mol.atoms[atom_idx];
    let ring_set: HashSet<usize> = ring.iter().copied().collect();
    let ring_bonds = mol.adjacency[atom_idx]
        .iter()
        .filter(|&&n| ring_set.contains(&n))
        .count();
    // Total neighbors (including explicit H/exocyclic substituents).
    // Used to distinguish pyrrole-like N (3 neighbors) from pyridine-like N (2 neighbors).
    let total_neighbors = mol.adjacency[atom_idx].len();

    match atom.atomic_number {
        6 => 1,
        7 => {
            if ring.len() == 5 && ring_bonds == 2 && total_neighbors >= 3 {
                // Pyrrole-like N in 5-membered ring (has H or substituent)
                2
            } else if ring_bonds == 3 && total_neighbors >= 3 {
                // Aniline-like N (3 ring bonds, has H or substituent)
                2
            } else {
                1
            }
        }
        8 | 16 => {
            if ring.len() == 5 && ring_bonds == 2 {
                // Furan/thiophene-like O/S always contributes lone pair
                2
            } else {
                1
            }
        }
        _ => 1,
    }
}

/// Check if atom is in an aromatic ring (ring membership + Huckel rule).
/// Uses RDKit-style aromaticity perception:
/// - All atoms in the ring must be aromatic candidates.
/// - Pi electrons are counted per atom based on donor type.
/// - Huckel's 4n+2 rule is applied to the total.
pub fn is_aromatic(atom_idx: usize, mol: &Molecule) -> bool {
    let rings = find_rings(mol);

    for ring in &rings {
        if !ring.contains(&atom_idx) {
            continue;
        }

        // All atoms in the ring must be aromatic candidates
        if !ring.iter().all(|&a| is_aromatic_candidate(a, ring, mol)) {
            continue;
        }

        let mut pi_electrons = 0i32;
        for &atom in ring {
            pi_electrons += count_pi_electrons(atom, ring, mol);
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

/// Find all 1-4 atom pairs (separated by exactly 3 bonds).
/// Used for scaled VDW and electrostatic interactions.
/// Includes pairs across non-rotatable bonds (double/triple bonds)
/// which are excluded from find_torsions().
pub fn find_one_four_pairs(mol: &Molecule) -> Vec<(usize, usize)> {
    let mut pairs = Vec::new();
    let mut seen = std::collections::HashSet::new();

    // Use every bond as a potential central bond in a 3-bond path
    for bond in &mol.bonds {
        let i = bond.atom1;
        let j = bond.atom2;
        let neighbors_i = get_neighbors(i, mol);
        let neighbors_j = get_neighbors(j, mol);

        for &k in neighbors_i {
            if k == j {
                continue;
            }
            for &l in neighbors_j {
                if l == i || l == k {
                    continue;
                }
                let (a, b) = (k.min(l), k.max(l));
                if seen.insert((a, b)) {
                    pairs.push((a, b));
                }
            }
        }
    }

    pairs
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
                ..Default::default()
            },
            Bond {
                atom1: 1,
                atom2: 2,
                bond_type: BondType::Single,
                ..Default::default()
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
                    ..Default::default()
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
                    ..Default::default()
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
                ..Default::default()
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
                ..Default::default()
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
                ..Default::default()
            },
            Bond {
                atom1: 0,
                atom2: 2,
                bond_type: BondType::Single,
                ..Default::default()
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
                    ..Default::default()
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
                    ..Default::default()
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
                ..Default::default()
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
                ..Default::default()
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

    #[test]
    fn test_is_aromatic_furan_kekule() {
        // Furan Kekule form: explicit single/double bonds
        // O at position 4, double bonds at 0-1 and 2-3
        let atoms: Vec<Atom> = vec![
            make_atom("C", 6, 0),
            make_atom("C", 6, 1),
            make_atom("C", 6, 2),
            make_atom("C", 6, 3),
            make_atom("O", 8, 4),
        ];
        let mut adjacency = vec![vec![]; 5];
        let mut bonds = Vec::new();

        for (i, j, bt) in [
            (0, 1, BondType::Double),
            (1, 2, BondType::Single),
            (2, 3, BondType::Double),
            (3, 4, BondType::Single),
            (4, 0, BondType::Single),
        ] {
            adjacency[i].push(j);
            adjacency[j].push(i);
            bonds.push(Bond {
                atom1: i,
                atom2: j,
                bond_type: bt,
                ..Default::default()
            });
        }

        let mol = Molecule {
            atoms,
            bonds,
            name: "Furan".to_string(),
            adjacency,
        };

        let aromatic = get_aromatic_atoms(&mol);
        assert_eq!(aromatic.len(), 5, "Furan should have 5 aromatic atoms");
    }

    #[test]
    fn test_is_aromatic_imidazole_kekule() {
        // Imidazole Kekule form: N0(pyrrole-like), C1=N2, N2-C3, C3=C4, C4-N0
        let mut atoms: Vec<Atom> = vec![
            make_atom("N", 7, 0),
            make_atom("C", 6, 1),
            make_atom("N", 7, 2),
            make_atom("C", 6, 3),
            make_atom("C", 6, 4),
        ];
        let mut adjacency = vec![vec![]; 6];
        let mut bonds = Vec::new();

        for (i, j, bt) in [
            (0, 1, BondType::Single),
            (1, 2, BondType::Double),
            (2, 3, BondType::Single),
            (3, 4, BondType::Double),
            (4, 0, BondType::Single),
        ] {
            adjacency[i].push(j);
            adjacency[j].push(i);
            bonds.push(Bond {
                atom1: i,
                atom2: j,
                bond_type: bt,
                ..Default::default()
            });
        }
        // Add explicit H on pyrrole-like N (atom 0)
        atoms.push(make_atom("H", 1, 5));
        adjacency[0].push(5);
        adjacency[5].push(0);
        bonds.push(Bond {
            atom1: 0,
            atom2: 5,
            bond_type: BondType::Single,
            ..Default::default()
        });

        let mol = Molecule {
            atoms,
            bonds,
            name: "Imidazole".to_string(),
            adjacency,
        };

        let aromatic = get_aromatic_atoms(&mol);
        assert_eq!(aromatic.len(), 5, "Imidazole should have 5 aromatic atoms");
    }

    #[test]
    fn test_is_aromatic_2_5_dihydrofuran() {
        // 2,5-dihydrofuran: 5-membered with O and 1 double bond — NOT aromatic
        let atoms: Vec<Atom> = vec![
            make_atom("C", 6, 0),
            make_atom("C", 6, 1),
            make_atom("C", 6, 2),
            make_atom("O", 8, 3),
            make_atom("C", 6, 4),
        ];
        let mut adjacency = vec![vec![]; 5];
        let mut bonds = Vec::new();

        for (i, j, bt) in [
            (0, 1, BondType::Double),
            (1, 2, BondType::Single),
            (2, 3, BondType::Single),
            (3, 4, BondType::Single),
            (4, 0, BondType::Single),
        ] {
            adjacency[i].push(j);
            adjacency[j].push(i);
            bonds.push(Bond {
                atom1: i,
                atom2: j,
                bond_type: bt,
                ..Default::default()
            });
        }

        let mol = Molecule {
            atoms,
            bonds,
            name: "2,5-dihydrofuran".to_string(),
            adjacency,
        };

        let aromatic = get_aromatic_atoms(&mol);
        assert!(aromatic.is_empty(), "2,5-dihydrofuran should not be aromatic");
    }


}
