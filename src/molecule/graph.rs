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

    let pi_bonds: f64 = mol
        .bonds
        .iter()
        .filter(|b| b.atom1 == atom_idx || b.atom2 == atom_idx)
        .map(|b| match b.bond_type {
            BondType::Single => 0.0,
            BondType::Double => 1.0,
            BondType::Triple => 2.0,
            BondType::Aromatic => 0.5,
        })
        .sum();
    // For the Sp1 (linear) rule, count only REAL multiple bonds (double/triple):
    // aromatic ring bonds contribute 0.5 each, so a carbonyl C in an aromatic
    // ring (C=O + 2 aromatic bonds = pi 2.0) would be wrongly typed Sp1 (linear)
    // once the ring is detected aromatic (e.g. xanthine's pyrimidinedione ring:
    // C1 typed Sp1 -> broken geometry, -34 kcal). A linear center needs two
    // genuine multiple bonds (C≡N, C≡C, C=C=C).
    let pi_multi: f64 = mol
        .bonds
        .iter()
        .filter(|b| b.atom1 == atom_idx || b.atom2 == atom_idx)
        .map(|b| match b.bond_type {
            BondType::Double => 1.0,
            BondType::Triple => 2.0,
            _ => 0.0,
        })
        .sum();

    // Hypervalent S (sulfone/sulfonate/sulfonamide...): the generic pi_bonds ->
    // Sp1/Sp2 rule is for C/N (linear/planar). S(=O)2 with >=3 neighbors is
    // tetrahedral (RDKit embeds sulfone S at ~109.5 deg, not linear/sp2). Kept
    // P out: WebMM's pipeline empirically embeds P(=O) compounds better with sp2.
    let hypervalent_s = symbol == "S" && num_bonds >= 3;
    if !hypervalent_s && pi_multi >= 2.0 {
        return Hybridization::Sp1;
    }
    if !hypervalent_s && pi_bonds >= 1.0 {
        return Hybridization::Sp2;
    }

    // Special case: N with 3 single bonds adjacent to aromatic ring (aniline, amide)
    // is sp2 due to conjugation with aromatic system
    if symbol == "N" && num_bonds == 3 && pi_bonds == 0.0 {
        let has_aromatic_neighbor = neighbors.iter().any(|&n| is_aromatic(n, mol));
        if has_aromatic_neighbor {
            return Hybridization::Sp2;
        }
        // N with 3 single bonds adjacent to a C with a double bond (C=O, C=N,
        // C=C) is sp2 (amide, amidine, enamine, guanidinium) — RDKit treats these
        // as planar. Excludes S=O (sulfinamide/sulfonamide) where N stays sp3.
        let neighbor_c_has_double = neighbors.iter().any(|&n| {
            mol.atoms[n].atomic_number == 6
                && mol.bonds.iter().any(|b| {
                    (b.atom1 == n || b.atom2 == n)
                        && matches!(b.bond_type, BondType::Double | BondType::Aromatic)
                })
        });
        if neighbor_c_has_double {
            return Hybridization::Sp2;
        }
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
            4 => Hybridization::Sp3,
            5 => Hybridization::Sp3D,
            _ if num_bonds > 5 => Hybridization::Sp3D2,
            _ => Hybridization::Sp3,
        },
        "S" => match num_bonds {
            // By the time we reach here pi_bonds == 0 (single bonds only):
            // 2-coordinate S (thiol R-S-H, thioether R-S-R) is bent sp3, not sp2.
            2 => Hybridization::Sp3,
            3 => Hybridization::Sp3,
            4 => Hybridization::Sp3D,
            5 => Hybridization::Sp3D,
            6 => Hybridization::Sp3D2,
            _ if num_bonds > 6 => Hybridization::Sp3D2,
            _ => Hybridization::Sp3,
        },
        _ => Hybridization::Sp3,
    }
}

/// Find all aromatic atoms in the molecule
pub fn get_aromatic_atoms(mol: &Molecule) -> HashSet<usize> {
    mmff_aromatic_atoms(mol)
        .into_iter()
        .enumerate()
        .filter(|(_, b)| *b)
        .map(|(i, _)| i)
        .collect()
}

/// Convert bonds in aromatic rings to BondType::Aromatic
/// This matches RDKit behavior where Kekulé structures are perceived as aromatic.
/// Iterated to a fixpoint because converting one ring's bonds to Aromatic can
/// change the pi-electron count enough for an adjacent fused ring to then pass
/// the aromaticity test (e.g. indole's pyrrole ring after the benzene ring is
/// perceived).
pub fn perceive_aromatic_bonds(mol: &mut Molecule) {
    loop {
        let aromatic_atoms = get_aromatic_atoms(mol);
        let rings = find_rings(mol);
        let mut changed = false;
        for ring in &rings {
            let ring_set: HashSet<usize> = ring.iter().copied().collect();
            if !ring.iter().all(|a| aromatic_atoms.contains(a)) {
                continue;
            }
            for bond in mol.bonds.iter_mut() {
                if ring_set.contains(&bond.atom1)
                    && ring_set.contains(&bond.atom2)
                    && bond.bond_type != BondType::Aromatic
                {
                    if bond.kekule_type.is_none() {
                        bond.kekule_type = Some(bond.bond_type);
                    }
                    bond.bond_type = BondType::Aromatic;
                    changed = true;
                }
            }
        }
        if !changed {
            break;
        }
    }
}

/// RDKit `MolOps::setMMFFAromaticity` (Aromaticity.cpp) — the SDM
/// ("simple delocalization model") used by MMFF atom typing.
/// `MMFFGetMoleculeProperties` re-perceives aromaticity with this model
/// internally, regardless of the molecule's incoming aromatic flags, so our
/// MMFF types (and therefore energies) must follow it — not the default
/// RDKit aromaticity model.
///
/// Per ring (iterated to a fixpoint so fused systems can propagate):
/// - +2 pi electrons per DOUBLE ring bond (kekulized view);
/// - a ring C, or N with total bond order 4, may bring contributions from
///   exocyclic multiple bonds: +1 if the exocyclic double partner is already
///   aromatic (fused-ring propagation), otherwise it sets `exoDouble`
///   (which blocks the lone-pair bonus);
/// - +2 if N/O/divalent-S is present, there is no exocyclic double, and the
///   ring size is odd (pyrrole-type lone-pair donation);
/// - every ring C/N must be sp2-like;
/// - aromatic iff pi_e > 2 and (pi_e - 2) % 4 == 0.
///
/// Atoms of rings that fail the rule are never marked, but "perceived" atoms
/// still unlock deferred rings in later passes (mirroring RDKit's
/// `aromBitVect` vs aromatic flags distinction).
/// Kekulize aromatic bonds that carry no recorded pre-aromatization order
/// (hand-built molecules, MOL2/aromatic-form molblocks). Computes an
/// alternating single/double assignment over each aromatic component so the
/// SDM electron counting sees a valid Kekule structure. Mirrors the spirit of
/// RDKit's `Kekulize` before `setMMFFAromaticity`. Bonds that were aromatized
/// by `perceive_aromatic_bonds` keep their recorded order and are skipped.
fn kekulize_unrecorded(mol: &Molecule) -> std::collections::HashMap<(usize, usize), i32> {
    use std::collections::HashMap;
    let mut assign: HashMap<(usize, usize), i32> = HashMap::new();
    let needs: Vec<(usize, usize)> = mol
        .bonds
        .iter()
        .filter(|b| b.bond_type == BondType::Aromatic && b.kekule_type.is_none())
        .map(|b| (b.atom1.min(b.atom2), b.atom1.max(b.atom2)))
        .collect();
    if needs.is_empty() {
        return assign;
    }

    let n = mol.atoms.len();
    // aromatic bond adjacency (indices into `needs`)
    let mut nbrs: Vec<Vec<usize>> = vec![vec![]; n];
    for (bi, &(a, b)) in needs.iter().enumerate() {
        nbrs[a].push(bi);
        nbrs[b].push(bi);
    }

    // Doubles each atom needs among its aromatic bonds (RDKit Kekulize
    // semantics): aromatic carbons want exactly one double unless they
    // already carry an explicit one; ring N wants one only when divalent
    // (pyridine-type); O and divalent S contribute lone pairs, not doubles.
    let want: Vec<i32> = (0..n)
        .map(|a| {
            let z = mol.atoms[a].atomic_number as i32;
            let explicit_doubles = mol
                .bonds
                .iter()
                .filter(|b| {
                    (b.atom1 == a || b.atom2 == a)
                        && b.bond_type != BondType::Aromatic
                        && matches!(b.bond_type, BondType::Double | BondType::Triple)
                })
                .count() as i32;
            match z {
                6 => (1 - explicit_doubles).max(0),
                7 if mol.adjacency[a].len() <= 2 && mol.atoms[a].charge == 0.0 => {
                    (1 - explicit_doubles).max(0)
                }
                _ => 0,
            }
        })
        .collect();

    // backtracking over unassigned aromatic bonds (components are tiny)
    let mut chosen: Vec<Option<bool>> = vec![None; needs.len()]; // true = double
    let mut got: Vec<i32> = vec![0; n];
    let mut budget: i64 = 200_000;
    fn solve(
        from: usize,
        needs: &[(usize, usize)],
        want: &[i32],
        got: &mut [i32],
        chosen: &mut [Option<bool>],
        budget: &mut i64,
    ) -> bool {
        if *budget <= 0 {
            return false;
        }
        *budget -= 1;
        let bi = match (from..needs.len()).find(|&i| chosen[i].is_none()) {
            Some(bi) => bi,
            None => return true,
        };
        for &is_double in &[true, false] {
            let (a, b) = needs[bi];
            if is_double && (got[a] + 1 > want[a] || got[b] + 1 > want[b]) {
                continue;
            }
            chosen[bi] = Some(is_double);
            if is_double {
                got[a] += 1;
                got[b] += 1;
            }
            if solve(bi + 1, needs, want, got, chosen, budget) {
                return true;
            }
            if is_double {
                got[a] -= 1;
                got[b] -= 1;
            }
            chosen[bi] = None;
        }
        false
    }
    if !solve(0, &needs, &want, &mut got, &mut chosen, &mut budget) {
        return assign; // give up: leave unrecorded (treated as double)
    }
    for (bi, &(a, b)) in needs.iter().enumerate() {
        assign.insert((a, b), if chosen[bi] == Some(true) { 2 } else { 1 });
    }
    assign
}

pub fn mmff_aromatic_atoms(mol: &Molecule) -> Vec<bool> {
    let rings = find_rings(mol);
    let n = mol.atoms.len();
    let mut aromatic = vec![false; n]; // final aromatic marking
    if rings.is_empty() {
        return aromatic;
    }

    // kekulized view of a bond order for the SDM counting: aromatized ring
    // bonds count with their recorded pre-aromatization (Kekule) order —
    // mirrors RDKit kekulizing the molecule before setMMFFAromaticity.
    // Aromatic bonds without a record are kekulized on the fly.
    let kek = kekulize_unrecorded(mol);
    let order = |b: &crate::molecule::Bond| -> i32 {
        let key = (b.atom1.min(b.atom2), b.atom1.max(b.atom2));
        if b.bond_type == BondType::Aromatic {
            if let Some(o) = kek.get(&key) {
                return *o;
            }
        }
        match b.kekule_type.unwrap_or(b.bond_type) {
            BondType::Single => 1,
            BondType::Double => 2,
            BondType::Triple => 3,
            BondType::Aromatic => 2, // defensive: kekulizer gave up
        }
    };
    let explicit_valence = |a: usize| -> i32 {
        mol.bonds
            .iter()
            .filter(|b| b.atom1 == a || b.atom2 == a)
            .map(&order)
            .sum()
    };
    let implicit_h = |a: usize| -> i32 {
        let dv = default_valence(mol.atoms[a].atomic_number as i32);
        if dv <= 0 {
            return 0;
        }
        (dv - explicit_valence(a)).max(0)
    };
    // sp2-like for the canBeAromatic check: incident multiple bond, hetero
    // pi-donor, or already marked aromatic
    let sp2_like = |a: usize, aromatic: &[bool]| -> bool {
        if aromatic[a] {
            return true;
        }
        let z = mol.atoms[a].atomic_number as i32;
        if z == 7 || z == 8 || z == 16 {
            return true;
        }
        mol.bonds
            .iter()
            .any(|b| (b.atom1 == a || b.atom2 == a) && order(b) >= 2)
    };

    let ring_members: Vec<std::collections::HashSet<usize>> =
        rings.iter().map(|r| r.iter().copied().collect()).collect();

    let mut perceived = vec![false; n]; // RDKit aromBitVect
    let mut old_n: i64 = -1;
    loop {
        for (ri, ring) in rings.iter().enumerate() {
            let len = ring.len();
            let mut pi_e: i32 = 0;
            let mut is_nos = false;
            let mut exo_double = false;
            let mut defer = false;
            for j in 0..len {
                let a = ring[j];
                let next = ring[(j + 1) % len];
                let z = mol.atoms[a].atomic_number as i32;
                if z == 7 || z == 8 || (z == 16 && mol.adjacency[a].len() == 2) {
                    is_nos = true;
                }
                let bond_to_next = mol.bonds.iter().find(|b| {
                    (b.atom1 == a && b.atom2 == next) || (b.atom2 == a && b.atom1 == next)
                });
                if bond_to_next.map(|b| order(b) == 2).unwrap_or(false) {
                    pi_e += 2;
                    continue;
                }
                // only C, or N with total bond order 4, can bring exocyclic
                // pi contributions
                let n_bo4 = z == 7 && (explicit_valence(a) + implicit_h(a)) == 4;
                if z != 6 && !n_bo4 {
                    continue;
                }
                for &nbr in &mol.adjacency[a] {
                    if ring_members[ri].contains(&nbr) {
                        continue; // only exocyclic neighbours
                    }
                    let b = mol.bonds.iter().find(|b| {
                        (b.atom1 == a && b.atom2 == nbr) || (b.atom2 == a && b.atom1 == nbr)
                    });
                    if b.map(|b| order(b) == 1).unwrap_or(true) {
                        continue;
                    }
                    // neighbour in an unprocessed ring: defer this ring
                    let nbr_in_ring = ring_members.iter().any(|rm| rm.contains(&nbr));
                    if nbr_in_ring && !perceived[nbr] {
                        defer = true;
                        break;
                    }
                    if b.map(|b| order(b) == 2).unwrap_or(false) {
                        if aromatic[nbr] {
                            pi_e += 1;
                        } else {
                            exo_double = true;
                        }
                    }
                }
                if defer {
                    break;
                }
            }
            if defer {
                continue;
            }
            let mut can_be_aromatic = true;
            for &a in ring.iter() {
                perceived[a] = true;
                let z = mol.atoms[a].atomic_number as i32;
                if (z == 6 || z == 7) && !sp2_like(a, &aromatic) {
                    can_be_aromatic = false;
                }
            }
            if !can_be_aromatic {
                continue;
            }
            if is_nos && !exo_double && len % 2 == 1 {
                pi_e += 2;
            }
            if pi_e > 2 && (pi_e - 2) % 4 == 0 {
                for &a in ring.iter() {
                    aromatic[a] = true;
                }
            }
        }
        // RDKit termination: stop when no new atoms were perceived
        let n_set: i64 = perceived.iter().filter(|p| **p).count() as i64;
        let all_done = rings.iter().flatten().all(|&a| perceived[a]);
        if all_done || n_set <= old_n {
            break;
        }
        old_n = n_set;
    }
    aromatic
}

/// Default valence for the organic subset (SDM N-bond-order check).
fn default_valence(z: i32) -> i32 {
    match z {
        5 => 3,
        6 => 4,
        7 => 3,
        8 => 2,
        14 => 4,
        15 => 3,
        16 => 2,
        33 => 3,
        34 => 2,
        _ => -1,
    }
}

/// Check if atom is in an aromatic ring (ring membership + Huckel rule).
/// Uses RDKit-style aromaticity perception:
/// - All atoms in the ring must be aromatic candidates.
/// - Pi electrons are counted per atom based on donor type.
/// - Huckel's 4n+2 rule is applied to the total.
pub fn is_aromatic(atom_idx: usize, mol: &Molecule) -> bool {
    mmff_aromatic_atoms(mol)[atom_idx]
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

/// Check if a bond between two atoms is in a ring
pub fn is_in_ring(atom1: usize, atom2: usize, mol: &Molecule) -> bool {
    let rings = find_rings(mol);
    for ring in &rings {
        let has_a = ring.contains(&atom1);
        let has_b = ring.contains(&atom2);
        if has_a && has_b && mol.adjacency[atom1].contains(&atom2) {
            return true;
        }
    }
    false
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
    Sp3D,
    Sp3D2,
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

    for j in 0..mol.atoms.len() {
        let neighbors: Vec<usize> = get_neighbors(j, mol).to_vec();
        if neighbors.len() < 2 {
            continue;
        }
        for idx1 in 0..neighbors.len() {
            for idx2 in (idx1 + 1)..neighbors.len() {
                let i = neighbors[idx1];
                let k = neighbors[idx2];
                angles.push(Angle {
                    atom1: i,
                    atom2: j,
                    atom3: k,
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

    for bond in &mol.bonds {
        let i = bond.atom1;
        let j = bond.atom2;
        let neighbors_i = get_neighbors(i, mol);
        let neighbors_j = get_neighbors(j, mol);

        for &k in neighbors_i {
            if k != j {
                for &l in neighbors_j {
                    if l != i && l != k {
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

    let aromatic_set = mmff_aromatic_atoms(mol);
    for (atom_idx, neighbors) in mol.adjacency.iter().enumerate() {
        let neighbors: Vec<usize> = neighbors.to_vec();

        // Only atoms with 3+ neighbors can have out-of-plane bending
        if neighbors.len() >= 3 {
            // Only sp2 and aromatic atoms typically have significant OOP
            let hybrid = determine_hybridization(atom_idx, mol);
            if hybrid == Hybridization::Sp2 || aromatic_set[atom_idx] {
                // RDKit creates 3 OOP terms per 3-neighbor combination: for each
                // choice of which neighbor is the "out-of-plane" atom (atom1),
                // the other two define the reference plane. This gives different
                // OOP angles and must all be summed.
                for i in 0..neighbors.len() {
                    for j in (i + 1)..neighbors.len() {
                        for k in (j + 1)..neighbors.len() {
                            // 3 cyclic permutations: each neighbor takes a turn
                            // as atom1 (the out-of-plane atom)
                            oops.push(OutOfPlane {
                                central: atom_idx,
                                atom1: neighbors[i],
                                atom2: neighbors[j],
                                atom3: neighbors[k],
                            });
                            oops.push(OutOfPlane {
                                central: atom_idx,
                                atom1: neighbors[j],
                                atom2: neighbors[i],
                                atom3: neighbors[k],
                            });
                            oops.push(OutOfPlane {
                                central: atom_idx,
                                atom1: neighbors[k],
                                atom2: neighbors[i],
                                atom3: neighbors[j],
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
    // --- MMFF SDM aromaticity regression tests (RDKit parity) ---

    fn mol_from_block(block: &str) -> Molecule {
        crate::molecule::parser::parse_sdf(block).unwrap()
    }

    const MALEIMIDE: &str = "maleimide\n     RDKit          3D\n\n 10 10  0  0  0  0  0  0  0  0999 V2000\n   -1.1905   -2.0775   -0.0430 O   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.5649   -1.0363   -0.0214 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.7738   -0.8719   -0.0146 N   0  0  0  0  0  0  0  0  0  0  0  0\n    1.0963    0.4376    0.0113 C   0  0  0  0  0  0  0  0  0  0  0  0\n    2.2047    0.9348    0.0238 O   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.1629    1.1890    0.0223 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.1614    0.3031    0.0027 C   0  0  0  0  0  0  0  0  0  0  0  0  0\n    1.4414   -1.6240   -0.0272 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.2148    2.2629    0.0427 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -2.2217    0.4823    0.0033 H   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  2  0\n  2  3  1  0\n  3  4  1  0\n  4  5  2  0\n  4  6  1  0\n  6  7  2  0\n  7  2  1  0\n  3  8  1  0\n  6  9  1  0\n  7 10  1  0\nM  END\n";

    #[test]
    fn mmff_aromaticity_maleimide_not_aromatic() {
        // The imide five-ring has only 4 in-cycle pi electrons (the two C=O
        // are exocyclic): RDKit's SDM demotes it and types the ring
        // C_2/N_AM/C_VIN. The old Huckel donor model wrongly aromatized it
        // (60 kcal/mol single-point error vs RDKit).
        let mol = mol_from_block(MALEIMIDE);
        let arom = mmff_aromatic_atoms(&mol);
        assert!(
            arom.iter().all(|b| !b),
            "maleimide ring must not be aromatic"
        );
    }

    const THIOUREA: &str = "thiourea\n     RDKit          3D\n\n  8  7  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n    1.3354    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    2.1018    1.1859    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n    1.8715   -1.5277    0.0000 S   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.4847    0.9344    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.4847   -0.9344    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    3.1118    1.0145    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    1.6936    2.1634    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0\n  2  3  1  0\n  2  4  2  0\n  1  5  1  0\n  1  6  1  0\n  3  7  1  0\n  3  8  1  0\nM  END\n";

    #[test]
    fn thioamide_typing_matches_rdkit() {
        // RDKit: N -> 10 (N_AM), C -> 3 (C_2), S -> 16 (S_2), N-H -> 28 (H_NAM).
        // The C(=S) partner must count as acyl for both the N and its H's.
        let mol = mol_from_block(THIOUREA);
        let ff = crate::mmff::MMFFForceField::new(&mol, crate::mmff::MMFFVariant::MMFF94s);
        let t = &ff.atom_types;
        use crate::mmff::MMFFAtomType::*;
        assert_eq!(t[0], N_AM);
        assert_eq!(t[1], C_2);
        assert_eq!(t[2], N_AM);
        assert_eq!(t[3], S_2);
        assert_eq!(t[4], H_NAM);
        assert_eq!(t[6], H_NAM);
    }

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
                stereo_parity: 0,
            },
            Atom {
                symbol: "C".to_string(),
                atomic_number: 6,
                mass: 12.011,
                charge: 0.0,
                position: [0.0; 3],
                index: 1,
                stereo_parity: 0,
            },
            Atom {
                symbol: "H".to_string(),
                atomic_number: 1,
                mass: 1.008,
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
            stereo_parity: 0,
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
        assert!(
            aromatic.is_empty(),
            "2,5-dihydrofuran should not be aromatic"
        );
    }
}
