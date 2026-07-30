//! MMFF94/MMFF94s force field implementation

use std::collections::{HashMap, HashSet};

use crate::molecule::Bond;
use crate::molecule::BondType;
use crate::molecule::Hybridization;
use crate::molecule::Molecule;

pub mod angle;
pub mod atom_types;
pub mod bond;
pub mod charges;
pub mod electrostatics;
pub mod estimation;
pub mod mmff_tables;
pub mod oop;
pub mod params;
pub mod stretch_bend;
pub mod torsion;
pub mod vdw;

pub use angle::*;
pub use atom_types::*;
pub use bond::*;
pub use charges::*;
pub use electrostatics::*;
pub use estimation::*;
pub use oop::*;
pub use params::*;
pub use stretch_bend::*;
pub use torsion::*;
pub use vdw::*;

pub use crate::MMFFVariant;

pub fn base_type(t: MMFFAtomType) -> MMFFAtomType {
    match t {
        MMFFAtomType::H_OH
        | MMFFAtomType::H_ONC
        | MMFFAtomType::H_COOH
        | MMFFAtomType::H_OAR
        | MMFFAtomType::H_N3
        | MMFFAtomType::H_NAM
        | MMFFAtomType::H_NIM
        | MMFFAtomType::HS => MMFFAtomType::H,
        MMFFAtomType::CR4R => MMFFAtomType::C_3,
        MMFFAtomType::CE4R => MMFFAtomType::C_2,
        MMFFAtomType::CR3R => MMFFAtomType::C_3,
        MMFFAtomType::HNRP => MMFFAtomType::H,
        MMFFAtomType::S2CM => MMFFAtomType::S_3,
        MMFFAtomType::S_O3 => MMFFAtomType::S_3,
        MMFFAtomType::S_CSO => MMFFAtomType::S_2,
        MMFFAtomType::P_ARM => MMFFAtomType::P_3,
        MMFFAtomType::HOS => MMFFAtomType::H,
        MMFFAtomType::H_OXP => MMFFAtomType::H,
        // 5-membered heteroaromatic ring atoms use same params as generic aromatic
        MMFFAtomType::C5A | MMFFAtomType::C5B => MMFFAtomType::C_AR,
        MMFFAtomType::C5A_M => MMFFAtomType::C_AR,
        MMFFAtomType::NPYL | MMFFAtomType::N5A | MMFFAtomType::N5B => MMFFAtomType::N_AR,
        MMFFAtomType::NPYL_M => MMFFAtomType::N_AR,
        MMFFAtomType::N_PYR => MMFFAtomType::N_AR,
        MMFFAtomType::N_T3 => MMFFAtomType::N_4,
        MMFFAtomType::N_POX2 => MMFFAtomType::N_POX,
        MMFFAtomType::N_RAD => MMFFAtomType::N_3,
        MMFFAtomType::OFUR => MMFFAtomType::O_3,
        MMFFAtomType::O_R => MMFFAtomType::O_3,
        MMFFAtomType::O_3_Z => MMFFAtomType::O_3,
        MMFFAtomType::O_3P => MMFFAtomType::O_3,
        MMFFAtomType::O_2P => MMFFAtomType::O_2,
        MMFFAtomType::N_NITROSO => MMFFAtomType::N_2,
        MMFFAtomType::CL4 => MMFFAtomType::Cl,
        MMFFAtomType::N_2Z | MMFFAtomType::N_1M => MMFFAtomType::N_1,
        MMFFAtomType::CID => MMFFAtomType::C_1,
        MMFFAtomType::NID => MMFFAtomType::N_1,
        MMFFAtomType::NCN_PLUS => MMFFAtomType::N_AM,
        MMFFAtomType::OXIDE => MMFFAtomType::O_3,
        MMFFAtomType::OH2 => MMFFAtomType::O_3,
        // SP3D/SP3D2 types fall back to base sp3 types for parameters
        MMFFAtomType::P_3D => MMFFAtomType::P_3,
        MMFFAtomType::S_3D | MMFFAtomType::S_3D2 => MMFFAtomType::S_3,
        other => other,
    }
}

/// MMFF atom type
#[allow(non_camel_case_types)]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum MMFFAtomType {
    // Hydrogen
    H,
    H_OH,   // H in water (H-O-H, all O neighbors are H)
    H_ONC,  // H bonded to O_3 where O bonded to C_3 (ethanol)
    H_COOH, // H bonded to O_3 where O bonded to C_2 (carboxylic acid -OH)
    H_OAR,  // H bonded to O_3 where O bonded to C_AR (phenol)
    H_N3,   // H bonded to N_3 (ammonia)
    H_NAM,  // H bonded to N_AM or N_AR (aniline, acetamide)
    H_NIM,  // H bonded to sp2 imine N with C=N double bond (MMFF 27)
    HS,     // H bonded to sulfur (thiol R-S-H) (MMFF 71)
    CR4R,   // sp3 carbon in a 4-membered ring (MMFF 20)
    CE4R,   // sp2 carbon in a 4-membered ring (MMFF 30)
    CR3R,   // sp3 carbon in a 3-membered ring (MMFF 22)
    HNRP,   // H on a positively-charged (ammonium) N (MMFF 36)
    S2CM,   // S in C(=S)=S / thiocarbonyl (MMFF 72)
    HOS,    // H on O bonded to S (sulfonic acid) (MMFF 33)
    H_OXP,  // H on oxonium O+ (MMFF 50)

    // Carbons
    C_3,
    C_2,
    C_VIN, // Generic sp2 alkene carbon, double bond(s) only to C (MMFF 2)
    C_CO2, // Carboxylate carbon, sp2 C with C=O + C-[O-] (MMFF 41)
    C_1,
    C_AR,
    C5A, // Alpha C in 5-membered heteroaromatic ring (MMFF 63)
    C5B, // Beta C in 5-membered heteroaromatic ring (MMFF 64)
    C5A_M, // Alpha C in pyrrolide anion ring (MMFF 78)
    C_CAT,
    C_AN,
    CID,   // Isonitrile carbon C- (MMFF 60)
    NID,   // Isonitrile/nitrile-oxide nitrogen N+ (MMFF 61)
    NCN_PLUS, // Guanidinium/amidinium N (MMFF 55)
    OXIDE,  // Nitrile oxide O- (MMFF 35)

    // Nitrogens
    N_3,
    N_2,
    N_1,
    N_AR, // Pyridine-type N in 6-membered ring (MMFF 38)
    NPYL, // Pyrrole-type N in 5-membered ring (MMFF 39)
    N_PL3,
    N_AM,
    N_4,
    N_2Z,
    N_1M,  // Terminal anionic N in azide/diazo (MMFF 47)
    N_SOM,
    N_NO2, // Nitro group nitrogen (MMFF 45)
    N_SO2, // Sulfonamide nitrogen, N bonded to SO2 sulfur (MMFF 43)
    N_NITROSO, // Nitroso nitrogen N=O (MMFF 46)
    N5A, // Alpha N in 5-membered heteroaromatic ring (MMFF 65)
    N5B, // Beta N in 5-membered heteroaromatic ring (MMFF 66)
    N_POX, // Pyridine N-oxide nitrogen (MMFF 69)
    N_RAD, // Nitrene / radical N (MMFF 62)
    NPYL_M, // Pyrrolide anion N (MMFF 76)
    N_PYR, // N-methylpyridinium N (MMFF 58)
    N_T3, // Trimethylamine N-oxide N (MMFF 68)
    N_POX2, // Pyridine N-oxide N (MMFF 69)

    // Oxygens
    O_3,
    O_2,
    O_R,  // Ether/alcohol O (MMFF 6) — note: water uses OH2 (70)
    OH2,  // Water oxygen (MMFF 70)
    OFUR, // Furan oxygen (MMFF 59)
    O_CO2,
    O_3_Z,
    O_3P, // Oxonium oxygen O+ (MMFF 49)
    O_2P, // Aromatic oxonium O+ in ring (MMFF 51)

    // Halogens
    F,
    Cl,
    CL4, // Perchlorate / hypervalent Cl (MMFF 77)
    Br,
    I,
    F_M,  // Fluoride anion (MMFF 89)
    CL_M, // Chloride anion (MMFF 90)
    BR_M, // Bromide anion (MMFF 91)

    // Sulfur, Phosphorus
    S_3,
    S_2,
    S_AR,
    S_OX, // Sulfoxide sulfur, one double bond to O (MMFF 17)
    S_O2, // Sulfone/sulfonamide sulfur, two double bonds to O (MMFF 18)
    S_O3, // Sulfite sulfur, 3-coord with 1 S=O (MMFF 73)
    S_CSO, // Sulfene sulfur, C=S=O (MMFF 74)
    P_ARM, // P in aromatic 3-ring (MMFF 75)
    P_3,
    P_4,
    Si, // Silicon (MMFF 19)

    // SP3D / SP3D2 (trigonal bipyramidal / octahedral)
    P_3D,  // Phosphorus in trigonal bipyramidal geometry (e.g., PF5)
    S_3D,  // Sulfur in trigonal bipyramidal geometry
    S_3D2, // Sulfur in octahedral geometry (e.g., SF6)

    // Ions
    Fe_P2,
    Fe_P3,
    Li,
    Na,
    K,
    Zn_P2,
    Ca_P2,
    Cu_P1,
    Cu_P2,
    Mg_P2,
}

/// MMFF force field
pub struct MMFFForceField {
    pub mol: Molecule,
    pub atom_types: Vec<MMFFAtomType>,
    pub type_ids: Vec<u8>,
    pub charges: Vec<f64>,
    pub variant: MMFFVariant,
    pub angles: Vec<crate::molecule::Angle>,
    pub torsions: Vec<crate::molecule::Torsion>,
    pub torsion_types: Vec<(u8, u8)>,
    pub oops: Vec<crate::molecule::OutOfPlane>,
    pub excluded_pairs: HashSet<(usize, usize)>,
    pub one_four_pairs: HashSet<(usize, usize)>,
    pub bond_map: HashMap<(usize, usize), Bond>,
    pub rings: Vec<Vec<usize>>,
}

impl MMFFForceField {
    pub fn new(mol: &Molecule, variant: MMFFVariant) -> Self {
        let atom_types = Self::assign_atom_types(mol);
        let charges = Self::calculate_charges(mol, &atom_types);

        let angles = crate::molecule::graph::find_angles(mol);
        let torsions = crate::molecule::graph::find_torsions(mol);
        let oops = crate::molecule::graph::find_out_of_planes(mol);

        let mut excluded_pairs = HashSet::new();
        let mut one_four_pairs = HashSet::new();
        let mut bond_map = HashMap::new();

        // 1-2 pairs (bonded)
        for bond in &mol.bonds {
            let (a, b) = (bond.atom1.min(bond.atom2), bond.atom1.max(bond.atom2));
            excluded_pairs.insert((a, b));
            bond_map.insert((a, b), bond.clone());
        }

        // 1-3 pairs (angle endpoints)
        for angle in &angles {
            let (a, c) = (angle.atom1.min(angle.atom3), angle.atom1.max(angle.atom3));
            excluded_pairs.insert((a, c));
        }

        // 1-4 pairs (torsion endpoints) — scaled by 0.75
        for torsion in &torsions {
            let (a, d) = (
                torsion.atom1.min(torsion.atom4),
                torsion.atom1.max(torsion.atom4),
            );
            one_four_pairs.insert((a, d));
        }

        // Also find 1-4 pairs across non-rotatable bonds (double/triple bonds)
        use crate::molecule::graph::find_one_four_pairs;
        for (a, d) in find_one_four_pairs(mol) {
            one_four_pairs.insert((a, d));
        }

        // Compute torsion types using RDKit's classification
        let rings = crate::molecule::graph::find_rings(mol);
        // type_ids report each atom's actual MMFF type (matching RDKit's
        // GetMMFFAtomType), and the equivalence-level protocol (EQ_LEVELS) handles
        // parameter-lookup fallback. A previous version collapsed subtypes via
        // base_type() (C5A->C_AR, NPYL->N_AR, ...) for energy correctness, but
        // once the 5-ring N classification / charge bugs were fixed the collapse
        // became energy-neutral, so all types now keep their real IDs. (C5A vs
        // C5B alpha/beta is the one remaining cosmetic gap — both EQ-fall to C_2.)
        let type_ids: Vec<u8> = atom_types.iter().map(|&at| mmff_type_id(at)).collect();

        let torsion_types: Vec<(u8, u8)> = torsions
            .iter()
            .map(|t| {
                let i = t.atom1;
                let j = t.atom2;
                let k = t.atom3;
                let l = t.atom4;

                let bond_ij_type = bond_map
                    .get(&(i.min(j), i.max(j)))
                    .map(|b| get_mmff_bond_type(b.bond_type, type_ids[i], type_ids[j]))
                    .unwrap_or(0);
                let bond_jk_type = bond_map
                    .get(&(j.min(k), j.max(k)))
                    .map(|b| get_mmff_bond_type(b.bond_type, type_ids[j], type_ids[k]))
                    .unwrap_or(0);
                let bond_kl_type = bond_map
                    .get(&(k.min(l), k.max(l)))
                    .map(|b| get_mmff_bond_type(b.bond_type, type_ids[k], type_ids[l]))
                    .unwrap_or(0);

                let bond_jk_actual = bond_map
                    .get(&(j.min(k), j.max(k)))
                    .map(|b| b.bond_type)
                    .unwrap_or(BondType::Single);

                let torsion_set: HashSet<usize> = [i, j, k, l].iter().copied().collect();

                let ring4 = rings
                    .iter()
                    .any(|r| r.len() == 4 && torsion_set.iter().all(|a| r.contains(a)));

                let ring5 = rings
                    .iter()
                    .any(|r| r.len() == 5 && torsion_set.iter().all(|a| r.contains(a)));

                get_mmff_torsion_type(
                    bond_ij_type,
                    bond_jk_type,
                    bond_kl_type,
                    bond_jk_actual,
                    ring4,
                    ring5,
                    type_ids[i],
                    type_ids[j],
                    type_ids[k],
                    type_ids[l],
                )
            })
            .collect();

        Self {
            mol: mol.clone(),
            atom_types,
            type_ids,
            charges,
            variant,
            angles,
            torsions,
            torsion_types,
            oops,
            excluded_pairs,
            one_four_pairs,
            bond_map,
            rings,
        }
    }

    pub fn angle_ring_size(&self, i: usize, j: usize, k: usize) -> u8 {
        let angle_set: HashSet<usize> = [i, j, k].iter().copied().collect();
        for ring in &self.rings {
            if (ring.len() == 3 || ring.len() == 4) && angle_set.iter().all(|a| ring.contains(a)) {
                return ring.len() as u8;
            }
        }
        0
    }

    fn assign_atom_types(mol: &Molecule) -> Vec<MMFFAtomType> {
        use crate::molecule::graph::{determine_hybridization, find_rings, get_aromatic_atoms};

        let aromatic_atoms = get_aromatic_atoms(mol);
        let rings = find_rings(mol);

        // Atoms in 3- and 4-membered rings (carbons get CR3R/CR4R/CE4R per RDKit MMFF)
        let in_3ring: std::collections::HashSet<usize> = rings
            .iter()
            .filter(|r| r.len() == 3)
            .flat_map(|r| r.iter().copied())
            .collect();
        let in_4ring: std::collections::HashSet<usize> = rings
            .iter()
            .filter(|r| r.len() == 4)
            .flat_map(|r| r.iter().copied())
            .collect();

        // Detect 5-membered heteroaromatic rings and classify their atoms
        let mut atom_5ring_type: std::collections::HashMap<usize, MMFFAtomType> =
            std::collections::HashMap::new();
        for ring in &rings {
            if ring.len() != 5 {
                continue;
            }
            let all_aromatic = ring.iter().all(|a| aromatic_atoms.contains(a));
            if !all_aromatic {
                continue;
            }
            let hetero_atoms: Vec<usize> = ring
                .iter()
                .filter(|&&a| {
                    let an = mol.atoms[a].atomic_number;
                    an == 7 || an == 8 || an == 16
                })
                .copied()
                .collect();
            if hetero_atoms.is_empty() {
                continue;
            }

            // Two passes: first classify heteroatoms (N/O/S), then classify
            // carbons from their ring hetero-neighbors' types. (RDKit's C5A vs
            // C5B keys on neighbor type: C adjacent to a pyrrole-type N (NPYL),
            // O (OFUR), or S (S_AR), or flanked by two heteroatoms, is C5A;
            // C adjacent only to a pyridine-type N (N5B) or no heteroatom is C5B.)
            use std::collections::HashSet;
            let ringset: HashSet<usize> = ring.iter().copied().collect();
            let mut tmp: std::collections::HashMap<usize, MMFFAtomType> =
                std::collections::HashMap::new();
            // Pass 1: heteroatoms
            for &atom_idx in ring {
                match mol.atoms[atom_idx].atomic_number {
                    7 => {
                        let n_neighbors = mol.adjacency[atom_idx].len();
                        if mol.atoms[atom_idx].charge < -0.5 {
                            tmp.insert(atom_idx, MMFFAtomType::NPYL_M);
                        } else {
                            tmp.insert(
                                atom_idx,
                                if n_neighbors >= 3 { MMFFAtomType::NPYL } else { MMFFAtomType::N5B },
                            );
                        }
                    }
                    8 => {
                        tmp.insert(atom_idx, MMFFAtomType::OFUR);
                    }
                    16 => {
                        tmp.insert(atom_idx, MMFFAtomType::S_AR);
                    }
                    _ => {}
                }
            }
            // Refine: a dicoordinate (pyridine-type) N adjacent in the ring to a
            // lone-pair heteroatom (NPYL/OFUR/S_AR) is the alpha N -> N5A (65);
            // otherwise beta -> N5B (66). E.g. pyrazole's =N- next to the pyrrole N.
            for (i, t) in tmp.clone() {
                if t != MMFFAtomType::N5B {
                    continue;
                }
                let adj_lp = mol.adjacency[i].iter().any(|&nbr| {
                    ringset.contains(&nbr)
                        && matches!(tmp.get(&nbr), Some(MMFFAtomType::NPYL) | Some(MMFFAtomType::OFUR) | Some(MMFFAtomType::S_AR))
                });
                if adj_lp {
                    tmp.insert(i, MMFFAtomType::N5A);
                }
            }
            // Pass 2: carbons
            for &atom_idx in ring {
                if mol.atoms[atom_idx].atomic_number != 6 {
                    continue;
                }
                let mut lonepair_het = 0; // NPYL / OFUR / S_AR neighbors
                let mut het_total = 0;
                for &nbr in &mol.adjacency[atom_idx] {
                    if !ringset.contains(&nbr) {
                        continue;
                    }
                    match tmp.get(&nbr) {
                        Some(MMFFAtomType::NPYL) | Some(MMFFAtomType::OFUR) | Some(MMFFAtomType::S_AR) => {
                            lonepair_het += 1;
                            het_total += 1;
                        }
                        Some(MMFFAtomType::N5B) => {
                            het_total += 1;
                        }
                        _ => {}
                    }
                }
                tmp.insert(
                    atom_idx,
                    if tmp.values().any(|&t| t == MMFFAtomType::NPYL_M) {
                        MMFFAtomType::C5A_M
                    } else if lonepair_het >= 1 || het_total >= 2 {
                        MMFFAtomType::C5A
                    } else {
                        MMFFAtomType::C5B
                    },
                );
            }
            for (&atom_idx, &typ) in tmp.iter() {
                atom_5ring_type.insert(atom_idx, typ);
            }
        }

        mol.atoms
            .iter()
            .map(|atom| {
                let idx = atom.index;

                // Check if atom is in a 5-membered heteroaromatic ring
                if let Some(&typ) = atom_5ring_type.get(&idx) {
                    return typ;
                }

                let hybrid = determine_hybridization(idx, mol);
                let aromatic = aromatic_atoms.contains(&idx);
                let num_bonds = crate::molecule::graph::get_neighbors(idx, mol).len();
                let charge = atom.charge;

                let has_c_o_neighbor = mol.bonds.iter().any(|b| {
                    let other = if b.atom1 == idx {
                        b.atom2
                    } else if b.atom2 == idx {
                        b.atom1
                    } else {
                        return false;
                    };
                    mol.atoms[other].atomic_number == 6
                        && mol.bonds.iter().any(|b2| {
                            (b2.atom1 == other || b2.atom2 == other)
                                && matches!(b2.bond_type, BondType::Double)
                                && mol.atoms[if b2.atom1 == other {
                                    b2.atom2
                                } else {
                                    b2.atom1
                                }]
                                .atomic_number
                                    == 8
                        })
                });

                // This N's own double bond to C (imine =N) — excludes it from the
                // enamine/amidine N_PL3 rules below (those are for the N single-
                // bonded to the C=X carbon; the =N itself is sp2 N_2).
                let n_owns_cn = mol.bonds
                    .iter()
                    .any(|b| {
                        b.bond_type == BondType::Double
                            && (b.atom1 == idx || b.atom2 == idx)
                            && mol.atoms[if b.atom1 == idx { b.atom2 } else { b.atom1 }].atomic_number
                                == 6
                    });
                // N bonded to a carbon that has a C=C double bond (enamine /
                // vinylamine) → planar N (MMFF N_PL3, 40), like aniline. (Aromatic
                // bonds are BondType::Aromatic, not Double, so this only catches
                // non-aromatic C=C neighbors; aniline is handled separately.)
                let has_c_c_neighbor = !n_owns_cn && mol.bonds.iter().any(|b| {
                    let other = if b.atom1 == idx {
                        b.atom2
                    } else if b.atom2 == idx {
                        b.atom1
                    } else {
                        return false;
                    };
                    mol.atoms[other].atomic_number == 6
                        && mol.bonds.iter().any(|b2| {
                            (b2.atom1 == other || b2.atom2 == other)
                                && matches!(b2.bond_type, BondType::Double)
                                && mol.atoms[if b2.atom1 == other {
                                    b2.atom2
                                } else {
                                    b2.atom1
                                }]
                                .atomic_number
                                    == 6
                        })
                });

                // Amidine / guanidinium N: this N is single-bonded to a carbon that
                // is double-bonded to a *different* nitrogen (R-N-C(=N')-). Such Ns
                // are planar → MMFF 40 (N_PL3). Excludes the =N' imine itself.
                let has_c_n_neighbor = !n_owns_cn
                    && mol.bonds.iter().any(|b| {
                        if b.bond_type != BondType::Single || (b.atom1 != idx && b.atom2 != idx) {
                            return false;
                        }
                        let c = if b.atom1 == idx { b.atom2 } else { b.atom1 };
                        mol.atoms[c].atomic_number == 6
                            && mol.bonds.iter().any(|b2| {
                                b2.bond_type == BondType::Double
                                    && (b2.atom1 == c || b2.atom2 == c)
                                    && mol.atoms[if b2.atom1 == c { b2.atom2 } else { b2.atom1 }]
                                        .atomic_number
                                        == 7
                            })
                    });

                let has_double_bond_to_c = mol.bonds.iter().any(|b| {
                    let other = if b.atom1 == idx {
                        b.atom2
                    } else if b.atom2 == idx {
                        b.atom1
                    } else {
                        return false;
                    };
                    mol.atoms[other].atomic_number == 6 && b.bond_type == BondType::Double
                });

                let carbon_neighbors: Vec<usize> = mol.adjacency[idx]
                    .iter()
                    .filter(|&&n| mol.atoms[n].atomic_number == 6)
                    .copied()
                    .collect();

                let oxygen_neighbors: Vec<usize> = mol.adjacency[idx]
                    .iter()
                    .filter(|&&n| mol.atoms[n].atomic_number == 8)
                    .copied()
                    .collect();

                let h_neighbor_count = mol.adjacency[idx]
                    .iter()
                    .filter(|&&n| mol.atoms[n].atomic_number == 1)
                    .count();

                // Counts double bonds from `center` to oxygen atoms
                let count_double_o = |center: usize| -> usize {
                    mol.bonds
                        .iter()
                        .filter(|b| b.bond_type == BondType::Double)
                        .filter(|b| {
                            let other = if b.atom1 == center {
                                b.atom2
                            } else if b.atom2 == center {
                                b.atom1
                            } else {
                                return false;
                            };
                            mol.atoms[other].atomic_number == 8
                        })
                        .count()
                };

                // Atomic numbers of this atom's double-bond partners
                let double_bond_partners: Vec<u8> = mol
                    .bonds
                    .iter()
                    .filter(|b| b.bond_type == BondType::Double)
                    .filter_map(|b| {
                        if b.atom1 == idx {
                            Some(b.atom2)
                        } else if b.atom2 == idx {
                            Some(b.atom1)
                        } else {
                            None
                        }
                    })
                    .map(|o| mol.atoms[o].atomic_number)
                    .collect();

                let double_o_count = double_bond_partners.iter().filter(|&&an| an == 8).count();
                let has_double_to_nos = double_bond_partners
                    .iter()
                    .any(|&an| an == 7 || an == 8 || an == 16);

                // Bond order to a neighboring SO2 sulfur (S with >=2 double-bonded O's)
                let so2_s_bond_order: Option<u8> = mol.adjacency[idx].iter().find_map(|&n| {
                    if mol.atoms[n].atomic_number != 16 || count_double_o(n) < 2 {
                        return None;
                    }
                    mol.bonds
                        .iter()
                        .find(|b| {
                            (b.atom1 == idx && b.atom2 == n) || (b.atom1 == n && b.atom2 == idx)
                        })
                        .map(|b| match b.bond_type {
                            BondType::Single => 1u8,
                            BondType::Double => 2u8,
                            _ => 0u8,
                        })
                });
                let bonded_to_so2_s = so2_s_bond_order.is_some();

                // Nitro N: exactly 2 O neighbors, >=1 double bond to O, no H neighbors
                let is_nitro_n = atom.atomic_number == 7
                    && oxygen_neighbors.len() == 2
                    && h_neighbor_count == 0
                    && double_o_count >= 1
                    && charge > 0.5;

                // Bonded (single or double) to a nitro-group N
                let bonded_to_nitro_n = atom.atomic_number == 8
                    && mol.adjacency[idx].iter().any(|&n| {
                        mol.atoms[n].atomic_number == 7
                            && mol.atoms[n].charge > 0.5
                            && mol.adjacency[n]
                                .iter()
                                .filter(|&&m| mol.atoms[m].atomic_number == 8)
                                .count()
                                == 2
                            && count_double_o(n) >= 1
                    });

                // Carboxylate C: sp2 C with >=1 C=O double bond AND >=1 single bond to O,
                // where at least one O neighbor has formal charge < -0.5 (RDKit-verified)
                let is_carboxylate_c = atom.atomic_number == 6
                    && double_o_count >= 1
                    && mol.bonds.iter().any(|b| {
                        b.bond_type == BondType::Single
                            && (b.atom1 == idx || b.atom2 == idx)
                            && {
                                let o_idx = if b.atom1 == idx { b.atom2 } else { b.atom1 };
                                mol.atoms[o_idx].atomic_number == 8
                            }
                    })
                    && mol.adjacency[idx].iter().any(|&n| {
                        mol.atoms[n].atomic_number == 8 && mol.atoms[n].charge < -0.5
                    });

                // O bonded to a carboxylate C (for carboxylate O typing)
                let o_bonded_to_carboxylate_c = atom.atomic_number == 8
                    && carbon_neighbors.iter().any(|&c| {
                        mol.atoms[c].atomic_number == 6
                            && count_double_o(c) >= 1
                            && mol.bonds.iter().any(|b| {
                                b.bond_type == BondType::Single
                                    && (b.atom1 == c || b.atom2 == c)
                                    && {
                                        let o_idx = if b.atom1 == c { b.atom2 } else { b.atom1 };
                                        mol.atoms[o_idx].atomic_number == 8
                                    }
                            })
                            && mol.adjacency[c].iter().any(|&n| {
                                mol.atoms[n].atomic_number == 8
                                    && mol.atoms[n].charge < -0.5
                            })
                    });

                // O bonded to an aromatic N (for N-oxide O typing)
                let bonded_to_aromatic_n = atom.atomic_number == 8
                    && mol.adjacency[idx].iter().any(|&n| {
                        mol.atoms[n].atomic_number == 7 && aromatic_atoms.contains(&n)
                    });

                // O⁻ bonded to S/Cl with multiple O neighbors (sulfite/perchlorate) → O_CO2
                let o_bonded_to_oxidized_center = atom.atomic_number == 8
                    && atom.charge < -0.5
                    && mol.adjacency[idx].iter().any(|&n| {
                        let z = mol.atoms[n].atomic_number;
                        if z == 7 { return true; } // N-oxide O⁻
                        if z != 16 && z != 17 { return false; }
                        // The center atom has multiple O neighbors
                        mol.adjacency[n].iter()
                            .filter(|&&m| mol.atoms[m].atomic_number == 8).count() >= 2
                    });

                // Aromatic N bonded to O (for N-oxide N typing)
                let is_noxide_n = atom.atomic_number == 7
                    && aromatic
                    && !oxygen_neighbors.is_empty();

                match (atom.atomic_number, hybrid, aromatic, num_bonds) {
                    // Hydrogen - context-dependent subtyping
                    (1, _, _, _) => {
                        let neighbor_idx = mol.adjacency[idx]
                            .first()
                            .expect("H must have exactly one neighbor");
                        let neighbor = &mol.atoms[*neighbor_idx];
                        Self::determine_h_subtype(idx, mol, *neighbor_idx, &aromatic_atoms)
                    }

                    // Isonitrile C- (sp, triple bond, negative charge) → CID (type 60)
                    (6, Hybridization::Sp1, _, _) if charge < -0.5 => MMFFAtomType::CID,
                    // Guanidinium/amidinium C (C=N+ with >=2 N neighbors) → CGD+ (type 57)
                    (6, _, _, _) if double_bond_partners.contains(&7)
                        && mol.adjacency[idx].iter()
                            .filter(|&&n| mol.atoms[n].atomic_number == 7).count() >= 2
                        && mol.bonds.iter().any(|b| {
                            b.bond_type == BondType::Double
                                && (b.atom1 == idx || b.atom2 == idx)
                                && {
                                    let other = if b.atom1 == idx { b.atom2 } else { b.atom1 };
                                    mol.atoms[other].charge > 0.5
                                }
                        }) => MMFFAtomType::C_AN,
                    // Carbon with formal charge
                    (6, _, _, _) if charge.abs() > 0.5 => {
                        if charge > 0.0 {
                            MMFFAtomType::C_CAT
                        } else {
                            MMFFAtomType::C_AN
                        }
                    }

                    // Carbon types — aromatic takes priority over hybridization
                    (6, _, true, _) => MMFFAtomType::C_AR,
                    // Carbons in small rings: RDKit uses CR3R (sp3, 22) for 3-rings,
                    // CR4R (sp3, 20) / CE4R (sp2, 30) for 4-rings.
                    (6, Hybridization::Sp3, false, _) if in_3ring.contains(&idx) => {
                        MMFFAtomType::CR3R
                    }
                    (6, Hybridization::Sp3, false, _) if in_4ring.contains(&idx) => {
                        MMFFAtomType::CR4R
                    }
                    (6, Hybridization::Sp2, false, _) if in_4ring.contains(&idx) => {
                        // C with double bond to heteroatom in 4-ring → C_2, not CE4R
                        let has_double_to_het = mol.bonds.iter().any(|b| {
                            b.bond_type == BondType::Double
                                && (b.atom1 == idx || b.atom2 == idx)
                                && {
                                    let other = if b.atom1 == idx { b.atom2 } else { b.atom1 };
                                    let z = mol.atoms[other].atomic_number;
                                    z != 6 && z != 1
                                }
                        });
                        if has_double_to_het { MMFFAtomType::C_2 } else { MMFFAtomType::CE4R }
                    }
                    (6, Hybridization::Sp3, false, 1..=4) => MMFFAtomType::C_3,
                    // sp2 C: carbonyl/imine/thiocarbonyl C (double bond to N/O/S) is MMFF 3;
                    // generic alkene C (double bond only to C) is MMFF 2 (RDKit-verified)
                    // Carboxylate C (C=O + C-[O-]) is MMFF 41 (CO2M) — must check first
                    (6, Hybridization::Sp2, false, 2..) if is_carboxylate_c => {
                        MMFFAtomType::C_CO2
                    }
                    (6, Hybridization::Sp2, false, 2..) => {
                        if has_double_to_nos {
                            MMFFAtomType::C_2
                        } else {
                            MMFFAtomType::C_VIN
                        }
                    }
                    (6, Hybridization::Sp1, _, 1..=2) => MMFFAtomType::C_1,

                    // Nitro N: exactly 2 O neighbors, >=1 N=O double bond, no H (MMFF 45)
                    (7, _, false, _) if is_nitro_n => MMFFAtomType::N_NO2,

                    // Sulfonamide N: bonded to SO2 sulfur (MMFF 43)
                    (7, _, false, _) if bonded_to_so2_s => MMFFAtomType::N_SO2,

                    // Trimethylamine N-oxide N+ (quaternary N with O neighbor) → type 68
                    (7, Hybridization::Sp3, _, 4) if charge > 0.5 && !oxygen_neighbors.is_empty() => MMFFAtomType::N_T3,
                    // Nitrogen with formal charge +1 (quaternary ammonium). Exclude
                    // acyl/amidine/enamine Ns (e.g. charged guanidinium N), which
                    // RDKit types N_PL3/N_AM, not N_4.
                    (7, Hybridization::Sp3, _, _)
                        if charge > 0.5
                            && !has_c_n_neighbor
                            && !has_c_o_neighbor
                            && !has_c_c_neighbor =>
                    {
                        MMFFAtomType::N_4
                    }

                    // Nitrogen: amide N (aromatic ring + bonded to C=O)
                    (7, _, true, _) if is_noxide_n => MMFFAtomType::N_POX,
                    // Pyridine N-oxide N (aromatic N with O- neighbor, charge +1) → N_POX2 (69)
                    (7, _, true, _) if charge > 0.5 && oxygen_neighbors.iter().any(|&o| mol.atoms[o].charge < -0.5) => MMFFAtomType::N_POX2,
                    // N-methylpyridinium N (aromatic N+, degree 3) → type 58
                    (7, _, true, 3) if charge > 0.5 => MMFFAtomType::N_PYR,
                    // Guanidinium/amidinium N+ (bonded to C that has C=N+ and >=2 N) → NCN+ (55)
                    (7, _, false, _) if mol.adjacency[idx].iter().any(|&c| {
                        mol.atoms[c].atomic_number == 6
                            && mol.adjacency[c].iter()
                                .filter(|&&n| mol.atoms[n].atomic_number == 7).count() >= 2
                            && mol.bonds.iter().any(|b| {
                                b.bond_type == BondType::Double
                                    && (b.atom1 == c || b.atom2 == c)
                                    && {
                                        let other = if b.atom1 == c { b.atom2 } else { b.atom1 };
                                        mol.atoms[other].atomic_number == 7
                                            && mol.atoms[other].charge > 0.5
                                    }
                            })
                    }) => MMFFAtomType::NCN_PLUS,
                    (7, _, true, _) if has_c_o_neighbor && !n_owns_cn => MMFFAtomType::N_AM,

                    // Amide / carbamate / urea N: non-aromatic N directly bonded to a
                    // carbonyl C is MMFF 10 (N_AM). Excludes N=C imines (n_owns_cn)
                    // like methyl isocyanate CH3-N=C=O where N has its own double bond.
                    (7, _, false, _) if has_c_o_neighbor && !n_owns_cn => MMFFAtomType::N_AM,

                    // Enamine / vinylamine N: non-aromatic N bonded to a C=C carbon
                    // is planar → MMFF 40 (N_PL3) (analogous to aniline).
                    (7, _, false, _) if has_c_c_neighbor => MMFFAtomType::N_PL3,

                    // Amidine / guanidinium N (bonded to C=N carbon) → MMFF 40 (N_PL3).
                    (7, _, false, _) if has_c_n_neighbor => MMFFAtomType::N_PL3,

                    // N_PL3: sp3 N bonded to aromatic carbon (aniline-like)
                    (7, Hybridization::Sp3, false, _)
                        if carbon_neighbors
                            .iter()
                            .any(|&c| aromatic_atoms.contains(&c)) =>
                    {
                        MMFFAtomType::N_PL3
                    }

                    // Nitrogen types — aromatic takes priority
                    (7, _, true, 2..=3) => MMFFAtomType::N_AR,
                    (7, _, true, _) => MMFFAtomType::N_AM,
                    (7, Hybridization::Sp3, false, 1..=3) => MMFFAtomType::N_3,
                    // N_PL3: sp2 N connected to aromatic carbon (aniline-like)
                    (7, Hybridization::Sp2, false, _)
                        if carbon_neighbors
                            .iter()
                            .any(|&c| aromatic_atoms.contains(&c)) =>
                    {
                        MMFFAtomType::N_PL3
                    }
                    // Nitroso N (N=O, no formal charge) → type 46
                    (7, Hybridization::Sp2, false, _)
                        if double_o_count == 1 && charge.abs() < 0.5 =>
                    {
                        MMFFAtomType::N_NITROSO
                    }
                    // Cumulated =N= (central N with 2+ double bonds) → type 53 (N_2Z)
                    (7, _, false, _) if double_bond_partners.len() >= 2 => {
                        MMFFAtomType::N_2Z
                    }
                    // Terminal anionic N (azide/diazo N-) → type 47 (N_1M)
                    (7, _, false, _) if charge < -0.5 => {
                        MMFFAtomType::N_1M
                    }
                    // Nitrene/radical N (type 62): sp2 N with degree 2 and no double bonds
                    (7, Hybridization::Sp2, false, 2) if double_bond_partners.is_empty() => MMFFAtomType::N_RAD,
                    (7, Hybridization::Sp2, false, 2..) => MMFFAtomType::N_2,
                    // Isonitrile/nitrile-oxide N+ (sp, triple bond, positive charge) → NID (type 61)
                    (7, Hybridization::Sp1, false, _) if charge > 0.5 => MMFFAtomType::NID,
                    (7, Hybridization::Sp1, false, 1..=2) => MMFFAtomType::N_1,

                    // Nitrile oxide O- (bonded to N+ with exactly 2 neighbors) → type 35
                    (8, _, _, _) if charge < -0.5
                        && mol.adjacency[idx].iter().any(|&n| {
                            mol.atoms[n].atomic_number == 7
                                && mol.atoms[n].charge > 0.5
                                && mol.adjacency[n].len() == 2
                        }) => MMFFAtomType::OXIDE,

                    // Terminal O on SO2 (sulfone/sulfonamide/sulfonate) or nitro groups
                    // is MMFF 32 (O2CM); bridging ester O's have carbon neighbors and
                    // are excluded (RDKit-verified)
                    (8, _, _, _)
                        if carbon_neighbors.is_empty()
                            && (so2_s_bond_order == Some(2)
                                || (so2_s_bond_order == Some(1) && charge < -0.5)
                                || bonded_to_nitro_n
                                || bonded_to_aromatic_n) =>
                    {
                        MMFFAtomType::O_CO2
                    }

                    // Carboxylate [O-]: single-bonded O on a carboxylate C → MMFF 32 (O2CM)
                    (8, _, _, _) if o_bonded_to_carboxylate_c => MMFFAtomType::O_CO2,
                    // O⁻ on sulfite/perchlorate center → O_CO2
                    (8, _, _, _) if o_bonded_to_oxidized_center => MMFFAtomType::O_CO2,

                    // O double-bonded to P (P=O), Cl (perchlorate), or S (sulfone/sulfite/sulfene) → MMFF 32
                    // Sulfoxide exception: S with degree 3, 2 C neighbors, 1 double O → O stays O_2
                    (8, _, _, _)
                        if mol.bonds.iter().any(|b| {
                            if b.bond_type != BondType::Double
                                || (b.atom1 != idx && b.atom2 != idx)
                            {
                                return false;
                            }
                            let other = if b.atom1 == idx { b.atom2 } else { b.atom1 };
                            match mol.atoms[other].atomic_number {
                                15 | 17 => true, // P=O, Cl=O
                                16 => { // S=O: type 32 unless sulfoxide/sulfinamide
                                    // Sulfoxide exception: S degree 3, 1 double O,
                                    // no other O neighbors (all non-=O neighbors are C/N)
                                    let s_deg = mol.adjacency[other].len();
                                    let s_has_o_single = mol.adjacency[other].iter()
                                        .any(|&n| n != idx && mol.atoms[n].atomic_number == 8);
                                    let s_double_o_total = mol.adjacency[other].iter()
                                        .filter(|&&n| mol.bonds.iter().any(|b2| {
                                            b2.bond_type == BondType::Double
                                                && (b2.atom1 == other && b2.atom2 == n
                                                    || b2.atom1 == n && b2.atom2 == other)
                                                && mol.atoms[n].atomic_number == 8
                                        })).count();
                                    // Sulfone (2+ double O) or sulfite (O single neighbor) → 32
                                    // Sulfoxide (degree 3, 1 double O, no O single) → 7
                                    s_double_o_total >= 2 || s_has_o_single
                                }
                                _ => false,
                            }
                        }) =>
                    {
                        MMFFAtomType::O_CO2
                    },

                    // Aromatic oxonium O+ (type 51): O+ with degree 2 (furanium)
                    (8, _, _, 2) if charge > 0.5 => MMFFAtomType::O_2P,
                    // O_CO2: O double-bonded to C that is also double-bonded to another O (CO2),
                    // OR C=O on a carboxylate C (where the other O has formal charge -1)
                    (8, _, _, _) if has_double_bond_to_c => {
                        let is_co2 = mol.bonds.iter().any(|b| {
                            if b.bond_type != BondType::Double {
                                return false;
                            }
                            let c_idx = if b.atom1 == idx {
                                b.atom2
                            } else if b.atom2 == idx {
                                b.atom1
                            } else {
                                return false;
                            };
                            if mol.atoms[c_idx].atomic_number != 6 {
                                return false;
                            }
                            // Check 1: C has another C=O double bond (e.g. CO2, neutral carboxyl)
                            let has_second_double_o = mol.bonds.iter().any(|b2| {
                                b2.bond_type == BondType::Double
                                    && b2.atom1 != idx
                                    && b2.atom2 != idx
                                    && (b2.atom1 == c_idx || b2.atom2 == c_idx)
                                    && (mol.atoms[b2.atom1].atomic_number == 8
                                        || mol.atoms[b2.atom2].atomic_number == 8)
                            });
                            // Check 2: C is a carboxylate C (has single-bonded O with charge -1)
                            let has_carboxylate_o = mol.bonds.iter().any(|b2| {
                                b2.bond_type == BondType::Single
                                    && (b2.atom1 == c_idx || b2.atom2 == c_idx)
                                    && {
                                        let o_idx = if b2.atom1 == c_idx {
                                            b2.atom2
                                        } else {
                                            b2.atom1
                                        };
                                        mol.atoms[o_idx].atomic_number == 8
                                            && mol.atoms[o_idx].charge < -0.5
                                    }
                            });
                            has_second_double_o || has_carboxylate_o
                        });
                        if is_co2 {
                            MMFFAtomType::O_CO2
                        } else {
                            MMFFAtomType::O_2
                        }
                    }

                    // Water O: both neighbors are H (true H2O). Not O-S/P-H hydroxyls.
                    (8, Hybridization::Sp3, _, 2)
                        if mol.adjacency[idx]
                            .iter()
                            .all(|&n| mol.atoms[n].atomic_number == 1) =>
                    {
                        MMFFAtomType::OH2
                    }

                    // O_R: sp3 O bonded to 2 carbons (ether)
                    (8, Hybridization::Sp3, _, 2) if carbon_neighbors.len() == 2 => {
                        MMFFAtomType::O_R
                    }

                    // Alcohol O: sp3 O bonded to 1 C and 1 H
                    (8, Hybridization::Sp3, _, 2) if carbon_neighbors.len() == 1 => {
                        MMFFAtomType::O_R
                    }

                    // Oxygen types
                    // Oxonium O+ (type 49): 3-coordinate sp3 O with positive charge
                    (8, Hybridization::Sp3, _, 3) if charge > 0.5 => MMFFAtomType::O_3P,
                    (8, Hybridization::Sp3, _, 3..) => MMFFAtomType::O_3,
                    (8, Hybridization::Sp2, _, _) => MMFFAtomType::O_2,
                    (8, Hybridization::Sp3, true, _) => MMFFAtomType::O_3_Z,

                    // Sulfur types — aromatic takes priority
                    (16, _, true, 2..=3) => MMFFAtomType::S_AR,
                    // Sulfene S (type 74): S with C=S=O (1 double to C + 1 double to O)
                    (16, _, _, 2) if double_o_count == 1
                        && double_bond_partners.iter().any(|&an| an == 6) =>
                    {
                        MMFFAtomType::S_CSO
                    }
                    // Oxidized sulfur: SO2 -> MMFF 18, S=O -> MMFF 17 (RDKit-verified)
                    (16, _, _, _) if double_o_count >= 2 => MMFFAtomType::S_O2,
                    // Sulfite S (type 73): S with 1 double O + an O single-bond neighbor
                    (16, _, _, 3) if double_o_count == 1
                        && mol.adjacency[idx].iter().any(|&n| {
                            mol.atoms[n].atomic_number == 8
                                && mol.bonds.iter().any(|b| {
                                    (b.atom1 == idx && b.atom2 == n
                                        || b.atom1 == n && b.atom2 == idx)
                                        && b.bond_type == BondType::Single
                                })
                        }) =>
                    {
                        MMFFAtomType::S_O3
                    }
                    (16, _, _, _) if double_o_count == 1 => MMFFAtomType::S_OX,
                    // CS2 / thiocarbonyl: S double-bonded to a C that is itself
                    // double-bonded to another S -> MMFF 72 (S2CM).
                    (16, _, _, _)
                        if mol.bonds.iter().any(|b| {
                            if b.bond_type != BondType::Double
                                || (b.atom1 != idx && b.atom2 != idx)
                            {
                                return false;
                            }
                            let c = if b.atom1 == idx { b.atom2 } else { b.atom1 };
                            mol.atoms[c].atomic_number == 6
                                && mol.bonds.iter().any(|b2| {
                                    b2.bond_type == BondType::Double
                                        && b2.atom1 != idx
                                        && b2.atom2 != idx
                                        && (b2.atom1 == c || b2.atom2 == c)
                                        && mol.atoms[if b2.atom1 == c {
                                            b2.atom2
                                        } else {
                                            b2.atom1
                                        }]
                                        .atomic_number
                                            == 16
                                })
                        }) =>
                    {
                        MMFFAtomType::S2CM
                    }
                    // Sulfonium S+ (3 bonds, positive charge) → type 17
                    (16, _, _, 3) if charge > 0.5 => MMFFAtomType::S_OX,
                    // Hypervalent S (4+ bonds) → type 18 (SX4) — RDKit-verified
                    (16, _, _, 4..) => MMFFAtomType::S_O2,
                    (16, Hybridization::Sp3, _, 2..) => MMFFAtomType::S_3,
                    (16, Hybridization::Sp2, _, _) => MMFFAtomType::S_2,

                    // Phosphorus types
                    // P in small ring with double bond → type 75
                    (15, _, _, 2) if double_bond_partners.len() >= 1 => MMFFAtomType::P_ARM,
                    (15, Hybridization::Sp3, _, 3..=4) => MMFFAtomType::P_3,
                    (15, Hybridization::Sp2, _, _) => MMFFAtomType::P_4,
                    (15, Hybridization::Sp3D, _, 5) => MMFFAtomType::P_3D,

                    // Silicon (MMFF 19) — always type 19 regardless of hybridization
                    (14, _, _, _) => MMFFAtomType::Si,

                    // Metal ions by element + formal charge (RDKit-verified)
                    (3, _, _, _) if charge > 0.5 => MMFFAtomType::Li,
                    (11, _, _, _) if charge > 0.5 => MMFFAtomType::Na,
                    (19, _, _, _) if charge > 0.5 => MMFFAtomType::K,
                    (12, _, _, _) if charge > 1.5 => MMFFAtomType::Mg_P2,
                    (20, _, _, _) if charge > 1.5 => MMFFAtomType::Ca_P2,
                    (30, _, _, _) if charge > 1.5 => MMFFAtomType::Zn_P2,
                    (26, _, _, _) if charge > 2.5 => MMFFAtomType::Fe_P3,
                    (26, _, _, _) if charge > 1.5 => MMFFAtomType::Fe_P2,
                    (29, _, _, _) if charge > 1.5 => MMFFAtomType::Cu_P2,
                    (29, _, _, _) if charge > 0.5 => MMFFAtomType::Cu_P1,

                    // Halide anions by formal charge (RDKit-verified; no I- in MMFF94)
                    (9, _, _, _) if charge < -0.5 => MMFFAtomType::F_M,
                    (17, _, _, _) if charge < -0.5 => MMFFAtomType::CL_M,
                    (35, _, _, _) if charge < -0.5 => MMFFAtomType::BR_M,

                    // Halogens
                    (9, _, _, _) => MMFFAtomType::F,
                    // Perchlorate Cl (type 77): 4+ coordinate Cl
                    (17, _, _, 4..) => MMFFAtomType::CL4,
                    (17, _, _, _) => MMFFAtomType::Cl,
                    (35, _, _, _) => MMFFAtomType::Br,
                    (53, _, _, _) => MMFFAtomType::I,

                    // Default fallback
                    _ => match atom.atomic_number {
                        6 => MMFFAtomType::C_3,
                        7 => MMFFAtomType::N_3,
                        8 => MMFFAtomType::O_3,
                        _ => MMFFAtomType::C_3,
                    },
                }
            })
            .collect()
    }

    fn calculate_charges(mol: &Molecule, atom_types: &[MMFFAtomType]) -> Vec<f64> {
        calculate_bci_charges(mol, atom_types)
    }

    fn determine_h_subtype(
        h_idx: usize,
        mol: &Molecule,
        neighbor_idx: usize,
        aromatic_atoms: &std::collections::HashSet<usize>,
    ) -> MMFFAtomType {
        let neighbor = &mol.atoms[neighbor_idx];
        match neighbor.atomic_number {
            8 => {
                // H on positively-charged oxonium O+ → MMFF 50 (H_OXP)
                if neighbor.charge > 0.5 {
                    return MMFFAtomType::H_OXP;
                }

                let o_neighbors: Vec<usize> = mol.adjacency[neighbor_idx]
                    .iter()
                    .filter(|&&n| n != h_idx)
                    .copied()
                    .collect();

                if o_neighbors.is_empty() {
                    return MMFFAtomType::H_OH;
                }

                let all_h = o_neighbors.iter().all(|&n| mol.atoms[n].atomic_number == 1);
                if all_h {
                    return MMFFAtomType::H_OH;
                }

                let non_h_neighbor = o_neighbors
                    .iter()
                    .find(|&&n| mol.atoms[n].atomic_number != 1);

                if let Some(&nn) = non_h_neighbor {
                    let nn_atom = &mol.atoms[nn];
                    match nn_atom.atomic_number {
                        // Sulfonic acid S-OH -> MMFF 33 (HOS)
                        16 => MMFFAtomType::HOS,
                        // Phosphonic / phosphoric acid P-OH -> MMFF 24 (H_COOH);
                        // RDKit reuses the carboxylic-acid H type for P-OH.
                        15 => MMFFAtomType::H_COOH,
                        6 => {
                            if aromatic_atoms.contains(&nn) {
                                MMFFAtomType::H_OAR
                            } else {
                                let has_double_bond = mol.bonds.iter().any(|b| {
                                    (b.atom1 == nn || b.atom2 == nn)
                                        && b.bond_type == BondType::Double
                                });
                                if has_double_bond {
                                    MMFFAtomType::H_COOH
                                } else {
                                    MMFFAtomType::H_ONC
                                }
                            }
                        }
                        _ => MMFFAtomType::H_ONC,
                    }
                } else {
                    MMFFAtomType::H_OH
                }
            }
            7 => {
                // Guanidinium/amidinium N-H → HNRP (36): N bonded to C with C=N+ and >=2 N
                let n_is_guanidinium = mol.adjacency[neighbor_idx].iter().any(|&c| {
                    mol.atoms[c].atomic_number == 6
                        && mol.adjacency[c].iter()
                            .filter(|&&n| mol.atoms[n].atomic_number == 7).count() >= 2
                        && mol.bonds.iter().any(|b| {
                            b.bond_type == BondType::Double
                                && (b.atom1 == c || b.atom2 == c)
                                && {
                                    let o = if b.atom1 == c { b.atom2 } else { b.atom1 };
                                    mol.atoms[o].atomic_number == 7 && mol.atoms[o].charge > 0.5
                                }
                        })
                });
                if n_is_guanidinium {
                    return MMFFAtomType::HNRP;
                }
                // Imine H: H on sp2 N with C=N double bond, not aromatic → MMFF 27 (H_NIM)
                let n_has_double_to_c = mol.bonds.iter().any(|b| {
                    b.bond_type == BondType::Double
                        && (b.atom1 == neighbor_idx || b.atom2 == neighbor_idx)
                        && {
                            let other = if b.atom1 == neighbor_idx {
                                b.atom2
                            } else {
                                b.atom1
                            };
                            mol.atoms[other].atomic_number == 6
                        }
                });
                let n_is_aromatic = aromatic_atoms.contains(&neighbor_idx);
                // Acyl/amidine N: N has a carbon neighbor that is double-bonded to O
                // (amide) or to another N (amidine/guanidinium).
                let n_is_acyl_or_amidine = mol.adjacency[neighbor_idx].iter().any(|&nbr| {
                    mol.atoms[nbr].atomic_number == 6
                        && mol.bonds.iter().any(|b2| {
                            b2.bond_type == BondType::Double
                                && (b2.atom1 == nbr || b2.atom2 == nbr)
                                && {
                                    let o = if b2.atom1 == nbr { b2.atom2 } else { b2.atom1 };
                                    let an = mol.atoms[o].atomic_number;
                                    an == 8 || an == 7 || an == 6
                                }
                        })
                });
                // Sulfonamide N: N bonded to an SO2 sulfur.
                let n_is_sulfonamide = mol.adjacency[neighbor_idx].iter().any(|&nbr| {
                    mol.atoms[nbr].atomic_number == 16
                        && mol.bonds.iter().filter(|b| {
                            b.bond_type == BondType::Double
                                && (b.atom1 == nbr || b.atom2 == nbr)
                        }).count() >= 2
                });
                // Aniline N: non-aromatic N bonded to an aromatic C.
                let n_is_aniline = !n_is_aromatic
                    && mol.adjacency[neighbor_idx]
                        .iter()
                        .any(|&nbr| aromatic_atoms.contains(&nbr) && mol.atoms[nbr].atomic_number == 6);

                if n_has_double_to_c && !n_is_aromatic {
                    // Imine N-H (C=N-H) → MMFF 27
                    MMFFAtomType::H_NIM
                } else if !n_is_aromatic
                    && (n_is_acyl_or_amidine || n_is_sulfonamide || n_is_aniline)
                {
                    // Amide / amidine / enamine / sulfonamide / aniline N-H → MMFF 28
                    // (covers charged guanidinium N-H, which is amidinium)
                    MMFFAtomType::H_NAM
                } else if mol.atoms[neighbor_idx].charge > 0.5 {
                    // Ammonium (sp3 N+) N-H → MMFF 36 (HNR+)
                    MMFFAtomType::HNRP
                } else {
                    // Generic amine or aromatic-ring N-H (pyrrole/indole) → MMFF 23
                    MMFFAtomType::H_N3
                }
            }
            // H bonded to S (thiol) or P → MMFF 71
            // (RDKit uses type 71 for both H-S and H-P)
            16 | 15 => MMFFAtomType::HS,
            _ => MMFFAtomType::H,
        }
    }

    pub fn calculate_energy_and_gradient(&self, coords: &[[f64; 3]]) -> (f64, Vec<[f64; 3]>) {
        let mut energy = 0.0;
        let mut gradient = vec![[0.0; 3]; self.mol.atoms.len()];

        // Bond stretching
        for bond in &self.mol.bonds {
            if let Some(params) = get_bond_params(
                self.atom_types[bond.atom1],
                self.atom_types[bond.atom2],
                bond.bond_type,
            ) {
                energy += bond_energy(coords, bond.atom1, bond.atom2, &params);
                let (gi, gj) = bond_gradient(coords, bond.atom1, bond.atom2, &params);
                gradient[bond.atom1][0] += gi[0];
                gradient[bond.atom1][1] += gi[1];
                gradient[bond.atom1][2] += gi[2];
                gradient[bond.atom2][0] += gj[0];
                gradient[bond.atom2][1] += gj[1];
                gradient[bond.atom2][2] += gj[2];
            }
        }

        // Angle bending
        for angle in &self.angles {
            let (i, j, k) = (angle.atom1, angle.atom2, angle.atom3);
            let bij_key = (i.min(j), i.max(j));
            let bkj_key = (k.min(j), k.max(j));
            let bt_ij = self
                .bond_map
                .get(&bij_key)
                .map(|b| get_mmff_bond_type(b.bond_type, self.type_ids[i], self.type_ids[j]))
                .unwrap_or(0);
            let bt_jk = self
                .bond_map
                .get(&bkj_key)
                .map(|b| get_mmff_bond_type(b.bond_type, self.type_ids[j], self.type_ids[k]))
                .unwrap_or(0);
            let ring_size = self.angle_ring_size(i, j, k);
            let r0_ij = self
                .bond_map
                .get(&bij_key)
                .and_then(|b| get_bond_params(self.atom_types[i], self.atom_types[j], b.bond_type))
                .map(|p| p.r0)
                .unwrap_or(1.5);
            let r0_jk = self
                .bond_map
                .get(&bkj_key)
                .and_then(|b| get_bond_params(self.atom_types[k], self.atom_types[j], b.bond_type))
                .map(|p| p.r0)
                .unwrap_or(1.5);
            if let Some(params) = get_angle_params_with_bond_info(
                self.atom_types[i],
                self.atom_types[j],
                self.atom_types[k],
                bt_ij,
                bt_jk,
                ring_size,
                r0_ij,
                r0_jk,
            ) {
                energy += angle_energy(coords, i, j, k, &params);
                let (g1, g2, g3) = angle_gradient(coords, i, j, k, &params);
                gradient[i][0] += g1[0];
                gradient[i][1] += g1[1];
                gradient[i][2] += g1[2];
                gradient[j][0] += g2[0];
                gradient[j][1] += g2[1];
                gradient[j][2] += g2[2];
                gradient[k][0] += g3[0];
                gradient[k][1] += g3[1];
                gradient[k][2] += g3[2];
            }
        }

        // Stretch-bend coupling
        for angle in &self.angles {
            let (i, j, k) = (angle.atom1, angle.atom2, angle.atom3);
            let bij_key = (i.min(j), i.max(j));
            let bkj_key = (k.min(j), k.max(j));
            let bond_ij = self.bond_map.get(&bij_key);
            let bond_kj = self.bond_map.get(&bkj_key);
            if let (Some(bij), Some(bkj)) = (bond_ij, bond_kj) {
                let bt_ij = get_mmff_bond_type(bij.bond_type, self.type_ids[i], self.type_ids[j]);
                let bt_jk = get_mmff_bond_type(bkj.bond_type, self.type_ids[j], self.type_ids[k]);
                let ring_size = self.angle_ring_size(i, j, k);
                let angle_type_val = mmff_tables::compute_angle_type(bt_ij, bt_jk, ring_size);
                if mmff_tables::is_linear_center(self.type_ids[j]) {
                    continue;
                }
                if let (
                    Some(sb_params),
                    Some(bond_params_ij),
                    Some(bond_params_kj),
                    Some(angle_params),
                ) = (
                    get_stretch_bend_params(
                        self.atom_types[i],
                        self.atom_types[j],
                        self.atom_types[k],
                        bt_ij,
                        bt_jk,
                        self.mol.atoms[i].atomic_number,
                        self.mol.atoms[j].atomic_number,
                        self.mol.atoms[k].atomic_number,
                        angle_type_val,
                    ),
                    get_bond_params(self.atom_types[i], self.atom_types[j], bij.bond_type),
                    get_bond_params(self.atom_types[k], self.atom_types[j], bkj.bond_type),
                    get_angle_params_with_bond_info(
                        self.atom_types[i],
                        self.atom_types[j],
                        self.atom_types[k],
                        bt_ij,
                        bt_jk,
                        ring_size,
                        get_bond_params(self.atom_types[i], self.atom_types[j], bij.bond_type)
                            .map(|p| p.r0)
                            .unwrap_or(1.5),
                        get_bond_params(self.atom_types[k], self.atom_types[j], bkj.bond_type)
                            .map(|p| p.r0)
                            .unwrap_or(1.5),
                    ),
                ) {
                    energy += stretch_bend_energy(
                        coords,
                        i,
                        j,
                        k,
                        bond_params_ij.r0,
                        bond_params_kj.r0,
                        angle_params.theta0.to_radians(),
                        &sb_params,
                    );
                    let (g1, g2, g3) = stretch_bend_gradient(
                        coords,
                        i,
                        j,
                        k,
                        bond_params_ij.r0,
                        bond_params_kj.r0,
                        angle_params.theta0.to_radians(),
                        &sb_params,
                    );
                    gradient[i][0] += g1[0];
                    gradient[i][1] += g1[1];
                    gradient[i][2] += g1[2];
                    gradient[j][0] += g2[0];
                    gradient[j][1] += g2[1];
                    gradient[j][2] += g2[2];
                    gradient[k][0] += g3[0];
                    gradient[k][1] += g3[1];
                    gradient[k][2] += g3[2];
                }
            }
        }

        // Torsion
        for (ti, torsion) in self.torsions.iter().enumerate() {
            let tor_type = self.torsion_types[ti];
            if let Some(params) = get_torsion_params(
                self.atom_types[torsion.atom1],
                self.atom_types[torsion.atom2],
                self.atom_types[torsion.atom3],
                self.atom_types[torsion.atom4],
                tor_type,
                self.variant,
            ) {
                energy += torsion_energy(
                    coords,
                    torsion.atom1,
                    torsion.atom2,
                    torsion.atom3,
                    torsion.atom4,
                    &params,
                );
                let (g1, g2, g3, g4) = torsion_gradient(
                    coords,
                    torsion.atom1,
                    torsion.atom2,
                    torsion.atom3,
                    torsion.atom4,
                    &params,
                );
                gradient[torsion.atom1][0] += g1[0];
                gradient[torsion.atom1][1] += g1[1];
                gradient[torsion.atom1][2] += g1[2];
                gradient[torsion.atom2][0] += g2[0];
                gradient[torsion.atom2][1] += g2[1];
                gradient[torsion.atom2][2] += g2[2];
                gradient[torsion.atom3][0] += g3[0];
                gradient[torsion.atom3][1] += g3[1];
                gradient[torsion.atom3][2] += g3[2];
                gradient[torsion.atom4][0] += g4[0];
                gradient[torsion.atom4][1] += g4[1];
                gradient[torsion.atom4][2] += g4[2];
            }
        }

        // Out-of-plane
        for oop in &self.oops {
            let params = get_oop_params(
                self.atom_types[oop.atom1],
                self.atom_types[oop.central],
                self.atom_types[oop.atom2],
                self.atom_types[oop.atom3],
                self.variant,
            );
            energy += oop_energy(
                coords,
                oop.central,
                oop.atom1,
                oop.atom2,
                oop.atom3,
                &params,
            );
            let (g_central, g1, g2, g3) = oop_gradient(
                coords,
                oop.central,
                oop.atom1,
                oop.atom2,
                oop.atom3,
                &params,
            );
            gradient[oop.central][0] += g_central[0];
            gradient[oop.central][1] += g_central[1];
            gradient[oop.central][2] += g_central[2];
            gradient[oop.atom1][0] += g1[0];
            gradient[oop.atom1][1] += g1[1];
            gradient[oop.atom1][2] += g1[2];
            gradient[oop.atom2][0] += g2[0];
            gradient[oop.atom2][1] += g2[1];
            gradient[oop.atom2][2] += g2[2];
            gradient[oop.atom3][0] += g3[0];
            gradient[oop.atom3][1] += g3[1];
            gradient[oop.atom3][2] += g3[2];
        }

        // Van der Waals (all nonbonded pairs, excluding 1-2 and 1-3)
        let n = self.mol.atoms.len();
        for i in 0..n {
            for j in (i + 1)..n {
                if !self.excluded_pairs.contains(&(i, j)) {
                    let params_i = get_vdw_params(self.atom_types[i]);
                    let params_j = get_vdw_params(self.atom_types[j]);
                    let is_14 = self.one_four_pairs.contains(&(i, j));
                    let (e, grad_i, grad_j) =
                        vdw_energy_and_gradient(coords, i, j, &params_i, &params_j, is_14);
                    energy += e;
                    gradient[i][0] += grad_i[0];
                    gradient[i][1] += grad_i[1];
                    gradient[i][2] += grad_i[2];
                    gradient[j][0] += grad_j[0];
                    gradient[j][1] += grad_j[1];
                    gradient[j][2] += grad_j[2];
                }
            }
        }

        // Electrostatics (all charged pairs, excluding 1-2 and 1-3)
        // 1-4 pairs are scaled by 0.75
        for i in 0..n {
            for j in (i + 1)..n {
                if !self.excluded_pairs.contains(&(i, j))
                    && (self.charges[i].abs() > 1e-6 || self.charges[j].abs() > 1e-6)
                {
                    let is_14 = self.one_four_pairs.contains(&(i, j));
                    let (e, grad_i, grad_j) =
                        electrostatic_energy_and_gradient(coords, &self.charges, i, j, 1.0, is_14);
                    energy += e;
                    gradient[i][0] += grad_i[0];
                    gradient[i][1] += grad_i[1];
                    gradient[i][2] += grad_i[2];
                    gradient[j][0] += grad_j[0];
                    gradient[j][1] += grad_j[1];
                    gradient[j][2] += grad_j[2];
                }
            }
        }

        (energy, gradient)
    }

    pub fn calculate_energy(&self, coords: &[[f64; 3]]) -> f64 {
        self.calculate_energy_and_gradient(coords).0
    }

    pub fn calculate_gradient(&self, coords: &[[f64; 3]]) -> Vec<[f64; 3]> {
        self.calculate_energy_and_gradient(coords).1
    }

    pub fn calculate_energy_breakdown(&self, coords: &[[f64; 3]]) -> EnergyBreakdown {
        let mut bd = EnergyBreakdown::default();
        let n = self.mol.atoms.len();

        for bond in &self.mol.bonds {
            if let Some(params) = get_bond_params(
                self.atom_types[bond.atom1],
                self.atom_types[bond.atom2],
                bond.bond_type,
            ) {
                bd.bond += bond_energy(coords, bond.atom1, bond.atom2, &params);
            }
        }

        for angle in &self.angles {
            let (i, j, k) = (angle.atom1, angle.atom2, angle.atom3);
            let bij_key = (i.min(j), i.max(j));
            let bkj_key = (k.min(j), k.max(j));
            let bt_ij = self
                .bond_map
                .get(&bij_key)
                .map(|b| get_mmff_bond_type(b.bond_type, self.type_ids[i], self.type_ids[j]))
                .unwrap_or(0);
            let bt_jk = self
                .bond_map
                .get(&bkj_key)
                .map(|b| get_mmff_bond_type(b.bond_type, self.type_ids[j], self.type_ids[k]))
                .unwrap_or(0);
            let ring_size = self.angle_ring_size(i, j, k);
            let r0_ij = self
                .bond_map
                .get(&bij_key)
                .and_then(|b| get_bond_params(self.atom_types[i], self.atom_types[j], b.bond_type))
                .map(|p| p.r0)
                .unwrap_or(1.5);
            let r0_jk = self
                .bond_map
                .get(&bkj_key)
                .and_then(|b| get_bond_params(self.atom_types[k], self.atom_types[j], b.bond_type))
                .map(|p| p.r0)
                .unwrap_or(1.5);
            if let Some(params) = get_angle_params_with_bond_info(
                self.atom_types[i],
                self.atom_types[j],
                self.atom_types[k],
                bt_ij,
                bt_jk,
                ring_size,
                r0_ij,
                r0_jk,
            ) {
                bd.angle += angle_energy(coords, i, j, k, &params);
            }
        }

        for angle in &self.angles {
            let (i, j, k) = (angle.atom1, angle.atom2, angle.atom3);
            let bij_key = (i.min(j), i.max(j));
            let bkj_key = (k.min(j), k.max(j));
            let bond_ij = self.bond_map.get(&bij_key);
            let bond_kj = self.bond_map.get(&bkj_key);
            if let (Some(bij), Some(bkj)) = (bond_ij, bond_kj) {
                let bt_ij = get_mmff_bond_type(bij.bond_type, self.type_ids[i], self.type_ids[j]);
                let bt_jk = get_mmff_bond_type(bkj.bond_type, self.type_ids[j], self.type_ids[k]);
                let ring_size = self.angle_ring_size(i, j, k);
                let angle_type_val = mmff_tables::compute_angle_type(bt_ij, bt_jk, ring_size);
                if mmff_tables::is_linear_center(self.type_ids[j]) {
                    continue;
                }
                if let (
                    Some(sb_params),
                    Some(bond_params_ij),
                    Some(bond_params_kj),
                    Some(angle_params),
                ) = (
                    get_stretch_bend_params(
                        self.atom_types[i],
                        self.atom_types[j],
                        self.atom_types[k],
                        bt_ij,
                        bt_jk,
                        self.mol.atoms[i].atomic_number,
                        self.mol.atoms[j].atomic_number,
                        self.mol.atoms[k].atomic_number,
                        angle_type_val,
                    ),
                    get_bond_params(self.atom_types[i], self.atom_types[j], bij.bond_type),
                    get_bond_params(self.atom_types[k], self.atom_types[j], bkj.bond_type),
                    get_angle_params_with_bond_info(
                        self.atom_types[i],
                        self.atom_types[j],
                        self.atom_types[k],
                        bt_ij,
                        bt_jk,
                        ring_size,
                        get_bond_params(self.atom_types[i], self.atom_types[j], bij.bond_type)
                            .map(|p| p.r0)
                            .unwrap_or(1.5),
                        get_bond_params(self.atom_types[k], self.atom_types[j], bkj.bond_type)
                            .map(|p| p.r0)
                            .unwrap_or(1.5),
                    ),
                ) {
                    bd.stretch_bend += stretch_bend_energy(
                        coords,
                        i,
                        j,
                        k,
                        bond_params_ij.r0,
                        bond_params_kj.r0,
                        angle_params.theta0.to_radians(),
                        &sb_params,
                    );
                }
            }
        }

        for (ti, torsion) in self.torsions.iter().enumerate() {
            let tor_type = self.torsion_types[ti];
            if let Some(params) = get_torsion_params(
                self.atom_types[torsion.atom1],
                self.atom_types[torsion.atom2],
                self.atom_types[torsion.atom3],
                self.atom_types[torsion.atom4],
                tor_type,
                self.variant,
            ) {
                bd.torsion += torsion_energy(
                    coords,
                    torsion.atom1,
                    torsion.atom2,
                    torsion.atom3,
                    torsion.atom4,
                    &params,
                );
            }
        }

        for oop in &self.oops {
            let params = get_oop_params(
                self.atom_types[oop.atom1],
                self.atom_types[oop.central],
                self.atom_types[oop.atom2],
                self.atom_types[oop.atom3],
                self.variant,
            );
            bd.oop += oop_energy(
                coords,
                oop.central,
                oop.atom1,
                oop.atom2,
                oop.atom3,
                &params,
            );
        }

        for i in 0..n {
            for j in (i + 1)..n {
                if !self.excluded_pairs.contains(&(i, j)) {
                    let params_i = get_vdw_params(self.atom_types[i]);
                    let params_j = get_vdw_params(self.atom_types[j]);
                    let is_14 = self.one_four_pairs.contains(&(i, j));
                    let (e, _, _) =
                        vdw_energy_and_gradient(coords, i, j, &params_i, &params_j, is_14);
                    bd.vdw += e;
                }
            }
        }

        for i in 0..n {
            for j in (i + 1)..n {
                if !self.excluded_pairs.contains(&(i, j))
                    && (self.charges[i].abs() > 1e-6 || self.charges[j].abs() > 1e-6)
                {
                    let is_14 = self.one_four_pairs.contains(&(i, j));
                    let (e, _, _) =
                        electrostatic_energy_and_gradient(coords, &self.charges, i, j, 1.0, is_14);
                    bd.electrostatic += e;
                }
            }
        }

        bd
    }
}

#[derive(Debug, Clone, Default)]
pub struct EnergyBreakdown {
    pub bond: f64,
    pub angle: f64,
    pub stretch_bend: f64,
    pub torsion: f64,
    pub oop: f64,
    pub vdw: f64,
    pub electrostatic: f64,
}

impl EnergyBreakdown {
    pub fn total(&self) -> f64 {
        self.bond
            + self.angle
            + self.stretch_bend
            + self.torsion
            + self.oop
            + self.vdw
            + self.electrostatic
    }
}
