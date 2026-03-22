//! MMFF94/MMFF94s force field implementation

use std::collections::HashSet;

use crate::molecule::BondType;
use crate::molecule::Hybridization;
use crate::molecule::Molecule;

pub mod angle;
pub mod atom_types;
pub mod bond;
pub mod electrostatics;
pub mod oop;
pub mod torsion;
pub mod vdw;

pub use angle::*;
pub use atom_types::*;
pub use bond::*;
pub use electrostatics::*;
pub use oop::*;
pub use torsion::*;
pub use vdw::*;

pub use crate::MMFFVariant;

/// MMFF atom type
#[allow(non_camel_case_types)]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum MMFFAtomType {
    // Hydrogen
    H,

    // Carbons
    C_3,
    C_2,
    C_1,
    C_AR,
    C_CAT,
    C_AN,

    // Nitrogens
    N_3,
    N_2,
    N_1,
    N_AR,
    N_PL3,
    N_AM,
    N_4,
    N_2Z,
    N_SOM,

    // Oxygens
    O_3,
    O_2,
    O_R,
    O_CO2,
    O_3_Z,

    // Halogens
    F,
    Cl,
    Br,
    I,

    // Sulfur, Phosphorus
    S_3,
    S_2,
    S_AR,
    P_3,
    P_4,

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
    pub charges: Vec<f64>,
    pub variant: MMFFVariant,
    angles: Vec<crate::molecule::Angle>,
    torsions: Vec<crate::molecule::Torsion>,
    oops: Vec<crate::molecule::OutOfPlane>,
    excluded_pairs: HashSet<(usize, usize)>,
}

impl MMFFForceField {
    pub fn new(mol: &Molecule, variant: MMFFVariant) -> Self {
        let atom_types = Self::assign_atom_types(mol);
        let charges = Self::calculate_charges(mol, &atom_types);

        let angles = crate::molecule::graph::find_angles(mol);
        let torsions = crate::molecule::graph::find_torsions(mol);
        let oops = crate::molecule::graph::find_out_of_planes(mol);

        let mut excluded_pairs = HashSet::new();

        // 1-2 pairs (bonded)
        for bond in &mol.bonds {
            let (a, b) = (bond.atom1.min(bond.atom2), bond.atom1.max(bond.atom2));
            excluded_pairs.insert((a, b));
        }

        // 1-3 pairs (angle endpoints)
        for angle in &angles {
            let (a, c) = (angle.atom1.min(angle.atom3), angle.atom1.max(angle.atom3));
            excluded_pairs.insert((a, c));
        }

        Self {
            mol: mol.clone(),
            atom_types,
            charges,
            variant,
            angles,
            torsions,
            oops,
            excluded_pairs,
        }
    }

    fn assign_atom_types(mol: &Molecule) -> Vec<MMFFAtomType> {
        use crate::molecule::graph::{determine_hybridization, find_rings, get_aromatic_atoms};

        let aromatic_atoms = get_aromatic_atoms(mol);
        let rings = find_rings(mol);

        mol.atoms
            .iter()
            .map(|atom| {
                let idx = atom.index;
                let hybrid = determine_hybridization(idx, mol);
                let aromatic = aromatic_atoms.contains(&idx);
                let num_bonds = crate::molecule::graph::get_neighbors(idx, mol).len();
                let charge = atom.charge;

                let _in_ring = rings.iter().any(|r| r.contains(&idx));

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

                match (atom.atomic_number, hybrid, aromatic, num_bonds) {
                    // Hydrogen
                    (1, _, _, _) => MMFFAtomType::H,

                    // Carbon with formal charge
                    (6, _, _, _) if charge.abs() > 0.5 => {
                        if charge > 0.0 {
                            MMFFAtomType::C_CAT
                        } else {
                            MMFFAtomType::C_AN
                        }
                    }

                    // Carbon types
                    (6, Hybridization::Sp3, false, 1..=4) => MMFFAtomType::C_3,
                    (6, Hybridization::Sp2, _, 2..) => MMFFAtomType::C_2,
                    (6, Hybridization::Sp1, _, 1..=2) => MMFFAtomType::C_1,
                    (6, _, true, _) => MMFFAtomType::C_AR,

                    // Nitrogen with formal charge +1
                    (7, Hybridization::Sp3, _, _) if charge.abs() > 0.5 && charge > 0.0 => {
                        MMFFAtomType::N_4
                    }

                    // Nitrogen: amide N (aromatic ring + bonded to C=O)
                    (7, _, true, _) if has_c_o_neighbor => MMFFAtomType::N_AM,

                    // N_PL3: sp3 N bonded to C=O
                    (7, Hybridization::Sp3, false, _) if has_c_o_neighbor => MMFFAtomType::N_PL3,

                    // Nitrogen types
                    (7, Hybridization::Sp3, false, 1..=3) => MMFFAtomType::N_3,
                    (7, Hybridization::Sp2, false, 2..) => MMFFAtomType::N_2,
                    (7, Hybridization::Sp1, false, 1..=2) => MMFFAtomType::N_1,
                    (7, _, true, 2..=3) => MMFFAtomType::N_AR,
                    (7, _, true, _) => MMFFAtomType::N_AM,

                    // O_CO2: sp2 O double-bonded to carbon
                    (8, _, _, _) if has_double_bond_to_c => MMFFAtomType::O_CO2,

                    // O_R: sp3 O bonded to 2 carbons (ether)
                    (8, Hybridization::Sp3, _, 2) if carbon_neighbors.len() == 2 => {
                        MMFFAtomType::O_R
                    }

                    // Oxygen types
                    (8, Hybridization::Sp3, _, 2..) => MMFFAtomType::O_3,
                    (8, Hybridization::Sp2, _, _) => MMFFAtomType::O_2,
                    (8, _, true, 2) => MMFFAtomType::O_R,
                    (8, Hybridization::Sp3, true, _) => MMFFAtomType::O_3_Z,

                    // Sulfur types
                    (16, Hybridization::Sp3, _, 2..) => MMFFAtomType::S_3,
                    (16, Hybridization::Sp2, _, _) => MMFFAtomType::S_2,
                    (16, _, true, 2..=3) => MMFFAtomType::S_AR,

                    // Phosphorus types
                    (15, Hybridization::Sp3, _, 3..=4) => MMFFAtomType::P_3,
                    (15, Hybridization::Sp2, _, _) => MMFFAtomType::P_4,

                    // Halogens
                    (9, _, _, _) => MMFFAtomType::F,
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

    fn calculate_charges(mol: &Molecule, _atom_types: &[MMFFAtomType]) -> Vec<f64> {
        // TODO: Implement bond charge increment method
        vec![0.0; mol.atoms.len()]
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
            if let Some(params) = get_angle_params(
                self.atom_types[angle.atom1],
                self.atom_types[angle.atom2],
                self.atom_types[angle.atom3],
            ) {
                energy += angle_energy(coords, angle.atom1, angle.atom2, angle.atom3, &params);
                let (g1, g2, g3) =
                    angle_gradient(coords, angle.atom1, angle.atom2, angle.atom3, &params);
                gradient[angle.atom1][0] += g1[0];
                gradient[angle.atom1][1] += g1[1];
                gradient[angle.atom1][2] += g1[2];
                gradient[angle.atom2][0] += g2[0];
                gradient[angle.atom2][1] += g2[1];
                gradient[angle.atom2][2] += g2[2];
                gradient[angle.atom3][0] += g3[0];
                gradient[angle.atom3][1] += g3[1];
                gradient[angle.atom3][2] += g3[2];
            }
        }

        // Torsion
        for torsion in &self.torsions {
            if let Some(params) = get_torsion_params(
                self.atom_types[torsion.atom1],
                self.atom_types[torsion.atom2],
                self.atom_types[torsion.atom3],
                self.atom_types[torsion.atom4],
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
            let params = get_oop_params(self.atom_types[oop.central], self.variant);
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
                    let (e, grad_i, grad_j) =
                        vdw_energy_and_gradient(coords, i, j, &params_i, &params_j);
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

        // Electrostatics (all charged pairs)
        for i in 0..n {
            for j in (i + 1)..n {
                if self.charges[i].abs() > 1e-6 || self.charges[j].abs() > 1e-6 {
                    let (e, grad_i, grad_j) =
                        electrostatic_energy_and_gradient(coords, &self.charges, i, j, 1.0);
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
}
