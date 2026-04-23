pub mod graph;
pub mod parser;

pub use graph::{Angle, Hybridization, OutOfPlane, Torsion};
pub use parser::*;

/// Atom in a molecule
#[derive(Debug, Clone)]
pub struct Atom {
    pub symbol: String,
    pub atomic_number: u8,
    pub mass: f64,
    pub charge: f64,
    pub position: [f64; 3],
    pub index: usize,
}

/// Bond type enumeration
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BondType {
    Single,
    Double,
    Triple,
    Aromatic,
}

/// Bond stereochemistry (from SDF bond block)
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum BondStereo {
    #[default]
    None,
    Cis,
    Trans,
}

/// Bond between two atoms
#[derive(Debug, Clone)]
pub struct Bond {
    pub atom1: usize,
    pub atom2: usize,
    pub bond_type: BondType,
    pub stereo: BondStereo,
}

impl Default for Bond {
    fn default() -> Self {
        Bond {
            atom1: 0,
            atom2: 0,
            bond_type: BondType::Single,
            stereo: BondStereo::None,
        }
    }
}

/// Molecule structure
#[derive(Debug, Clone)]
pub struct Molecule {
    pub atoms: Vec<Atom>,
    pub bonds: Vec<Bond>,
    pub name: String,
    pub adjacency: Vec<Vec<usize>>,
}
