pub mod graph;
pub mod parser;

pub use graph::*;
pub use parser::*;

pub use graph::Angle;
pub use graph::Hybridization;
pub use graph::OutOfPlane;
pub use graph::Torsion;

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

/// Bond between two atoms
#[derive(Debug, Clone)]
pub struct Bond {
    pub atom1: usize,
    pub atom2: usize,
    pub bond_type: BondType,
}

/// Molecule structure
#[derive(Debug, Clone)]
pub struct Molecule {
    pub atoms: Vec<Atom>,
    pub bonds: Vec<Bond>,
    pub name: String,
}
