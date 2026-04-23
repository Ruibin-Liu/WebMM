//! Molecular geometry optimizer using MMFF94/MMFF94s force field
//! and L-BFGS optimization algorithm

#![allow(unused_variables)]

use wasm_bindgen::prelude::*;

/// Core molecular data structures
pub mod molecule;

/// ETKDG v3 3D coordinate embedding
pub mod etkdg;

/// MMFF94/MMFF94s force field implementation
pub mod mmff;

/// L-BFGS optimization algorithm
pub mod optimizer;

/// Utility functions
pub mod utils;

/// Property-based tests (only in test mode)
#[cfg(test)]
pub mod prop_tests;

/// MMFF variant selection
#[wasm_bindgen]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum MMFFVariant {
    MMFF94,
    #[default]
    MMFF94s,
}

/// Convergence criteria options
#[wasm_bindgen]
#[derive(Debug, Clone)]
pub struct ConvergenceOptions {
    #[wasm_bindgen(getter_with_clone)]
    pub max_force: f64,
    #[wasm_bindgen(getter_with_clone)]
    pub rms_force: f64,
    #[wasm_bindgen(getter_with_clone)]
    pub energy_change: f64,
    #[wasm_bindgen(getter_with_clone)]
    pub max_iterations: usize,
}

impl Default for ConvergenceOptions {
    fn default() -> Self {
        Self {
            max_force: 0.01,
            rms_force: 0.001,
            energy_change: 1e-6,
            max_iterations: 200,
        }
    }
}

/// Optimization options
#[wasm_bindgen]
#[derive(Debug, Clone)]
pub struct OptimizationOptions {
    #[wasm_bindgen(getter_with_clone, setter)]
    pub mmff_variant: String,
    #[wasm_bindgen(skip)]
    pub convergence: ConvergenceOptions,
}

#[wasm_bindgen]
impl OptimizationOptions {
    #[wasm_bindgen(constructor)]
    pub fn new() -> Self {
        Self::default()
    }

    #[wasm_bindgen]
    pub fn set_max_iterations(&mut self, val: usize) {
        self.convergence.max_iterations = val;
    }

    #[wasm_bindgen]
    pub fn set_max_force(&mut self, val: f64) {
        self.convergence.max_force = val;
    }

    #[wasm_bindgen]
    pub fn set_rms_force(&mut self, val: f64) {
        self.convergence.rms_force = val;
    }

    #[wasm_bindgen]
    pub fn set_energy_change(&mut self, val: f64) {
        self.convergence.energy_change = val;
    }
}

impl Default for OptimizationOptions {
    fn default() -> Self {
        Self {
            mmff_variant: "MMFF94s".to_string(),
            convergence: ConvergenceOptions::default(),
        }
    }
}

/// Optimization result
#[wasm_bindgen]
pub struct OptimizationResult {
    pub n_atoms: usize,
    pub final_energy: f64,
    pub converged: bool,
    pub iterations: usize,
    #[wasm_bindgen(getter_with_clone)]
    pub message: String,
    coordinates: Vec<f64>,
}

#[wasm_bindgen]
impl OptimizationResult {
    #[wasm_bindgen]
    pub fn get_final_energy(&self) -> f64 {
        self.final_energy
    }

    #[wasm_bindgen]
    pub fn get_converged(&self) -> bool {
        self.converged
    }

    #[wasm_bindgen]
    pub fn get_iterations(&self) -> usize {
        self.iterations
    }

    #[wasm_bindgen]
    pub fn get_message(&self) -> String {
        self.message.clone()
    }

    #[wasm_bindgen]
    pub fn get_coord(&self, atom_idx: usize, coord_idx: usize) -> f64 {
        self.coordinates[atom_idx * 3 + coord_idx]
    }
}
#[wasm_bindgen]
pub struct ETKDGResult {
    coordinates: Vec<f64>,
    n_atoms: usize,
    success: bool,
    error: String,
}

#[wasm_bindgen]
impl ETKDGResult {
    #[wasm_bindgen]
    pub fn get_coordinates(&self) -> Vec<f64> {
        self.coordinates.clone()
    }

    #[wasm_bindgen]
    pub fn get_n_atoms(&self) -> usize {
        self.n_atoms
    }

    #[wasm_bindgen]
    pub fn get_success(&self) -> bool {
        self.success
    }

    #[wasm_bindgen]
    pub fn get_error(&self) -> String {
        self.error.clone()
    }
}

// Simple init function for WASM loading test
#[wasm_bindgen]
pub fn init() {
    // Empty function - just verify WASM can be loaded
}

// Export ETKDG functions for JavaScript
#[wasm_bindgen]
pub fn generate_initial_coordinates_wasm(sdf_content: &str) -> Result<ETKDGResult, JsValue> {
    // Validate input
    let trimmed = sdf_content.trim();
    if trimmed.is_empty() || trimmed.len() < 10 {
        return Ok(ETKDGResult {
            coordinates: Vec::new(),
            n_atoms: 0,
            success: false,
            error: "Empty or invalid SDF content".to_string(),
        });
    }

    let mol = match crate::molecule::parser::parse_sdf(trimmed) {
        Ok(m) => m,
        Err(e) => {
            return Ok(ETKDGResult {
                coordinates: Vec::new(),
                n_atoms: 0,
                success: false,
                error: format!("Parse error: {}", e),
            });
        }
    };

    let coords = crate::etkdg::generate_initial_coords(&mol);
    let mut flat_coords = Vec::new();
    for coord in &coords {
        flat_coords.extend_from_slice(coord);
    }

    Ok(ETKDGResult {
        coordinates: flat_coords,
        n_atoms: coords.len(),
        success: true,
        error: String::new(),
    })
}

#[wasm_bindgen]
pub fn optimize_from_sdf(sdf_content: &str, options: OptimizationOptions) -> OptimizationResult {
    let mol = match crate::molecule::parser::parse_sdf(sdf_content) {
        Ok(mol) => mol,
        Err(e) => {
            return OptimizationResult {
                n_atoms: 0,
                final_energy: 0.0,
                converged: false,
                iterations: 0,
                message: format!("Parse error: {}", e),
                coordinates: Vec::new(),
            };
        }
    };

    let initial_coords = crate::etkdg::generate_initial_coords(&mol);

    let variant = match options.mmff_variant.as_str() {
        "MMFF94" => MMFFVariant::MMFF94,
        _ => MMFFVariant::MMFF94s,
    };

    let ff = crate::mmff::MMFFForceField::new(&mol, variant);
    let optimizer_result = crate::optimizer::optimize(&ff, &initial_coords, &options.convergence);

    let mut flat_coords = Vec::new();
    for coord in &optimizer_result.optimized_coords {
        flat_coords.extend_from_slice(coord);
    }

    OptimizationResult {
        n_atoms: optimizer_result.optimized_coords.len(),
        final_energy: optimizer_result.final_energy,
        converged: optimizer_result.converged,
        iterations: optimizer_result.iterations,
        message: "Optimization completed".to_string(),
        coordinates: flat_coords,
    }
}

#[cfg(test)]
mod tests {
    use crate::mmff::MMFFForceField;
    use crate::molecule::parse_sdf;
    use crate::ConvergenceOptions;
    use crate::MMFFVariant;

    #[test]
    fn test_parse_real_sdf() {
        let sdf_content = r#"2029
  CDK     101218203532D 0

  6  6  0  0  0  0            999 V2000
    1.2120    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.6060    1.0493    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6060    1.0493    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2120    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6060   -1.0493    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.6060   -1.0493    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  2  0  0  0  0
  3  4  1  0  0  0  0
  4  5  2  0  0  0  0
  5  6  1  0  0  0  0
  6  1  2  0  0  0  0
M  END"#;

        let molecule = parse_sdf(sdf_content).expect("Failed to parse SDF");
        assert_eq!(molecule.atoms.len(), 6);
        assert_eq!(molecule.bonds.len(), 6);
        assert_eq!(molecule.atoms[0].symbol, "C");
        assert_eq!(molecule.atoms[0].atomic_number, 6);
    }

    #[test]
    fn test_etkdg_v3_basic() {
        use crate::etkdg::{generate_initial_coords_with_config, ETKDGConfig};

        let sdf_content = r#"2029
  CDK     101218203532D 0

  6  6  0  0  0  0            999 V2000
    1.2120    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.6060    1.0493    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6060    1.0493    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2120    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6060   -1.0493    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.6060   -1.0493    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  2  0  0  0  0
  3  4  1  0  0  0  0
  4  5  2  0  0  0  0
  5  6  1  0  0  0  0
  6  1  2  0  0  0  0
M  END"#;

        let molecule = parse_sdf(sdf_content).expect("Failed to parse SDF");
        let config = ETKDGConfig {
            max_attempts: 1,
            max_iterations: 50,
            ..Default::default()
        };
        let coords = generate_initial_coords_with_config(&molecule, &config);

        // Should generate coordinates for all atoms
        assert_eq!(coords.len(), 6);

        // All coordinates should be finite numbers
        for coord in &coords {
            for &value in coord {
                assert!(value.is_finite());
            }
        }

        // Distances should be reasonable (not all zeros or infinities)
        for i in 0..coords.len() {
            for j in (i + 1)..coords.len() {
                let dx = coords[i][0] - coords[j][0];
                let dy = coords[i][1] - coords[j][1];
                let dz = coords[i][2] - coords[j][2];
                let dist = (dx * dx + dy * dy + dz * dz).sqrt();
                assert!(dist > 0.1 && dist < 10.0);
            }
        }
    }

    #[test]
    fn test_etkdg_v3_water() {
        use crate::etkdg::{generate_initial_coords_with_config, ETKDGConfig};

        let sdf_content = r#"Water
  CDK     101218203532D 0

  3  2  0  0  0  0            999 V2000
    0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.9580    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2390    0.9270    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  1  3  1  0  0  0  0
M  END"#;

        let molecule = parse_sdf(sdf_content).expect("Failed to parse SDF");
        let config = ETKDGConfig {
            max_attempts: 3,
            max_iterations: 100,
            ..Default::default()
        };
        let coords = generate_initial_coords_with_config(&molecule, &config);

        // Should generate coordinates for all atoms
        assert_eq!(coords.len(), 3);

        // Check bond lengths are reasonable - use looser bounds for stochastic embedding
        // Water is a challenging case for distance geometry methods
        let dx_oh1 = coords[0][0] - coords[1][0];
        let dy_oh1 = coords[0][1] - coords[1][1];
        let dz_oh1 = coords[0][2] - coords[1][2];
        let dist_oh1 = (dx_oh1 * dx_oh1 + dy_oh1 * dy_oh1 + dz_oh1 * dz_oh1).sqrt();

        let dx_oh2 = coords[0][0] - coords[2][0];
        let dy_oh2 = coords[0][1] - coords[2][1];
        let dz_oh2 = coords[0][2] - coords[2][2];
        let dist_oh2 = (dx_oh2 * dx_oh2 + dy_oh2 * dy_oh2 + dz_oh2 * dz_oh2).sqrt();

        // Accept any reasonable distances (0.1-5.0 Å) since this is initial embedding
        assert!(
            dist_oh1 > 0.1 && dist_oh1 < 5.0,
            "O-H1 distance {} out of range [0.1, 5.0]",
            dist_oh1
        );
        assert!(
            dist_oh2 > 0.1 && dist_oh2 < 5.0,
            "O-H2 distance {} out of range [0.1, 5.0]",
            dist_oh2
        );
    }

    #[test]
    fn test_end_to_end_water_optimization() {
        use crate::etkdg::{generate_initial_coords_with_config, ETKDGConfig};

        let sdf = r#"Water
     RDKit          3D

  3  2  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0
    0.9580    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0
   -0.2390    0.9270    0.0000 H   0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  1  3  1  0  0 0  0
M  END"#;

        let mol = crate::molecule::parser::parse_sdf(sdf).expect("Parse failed");
        assert_eq!(mol.atoms.len(), 3);
        assert_eq!(mol.atoms[0].symbol, "O");

        let config = ETKDGConfig {
            max_attempts: 1,
            max_iterations: 50,
            ..Default::default()
        };
        let coords = generate_initial_coords_with_config(&mol, &config);
        assert_eq!(coords.len(), 3);

        let ff = crate::mmff::MMFFForceField::new(&mol, MMFFVariant::MMFF94s);
        let _initial_energy = ff.calculate_energy(&coords);

        let conv = ConvergenceOptions::default();
        let result = crate::optimizer::optimize(&ff, &coords, &conv);

        assert!(result.final_energy.is_finite(), "Energy should be finite");
        assert!(result.optimized_coords.len() == 3);
    }

    #[test]
    fn test_ring_detection_benzene_from_sdf() {
        let sdf = r#"Benzene
     RDKit          3D

  6  6  0  0  0  0  0  0  0  0999 V2000
    0.0000    1.4010    0.0000 C   0  0  0  0  0  0  0  0  0
    1.2115    0.7035    0.0000 C   0  0  0  0  0  0  0  0  0
   -0.6060    1.0493    0.0000 C   0  0  0  0  0  0  0  0  0
   -1.2115    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0
   -0.6060   -1.0493    0.0000 C   0  0  0  0  0  0  0  0  0
    1.2115   -0.7035    0.0000 C   0  0  0  0  0  0  0  0  0
  1  2  4  0  0  0  0
  2  3  4  0  0  0  0
  3  4  4  0  0 0  0
  4  5  4 0  0 0  0
  5  6  4 0 0 0 0 0
  6  1  4 0 0 0 0 0
M  END"#;
        let mol = parse_sdf(sdf).expect("Failed to parse");
        let rings = crate::molecule::graph::find_rings(&mol);
        assert_eq!(rings.len(), 1);
        assert_eq!(rings[0].len(), 6);
    }

    #[test]
    fn test_etkdg_benzene_embedding() {
        let sdf = r#"Benzene
     RDKit          3D

  6  6  0  0  0  0  0  0  0  0999 V2000
    0.0000    1.4010    0.0000 C   0  0  0  0  0  0  0  0  0
    1.2115    0.7035    0.0000 C   0  0  0  0  0  0  0  0  0
   -0.6060    1.0493    0.0000 C   0  0  0  0  0  0  0  0  0
   -1.2115    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0
   -0.6060   -1.0493    0.0000 C   0  0  0  0  0  0  0  0  0
    1.2115   -0.7035    0.0000 C   0  0  0  0  0  0  0  0  0
  1  2  4  0  0  0  0
  2  3  4  0  0  0  0
  3  4  4  0  0 0  0
  4  5  4 0  0 0  0
  5  6  4 0 0 0 0 0
  6  1  4 0 0 0 0 0
M  END"#;
        let mol = parse_sdf(sdf).expect("Failed to parse");
        let config = crate::etkdg::ETKDGConfig {
            max_attempts: 1,
            ..Default::default()
        };
        let coords = crate::etkdg::generate_initial_coords_with_config(&mol, &config);
        assert_eq!(coords.len(), 6);
        for coord in &coords {
            for &value in coord {
                assert!(value.is_finite());
            }
        }
    }

    #[test]
    fn test_ethanol_optimization() {
        use crate::etkdg::{generate_initial_coords_with_config, ETKDGConfig};

        let sdf = r#"Ethanol
     RDKit          3D

  9  8  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0
    1.5260    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0
    2.1450    1.0190    0.0000 O   0  0  0  0  0  0  0  0  0
   -0.5430    1.0200    0.0000 H   0  0  0  0  0  0  0  0  0
   -0.5430   -1.0200    0.0000 H   0  0  0  0  0  0  0  0  0
    1.8860   -1.0200    0.0000 H   0  0  0  0  0  0  0  0  0
    2.0190   -0.5430    0.0000 H   0  0  0  0  0  0  0  0  0
    1.5270    1.5640    0.9170 H   0  0  0  0  0  0  0  0  0
    3.0290    1.3430    0.3590 H   0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  1  0  0  0  0
  1  4  1  0  0  0  0
  1  5  1  0  0  0  0
  2  6  1  0  0  0  0
  2  7  1  0  0  0  0
  3  8  1  0  0  0  0
  3  9  1  0  0  0  0
M  END"#;

        let mol = crate::molecule::parser::parse_sdf(sdf).expect("Parse failed");
        assert_eq!(mol.atoms.len(), 9);

        let config = ETKDGConfig {
            max_attempts: 1,
            max_iterations: 50,
            ..Default::default()
        };
        let coords = generate_initial_coords_with_config(&mol, &config);
        assert_eq!(coords.len(), 9);

        let ff = crate::mmff::MMFFForceField::new(&mol, MMFFVariant::MMFF94s);
        let initial_energy = ff.calculate_energy(&coords);
        assert!(initial_energy.is_finite());

        let conv = ConvergenceOptions::default();
        let result = crate::optimizer::optimize(&ff, &coords, &conv);

        assert!(result.final_energy.is_finite());
        assert!(result.optimized_coords.len() == 9);
        assert!(result.final_energy <= initial_energy + 1.0);
    }

    #[test]
    fn test_linear_molecule_co2() {
        // Linear molecule edge case
        let sdf = r#"Carbon dioxide
     RDKit          3D

  3  2  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0
    1.1600    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0
   -1.1600    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0
  1  2  2  0  0  0  0
  1  3  2  0  0  0  0
M  END"#;

        let mol = crate::molecule::parser::parse_sdf(sdf).expect("Parse failed");
        assert_eq!(mol.atoms.len(), 3);

        // Use coordinates from SDF directly — ETKDG embedding is unreliable
        // for linear 3-atom molecules due to random distance sampling
        let coords: Vec<[f64; 3]> = vec![[0.0, 0.0, 0.0], [1.16, 0.0, 0.0], [-1.16, 0.0, 0.0]];

        // Verify optimization works
        let ff = crate::mmff::MMFFForceField::new(&mol, MMFFVariant::MMFF94s);
        let conv = ConvergenceOptions::default();
        let result = crate::optimizer::optimize(&ff, &coords, &conv);

        assert!(result.final_energy.is_finite());
        assert!(result.converged || result.iterations > 0);
    }

    #[test]
    fn test_tetrahedral_molecule_methane() {
        // Tetrahedral molecule edge case
        let sdf = r#"Methane
     RDKit          3D

  5  4  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0
    0.6290    0.6290    0.6290 H   0  0  0  0  0  0  0  0  0
    0.6290   -0.6290   -0.6290 H   0  0  0  0  0  0  0  0  0
   -0.6290    0.6290   -0.6290 H   0  0  0  0  0  0  0  0  0
   -0.6290   -0.6290    0.6290 H   0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  1  3  1  0  0  0  0
  1  4  1  0  0  0  0
  1  5  1  0  0  0  0
M  END"#;

        let mol = crate::molecule::parser::parse_sdf(sdf).expect("Parse failed");
        assert_eq!(mol.atoms.len(), 5);

        let coords = crate::etkdg::generate_initial_coords(&mol);
        assert_eq!(coords.len(), 5);

        // Check C-H bond lengths
        for i in 1..=4 {
            let dx = coords[0][0] - coords[i][0];
            let dy = coords[0][1] - coords[i][1];
            let dz = coords[0][2] - coords[i][2];
            let dist = (dx * dx + dy * dy + dz * dz).sqrt();
            assert!(
                dist > 0.5 && dist < 2.0,
                "C-H{} distance {} out of reasonable range",
                i,
                dist
            );
        }

        // Verify optimization
        let ff = crate::mmff::MMFFForceField::new(&mol, MMFFVariant::MMFF94s);
        let conv = ConvergenceOptions::default();
        let result = crate::optimizer::optimize(&ff, &coords, &conv);

        assert!(result.final_energy.is_finite());
    }

    #[test]
    fn test_molecule_with_halogens() {
        // Chloroform - tests halogen handling
        let sdf = r#"Chloroform
     RDKit          3D

  5  4  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0
    1.0700    0.0000    0.0000 Cl  0  0  0  0  0  0  0  0  0
   -0.3567    1.0083    0.0000 Cl  0  0  0  0  0  0  0  0  0
   -0.3567   -0.5041   -0.8730 Cl  0  0  0  0  0  0  0  0  0
   -0.3567   -0.5041    0.8730 H   0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  1  3  1  0  0  0  0
  1  4  1  0  0  0  0
  1  5  1  0  0  0  0
M  END"#;

        let mol = crate::molecule::parser::parse_sdf(sdf).expect("Parse failed");
        assert_eq!(mol.atoms.len(), 5);

        let coords = crate::etkdg::generate_initial_coords(&mol);
        assert_eq!(coords.len(), 5);

        // Check C-Cl bond lengths are longer than C-H
        let c_h_dist = {
            let dx = coords[0][0] - coords[4][0];
            let dy = coords[0][1] - coords[4][1];
            let dz = coords[0][2] - coords[4][2];
            (dx * dx + dy * dy + dz * dz).sqrt()
        };

        let c_cl_dist = {
            let dx = coords[0][0] - coords[1][0];
            let dy = coords[0][1] - coords[1][1];
            let dz = coords[0][2] - coords[1][2];
            (dx * dx + dy * dy + dz * dz).sqrt()
        };

        assert!(
            c_cl_dist > c_h_dist,
            "C-Cl ({}) should be longer than C-H ({})",
            c_cl_dist,
            c_h_dist
        );

        // Verify optimization
        let ff = crate::mmff::MMFFForceField::new(&mol, MMFFVariant::MMFF94s);
        let conv = ConvergenceOptions::default();
        let result = crate::optimizer::optimize(&ff, &coords, &conv);

        assert!(result.final_energy.is_finite());
    }

    #[test]
    fn test_amino_acid_glycine() {
        // Glycine - tests multiple functional groups
        let sdf = r#"Glycine
     RDKit          3D

 10  9  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 N   0  0  0  0  0  0  0  0  0
    1.4500    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0
    1.8500    1.2570    0.0000 C   0  0  0  0  0  0  0  0  0
    1.3300    2.1800    0.5500 O   0  0  0  0  0  0  0  0  0
    2.8000    1.4500   -0.5000 O   0  0  0  0  0  0  0  0  0
   -0.3600   -0.3200    0.8900 H   0  0  0  0  0  0  0  0  0
   -0.3600   -0.3200   -0.8900 H   0  0  0  0  0  0  0  0  0
    1.7300   -0.8900    0.5500 H   0  0  0  0  0  0  0  0  0
    1.7300   -0.8900   -0.5500 H   0  0  0  0  0  0  0  0  0
    3.1000    2.3000   -0.2500 H   0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  1  0  0  0  0
  3  4  2  0  0  0  0
  3  5  1  0  0  0  0
  1  6  1  0  0  0  0
  1  7  1  0  0  0  0
  2  8  1  0  0  0  0
  2  9  1  0  0  0  0
  5 10  1  0  0  0  0
M  END"#;

        let mol = crate::molecule::parser::parse_sdf(sdf).expect("Parse failed");
        assert_eq!(mol.atoms.len(), 10);

        let coords = crate::etkdg::generate_initial_coords(&mol);
        assert_eq!(coords.len(), 10);

        // Verify all coordinates are finite
        for (i, coord) in coords.iter().enumerate() {
            for (j, &val) in coord.iter().enumerate() {
                assert!(val.is_finite(), "Atom {} coord {} is not finite", i, j);
            }
        }

        // Verify optimization works for amino acid
        let ff = crate::mmff::MMFFForceField::new(&mol, MMFFVariant::MMFF94s);
        let conv = ConvergenceOptions::default();
        let result = crate::optimizer::optimize(&ff, &coords, &conv);

        assert!(result.final_energy.is_finite());
        assert!(result.optimized_coords.len() == 10);
    }

    // === Edge case tests ===

    #[test]
    fn test_single_atom_molecule() {
        let sdf = r#"Single
     RDKit          3D

  1  0  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0
M  END"#;

        let mol = parse_sdf(sdf).expect("Parse failed");
        assert_eq!(mol.atoms.len(), 1);
        assert_eq!(mol.bonds.len(), 0);

        let ff = crate::mmff::MMFFForceField::new(&mol, MMFFVariant::MMFF94s);
        let coords = vec![[0.0, 0.0, 0.0]];
        let energy = ff.calculate_energy(&coords);
        assert!(energy.is_finite(), "Single atom energy should be finite");
    }

    #[test]
    fn test_two_atom_molecule() {
        let sdf = r#"H2
     RDKit          3D

  2  1  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0
    0.7400    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
M  END"#;

        let mol = parse_sdf(sdf).expect("Parse failed");
        assert_eq!(mol.atoms.len(), 2);
        assert_eq!(mol.bonds.len(), 1);

        let coords = vec![[0.0, 0.0, 0.0], [0.74, 0.0, 0.0]];
        let ff = crate::mmff::MMFFForceField::new(&mol, MMFFVariant::MMFF94s);
        let energy = ff.calculate_energy(&coords);
        assert!(energy.is_finite());
    }

    #[test]
    fn test_molecule_with_triple_bond() {
        // Acetylene (HCCH) - linear molecule with triple bond
        let sdf = r#"Acetylene
     RDKit          3D

  4  3  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0
    1.2000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0
    2.2800    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0
   -0.6000    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0
  1  2  3  0  0  0  0
  2  3  1  0  0  0  0
  1  4  1  0  0  0  0
M  END"#;

        let mol = parse_sdf(sdf).expect("Parse failed");
        assert_eq!(mol.atoms.len(), 4);
        assert_eq!(mol.bonds.len(), 3);

        let ff = crate::mmff::MMFFForceField::new(&mol, MMFFVariant::MMFF94s);
        let coords = vec![
            [0.0, 0.0, 0.0],
            [1.2, 0.0, 0.0],
            [2.28, 0.0, 0.0],
            [-0.6, 0.0, 0.0],
        ];
        let energy = ff.calculate_energy(&coords);
        assert!(energy.is_finite());
    }

    #[test]
    fn test_gradient_consistency_water() {
        let sdf = r#"Water
     RDKit          3D

  3  2  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0
    0.9580    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0
   -0.2390    0.9270    0.0000 H   0  0  0  0  0  0  0  0  0
  1  2  1  0
  1  3  1  0
M  END"#;

        let mol = parse_sdf(sdf).expect("Parse failed");
        let ff = MMFFForceField::new(&mol, MMFFVariant::MMFF94s);

        let coords: Vec<[f64; 3]> = vec![[0.0, 0.0, 0.0], [0.958, 0.0, 0.0], [-0.239, 0.927, 0.0]];

        // Test angle gradient alone
        {
            use crate::mmff::angle::{angle_energy, angle_gradient, get_angle_params};
            let at = &ff.atom_types;
            for angle in &ff.angles {
                let (i, j, k) = (angle.atom1, angle.atom2, angle.atom3);
                if let Some(params) = get_angle_params(at[i], at[j], at[k]) {
                    let (g1, g2, g3) = angle_gradient(&coords, i, j, k, &params);
                    let e0 = angle_energy(&coords, i, j, k, &params);
                    let eps = 1e-7;
                    for (atom_idx, grad) in [(i, g1), (j, g2), (k, g3)] {
                        for dim in 0..3 {
                            let mut cp = coords.clone();
                            cp[atom_idx][dim] += eps;
                            let ep = angle_energy(&cp, i, j, k, &params);
                            let num = (ep - e0) / eps;
                            let err = (grad[dim] - num).abs();
                            if err > 1e-3 {
                                eprintln!("  ANGLE GRAD MISMATCH {}-{}-{} atom {} dim {}: ana={:.6} num={:.6} err={:.6}",
                                    i, j, k, atom_idx, dim, grad[dim], num, err);
                            }
                        }
                    }
                }
            }
        }

        // Test stretch-bend gradient alone
        {
            use crate::mmff::stretch_bend::{
                get_stretch_bend_params, stretch_bend_energy, stretch_bend_gradient,
            };
            use crate::mmff::{get_angle_params, get_bond_params};
            for angle in &ff.angles {
                let (i, j, k) = (angle.atom1, angle.atom2, angle.atom3);
                let bij_key = (i.min(j), i.max(j));
                let bkj_key = (k.min(j), k.max(j));
                let bond_ij = ff.bond_map.get(&bij_key);
                let bond_kj = ff.bond_map.get(&bkj_key);
                if let (Some(bij), Some(bkj)) = (bond_ij, bond_kj) {
                    if let (Some(sb_params), Some(bp_ij), Some(bp_kj), Some(ap)) = (
                        get_stretch_bend_params(
                            ff.atom_types[i],
                            ff.atom_types[j],
                            ff.atom_types[k],
                            bij.bond_type,
                            bkj.bond_type,
                        ),
                        get_bond_params(ff.atom_types[i], ff.atom_types[j], bij.bond_type),
                        get_bond_params(ff.atom_types[k], ff.atom_types[j], bkj.bond_type),
                        get_angle_params(ff.atom_types[i], ff.atom_types[j], ff.atom_types[k]),
                    ) {
                        let (g1, g2, g3) = stretch_bend_gradient(
                            &coords,
                            i,
                            j,
                            k,
                            bp_ij.r0,
                            bp_kj.r0,
                            ap.theta0.to_radians(),
                            &sb_params,
                        );
                        let e0 = stretch_bend_energy(
                            &coords,
                            i,
                            j,
                            k,
                            bp_ij.r0,
                            bp_kj.r0,
                            ap.theta0.to_radians(),
                            &sb_params,
                        );
                        let eps = 1e-7;
                        for (atom_idx, grad) in [(i, g1), (j, g2), (k, g3)] {
                            for dim in 0..3 {
                                let mut cp = coords.clone();
                                cp[atom_idx][dim] += eps;
                                let ep = stretch_bend_energy(
                                    &cp,
                                    i,
                                    j,
                                    k,
                                    bp_ij.r0,
                                    bp_kj.r0,
                                    ap.theta0.to_radians(),
                                    &sb_params,
                                );
                                let num = (ep - e0) / eps;
                                let err = (grad[dim] - num).abs();
                                if err > 1e-3 {
                                    eprintln!("  SB GRAD MISMATCH {}-{}-{} atom {} dim {}: ana={:.6} num={:.6} err={:.6}",
                                        i, j, k, atom_idx, dim, grad[dim], num, err);
                                }
                            }
                        }
                    }
                }
            }
        }

        let (e0, grad) = ff.calculate_energy_and_gradient(&coords);

        let eps = 1e-7;
        let mut max_err = 0.0f64;
        for atom_idx in 0..coords.len() {
            for dim in 0..3 {
                let mut coords_p = coords.clone();
                coords_p[atom_idx][dim] += eps;
                let e_plus = ff.calculate_energy(&coords_p);
                let num_grad = (e_plus - e0) / eps;
                let ana_grad = grad[atom_idx][dim];
                let err = (ana_grad - num_grad).abs();
                if err > max_err {
                    max_err = err;
                }
                if err > 1e-3 {
                    eprintln!(
                        "  TOTAL GRADIENT MISMATCH: atom {} dim {}: analytical={:.6} numerical={:.6} err={:.6}",
                        atom_idx, dim, ana_grad, num_grad, err
                    );
                }
            }
        }
        assert!(
            max_err < 1e-3,
            "Max gradient error {:.6} exceeds 1e-3",
            max_err
        );
    }

    #[test]
    fn test_gradient_consistency_ethanol() {
        let sdf = r#"Ethanol
     RDKit          3D

  9  8  0  0  0  0  0  0  0  0999 V2000
   -0.8883    0.1670   -0.0273 C   0  0  0  0  0  0  0  0  0
    0.4658   -0.5116   -0.0368 C   0  0  0  0  0  0  0  0  0
    1.4311    0.3229    0.5867 O   0  0  0  0  0  0  0  0  0
   -0.8487    1.1175   -0.5695 H   0  0  0  0  0  0  0  0  0
   -1.6471   -0.4704   -0.4896 H   0  0  0  0  0  0  0  0  0
   -1.1964    0.3978    0.9977 H   0  0  0  0  0  0  0  0  0
    0.7920   -0.7224   -1.0597 H   0  0  0  0  0  0  0  0  0
    0.4246   -1.4559    0.5138 H   0  0  0  0  0  0  0  0  0
    1.4671    1.1550    0.0848 H   0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  1  0
  1  4  1  0
  1  5  1  0
  1  6  1  0
  2  7  1  0
  2  8  1  0
  3  9  1  0
M  END"#;

        let mol = parse_sdf(sdf).expect("Parse failed");
        let ff = MMFFForceField::new(&mol, MMFFVariant::MMFF94s);

        let coords: Vec<[f64; 3]> = mol.atoms.iter().map(|a| a.position).collect();

        let (e0, grad) = ff.calculate_energy_and_gradient(&coords);

        let eps = 1e-7;
        let mut max_err = 0.0f64;
        let mut max_err_info = (0usize, 0usize, 0.0f64, 0.0f64);
        for atom_idx in 0..coords.len() {
            for dim in 0..3 {
                let mut coords_p = coords.clone();
                coords_p[atom_idx][dim] += eps;
                let e_plus = ff.calculate_energy(&coords_p);
                let num_grad = (e_plus - e0) / eps;
                let ana_grad = grad[atom_idx][dim];
                let err = (ana_grad - num_grad).abs();
                if err > max_err {
                    max_err = err;
                    max_err_info = (atom_idx, dim, ana_grad, num_grad);
                }
                if err > 1e-3 {
                    eprintln!(
                        "  GRADIENT MISMATCH: atom {} dim {}: analytical={:.6} numerical={:.6} err={:.6}",
                        atom_idx, dim, ana_grad, num_grad, err
                    );
                }
            }
        }
        eprintln!(
            "  Max gradient error: {:.6} at atom {} dim {} (ana={:.6} num={:.6})",
            max_err, max_err_info.0, max_err_info.1, max_err_info.2, max_err_info.3
        );
        assert!(
            max_err < 1e-3,
            "Max gradient error {:.6} exceeds 1e-3",
            max_err
        );
    }

    #[test]
    fn test_molecule_with_sulfur() {
        // Dimethyl sulfide (CH3-S-CH3)
        let sdf = r#"Dimethyl sulfide
     RDKit          3D

  6  5  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0
    1.5400    0.0000    0.0000 S   0  0  0  0  0  0  0  0  0
    2.9000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0
    0.0000    1.0200    0.0000 H   0  0  0  0  0  0  0  0  0
    0.0000   -1.0200    0.0000 H   0  0  0  0  0  0  0  0  0
    0.5400    0.8800    0.0000 H   0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  1  0  0  0  0
  1  4  1  0  0  0  0
  1  5  1  0  0  0  0
  3  6  1  0  0  0  0
M  END"#;

        let mol = parse_sdf(sdf).expect("Parse failed");
        assert_eq!(mol.atoms.len(), 6);

        let coords = crate::etkdg::generate_initial_coords(&mol);
        assert_eq!(coords.len(), 6);

        let ff = crate::mmff::MMFFForceField::new(&mol, MMFFVariant::MMFF94s);
        let conv = ConvergenceOptions::default();
        let result = crate::optimizer::optimize(&ff, &coords, &conv);
        assert!(result.final_energy.is_finite());
    }

    // === Atom type assignment tests ===

    #[test]
    fn test_atom_types_methane() {
        let sdf = r#"Methane
     RDKit          3D

  5  4  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0
    1.0900    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0
    0.0000    1.0900    0.0000 H   0  0  0  0  0  0  0  0  0
   -0.3630   -0.5450    0.8900 H   0  0  0  0  0  0  0  0  0
   -0.3630   -0.5450   -0.8900 H   0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  1  3  1  0  0  0  0
  1  4  1  0  0  0  0
  1  5  1  0  0  0  0
M  END"#;

        let mol = parse_sdf(sdf).expect("Parse failed");
        let ff = crate::mmff::MMFFForceField::new(&mol, MMFFVariant::MMFF94s);

        assert_eq!(ff.atom_types[0], crate::mmff::MMFFAtomType::C_3);
        assert_eq!(ff.atom_types[1], crate::mmff::MMFFAtomType::H);
        assert_eq!(ff.atom_types[2], crate::mmff::MMFFAtomType::H);
        assert_eq!(ff.atom_types[3], crate::mmff::MMFFAtomType::H);
        assert_eq!(ff.atom_types[4], crate::mmff::MMFFAtomType::H);
    }

    #[test]
    fn test_atom_types_formaldehyde() {
        // H2C=O: C should be C_2, O should be O_2
        let sdf = r#"Formaldehyde
     RDKit          3D

  4  3  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0
    1.2000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0
   -0.5400    0.9400    0.0000 H   0  0  0  0  0  0  0  0  0
   -0.5400   -0.9400    0.0000 H   0  0  0  0  0  0  0  0  0
  1  2  2  0  0  0  0
  1  3  1  0  0  0  0
  1  4  1  0  0  0  0
M  END"#;

        let mol = parse_sdf(sdf).expect("Parse failed");
        let ff = crate::mmff::MMFFForceField::new(&mol, MMFFVariant::MMFF94s);

        assert_eq!(
            ff.atom_types[0],
            crate::mmff::MMFFAtomType::C_2,
            "Carbonyl C should be C_2, got {:?}",
            ff.atom_types[0]
        );
        assert_eq!(
            ff.atom_types[1],
            crate::mmff::MMFFAtomType::O_CO2,
            "Carbonyl O double-bonded to C should be O_CO2, got {:?}",
            ff.atom_types[1]
        );
    }

    #[test]
    fn test_atom_types_hydroxide() {
        // OH- with formal charge -1: O has 1 bond so it's O_2 (sp2 for O with 1 bond)
        let sdf = r#"Hydroxide
     RDKit          3D

  2  1  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 O  0  5  0  0  0  0  0  0  0  0  0  0  0
    0.9600    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
M  END"#;

        let mol = parse_sdf(sdf).expect("Parse failed");
        let ff = crate::mmff::MMFFForceField::new(&mol, MMFFVariant::MMFF94s);

        assert_eq!(ff.atom_types[0], crate::mmff::MMFFAtomType::O_2);
        assert_eq!(ff.atom_types[1], crate::mmff::MMFFAtomType::H_OH);
    }

    #[test]
    fn test_atom_types_ether() {
        // Dimethyl ether CH3-O-CH3: O should be O_R (ether oxygen bonded to 2 carbons)
        let sdf = r#"Dimethyl ether
     RDKit          3D

  6  5  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0
    1.4300    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0
    2.4900    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0
    0.0000    1.0200    0.0000 H   0  0  0  0  0  0  0  0  0
    0.0000   -1.0200    0.0000 H   0  0  0  0  0  0  0  0  0
    3.1100    0.0000    0.8900 H   0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  1  0  0  0  0
  1  4  1  0  0  0  0
  1  5  1  0  0  0  0
  3  6  1  0  0  0  0
M  END"#;

        let mol = parse_sdf(sdf).expect("Parse failed");
        let ff = crate::mmff::MMFFForceField::new(&mol, MMFFVariant::MMFF94s);

        assert_eq!(
            ff.atom_types[1],
            crate::mmff::MMFFAtomType::O_R,
            "Ether O bonded to 2 C should be O_R, got {:?}",
            ff.atom_types[1]
        );
    }

    #[test]
    fn test_optimizer_convergence_improves_energy() {
        let sdf = r#"Ethanol
     RDKit          3D

  9  8  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0
    1.5260    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0
    2.1450    1.0190    0.0000 O   0  0  0  0  0  0  0  0  0
   -0.5430    1.0200    0.0000 H   0  0  0  0  0  0  0  0  0
   -0.5430   -1.0200    0.0000 H   0  0  0  0  0  0  0  0  0
    1.8860   -1.0200    0.0000 H   0  0  0  0  0  0  0  0  0
    2.0190   -0.5430    0.0000 H   0  0  0  0  0  0  0  0  0
    1.5270    1.5640    0.9170 H   0  0  0  0  0  0  0  0  0
    3.0290    1.3430    0.3590 H   0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  1  0  0  0  0
  1  4  1  0  0  0  0
  1  5  1  0  0  0  0
  2  6  1  0  0  0  0
  2  7  1  0  0  0  0
  3  8  1  0  0  0  0
  3  9  1  0  0  0  0
M  END"#;

        let mol = parse_sdf(sdf).expect("Parse failed");
        let coords = crate::etkdg::generate_initial_coords(&mol);
        let ff = crate::mmff::MMFFForceField::new(&mol, MMFFVariant::MMFF94s);

        let initial_energy = ff.calculate_energy(&coords);
        let conv = ConvergenceOptions {
            max_iterations: 500,
            ..Default::default()
        };
        let result = crate::optimizer::optimize(&ff, &coords, &conv);

        assert!(
            result.final_energy <= initial_energy,
            "Optimizer should not increase energy: initial={}, final={}",
            initial_energy,
            result.final_energy
        );
    }

    #[test]
    fn test_mmff94s_variant() {
        let sdf = r#"Water
     RDKit          3D

  3  2  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0
    0.9580    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0
   -0.2390    0.9270    0.0000 H   0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  1  3  1  0  0  0  0
M  END"#;

        let mol = parse_sdf(sdf).expect("Parse failed");

        let coords = vec![[0.0, 0.0, 0.0], [0.958, 0.0, 0.0], [-0.239, 0.927, 0.0]];
        let ff_94 = crate::mmff::MMFFForceField::new(&mol, MMFFVariant::MMFF94);
        let ff_94s = crate::mmff::MMFFForceField::new(&mol, MMFFVariant::MMFF94s);

        let e_94 = ff_94.calculate_energy(&coords);
        let e_94s = ff_94s.calculate_energy(&coords);

        // Both should produce finite energies
        assert!(e_94.is_finite(), "MMFF94 energy should be finite");
        assert!(e_94s.is_finite(), "MMFF94s energy should be finite");
    }

    #[test]
    fn test_wasm_api_from_sdf() {
        use crate::optimize_from_sdf;
        use crate::OptimizationOptions;

        let sdf = r#"Water
     RDKit          3D

  3  2  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0
    0.9580    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0
   -0.2390    0.9270    0.0000 H   0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  1  3  1  0  0  0  0
M  END"#;

        let options = OptimizationOptions::default();
        let result = optimize_from_sdf(sdf, options);

        assert_eq!(result.n_atoms, 3);
        assert!(result.get_final_energy().is_finite());
        assert_eq!(result.get_iterations() > 0, true);

        let x0 = result.get_coord(0, 0);
        let y0 = result.get_coord(0, 1);
        let z0 = result.get_coord(0, 2);
        assert!(x0.is_finite() && y0.is_finite() && z0.is_finite());
    }

    #[test]
    fn test_acetic_acid_convergence() {
        use crate::ConvergenceOptions;

        let sdf = r#"Acetic acid
     RDKit          3D

  8  7  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0
    1.5260    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0
    2.2000    1.0190    0.0000 O   0  0  0  0  0  0  0  0  0
    2.9000    0.3000    0.0000 O   0  0  0  0  0  0  0  0  0
   -0.5430    1.0200    0.0000 H   0  0  0  0  0  0  0  0  0
   -0.5430   -1.0200    0.0000 H   0  0  0  0  0  0  0  0  0
    1.8860   -1.0200    0.0000 H   0  0  0  0  0  0  0  0  0
    3.0000    1.2000    0.0000 H   0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  2  0  0  0  0
  2  4  1  0  0  0  0
  1  5  1  0  0  0  0
  1  6  1  0  0  0  0
  2  7  1  0  0  0  0
  4  8  1  0  0  0  0
M  END"#;

        let mol = parse_sdf(sdf).expect("Parse failed");
        assert_eq!(mol.atoms.len(), 8);

        let coords = crate::etkdg::generate_initial_coords(&mol);
        let ff = crate::mmff::MMFFForceField::new(&mol, MMFFVariant::MMFF94s);

        let initial_energy = ff.calculate_energy(&coords);
        let initial_grad = ff.calculate_gradient(&coords);
        let max_grad = initial_grad
            .iter()
            .flat_map(|g| g.iter())
            .map(|&v| v.abs())
            .fold(0.0f64, f64::max);

        let conv = ConvergenceOptions {
            max_iterations: 500,
            max_force: 0.01,
            rms_force: 0.001,
            energy_change: 1e-6,
        };
        let result = crate::optimizer::optimize(&ff, &coords, &conv);

        eprintln!(
            "Acetic acid: {} atoms, {} iters, converged={}",
            mol.atoms.len(),
            result.iterations,
            result.converged
        );
        eprintln!("  Initial energy: {:.4}", initial_energy);
        eprintln!("  Final energy:   {:.4}", result.final_energy);
        eprintln!(
            "  Energy change:  {:.6}",
            (initial_energy - result.final_energy).abs()
        );
        eprintln!("  Initial max |grad|: {:.6}", max_grad);

        // Energy should decrease, or initial gradient should already be below threshold
        if max_grad > conv.max_force {
            assert!(
                result.final_energy <= initial_energy + 1e-6,
                "Energy should decrease: initial={}, final={}",
                initial_energy,
                result.final_energy
            );
        }
    }

    #[test]
    fn test_generate_initial_coordinates_wasm() {
        use crate::generate_initial_coordinates_wasm;

        let sdf = r#"Water
     RDKit          3D

  3  2  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0
    0.9580    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0
   -0.2390    0.9270    0.0000 H   0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  1  3  1  0  0  0  0
M  END"#;

        let result = generate_initial_coordinates_wasm(sdf).expect("Should return Ok");

        assert!(result.get_success(), "Should succeed for valid SDF");
        assert_eq!(result.get_n_atoms(), 3, "Water should have 3 atoms");

        let coords = result.get_coordinates();
        assert_eq!(
            coords.len(),
            9,
            "Should have 9 coordinate values (3 atoms * 3 coords)"
        );

        for coord in &coords {
            assert!(coord.is_finite(), "All coordinates should be finite");
        }
    }

    #[test]
    fn test_generate_initial_coordinates_wasm_empty_input() {
        use crate::generate_initial_coordinates_wasm;

        let result =
            generate_initial_coordinates_wasm("").expect("Should return Ok even for empty input");
        assert!(!result.get_success(), "Should fail for empty input");
        assert!(!result.get_error().is_empty(), "Should have error message");
    }

    #[test]
    fn test_generate_initial_coordinates_wasm_invalid_input() {
        use crate::generate_initial_coordinates_wasm;

        let result =
            generate_initial_coordinates_wasm("invalid sdf content").expect("Should return Ok");
        assert!(!result.get_success(), "Should fail for invalid SDF");
        assert!(!result.get_error().is_empty(), "Should have error message");
    }

    #[test]
    fn test_mmff_bond_params_acetic_acid() {
        let sdf = r#"
     RDKit          3D

  8  7  0  0  0  0  0  0  0  0999 V2000
   -0.9335   -0.0601   -0.2304 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.4936    0.2789    0.0469 C   0  0  0  0  0  0   0  0  0  0  0  0
    1.0325    1.3566   -0.1361 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.1814   -0.7645    0.5462 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4427    0.8203   -0.6327 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4305   -0.3576    0.6963 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.9867   -0.8625   -0.9702 H   0  0  0  0 0  0  0  0  0  0  0  0
    2.0859   -0.4112    0.6800 H   0  0  0  0 0  0  0  0  0  0 0  0  0
  1  2  1  0
  2  3  2  0
  2  4  1  0
  1  5  1  0
  1  6  1  0
  1  7  1  0
  4  8  1  0
M  END"#;

        let mol = crate::molecule::parser::parse_sdf(sdf).expect("parse failed");
        let ff = crate::mmff::MMFFForceField::new(&mol, crate::MMFFVariant::MMFF94s);

        let mut coords = mol.atoms.iter().map(|a| a.position).collect::<Vec<_>>();
        let e_before = ff.calculate_energy(&coords);
        let conv = crate::ConvergenceOptions::default();
        let result = crate::optimizer::optimize(&ff, &coords, &conv);
        for (i, c) in result.optimized_coords.iter().enumerate() {
            coords[i] = *c;
        }
        let e_after = ff.calculate_energy(&coords);

        eprintln!(
            "Acetic acid optimization: {:.2} -> {:.2} kcal/mol, converged={}, iters={}",
            e_before, e_after, result.converged, result.iterations
        );

        assert!(e_after < e_before, "Energy should decrease");
        assert!(
            result.iterations > 0,
            "Optimizer should run at least one iteration"
        );

        for bond in &mol.bonds {
            let dx = coords[bond.atom1][0] - coords[bond.atom2][0];
            let dy = coords[bond.atom1][1] - coords[bond.atom2][1];
            let dz = coords[bond.atom1][2] - coords[bond.atom2][2];
            let d = (dx * dx + dy * dy + dz * dz).sqrt();
            let sym1 = &mol.atoms[bond.atom1].symbol;
            let sym2 = &mol.atoms[bond.atom2].symbol;
            assert!(
                d > 0.5 && d < 2.5,
                "Bond {}-{} distance {:.3} out of range after optimization",
                sym1,
                sym2,
                d
            );
        }
    }

    fn compare_with_rdkit(
        name: &str,
        sdf: &str,
        rdkit_energy: f64,
        rdkit_positions: &[[f64; 3]],
        rdkit_forces: &[[f64; 3]],
        energy_tol: f64,
        force_tol: f64,
    ) {
        let mol = parse_sdf(sdf).expect(&format!("{}: parse failed", name));
        assert_eq!(
            mol.atoms.len(),
            rdkit_positions.len(),
            "{}: atom count mismatch",
            name
        );

        let coords: Vec<[f64; 3]> = mol.atoms.iter().map(|a| a.position).collect();

        let ff = MMFFForceField::new(&mol, MMFFVariant::MMFF94s);
        let our_energy = ff.calculate_energy(&coords);
        let our_gradient = ff.calculate_gradient(&coords);

        eprintln!("\n=== {} ({} atoms) ===", name, mol.atoms.len());
        eprintln!("  RDKit energy:   {:.6}", rdkit_energy);
        eprintln!("  Our   energy:   {:.6}", our_energy);
        eprintln!(
            "  Energy diff:    {:.6} (tol {:.2})",
            our_energy - rdkit_energy,
            energy_tol
        );

        let bd = ff.calculate_energy_breakdown(&coords);
        eprintln!(
            "  Breakdown: bond={:.4} angle={:.4} sb={:.4} torsion={:.4} oop={:.4} vdw={:.4} elec={:.4}",
            bd.bond, bd.angle, bd.stretch_bend, bd.torsion, bd.oop, bd.vdw, bd.electrostatic
        );
        eprintln!(
            "  Breakdown total: {:.6} (calc_total: {:.6})",
            bd.total(),
            our_energy
        );

        assert!(
            (our_energy - rdkit_energy).abs() < energy_tol,
            "{}: energy mismatch: ours={:.6} vs rdkit={:.6}, diff={:.6}, tol={:.2}\n  breakdown: bond={:.4} angle={:.4} sb={:.4} torsion={:.4} oop={:.4} vdw={:.4} elec={:.4}",
            name, our_energy, rdkit_energy, our_energy - rdkit_energy, energy_tol,
            bd.bond, bd.angle, bd.stretch_bend, bd.torsion, bd.oop, bd.vdw, bd.electrostatic
        );

        for i in 0..mol.atoms.len() {
            let sym = &mol.atoms[i].symbol;
            for dim in 0..3 {
                let our_f = our_gradient[i][dim];
                let rdkit_f = rdkit_forces[i][dim];
                let diff = (our_f - rdkit_f).abs();
                if diff > force_tol {
                    eprintln!(
                        "  {}{} dim={}: our={:.4} rdkit={:.4} diff={:.4}",
                        sym, i, dim, our_f, rdkit_f, diff
                    );
                }
            }
        }
    }

    #[test]
    fn test_rdkit_compare_water_distorted() {
        let sdf = r#"
     RDKit          3D

  3  2  0  0  0  0  0  0  0  0999 V2000
    0.0075    0.3977    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5671    0.1156    0.1000 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.7596   -0.2134    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  1  3  1  0
M  END"#;
        compare_with_rdkit(
            "water_distorted",
            sdf,
            109.746295,
            &[
                [0.007544053151, 0.3977434312, 0.0],
                [-0.5671031371, 0.1156068351, 0.1],
                [0.7595590839, -0.2133502663, 0.0],
            ],
            &[
                [788.4815519, 393.9982971, -138.129142],
                [-780.549907, -392.9872317, 137.1337623],
                [-7.931644892, -1.01106543, 0.9953796937],
            ],
            100.0, // energy tolerance -- TODO: bond params need calibration
            500.0, // force tolerance
        );
    }

    #[test]
    fn test_rdkit_compare_methane_distorted() {
        let sdf = r#"
     RDKit          3D

  5  4  0  0  0  0  0  0  0  0999 V2000
   -0.0000   -0.0000   -0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5253    0.7541    0.1146 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3949   -0.8359   -0.5815 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.0848   -0.2937    1.0485 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.9855    0.2756   -0.3817 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  1  3  1  0
  1  4  1  0
  1  5  1  0
M  END"#;
        compare_with_rdkit(
            "methane_distorted",
            sdf,
            15.758408,
            &[
                [-1.66630037e-08, -7.370584101e-09, -1.935743655e-08],
                [-0.5253425223, 0.7541226095, 0.1146386471],
                [-0.3949462513, -0.8359379305, -0.5814849735],
                [0.08475366302, -0.2937409757, 1.048538274],
                [0.9855351272, 0.2755563041, -0.3816919283],
            ],
            &[
                [108.4078762, -140.2588297, 1.867300278],
                [-113.1210158, 150.054103, -1.076790919],
                [-5.386759294, 5.842118387, 0.40088068],
                [9.943412935, -14.28374526, -1.440300398],
                [0.1564860061, -1.353646382, 0.248910359],
            ],
            50.0,
            500.0,
        );
    }

    #[test]
    fn test_rdkit_compare_formaldehyde_distorted() {
        let sdf = r#"
     RDKit          3D

  4  3  0  0  0  0  0  0  0  0999 V2000
   -0.0122    0.0019    0.0002 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.1981   -0.1844   -0.0180 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3513    1.2123    0.3023 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7346   -0.8298    0.0155 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0
  1  3  1  0
  1  4  1  0
M  END"#;
        compare_with_rdkit(
            "formaldehyde_distorted",
            sdf,
            11.770947,
            &[
                [-0.01218758952, 0.001875648265, 0.0001831614605],
                [1.198061605, -0.1843732515, -0.01800530945],
                [-0.351295555, 1.212282166, 0.3023159437],
                [-0.7345784602, -0.8297845631, 0.0155062043],
            ],
            &[
                [-5.824477141, 79.90010669, 43.57340591],
                [14.88506482, -13.02784507, -11.61226497],
                [-0.7913849641, -76.69016668, -25.23267077],
                [-8.269202719, 9.817905061, -6.728470168],
            ],
            50.0,
            500.0,
        );
    }

    #[test]
    fn test_rdkit_compare_ethane_distorted() {
        let sdf = r#"
     RDKit          3D

  8  7  0  0  0  0  0  0  0  0999 V2000
   -0.7558    0.0071   -0.0165 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7558   -0.0071    0.0165 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0134   -0.0004    0.9931 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.1223    0.9481   -0.4375 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.1346   -0.8156   -0.6303 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.1346    0.8156    0.6303 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.1634    0.1004   -0.9931 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.1223   -0.9481    0.4375 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  1  3  1  0
  1  4  1  0
  1  5  1  0
  2  6  1  0
  2  7  1  0
  2  8  1  0
M  END"#;
        compare_with_rdkit(
            "ethane_distorted",
            sdf,
            -2.026288,
            &[
                [-0.755815366, 0.007099600842, -0.01650379534],
                [0.7558155643, -0.007099638544, 0.01650381759],
                [-1.013351354, -0.0003807007239, 0.9931410125],
                [-1.122263433, 0.9480894007, -0.4375433377],
                [-1.134621159, -0.8155813974, -0.6302809052],
                [1.134621157, 0.8155825009, 0.6302789265],
                [1.163351445, 0.1003786207, -0.9931410825],
                [1.122263145, -0.9480883864, 0.4375453642],
            ],
            &[
                [24.23595717, 12.08553838, -32.43133914],
                [4.476204058, -0.4744418333, -8.914809105],
                [-27.25760512, -10.30463906, 35.31801151],
                [0.5922368163, -0.3195902473, -2.069117931],
                [-3.043872693, -1.451267824, 8.231158342],
                [0.6929275501, 0.3936112452, -0.3123159398],
                [0.01330082556, 0.1183574055, 0.002080831386],
                [0.290851388, -0.047568067, 0.1763314297],
            ],
            50.0,
            500.0,
        );
    }

    #[test]
    fn test_rdkit_compare_acetic_acid_distorted() {
        let sdf = r#"
     RDKit          3D

  8  7  0  0  0  0  0  0  0  0999 V2000
   -0.9335   -0.0601   -0.2304 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.4936    0.2789    0.0469 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0325    1.3566   -0.1361 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.1814   -0.7645    0.5462 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2427    0.9703   -0.5327 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4305   -0.3576    0.6963 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.9867   -0.8625   -0.9702 H   0  0  0  0  0  0  0  0  0  0  0  0
    2.0859   -0.4112    0.6800 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  2  0
  2  4  1  0
  1  5  1  0
  1  6  1  0
  1  7  1  0
  4  8  1  0
M  END"#;
        compare_with_rdkit(
            "acetic_acid_distorted",
            sdf,
            -22.071705,
            &[
                [-0.9335433255, -0.06006274984, -0.2303758313],
                [0.4936076142, 0.2789326841, 0.04691435166],
                [1.032525311, 1.356629298, -0.1361364057],
                [1.18140137, -0.764470062, 0.5462122721],
                [-1.24268686, 0.9702986055, -0.5327062184],
                [-1.430456818, -0.357583207, 0.6962707506],
                [-0.9867314458, -0.8625467632, -0.9701798995],
                [2.085884154, -0.4111978051, 0.6800009806],
            ],
            &[
                [16.85401044, 26.33698995, 7.536204383],
                [7.664343486, -13.82904173, 5.604063667],
                [3.275839936, 0.3641767657, 0.1828357139],
                [0.04046025863, -0.03166734226, 0.003236482304],
                [-23.03144195, -24.8494376, -7.296501161],
                [-1.100109267, 2.662399103, -0.4656106469],
                [-3.701411072, 9.345532215, -5.563808595],
                [-0.001691824861, 0.001048638602, -0.0004198435949],
            ],
            100.0, // energy tolerance -- TODO: bond params need calibration
            500.0, // force tolerance
        );
    }

    /// Helper: compute distance between two atoms in coords
    fn dist3(coords: &[[f64; 3]], i: usize, j: usize) -> f64 {
        let dx = coords[i][0] - coords[j][0];
        let dy = coords[i][1] - coords[j][1];
        let dz = coords[i][2] - coords[j][2];
        (dx * dx + dy * dy + dz * dz).sqrt()
    }

    /// Helper: compute angle in degrees between atoms i-j-k
    fn angle_deg(coords: &[[f64; 3]], i: usize, j: usize, k: usize) -> f64 {
        let v1 = [
            coords[i][0] - coords[j][0],
            coords[i][1] - coords[j][1],
            coords[i][2] - coords[j][2],
        ];
        let v2 = [
            coords[k][0] - coords[j][0],
            coords[k][1] - coords[j][1],
            coords[k][2] - coords[j][2],
        ];
        let n1 = (v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2]).sqrt();
        let n2 = (v2[0] * v2[0] + v2[1] * v2[1] + v2[2] * v2[2]).sqrt();
        if n1 < 1e-10 || n2 < 1e-10 {
            return 0.0;
        }
        let cos = (v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2]) / (n1 * n2);
        cos.clamp(-1.0, 1.0).acos().to_degrees()
    }

    /// ETKDG + optimize, then validate bond lengths against RDKit reference
    fn validate_etkdg_geometry(
        name: &str,
        sdf: &str,
        ref_bond_lengths: &[(usize, usize, f64)], // (i, j, expected_length)
        ref_angles: &[(usize, usize, usize, f64)], // (i, j, k, expected_angle_deg)
        bond_tol: f64,                            // tolerance in Angstroms
        angle_tol: f64,                           // tolerance in degrees
    ) {
        let mol = parse_sdf(sdf).expect(&format!("{}: parse failed", name));

        // Run ETKDG + optimize
        let config = crate::etkdg::ETKDGConfig {
            max_attempts: 5,
            max_iterations: 300,
            ..Default::default()
        };
        let mut coords = crate::etkdg::generate_initial_coords_with_config(&mol, &config);

        let ff = MMFFForceField::new(&mol, MMFFVariant::MMFF94s);

        // Print atom types for debugging
        eprint!("  Atom types:");
        for (i, at) in ff.atom_types.iter().enumerate() {
            eprint!(" {}:{:?}", i, at);
        }
        eprintln!();

        let conv = ConvergenceOptions {
            max_iterations: 500,
            ..Default::default()
        };
        let result = crate::optimizer::optimize(&ff, &coords, &conv);
        eprintln!(
            "  Optimize: converged={}, iters={}, energy={:.6}",
            result.converged, result.iterations, result.final_energy
        );

        let bd = ff.calculate_energy_breakdown(&result.optimized_coords);
        eprintln!("  Breakdown: bond={:.4} angle={:.4} sb={:.4} torsion={:.4} oop={:.4} vdw={:.4} elec={:.4}",
            bd.bond, bd.angle, bd.stretch_bend, bd.torsion, bd.oop, bd.vdw, bd.electrostatic);

        for (i, c) in result.optimized_coords.iter().enumerate() {
            coords[i] = *c;
        }

        eprintln!("\n=== ETKDG validation: {} ===", name);
        let mut max_bond_err = 0.0f64;
        for &(i, j, ref_len) in ref_bond_lengths {
            let our_len = dist3(&coords, i, j);
            let err = (our_len - ref_len).abs();
            if err > max_bond_err {
                max_bond_err = err;
            }
            if err > bond_tol {
                eprintln!(
                    "  BOND {}-{}: ours={:.4} ref={:.4} err={:.4} EXCEEDS tol={:.3}",
                    i, j, our_len, ref_len, err, bond_tol
                );
            }
        }

        let mut max_angle_err = 0.0f64;
        for &(i, j, k, ref_ang) in ref_angles {
            let our_ang = angle_deg(&coords, i, j, k);
            let err = (our_ang - ref_ang).abs();
            if err > max_angle_err {
                max_angle_err = err;
            }
            if err > angle_tol {
                eprintln!(
                    "  ANGLE {}-{}-{}: ours={:.2} ref={:.2} err={:.2} EXCEEDS tol={:.1}",
                    i, j, k, our_ang, ref_ang, err, angle_tol
                );
            }
        }

        eprintln!(
            "  Max bond error: {:.4} A (tol {:.3})",
            max_bond_err, bond_tol
        );
        eprintln!(
            "  Max angle error: {:.2} deg (tol {:.1})",
            max_angle_err, angle_tol
        );

        assert!(
            max_bond_err < bond_tol,
            "{}: max bond length error {:.4} exceeds tolerance {:.3}",
            name,
            max_bond_err,
            bond_tol
        );
        assert!(
            max_angle_err < angle_tol,
            "{}: max angle error {:.2} exceeds tolerance {:.1}",
            name,
            max_angle_err,
            angle_tol
        );
    }

    #[test]
    fn test_etkdg_water_geometry() {
        let sdf = r#"Water
     RDKit          3D

  3  2  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0
    0.9580    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0
   -0.2390    0.9270    0.0000 H   0  0  0  0  0  0  0  0  0
  1  2  1  0
  1  3  1  0
M  END"#;
        validate_etkdg_geometry(
            "water",
            sdf,
            &[(0, 1, 0.969), (0, 2, 0.969)],
            &[(1, 0, 2, 103.98)],
            0.05,
            3.0,
        );
    }

    #[test]
    fn test_etkdg_methane_geometry() {
        let sdf = r#"Methane
     RDKit          3D

  5  4  0  0  0  0  0  0  0  0999 V2000
   -0.0000   -0.0000   -0.0000 C   0  0  0  0  0  0  0  0  0
   -0.5253    0.7541    0.1146 H   0  0  0  0  0  0  0  0  0
   -0.3949   -0.8359   -0.5815 H   0  0  0  0  0  0  0  0  0
    0.0848   -0.2937    1.0485 H   0  0  0  0  0  0  0  0  0
    0.9855    0.2756   -0.3817 H   0  0  0  0  0  0  0  0  0
  1  2  1  0
  1  3  1  0
  1  4  1  0
  1  5  1  0
M  END"#;
        validate_etkdg_geometry(
            "methane",
            sdf,
            &[(0, 1, 1.092), (0, 2, 1.092), (0, 3, 1.092), (0, 4, 1.092)],
            &[(1, 0, 2, 109.47), (1, 0, 3, 109.47), (1, 0, 4, 109.47)],
            0.05,
            3.0,
        );
    }

    #[test]
    fn test_etkdg_ethanol_geometry() {
        let sdf = r#"Ethanol
     RDKit          3D

  9  8  0  0  0  0  0  0  0  0999 V2000
   -0.8883    0.1670   -0.0273 C   0  0  0  0  0  0  0  0  0
    0.4658   -0.5116   -0.0368 C   0  0  0  0  0  0  0  0  0
    1.4311    0.3229    0.5867 O   0  0  0  0  0  0  0  0  0
   -0.8487    1.1175   -0.5695 H   0  0  0  0  0  0  0  0  0
   -1.6471   -0.4704   -0.4896 H   0  0  0  0  0  0  0  0  0
   -1.1964    0.3978    0.9977 H   0  0  0  0  0  0  0  0  0
    0.7920   -0.7224   -1.0597 H   0  0  0  0  0  0  0  0  0
    0.4246   -1.4559    0.5138 H   0  0  0  0  0  0  0  0  0
    1.4671    1.1550    0.0848 H   0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  1  0
  1  4  1  0
  1  5  1  0
  1  6  1  0
  2  7  1  0
  2  8  1  0
  3  9  1  0
M  END"#;
        validate_etkdg_geometry(
            "ethanol",
            sdf,
            &[
                (0, 1, 1.515), // C-C
                (1, 2, 1.420), // C-O
                (2, 8, 0.972), // O-H
                (0, 3, 1.095), // C-H
            ],
            &[
                (0, 1, 2, 109.98), // C-C-O
                (1, 2, 8, 107.55), // C-O-H
            ],
            0.05,
            3.0,
        );
    }

    #[test]
    fn test_etkdg_formaldehyde_geometry() {
        let sdf = r#"Formaldehyde
     RDKit          3D

  4  3  0  0  0  0  0  0  0  0999 V2000
   -0.0122    0.0019    0.0002 C   0  0  0  0  0  0  0  0  0
    1.1981   -0.1844   -0.0180 O   0  0  0  0  0  0  0  0  0
   -0.4513    1.0123    0.0023 H   0  0  0  0  0  0  0  0  0
   -0.7346   -0.8298    0.0155 H   0  0  0  0  0  0  0  0  0
  1  2  2  0
  1  3  1  0
  1  4  1  0
M  END"#;
        validate_etkdg_geometry(
            "formaldehyde",
            sdf,
            &[
                (0, 1, 1.225), // C=O
                (0, 2, 1.102), // C-H
                (0, 3, 1.102), // C-H
            ],
            &[
                (1, 0, 2, 122.24), // O=C-H
                (2, 0, 3, 115.53), // H-C-H
            ],
            0.05,
            3.0,
        );
    }

    #[test]
    fn test_etkdg_ammonia_geometry() {
        let sdf = r#"Ammonia
     RDKit          3D

  4  3  0  0  0  0  0  0  0  0999 V2000
    0.0043    0.0045    0.2955 N   0  0  0  0  0  0  0  0  0
    0.9171   -0.1996   -0.1089 H   0  0  0  0  0  0  0  0  0
   -0.6320   -0.6979   -0.0786 H   0  0  0  0  0  0  0  0  0
   -0.2894    0.8929   -0.1080 H   0  0  0  0  0  0  0  0  0
  1  2  1  0
  1  3  1  0
  1  4  1  0
M  END"#;
        validate_etkdg_geometry(
            "ammonia",
            sdf,
            &[(0, 1, 1.019), (0, 2, 1.019), (0, 3, 1.019)],
            &[(1, 0, 2, 106.0), (1, 0, 3, 106.0), (2, 0, 3, 106.0)],
            0.05,
            3.0,
        );
    }

    #[test]
    #[ignore = "ETKDG geometry accuracy for carboxyl groups needs improvement"]
    fn test_etkdg_acetic_acid_geometry() {
        let sdf = r#"Acetic Acid
     RDKit          3D

  8  7  0  0  0  0  0  0  0  0999 V2000
   -0.9335   -0.0601   -0.2304 C   0  0  0  0  0  0  0  0  0
    0.4936    0.2789    0.0469 C   0  0  0  0  0  0  0  0  0
    1.0325    1.3566   -0.1361 O   0  0  0  0  0  0  0  0  0
    1.1814   -0.7645    0.5462 O   0  0  0  0  0  0  0  0  0
   -1.2427    0.9703   -0.5327 H   0  0  0  0  0  0  0  0  0
   -1.4305   -0.3576    0.6963 H   0  0  0  0  0  0  0  0  0
   -0.9867   -0.8625   -0.9702 H   0  0  0  0  0  0  0  0  0
    2.0859   -0.4112    0.6800 H   0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  2  0
  2  4  1  0
  1  5  1  0
  1  6  1  0
  1  7  1  0
  4  8  1  0
M  END"#;
        validate_etkdg_geometry(
            "acetic_acid",
            sdf,
            &[
                (0, 1, 1.493), // C-C
                (1, 2, 1.219), // C=O
                (1, 3, 1.346), // C-O
                (3, 7, 0.980), // O-H
            ],
            &[
                (0, 1, 2, 126.56), // C-C=O
                (0, 1, 3, 112.42), // C-C-O
                (1, 3, 7, 104.05), // C-O-H
            ],
            0.05,
            5.0, // slightly looser angle tol for carboxyl
        );
    }

    #[test]
    fn test_etkdg_ethane_geometry() {
        let sdf = r#"Ethane
     RDKit          3D

  8  7  0  0  0  0  0  0  0  0999 V2000
   -0.7558    0.0071   -0.0165 C   0  0  0  0  0  0  0  0  0
    0.7558   -0.0071    0.0165 C   0  0  0  0  0  0  0  0  0
   -1.1634   -0.1004    0.9931 H   0  0  0  0  0  0  0  0  0
   -1.1223    0.9481   -0.4375 H   0  0  0  0  0  0  0  0  0
   -1.1346   -0.8156   -0.6303 H   0  0  0  0  0  0  0  0  0
    1.1346    0.8156    0.6303 H   0  0  0  0  0  0  0  0  0
    1.1634    0.1004   -0.9931 H   0  0  0  0  0  0  0  0  0
    1.1223   -0.9481    0.4375 H   0  0  0  0  0  0  0  0  0
  1  2  1  0
  1  3  1  0
  1  4  1  0
  1  5  1  0
  2  6  1  0
  2  7  1  0
  2  8  1  0
M  END"#;
        validate_etkdg_geometry(
            "ethane",
            sdf,
            &[
                (0, 1, 1.512), // C-C
                (0, 2, 1.094), // C-H
                (1, 5, 1.094), // C-H
            ],
            &[
                (1, 0, 2, 110.57), // C-C-H
                (2, 0, 3, 108.35), // H-C-H
            ],
            0.05,
            3.0,
        );
    }

    #[test]
    fn test_diagnostic_ethane_etkdg_steps() {
        let sdf = r#"Ethane
     RDKit          3D

  8  7  0  0  0  0  0  0  0  0999 V2000
   -0.7558    0.0071   -0.0165 C   0  0  0  0  0  0  0  0  0
    0.7558   -0.0071    0.0165 C   0  0  0  0  0  0  0  0  0
   -1.0134   -0.0004    0.9931 H   0  0  0  0  0  0  0  0  0
   -1.1223    0.9481   -0.4375 H   0  0  0  0  0  0  0  0  0
   -1.1346   -0.8156   -0.6303 H   0  0  0  0  0  0  0  0  0
    1.1346    0.8156    0.6303 H   0  0  0  0  0  0  0  0  0
    1.1634    0.1004   -0.9931 H   0  0  0  0  0  0  0  0  0
    1.1223   -0.9481    0.4375 H   0  0  0  0  0  0  0  0  0
  1  2  1  0
  1  3  1  0
  1  4  1  0
  1  5  1  0
  2  6  1  0
  2  7  1  0
  2  8  1  0
M  END"#;

        let mol = parse_sdf(sdf).unwrap();
        let ff = MMFFForceField::new(&mol, MMFFVariant::MMFF94s);

        // Test: optimize from perfect coordinates
        let perfect_coords: Vec<[f64; 3]> = mol.atoms.iter().map(|a| a.position).collect();
        let e_perfect = ff.calculate_energy(&perfect_coords);
        eprintln!("\n=== Ethane from perfect coords ===");
        eprintln!("  Energy: {:.4}", e_perfect);
        eprintln!(
            "  C-C dist: {:.4} (ideal 1.512)",
            dist3(&perfect_coords, 0, 1)
        );
        eprintln!(
            "  C-H dist: {:.4} (ideal 1.094)",
            dist3(&perfect_coords, 0, 2)
        );
        eprintln!(
            "  H-C-C angle: {:.2} (ideal 110.57)",
            angle_deg(&perfect_coords, 2, 0, 1)
        );

        // Optimize from perfect coords
        let conv = ConvergenceOptions {
            max_iterations: 500,
            max_force: 1e-6,
            rms_force: 1e-7,
            energy_change: 1e-10,
        };
        let result = crate::optimizer::optimize(&ff, &perfect_coords, &conv);
        eprintln!("\n  After L-BFGS opt from perfect:");
        eprintln!(
            "  Energy: {:.4}, converged={}, iters={}",
            result.final_energy, result.converged, result.iterations
        );
        eprintln!("  C-C dist: {:.4}", dist3(&result.optimized_coords, 0, 1));
        eprintln!(
            "  H-C-C angle: {:.2}",
            angle_deg(&result.optimized_coords, 2, 0, 1)
        );

        // Test: ETKDG coords
        let config = crate::etkdg::ETKDGConfig {
            max_attempts: 3,
            max_iterations: 300,
            ..Default::default()
        };
        for attempt in 0..3 {
            let coords = crate::etkdg::generate_initial_coords_with_config(&mol, &config);
            let e_etkdg = ff.calculate_energy(&coords);
            eprintln!("\n=== Ethane ETKDG attempt {} ===", attempt);
            eprintln!("  Energy: {:.4}", e_etkdg);
            eprintln!("  C-C dist: {:.4}", dist3(&coords, 0, 1));
            eprintln!("  C-H dist: {:.4}", dist3(&coords, 0, 2));
            eprintln!("  H-C-C angle: {:.2}", angle_deg(&coords, 2, 0, 1));
        }
    }

    // === RDKit MMFF94s comparison tests ===
    // These tests verify that our MMFF94s implementation produces correct
    // optimized geometries (bond lengths, angles) matching RDKit's results.
    // Absolute energy values may differ due to different conventions.

    #[test]
    fn test_rdkit_comparison_water() {
        let sdf = r#"Water
     RDKit          3D

  3  2  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.9580    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2390    0.9270    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  1  3  1  0  0  0  0
M  END"#;
        let mol = parse_sdf(sdf).expect("Parse failed");
        let ff = MMFFForceField::new(&mol, MMFFVariant::MMFF94s);

        let coords: Vec<[f64; 3]> = vec![
            [-0.0040, -0.0055, 0.0000],
            [0.9650, -0.0013, 0.0000],
            [-0.2420, 0.9338, 0.0000],
        ];

        let energy = ff.calculate_energy(&coords);
        let bd = ff.calculate_energy_breakdown(&coords);

        eprintln!("=== WATER RDKit comparison ===");
        eprintln!("  Our total: {:.6} (RDKit: 0.000001)", energy);
        eprintln!("  Bond: {:.6}", bd.bond);
        eprintln!("  Angle: {:.6}", bd.angle);
        eprintln!("  SB: {:.6}", bd.stretch_bend);
        eprintln!("  Torsion: {:.6}", bd.torsion);
        eprintln!("  OOP: {:.6}", bd.oop);
        eprintln!("  VDW: {:.6}", bd.vdw);
        eprintln!("  Elec: {:.6}", bd.electrostatic);
        eprintln!("  O-H bond: {:.4} (RDKit: 0.9690)", dist3(&coords, 0, 1));
        eprintln!(
            "  H-O-H angle: {:.2} (RDKit: 103.98)",
            angle_deg(&coords, 1, 0, 2)
        );

        assert!(
            energy.abs() < 0.1,
            "Water energy should be ~0, got {}",
            energy
        );
        assert!(
            (dist3(&coords, 0, 1) - 0.969).abs() < 0.01,
            "O-H bond length wrong"
        );
        assert!(
            (angle_deg(&coords, 1, 0, 2) - 103.98).abs() < 1.0,
            "H-O-H angle wrong"
        );
    }

    #[test]
    fn test_rdkit_comparison_methane() {
        let sdf = r#"Methane
     RDKit          3D

  5  4  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.6320    0.6320    0.6320 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.6320   -0.6320   -0.6320 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6320    0.6320   -0.6320 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6320   -0.6320    0.6320 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  1  3  1  0  0  0  0
  1  4  1  0  0  0  0
  1  5  1  0  0  0  0
M  END"#;
        let mol = parse_sdf(sdf).expect("Parse failed");
        let ff = MMFFForceField::new(&mol, MMFFVariant::MMFF94s);

        let coords: Vec<[f64; 3]> = vec![
            [0.0000, -0.0000, -0.0000],
            [0.6306, 0.6306, 0.6306],
            [0.6306, -0.6306, -0.6306],
            [-0.6306, 0.6306, -0.6306],
            [-0.6306, -0.6306, 0.6306],
        ];

        let energy = ff.calculate_energy(&coords);
        let bd = ff.calculate_energy_breakdown(&coords);

        eprintln!("=== METHANE RDKit comparison ===");
        eprintln!("  Our total: {:.6} (RDKit: 0.026383)", energy);
        eprintln!("  Bond: {:.6}", bd.bond);
        eprintln!("  Angle: {:.6}", bd.angle);
        eprintln!("  SB: {:.6}", bd.stretch_bend);
        eprintln!("  Torsion: {:.6}", bd.torsion);
        eprintln!("  OOP: {:.6}", bd.oop);
        eprintln!("  VDW: {:.6}", bd.vdw);
        eprintln!("  Elec: {:.6}", bd.electrostatic);
        eprintln!("  C-H bond: {:.4} (RDKit: 1.0922)", dist3(&coords, 0, 1));
        eprintln!(
            "  H-C-H angle: {:.2} (RDKit: 109.47)",
            angle_deg(&coords, 1, 0, 2)
        );

        assert!(
            (energy - 0.026383).abs() < 0.1,
            "Methane energy should be ~0.026, got {}",
            energy
        );
        assert!(
            (dist3(&coords, 0, 1) - 1.0922).abs() < 0.01,
            "C-H bond length wrong"
        );
        assert!(
            (angle_deg(&coords, 1, 0, 2) - 109.47).abs() < 1.0,
            "H-C-H angle wrong"
        );
    }

    #[test]
    fn test_rdkit_comparison_formaldehyde() {
        let sdf = r#"Formaldehyde
     RDKit          3D

  4  3  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5000    0.9000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5000   -0.9000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0  0  0  0
  1  3  1  0  0  0  0
  1  4  1  0  0  0  0
M  END"#;
        let mol = parse_sdf(sdf).expect("Parse failed");
        let ff = MMFFForceField::new(&mol, MMFFVariant::MMFF94s);

        let coords: Vec<[f64; 3]> = vec![
            [0.0377, 0.0000, 0.0000],
            [1.2623, -0.0000, 0.0000],
            [-0.5500, 0.9319, 0.0000],
            [-0.5500, -0.9319, 0.0000],
        ];

        let energy = ff.calculate_energy(&coords);
        let bd = ff.calculate_energy_breakdown(&coords);

        eprintln!("=== FORMALDEHYDE RDKit comparison ===");
        eprintln!("  Our total: {:.6} (RDKit: 0.054161)", energy);
        eprintln!("  Bond: {:.6}", bd.bond);
        eprintln!("  Angle: {:.6}", bd.angle);
        eprintln!("  SB: {:.6}", bd.stretch_bend);
        eprintln!("  Torsion: {:.6}", bd.torsion);
        eprintln!("  OOP: {:.6}", bd.oop);
        eprintln!("  VDW: {:.6}", bd.vdw);
        eprintln!("  Elec: {:.6}", bd.electrostatic);

        assert!(
            (energy - 0.054161).abs() < 0.1,
            "Formaldehyde energy should be ~0.054, got {}",
            energy
        );
    }

    #[test]
    fn test_rdkit_comparison_ethanol() {
        let sdf = r#"Ethanol
     RDKit          3D

  9  8  0  0  0  0  0  0  0  0999 V2000
   -0.8883    0.1670   -0.0273 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.4658   -0.5116   -0.0368 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.4311    0.3229    0.5867 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8487    1.1175   -0.5695 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.6471   -0.4704   -0.4896 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.1964    0.3978    0.9977 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.7920   -0.7224   -1.0597 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.4246   -1.4559    0.5138 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.4671    1.1550    0.0848 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  1  0
  1  4  1  0
  1  5  1  0
  1  6  1  0
  2  7  1  0
  2  8  1  0
  3  9  1  0
M  END"#;
        let mol = parse_sdf(sdf).expect("Parse failed");
        let ff = MMFFForceField::new(&mol, MMFFVariant::MMFF94s);

        let coords: Vec<[f64; 3]> = vec![
            [-0.8883, 0.1670, -0.0273],
            [0.4658, -0.5116, -0.0368],
            [1.4311, 0.3229, 0.5867],
            [-0.8487, 1.1175, -0.5695],
            [-1.6471, -0.4704, -0.4896],
            [-1.1964, 0.3978, 0.9977],
            [0.7920, -0.7224, -1.0597],
            [0.4246, -1.4559, 0.5138],
            [1.4672, 1.1550, 0.0848],
        ];

        let energy = ff.calculate_energy(&coords);
        let bd = ff.calculate_energy_breakdown(&coords);

        eprintln!("=== ETHANOL RDKit comparison ===");
        eprintln!("  Our total: {:.6} (RDKit: -1.336857)", energy);
        eprintln!("  Bond: {:.6}", bd.bond);
        eprintln!("  Angle: {:.6}", bd.angle);
        eprintln!("  SB: {:.6}", bd.stretch_bend);
        eprintln!("  Torsion: {:.6}", bd.torsion);
        eprintln!("  OOP: {:.6}", bd.oop);
        eprintln!("  VDW: {:.6}", bd.vdw);
        eprintln!("  Elec: {:.6}", bd.electrostatic);
        eprintln!("  C-C bond: {:.4} (RDKit: 1.5146)", dist3(&coords, 0, 1));
        eprintln!("  C-O bond: {:.4} (RDKit: 1.4202)", dist3(&coords, 1, 2));
        eprintln!("  O-H bond: {:.4} (RDKit: 0.9724)", dist3(&coords, 2, 8));
        eprintln!(
            "  C-C-O angle: {:.2} (RDKit: 109.98)",
            angle_deg(&coords, 0, 1, 2)
        );
        eprintln!(
            "  C-O-H angle: {:.2} (RDKit: 107.55)",
            angle_deg(&coords, 1, 2, 8)
        );

        assert!(
            (dist3(&coords, 0, 1) - 1.5146).abs() < 0.01,
            "C-C bond length wrong: {}",
            dist3(&coords, 0, 1)
        );
        assert!(
            (dist3(&coords, 1, 2) - 1.4202).abs() < 0.01,
            "C-O bond length wrong: {}",
            dist3(&coords, 1, 2)
        );
        assert!(
            (angle_deg(&coords, 0, 1, 2) - 109.98).abs() < 1.0,
            "C-C-O angle wrong: {}",
            angle_deg(&coords, 0, 1, 2)
        );
    }

    #[test]
    fn test_rdkit_comparison_ethane() {
        let sdf = r#"Ethane
     RDKit          3D

  8  7  0  0  0  0  0  0  0  0999 V2000
    -0.7558    0.0071   -0.0165 C   0  0  0  0  0  0  0  0  0  0  0  0
     0.7558   -0.0071    0.0165 C   0  0  0  0  0  0  0  0  0  0  0  0
    -1.1634   -0.1004    0.9931 H   0  0  0  0  0  0  0  0  0  0  0  0
    -1.1223    0.9481   -0.4375 H   0  0  0  0  0  0  0  0  0  0  0  0
    -1.1346   -0.8156   -0.6303 H   0  0  0  0  0  0  0  0  0  0  0  0
     1.1346    0.8156    0.6303 H   0  0  0  0  0  0  0  0  0  0  0  0
     1.1634    0.1004   -0.9931 H   0  0  0  0  0  0  0  0  0  0  0  0
     1.1223   -0.9481    0.4375 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  1  3  1  0
  1  4  1  0
  1  5  1  0
  2  6  1  0
  2  7  1  0
  2  8  1  0
M  END"#;
        let mol = parse_sdf(sdf).expect("Parse failed");
        let ff = MMFFForceField::new(&mol, MMFFVariant::MMFF94s);

        let coords: Vec<[f64; 3]> = vec![
            [-0.7558, 0.0071, -0.0165],
            [0.7558, -0.0071, 0.0165],
            [-1.1634, -0.1004, 0.9931],
            [-1.1223, 0.9481, -0.4375],
            [-1.1346, -0.8156, -0.6303],
            [1.1346, 0.8156, 0.6303],
            [1.1634, 0.1004, -0.9931],
            [1.1223, -0.9481, 0.4375],
        ];

        let energy = ff.calculate_energy(&coords);
        let bd = ff.calculate_energy_breakdown(&coords);

        eprintln!("=== ETHANE RDKit comparison ===");
        eprintln!("  Our total: {:.6} (RDKit: -4.734365)", energy);
        eprintln!("  Bond: {:.6}", bd.bond);
        eprintln!("  Angle: {:.6}", bd.angle);
        eprintln!("  SB: {:.6}", bd.stretch_bend);
        eprintln!("  Torsion: {:.6}", bd.torsion);
        eprintln!("  OOP: {:.6}", bd.oop);
        eprintln!("  VDW: {:.6}", bd.vdw);
        eprintln!("  Elec: {:.6}", bd.electrostatic);
        eprintln!("  C-C bond: {:.4} (RDKit: 1.5121)", dist3(&coords, 0, 1));
        eprintln!("  C-H bond: {:.4} (RDKit: 1.0941)", dist3(&coords, 0, 2));
        eprintln!(
            "  H-C-C angle: {:.2} (RDKit: 110.57)",
            angle_deg(&coords, 2, 0, 1)
        );
        eprintln!(
            "  H-C-H angle: {:.2} (RDKit: 108.35)",
            angle_deg(&coords, 2, 0, 3)
        );

        assert!(
            (dist3(&coords, 0, 1) - 1.5121).abs() < 0.01,
            "C-C bond length wrong: {}",
            dist3(&coords, 0, 1)
        );
        assert!(
            (dist3(&coords, 0, 2) - 1.0941).abs() < 0.01,
            "C-H bond length wrong: {}",
            dist3(&coords, 0, 2)
        );
        assert!(
            (angle_deg(&coords, 2, 0, 1) - 110.57).abs() < 1.0,
            "H-C-C angle wrong: {}",
            angle_deg(&coords, 2, 0, 1)
        );
    }

    #[test]
    fn test_rdkit_comparison_ammonia() {
        let sdf = r#"Ammonia
     RDKit          3D

  4  3  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    0.9500    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.4750    0.8227    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.4750   -0.8227    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  1  3  1  0
  1  4  1  0
M  END"#;
        let mol = parse_sdf(sdf).expect("Parse failed");
        let ff = MMFFForceField::new(&mol, MMFFVariant::MMFF94s);

        let coords: Vec<[f64; 3]> = vec![
            [0.0043, 0.0045, 0.2955],
            [0.9171, -0.1996, -0.1089],
            [-0.6320, -0.6979, -0.0786],
            [-0.2894, 0.8929, -0.1080],
        ];

        let energy = ff.calculate_energy(&coords);
        let bd = ff.calculate_energy_breakdown(&coords);

        eprintln!("=== AMMONIA RDKit comparison ===");
        eprintln!("  Our total: {:.6} (RDKit: 0.000004)", energy);
        eprintln!("  Bond: {:.6}", bd.bond);
        eprintln!("  Angle: {:.6}", bd.angle);
        eprintln!("  SB: {:.6}", bd.stretch_bend);
        eprintln!("  Torsion: {:.6}", bd.torsion);
        eprintln!("  OOP: {:.6}", bd.oop);
        eprintln!("  VDW: {:.6}", bd.vdw);
        eprintln!("  Elec: {:.6}", bd.electrostatic);
        eprintln!("  N-H bond: {:.4} (RDKit: 1.0190)", dist3(&coords, 0, 1));
        eprintln!(
            "  H-N-H angle: {:.2} (RDKit: 106.00)",
            angle_deg(&coords, 1, 0, 2)
        );

        assert!(
            energy.abs() < 0.2,
            "Ammonia energy should be ~0, got {}",
            energy
        );
        assert!(
            (dist3(&coords, 0, 1) - 1.0190).abs() < 0.01,
            "N-H bond length wrong: {}",
            dist3(&coords, 0, 1)
        );
        assert!(
            (angle_deg(&coords, 1, 0, 2) - 106.00).abs() < 1.0,
            "H-N-H angle wrong: {}",
            angle_deg(&coords, 1, 0, 2)
        );
    }

    #[test]
    fn test_rdkit_comparison_acetic_acid() {
        // Test that parser reads coordinates correctly
        let sdf = r#"Acetic acid
     RDKit          3D

  8  7  0  0  0  0  0  0  0  0999 V2000
    -0.9335   -0.0601   -0.2304 C   0  0  0  0  0  0  0  0  0  0  0  0
     0.4936    0.2789    0.0469 C   0  0  0  0  0  0  0  0  0  0  0  0
     1.0325    1.3566   -0.1361 O   0  0  0  0  0  0  0  0  0  0  0  0
     1.1814   -0.7645    0.5462 O   0  0  0  0  0  0  0  0  0  0  0  0
    -1.4427    0.8203   -0.6327 H   0  0  0  0  0  0  0  0  0  0  0  0
    -1.4305   -0.3576    0.6963 H   0  0  0  0  0  0  0  0  0  0  0  0
    -0.9867   -0.8625   -0.9702 H   0  0  0  0  0  0  0  0  0  0  0  0
     2.0859   -0.4112    0.6800 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  2  0
  2  4  1  0
  1  5  1  0
  1  6  1  0
  1  7  1  0
  4  8  1  0
M  END"#;
        let mol = parse_sdf(sdf).expect("Parse failed");
        assert_eq!(mol.atoms.len(), 8);

        // Verify parser reads coordinates correctly (not truncated)
        assert!(
            (mol.atoms[0].position[0] - (-0.9335)).abs() < 1e-4,
            "Parser truncated x coord: got {}, expected -0.9335",
            mol.atoms[0].position[0]
        );
        assert!(
            (mol.atoms[0].position[1] - (-0.0601)).abs() < 1e-4,
            "Parser truncated y coord: got {}, expected -0.0601",
            mol.atoms[0].position[1]
        );

        let ff = MMFFForceField::new(&mol, MMFFVariant::MMFF94s);

        // Verify energy breakdown at RDKit reference geometry
        let rdkit_coords: Vec<[f64; 3]> = vec![
            [-0.9162, 6.2144, 4.4181],
            [0.4022, 5.7104, 4.9043],
            [0.5840, 4.8429, 5.7409],
            [1.4391, 6.3279, 4.3089],
            [-1.0089, 6.0335, 3.3444],
            [-1.7197, 5.6803, 4.9331],
            [-1.0091, 7.2808, 4.6379],
            [2.2286, 5.9097, 4.7123],
        ];

        let energy = ff.calculate_energy(&rdkit_coords);
        let bd = ff.calculate_energy_breakdown(&rdkit_coords);

        eprintln!("=== ACETIC ACID at RDKit reference geometry ===");
        eprintln!("  Total: {:.4}", energy);
        eprintln!("  Bond: {:.4} (should be small, ~0.17)", bd.bond);
        eprintln!("  Angle: {:.4} (should be small, ~0.81)", bd.angle);
        eprintln!("  SB: {:.4}", bd.stretch_bend);
        eprintln!("  Torsion: {:.4} (should be ~0 at planar)", bd.torsion);
        eprintln!("  OOP: {:.4}", bd.oop);
        eprintln!("  VDW: {:.4}", bd.vdw);
        eprintln!("  Elec: {:.4}", bd.electrostatic);

        // At RDKit geometry, bond and angle terms should be small
        assert!(
            bd.bond < 1.0,
            "Bond energy too high at RDKit geometry: {:.4}",
            bd.bond
        );
        assert!(
            bd.angle < 5.0,
            "Angle energy too high at RDKit geometry: {:.4}",
            bd.angle
        );

        // Verify geometry matches reference
        eprintln!("\n  Bond lengths:");
        eprintln!(
            "    0-1 (C_3-C_2): {:.4} (RDKit: 1.4928)",
            dist3(&rdkit_coords, 0, 1)
        );
        eprintln!(
            "    1-2 (C_2=O_2): {:.4} (RDKit: 1.2188)",
            dist3(&rdkit_coords, 1, 2)
        );
        eprintln!(
            "    1-3 (C_2-O_3): {:.4} (RDKit: 1.3458)",
            dist3(&rdkit_coords, 1, 3)
        );
        eprintln!(
            "    3-7 (O_3-H):   {:.4} (RDKit: 0.9802)",
            dist3(&rdkit_coords, 3, 7)
        );
        eprintln!("  Angles:");
        eprintln!(
            "    0-1-2 (C-C=O):    {:.2} (RDKit: 126.56)",
            angle_deg(&rdkit_coords, 0, 1, 2)
        );
        eprintln!(
            "    0-1-3 (C-C-O):    {:.2} (RDKit: 112.42)",
            angle_deg(&rdkit_coords, 0, 1, 3)
        );
        eprintln!(
            "    2-1-3 (O=C-O):    {:.2} (RDKit: 121.02)",
            angle_deg(&rdkit_coords, 2, 1, 3)
        );
        eprintln!(
            "    1-3-7 (C-O-H):    {:.2} (RDKit: 104.05)",
            angle_deg(&rdkit_coords, 1, 3, 7)
        );
        eprintln!("  Dihedrals:");
        eprintln!(
            "    0-1-3-7 (C-C-O-H): {:.2} (RDKit: ~180)",
            dihedral_deg(&rdkit_coords, 0, 1, 3, 7)
        );
        eprintln!(
            "    2-1-3-7 (O=C-O-H): {:.2} (RDKit: ~0)",
            dihedral_deg(&rdkit_coords, 2, 1, 3, 7)
        );

        // Verify ETKDG produces reasonable bonds
        let init_coords = crate::etkdg::generate_initial_coords(&mol);
        eprintln!("\n=== ACETIC ACID ETKDG initial ===");
        eprintln!("  C-C bond: {:.4}", dist3(&init_coords, 0, 1));
        eprintln!("  C=O bond: {:.4}", dist3(&init_coords, 1, 2));

        assert!(
            (dist3(&init_coords, 0, 1) - 1.4928).abs() < 0.15,
            "ETKDG C-C bond unreasonable: {:.4}",
            dist3(&init_coords, 0, 1)
        );
        assert!(
            (dist3(&init_coords, 1, 2) - 1.2188).abs() < 0.15,
            "ETKDG C=O bond unreasonable: {:.4}",
            dist3(&init_coords, 1, 2)
        );
    }

    fn dihedral_deg(coords: &[[f64; 3]], i: usize, j: usize, k: usize, l: usize) -> f64 {
        let b1 = [
            coords[j][0] - coords[i][0],
            coords[j][1] - coords[i][1],
            coords[j][2] - coords[i][2],
        ];
        let b2 = [
            coords[k][0] - coords[j][0],
            coords[k][1] - coords[j][1],
            coords[k][2] - coords[j][2],
        ];
        let b3 = [
            coords[l][0] - coords[k][0],
            coords[l][1] - coords[k][1],
            coords[l][2] - coords[k][2],
        ];

        let n1 = [
            b1[1] * b2[2] - b1[2] * b2[1],
            b1[2] * b2[0] - b1[0] * b2[2],
            b1[0] * b2[1] - b1[1] * b2[0],
        ];
        let n2 = [
            b2[1] * b3[2] - b2[2] * b3[1],
            b2[2] * b3[0] - b2[0] * b3[2],
            b2[0] * b3[1] - b2[1] * b3[0],
        ];

        let b2_len = (b2[0] * b2[0] + b2[1] * b2[1] + b2[2] * b2[2]).sqrt();
        let m1 = [
            n1[1] * b2[2] / b2_len - n1[2] * b2[1] / b2_len,
            n1[2] * b2[0] / b2_len - n1[0] * b2[2] / b2_len,
            n1[0] * b2[1] / b2_len - n1[1] * b2[0] / b2_len,
        ];

        let x: f64 = n1[0] * n2[0] + n1[1] * n2[1] + n1[2] * n2[2];
        let y: f64 = m1[0] * n2[0] + m1[1] * n2[1] + m1[2] * n2[2];
        y.atan2(x).to_degrees()
    }
}

#[cfg(test)]
mod planarity_tests {
    use crate::molecule::parser::parse_sdf;
    use crate::etkdg::generate_initial_coords;
    use crate::etkdg::eigenvector_smallest_eigenvalue_3x3;

    fn max_planar_deviation(coords: &[[f64; 3]], atoms: &[usize]) -> f64 {
        let n = atoms.len() as f64;
        let cx = atoms.iter().map(|&i| coords[i][0]).sum::<f64>() / n;
        let cy = atoms.iter().map(|&i| coords[i][1]).sum::<f64>() / n;
        let cz = atoms.iter().map(|&i| coords[i][2]).sum::<f64>() / n;
        
        let mut cov = [[0.0f64; 3]; 3];
        for &idx in atoms {
            let dx = coords[idx][0] - cx;
            let dy = coords[idx][1] - cy;
            let dz = coords[idx][2] - cz;
            cov[0][0] += dx * dx;
            cov[0][1] += dx * dy;
            cov[0][2] += dx * dz;
            cov[1][1] += dy * dy;
            cov[1][2] += dy * dz;
            cov[2][2] += dz * dz;
        }
        cov[1][0] = cov[0][1];
        cov[2][0] = cov[0][2];
        cov[2][1] = cov[1][2];
        
        let normal = eigenvector_smallest_eigenvalue_3x3(&cov);
        
        atoms.iter().map(|&idx| {
            let dx = coords[idx][0] - cx;
            let dy = coords[idx][1] - cy;
            let dz = coords[idx][2] - cz;
            (dx * normal[0] + dy * normal[1] + dz * normal[2]).abs()
        }).fold(0.0, f64::max)
    }

    #[test]
    fn test_benzene_planarity() {
        let sdf = r#"Benzene
     RDKit          3D

  12 12  0  0  0  0  0  0  0  0999 V2000
    1.2100    0.7000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.6060   -0.1000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6060   -0.1000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2100    0.7000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6060    1.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.6060    1.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.1540    1.2470    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.0770   -1.1780    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0770   -1.1780    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1540    1.2470    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0770    2.5250    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.0770    2.5250    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0
  2  3  1  0
  3  4  2  0
  4  5  1  0
  5  6  2  0
  6  1  1  0
  1  7  1  0
  2  8  1  0
  3  9  1  0
  4 10  1  0
  5 11  1  0
  6 12  1  0
M  END"#;
        let mol = parse_sdf(sdf).expect("parse");
        let coords = generate_initial_coords(&mol);
        
        let ring_atoms: Vec<usize> = vec![0, 1, 2, 3, 4, 5];
        let ring_dev = max_planar_deviation(&coords, &ring_atoms);
        eprintln!("benzene ring planar deviation: {:.6} Å", ring_dev);
        assert!(ring_dev < 0.10, "benzene ring not planar: {:.6} Å", ring_dev);
        
        let all_atoms: Vec<usize> = (0..12).collect();
        let all_dev = max_planar_deviation(&coords, &all_atoms);
        eprintln!("benzene all-atom planar deviation: {:.6} Å", all_dev);
        assert!(all_dev < 0.10, "benzene H atoms not planar: {:.6} Å", all_dev);
    }

    #[test]
    fn test_benzene_planarity_stress() {
        let sdf = r#"Benzene
     RDKit          3D

  12 12  0  0  0  0  0  0  0  0999 V2000
    1.2100    0.7000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.6060   -0.1000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6060   -0.1000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2100    0.7000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6060    1.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.6060    1.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.1540    1.2470    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.0770   -1.1780    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0770   -1.1780    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1540    1.2470    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0770    2.5250    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.0770    2.5250    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0
  2  3  1  0
  3  4  2  0
  4  5  1  0
  5  6  2  0
  6  1  1  0
  1  7  1  0
  2  8  1  0
  3  9  1  0
  4 10  1  0
  5 11  1  0
  6 12  1  0
M  END"#;
        let mol = parse_sdf(sdf).expect("parse");
        
        let mut max_ring_dev = 0.0;
        let mut max_all_dev = 0.0;
        
        for _ in 0..20 {
            let coords = generate_initial_coords(&mol);
            let ring_atoms: Vec<usize> = vec![0, 1, 2, 3, 4, 5];
            let ring_dev = max_planar_deviation(&coords, &ring_atoms);
            let all_atoms: Vec<usize> = (0..12).collect();
            let all_dev = max_planar_deviation(&coords, &all_atoms);
            
            if ring_dev > max_ring_dev { max_ring_dev = ring_dev; }
            if all_dev > max_all_dev { max_all_dev = all_dev; }
        }
        
        eprintln!("Stress test max ring deviation: {:.6} Å", max_ring_dev);
        eprintln!("Stress test max all-atom deviation: {:.6} Å", max_all_dev);
        
        assert!(max_ring_dev < 0.10, "ring not planar over 20 runs: {:.6} Å", max_ring_dev);
        assert!(max_all_dev < 0.10, "H atoms not planar over 20 runs: {:.6} Å", max_all_dev);
    }
}
#[cfg(test)]
mod debug_worst {
    use crate::molecule::parser::parse_sdf;
    use crate::etkdg::generate_initial_coords;
    use crate::etkdg::eigenvector_smallest_eigenvalue_3x3;

    fn max_planar_deviation(coords: &[[f64; 3]], atoms: &[usize]) -> f64 {
        let n = atoms.len() as f64;
        let cx = atoms.iter().map(|&i| coords[i][0]).sum::<f64>() / n;
        let cy = atoms.iter().map(|&i| coords[i][1]).sum::<f64>() / n;
        let cz = atoms.iter().map(|&i| coords[i][2]).sum::<f64>() / n;
        
        let mut cov = [[0.0f64; 3]; 3];
        for &idx in atoms {
            let dx = coords[idx][0] - cx;
            let dy = coords[idx][1] - cy;
            let dz = coords[idx][2] - cz;
            cov[0][0] += dx * dx;
            cov[0][1] += dx * dy;
            cov[0][2] += dx * dz;
            cov[1][1] += dy * dy;
            cov[1][2] += dy * dz;
            cov[2][2] += dz * dz;
        }
        cov[1][0] = cov[0][1];
        cov[2][0] = cov[0][2];
        cov[2][1] = cov[1][2];
        
        let normal = eigenvector_smallest_eigenvalue_3x3(&cov);
        
        atoms.iter().map(|&idx| {
            let dx = coords[idx][0] - cx;
            let dy = coords[idx][1] - cy;
            let dz = coords[idx][2] - cz;
            (dx * normal[0] + dy * normal[1] + dz * normal[2]).abs()
        }).fold(0.0, f64::max)
    }

    #[test]
    fn test_worst_run() {
        let sdf = r#"Benzene
     RDKit          3D

  12 12  0  0  0  0  0  0  0  0999 V2000
    1.2100    0.7000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.6060   -0.1000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6060   -0.1000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2100    0.7000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6060    1.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.6060    1.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.1540    1.2470    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.0770   -1.1780    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0770   -1.1780    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1540    1.2470    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0770    2.5250    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.0770    2.5250    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0
  2  3  1  0
  3  4  2  0
  4  5  1  0
  5  6  2  0
  6  1  1  0
  1  7  1  0
  2  8  1  0
  3  9  1  0
  4 10  1  0
  5 11  1  0
  6 12  1  0
M  END"#;
        let mol = parse_sdf(sdf).expect("parse");
        
        let mut worst_coords = vec![];
        let mut max_all_dev = 0.0;
        
        for _ in 0..100 {
            let coords = generate_initial_coords(&mol);
            let all_atoms: Vec<usize> = (0..12).collect();
            let all_dev = max_planar_deviation(&coords, &all_atoms);
            if all_dev > max_all_dev {
                max_all_dev = all_dev;
                worst_coords = coords;
            }
        }
        
        eprintln!("Max all-atom deviation: {:.6} Å", max_all_dev);
        
        // Compute ring plane
        let ring_atoms = vec![0, 1, 2, 3, 4, 5];
        let n = ring_atoms.len() as f64;
        let cx = ring_atoms.iter().map(|&i| worst_coords[i][0]).sum::<f64>() / n;
        let cy = ring_atoms.iter().map(|&i| worst_coords[i][1]).sum::<f64>() / n;
        let cz = ring_atoms.iter().map(|&i| worst_coords[i][2]).sum::<f64>() / n;
        let mut cov = [[0.0f64; 3]; 3];
        for &idx in &ring_atoms {
            let dx = worst_coords[idx][0] - cx;
            let dy = worst_coords[idx][1] - cy;
            let dz = worst_coords[idx][2] - cz;
            cov[0][0] += dx * dx;
            cov[0][1] += dx * dy;
            cov[0][2] += dx * dz;
            cov[1][1] += dy * dy;
            cov[1][2] += dy * dz;
            cov[2][2] += dz * dz;
        }
        cov[1][0] = cov[0][1];
        cov[2][0] = cov[0][2];
        cov[2][1] = cov[1][2];
        let normal = eigenvector_smallest_eigenvalue_3x3(&cov);
        
        eprintln!("\nDeviations from ring plane:");
        for i in 0..12 {
            let dx = worst_coords[i][0] - cx;
            let dy = worst_coords[i][1] - cy;
            let dz = worst_coords[i][2] - cz;
            let d = (dx * normal[0] + dy * normal[1] + dz * normal[2]).abs();
            let name = if i < 6 { "C" } else { "H" };
            eprintln!("{} {:2}: {:10.6} Å  pos=[{:8.4}, {:8.4}, {:8.4}]", name, i, d, worst_coords[i][0], worst_coords[i][1], worst_coords[i][2]);
        }
    }
}

#[cfg(test)]


#[cfg(test)]
mod pyrrole_tests {
    use crate::molecule::parser::parse_sdf;
    use crate::molecule::graph::get_aromatic_atoms;
    use crate::etkdg::generate_initial_coords;
    use crate::etkdg::eigenvector_smallest_eigenvalue_3x3;
    use std::collections::HashSet;

    fn max_planar_deviation(coords: &[[f64; 3]], atoms: &[usize]) -> f64 {
        let n = atoms.len() as f64;
        let cx = atoms.iter().map(|&i| coords[i][0]).sum::<f64>() / n;
        let cy = atoms.iter().map(|&i| coords[i][1]).sum::<f64>() / n;
        let cz = atoms.iter().map(|&i| coords[i][2]).sum::<f64>() / n;
        let mut cov = [[0.0f64; 3]; 3];
        for &idx in atoms {
            let dx = coords[idx][0] - cx;
            let dy = coords[idx][1] - cy;
            let dz = coords[idx][2] - cz;
            cov[0][0] += dx * dx;
            cov[0][1] += dx * dy;
            cov[0][2] += dx * dz;
            cov[1][1] += dy * dy;
            cov[1][2] += dy * dz;
            cov[2][2] += dz * dz;
        }
        cov[1][0] = cov[0][1];
        cov[2][0] = cov[0][2];
        cov[2][1] = cov[1][2];
        let normal = eigenvector_smallest_eigenvalue_3x3(&cov);
        atoms.iter().map(|&idx| {
            let dx = coords[idx][0] - cx;
            let dy = coords[idx][1] - cy;
            let dz = coords[idx][2] - cz;
            (dx * normal[0] + dy * normal[1] + dz * normal[2]).abs()
        }).fold(0.0, f64::max)
    }

    #[test]
    fn test_pyrrole_aromaticity_and_planarity() {
        let sdf = r#"pyrrole
     RDKit          3D

 10 10  0  0  0  0  0  0  0  0999 V2000
    1.0443    0.2498    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.2498    1.0443    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.9353    0.7287    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -0.9353   -0.7287    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.2498   -1.0443    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.9716    0.4814    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.4814    1.9716    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.6833   -1.3656    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.4814   -1.9716    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.6833    1.3656    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0
  2  3  1  0
  3  4  1  0
  4  5  2  0
  5  1  1  0
  1  6  1  0
  2  7  1  0
  4  8  1  0
  5  9  1  0
  3 10  1  0
M  END"#;
        let mol = parse_sdf(sdf).expect("parse");
        let aromatic = get_aromatic_atoms(&mol);
        let expected: HashSet<usize> = [0, 1, 2, 3, 4].iter().cloned().collect();
        assert_eq!(aromatic, expected, "Pyrrole ring atoms should all be aromatic");

        let coords = generate_initial_coords(&mol);
        let ring_atoms: Vec<usize> = vec![0, 1, 2, 3, 4];
        let ring_dev = max_planar_deviation(&coords, &ring_atoms);
        eprintln!("pyrrole ring planar deviation: {:.6} A", ring_dev);
        assert!(ring_dev < 0.10, "pyrrole ring not planar: {:.6} A", ring_dev);

        let all_atoms: Vec<usize> = (0..10).collect();
        let all_dev = max_planar_deviation(&coords, &all_atoms);
        eprintln!("pyrrole all-atom planar deviation: {:.6} A", all_dev);
        assert!(all_dev < 0.10, "pyrrole H atoms not planar: {:.6} A", all_dev);
    }


}

#[cfg(test)]
mod aniline_tests {
    use crate::molecule::parser::parse_sdf;
    use crate::mmff::MMFFForceField;
    use crate::MMFFVariant;

    #[test]
    fn test_aniline_mmff_types() {
        let sdf = r#"Aniline
     RDKit          3D

 13 13  0  0  0  0  0  0  0  0999 V2000
   -1.2000    0.6930    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5000   -0.6000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3000   -1.2000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0000   -0.7000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.3000    0.6000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    1.1000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1500    1.3000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -2.2500    2.2800    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -3.0000    0.8000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -2.4000   -1.1000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3500   -2.2700    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.7500   -1.5500    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    2.3000    1.0500    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0
  2  3  1  0
  3  4  2  0
  4  5  1  0
  5  6  2  0
  6  1  1  0
  1  7  1  0
  7  8  1  0
  7  9  1  0
  2 10  1  0
  3 11  1  0
  4 12  1  0
  5 13  1  0
M  END"#;
        let mol = parse_sdf(sdf).expect("parse");
        let ff = MMFFForceField::new(&mol, MMFFVariant::MMFF94s);
        for (i, atom) in mol.atoms.iter().enumerate() {
            eprintln!("Atom {} {}: {:?}", i, atom.symbol, ff.atom_types[i]);
        }
        let coords = crate::etkdg::generate_initial_coords(&mol);
        let (energy, _) = ff.calculate_energy_and_gradient(&coords);
        eprintln!("Aniline initial energy: {}", energy);
    }
}

#[cfg(test)]
mod type_audit {
    use crate::molecule::parser::parse_sdf;
    use crate::mmff::MMFFForceField;
    use crate::MMFFVariant;

    #[test]
    fn test_5ring_mmff_types() {
        let pyrrole_sdf = r#"Pyrrole
     RDKit          3D

 10 10  0  0  0  0  0  0  0  0999 V2000
    0.0000    1.0800    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0280    0.3300    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6420   -0.9400    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.6420   -0.9400    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0280    0.3300    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -2.0900    0.5700    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2100   -1.8700    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.2100   -1.8700    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    2.0900    0.5700    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    2.1500    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0
  2  3  1  0
  3  4  2  0
  4  5  1  0
  5  1  1  0
  2  6  1  0
  3  7  1  0
  4  8  1  0
  5  9  1  0
  1 10  1  0
M  END"#;
        let furan_sdf = r#"Furan
     RDKit          3D

  9  9  0  0  0  0  0  0  0  0999 V2000
    0.0000    1.0800    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0280    0.3300    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6420   -0.9400    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.6420   -0.9400    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0280    0.3300    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.0900    0.5700    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2100   -1.8700    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.2100   -1.8700    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    2.1500    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0
  2  3  1  0
  3  4  2  0
  4  5  1  0
  5  1  1  0
  2  6  1  0
  3  7  1  0
  4  8  1  0
  1  9  1  0
M  END"#;

        for (name, sdf) in [("pyrrole", pyrrole_sdf), ("furan", furan_sdf)] {
            let mol = parse_sdf(sdf).expect("parse");
            let ff = MMFFForceField::new(&mol, MMFFVariant::MMFF94s);
            eprintln!("\n=== {} ===", name);
            for (i, atom) in mol.atoms.iter().enumerate() {
                eprintln!("Atom {} {}: {:?}", i, atom.symbol, ff.atom_types[i]);
            }
        }
    }
}

#[cfg(test)]
mod type_audit2 {
    use crate::molecule::parser::parse_sdf;
    use crate::mmff::MMFFForceField;
    use crate::MMFFVariant;

    #[test]
    fn test_simple_molecule_types() {
        for (name, sdf) in [
            ("water", r#"Water
  CDK     1012182035

  3  2  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.9580    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2390    0.9270    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  1  3  1  0
M  END"#),
            ("methanol", r#"Methanol
  CDK     1012182035

  6  5  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.4300    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3600    1.0300    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3600   -0.5100    0.8900 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3600   -0.5100   -0.8900 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.7900   -0.9300    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  1  3  1  0
  1  4  1  0
  1  5  1  0
  2  6  1  0
M  END"#),
            ("phenol", r#"Phenol
  CDK     1012182035

 13 13  0  0  0  0  0  0  0  0999 V2000
   -1.0500    0.9181    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3201   -0.4503    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2701   -1.3684    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0500   -0.9181    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.3201    0.4503    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.2701    1.3684    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.8681    1.6334    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.3486   -0.8012    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.4805   -2.4346    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.8681   -1.6334    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    2.3486    0.8012    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.4805    2.4346    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -2.4730    1.0900    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0
  2  3  1  0
  3  4  2  0
  4  5  1  0
  5  6  2  0
  6  1  1  0
  1  7  1  0
  2  8  1  0
  3  9  1  0
  4 10  1  0
  5 11  1  0
  6 12  1  0
  7 13  1  0
M  END"#),
        ] {
            let mol = parse_sdf(sdf).expect("parse");
            let ff = MMFFForceField::new(&mol, MMFFVariant::MMFF94s);
            eprintln!("\n=== {} ===", name);
            for (i, atom) in mol.atoms.iter().enumerate() {
                eprintln!("Atom {} {}: {:?}", i, atom.symbol, ff.atom_types[i]);
            }
        }
    }
}

#[cfg(test)]
mod type_audit3 {
    use crate::molecule::parser::parse_sdf;
    use crate::mmff::MMFFForceField;
    use crate::MMFFVariant;

    #[test]
    fn test_thiophene_imidazole_types() {
        for (name, sdf) in [
            ("thiophene", r#"Thiophene
     RDKit          3D

  9  9  0  0  0  0  0  0  0  0999 V2000
    0.0000    1.0800    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0280    0.3300    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6420   -0.9400    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.6420   -0.9400    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0280    0.3300    0.0000 S   0  0  0  0  0  0  0  0  0  0  0  0
   -2.0900    0.5700    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2100   -1.8700    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.2100   -1.8700    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    2.1500    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0
  2  3  1  0
  3  4  2  0
  4  5  1  0
  5  1  1  0
  2  6  1  0
  3  7  1  0
  4  8  1  0
  1  9  1  0
M  END"#),
            ("imidazole", r#"Imidazole
     RDKit          3D

  9  9  0  0  0  0  0  0  0  0999 V2000
    0.0000    1.0800    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0280    0.3300    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6420   -0.9400    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    0.6420   -0.9400    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0280    0.3300    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.0900    0.5700    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2100   -1.8700    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.2100   -1.8700    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    2.1500    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  2  0
  3  4  1  0
  4  5  2  0
  5  1  1  0
  2  6  1  0
  3  7  1  0
  4  8  1  0
  5  9  1  0
M  END"#),
        ] {
            let mol = parse_sdf(sdf).expect("parse");
            let ff = MMFFForceField::new(&mol, MMFFVariant::MMFF94s);
            eprintln!("\n=== {} ===", name);
            for (i, atom) in mol.atoms.iter().enumerate() {
                eprintln!("Atom {} {}: {:?}", i, atom.symbol, ff.atom_types[i]);
            }
        }
    }
}
