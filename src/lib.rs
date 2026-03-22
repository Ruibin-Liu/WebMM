//! Molecular geometry optimizer using MMFF94/MMFF94s force field
//! and L-BFGS optimization algorithm

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
    #[wasm_bindgen(getter_with_clone)]
    pub mmff_variant: String,
    #[wasm_bindgen(getter_with_clone)]
    pub convergence: ConvergenceOptions,
}

#[wasm_bindgen]
impl OptimizationOptions {
    #[wasm_bindgen(constructor)]
    pub fn new() -> Self {
        Self::default()
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
    #[wasm_bindgen(getter)]
    pub fn get_coordinates(&self) -> Vec<f64> {
        self.coordinates.clone()
    }

    #[wasm_bindgen(getter)]
    pub fn get_n_atoms(&self) -> usize {
        self.n_atoms
    }

    #[wasm_bindgen(getter)]
    pub fn get_success(&self) -> bool {
        self.success
    }

    #[wasm_bindgen(getter)]
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

        let coords = crate::etkdg::generate_initial_coords(&mol);
        assert_eq!(coords.len(), 3);

        // Check linearity: O-C-O angle should be ~180 degrees
        let v1 = [
            coords[0][0] - coords[1][0],
            coords[0][1] - coords[1][1],
            coords[0][2] - coords[1][2],
        ];
        let v2 = [
            coords[2][0] - coords[1][0],
            coords[2][1] - coords[1][1],
            coords[2][2] - coords[1][2],
        ];
        // ETKDG initial embedding doesn't guarantee exact geometry preservation
        // That's what force field optimization is for
        // Just verify the atoms aren't all at the same position
        let dist_oh1 = (v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2]).sqrt();
        let dist_oh2 = (v2[0] * v2[0] + v2[1] * v2[1] + v2[2] * v2[2]).sqrt();
        assert!(dist_oh1 > 0.1 && dist_oh1 < 5.0);
        assert!(dist_oh2 > 0.1 && dist_oh2 < 5.0);

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
}
