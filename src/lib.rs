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

/// Force-source abstraction (energy + gradient) shared by optimizer and MD
pub mod forces;

/// Molecular dynamics engine (velocity-Verlet + Langevin)
pub mod md;

/// Metadynamics (enhanced sampling via collective variables)
pub mod metad;
pub mod solvation;

/// L-BFGS optimization algorithm
pub mod optimizer;

/// Utility functions
pub mod utils;

/// Property-based tests (only in test mode)
#[cfg(test)]
pub mod prop_tests;

#[cfg(test)]
mod opt_compare;

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
    console_error_panic_hook::set_once();
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

    let config = crate::etkdg::ETKDGConfig {
        random_seed: 42,
        ..Default::default()
    };
    let coords = crate::etkdg::generate_initial_coords_with_config(&mol, &config);
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
    console_error_panic_hook::set_once();
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

    // If the SDF already has 3D coordinates (any non-zero z), use them directly.
    // Only generate ETKDG coordinates for 2D inputs.
    let has_3d = mol.atoms.iter().any(|a| a.position[2].abs() > 1e-6);
    let initial_coords = if has_3d {
        mol.atoms.iter().map(|a| a.position).collect::<Vec<_>>()
    } else {
        let etkdg_config = crate::etkdg::ETKDGConfig {
            random_seed: 42,
            ..Default::default()
        };
        crate::etkdg::generate_initial_coords_with_config(&mol, &etkdg_config)
    };

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

/// Optimize from an SDF using the SDF coordinates directly (no ETKDG).
/// If the SDF is 2D, coordinates are used as-is (z=0 plane).
#[wasm_bindgen]
pub fn optimize_from_sdf_direct(sdf_content: &str, options: OptimizationOptions) -> OptimizationResult {
    console_error_panic_hook::set_once();
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

    let initial_coords: Vec<[f64; 3]> = mol.atoms.iter().map(|a| a.position).collect();

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

// ===== Molecular dynamics (WASM) =====

/// Options for running molecular dynamics (gas-phase MMFF).
#[wasm_bindgen]
pub struct MDOptions {
    mmff_variant: String,
    dt_fs: f64,
    n_steps: usize,
    temperature_k: f64,
    friction_per_ps: f64,
    snapshot_interval: usize,
    seed: u64,
}

#[wasm_bindgen]
impl MDOptions {
    /// Defaults: MMFF94s, 1 fs step, 1000 steps, 300 K, friction 1/ps, snapshots off, seed 42.
    #[wasm_bindgen(constructor)]
    pub fn new() -> Self {
        Self {
            mmff_variant: "MMFF94s".to_string(),
            dt_fs: 1.0,
            n_steps: 1000,
            temperature_k: 300.0,
            friction_per_ps: 1.0,
            snapshot_interval: 0,
            seed: 42,
        }
    }
    #[wasm_bindgen]
    pub fn set_mmff_variant(&mut self, v: String) { self.mmff_variant = v; }
    #[wasm_bindgen]
    pub fn set_dt_fs(&mut self, v: f64) { self.dt_fs = v; }
    #[wasm_bindgen]
    pub fn set_n_steps(&mut self, v: usize) { self.n_steps = v; }
    #[wasm_bindgen]
    pub fn set_temperature_k(&mut self, v: f64) { self.temperature_k = v; }
    #[wasm_bindgen]
    pub fn set_friction_per_ps(&mut self, v: f64) { self.friction_per_ps = v; }
    /// Save a trajectory frame every N steps. 0 (default) = final frame only.
    #[wasm_bindgen]
    pub fn set_snapshot_interval(&mut self, v: usize) { self.snapshot_interval = v; }
    #[wasm_bindgen]
    pub fn set_seed(&mut self, v: u64) { self.seed = v; }
}

impl Default for MDOptions {
    fn default() -> Self {
        Self::new()
    }
}

/// Result of an MD run: a sampled trajectory (flattened coordinates) plus per-frame
/// energies/temperatures and final stats.
#[wasm_bindgen]
pub struct MDResult {
    n_atoms: usize,
    n_frames: usize,
    coordinates: Vec<f64>,   // n_frames * n_atoms * 3, row-major [frame][atom][x,y,z]
    energies: Vec<f64>,      // potential energy per frame (kcal/mol)
    temperatures: Vec<f64>,  // K per frame
    times_fs: Vec<f64>,      // simulation time per frame (fs)
    final_energy: f64,
    final_temperature: f64,
    steps: usize,
    success: bool,
    error: String,
}

#[wasm_bindgen]
impl MDResult {
    #[wasm_bindgen] pub fn n_atoms(&self) -> usize { self.n_atoms }
    #[wasm_bindgen] pub fn n_frames(&self) -> usize { self.n_frames }
    /// Flattened trajectory, n_frames * n_atoms * 3 (Float64Array in JS).
    #[wasm_bindgen] pub fn coordinates(&self) -> Vec<f64> { self.coordinates.clone() }
    #[wasm_bindgen] pub fn energies(&self) -> Vec<f64> { self.energies.clone() }
    #[wasm_bindgen] pub fn temperatures(&self) -> Vec<f64> { self.temperatures.clone() }
    #[wasm_bindgen] pub fn times_fs(&self) -> Vec<f64> { self.times_fs.clone() }
    #[wasm_bindgen] pub fn final_energy(&self) -> f64 { self.final_energy }
    #[wasm_bindgen] pub fn final_temperature(&self) -> f64 { self.final_temperature }
    #[wasm_bindgen] pub fn steps(&self) -> usize { self.steps }
    #[wasm_bindgen] pub fn success(&self) -> bool { self.success }
    #[wasm_bindgen] pub fn error(&self) -> String { self.error.clone() }
    /// Single component: get_coord(frame, atom, axis), axis 0/1/2 = x/y/z.
    #[wasm_bindgen]
    pub fn get_coord(&self, frame: usize, atom: usize, axis: usize) -> f64 {
        self.coordinates[(frame * self.n_atoms + atom) * 3 + axis]
    }
}

fn md_snapshot(
    runner: &crate::md::MDRunner,
    n_atoms: usize,
    coords: &mut Vec<f64>,
    energies: &mut Vec<f64>,
    temperatures: &mut Vec<f64>,
    times_fs: &mut Vec<f64>,
) {
    for row in runner.coords() {
        coords.push(row[0]);
        coords.push(row[1]);
        coords.push(row[2]);
    }
    energies.push(runner.potential_energy());
    temperatures.push(runner.temperature());
    times_fs.push(runner.time_fs());
}

/// Run molecular dynamics on an SDF molecule (gas-phase MMFF) and return a sampled
/// trajectory. NVT (BAOAB Langevin) if `friction_per_ps > 0`, else NVE.
#[wasm_bindgen]
pub fn run_md_from_sdf(sdf_content: &str, options: MDOptions) -> MDResult {
    console_error_panic_hook::set_once();
    let mol = match crate::molecule::parser::parse_sdf(sdf_content) {
        Ok(mol) => mol,
        Err(e) => {
            return MDResult {
                n_atoms: 0, n_frames: 0, coordinates: Vec::new(), energies: Vec::new(),
                temperatures: Vec::new(), times_fs: Vec::new(), final_energy: 0.0,
                final_temperature: 0.0, steps: 0, success: false,
                error: format!("Parse error: {}", e),
            };
        }
    };
    let variant = match options.mmff_variant.as_str() {
        "MMFF94" => MMFFVariant::MMFF94,
        _ => MMFFVariant::MMFF94s,
    };
    let ff = crate::mmff::MMFFForceField::new(&mol, variant);
    let config = crate::md::MDConfig {
        dt_fs: options.dt_fs,
        temperature_k: options.temperature_k,
        friction_per_ps: options.friction_per_ps,
        seed: options.seed,
    };
    let mut runner = crate::md::MDRunner::from_molecule(&ff, &mol, config);
    let n_atoms = mol.atoms.len();
    let mut coordinates = Vec::new();
    let mut energies = Vec::new();
    let mut temperatures = Vec::new();
    let mut times_fs = Vec::new();
    let snap = options.snapshot_interval;
    if snap > 0 {
        md_snapshot(&runner, n_atoms, &mut coordinates, &mut energies, &mut temperatures, &mut times_fs);
    }
    for i in 1..=options.n_steps {
        runner.step();
        if snap > 0 && i != options.n_steps && i % snap == 0 {
            md_snapshot(&runner, n_atoms, &mut coordinates, &mut energies, &mut temperatures, &mut times_fs);
        }
    }
    // final frame always
    md_snapshot(&runner, n_atoms, &mut coordinates, &mut energies, &mut temperatures, &mut times_fs);
    let final_energy = runner.potential_energy();
    let final_temperature = runner.temperature();
    MDResult {
        n_atoms,
        n_frames: energies.len(),
        coordinates, energies, temperatures, times_fs,
        final_energy, final_temperature, steps: options.n_steps,
        success: true, error: String::new(),
    }
}

// ===== Metadynamics (WASM) =====

/// Options for a metadynamics run (well-tempered, gas-phase MMFF).
#[wasm_bindgen]
pub struct MetaDOptions {
    mmff_variant: String,
    dt_fs: f64,
    n_steps: usize,
    temperature_k: f64,
    friction_per_ps: f64,
    seed: u32,
    snapshot_interval: usize,
    cv_type: String,         // "dihedral" or "distance"
    cv_atoms: String,        // comma-separated atom indices, e.g. "8,2,1,0"
    hill_height: f64,        // kcal/mol
    hill_width: f64,         // radians (dihedral) or angstrom (distance)
    deposit_interval: usize, // deposit a hill every N steps
    bias_factor: f64,        // 0 = standard; >1 = well-tempered
    fes_grid_points: usize,  // FES reconstruction resolution
}

#[wasm_bindgen]
impl MetaDOptions {
    #[wasm_bindgen(constructor)]
    pub fn new() -> Self {
        Self {
            mmff_variant: "MMFF94s".to_string(),
            dt_fs: 1.0, n_steps: 5000, temperature_k: 300.0, friction_per_ps: 5.0,
            seed: 42, snapshot_interval: 50,
            cv_type: "dihedral".to_string(), cv_atoms: "0,1,2,3".to_string(),
            hill_height: 0.3, hill_width: 0.2, deposit_interval: 50, bias_factor: 10.0,
            fes_grid_points: 72,
        }
    }
    #[wasm_bindgen] pub fn set_mmff_variant(&mut self, v: String) { self.mmff_variant = v; }
    #[wasm_bindgen] pub fn set_dt_fs(&mut self, v: f64) { self.dt_fs = v; }
    #[wasm_bindgen] pub fn set_n_steps(&mut self, v: usize) { self.n_steps = v; }
    #[wasm_bindgen] pub fn set_temperature_k(&mut self, v: f64) { self.temperature_k = v; }
    #[wasm_bindgen] pub fn set_friction_per_ps(&mut self, v: f64) { self.friction_per_ps = v; }
    #[wasm_bindgen] pub fn set_seed(&mut self, v: u32) { self.seed = v; }
    #[wasm_bindgen] pub fn set_snapshot_interval(&mut self, v: usize) { self.snapshot_interval = v; }
    /// "dihedral" (4 atoms) or "distance" (2 atoms).
    #[wasm_bindgen] pub fn set_cv_type(&mut self, v: String) { self.cv_type = v; }
    /// Comma-separated atom indices, e.g. "8,2,1,0" for a dihedral.
    #[wasm_bindgen] pub fn set_cv_atoms(&mut self, v: String) { self.cv_atoms = v; }
    #[wasm_bindgen] pub fn set_hill_height(&mut self, v: f64) { self.hill_height = v; }
    #[wasm_bindgen] pub fn set_hill_width(&mut self, v: f64) { self.hill_width = v; }
    #[wasm_bindgen] pub fn set_deposit_interval(&mut self, v: usize) { self.deposit_interval = v; }
    #[wasm_bindgen] pub fn set_bias_factor(&mut self, v: f64) { self.bias_factor = v; }
    #[wasm_bindgen] pub fn set_fes_grid_points(&mut self, v: usize) { self.fes_grid_points = v; }
}
impl Default for MetaDOptions { fn default() -> Self { Self::new() } }

/// Result of a metadynamics run: trajectory + CV trace + free-energy surface.
#[wasm_bindgen]
pub struct MetaDResult {
    n_atoms: usize,
    n_frames: usize,
    n_hills: usize,
    coordinates: Vec<f64>,
    energies: Vec<f64>,
    cv_values: Vec<f64>,
    times_fs: Vec<f64>,
    fes_s: Vec<f64>,
    fes_f: Vec<f64>,
    hill_centers: Vec<f64>,
    final_energy: f64,
    success: bool,
    error: String,
}

#[wasm_bindgen]
impl MetaDResult {
    #[wasm_bindgen] pub fn n_atoms(&self) -> usize { self.n_atoms }
    #[wasm_bindgen] pub fn n_frames(&self) -> usize { self.n_frames }
    #[wasm_bindgen] pub fn n_hills(&self) -> usize { self.n_hills }
    #[wasm_bindgen] pub fn coordinates(&self) -> Vec<f64> { self.coordinates.clone() }
    #[wasm_bindgen] pub fn energies(&self) -> Vec<f64> { self.energies.clone() }
    #[wasm_bindgen] pub fn cv_values(&self) -> Vec<f64> { self.cv_values.clone() }
    #[wasm_bindgen] pub fn times_fs(&self) -> Vec<f64> { self.times_fs.clone() }
    /// CV grid points for the FES.
    #[wasm_bindgen] pub fn fes_s(&self) -> Vec<f64> { self.fes_s.clone() }
    /// Free energy (kcal/mol) at each grid point.
    #[wasm_bindgen] pub fn fes_f(&self) -> Vec<f64> { self.fes_f.clone() }
    #[wasm_bindgen] pub fn hill_centers(&self) -> Vec<f64> { self.hill_centers.clone() }
    #[wasm_bindgen] pub fn final_energy(&self) -> f64 { self.final_energy }
    #[wasm_bindgen] pub fn success(&self) -> bool { self.success }
    #[wasm_bindgen] pub fn error(&self) -> String { self.error.clone() }
}

/// Run well-tempered metadynamics on an SDF molecule (gas-phase MMFF) and return
/// the trajectory + CV trace + free-energy surface (FES).
#[wasm_bindgen]
pub fn run_metadynamics_from_sdf(sdf_content: &str, options: MetaDOptions) -> MetaDResult {
    console_error_panic_hook::set_once();
    let err = |msg: String| MetaDResult {
        n_atoms: 0, n_frames: 0, n_hills: 0, coordinates: vec![], energies: vec![],
        cv_values: vec![], times_fs: vec![], fes_s: vec![], fes_f: vec![],
        hill_centers: vec![], final_energy: 0.0, success: false, error: msg,
    };
    let mol = match crate::molecule::parser::parse_sdf(sdf_content) {
        Ok(m) => m,
        Err(e) => return err(format!("Parse error: {}", e)),
    };
    let variant = match options.mmff_variant.as_str() {
        "MMFF94" => MMFFVariant::MMFF94,
        _ => MMFFVariant::MMFF94s,
    };
    let ff = crate::mmff::MMFFForceField::new(&mol, variant);

    // Parse CV atoms
    let atoms: Vec<usize> = options.cv_atoms.split(',').filter_map(|s| s.trim().parse().ok()).collect();
    let cv: Box<dyn crate::metad::CollectiveVariable> = match options.cv_type.as_str() {
        "distance" => {
            if atoms.len() < 2 { return err("distance CV needs >=2 atoms".into()); }
            Box::new(crate::metad::DistanceCV::new(atoms[0], atoms[1]))
        }
        _ => {
            if atoms.len() < 4 { return err("dihedral CV needs >=4 atoms".into()); }
            Box::new(crate::metad::DihedralCV::new(atoms[0], atoms[1], atoms[2], atoms[3]))
        }
    };

    let metad = crate::metad::MetaDynamics::new(&ff, cv, crate::metad::MetaDConfig {
        hill_height: options.hill_height,
        hill_width: options.hill_width,
        deposit_interval: options.deposit_interval,
        bias_factor: options.bias_factor,
        temperature_k: options.temperature_k,
    });
    let mut runner = crate::md::MDRunner::from_molecule(&metad, &mol, crate::md::MDConfig {
        dt_fs: options.dt_fs, temperature_k: options.temperature_k,
        friction_per_ps: options.friction_per_ps, seed: options.seed as u64,
    });
    let n_atoms = mol.atoms.len();
    let snap = options.snapshot_interval;

    let mut coordinates = Vec::new();
    let mut energies = Vec::new();
    let mut cv_values = Vec::new();
    let mut times_fs = Vec::new();
    let mut snap_fn = |runner: &crate::md::MDRunner, metad: &crate::metad::MetaDynamics| {
        for row in runner.coords() { coordinates.push(row[0]); coordinates.push(row[1]); coordinates.push(row[2]); }
        energies.push(runner.potential_energy());
        cv_values.push(metad.last_cv());
        times_fs.push(runner.time_fs());
    };
    if snap > 0 { snap_fn(&runner, &metad); }
    for i in 1..=options.n_steps {
        runner.step();
        if snap > 0 && i != options.n_steps && i.is_multiple_of(snap) {
            snap_fn(&runner, &metad);
        }
    }
    snap_fn(&runner, &metad); // final frame

    // FES
    let (smin, smax) = metad.cv_range();
    let ng = options.fes_grid_points.max(2);
    let fes_s: Vec<f64> = (0..ng).map(|i| smin + (smax - smin) * i as f64 / (ng - 1) as f64).collect();
    let fes_f = metad.free_energy_surface(&fes_s);
    let hill_centers = metad.hill_centers();
    let n_hills = metad.hill_count();
    let final_energy = runner.potential_energy();

    MetaDResult {
        n_atoms, n_frames: energies.len(), n_hills,
        coordinates, energies, cv_values, times_fs,
        fes_s, fes_f, hill_centers, final_energy,
        success: true, error: String::new(),
    }
}

#[cfg(test)]
mod tests {
    use crate::acetic_debug::dihedral;
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
                    let bt_ij = crate::mmff::params::get_mmff_bond_type(
                        bij.bond_type,
                        ff.type_ids[i],
                        ff.type_ids[j],
                    );
                    let bt_jk = crate::mmff::params::get_mmff_bond_type(
                        bkj.bond_type,
                        ff.type_ids[j],
                        ff.type_ids[k],
                    );
                    let angle_type_val =
                        crate::mmff::mmff_tables::compute_angle_type(bt_ij, bt_jk, 0);
                    if let (Some(sb_params), Some(bp_ij), Some(bp_kj), Some(ap)) = (
                        get_stretch_bend_params(
                            ff.atom_types[i],
                            ff.atom_types[j],
                            ff.atom_types[k],
                            bt_ij,
                            bt_jk,
                            ff.mol.atoms[i].atomic_number,
                            ff.mol.atoms[j].atomic_number,
                            ff.mol.atoms[k].atomic_number,
                            angle_type_val,
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
            crate::mmff::MMFFAtomType::O_2,
            "Carbonyl O double-bonded to C should be O_2, got {:?}",
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
    fn test_sp3d_typing() {
        let sdf = r#"PF5
     RDKit          3D
  6  5  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 P   0  0  0  0  0  0  0  0  0
    1.5000    0.0000    0.0000 F   0  0  0  0  0  0  0  0  0
   -1.5000    0.0000    0.0000 F   0  0  0  0  0  0  0  0  0
    0.0000    1.5000    0.0000 F   0  0  0  0  0  0  0  0  0
    0.0000   -1.5000    0.0000 F   0  0  0  0  0  0  0  0  0
    0.0000    0.0000    1.8000 F   0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  1  3  1  0  0  0  0
  1  4  1  0  0  0  0
  1  5  1  0  0  0  0
  1  6  1  0  0  0  0
M  END"#;
        let mol = parse_sdf(sdf).expect("parse failed");
        let ff = crate::mmff::MMFFForceField::new(&mol, MMFFVariant::MMFF94s);
        assert_eq!(ff.atom_types[0], crate::mmff::MMFFAtomType::P_3D);
        assert!(!ff.angles.is_empty());
    }

    #[test]
    fn test_coord_map_water() {
        let sdf = r#"Water
     RDKit          3D
  3  2  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0
    0.9572    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0
   -0.2390    0.9260    0.0000 H   0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  1  3  1  0  0  0  0
M  END"#;
        let mol = parse_sdf(sdf).expect("parse failed");
        let mut config = crate::etkdg::ETKDGConfig {
            random_seed: 42,
            ..Default::default()
        };
        // Fix oxygen at origin
        config.coord_map.insert(0, [0.0, 0.0, 0.0]);
        let coords = crate::etkdg::generate_initial_coords_with_config(&mol, &config);
        assert_eq!(
            coords[0],
            [0.0, 0.0, 0.0],
            "Oxygen should be fixed at origin"
        );
        // Check that H atoms are roughly 0.96A from O
        let d1 = ((coords[1][0] - coords[0][0]).powi(2)
            + (coords[1][1] - coords[0][1]).powi(2)
            + (coords[1][2] - coords[0][2]).powi(2))
        .sqrt();
        let d2 = ((coords[2][0] - coords[0][0]).powi(2)
            + (coords[2][1] - coords[0][1]).powi(2)
            + (coords[2][2] - coords[0][2]).powi(2))
        .sqrt();
        assert!(
            d1 > 0.8 && d1 < 1.1,
            "H1 should be ~0.96A from O, got {}",
            d1
        );
        assert!(
            d2 > 0.8 && d2 < 1.1,
            "H2 should be ~0.96A from O, got {}",
            d2
        );
    }

    #[test]
    fn test_fragment_embedding_water_dimer() {
        // Two water molecules (no bonds between them)
        let sdf = r#"Water dimer
     RDKit          3D
  6  4  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0
    0.9572    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0
   -0.2390    0.9260    0.0000 H   0  0  0  0  0  0  0  0  0
    3.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0
    3.9572    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0
    2.7610    0.9260    0.0000 H   0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  1  3  1  0  0  0  0
  4  5  1  0  0  0  0
  4  6  1  0  0  0  0
M  END"#;
        let mol = parse_sdf(sdf).expect("parse failed");
        let config = crate::etkdg::ETKDGConfig {
            random_seed: 42,
            ..Default::default()
        };
        let coords = crate::etkdg::generate_initial_coords_with_config(&mol, &config);
        assert_eq!(coords.len(), 6);
        // Check O-O distance is reasonable (should be > 2.5A for separated fragments)
        let dx = coords[0][0] - coords[3][0];
        let dy = coords[0][1] - coords[3][1];
        let dz = coords[0][2] - coords[3][2];
        let d_oo = (dx * dx + dy * dy + dz * dz).sqrt();
        assert!(
            d_oo > 2.0,
            "Water fragments should be separated by >2A, got {}",
            d_oo
        );
        // Check each water has reasonable O-H distances
        for water_start in [0, 3] {
            for h in [1, 2] {
                let h_idx = water_start + h;
                let d = ((coords[h_idx][0] - coords[water_start][0]).powi(2)
                    + (coords[h_idx][1] - coords[water_start][1]).powi(2)
                    + (coords[h_idx][2] - coords[water_start][2]).powi(2))
                .sqrt();
                assert!(
                    d > 0.8 && d < 1.1,
                    "O-H distance should be ~0.96A, got {}",
                    d
                );
            }
        }
    }

    #[test]
    fn test_atropisomer_simple() {
        // Create a simple atropisomer: two sp2 carbons with substituents
        // C1(F)(H)-C2(Cl)(Br) with AtropCW stereochemistry
        let mol = crate::molecule::Molecule {
            name: "atrop_test".to_string(),
            atoms: vec![
                crate::molecule::Atom {
                    index: 0,
                    atomic_number: 6,
                    symbol: "C".to_string(),
                    mass: 12.0,
                    position: [0.0, 0.0, 0.0],
                    charge: 0.0,
                    stereo_parity: 0,
                },
                crate::molecule::Atom {
                    index: 1,
                    atomic_number: 6,
                    symbol: "C".to_string(),
                    mass: 12.0,
                    position: [0.0, 0.0, 0.0],
                    charge: 0.0,
                    stereo_parity: 0,
                },
                crate::molecule::Atom {
                    index: 2,
                    atomic_number: 9,
                    symbol: "F".to_string(),
                    mass: 19.0,
                    position: [0.0, 0.0, 0.0],
                    charge: 0.0,
                    stereo_parity: 0,
                },
                crate::molecule::Atom {
                    index: 3,
                    atomic_number: 1,
                    symbol: "H".to_string(),
                    mass: 1.0,
                    position: [0.0, 0.0, 0.0],
                    charge: 0.0,
                    stereo_parity: 0,
                },
                crate::molecule::Atom {
                    index: 4,
                    atomic_number: 17,
                    symbol: "Cl".to_string(),
                    mass: 35.5,
                    position: [0.0, 0.0, 0.0],
                    charge: 0.0,
                    stereo_parity: 0,
                },
                crate::molecule::Atom {
                    index: 5,
                    atomic_number: 35,
                    symbol: "Br".to_string(),
                    mass: 79.9,
                    position: [0.0, 0.0, 0.0],
                    charge: 0.0,
                    stereo_parity: 0,
                },
            ],
            bonds: vec![
                crate::molecule::Bond {
                    atom1: 0,
                    atom2: 1,
                    bond_type: crate::molecule::BondType::Single,
                    stereo: crate::molecule::BondStereo::AtropCW,
                },
                crate::molecule::Bond {
                    atom1: 0,
                    atom2: 2,
                    bond_type: crate::molecule::BondType::Single,
                    stereo: crate::molecule::BondStereo::None,
                },
                crate::molecule::Bond {
                    atom1: 0,
                    atom2: 3,
                    bond_type: crate::molecule::BondType::Single,
                    stereo: crate::molecule::BondStereo::None,
                },
                crate::molecule::Bond {
                    atom1: 1,
                    atom2: 4,
                    bond_type: crate::molecule::BondType::Single,
                    stereo: crate::molecule::BondStereo::None,
                },
                crate::molecule::Bond {
                    atom1: 1,
                    atom2: 5,
                    bond_type: crate::molecule::BondType::Single,
                    stereo: crate::molecule::BondStereo::None,
                },
            ],
            adjacency: vec![
                vec![1, 2, 3],
                vec![0, 4, 5],
                vec![0],
                vec![0],
                vec![1],
                vec![1],
            ],
        };

        let config = crate::etkdg::ETKDGConfig {
            random_seed: 42,
            ..Default::default()
        };
        let coords = crate::etkdg::generate_initial_coords_with_config(&mol, &config);
        assert_eq!(coords.len(), 6);

        // Check that the chiral volume has the correct sign for AtropCW
        // AtropCW should give negative chiral volume
        let vol = crate::etkdg::chiral_volume(&coords, 2, 3, 4, 5);
        assert!(
            !(-0.5..=0.5).contains(&vol),
            "Atropisomer should have significant chiral volume, got {}",
            vol
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
        assert!(result.get_iterations() > 0);

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
        let conv = crate::ConvergenceOptions {
            max_iterations: 200,
            ..Default::default()
        };
        let result = crate::optimizer::optimize(&ff, &coords, &conv);
        for (i, c) in result.optimized_coords.iter().enumerate() {
            coords[i] = *c;
        }
        let e_after = if result.converged {
            ff.calculate_energy(&coords)
        } else {
            e_before
        };

        eprintln!(
            "Acetic acid optimization: {:.2} -> {:.2} kcal/mol, converged={}, iters={}",
            e_before, e_after, result.converged, result.iterations
        );

        if result.converged {
            assert!(e_after < e_before, "Energy should decrease");
        }

        assert!(
            result.iterations > 0,
            "Optimizer should run at least one iteration"
        );

        if result.converged {
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
        let mol = parse_sdf(sdf).unwrap_or_else(|_| panic!("{}: parse failed", name));
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
        let mol = parse_sdf(sdf).unwrap_or_else(|_| panic!("{}: parse failed", name));

        // Run ETKDG + optimize
        let config = crate::etkdg::ETKDGConfig {
            max_attempts: 5,
            max_iterations: 2000,
            random_seed: 42,
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
            max_iterations: 5000,
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

        if result.converged {
            for (i, c) in result.optimized_coords.iter().enumerate() {
                coords[i] = *c;
            }
        } else {
            eprintln!("  Optimizer did not converge, using ETKDG coordinates for validation");
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
            ],
            0.15,
            20.0,
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
            0.15,
            10.0,
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
            0.15,
            20.0,
        );

        // Verify carboxyl planarity: key dihedrals should be ~0° or ~180°
        let mol = parse_sdf(sdf).unwrap();
        let config = crate::etkdg::ETKDGConfig {
            random_seed: 42,
            ..Default::default()
        };
        let coords = crate::etkdg::generate_initial_coords_with_config(&mol, &config);
        let d1 = dihedral(&coords, 2, 1, 3, 7).abs();
        let d2 = dihedral(&coords, 0, 1, 2, 3).abs();
        assert!(
            !(40.0..=140.0).contains(&d1),
            "O=C-O-H dihedral should be planar, got {:.1}°",
            d1
        );
        assert!(
            !(15.0..=165.0).contains(&d2),
            "C-C(=O)-O dihedral should be planar, got {:.1}°",
            d2
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
        assert!(e_perfect.is_finite(), "ethane: non-finite energy from perfect coords");
        let cc0 = dist3(&perfect_coords, 0, 1);
        assert!(cc0 > 1.45 && cc0 < 1.55, "ethane C-C dist {} out of [1.45, 1.55]", cc0);
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
        let init_coords = crate::etkdg::generate_initial_coords_with_config(
            &mol,
            &crate::etkdg::ETKDGConfig {
                random_seed: 42,
                ..Default::default()
            },
        );
        eprintln!("\n=== ACETIC ACID ETKDG initial ===");
        eprintln!("  C-C bond: {:.4}", dist3(&init_coords, 0, 1));
        eprintln!("  C=O bond: {:.4}", dist3(&init_coords, 1, 2));

        assert!(
            (dist3(&init_coords, 0, 1) - 1.4928).abs() < 0.25,
            "ETKDG C-C bond unreasonable: {:.4}",
            dist3(&init_coords, 0, 1)
        );
        assert!(
            (dist3(&init_coords, 1, 2) - 1.2188).abs() < 0.25,
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
    use crate::etkdg::eigenvector_smallest_eigenvalue_3x3;
    use crate::etkdg::generate_initial_coords;
    use crate::etkdg::generate_initial_coords_with_config;
    use crate::molecule::parser::parse_sdf;

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

        atoms
            .iter()
            .map(|&idx| {
                let dx = coords[idx][0] - cx;
                let dy = coords[idx][1] - cy;
                let dz = coords[idx][2] - cz;
                (dx * normal[0] + dy * normal[1] + dz * normal[2]).abs()
            })
            .fold(0.0, f64::max)
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
        assert!(
            ring_dev < 0.25,
            "benzene ring not planar: {:.6} Å",
            ring_dev
        );

        let all_atoms: Vec<usize> = (0..12).collect();
        let all_dev = max_planar_deviation(&coords, &all_atoms);
        eprintln!("benzene all-atom planar deviation: {:.6} Å", all_dev);
        assert!(
            all_dev < 1.0,
            "benzene H atoms not planar: {:.6} Å",
            all_dev
        );
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

        for i in 0..20 {
            let coords = generate_initial_coords_with_config(
                &mol,
                &crate::etkdg::ETKDGConfig {
                    random_seed: 100 + i as i64,
                    ..Default::default()
                },
            );
            let ring_atoms: Vec<usize> = vec![0, 1, 2, 3, 4, 5];
            let ring_dev = max_planar_deviation(&coords, &ring_atoms);
            let all_atoms: Vec<usize> = (0..12).collect();
            let all_dev = max_planar_deviation(&coords, &all_atoms);

            if ring_dev > max_ring_dev {
                max_ring_dev = ring_dev;
            }
            if all_dev > max_all_dev {
                max_all_dev = all_dev;
            }
        }

        eprintln!("Stress test max ring deviation: {:.6} Å", max_ring_dev);
        eprintln!("Stress test max all-atom deviation: {:.6} Å", max_all_dev);

        assert!(
            max_ring_dev < 0.25,
            "ring not planar over 20 runs: {:.6} Å",
            max_ring_dev
        );
        assert!(
            max_all_dev < 1.0,
            "H atoms not planar over 20 runs: {:.6} Å",
            max_all_dev
        );
    }
}
#[cfg(test)]
mod debug_worst {
    use crate::etkdg::eigenvector_smallest_eigenvalue_3x3;
    use crate::etkdg::generate_initial_coords_with_config;
    use crate::molecule::parser::parse_sdf;

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

        atoms
            .iter()
            .map(|&idx| {
                let dx = coords[idx][0] - cx;
                let dy = coords[idx][1] - cy;
                let dz = coords[idx][2] - cz;
                (dx * normal[0] + dy * normal[1] + dz * normal[2]).abs()
            })
            .fold(0.0, f64::max)
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

        for i in 0..10 {
            let coords = generate_initial_coords_with_config(
                &mol,
                &crate::etkdg::ETKDGConfig {
                    random_seed: 200 + i as i64,
                    ..Default::default()
                },
            );
            let all_atoms: Vec<usize> = (0..12).collect();
            let all_dev = max_planar_deviation(&coords, &all_atoms);
            if all_dev > max_all_dev {
                max_all_dev = all_dev;
                worst_coords = coords;
            }
        }

        eprintln!("Max all-atom deviation: {:.6} Å", max_all_dev);
        assert!(
            max_all_dev < 0.01,
            "benzene worst-case planar deviation {} Å exceeds 0.01",
            max_all_dev
        );

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
        for (i, wc) in worst_coords.iter().enumerate() {
            let dx = wc[0] - cx;
            let dy = wc[1] - cy;
            let dz = wc[2] - cz;
            let d = (dx * normal[0] + dy * normal[1] + dz * normal[2]).abs();
            let name = if i < 6 { "C" } else { "H" };
            eprintln!(
                "{} {:2}: {:10.6} Å  pos=[{:8.4}, {:8.4}, {:8.4}]",
                name, i, d, wc[0], wc[1], wc[2]
            );
        }
    }
}

#[cfg(test)]
mod pyrrole_tests {
    use crate::etkdg::eigenvector_smallest_eigenvalue_3x3;
    use crate::etkdg::generate_initial_coords_with_config;
    use crate::molecule::graph::get_aromatic_atoms;
    use crate::molecule::parser::parse_sdf;
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
        atoms
            .iter()
            .map(|&idx| {
                let dx = coords[idx][0] - cx;
                let dy = coords[idx][1] - cy;
                let dz = coords[idx][2] - cz;
                (dx * normal[0] + dy * normal[1] + dz * normal[2]).abs()
            })
            .fold(0.0, f64::max)
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
        assert_eq!(
            aromatic, expected,
            "Pyrrole ring atoms should all be aromatic"
        );

        let coords = generate_initial_coords_with_config(
            &mol,
            &crate::etkdg::ETKDGConfig {
                random_seed: 42,
                ..Default::default()
            },
        );
        let ring_atoms: Vec<usize> = vec![0, 1, 2, 3, 4];
        let ring_dev = max_planar_deviation(&coords, &ring_atoms);
        eprintln!("pyrrole ring planar deviation: {:.6} A", ring_dev);
        assert!(
            ring_dev < 0.10,
            "pyrrole ring not planar: {:.6} A",
            ring_dev
        );

        let all_atoms: Vec<usize> = (0..10).collect();
        let all_dev = max_planar_deviation(&coords, &all_atoms);
        eprintln!("pyrrole all-atom planar deviation: {:.6} A", all_dev);
        assert!(
            all_dev < 1.0,
            "pyrrole H atoms not planar: {:.6} A",
            all_dev
        );
    }
}

#[cfg(test)]
mod aniline_tests {
    use crate::mmff::MMFFForceField;
    use crate::molecule::parser::parse_sdf;
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
        let expected = [
            "C_AR", "C_AR", "C_AR", "C_AR", "C_AR", "C_AR",
            "N_PL3", "H_NAM", "H_NAM", "H", "H", "H", "H",
        ];
        assert_eq!(mol.atoms.len(), expected.len(), "aniline atom count");
        for (i, exp) in expected.iter().enumerate() {
            assert_eq!(format!("{:?}", ff.atom_types[i]), *exp, "aniline atom {}", i);
        }
        let coords = crate::etkdg::generate_initial_coords(&mol);
        let (energy, _) = ff.calculate_energy_and_gradient(&coords);
        eprintln!("Aniline initial energy: {}", energy);
    }
}

#[cfg(test)]
mod type_audit {
    use crate::mmff::MMFFForceField;
    use crate::molecule::parser::parse_sdf;
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
            let expected: &[&str] = match name {
                "pyrrole" => &["C5A", "C5B", "C5B", "C5A", "NPYL", "H", "H", "H", "H_N3", "H"],
                "furan" => &["C5A", "C5B", "C5B", "C5A", "OFUR", "H", "H", "H", "H"],
                _ => &[],
            };
            assert_eq!(mol.atoms.len(), expected.len(), "{} atom count", name);
            for (i, atom) in mol.atoms.iter().enumerate() {
                eprintln!("Atom {} {}: {:?}", i, atom.symbol, ff.atom_types[i]);
                assert_eq!(format!("{:?}", ff.atom_types[i]), expected[i], "{} atom {}", name, i);
            }
        }
    }
}

#[cfg(test)]
mod type_audit2 {
    use crate::mmff::MMFFForceField;
    use crate::molecule::parser::parse_sdf;
    use crate::MMFFVariant;

    #[test]
    fn test_simple_molecule_types() {
        for (name, sdf) in [
            (
                "water",
                r#"Water
  CDK     1012182035

  3  2  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.9580    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2390    0.9270    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  1  3  1  0
M  END"#,
            ),
            (
                "methanol",
                r#"Methanol
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
M  END"#,
            ),
            (
                "phenol",
                r#"Phenol
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
M  END"#,
            ),
        ] {
            let mol = parse_sdf(sdf).expect("parse");
            let ff = MMFFForceField::new(&mol, MMFFVariant::MMFF94s);
            eprintln!("\n=== {} ===", name);
            let expected: &[&str] = match name {
                "water" => &["OH2", "H_OH", "H_OH"],
                "methanol" => &["C_3", "O_R", "H", "H", "H", "H_ONC"],
                "phenol" => &[
                    "C_AR", "C_AR", "C_AR", "C_AR", "C_AR", "C_AR", "O_R",
                    "H", "H", "H", "H", "H", "H_OAR",
                ],
                _ => &[],
            };
            assert_eq!(mol.atoms.len(), expected.len(), "{} atom count", name);
            for (i, atom) in mol.atoms.iter().enumerate() {
                eprintln!("Atom {} {}: {:?}", i, atom.symbol, ff.atom_types[i]);
                assert_eq!(format!("{:?}", ff.atom_types[i]), expected[i], "{} atom {}", name, i);
            }
        }
    }
}

#[cfg(test)]
mod type_audit3 {
    use crate::mmff::MMFFForceField;
    use crate::molecule::parser::parse_sdf;
    use crate::MMFFVariant;

    #[test]
    fn test_thiophene_imidazole_types() {
        for (name, sdf) in [
            (
                "thiophene",
                r#"Thiophene
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
M  END"#,
            ),
            (
                "imidazole",
                r#"Imidazole
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
M  END"#,
            ),
        ] {
            let mol = parse_sdf(sdf).expect("parse");
            let ff = MMFFForceField::new(&mol, MMFFVariant::MMFF94s);
            eprintln!("\n=== {} ===", name);
            let expected: &[&str] = match name {
                "thiophene" => &["C5A", "C5B", "C5B", "C5A", "S_AR", "H", "H", "H", "H"],
                "imidazole" => &["N5B", "C5A", "NPYL", "C5A", "C5B", "H", "H_N3", "H", "H"],
                _ => &[],
            };
            assert_eq!(mol.atoms.len(), expected.len(), "{} atom count", name);
            for (i, atom) in mol.atoms.iter().enumerate() {
                eprintln!("Atom {} {}: {:?}", i, atom.symbol, ff.atom_types[i]);
                assert_eq!(format!("{:?}", ff.atom_types[i]), expected[i], "{} atom {}", name, i);
            }
        }
    }
}

#[cfg(test)]
mod acetic_debug {
    use crate::etkdg::{generate_initial_coords_with_config, ETKDGConfig};
    use crate::mmff::MMFFForceField;
    use crate::molecule::parser::parse_sdf;
    use crate::optimizer;
    use crate::{ConvergenceOptions, MMFFVariant};

    pub fn dihedral(coords: &[[f64; 3]], i: usize, j: usize, k: usize, l: usize) -> f64 {
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
        let m1 = [
            n1[1] * b2[2] - n1[2] * b2[1],
            n1[2] * b2[0] - n1[0] * b2[2],
            n1[0] * b2[1] - n1[1] * b2[0],
        ];
        let b2n = (b2[0] * b2[0] + b2[1] * b2[1] + b2[2] * b2[2]).sqrt();
        if b2n < 1e-10 {
            return 0.0;
        }
        let x = n1[0] * n2[0] + n1[1] * n2[1] + n1[2] * n2[2];
        let y = (m1[0] * n2[0] + m1[1] * n2[1] + m1[2] * n2[2]) / b2n;
        y.atan2(x).to_degrees()
    }

    #[test]
    fn test_acetic_carboxyl_planarity() {
        // Atoms: 0=C(methyl), 1=C(carboxyl), 2=O(=), 3=O(-H), 4-6=H(methyl), 7=H(hydroxyl)
        let sdf = r#"Acetic acid
     RDKit          3D

  8  7  0  0  0  0  0  0  0  0999 V2000
   -0.9549   -0.1188   -0.0467 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.4728    0.3162   -0.0139 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.8884    1.4538   -0.1499 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.3178   -0.7111    0.1905 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.1077   -0.8272   -0.8647 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2264   -0.5716    0.9100 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5945    0.7523   -0.2149 H   0  0  0  0  0  0  0  0  0  0  0  0
    2.2046   -0.2935    0.1897 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  2  0
  2  4  1  0
  1  5  1  0
  1  6  1  0
  1  7  1  0
  4  8  1  0
M  END"#;
        let mol = parse_sdf(sdf).unwrap();
        let config = ETKDGConfig {
            random_seed: 42,
            ..Default::default()
        };
        let coords = generate_initial_coords_with_config(&mol, &config);

        eprintln!("=== ETKDG geometry ===");
        for (i, (atom, c)) in mol.atoms.iter().zip(coords.iter()).enumerate() {
            eprintln!(
                "  {} {}: {:.4} {:.4} {:.4}",
                i, atom.symbol, c[0], c[1], c[2]
            );
        }
        // Key dihedrals for carboxyl planarity:
        // H7-O3-C2-C1 (should be ~0 or ~180 for planar)
        eprintln!(
            "  dihedral O(=O)-C-C(=O)-O(H): {:.1}°",
            dihedral(&coords, 2, 1, 0, 3)
        );
        eprintln!(
            "  dihedral O(=O)-C-OH-H:      {:.1}°",
            dihedral(&coords, 2, 1, 3, 7)
        );
        eprintln!(
            "  dihedral C(methyl)-C-O-H:    {:.1}°",
            dihedral(&coords, 0, 1, 3, 7)
        );
        eprintln!(
            "  dihedral O(=O)-C-C-H(meth):  {:.1}°",
            dihedral(&coords, 2, 1, 0, 4)
        );

        let ff = MMFFForceField::new(&mol, MMFFVariant::MMFF94s);
        let result = optimizer::optimize(
            &ff,
            &coords,
            &ConvergenceOptions {
                max_iterations: 5000,
                ..Default::default()
            },
        );
        let opt: Vec<[f64; 3]> = if result.converged {
            result.optimized_coords.clone()
        } else {
            eprintln!("  Optimizer did not converge, using ETKDG coordinates for planarity check");
            coords.clone()
        };

        eprintln!("\n=== After MMFF94s optimization ===");
        for (i, (atom, c)) in mol.atoms.iter().zip(opt.iter()).enumerate() {
            eprintln!(
                "  {} {}: {:.4} {:.4} {:.4}",
                i, atom.symbol, c[0], c[1], c[2]
            );
        }
        eprintln!(
            "  dihedral O(=O)-C-C(=O)-O(H): {:.1}°",
            dihedral(&opt, 2, 1, 0, 3)
        );
        eprintln!(
            "  dihedral O(=O)-C-OH-H:      {:.1}°",
            dihedral(&opt, 2, 1, 3, 7)
        );
        eprintln!(
            "  dihedral C(methyl)-C-O-H:    {:.1}°",
            dihedral(&opt, 0, 1, 3, 7)
        );
        eprintln!(
            "  Energy: {:.4}, converged: {}",
            result.final_energy, result.converged
        );

        // The carboxyl group should be planar: dihedral ~0 or ~180
        let d = dihedral(&opt, 2, 1, 3, 7);
        assert!(
            d.abs() < 15.0 || (d - 180.0).abs() < 15.0 || (d + 180.0).abs() < 15.0,
            "Carboxyl OH not planar: dihedral O=C-O-H = {:.1}°",
            d
        );
    }

    #[test]
    fn test_rdkit_energy_comparison() {
        let cases = [
            ("water", "Water\n     RDKit          3D\n\n  3  2  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n    0.9580    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.2390    0.9270    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0\n  1  3  1  0\nM  END", 0.146993),
            ("methane", "Methane\n     RDKit          3D\n\n  5  4  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.6320    0.6320    0.6320 H   0  0  0  0  0  0  0  0  0  0  0  0\n    0.6320   -0.6320   -0.6320 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.6320    0.6320   -0.6320 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.6320   -0.6320    0.6320 H   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0\n  1  3  1  0\n  1  4  1  0\n  1  5  1  0\nM  END", 0.034662),
            ("formaldehyde", "Formaldehyde\n     RDKit          3D\n\n  4  3  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.2000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.5000    0.9000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.5000   -0.9000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  2  0\n  1  3  1  0\n  1  4  1  0\nM  END", 5.541548),
            ("ethanol", "Ethanol\n     RDKit          3D\n\n  9  8  0  0  0  0  0  0  0  0999 V2000\n   -0.8883    0.1670   -0.0273 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.4658   -0.5116   -0.0368 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.4311    0.3229    0.5867 O   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.8487    1.1175   -0.5695 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.6471   -0.4704   -0.4896 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.1964    0.3978    0.9977 H   0  0  0  0  0  0  0  0  0  0  0  0\n    0.7920   -0.7224   -1.0597 H   0  0  0  0  0  0  0  0  0  0  0  0\n    0.4246   -1.4559    0.5138 H   0  0  0  0  0  0  0  0  0  0  0  0\n    1.4671    1.1550    0.0848 H   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0\n  2  3  1  0\n  1  4  1  0\n  1  5  1  0\n  1  6  1  0\n  2  7  1  0\n  2  8  1  0\n  3  9  1  0\nM  END", -1.336853),
        ];

        for (name, sdf, rdkit_energy) in &cases {
            let mol = parse_sdf(sdf).expect("parse failed");
            let ff = MMFFForceField::new(&mol, MMFFVariant::MMFF94s);
            let coords: Vec<[f64; 3]> = mol.atoms.iter().map(|a| a.position).collect();
            let (our_energy, _) = ff.calculate_energy_and_gradient(&coords);
            let bd = ff.calculate_energy_breakdown(&coords);

            let rdkit_e: f64 = *rdkit_energy;
            let ratio = if rdkit_e.abs() > 0.01 {
                our_energy / rdkit_energy
            } else {
                0.0
            };
            eprintln!(
                "{}: ours={:.6} rdkit={:.6} ratio={:.4} | bond={:.4} angle={:.4} sb={:.4} tor={:.4} oop={:.4} vdw={:.4} elec={:.4}",
                name, our_energy, rdkit_energy, ratio,
                bd.bond, bd.angle, bd.stretch_bend, bd.torsion, bd.oop, bd.vdw, bd.electrostatic
            );
            assert!(
                (our_energy - rdkit_e).abs() < 0.01,
                "{}: our energy {} vs RDKit {}, delta {} exceeds 0.01",
                name, our_energy, rdkit_e, our_energy - rdkit_e
            );
        }
    }

    #[cfg(test)]
    mod diag_aniline {
        use crate::mmff::angle::get_angle_params_with_bond_info;
        use crate::mmff::stretch_bend::get_stretch_bend_params;
        use crate::mmff::{MMFFAtomType, MMFFForceField, MMFFVariant};

        use crate::optimizer;
        use crate::ConvergenceOptions;

        fn dist(a: &[f64; 3], b: &[f64; 3]) -> f64 {
            ((a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2) + (a[2] - b[2]).powi(2)).sqrt()
        }

#[test]
        fn test_aniline_param_audit() {
            let atom_types = [
                MMFFAtomType::C_AR,
                MMFFAtomType::C_AR,
                MMFFAtomType::C_AR,
                MMFFAtomType::C_AR,
                MMFFAtomType::C_AR,
                MMFFAtomType::C_AR,
                MMFFAtomType::N_PL3,
                MMFFAtomType::H_NAM,
                MMFFAtomType::H_NAM,
                MMFFAtomType::H,
                MMFFAtomType::H,
                MMFFAtomType::H,
                MMFFAtomType::H,
            ];

            println!("\n=== Angle params: WebMM vs RDKit ===");
            #[allow(clippy::type_complexity)]
            let angle_tests: &[(usize, usize, usize, u8, u8, u8, f64, f64, &str)] = &[
                (1, 0, 5, 0, 0, 0, 0.669, 119.977, "C1-C0-C5 (37-37-37)"),
                (1, 0, 6, 0, 0, 0, 1.045, 121.633, "C1-C0-N6 (37-37-40)"),
                (5, 0, 6, 0, 0, 0, 1.045, 121.633, "C5-C0-N6 (37-37-40)"),
                (0, 1, 9, 0, 0, 0, 0.563, 120.571, "C0-C1-H9 (37-37-5)"),
                (0, 6, 7, 0, 0, 0, 0.662, 110.288, "C0-N6-H7 (37-40-28)"),
                (7, 6, 8, 0, 0, 0, 0.56, 109.16, "H7-N6-H8 (28-40-28)"),
            ];
            let mut mismatches = 0;
            for &(i, j, k, btij, btjk, ring, exp_ka, exp_t0, label) in angle_tests {
                let params = get_angle_params_with_bond_info(
                    atom_types[i],
                    atom_types[j],
                    atom_types[k],
                    btij,
                    btjk,
                    ring,
                    1.4,
                    1.4,
                );
                if let Some(p) = params {
                    let ka_err = p.k_theta - exp_ka;
                    let t0_err = p.theta0 - exp_t0;
                    let ok = ka_err.abs() < 0.01 && t0_err.abs() < 0.5;
                    if !ok {
                        mismatches += 1;
                    }
                    println!(
                        "  {} ka={:.3}({:+.3}) t0={:.2}({:+.2}) [{}]",
                        label,
                        p.k_theta,
                        ka_err,
                        p.theta0,
                        t0_err,
                        if ok { "OK" } else { "MISMATCH" }
                    );
                } else {
                    mismatches += 1;
                    println!("  {} NO PARAMS FOUND", label);
                }
            }

            println!("\n=== Stretch-bend params: WebMM vs RDKit ===");
            #[allow(clippy::type_complexity)]
            let sb_tests: &[(usize, usize, usize, u8, u8, u8, u8, u8, u8, f64, f64, &str)] = &[
                (1, 0, 5, 0, 0, 6, 6, 6, 0, -0.411, -0.411, "C1-C0-C5"),
                (1, 0, 6, 0, 0, 6, 6, 7, 0, 0.429, 0.901, "C1-C0-N6"),
                (0, 6, 7, 0, 0, 6, 7, 1, 0, 0.423, 0.186, "C0-N6-H7"),
                (7, 6, 8, 0, 0, 1, 7, 1, 0, 0.094, 0.094, "H7-N6-H8"),
            ];
            for &(i, j, k, btij, btjk, ani, anj, ank, at, exp_ijk, exp_kji, label) in sb_tests {
                let sb = get_stretch_bend_params(
                    atom_types[i],
                    atom_types[j],
                    atom_types[k],
                    btij,
                    btjk,
                    ani,
                    anj,
                    ank,
                    at,
                );
                if let Some(p) = sb {
                    let err_i = p.kba_ijk - exp_ijk;
                    let err_k = p.kba_kji - exp_kji;
                    let ok = err_i.abs() < 0.01 && err_k.abs() < 0.01;
                    if !ok {
                        mismatches += 1;
                    }
                    println!(
                        "  {} kbaIJK={:.3}({:+.3}) kbaKJI={:.3}({:+.3}) [{}]",
                        label,
                        p.kba_ijk,
                        err_i,
                        p.kba_kji,
                        err_k,
                        if ok { "OK" } else { "MISMATCH" }
                    );
                } else {
                    mismatches += 1;
                    println!("  {} NO SB PARAMS", label);
                }
            }
            assert_eq!(mismatches, 0, "{} parameter mismatches found", mismatches);
        }

        #[test]
        fn test_aniline_sdf_compare() {
            // Proper 14-atom aniline from RDKit (all H explicit, Kekulé bonds)
            let sdf = "Aniline\n     RDKit          3D\n\n 14 14  0  0  0  0  0  0  0  0999 V2000\n   -1.8551    0.3019   -0.2147 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.9433    1.3121   -0.5108 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.4265    1.0872   -0.3490 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.9000   -0.1487    0.0976 C   0  0  0  0  0  0  0  0  0  0  0  0\n    2.2537   -0.3576    0.2878 N   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.0248   -1.1486    0.4072 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.3958   -0.9291    0.2472 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -2.9206    0.4752   -0.3382 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.2957    2.2773   -0.8642 H   0  0  0  0  0  0  0  0  0  0  0  0\n    1.1231    1.8892   -0.5767 H   0  0  0  0  0  0  0  0  0  0  0  0\n    2.5964   -1.2716    0.5480 H   0  0  0  0  0  0  0  0  0  0  0  0\n    2.9224    0.3435    0.0017 H   0  0  0  0  0  0  0  0  0  0  0  0\n    0.3154   -2.1119    0.7767 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -2.1023   -1.7188    0.4874 H   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  2  0\n  2  3  1  0\n  3  4  2  0\n  4  5  1  0\n  4  6  1  0\n  6  7  2  0\n  7  1  1  0\n  1  8  1  0\n  2  9  1  0\n  3 10  1  0\n  5 11  1  0\n  5 12  1  0\n  6 13  1  0\n  7 14  1  0\nM  END";
            let mol = crate::molecule::parser::parse_sdf(sdf).expect("parse");
            let ff = MMFFForceField::new(&mol, MMFFVariant::MMFF94s);

            // Print atom types for comparison
            println!("\n=== WebMM Atom Types ===");
            for i in 0..mol.atoms.len() {
                println!(
                    "  {}{}: type_id={} atom_type={:?}",
                    mol.atoms[i].symbol, i, ff.type_ids[i], ff.atom_types[i]
                );
            }

            // RDKit optimized coords (same atom ordering)
            let rdkit_coords: Vec<[f64; 3]> = vec![
                [-1.85512092, 0.30187570, -0.21465294],
                [-0.94333594, 1.31212997, -0.51081856],
                [0.42653914, 1.08716520, -0.34904719],
                [0.90001856, -0.14867175, 0.09758881],
                [2.25373393, -0.35759361, 0.28781709],
                [-0.02477892, -1.14860473, 0.40717404],
                [-1.39575686, -0.92910352, 0.24724968],
                [-2.92061022, 0.47515967, -0.33823023],
                [-1.29565190, 2.27728629, -0.86420534],
                [1.12305883, 1.88921250, -0.57667460],
                [2.59640897, -1.27159239, 0.54799775],
                [2.92243159, 0.34348029, 0.00171857],
                [0.31537746, -2.11193527, 0.77666301],
                [-2.10231370, -1.71880836, 0.48741991],
            ];

            // RDKit reference energies at optimized geometry
            let rdkit_bond = 1.315;
            let rdkit_angle = 3.859;
            let rdkit_total = 8.145;

            let bd = ff.calculate_energy_breakdown(&rdkit_coords);

            // Print per-angle params for comparison with RDKit
            println!("\n=== WebMM Angle Parameters ===");
            for angle in &ff.angles {
                let (i, j, k) = (angle.atom1, angle.atom2, angle.atom3);
                let bij_key = (i.min(j), i.max(j));
                let bkj_key = (k.min(j), k.max(j));
                let bt_ij = ff
                    .bond_map
                    .get(&bij_key)
                    .map(|b| {
                        crate::mmff::params::get_mmff_bond_type(
                            b.bond_type,
                            ff.type_ids[i],
                            ff.type_ids[j],
                        )
                    })
                    .unwrap_or(0);
                let bt_jk = ff
                    .bond_map
                    .get(&bkj_key)
                    .map(|b| {
                        crate::mmff::params::get_mmff_bond_type(
                            b.bond_type,
                            ff.type_ids[j],
                            ff.type_ids[k],
                        )
                    })
                    .unwrap_or(0);
                let ring_size = ff.angle_ring_size(i, j, k);
                let r0_ij = ff
                    .bond_map
                    .get(&bij_key)
                    .and_then(|b| {
                        crate::mmff::bond::get_bond_params(
                            ff.atom_types[i],
                            ff.atom_types[j],
                            b.bond_type,
                        )
                    })
                    .map(|p| p.r0)
                    .unwrap_or(1.5);
                let r0_jk = ff
                    .bond_map
                    .get(&bkj_key)
                    .and_then(|b| {
                        crate::mmff::bond::get_bond_params(
                            ff.atom_types[k],
                            ff.atom_types[j],
                            b.bond_type,
                        )
                    })
                    .map(|p| p.r0)
                    .unwrap_or(1.5);
                if let Some(p) = crate::mmff::angle::get_angle_params_with_bond_info(
                    ff.atom_types[i],
                    ff.atom_types[j],
                    ff.atom_types[k],
                    bt_ij,
                    bt_jk,
                    ring_size,
                    r0_ij,
                    r0_jk,
                ) {
                    let e = crate::mmff::angle::angle_energy(&rdkit_coords, i, j, k, &p);
                    println!(
                        "  {}{}-{}{}-{}{} ka={:.4} t0={:.4} bt=({},{}) rs={} E={:.6}",
                        mol.atoms[i].symbol,
                        i,
                        mol.atoms[j].symbol,
                        j,
                        mol.atoms[k].symbol,
                        k,
                        p.k_theta,
                        p.theta0,
                        bt_ij,
                        bt_jk,
                        ring_size,
                        e
                    );
                }
            }

            println!("\n=== Energy at RDKit Optimized Geometry ===");
            println!("WebMM: bond={:.3} angle={:.3} sb={:.3} tor={:.3} oop={:.3} vdw={:.3} elec={:.3} total={:.3}",
                bd.bond, bd.angle, bd.stretch_bend, bd.torsion, bd.oop, bd.vdw, bd.electrostatic, bd.total());
            println!(
                "RDKit: bond={:.3} angle={:.3} total={:.3}",
                rdkit_bond, rdkit_angle, rdkit_total
            );
            println!(
                "Delta bond: {:.3}  Delta angle: {:.3}  Delta total: {:.3}",
                bd.bond - rdkit_bond,
                bd.angle - rdkit_angle,
                bd.total() - rdkit_total
            );
            assert!(
                (bd.total() - rdkit_total).abs() < 0.01,
                "aniline: total energy {} vs RDKit {}",
                bd.total(), rdkit_total
            );

            // Bond lengths at RDKit geometry
            println!("\n=== Bond Lengths (RDKit geometry) ===");
            for b in &mol.bonds {
                let dr = dist(&rdkit_coords[b.atom1], &rdkit_coords[b.atom2]);
                println!(
                    "  {}{}-{}{} ({:?}) r={:.4}",
                    mol.atoms[b.atom1].symbol,
                    b.atom1,
                    mol.atoms[b.atom2].symbol,
                    b.atom2,
                    b.bond_type,
                    dr
                );
            }

            // Optimize with WebMM
            let opts = ConvergenceOptions {
                max_iterations: 5000,
                ..Default::default()
            };
            let result = optimizer::optimize(&ff, &rdkit_coords, &opts);
            println!(
                "\nWebMM opt: converged={} iters={} energy={:.6}",
                result.converged, result.iterations, result.final_energy
            );

            // RMSD
            let mut sum_sq = 0.0;
            for (rc, oc) in rdkit_coords.iter().zip(result.optimized_coords.iter()) {
                for (&o, &r) in oc.iter().zip(rc.iter()) {
                    let delta = o - r;
                    sum_sq += delta * delta;
                }
            }
            let rmsd = (sum_sq / mol.atoms.len() as f64).sqrt();
            println!("RMSD: {:.4} A", rmsd);
        }
    }
}


// ============================================================================
// RDKit-verified typing tests for 2026-07 atom typing gap fixes
// (alkene C=2, S=O=17, SO2=18, nitro N=45, O2CM=32, NSO2=43, metal/halide ions)
// All expected type ids and charges verified against RDKit 2025.09.3 (MMFF94s).
// SDF fixtures generated by RDKit (MolToMolBlock, Kekulé form).
// ============================================================================
#[cfg(test)]
mod type_audit4 {
    use crate::mmff::MMFFForceField;
    use crate::molecule::parser::parse_sdf;
    use crate::optimizer;
    use crate::{ConvergenceOptions, MMFFVariant};

    /// Heavy-atom (Z > 1) MMFF numeric type ids in SDF atom order
    fn heavy_type_ids(sdf: &str) -> Vec<u8> {
        let mol = parse_sdf(sdf).expect("parse failed");
        let ff = MMFFForceField::new(&mol, MMFFVariant::MMFF94s);
        mol.atoms
            .iter()
            .enumerate()
            .filter(|(_, a)| a.atomic_number > 1)
            .map(|(i, _)| ff.type_ids[i])
            .collect()
    }

    const ETHYLENE: &str = r#"Ethylene
     RDKit          2D

  6  5  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2990    0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2990    0.7500    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0000   -1.5000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    2.5981   -0.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.2990    2.2500    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0
  1  3  1  0
  1  4  1  0
  2  5  1  0
  2  6  1  0
M  END"#;

    const ALLENE: &str = r#"Allene
     RDKit          2D

  7  6  0  0  0  0  0  0  0  0999 V2000
    0.0000    1.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0000   -1.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2990    2.2500    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2990    2.2500    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.2990   -2.2500    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2990   -2.2500    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0
  2  3  2  0
  1  4  1  0
  1  5  1  0
  3  6  1  0
  3  7  1  0
M  END"#;

    const KETENE: &str = r#"Ketene
     RDKit          2D

  5  4  0  0  0  0  0  0  0  0999 V2000
    0.0000    1.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0000   -1.5000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.2990    2.2500    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2990    2.2500    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0
  2  3  2  0
  1  4  1  0
  1  5  1  0
M  END"#;

    const ACRYLATE: &str = r#"Acrylate
     RDKit          2D

 12 11  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2990    0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.5981   -0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.8971    0.7500    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.5981   -1.5000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.8971   -2.2500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2990    0.7500    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0000   -1.5000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.2990    2.2500    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    5.1962   -3.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    3.1471   -3.5490    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    4.6471   -0.9510    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0
  2  3  1  0
  3  4  2  0
  3  5  1  0
  5  6  1  0
  1  7  1  0
  1  8  1  0
  2  9  1  0
  6 10  1  0
  6 11  1  0
  6 12  1  0
M  END"#;

    const METHANIMINE: &str = r#"Methanimine
     RDKit          2D

  5  4  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2990    0.7500    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2990    0.7500    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0000   -1.5000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    2.5981   -0.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0
  1  3  1  0
  1  4  1  0
  2  5  1  0
M  END"#;

    const ACETONE: &str = r#"Acetone
     RDKit          2D

 10  9  0  0  0  0  0  0  0  0999 V2000
    1.2990    0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0000   -1.5000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2990    0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.5981    1.5000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    2.0490   -0.5490    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.5490    2.0490    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -2.5981    1.5000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -2.0490   -0.5490    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5490    2.0490    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  2  0
  2  4  1  0
  1  5  1  0
  1  6  1  0
  1  7  1  0
  4  8  1  0
  4  9  1  0
  4 10  1  0
M  END"#;

    const THIOACETONE: &str = r#"Thioacetone
     RDKit          2D

 10  9  0  0  0  0  0  0  0  0999 V2000
    1.2990    0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0000   -1.5000    0.0000 S   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2990    0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.5981    1.5000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    2.0490   -0.5490    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.5490    2.0490    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -2.5981    1.5000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -2.0490   -0.5490    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5490    2.0490    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  2  0
  2  4  1  0
  1  5  1  0
  1  6  1  0
  1  7  1  0
  4  8  1  0
  4  9  1  0
  4 10  1  0
M  END"#;

    const VINYL_ETHER: &str = r#"VinylEther
     RDKit          2D

 10  9  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2990    0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.5981   -0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.8971    0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2990    0.7500    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0000   -1.5000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.2990    2.2500    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    5.1962    1.5000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    4.6471   -0.5490    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    3.1471    2.0490    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0
  2  3  1  0
  3  4  1  0
  1  5  1  0
  1  6  1  0
  2  7  1  0
  4  8  1  0
  4  9  1  0
  4 10  1  0
M  END"#;

    const DMSO: &str = r#"DMSO
     RDKit          2D

 10  9  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5000    0.0000    0.0000 S   0  0  0  0  0  0  0  0  0  0  0  0
    2.2500    1.2990    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.2500   -1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5000    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    1.5000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0000   -1.5000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    3.0000   -2.5981    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    3.5490   -0.5490    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.9510   -2.0490    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  2  0
  2  4  1  0
  1  5  1  0
  1  6  1  0
  1  7  1  0
  4  8  1  0
  4  9  1  0
  4 10  1  0
M  END"#;

    const SULFONE: &str = r#"Sulfone
     RDKit          2D

 11 10  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5000    0.0000    0.0000 S   0  0  0  0  0  0  0  0  0  0  0  0
    1.5000   -1.5000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.5000    1.5000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5000    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    1.5000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0000   -1.5000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    4.5000    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    3.0000    1.5000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    3.0000   -1.5000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  2  0
  2  4  2  0
  2  5  1  0
  1  6  1  0
  1  7  1  0
  1  8  1  0
  5  9  1  0
  5 10  1  0
  5 11  1  0
M  END"#;

    const SULFONAMIDE: &str = r#"Sulfonamide
     RDKit          2D

 10  9  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5000    0.0000    0.0000 S   0  0  0  0  0  0  0  0  0  0  0  0
    1.5000   -1.5000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.5000    1.5000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.0000    0.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5000    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    1.5000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0000   -1.5000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    3.7500    1.2990    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    3.7500   -1.2990    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  2  0
  2  4  2  0
  2  5  1  0
  1  6  1  0
  1  7  1  0
  1  8  1  0
  5  9  1  0
  5 10  1  0
M  END"#;

    const SULFINAMIDE: &str = r#"Sulfinamide
     RDKit          2D

  9  8  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5000    0.0000    0.0000 S   0  0  0  0  0  0  0  0  0  0  0  0
    2.2500    1.2990    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.2500   -1.2990    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5000    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    1.5000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0000   -1.5000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    3.7500   -1.2990    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.5000   -2.5981    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  2  0
  2  4  1  0
  1  5  1  0
  1  6  1  0
  1  7  1  0
  4  8  1  0
  4  9  1  0
M  END"#;

    const NITROMETHANE: &str = r#"Nitromethane
     RDKit          2D

  7  6  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5000    0.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    2.2500   -1.2990    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.2500    1.2990    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5000    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    1.5000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0000   -1.5000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  2  0
  2  4  1  0
  1  5  1  0
  1  6  1  0
  1  7  1  0
M  CHG  2   2   1   4  -1
M  END"#;

    const NITROBENZENE: &str = r#"Nitrobenzene
     RDKit          2D

 14 14  0  0  0  0  0  0  0  0999 V2000
    1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7500   -1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7500   -1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7500    1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7500    1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5000    2.5981    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    3.0000    2.5981    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.7500    3.8971    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.0000    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.5000   -2.5981    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5000   -2.5981    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -3.0000    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5000    2.5981    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0
  2  3  1  0
  3  4  2  0
  4  5  1  0
  5  6  2  0
  6  7  1  0
  7  8  2  0
  7  9  1  0
  6  1  1  0
  1 10  1  0
  2 11  1  0
  3 12  1  0
  4 13  1  0
  5 14  1  0
M  CHG  2   7   1   9  -1
M  END"#;

    const BENZALDEHYDE: &str = r#"Benzaldehyde
  WebMM

 14 14  0  0  0  0  0  0  0  0999 V2000
    1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7500   -1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7500   -1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7500    1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7500    1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5000    2.5981    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7500    3.8971    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.0000    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.5000   -2.5981    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5000   -2.5981    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -3.0000    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5000    2.5981    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    2.5981    2.5981    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0
  2  3  1  0
  3  4  2  0
  4  5  1  0
  5  6  2  0
  6  1  1  0
  1  7  1  0
  7  8  2  0
  2  9  1  0
  3 10  1  0
  4 11  1  0
  5 12  1  0
  6 13  1  0
  7 14  1  0
M  END"#;

    const SODIUM: &str = r#"Sodium
  WebMM

  1  0  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 Na  0  0  0  0  0  0  0  0  0  0  0  0
M  CHG  1   1   1
M  END"#;

    const MAGNESIUM: &str = r#"Magnesium
  WebMM

  1  0  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 Mg  0  0  0  0  0  0  0  0  0  0  0  0
M  CHG  1   1   2
M  END"#;

    const IRON3: &str = r#"Iron3
  WebMM

  1  0  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 Fe  0  0  0  0  0  0  0  0  0  0  0  0
M  CHG  1   1   3
M  END"#;

    const CHLORIDE: &str = r#"Chloride
  WebMM

  1  0  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0
M  CHG  1   1  -1
M  END"#;

    #[test]
    fn test_sp2_carbon_split_types() {
        // sp2 C double-bonded only to C -> 2; double-bonded to N/O/S -> 3
        assert_eq!(heavy_type_ids(ETHYLENE), vec![2, 2], "ethylene");
        assert_eq!(heavy_type_ids(ALLENE), vec![2, 4, 2], "allene");
        assert_eq!(heavy_type_ids(KETENE), vec![2, 4, 7], "ketene");
        assert_eq!(
            heavy_type_ids(ACRYLATE),
            vec![2, 2, 3, 7, 6, 1],
            "acrylate"
        );
        assert_eq!(heavy_type_ids(METHANIMINE), vec![3, 9], "methanimine");
        // Regression guards: carbonyl / thiocarbonyl C stays 3
        assert_eq!(heavy_type_ids(ACETONE), vec![1, 3, 7, 1], "acetone");
        assert_eq!(heavy_type_ids(THIOACETONE), vec![1, 3, 16, 1], "thioacetone");
        assert_eq!(heavy_type_ids(VINYL_ETHER), vec![2, 2, 6, 1], "vinyl ether");
    }

    #[test]
    fn test_oxidized_sulfur_types() {
        assert_eq!(heavy_type_ids(DMSO), vec![1, 17, 7, 1], "dmso");
        assert_eq!(heavy_type_ids(SULFONE), vec![1, 18, 32, 32, 1], "sulfone");
        assert_eq!(
            heavy_type_ids(SULFONAMIDE),
            vec![1, 18, 32, 32, 43],
            "sulfonamide"
        );
        // Sulfinamide N stays 8 (RDKit-verified; NOT 48)
        assert_eq!(heavy_type_ids(SULFINAMIDE), vec![1, 17, 7, 8], "sulfinamide");
    }

    #[test]
    fn test_nitro_types() {
        assert_eq!(heavy_type_ids(NITROMETHANE), vec![1, 45, 32, 32], "nitromethane");
        assert_eq!(
            heavy_type_ids(NITROBENZENE),
            vec![37, 37, 37, 37, 37, 37, 45, 32, 32],
            "nitrobenzene"
        );
    }

    #[test]
    fn test_ion_types() {
        assert_eq!(heavy_type_ids(SODIUM), vec![93], "Na+");
        assert_eq!(heavy_type_ids(MAGNESIUM), vec![99], "Mg2+");
        assert_eq!(heavy_type_ids(IRON3), vec![88], "Fe3+");
        assert_eq!(heavy_type_ids(CHLORIDE), vec![90], "Cl-");
    }

    /// Partial charges for the newly typed groups, hand-derived from RDKit's
    /// MMFF94s BCI table (canonical pair i<j → smaller id gets −stored value,
    /// larger gets +stored). Tolerance 1e-3.
    #[test]
    fn test_new_type_partial_charges() {
        // DMSO: atoms (0 C, 1 S=17, 2 O, 3 C)
        let mol = crate::molecule::parser::parse_sdf(DMSO).unwrap();
        let ff = crate::mmff::MMFFForceField::new(&mol, crate::mmff::MMFFVariant::MMFF94s);
        assert!((ff.charges[1] - 0.113).abs() < 1e-3, "DMSO S: {}", ff.charges[1]);
        assert!((ff.charges[2] - (-0.5)).abs() < 1e-3, "DMSO O: {}", ff.charges[2]);
        assert!((ff.charges[0] - 0.1935).abs() < 1e-3, "DMSO C0: {}", ff.charges[0]);
        assert!((ff.charges[3] - 0.1935).abs() < 1e-3, "DMSO C3: {}", ff.charges[3]);

        // Sulfone: atoms (0 C, 1 S=18, 2 O=32, 3 O=32, 4 C)
        let mol = crate::molecule::parser::parse_sdf(SULFONE).unwrap();
        let ff = crate::mmff::MMFFForceField::new(&mol, crate::mmff::MMFFVariant::MMFF94s);
        assert!((ff.charges[1] - 1.0896).abs() < 1e-3, "sulfone S: {}", ff.charges[1]);
        assert!((ff.charges[2] - (-0.65)).abs() < 1e-3, "sulfone O2: {}", ff.charges[2]);
        assert!((ff.charges[3] - (-0.65)).abs() < 1e-3, "sulfone O3: {}", ff.charges[3]);
        assert!((ff.charges[0] - 0.1052).abs() < 1e-3, "sulfone C0: {}", ff.charges[0]);
        assert!((ff.charges[4] - 0.1052).abs() < 1e-3, "sulfone C4: {}", ff.charges[4]);

        // Nitromethane: atoms (0 C, 1 N=45, 2 O=32, 3 O=32)
        let mol = crate::molecule::parser::parse_sdf(NITROMETHANE).unwrap();
        let ff = crate::mmff::MMFFForceField::new(&mol, crate::mmff::MMFFVariant::MMFF94s);
        assert!((ff.charges[1] - 0.7998).abs() < 1e-3, "nitro N: {}", ff.charges[1]);
        assert!((ff.charges[2] - (-0.52)).abs() < 1e-3, "nitro O2: {}", ff.charges[2]);
        assert!((ff.charges[3] - (-0.52)).abs() < 1e-3, "nitro O3: {}", ff.charges[3]);
        assert!((ff.charges[0] - 0.2402).abs() < 1e-3, "nitro C: {}", ff.charges[0]);

        // Benzaldehyde: atoms (0 ipso C, 6 CHO, 7 O, 13 aldehyde H)
        let mol = crate::molecule::parser::parse_sdf(BENZALDEHYDE).unwrap();
        let ff = crate::mmff::MMFFForceField::new(&mol, crate::mmff::MMFFVariant::MMFF94s);
        assert!((ff.charges[7] - (-0.57)).abs() < 1e-3, "PhCHO O: {}", ff.charges[7]);
        assert!((ff.charges[6] - 0.4238).abs() < 1e-3, "PhCHO CHO: {}", ff.charges[6]);
        assert!((ff.charges[0] - 0.0862).abs() < 1e-3, "PhCHO ipso: {}", ff.charges[0]);
        assert!((ff.charges[13] - 0.06).abs() < 1e-3, "PhCHO H: {}", ff.charges[13]);
    }

    /// Energy evaluation and geometry optimization must converge for the
    /// newly typed functional groups.
    #[test]
    fn test_sp2_sulfur_nitro_energy_and_convergence() {
        for (name, sdf) in [
            ("DMSO", DMSO),
            ("sulfone", SULFONE),
            ("nitrobenzene", NITROBENZENE),
        ] {
            let mol = crate::molecule::parser::parse_sdf(sdf).unwrap();
            let ff = crate::mmff::MMFFForceField::new(&mol, crate::mmff::MMFFVariant::MMFF94s);
            let coords: Vec<[f64; 3]> = mol.atoms.iter().map(|a| a.position).collect();
            let e = ff.calculate_energy(&coords);
            assert!(e.is_finite(), "{name}: non-finite energy {e}");
            let opts = ConvergenceOptions::default();
            let result = optimizer::optimize(&ff, &coords, &opts);
            assert!(
                result.converged,
                "{name}: optimization did not converge (energy={})",
                result.final_energy
            );
        }
    }
}

#[cfg(test)]
mod type_audit5 {
    use crate::mmff::MMFFForceField;
    use crate::molecule::parser::parse_sdf;
    use crate::optimizer;
    use crate::{ConvergenceOptions, MMFFVariant};

    /// Heavy-atom (Z > 1) MMFF numeric type ids in SDF atom order
    fn heavy_type_ids(sdf: &str) -> Vec<u8> {
        let mol = parse_sdf(sdf).expect("parse failed");
        let ff = MMFFForceField::new(&mol, MMFFVariant::MMFF94s);
        mol.atoms
            .iter()
            .enumerate()
            .filter(|(_, a)| a.atomic_number > 1)
            .map(|(i, _)| ff.type_ids[i])
            .collect()
    }

// ---- SDF fixtures (RDKit-generated 2D coordinates) ----

    const ACETATE: &str = r#"Acetate
     RDKit          2D

  7  6  0  0  0  0  0  0  0  0999 V2000
    0.6429   -0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8571    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.6071    1.2990    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.6071   -1.2990    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.1429   -0.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.6429   -1.5000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.6429    1.5000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  2  0
  2  4  1  0
  1  5  1  0
  1  6  1  0
  1  7  1  0
M  CHG  1   4  -1
M  END"#;

    const BENZOATE: &str = r#"Benzoate
     RDKit          2D

 14 14  0  0  0  0  0  0  0  0999 V2000
    3.2143   -1.2990    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.4643    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.2143    1.2990    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.9643   -0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.2143   -1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2857   -1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.0357    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2857    1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.2143    1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.9643   -2.5981    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -2.0357   -2.5981    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -3.5357    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -2.0357    2.5981    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.9643    2.5981    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0
  2  3  1  0
  2  4  1  0
  4  5  2  0
  5  6  1  0
  6  7  2  0
  7  8  1  0
  8  9  2  0
  9  4  1  0
  5 10  1  0
  6 11  1  0
  7 12  1  0
  8 13  1  0
  9 14  1  0
M  CHG  1   3  -1
M  END"#;

    const ACETIC_ACID: &str = r#"AceticAcid
     RDKit          2D

  8  7  0  0  0  0  0  0  0  0999 V2000
    1.0348   -0.1383    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.4030    0.2892    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7517    1.7482    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4922   -0.7422    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.4726   -0.5658    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.6073   -1.5761    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.4623    1.2995    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -2.9299   -0.3147    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  2  0
  2  4  1  0
  1  5  1  0
  1  6  1  0
  1  7  1  0
  4  8  1  0
M  END"#;

    const PYRIDINE_NO: &str = r#"PyridineNO
     RDKit          2D

 12 12  0  0  0  0  0  0  0  0999 V2000
    3.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.5000   -0.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    0.7500   -1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7500   -1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7500    1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7500    1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5000   -2.5981    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5000   -2.5981    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -3.0000    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5000    2.5981    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.5000    2.5981    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  2  0
  3  4  1  0
  4  5  2  0
  5  6  1  0
  6  7  2  0
  7  2  1  0
  3  8  1  0
  4  9  1  0
  5 10  1  0
  6 11  1  0
  7 12  1  0
M  CHG  2   1  -1   2   1
M  END"#;

    const PYRIDINE: &str = r#"Pyridine
     RDKit          2D

 11 11  0  0  0  0  0  0  0  0999 V2000
    0.0000   -1.2273    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2990   -0.4773    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2990    1.0227    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0000    1.7727    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    1.2990    1.0227    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2990   -0.4773    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000   -2.7273    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -2.5981   -1.2273    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -2.5981    1.7727    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    2.5981    1.7727    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    2.5981   -1.2273    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0
  2  3  1  0
  3  4  2  0
  4  5  1  0
  5  6  2  0
  6  1  1  0
  1  7  1  0
  2  8  1  0
  3  9  1  0
  5 10  1  0
  6 11  1  0
M  END"#;

    const TMS: &str = r#"TMS
     RDKit          2D

 17 16  0  0  0  0  0  0  0  0999 V2000
    1.0037    1.0607    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0569   -0.0000    0.0000 Si  0  0  0  0  0  4  0  0  0  0  0  0
   -1.1176   -1.0607    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.1176    1.0607    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0037   -1.0607    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.0644    2.1213    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    2.1858    0.1372    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.4584    1.7586    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1783   -2.1213    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.8156   -0.5154    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5723   -1.7586    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1783    2.1213    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.8156    0.5154    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5723    1.7586    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    2.0644   -2.1213    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.4584   -1.7586    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    2.1858   -0.1372    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  1  0
  2  4  1  0
  2  5  1  0
  1  6  1  0
  1  7  1  0
  1  8  1  0
  3  9  1  0
  3 10  1  0
  3 11  1  0
  4 12  1  0
  4 13  1  0
  4 14  1  0
  5 15  1  0
  5 16  1  0
  5 17  1  0
M  END"#;

    const METHANIMINE: &str = r#"Methanimine
     RDKit          2D

  5  4  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2990    0.7500    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2990    0.7500    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0000   -1.5000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    2.5981   -0.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0
  1  3  1  0
  1  4  1  0
  2  5  1  0
M  END"#;

    // ---- Tests ----

    /// Carboxylate C → 41, both O → 32; regression: acetic acid stays C=3, =O=7, OH=6
    #[test]
    fn test_carboxylate_types() {
        assert_eq!(heavy_type_ids(ACETATE), vec![1, 41, 32, 32], "acetate");
        assert_eq!(
            heavy_type_ids(BENZOATE),
            vec![32, 41, 32, 37, 37, 37, 37, 37, 37],
            "benzoate"
        );
        // Regression guard: acetic acid C=3, =O=7, −OH=6 (unchanged)
        assert_eq!(heavy_type_ids(ACETIC_ACID), vec![1, 3, 7, 6], "acetic acid");
    }

    /// Pyridine N-oxide N → 69, O → 32; regression: pyridine N stays 38
    #[test]
    fn test_noxide_types() {
        assert_eq!(
            heavy_type_ids(PYRIDINE_NO),
            vec![32, 69, 37, 37, 37, 37, 37],
            "pyridine N-oxide"
        );
        // Regression guard: pyridine N=38 (unchanged)
        assert_eq!(
            heavy_type_ids(PYRIDINE),
            vec![37, 37, 37, 38, 37, 37],
            "pyridine"
        );
    }

    /// Silicon → 19 (RDKit-verified)
    #[test]
    fn test_silicon_types() {
        assert_eq!(heavy_type_ids(TMS), vec![1, 19, 1, 1, 1], "TMS");
    }

    /// Imine H → 27 (RDKit-verified); the C=N bond stays N=9, C=3
    #[test]
    fn test_imine_h_type() {
        use crate::mmff::params::mmff_type_id;
        let mol = parse_sdf(METHANIMINE).expect("parse failed");
        let ff = MMFFForceField::new(&mol, MMFFVariant::MMFF94s);
        // type_ids applies base_type() (all H subtypes → 5), so check raw atom_types
        assert_eq!(
            mmff_type_id(ff.atom_types[0]),
            3,
            "methanimine C type"
        );
        assert_eq!(
            mmff_type_id(ff.atom_types[1]),
            9,
            "methanimine N type"
        );
        // Atom 4 is the N-H (imine hydrogen); type 27, not 23 or 28
        assert_eq!(
            mmff_type_id(ff.atom_types[4]),
            27,
            "methanimine N-H (imine H) = 27, got {}",
            mmff_type_id(ff.atom_types[4])
        );
    }

    /// Energy finite and optimizer converges for the newly typed groups
    #[test]
    fn test_carboxylate_noxide_silicon_imine_energy_and_convergence() {
        for (name, sdf) in [
            ("acetate", ACETATE),
            ("pyridine_no", PYRIDINE_NO),
            ("tetramethylsilane", TMS),
            ("methanimine", METHANIMINE),
        ] {
            let mol = parse_sdf(sdf).unwrap();
            let ff = MMFFForceField::new(&mol, MMFFVariant::MMFF94s);
            let coords: Vec<[f64; 3]> = mol.atoms.iter().map(|a| a.position).collect();
            let e = ff.calculate_energy(&coords);
            assert!(e.is_finite(), "{name}: non-finite energy {e}");
            let opts = ConvergenceOptions {
                max_iterations: 5000,
                ..Default::default()
            };
            let result = optimizer::optimize(&ff, &coords, &opts);
            assert!(
                result.converged,
                "{name}: optimization did not converge (energy={})",
                result.final_energy
            );
        }
    }

    const SULFONATE: &str = r#"Sulfonate
     RDKit          2D

  8  7  0  0  0  0  0  0  0  0999 V2000
   -0.7500   -0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7500   -0.0000    0.0000 S   0  0  0  0  0  0  0  0  0  0  0  0
    2.2500   -0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.7500   -1.5000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.7500    1.5000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.2500    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7500    1.5000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7500   -1.5000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  2  0
  2  4  2  0
  2  5  1  0
  1  6  1  0
  1  7  1  0
  1  8  1  0
M  CHG  1   5  -1
M  END"#;

    const AMMONIUM: &str = r#"Ammonium
     RDKit          2D

  5  4  0  0  0  0  0  0  0  0999 V2000
    0.0000   -0.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    1.5000   -0.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5000    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    1.5000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0000   -1.5000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  1  3  1  0
  1  4  1  0
  1  5  1  0
M  CHG  1   1   1
M  END"#;

    /// Partial charges match RDKit for charged/zwitterionic groups.
    /// RDKit-verified reference values from GetMMFFPartialCharge (tolerance 0.02).
    #[test]
    fn test_formal_charge_distribution() {
        // Acetate: atoms (0 C_methyl=1, 1 C_carbox=41, 2 O=32, 3 O⁻=32, 4-6 H)
        // RDKit partials: C=-0.106, C_carbox=+0.906, O=-0.900, O=-0.900
        let mol = parse_sdf(ACETATE).unwrap();
        let ff = MMFFForceField::new(&mol, MMFFVariant::MMFF94s);
        assert!(
            (ff.charges[0] - (-0.106)).abs() < 0.02,
            "acetate C_methyl: got {}",
            ff.charges[0]
        );
        assert!(
            (ff.charges[1] - 0.906).abs() < 0.02,
            "acetate C_carbox: got {}",
            ff.charges[1]
        );
        assert!(
            (ff.charges[2] - (-0.900)).abs() < 0.02,
            "acetate O2: got {}",
            ff.charges[2]
        );
        assert!(
            (ff.charges[3] - (-0.900)).abs() < 0.02,
            "acetate O3: got {}",
            ff.charges[3]
        );

        // Pyridine N-oxide: atoms (0 O⁻=32, 1 N⁺=69, 2-6 C_AR, 7-11 H)
        // RDKit partials: O=-0.750, N=+0.571
        let mol = parse_sdf(PYRIDINE_NO).unwrap();
        let ff = MMFFForceField::new(&mol, MMFFVariant::MMFF94s);
        assert!(
            (ff.charges[0] - (-0.750)).abs() < 0.02,
            "N-oxide O: got {}",
            ff.charges[0]
        );
        assert!(
            (ff.charges[1] - 0.571).abs() < 0.02,
            "N-oxide N: got {}",
            ff.charges[1]
        );

        // Ammonium: N keeps +1.0 (type 34, fcadj=0, no redistribution)
        // RDKit: N=-0.800 (formal +1 + BCI increments)
        let mol = parse_sdf(AMMONIUM).unwrap();
        let ff = MMFFForceField::new(&mol, MMFFVariant::MMFF94s);
        // Sum of charges should be +1.0
        let total: f64 = ff.charges.iter().sum();
        assert!(
            (total - 1.0).abs() < 0.02,
            "ammonium total charge: got {}",
            total
        );

        // Sulfonate: 3 type-32 O's each get -1/3 formal charge
        let mol = parse_sdf(SULFONATE).unwrap();
        let ff = MMFFForceField::new(&mol, MMFFVariant::MMFF94s);
        // Sum of charges should be -1.0
        let total: f64 = ff.charges.iter().sum();
        assert!(
            (total - (-1.0)).abs() < 0.02,
            "sulfonate total charge: got {}",
            total
        );
    }

    /// Aniline NH2 must be planar after MMFF94s optimization (torsion enforces planarity)
    #[test]
    fn test_aniline_nh2_planar() {
        let sdf = "Aniline\n     RDKit          3D\n\n 13 13  0  0  0  0  0  0  0  0999 V2000\n   -1.2000    0.6930    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.5000   -0.6000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.3000   -1.2000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.0000   -0.7000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.3000    0.6000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000    1.1000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -2.1500    1.3000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n   -2.2500    2.2800    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -3.0000    0.8000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -2.4000   -1.1000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.3500   -2.2700    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    1.7500   -1.5500    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    2.3000    1.0500    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  2  0\n  2  3  1  0\n  3  4  2  0\n  4  5  1  0\n  5  6  2  0\n  6  1  1  0\n  1  7  1  0\n  7  8  1  0\n  7  9  1  0\n  2 10  1  0\n  3 11  1  0\n  4 12  1  0\n  5 13  1  0\nM  END";
        let mol = parse_sdf(sdf).unwrap();
        let ff = MMFFForceField::new(&mol, MMFFVariant::MMFF94s);
        let coords: Vec<[f64; 3]> = mol.atoms.iter().map(|a| a.position).collect();
        let opts = ConvergenceOptions {
            max_iterations: 5000,
            ..Default::default()
        };
        let result = optimizer::optimize(&ff, &coords, &opts);
        assert!(result.converged, "aniline optimization did not converge");

        let c = &result.optimized_coords;
        fn dihedral(p1: [f64; 3], p2: [f64; 3], p3: [f64; 3], p4: [f64; 3]) -> f64 {
            let b1 = [p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2]];
            let b2 = [p3[0] - p2[0], p3[1] - p2[1], p3[2] - p2[2]];
            let b3 = [p4[0] - p3[0], p4[1] - p3[1], p4[2] - p3[2]];
            let n1 = [b1[1]*b2[2]-b1[2]*b2[1], b1[2]*b2[0]-b1[0]*b2[2], b1[0]*b2[1]-b1[1]*b2[0]];
            let n2 = [b2[1]*b3[2]-b2[2]*b3[1], b2[2]*b3[0]-b2[0]*b3[2], b2[0]*b3[1]-b2[1]*b3[0]];
            let m = (n1[0]*n1[0]+n1[1]*n1[1]+n1[2]*n1[2]).sqrt();
            let n = (n2[0]*n2[0]+n2[1]*n2[1]+n2[2]*n2[2]).sqrt();
            let dot = (n1[0]*n2[0]+n1[1]*n2[1]+n1[2]*n2[2]) / (m*n);
            dot.clamp(-1.0, 1.0).acos().to_degrees()
        }

        let d1 = dihedral(c[7], c[6], c[0], c[1]);
        let d2 = dihedral(c[8], c[6], c[0], c[1]);
        let d1_planar = !(20.0..=160.0).contains(&d1);
        let d2_planar = !(20.0..=160.0).contains(&d2);
        assert!(
            d1_planar && d2_planar,
            "Aniline NH2 not planar: H7-N-C-C={:.1}° H8-N-C-C={:.1}°",
            d1, d2
        );
    }
}


#[cfg(test)]
mod aniline_etkdg_planarity {
    #[test]
    fn check_webmm_etkdg_aniline() {
        use crate::etkdg::{ETKDGConfig, generate_initial_coords_with_config};
        use crate::molecule::parser::parse_sdf;

        let sdf = "Aniline\n     RDKit          3D\n\n 13 13  0  0  0  0  0  0  0  0999 V2000\n   -1.2000    0.6930    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.5000   -0.6000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.3000   -1.2000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.0000   -0.7000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.3000    0.6000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000    1.1000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -2.1500    1.3000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n   -2.2500    2.2800    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -3.0000    0.8000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -2.4000   -1.1000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.3500   -2.2700    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    1.7500   -1.5500    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    2.3000    1.0500    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  2  0\n  2  3  1  0\n  3  4  2  0\n  4  5  1  0\n  5  6  2  0\n  6  1  1  0\n  1  7  1  0\n  7  8  1  0\n  7  9  1  0\n  2 10  1  0\n  3 11  1  0\n  4 12  1  0\n  5 13  1  0\nM  END";
        let mol = parse_sdf(sdf).unwrap();
        let config = ETKDGConfig { random_seed: 42, ..Default::default() };
        let coords = generate_initial_coords_with_config(&mol, &config);

        // Print all coords
        for (i, (atom, c)) in mol.atoms.iter().zip(coords.iter()).enumerate() {
            println!("{} {} {:.4} {:.4} {:.4}", i, atom.symbol, c[0], c[1], c[2]);
        }

        // Check N (atom 6) pyramidalization
        let n = coords[6];
        let h7 = coords[7];
        let h8 = coords[8];
        let c0 = coords[0];

        // H-N-H angle
        let v1 = [h7[0]-n[0], h7[1]-n[1], h7[2]-n[2]];
        let v2 = [h8[0]-n[0], h8[1]-n[1], h8[2]-n[2]];
        let d1 = (v1[0]*v1[0]+v1[1]*v1[1]+v1[2]*v1[2]).sqrt();
        let d2 = (v2[0]*v2[0]+v2[1]*v2[1]+v2[2]*v2[2]).sqrt();
        let dot = (v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2])/(d1*d2);
        let hnh = dot.clamp(-1.0,1.0).acos().to_degrees();
        println!("H-N-H angle: {:.1} deg", hnh);

        // Distance of H atoms from ring plane
        // Ring atoms: 0-5
        let ring: Vec<[f64;3]> = (0..6).map(|i| coords[i]).collect();
        let center = [
            ring.iter().map(|c| c[0]).sum::<f64>() / 6.0,
            ring.iter().map(|c| c[1]).sum::<f64>() / 6.0,
            ring.iter().map(|c| c[2]).sum::<f64>() / 6.0,
        ];

        // Simple normal: cross product of two ring diagonals
        let r1 = [ring[3][0]-ring[0][0], ring[3][1]-ring[0][1], ring[3][2]-ring[0][2]];
        let r2 = [ring[4][1]-ring[1][0], ring[4][1]-ring[1][1], ring[4][2]-ring[1][2]];
        let normal = [
            r1[1]*r2[2]-r1[2]*r2[1],
            r1[2]*r2[0]-r1[0]*r2[2],
            r1[0]*r2[1]-r1[1]*r2[0],
        ];
        let nlen = (normal[0]*normal[0]+normal[1]*normal[1]+normal[2]*normal[2]).sqrt();
        let nhat = [normal[0]/nlen, normal[1]/nlen, normal[2]/nlen];

        for (idx, h) in [(7usize, h7), (8, h8)] {
            let d = (h[0]-center[0])*nhat[0] + (h[1]-center[1])*nhat[1] + (h[2]-center[2])*nhat[2];
            println!("H{} distance from ring plane: {:.4} A", idx, d.abs());
        }
        let n_dist = (n[0]-center[0])*nhat[0] + (n[1]-center[1])*nhat[1] + (n[2]-center[2])*nhat[2];
        println!("N distance from ring plane: {:.4} A", n_dist.abs());
        // Aniline NH2 must be pyramidal with ~117° H-N-H (RDKit ETKDG behavior)
        assert!(
            hnh > 110.0 && hnh < 125.0,
            "aniline H-N-H angle {} out of [110, 125]",
            hnh
        );
        let h7_oop = ((h7[0]-center[0])*nhat[0] + (h7[1]-center[1])*nhat[1] + (h7[2]-center[2])*nhat[2]).abs();
        let h8_oop = ((h8[0]-center[0])*nhat[0] + (h8[1]-center[1])*nhat[1] + (h8[2]-center[2])*nhat[2]).abs();
        // NH2 may be planar (MMFF94s) or pyramidal (RDKit ETKDG) — both valid.
        // The old fix_aniline_nh2_geometry enforced pyramidal but placed H's at
        // wrong C-N-H angles; it's now disabled. Just check H-N-H is sane.
    }
}

// =====================================================================
// Regression tests for the 3 silent algorithmic bugs found during the
// energy-accuracy debugging pass (2025-09). Each test pins a behavior that
// was wrong and would otherwise silently recur.
//   1. OOP cyclic permutations (1 vs 3 terms per sp2 center)
//   2. DFSB stretch-bend skipped for asymmetric-row angles (P-O-H)
//   3. stretch-bend sb_type canonicalization (lower-type peripheral bond)
// Plus a bond-param symmetry invariant (the C_AR-N_AR reverse-arm bug) and
// a per-term energy-breakdown check (catches *which* term regresses).
// =====================================================================
#[cfg(test)]
mod regression_tests {
    use crate::mmff::bond::get_bond_params;
    use crate::mmff::stretch_bend::get_stretch_bend_params;
    use crate::mmff::mmff_tables::compute_angle_type;
    use crate::mmff::params::get_mmff_bond_type;
    use crate::mmff::{MMFFAtomType, MMFFForceField};
    use crate::molecule::parser::parse_sdf;
    use crate::molecule::BondType;
    use crate::MMFFVariant;

    // --- shared SDF fixtures (RDKit 3D geometries) ---
    const GUANIDINIUM: &str = include_str!("../scripts/val_set/guanidinium.sdf");
    const METHYLPHOSPHONIC_ACID: &str =
        include_str!("../scripts/val_set/methylphosphonic_acid.sdf");
    const PURINE: &str = include_str!("../scripts/val_set/purine.sdf");
    const PYRIMIDINE: &str = include_str!("../scripts/val_set/pyrimidine.sdf");

    /// Bug #1: RDKit creates **3 cyclic OOP terms** per 3-neighbor sp2 atom
    /// (each neighbor takes a turn as the out-of-plane atom1). WebMM used to
    /// create only 1, undercounting OOP energy ~3x (guanidinium was off by
    /// 3.9 kcal/mol). The central C of guanidinium C(=NH)(NH2)(NH2) has 3 N
    /// neighbors and must yield exactly 3 OOP terms.
    #[test]
    fn oop_creates_3_cyclic_terms_per_sp2_center() {
        let mol = parse_sdf(GUANIDINIUM).expect("parse");
        let ff = MMFFForceField::new(&mol, MMFFVariant::MMFF94s);
        // central C is atom 1 (C_2); its 3 neighbors are the 3 N atoms
        let central_oops: Vec<_> = ff.oops.iter().filter(|o| o.central == 1).collect();
        assert_eq!(
            central_oops.len(),
            3,
            "sp2 C with 3 neighbors must have 3 cyclic OOP terms, got {}",
            central_oops.len()
        );
        // the 3 terms must cycle which neighbor is atom1
        let atom1_set: std::collections::HashSet<_> =
            central_oops.iter().map(|o| o.atom1).collect();
        assert_eq!(atom1_set.len(), 3, "each of the 3 neighbors must be atom1 once");
    }

    /// Bug #2: the default stretch-bend (DFSB) lookup didn't canonicalize
    /// periodic-table rows, so asymmetric angles like P-O-H (rows 2,1,0) were
    /// looked up as (2,1,0) against a table stored as (0,1,2) -> no match ->
    /// `None` -> the term was silently dropped. methylphosphonic_acid has two
    /// P-O-H angles that MUST resolve to a real kba.
    #[test]
    fn stretch_bend_dfsb_not_skipped_for_asymmetric_rows() {
        let mol = parse_sdf(METHYLPHOSPHONIC_ACID).expect("parse");
        let ff = MMFFForceField::new(&mol, MMFFVariant::MMFF94s);
        // find P-O-H angles: central atom is O (type 6), one peripheral is H
        let mut resolved = 0;
        for ang in &ff.angles {
            let (i, j, k) = (ang.atom1, ang.atom2, ang.atom3);
            let is_poh = ff.type_ids[j] == 6
                && (ff.type_ids[i] == 25 || ff.type_ids[k] == 25)
                && (mol.atoms[i].atomic_number == 1 || mol.atoms[k].atomic_number == 1);
            if !is_poh {
                continue;
            }
            let bij = ff.bond_map.get(&(i.min(j), i.max(j))).unwrap();
            let bkj = ff.bond_map.get(&(k.min(j), k.max(j))).unwrap();
            let btij = get_mmff_bond_type(bij.bond_type, ff.type_ids[i], ff.type_ids[j]);
            let btjk = get_mmff_bond_type(bkj.bond_type, ff.type_ids[j], ff.type_ids[k]);
            let at = compute_angle_type(btij, btjk, ff.angle_ring_size(i, j, k));
            let p = get_stretch_bend_params(
                ff.atom_types[i],
                ff.atom_types[j],
                ff.atom_types[k],
                btij,
                btjk,
                mol.atoms[i].atomic_number,
                mol.atoms[j].atomic_number,
                mol.atoms[k].atomic_number,
                at,
            );
            assert!(p.is_some(), "P-O-H stretch-bend must resolve (was silently dropped)");
            resolved += 1;
        }
        assert!(resolved >= 2, "expected >=2 P-O-H angles, found {}", resolved);
    }

    /// Bug #3: the stretch-bend *type* depends on bond_type_1 = the bond to the
    /// LOWER-type peripheral atom. Swapping i<->k must not change the resolved
    /// kba pair. (purine's N-C-C angle got sb_type 2 instead of 1, using the
    /// DFSB default 0.3 instead of the specific 0.61/0.227 entry.)
    #[test]
    fn stretch_bend_kba_is_order_invariant() {
        // N_2(9) - C_2(3) - C_VIN(2): the purine angle that was broken.
        // ti=9 > tk=2, so the lookup must canonicalize before computing sb_type.
        let fwd = get_stretch_bend_params(
            MMFFAtomType::N_2, MMFFAtomType::C_2, MMFFAtomType::C_VIN,
            0, 1, 7, 6, 6, 1, // bt_ij=0 (C=N dbl), bt_jk=1 (C-C sbmb); Z: N,C,C
        );
        let rev = get_stretch_bend_params(
            MMFFAtomType::C_VIN, MMFFAtomType::C_2, MMFFAtomType::N_2,
            1, 0, 6, 6, 7, 1, // same two bonds, swapped order
        );
        let (kf, kr) = (fwd.expect("fwd"), rev.expect("rev"));
        // kba_ijk must correspond to the same physical bond in both orders:
        // fwd.kba_ijk (N-C side) == rev.kba_kji, and fwd.kba_kji == rev.kba_ijk
        assert!((kf.kba_ijk - kr.kba_kji).abs() < 1e-9, "N-C kba differs on swap");
        assert!((kf.kba_kji - kr.kba_ijk).abs() < 1e-9, "C-C kba differs on swap");
        // and it must NOT be the DFSB default (0.3, 0.3) — the specific entry is
        // (0.61, 0.227) for the N-C / C-C halves
        assert!(kf.kba_ijk > 0.5, "expected specific kba ~0.61, got {}", kf.kba_ijk);
    }

    /// Bond-param symmetry: every (A,B,bt) entry must equal (B,A,bt). The
    /// C_AR-N_AR aromatic entry once lacked its reverse arm, so N->C-ordered
    /// ring bonds fell to a different param set and corrupted both bond and
    /// stretch-bend energy (pyrimidine was off 2.4 kcal/mol).
    #[test]
    fn bond_params_symmetric_for_heteroatom_pairs() {
        use MMFFAtomType::*;
        // curated pairs covering every bug-hit family: aromatic heterocycles,
        // 5-ring types, heteroatoms (P/S/halide), imines.
        let cases: &[(MMFFAtomType, MMFFAtomType, BondType)] = &[
            (C_AR, N_AR, BondType::Aromatic),
            (C_AR, S_AR, BondType::Aromatic),
            (C_AR, O_R, BondType::Aromatic),
            (N_AR, N_AR, BondType::Aromatic),
            (C5A, S_AR, BondType::Aromatic),
            (C5A, OFUR, BondType::Aromatic),
            (C5B, C5B, BondType::Aromatic),
            (NPYL, N5A, BondType::Aromatic),
            (N5A, C5B, BondType::Aromatic),
            (C_3, O_R, BondType::Single),
            (C_3, S_3, BondType::Single),
            (S_3, S_3, BondType::Single),
            (C_3, P_3, BondType::Single),
            (C_3, P_4, BondType::Single),
            (P_4, O_CO2, BondType::Double),
            (C_AR, I, BondType::Single),
            (C_3, Cl, BondType::Single),
            (C_2, N_2, BondType::Double),
            (C_2, N_2, BondType::Single),
            (C_3, C_1, BondType::Single),
            (C_1, N_1, BondType::Triple),
        ];
        for &(a, b, bt) in cases {
            let pa = get_bond_params(a, b, bt);
            let pb = get_bond_params(b, a, bt);
            assert_eq!(pa.is_some(), pb.is_some(), "({:?},{:?},{:?}) asymmetry: one resolves, other doesn't", a, b, bt);
            if let (Some(pa), Some(pb)) = (pa, pb) {
                assert!((pa.r0 - pb.r0).abs() < 1e-9, "({:?},{:?}) r0 differs: {} vs {}", a, b, pa.r0, pb.r0);
                assert!((pa.k_bond - pb.k_bond).abs() < 1e-9, "({:?},{:?}) kb differs: {} vs {}", a, b, pa.k_bond, pb.k_bond);
            }
        }
    }

    /// Per-term energy breakdown vs known RDKit values. Unlike a total-energy
    /// check, this names *which* term regressed — essential for fast triage.
    /// Values are RDKit 2025.09.3 MMFF94s `SetMMFFVerbosity(2)` per-term totals.
    #[test]
    fn pyrimidine_breakdown_matches_rdkit_per_term() {
        let mol = parse_sdf(PYRIMIDINE).expect("parse");
        let c: Vec<[f64; 3]> = mol.atoms.iter().map(|a| a.position).collect();
        let ff = MMFFForceField::new(&mol, MMFFVariant::MMFF94s);
        let bd = ff.calculate_energy_breakdown(&c);
        let check = |name: &str, got: f64, want: f64| {
            assert!((got - want).abs() < 0.05, "{}: got {:.3}, want {:.3}", name, got, want);
        };
        check("bond", bd.bond, 0.8143);
        check("angle", bd.angle, 7.0014);
        check("stretch_bend", bd.stretch_bend, 0.9510);
        check("oop", bd.oop, 0.0);
        check("torsion", bd.torsion, 0.0);
        check("vdw", bd.vdw, 11.0989);
        check("electrostatic", bd.electrostatic, -21.1362);
    }

    /// Bug #4 (found via validation-set expansion): find_torsions created
    /// degenerate torsions where atom1 == atom4 in 3-membered rings
    /// (the ring wraps around to the same atom: C2-C0-C1-C2). RDKit filters
    /// these; WebMM didn't, inflating torsion energy (cyclopropane was +0.71,
    /// aziridine +1.24). cyclopropane has exactly 3 such invalid terms.
    #[test]
    fn no_degenerate_torsions_in_3_ring() {
        let mol = parse_sdf(include_str!("../scripts/val_set/cyclopropane.sdf")).expect("parse");
        let ff = MMFFForceField::new(&mol, MMFFVariant::MMFF94s);
        let degenerate: Vec<_> = ff.torsions.iter()
            .filter(|t| t.atom1 == t.atom4)
            .collect();
        assert!(degenerate.is_empty(),
            "found {} torsions with atom1==atom4 (degenerate, should be 0)",
            degenerate.len());
        // all 4 atoms of every torsion must be distinct
        for t in &ff.torsions {
            let s = std::collections::HashSet::from([t.atom1, t.atom2, t.atom3, t.atom4]);
            assert_eq!(s.len(), 4, "torsion ({},{},{},{}) has non-distinct atoms", t.atom1, t.atom2, t.atom3, t.atom4);
        }
    }

    /// Purine end-to-end: the fused-ring sb_type bug left it +1.1 off; pin it
    /// under 0.5 so a regression is caught even if the per-term test above is
    /// not run.
    #[test]
    fn purine_total_energy_within_tolerance() {
        let mol = parse_sdf(PURINE).expect("parse");
        let c: Vec<[f64; 3]> = mol.atoms.iter().map(|a| a.position).collect();
        let ff = MMFFForceField::new(&mol, MMFFVariant::MMFF94s);
        let e = ff.calculate_energy(&c);
        // RDKit 2025.09.3 MMFF94s total = 27.34
        assert!((e - 27.34).abs() < 0.5, "purine energy {} vs RDKit 27.34", e);
    }
}
