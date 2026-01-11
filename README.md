# MolGopt - MMFF94/MMFF94s Geometry Optimizer

WASM-based molecular geometry optimizer for drug-like compounds using ETKDG v3 embedding and MMFF94/MMFF94s force field.

## Project Status

This is a **work in progress** project. Core structure is set up but most functionality is TODO.

## Project Structure

```
molgopt/
├── Cargo.toml                 # Rust project configuration
├── data/
│   └── mmff94_sample_parameters.json  # MMFF94 parameters (sample subset)
├── src/
│   ├── lib.rs              # WASM entry point & API
│   ├── molecule/
│   │   ├── mod.rs         # Module exports
│   │   ├── parser.rs      # SDF/MOL file parsing (✓ Complete)
│   │   └── graph.rs       # Molecular graph analysis (✓ Complete)
│   ├── etkdg/
│   │   └── mod.rs         # ETKDG v3 embedding (⚠️ Placeholder)
│   ├── mmff/
│   │   ├── mod.rs         # MMFF force field orchestrator (⚠️ Partial)
│   │   ├── bond.rs       # Bond stretching (✓ Complete)
│   │   ├── angle.rs      # Angle bending (✓ Complete)
│   │   ├── torsion.rs    # Torsion (✓ Complete)
│   │   ├── oop.rs        # Out-of-plane (✓ Complete)
│   │   ├── vdw.rs        # van der Waals (✓ Complete)
│   │   └── electrostatics.rs  # Electrostatics (✓ Complete)
│   ├── optimizer/
│   │   └── mod.rs         # L-BFGS optimizer (⚠️ Implemented, untested)
│   └── utils/
│       └── mod.rs         # Utilities (⚠️ Minimal)
└── scripts/
    └── generate_test_molecules.py  # Test data generator (✓ Works with RDKit)
```

## Progress

### ✓ Completed
- **Molecule parsing**: SDF/MOL file parser with full support for atoms, bonds, properties
- **Graph analysis**: Adjacency lists, hybridization, aromaticity, angles, torsions, OOP groups
- **MMFF energy terms**: All 6 energy terms implemented with sample parameters:
  - Bond stretching
  - Angle bending
  - Torsion (Fourier series)
  - Out-of-plane bending (Wilson Fourier)
  - Van der Waals (14-7 buffered potential)
  - Electrostatics (Coulomb with dielectric)
- **L-BFGS optimizer**: Full implementation with two-loop recursion and Armijo line search

### ⚠️ Partial / TODO
- **ETKDG v3**: Only placeholder - needs full implementation:
  - Distance bounds matrix
  - Triangle smoothing
  - Coordinate embedding
  - Distance geometry force field
  - Experimental torsions (ET)
  - Knowledge terms (K)
- **MMFF atom typing**: Simple heuristic only - needs complete MMFF94 rules
- **MMFF parameter tables**: Only sample values - need full MMFF94 parameter set
- **MMFF gradients**: Angle, torsion, OOP gradients not yet implemented
- **WASM build**: Not yet set up (requires `cargo install`)

### ❌ Not Started
- **Testing**: No unit tests written yet
- **Validation**: No comparison against RDKit MMFF94 energies
- **Documentation**: API documentation
- **Web interface**: No HTML/JS example yet

## API (Planned)

```rust
// Main WASM entry point
pub fn optimize_from_sdf(
    sdf_content: &str,
    options: OptimizationOptions,
) -> Result<OptimizationResult, String>

// Options
pub struct OptimizationOptions {
    pub mmff_variant: String,  // "MMFF94" or "MMFF94s"
    pub convergence: ConvergenceOptions,
}

// Result
pub struct OptimizationResult {
    pub optimized_coords: Vec<Vec<f64>>>,  // 3D coordinates
    pub final_energy: f64,
    pub converged: bool,
    pub iterations: usize,
    pub message: String,
}
```

## Build Instructions

### Prerequisites

1. **Install Rust** (M2 Mac compatible):
   ```bash
   # Using Homebrew
   brew install rust
   
   # Or use rustup
   curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
   ```

2. **Add WASM target**:
   ```bash
   rustup target add wasm32-unknown-unknown
   cargo install wasm-bindgen-cli
   ```

### Build

```bash
# Build WASM library
cargo build --release --target wasm32-unknown-unknown

# Build with wasm-bindgen (for JS integration)
wasm-pack build --release
```

### Test

```bash
# Run unit tests
cargo test

# Check compilation
cargo check
```

3. Implement coordinate embedding (random + DGFF minimization)
4. Implement DGFF (simple force field)
5. Add experimental torsions (optional initially)
6. Add knowledge terms (optional initially)

### Priority 3: Complete MMFF implementation
1. Implement full MMFF94 atom typing (100+ types)
2. Load parameters from JSON (expand data/mmff94_sample_parameters.json)
3. Implement remaining gradients (angle, torsion, OOP)
4. Add MMFF94s variant parameters

### Priority 4: Testing
1. Write unit tests for each module
2. Test with molecules from `test_molecules/` (from Python script)
3. Validate against RDKit MMFF94 energies

### Priority 5: WASM Integration
1. Ensure wasm-bindgen works
2. Create JavaScript wrapper examples
3. Create simple HTML demo page
4. Test in browser

## Resources

- **Original MMFF94 paper**: Halgren, T.A. J. Comput. Chem. 17, 490-519 (1996)
- **RDKit MMFF implementation**: https://github.com/rdkit/rdkit/tree/master/Code/GraphMol/ForceFieldHelpers/MMFF
- **ETKDG v3 paper**: Wang et al. J. Chem. Inf. Model. 61, 6598-6607 (2020)
- **RDKit ETKDG blog**: https://greglandrum.github.io/rdkit-blog/

## Development Notes
 
- This is a **local development** project - not intended for production use
- Parameter tables are stored locally in JSON (not embedded in binary)
- Target: Drug-like molecules (typical 5-50 heavy atoms)
- Platform optimized for M2 MacBook Air (Apple Silicon)
- Default convergence: RMS force < 0.001, Max force < 0.01 kcal/mol/Å
 
## WASM Integration

### Build Commands

```bash
# Build release version
cargo build --release --target wasm32-unknown-unknown

# Generate JavaScript bindings
wasm-bindgen --out-dir pkg target/wasm32-unknown-unknown/release/molgopt.wasm

# Generated files: pkg/molgopt.js, pkg/molgopt_bg.js, pkg/molgopt.d.ts, pkg/molgopt.wasm
```

### Browser Usage

```javascript
import init, { optimize_from_sdf, OptimizationOptions, MMFFVariant } from './molgopt.js';

// Load WASM module
await init('./molgopt_bg.wasm');

// Configure optimization
const options = new OptimizationOptions();
options.convergence.max_force = 0.01;  // Maximum force component
options.convergence.rms_force = 0.001;  // RMS force
options.convergence.energy_change = 1e-6;  // Energy change threshold
options.convergence.max_iterations = 200;  // Maximum iterations
options.mmff_variant = 'MMFF94s';  // MMFF variant (MMFF94 or MMFF94s)

// Optimize molecule
const result = optimize_from_sdf(sdfContent, options);

if (result.converged()) {
    console.log('✓ Optimization converged!');
    console.log('Final energy:', result.get_final_energy());
    console.log('Iterations:', result.get_iterations());
    console.log('Atoms:', result.get_n_atoms());
    
    // Access optimized coordinates (atom index 0-2, coord index 0-2)
    for (let i = 0; i < result.get_n_atoms(); i++) {
        const x = result.get_coord(i, 0);
        const y = result.get_coord(i, 1);
        const z = result.get_coord(i, 2);
        console.log(`Atom ${i}:`, x, y, z);
    }
} else {
    console.log('✗ Optimization did not converge');
    console.log('Final energy:', result.get_final_energy());
    console.log('Iterations:', result.get_iterations());
    console.log('Message:', result.get_message());
}
```

### Testing

A test page is provided in `pkg/index.html` that can be opened in a browser to test the WASM module.

### API Reference

#### `optimize_from_sdf(sdfContent, options) -> OptimizationResult`

Optimizes a 3D molecular structure from SDF/MOL file content.

**Parameters:**
- `sdfContent`: SDF/MOL file content as string
- `options`: OptimizationOptions configuration

**Returns:**
- `n_atoms`: Number of atoms
- `final_energy`: Final MMFF energy (kcal/mol)
- `converged()`: Whether optimization converged (getter method)
- `iterations()`: Number of optimization iterations (getter method)
- `message`: Status/error message (getter method)

#### `get_coord(atom_idx, coord_idx) -> f64`

Access individual coordinate of an optimized atom.
- `atom_idx`: Atom index (0-based)
- `coord_idx`: Coordinate index (0=x, 1=y, 2=z)

#### `OptimizationOptions`

**Convergence Settings:**
- `max_force`: Maximum force component (default: 0.01)
- `rms_force`: RMS force (default: 0.001)
- `energy_change`: Energy change threshold (default: 1e-6)
- `max_iterations`: Maximum iterations (default: 200)
- `mmff_variant`: MMFF variant ("MMFF94" or "MMFF94s")

#### `MMFFVariant`

Enum for MMFF force field variant:
- `MMFF94`: Standard MMFF94 force field
- `MMFF94s`: Scaled MMFF94s force field (modified for better performance)

## License
 
TBD (decide when ready to share)
