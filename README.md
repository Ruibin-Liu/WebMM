# MolGopt - MMFF94/MMFF94s Geometry Optimizer

WASM-based molecular geometry optimizer for drug-like compounds using ETKDG v3 embedding and MMFF94/MMFF94s force field.

## Project Structure

```
molgopt/
├── Cargo.toml                 # Rust project configuration
├── data/
│   └── mmff94_sample_parameters.json  # MMFF94 parameters (embedded at compile time)
├── src/
│   ├── lib.rs              # WASM entry point & API
│   ├── molecule/
│   │   ├── mod.rs         # Molecule types with cached adjacency
│   │   ├── parser.rs      # SDF/MOL V2000 file parser
│   │   └── graph.rs       # Molecular graph analysis (hybridization, angles, torsions, OOP)
│   ├── etkdg/
│   │   └── mod.rs         # ETKDG v3 3D coordinate embedding
│   ├── mmff/
│   │   ├── mod.rs         # MMFF force field orchestrator (single-pass energy+gradient)
│   │   ├── bond.rs        # Bond stretching (energy + gradient)
│   │   ├── angle.rs       # Angle bending (energy + gradient)
│   │   ├── torsion.rs     # Torsion (energy + numerical gradient)
│   │   ├── oop.rs         # Out-of-plane (energy + numerical gradient)
│   │   ├── vdw.rs         # van der Waals buffered 14-7 (energy + numerical gradient)
│   │   └── electrostatics.rs  # Electrostatics (energy + gradient)
│   ├── optimizer/
│   │   └── mod.rs         # L-BFGS optimizer with correct two-loop recursion
│   └── utils/
│       └── mod.rs         # Parameter loading from embedded JSON
├── pkg/
│   └── index.html         # Browser test page
├── scripts/
│   └── generate_test_molecules.py  # Test data generator (requires RDKit)
└── test_molecules/         # 43 test SDF files
```

## Progress

### Completed
- **Molecule parsing**: Full SDF/MOL V2000 parser with correct column layout, V2000 charge encoding, multi-element support
- **Graph analysis**: Cached adjacency lists, bond-order-aware hybridization, aromaticity, angle/torsion/OOP detection
- **MMFF energy terms**: All 6 terms with correct energy formulas:
  - Bond stretching (harmonic)
  - Angle bending (harmonic with cubic/quartic corrections)
  - Torsion (Fourier series)
  - Out-of-plane bending (Wilson Fourier)
  - Van der Waals (buffered 14-7 potential with attractive well)
  - Electrostatics (Coulomb with dielectric)
- **MMFF gradients**: All 6 terms with verified gradients (numerical for torsion/OOP/VDW, analytical for bond/angle/electrostatics)
- **L-BFGS optimizer**: Correct two-loop recursion, H0 scaling, Armijo line search, energy change convergence
- **ETKDG v3**: Distance bounds, triangle smoothing, 4D stochastic embedding, DGFF minimization
- **WASM API**: Full JavaScript interface with coordinate access, constructor support
- **Parameter loading**: MMFF parameters embedded from JSON at compile time
- **Testing**: 21 tests including numerical gradient verification and end-to-end optimization

### Partial / TODO
- **MMFF atom typing**: Heuristic-based (C, H, N, O, S, P, halogens) -- needs full MMFF94 decision tree
- **MMFF parameters**: Sample subset from JSON -- needs full MMFF94 parameter set
- **Partial charges**: Currently zero (electrostatics disabled) -- needs bond charge increment method
- **ETKDG v3 advanced**: Missing experimental torsions and knowledge-based terms

## Build Instructions

### Prerequisites

1. **Install Rust**:
   ```bash
   curl --proto '=https' --tlsv1.2 -sSf https://shustart rustup.rs | sh
   ```

2. **Add WASM target**:
   ```bash
   rustup target add wasm32-unknown-unknown
   cargo install wasm-bindgen-cli
   ```

### Build

```bash
# Build native library
cargo build --release

# Build WASM library
cargo build --release --target wasm32-unknown-unknown
wasm-bindgen --out-dir pkg target/wasm32-unknown-unknown/release/molgopt.wasm
```

### Test

```bash
cargo test
```

## WASM API

### Browser Usage

```javascript
import { optimize_from_sdf, OptimizationOptions, MMFFVariant } from './molgopt.js';

// Configure optimization
const options = new OptimizationOptions();
options.convergence.max_force = 0.01;
options.convergence.rms_force = 0.001;
options.convergence.energy_change = 1e-6;
options.convergence.max_iterations = 200;
options.mmff_variant = 'MMFF94s';

// Optimize molecule
const result = optimize_from_sdf(sdfContent, options);

if (result.get_converged()) {
    console.log('Final energy:', result.get_final_energy());
    console.log('Iterations:', result.get_iterations());

    // Access optimized coordinates
    for (let i = 0; i < result.n_atoms; i++) {
        const x = result.get_coord(i, 0);
        const y = result.get_coord(i, 1);
        const z = result.get_coord(i, 2);
    }
}
```

### API Reference

#### `optimize_from_sdf(sdfContent, options) -> OptimizationResult`

Optimizes a 3D molecular structure from SDF/MOL file content.

#### `OptimizationResult`

| Field / Method | Type | Description |
|---|---|---|
| `n_atoms` | `usize` | Number of atoms |
| `final_energy` | `f64` | Final MMFF energy (kcal/mol) |
| `get_converged()` | `bool` | Whether optimization converged |
| `get_iterations()` | `usize` | Number of iterations |
| `get_message()` | `String` | Status message |
| `get_coord(atom, dim)` | `f64` | Coordinate (atom 0-based, dim: 0=x, 1=y, 2=z) |
| `get_coordinates()` | `Vec<f64>` | Flat coordinate array [x0,y0,z0,x1,...] |

#### `OptimizationOptions`

| Field | Default | Description |
|---|---|---|
| `mmff_variant` | `"MMFF94s"` | `"MMFF94"` or `"MMFF94s"` |
| `convergence.max_force` | `0.01` | Max force component (kcal/mol/A) |
| `convergence.rms_force` | `0.001` | RMS force |
| `convergence.energy_change` | `1e-6` | Energy change threshold |
| `convergence.max_iterations` | `200` | Max iterations |

## Resources

- **MMFF94 paper**: Halgren, T.A. J. Comput. Chem. 17, 490-519 (1996)
- **ETKDG v3 paper**: Wang et al. J. Chem. Inf. Model. 61, 6598-6607 (2020)
- **RDKit MMFF**: https://github.com/rdkit/rdkit/tree/master/Code/GraphMol/ForceFieldHelpers/MMFF

## License

TBD
