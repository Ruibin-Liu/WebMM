# WebMM — MMFF94/MMFF94s Geometry Optimizer

WASM-based molecular modeling toolkit for drug-like compounds: ETKDG v3
conformer embedding, MMFF94/MMFF94s force field with L-BFGS optimization,
gas-phase molecular dynamics, well-tempered metadynamics, and GBSA implicit
solvation. Runs fully in the browser via WebAssembly.

## Project Structure

```
webmm/
├── Cargo.toml                 # Rust project configuration
├── data/
│   └── mmff94_sample_parameters.json  # MMFF94 parameters (embedded at compile time)
├── src/
│   ├── lib.rs              # WASM entry point & public API
│   ├── forces.rs           # ForceField trait (composable force sources)
│   ├── molecule/
│   │   ├── mod.rs         # Molecule types with cached adjacency
│   │   ├── parser.rs      # SDF/MOL V2000 file parser
│   │   └── graph.rs       # Molecular graph analysis (hybridization, angles, torsions, OOP)
│   ├── etkdg/
│   │   └── mod.rs         # ETKDG v3 3D coordinate embedding
│   ├── mmff/
│   │   ├── mod.rs         # MMFF force field orchestrator (cached term params)
│   │   ├── atom_types.rs  # Atom typing (ring/aromaticity/charge-aware)
│   │   ├── bond.rs        # Bond stretching (energy + gradient)
│   │   ├── angle.rs       # Angle bending (energy + gradient)
│   │   ├── torsion.rs     # Torsion (energy + gradient)
│   │   ├── oop.rs         # Out-of-plane (energy + gradient)
│   │   ├── stretch_bend.rs# Stretch-bend (energy + gradient)
│   │   ├── vdw.rs         # van der Waals buffered 14-7 (energy + gradient)
│   │   ├── electrostatics.rs  # Electrostatics (energy + gradient)
│   │   ├── charges.rs     # MMFF eq. 15 partial charges (BCI)
│   │   ├── estimation.rs  # Parameter estimation fallbacks
│   │   ├── params.rs      # Resolved parameter cache
│   │   └── mmff_tables.rs # Primary parameter tables
│   ├── md/
│   │   └── mod.rs         # Molecular dynamics (velocity-Verlet NVE / BAOAB Langevin NVT)
│   ├── metad/
│   │   └── mod.rs         # Well-tempered metadynamics + FES reconstruction
│   ├── solvation/
│   │   └── mod.rs         # GBSA implicit solvation (OBC2 + LCPO surface area)
│   ├── optimizer/
│   │   └── mod.rs         # L-BFGS optimizer with correct two-loop recursion
│   └── utils/
│       └── mod.rs         # Parameter loading from embedded JSON
├── site/
│   └── index.html         # Committed demo landing page (staged into pkg/ at build)
├── pkg/                   # Gitignored WASM build output (wasm-pack + staged index.html)
├── scripts/               # Python validation/benchmark tooling (RDKit only in scripts)
│   ├── benchmark_mmff.py  # 230-molecule MMFF validation gate vs RDKit
│   ├── validate_etkdg.py  # Multi-seed ETKDG harness vs RDKit
│   ├── val_set*/          # Validation SDF sets + RDKit reference JSONs
│   └── ...                # parameter extraction / audit / diagnostics
├── examples/              # Native diagnostic & benchmark examples
└── docs/                  # Validation and coverage notes
```

## Progress

### Completed
- **Molecule parsing**: Full SDF/MOL V2000 parser with correct column layout, V2000 charge encoding, multi-element support
- **Graph analysis**: Cached adjacency lists, bond-order-aware hybridization, aromaticity, angle/torsion/OOP detection
- **SSSR ring detection**: BFS-based smallest set of smallest rings with canonical deduplication
- **MMFF atom typing**: Context-sensitive assignment using ring membership, aromaticity, formal charge, neighbor C=O/ether context (C_3, C_2, C_1, C_AR, C_CAT, C_AN, N_3, N_2, N_1, N_AR, N_PL3, N_AM, N_4, O_3, O_2, O_R, O_CO2, S_3, S_2, S_AR, P_3, P_4, halogens)
- **Atom type property table**: Full Halgren 1996 Table II properties (cr, phi, Z, anc)
- **MMFF energy terms**: All MMFF terms with correct energy formulas:
  - Bond stretching (harmonic)
  - Angle bending (harmonic with cubic/quartic corrections)
  - Torsion (Fourier series)
  - Out-of-plane bending (Wilson Fourier)
  - Stretch-bend coupling
  - Van der Waals (buffered 14-7 potential with attractive well)
  - Electrostatics (Coulomb with dielectric)
- **MMFF gradients**: All terms with verified analytical gradients (finite-difference-validated)
- **MMFF validation vs RDKit**: **230/230 molecules match RDKit to <0.01 kcal/mol** (atom types, charges, energies — `scripts/benchmark_mmff.py`, a regression gate)
- **L-BFGS optimizer**: Correct two-loop recursion, H0 scaling, Armijo line search, energy change convergence
- **ETKDG v3**: Distance bounds (bond + angle 1-3 + torsion 1-4 + ring closure), triangle smoothing, 4D stochastic embedding, eigenvector 4D-to-3D projection, FF-based refinement with L-BFGS, multi-conformer selection; multi-seed validated vs RDKit (`scripts/validate_etkdg.py`)
- **Molecular dynamics**: Allocation-free force evaluation, velocity-Verlet (NVE) + BAOAB Langevin (NVT) integrators, Maxwell–Boltzmann initialization, deterministic seeded PRNG; live-steppable `MDLive` WASM handle for in-browser trajectory animation
- **Metadynamics**: Well-tempered Gaussian bias on dihedral/distance collective variables, free-energy surface reconstruction; live-steppable `MetaDLive` WASM handle (CV/hills/FES queried between animation frames)
- **GBSA implicit solvation**: Onufriev–Bashford–Case (OBC2) Born radii via exact HCT desolvation integrals, analytical gradient, LCPO surface-area (SA) nonpolar term
- **WASM API**: Full JavaScript interface — optimization, embedding, MD, and metadynamics with trajectory/FES results
- **Parameter loading**: MMFF parameters embedded from JSON at compile time with fallback lookup
- **Testing**: 199 tests including numerical gradient verification, end-to-end optimization, ring detection, V3000 parsing, property-based invariants, atom type assignment, NVE/NVT stability, and edge cases

## Validation

The library is validated against RDKit 2025.09.3 (dev-time tooling only; no
Python runtime dependency in the library):

- **MMFF**: `python3 scripts/benchmark_mmff.py --no-speed` — 230/230 molecules
  match RDKit atom types, charges, and energies to <0.01 kcal/mol. This is the
  regression gate; exit code 0 required.
- **ETKDG**: `python3 scripts/gen_etkdg_ref.py` + `scripts/validate_etkdg.py` —
  multi-seed embedding harness vs RDKit conformers.
- See `docs/atom-type-coverage.md` and `docs/validation-energy-analysis.md`.

## Build Instructions

### Prerequisites

1. **Install Rust**:
   ```bash
   curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
   ```

2. **Add WASM target and install wasm-pack**:
   ```bash
   rustup target add wasm32-unknown-unknown
   npm install    # installs the pinned wasm-pack devDependency (see package.json)
   ```

### Build

**Option 1: npm + wasm-pack (recommended)**

```bash
npm run build   # wasm-pack build --target web --out-dir pkg
```

**Option 2: Manual (uses wasm-bindgen-cli, must match the wasm-bindgen version in Cargo.toml)**

```bash
# Build native library
cargo build --release

# Build WASM library
cargo build --release --target wasm32-unknown-unknown
wasm-bindgen --out-dir pkg --target web target/wasm32-unknown-unknown/release/webmm.wasm
```

The WASM release profile enables `wasm32-simd128` via `.cargo/config.toml`
(scoped to the wasm target — native builds/tests are unaffected). Safari
16.4+ / Chrome 91+ / Firefox 89+ required.

### Serve the demo locally

The WASM bindings are written to the gitignored `pkg/` directory, while the
demo landing page is the committed `site/index.html`. To serve the demo,
build the WASM, stage the landing page into `pkg/`, then serve `pkg/`:

```bash
npm run build                              # wasm-pack build --target web --out-dir pkg
cp site/index.html pkg/index.html          # stage the landing page
cp site/playground.html pkg/playground.html  # stage the playground page
python3 -m http.server 8000 --directory pkg   # serve on http://localhost:8000
```

Then open **http://localhost:8000** in your browser. (This mirrors what the
GitHub Pages workflow does — see `.github/workflows/pages.yml`.)

The **playground** (`/playground.html`) is an interactive physics toy on top of
the live MD engine: drag atoms with the pointer (the force field fights back),
twist bonds to rotate dihedrals, crank the temperature, watch per-atom force
glow, and run metadynamics with a live free-energy surface. It uses the
`MDLive` perturbation exports (`set_atom_position` / `rescale_temperature` /
`force_magnitudes`).

> **Note:** the demo loads `3Dmol.js` (the 3D viewer) from a CDN, so you need
> internet access for the viewer. The MMFF/ETKDG/MD/metadynamics engines run
> fully in-browser via WASM — no network calls.

#### Using the demo

- **Molecule buttons** — load any built-in RDKit-generated 3D structure
  (caffeine, aspirin, benzene, …) into the viewer.
- **Render from input** — re-render whatever SDF is in the textarea.
- **🌐 Embed 3D (ETKDG)** — regenerate coordinates from scratch using the
  ETKDG v3 embedding algorithm.
- **⚙ Optimize (MMFF94s)** — run the full MMFF94s minimization. The readout
  panel shows final energy (kcal/mol), iteration count, atom count, and time.
- **🔥 Run Molecular Dynamics** — gas-phase MMFF MD; NVT (BAOAB Langevin) with
  the thermostat on, NVE otherwise. Trajectory playback with energy/temperature
  charts.
- **🔬 Run Metadynamics** — well-tempered metadynamics on a dihedral or
  distance collective variable, with a free-energy surface (FES) plot plus
  deposited-hill markers and CV-over-time overlay.
- **Viewer style buttons** — switch between Stick / Ball+Stick / Space Fill /
  Wireframe.

#### Stopping the server

```bash
lsof -ti tcp:8000 | xargs kill
```

#### Rebuilding after edits

After changing Rust source or `site/index.html`, re-run the three commands
above (build → stage → serve).

### Test

```bash
cargo test          # 199 tests
cargo clippy --all-targets   # must stay at 0 warnings
```

## WASM API

### Browser Usage

```javascript
import { optimize_from_sdf, run_md_from_sdf, run_metadynamics_from_sdf,
         OptimizationOptions, MDOptions, MetaDOptions, MMFFVariant } from './webmm.js';

// Optimize a molecule (MMFF94s minimization)
const options = new OptimizationOptions();
options.convergence.max_force = 0.01;
options.convergence.rms_force = 0.001;
options.convergence.energy_change = 1e-6;
options.convergence.max_iterations = 200;
options.mmff_variant = 'MMFF94s';

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

Optimizes a 3D molecular structure from SDF/MOL file content. Uses the SDF
coordinates directly when they are already 3D, otherwise ETKDG-embeds (seed 42)
first.
`optimize_from_sdf_direct(sdfContent, options)` skips ETKDG entirely and
optimizes the SDF coordinates as-is (2D SDFs are optimized from the flat
z=0 plane). `generate_initial_coordinates_wasm(sdfContent) -> ETKDGResult` runs
ETKDG v3 embedding only (returns coordinates without MMFF minimization).

#### `run_md_from_sdf(sdfContent, options) -> MDResult`

Runs gas-phase MMFF molecular dynamics and returns a sampled trajectory
(NVT BAOAB Langevin if `friction_per_ps > 0`, else NVE).

#### `run_metadynamics_from_sdf(sdfContent, options) -> MetaDResult`

Runs well-tempered metadynamics and returns the trajectory, CV trace, hill
centers, and a reconstructed free-energy surface.

#### Live handles: `MDLive` / `MetaDLive` (stateful stepping)

For live trajectory animation (the demo's MD / metadynamics buttons), construct
`new MDLive(sdf, options)` / `new MetaDLive(sdf, options)`, then advance the
simulation in small chunks — `step(nSteps)` — and read
`coords()`, `potential_energy()`, `temperature()`, `time_fs()` between
animation frames. `MetaDLive` additionally exposes `last_cv()`, `hill_count()`,
`hill_centers()`, and `fes_s(gridPoints)` / `fes_f(gridPoints)` for the
free-energy surface. Stepping is deterministic and identical to the
corresponding batch API.

`MDLive` also exposes interactive-perturbation methods (for the playground
page): `set_atom_position(i, x, y, z)` (kinematic dragging — sets the
position, zeroes the atom's velocity, and refreshes the cached energy/forces),
`rescale_temperature(t)` (exact instantaneous velocity rescale + thermostat
retarget; re-initializes from Maxwell–Boltzmann at ~0 K; `t <= 0` freezes),
and `force_magnitudes()` (per-atom |F| in kcal/mol/Å from the last force
evaluation).

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
| `get_success()` / `get_error()` | `bool` / `String` | Failure reporting |

#### `OptimizationOptions`

| Field | Default | Description |
|---|---|---|
| `mmff_variant` | `"MMFF94s"` | `"MMFF94"` or `"MMFF94s"` |
| `convergence.max_force` | `0.01` | Max force component (kcal/mol/A) |
| `convergence.rms_force` | `0.001` | RMS force |
| `convergence.energy_change` | `1e-6` | Energy change threshold |
| `convergence.max_iterations` | `200` | Max iterations |

#### `MDOptions` / `MetaDOptions`

| Field | Description |
|---|---|
| `mmff_variant` | `"MMFF94"` or `"MMFF94s"` |
| `dt_fs` | Integration timestep in fs |
| `n_steps` | Total integration steps |
| `temperature_k` | Target temperature (K) |
| `friction_per_ps` | Langevin friction (1/ps); 0 → NVE |
| `seed` | PRNG seed (deterministic trajectories) |
| `snapshot_interval` | Frames per recorded snapshot |
| `cv_type` / `cv_atoms` (metad) | `"dihedral"` or `"distance"` CV over atom indices |
| `hill_height` / `hill_width` / `deposit_interval` / `bias_factor` (metad) | Well-tempered hill parameters |
| `fes_grid_points` (metad) | FES grid resolution |

#### `MDResult`

| Field / Method | Type | Description |
|---|---|---|
| `n_atoms` / `n_frames` | `usize` | Trajectory shape |
| `coordinates()` | `Vec<f64>` | Flattened frames [frame, atom, xyz] |
| `energies()` / `temperatures()` / `times_fs()` | `Vec<f64>` | Per-frame properties |
| `final_energy()` / `final_temperature()` | `f64` | Final state |
| `get_coord(frame, atom, axis)` | `f64` | Single coordinate access |
| `success()` / `error()` | `bool` / `String` | Failure reporting |

`MetaDResult` additionally exposes `cv_values()`, `n_hills()`, `fes_s()`,
and `fes_f()` (CV trace and reconstructed free-energy surface).

## Resources

- **MMFF94 paper**: Halgren, T.A. J. Comput. Chem. 17, 490-519 (1996)
- **ETKDG v3 paper**: Wang et al. J. Chem. Inf. Model. 61, 6598-6607 (2020)
- **RDKit MMFF**: https://github.com/rdkit/rdkit/tree/master/Code/GraphMol/ForceFieldHelpers/MMFF
- **Validation notes**: `docs/atom-type-coverage.md`, `docs/validation-energy-analysis.md`

## License

TBD
