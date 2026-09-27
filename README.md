# WebMM

**Molecular mechanics in the browser — no server, no install.**

WebMM is a Rust molecular-modeling engine compiled to WebAssembly. It embeds
3D conformers (ETKDG v3), optimizes geometries (MMFF94/MMFF94s, GFN-FF), and
runs molecular dynamics and well-tempered metadynamics entirely client-side.
Structures never leave the user's machine.

<!-- TODO: replace with a real GIF of the playground (drag atoms → optimize → MD → metad FES) -->
<!-- ![WebMM demo](docs/assets/demo.gif) -->

🔗 **Live demo:** <https://ruibin-liu.github.io/WebMM/>

---

## What it does

| Task | Method | Scope & notes |
|---|---|---|
| 3D embedding | ETKDG v3 distance geometry | From SDF/MOL connectivity; stereo-aware, seeded |
| Conformer ensembles | Embed + optimize + RMSD prune | Batch single-call pipeline (30 conformers of a 33-atom drug in ~2 s in-browser) |
| Geometry optimization | MMFF94/MMFF94s + L-BFGS/BFGS | Validated 230/230 vs RDKit to <0.01 kcal/mol |
| Broad-coverage force field | GFN-FF | For elements/patterns MMFF94 does not parameterize (incl. metal coordination complexes); validated vs xtb |
| Molecular dynamics | Velocity-Verlet (NVE) / BAOAB Langevin (NVT) | Live-steppable for in-browser trajectory animation |
| Enhanced sampling | Well-tempered metadynamics | Dihedral and distance collective variables, FES reconstruction |
| Implicit solvation | GBSA (OBC2 + LCPO) | Optional add-on; gas phase by default |
| Ligand workbench | `app/` | Draw (JSME) → descriptors (RDKit-js) → embed/optimize/rank conformers → export; offline-capable |

**Not covered:** proteins, periodic systems, QM beyond GFN-FF, explicit
solvent. MMFF typing refuses metal-bonded systems exactly like RDKit's MMFF
(NULL, never garbage-typed).

## Architecture

```
┌──────────────────────────────────────────────┐
│  app/   ligand CADD workbench (JSME + RDKit-js)
│  site/  engine demo + interactive playground
├──────────────────────────────────────────────┤
│  pkg/   WASM bindings (wasm-bindgen, --target web)
├──────────────────────────────────────────────┤
│  src/   Rust core: parsing, graph analysis, typing,
│         MMFF + GFN-FF, ETKDG, optimizers, MD/metad
└──────────────────────────────────────────────┘
```

`pkg/` is generated (gitignored). Python appears only in dev-time validation
scripts (`scripts/`); the shipped library has no Python dependency.


## Quick start

### Browser support

WebAssembly with `simd128`: **Safari 16.4+, Chrome 91+, Firefox 89+.**

### Option 0 — use the compiled engine (no toolchain, no build)

The Pages deployment serves the engine directly, so you can import it from
any web page:

```html
<script type="module">
  import init, { OptimizationOptions, optimize_from_sdf } from
    "https://ruibin-liu.github.io/WebMM/webmm.js";
  await init();   // fetches + instantiates webmm_bg.wasm
  // ... same API as below
</script>
```

Requirements: a modern browser and network access **at load time** (the
engine itself runs entirely locally afterwards). Caveat: this URL always
serves the **latest** deployed build — there is no version pinning. For a
pinned version, download `webmm.js` + `webmm_bg.wasm` once and self-host,
or build from source.

### Building from source

The remaining paths share these prerequisites:

- [Rust](https://rustup.rs/) stable + `wasm32-unknown-unknown` target
- [wasm-pack](https://rustwasm.github.io/wasm-pack/) — pinned via `npm install` (devDependency)
- Node.js ≥ 18 (only for the pinned toolchain)
- Python 3 (only as a static file server)

```bash
git clone https://github.com/Ruibin-Liu/WebMM.git
cd WebMM
rustup target add wasm32-unknown-unknown
npm install
```

### Option A — the workbench app

```bash
npm run build                      # wasm-pack -> pkg/
python3 -m http.server 8000        # from the repo root
# open http://localhost:8000/app/index.html
```

Draw a molecule or paste SMILES / a multi-record SDF; the batch mode
enumerates, embeds, optimizes and exports 3D SDF + CSV.

### Option B — the engine demo / playground

```bash
npm run build
cp site/index.html site/playground.html pkg/
python3 -m http.server 8000 --directory pkg
# open http://localhost:8000/        (demo)
#      http://localhost:8000/playground.html   (interactive physics toy)
```

The playground adds live MD with atom dragging, dihedral twisting, force
visualization, and metadynamics with a live FES.

> Serve over HTTP — browsers block WASM from `file://` URLs. The 3D viewer
> (3Dmol.js) loads from a CDN; the engine itself makes no network calls.

### Option C — native Rust library

```bash
cargo build --release
cargo test --release        # 281 tests
```

## Usage (JavaScript)

Real exports from `pkg/webmm.d.ts` — all molecule input is SDF/MOL text
(V2000; SMILES input is an app-layer concern via JSME + RDKit-js):

```js
import init, {
  OptimizationOptions, optimize_from_sdf,
  generate_optimized_conformers_wasm,   // batch: embed -> attachH -> optimize
  run_md_from_sdf, run_metadynamics_from_sdf,
  energy_terms_wasm,
} from "./pkg/webmm.js";

await init();

// 1. Optimize (MMFF94s; the engine ETKDG-embeds 2D input first)
const opt = new OptimizationOptions();
opt.mmff_variant = "MMFF94s";
opt.convergence.max_iterations = 1000;
const r = optimize_from_sdf(sdfText, opt);
if (r.get_converged()) console.log(r.get_final_energy(), r.get_iterations());

// 2. Conformer ensemble in one call (30 conformers, MMFF94s, 100-iter protocol)
const batch = generate_optimized_conformers_wasm(sdfText, 30, 42n, "MMFF94s", 100);
const energies = batch.get_energies();       // Float64Array, kcal/mol
const coords   = batch.get_coordinates();    // flat [conf][atom][xyz]

// 3. MD (BAOAB Langevin NVT when friction_per_ps > 0, else NVE)
const md = run_md_from_sdf(sdfText, mdOptions /* see below */);
// 4. Metadynamics (dihedral CV over 4 atoms, well-tempered)
const meta = run_metadynamics_from_sdf(sdfText, metaOptions /* see below */);
```

Key option tables (`MDOptions` / `MetaDOptions`): `dt_fs`, `n_steps`,
`temperature_k`, `friction_per_ps`, `seed`, `snapshot_interval`; metad adds
`cv_type` (`"dihedral"` | `"distance"`), `cv_atoms`, `hill_height`,
`hill_width`, `deposit_interval`, `bias_factor`, `fes_grid_points`.

**Live handles** for animation: `new MDLive(sdf, opts)` / `new MetaDLive(...)`
— `step(n)` between frames, read `coords()`, `temperature()`, and (metad)
`last_cv()`, `hill_count()`, `fes_s(n)`. `MDLive` also supports interactive
perturbation (`set_atom_position`, `rescale_temperature`,
`force_magnitudes`).

Full field-level API reference: see `pkg/webmm.d.ts` after a build.

---

## Validation

All references reproduced identically by RDKit 2025.09.3 / 2026.03.6 and
xtb 6.7.1 (dev-time tooling only).

| Check | Reference | Result |
|---|---|---|
| MMFF94s atom typing, charges, energies | RDKit | **230/230 molecules < 0.01 kcal/mol** (regression gate: `scripts/benchmark_mmff.py`) |
| MMFF94s single-point locks | RDKit | 97 molecules at 0.00000 kcal/mol (90 neutral + 7 charged) |
| GFN-FF per-term energies | xtb | 32/32 organics **< 4×10⁻⁶ Eh**; 5 metal complexes (Ni(CO)₄ exact … Zn(NH₃)₄²⁺ 5.5×10⁻² Eh) |
| ETKDG + MMFF ensembles (30 seeds × 6 molecules) | RDKit | 6/6 global minima bit-identical |
| Gradient audit (analytic vs FD, perturbed geometries) | internal | 230-molecule corpus + GFN-FF fixtures: 0 failures (worst 2×10⁻⁵ rel) |

```bash
cargo test --release            # the full suite
python3 scripts/benchmark_mmff.py --no-speed    # the RDKit gate (needs RDKit)
```

Details: `docs/atom-type-coverage.md`, `docs/validation-energy-analysis.md`,
`docs/gfnff-porting-notes.md`.

## Performance (measured, in-browser)

WASM vs **native** RDKit C++ / xtb on an M-series laptop (2026-09):

| Task (ibuprofen, 33 atoms) | WebMM WASM | Reference (native) |
|---|---|---|
| MMFF optimization | 7.7 ms | RDKit 9.9 ms |
| Conformer pipeline (per conformer) | ~14 ms | RDKit ~17 ms |
| ETKDG embedding | ~11 ms | RDKit ~12 ms |
| MD step (MMFF, NVT) | ~27 µs (≈ 37k steps/s) | — |

Metadynamics per-step cost is flat in the number of deposited hills
(4σ-truncated bias); long runs do not degrade.

## Limitations

- **Size.** MD is practical for small molecules (tens of atoms; ~37k steps/s
  at 33 atoms MMFF, ~4k steps/s GFN-FF). Not for proteins or production
  sampling.
- **Force-field coverage.** MMFF94 parameterizes main-group organic
  chemistry; unusual bonding is out of scope. GFN-FF extends coverage
  (including many metal complexes) with reduced accuracy — treat results
  with caution.
- **Implicit solvent only.** GBSA is an approximation, not a substitute for
  explicit solvent.
- **Physics scope.** Classical force fields; no QM, no PBC, no reactions.
- **Maturity.** Aimed at education, quick exploration, and prototyping.
  Validate important results with established packages.

---

<details>
<summary><b>Project structure</b></summary>

```
WebMM/
├── src/         Rust core (molecule/, mmff/, gfnff/, etkdg/, optimizer/, md/, metad/, solvation/)
├── pkg/         wasm-pack output (generated, gitignored)
├── site/        demo + playground pages
├── app/         ligand workbench
├── tests/       fixtures + reference data
├── scripts/     Python validation tooling (RDKit/xtb refs; dev-time only)
├── examples/    native benchmarks & audits (bench_mmff, grad_audit, md_audit, …)
└── docs/        validation & porting notes
```

</details>

## Contributing

Issues and pull requests welcome. Before submitting:

```bash
cargo fmt && cargo clippy --all-targets && cargo test
```

(`clippy` must stay at 0 warnings; the RDKit benchmark gate must pass when
MMFF code changes.)

## References

1. Halgren, T. A. Merck molecular force field. I–V. *J. Comput. Chem.* **1996**, *17*, 490–641.
2. Halgren, T. A. MMFF VI. *J. Comput. Chem.* **1999**, *20*, 720–729.
3. Wang, W. et al. ETKDGv3. *J. Chem. Inf. Model.* **2020**, *61*, 6598–6607.
4. Spicher, S.; Grimme, S. Robust atomistic modeling of materials, organometallic, and biochemical systems. *Angew. Chem. Int. Ed.* **2020**, *59*, 15665–15673. (GFN-FF)
5. Onufriev, A.; Bashford, D.; Case, D. A. *Proteins* **2004**, *55*, 383–394. (GBSA OBC)
6. Leimkuhler, B.; Matthews, C. *Appl. Math. Res. Express* **2013**, 34–56. (BAOAB)
7. Laio, A.; Parrinello, M. *PNAS* **2002**, *99*, 12562–12566. (Metadynamics)
8. Liu, D. C.; Nocedal, J. *Math. Program.* **1989**, *45*, 503–528. (L-BFGS)

## Citation

```bibtex
@software{webmm,
  author = {Liu, Ruibin},
  title  = {WebMM: Molecular mechanics in the browser},
  url    = {https://github.com/Ruibin-Liu/WebMM},
  year   = {2026}
}
```

## License

[MIT](LICENSE) © 2026 Ruibin Liu
