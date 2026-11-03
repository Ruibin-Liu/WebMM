// m22_api_surface.test.js — interface stabilization guard (feature freeze).
// Pure Node: parses pkg/webmm.d.ts, the three workers, the task contract and
// the IMTD driver, and compares every public surface against the committed
// GOLDEN snapshots below. ANY change (add, remove, rename) fails — updating
// the golden is a deliberate act that must come with a version bump per
// docs/STABILITY.md.
const fs = require('fs');
const path = require('path');

const ROOT = path.join(__dirname, '..', '..');
let pass = 0, fail = 0;
function check(name, cond, detail) {
  if (cond) { pass++; console.log('  ✓ ' + name); }
  else { fail++; console.log('  ✗ ' + name + (detail ? ' — ' + detail : '')); }
}

// ── GOLDEN: wasm export surface (v1.12.0 freeze line) ──
const GOLDEN_FUNCTIONS = [
  
  "add_hydrogens_wasm",
  "attach_hydrogens_3d_wasm",
  "ecfp4_fingerprint_wasm",
  "energy_terms_wasm",
  "generate_conformers_wasm",
  "generate_initial_coordinates_wasm",
  "generate_optimized_conformers_wasm",
  "imtd_boltzmann_populations_g_json",
  "imtd_boltzmann_populations_json",
  "kabsch_align_wasm",
  "optimize_from_sdf",
  "optimize_from_sdf_direct",
  "run_imtd_conformer_search",
  "run_imtd_metad_run_json",
  "run_imtd_neb_json",
  "run_md_from_sdf",
  "run_metadynamics_from_sdf",
  "sa_load_table_wasm",
  "sa_score_wasm",
  "shape_align_color_wasm",
  "shape_align_wasm",
  "tanimoto_wasm",
  "webmm_version"
];
const GOLDEN_CLASSES = {
  "ConformerBatch": [
    "get_coordinates",
    "get_error",
    "get_n_atoms",
    "get_n_confs",
    "get_success"
  ],
  "ConvergenceOptions": [],
  "ETKDGResult": [
    "get_coordinates",
    "get_error",
    "get_n_atoms",
    "get_success"
  ],
  "ImtdOptions": [
    "constructor",
    "set_bias_factor",
    "set_convergence_new_minima",
    "set_crossing_pairs",
    "set_deposit_interval",
    "set_dt_fs",
    "set_energy_window",
    "set_friction_per_ps",
    "set_hill_height",
    "set_hill_width",
    "set_hill_width_distance",
    "set_max_conformers",
    "set_max_rounds",
    "set_mmff_variant",
    "set_quench_interval",
    "set_quench_max_iterations",
    "set_rmsd_tol",
    "set_runs_per_round",
    "set_seed",
    "set_steps_per_run",
    "set_temperature_k"
  ],
  "ImtdSearch": [
    "constructor",
    "done",
    "error",
    "gather",
    "gather_refine",
    "gather_verify",
    "landscape_json",
    "plan_json",
    "plan_refine_json",
    "plan_verify_json",
    "result",
    "round",
    "success",
    "thermo_json"
  ],
  "ImtdSearchResult": [
    "converged",
    "coordinates",
    "energies",
    "error",
    "n_atoms",
    "n_confs",
    "n_rotatable",
    "new_minima_per_round",
    "pool_size",
    "round_found",
    "rounds_executed",
    "success"
  ],
  "MDLive": [
    "constructor",
    "coords",
    "error",
    "force_magnitudes",
    "masses",
    "n_atoms",
    "potential_energy",
    "rescale_temperature",
    "set_atom_position",
    "step",
    "steps_done",
    "success",
    "temperature",
    "time_fs"
  ],
  "MDOptions": [
    "constructor",
    "set_dt_fs",
    "set_friction_per_ps",
    "set_mmff_variant",
    "set_n_steps",
    "set_seed",
    "set_snapshot_interval",
    "set_temperature_k"
  ],
  "MDResult": [
    "coordinates",
    "energies",
    "error",
    "final_energy",
    "final_temperature",
    "get_coord",
    "n_atoms",
    "n_frames",
    "steps",
    "success",
    "temperatures",
    "times_fs"
  ],
  "MetaDLive": [
    "constructor",
    "coords",
    "error",
    "fes_f",
    "fes_s",
    "hill_centers",
    "hill_count",
    "last_cv",
    "masses",
    "n_atoms",
    "potential_energy",
    "step",
    "steps_done",
    "success",
    "temperature",
    "time_fs"
  ],
  "MetaDOptions": [
    "constructor",
    "set_bias_factor",
    "set_cv_atoms",
    "set_cv_type",
    "set_deposit_interval",
    "set_dt_fs",
    "set_fes_grid_points",
    "set_friction_per_ps",
    "set_hill_height",
    "set_hill_width",
    "set_mmff_variant",
    "set_n_steps",
    "set_seed",
    "set_snapshot_interval",
    "set_temperature_k"
  ],
  "MetaDResult": [
    "coordinates",
    "cv_values",
    "energies",
    "error",
    "fes_f",
    "fes_s",
    "final_energy",
    "hill_centers",
    "n_atoms",
    "n_frames",
    "n_hills",
    "success",
    "times_fs"
  ],
  "OptimizationOptions": [
    "constructor",
    "set_energy_change",
    "set_max_force",
    "set_max_iterations",
    "set_restraints",
    "set_rms_force"
  ],
  "OptimizationResult": [
    "get_charges",
    "get_converged",
    "get_coord",
    "get_energy_terms_json",
    "get_engine",
    "get_final_energy",
    "get_iterations",
    "get_message"
  ],
  "OptimizedConformers": [
    "get_converged",
    "get_coordinates",
    "get_energies",
    "get_error",
    "get_iterations",
    "get_n_atoms",
    "get_n_confs",
    "get_n_heavy",
    "get_seeds",
    "get_success",
    "get_template_sdf"
  ]
};
// ── GOLDEN: worker message protocol ──
const GOLDEN_WORKERS = {
  "conf.worker.js": [
    "imtdNeb",
    "imtdRun",
    "init",
    "optimize",
    "run"
  ],
  "shape.worker.js": [
    "full",
    "initq",
    "prep",
    "screen"
  ],
  "analog.worker.js": [
    "heap",
    "init",
    "score"
  ]
};
// ── GOLDEN: task contract kinds + envelope ──
const GOLDEN_KINDS = ['analogScan', 'conformers', 'imtd', 'optimize', 'shapeSearch'];
// ── GOLDEN: IMTD driver surface ──
const GOLDEN_IMTD_DRIVER = ['DEFAULTS', 'optsToWasm', 'runParallel', 'start'];

const dts = path.join(ROOT, 'pkg', 'webmm.d.ts');
if (!fs.existsSync(dts)) {
  console.log('RESULT: 0 passed, 1 failed — pkg/webmm.d.ts missing (build first)');
  process.exit(1);
}
const src = fs.readFileSync(dts, 'utf8');

// extract exported functions and class methods
const INTERNALS = new Set(['free', 'toJSON', 'toString', 'clone', 'delete', 'init', 'initSync']);
const fns = [...src.matchAll(/^export function ([A-Za-z_][A-Za-z0-9_]*)/gm)].map(m => m[1])
  .filter(n => !INTERNALS.has(n)).sort();
const classRe = /^export class ([A-Za-z_][A-Za-z0-9_]*) \{([\s\S]*?)^\}/gm;
const classes = {};
let cm;
while ((cm = classRe.exec(src)) !== null) {
  const methods = [...cm[2].matchAll(/^\s{2}([A-Za-z_][A-Za-z0-9_]*)\(/gm)].map(m => m[1]);
  classes[cm[1]] = [...new Set(methods)].filter(n => !INTERNALS.has(n)).sort();
}
const classSort = o => Object.fromEntries(Object.entries(o).map(([k, v]) => [k, [...v].sort()]));

check('wasm export functions == golden',
  JSON.stringify(fns) === JSON.stringify([...GOLDEN_FUNCTIONS].sort()),
  diffList(fns, GOLDEN_FUNCTIONS));
const gClass = classSort(GOLDEN_CLASSES);
for (const [name, methods] of Object.entries(classes).sort()) {
  const want = gClass[name];
  if (!want) { check(`class ${name} known`, false, 'not in golden'); continue; }
  check(`class ${name} methods == golden`, JSON.stringify(methods) === JSON.stringify(want),
    diffList(methods, want));
}
for (const name of Object.keys(gClass)) {
  if (!classes[name]) check(`golden class ${name} present`, false, 'missing from d.ts');
}

// workers
for (const [wf, golden] of Object.entries(GOLDEN_WORKERS)) {
  const wsrc = fs.readFileSync(path.join(ROOT, 'app', wf), 'utf8');
  const types = [...wsrc.matchAll(/msg\.type\s*(?:[!=]==?)\s*'([A-Za-z]+)'/g)].map(m => m[1]);
  const uniq = [...new Set(types)].sort();
  check(`worker ${wf} message types == golden`, JSON.stringify(uniq) === JSON.stringify([...golden].sort()),
    diffList(uniq, golden));
}

// task contract
const C = require(path.join(ROOT, 'app', 'tasks', 'contract.js'));
check('task KINDS == golden', JSON.stringify([...C.KINDS].sort()) === JSON.stringify([...GOLDEN_KINDS].sort()),
  diffList(C.KINDS, GOLDEN_KINDS));

// IMTD driver surface (parse the window.IMTD assignment)
const drv = fs.readFileSync(path.join(ROOT, 'app', 'tasks', 'imtd.js'), 'utf8');
const im = drv.match(/window\.IMTD = \{([^}]*)\}/);
const imtdKeys = im ? im[1].split(',').map(s => s.trim().split(/:\s|\s/)[0]).filter(Boolean).sort() : [];
check('window.IMTD surface == golden', JSON.stringify(imtdKeys) === JSON.stringify([...GOLDEN_IMTD_DRIVER].sort()),
  diffList(imtdKeys, GOLDEN_IMTD_DRIVER));
// handle surface (returned object from start())
const handleOk = /return \{[\s\S]{0,200}?promise,[\s\S]{0,80}?verify,[\s\S]{0,80}?refine,[\s\S]{0,120}?thermoJson:/.test(drv);
check('driver handle surface (promise/verify/refine/thermoJson)', handleOk);

// version consistency
const cargoV = fs.readFileSync(path.join(ROOT, 'Cargo.toml'), 'utf8').match(/^version = "(.+)"$/m)[1];
const pkgV = JSON.parse(fs.readFileSync(path.join(ROOT, 'package.json'), 'utf8')).version;
check('Cargo.toml == package.json version', cargoV === pkgV, `${cargoV} vs ${pkgV}`);
check('d.ts exposes webmm_version', /webmm_version/.test(src));

function diffList(a, b) {
  const A = new Set(a), B = new Set(b);
  const miss = [...B].filter(x => !A.has(x)), extra = [...A].filter(x => !B.has(x));
  return (miss.length ? 'missing: ' + miss.join(',') + ' ' : '') + (extra.length ? 'unexpected: ' + extra.join(',') : '');
}

console.log(`\nRESULT: ${pass} passed, ${fail} failed`);
process.exit(fail ? 1 : 0);
