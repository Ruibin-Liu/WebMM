// Conformer farm worker (ES module): embeds its seed share with ETKDG,
// attaches hydrogens in 3D, optimizes with the selected engine, and streams
// each finished conformer back to the main thread.
//
// Self-contained copies of sdfLines/buildSdfFromCoords (the main thread keeps
// its own; ~30 lines duplicated to keep the worker a single import).
// The wasm module URL is layout-dependent (repo root vs deployed pkg root);
// the main thread resolves it and hands it over in the init message.
let wasm = null;
let ready = false;
let initPromise = null;   // serializes init vs run (async onmessage doesn't queue)

function sdfLines(sdf) {
  const L = sdf.replace(/\r/g, '').split('\n');
  const ci = L.findIndex(l => l.includes('V2000'));
  if (ci === 3) return L;
  if (ci > 3) return L.slice(ci - 3);
  for (let k = 3 - Math.max(ci, 0); k > 0; k--) L.unshift('');
  return L;
}

function buildSdfFromCoords(coords, sdf) {
  const lines = sdfLines(sdf);
  const na = parseInt(lines[3].substring(0, 3));
  const nb = parseInt(lines[3].substring(3, 6));
  const header = lines.slice(0, 4);
  const atomLines = [];
  for (let i = 0; i < na; i++) {
    const x = coords[i * 3].toFixed(4).padStart(10);
    const y = coords[i * 3 + 1].toFixed(4).padStart(10);
    const z = coords[i * 3 + 2].toFixed(4).padStart(10);
    atomLines.push(`${x}${y}${z}${lines[4 + i].substring(30)}`);
  }
  const bondLines = lines.slice(4 + na, 4 + na + nb);
  // Preserve everything after the bond block — M CHG (formal charges!),
  // other M properties and M END. Dropping this tail silently neutralized
  // charged molecules at the 3D stage (e.g. piperazinium [NH+]).
  const tail = lines.slice(4 + na + nb);
  return [...header, ...atomLines, ...bondLines, ...(tail.length ? tail : ['M  END'])].join('\n');
}

self.onmessage = async (e) => {
  const msg = e.data;
  if (msg.type === 'init') {
    initPromise = (async () => {
      wasm = await import(msg.wasmUrl);
      await wasm.default();   // default export = the loader
      ready = true;
    })();
    await initPromise;
    self.postMessage({ type: 'ready', seedBase: msg.seedBase });
    return;
  }
  // Single-structure optimization off the main thread: keeps the UI
  // responsive for large molecules (remdesivir ~19 s on the main thread).
  if (msg.type === 'optimize') {
    try {
      await initPromise;
      const { sdf, engine, maxIter } = msg;
      const opts = new wasm.OptimizationOptions();
      opts.engine = engine;
      opts.set_max_iterations(maxIter);
      const res = wasm.optimize_from_sdf(sdf, opts);
      const coords = [];
      for (let a = 0; a < res.n_atoms; a++)
        for (let d = 0; d < 3; d++) coords.push(res.get_coord(a, d));
      const charges = [];
      for (let a = 0; a < res.n_atoms; a++) charges.push(res.get_charges()[a]);
      self.postMessage({
        type: 'optimized',
        E: res.final_energy,
        converged: res.converged,
        iterations: res.iterations,
        message: res.get_message(),
        engine: res.get_engine(),
        nAtoms: res.n_atoms,
        coords,
        charges,
        termsJson: res.get_energy_terms_json(),
      });
    } catch (err) {
      self.postMessage({ type: 'error', message: String(err && err.message ? err.message : err) });
    }
    return;
  }
  if (msg.type !== 'run') return;
  try {
    await initPromise;   // 'run' waits for init (async handlers don't queue)
    const { sdfHeavy, seedBase, count, engine, maxIter } = msg;
    const t0 = Date.now();
    // Native batch pipeline (v1.2.0): embed -> attach H -> optimize runs
    // entirely inside one WASM call per chunk — flat arrays cross the
    // boundary instead of two SDF string round-trips per conformer, and the
    // MMFF force field is built once and reused across all conformers.
    // Chunks keep seeds continuous (Rust validates 1..=500 per call) and
    // preserve the streaming progress updates.
    const CHUNK = 500;
    let template = null;
    for (let off = 0; off < count; off += CHUNK) {
      const c = Math.min(CHUNK, count - off);
      const batch = wasm.generate_optimized_conformers_wasm(
        sdfHeavy, c, BigInt(seedBase + off), engine, maxIter);
      if (!batch.get_success()) throw new Error(batch.get_error());
      const nc = batch.get_n_confs();
      if (!template && nc > 0) template = batch.get_template_sdf();
      self.postMessage({ type: 'meta', seedBase, nConfs: off + nc, template });
      const flat = batch.get_coordinates();
      const na = batch.get_n_atoms();
      const energies = batch.get_energies();
      const conv = batch.get_converged();
      const iters = batch.get_iterations();
      const seeds = batch.get_seeds();
      for (let i = 0; i < nc; i++) {
        const coords = Array.from(flat.slice(i * na * 3, (i + 1) * na * 3));
        self.postMessage({
          type: 'conf', seed: seeds[i], E: energies[i],
          converged: !!conv[i], iterations: iters[i],
          termsJson: null,
          sdf: buildSdfFromCoords(coords, template), coords,
        });
      }
    }
    console.log('[worker] embedded+optimized ' + count + ' confs in ' + (Date.now() - t0) + 'ms');
    self.postMessage({ type: 'done', seedBase });
  } catch (err) {
    self.postMessage({ type: 'error', message: String(err && err.message ? err.message : err) });
  }
};
