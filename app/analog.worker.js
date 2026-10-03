// Analog-explorer farm worker (M3): for one candidate (molblock + nConfs +
// seed), generate its conformer ensemble and rigid-align every conformer
// against the parent's reference conformer (shape-only in the worker —
// color rescoring of winners happens main-thread where RDKit lives).
// Messages chain through one promise queue (same serialization lesson as
// shape.worker.js).
let wasm = null;
let initPromise = null;
let queue = Promise.resolve();

self.onmessage = (e) => {
  const msg = e.data;
  queue = queue.then(() => handle(msg)).catch(err => {
    self.postMessage({ type: 'done', idx: msg.idx, failed: String(err && err.message || err).slice(0, 120) });
  });
};

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
  const tail = lines.slice(4 + na + nb).filter(l => l.startsWith('M  ') || l.startsWith('A  ') || l.trim() === '');
  return header.concat(atomLines, bondLines, tail).join('\n');
}

async function handle(msg) {
  if (msg.type === 'init') {
    if (!initPromise) initPromise = (async () => {
      wasm = await import(msg.wasmUrl);
      await wasm.default();
    })();
    await initPromise;
    self.postMessage({ type: 'ready' });
    return;
  }
  if (msg.type === 'score') {
    const r = wasm.generate_optimized_conformers_wasm(msg.molblock, msg.nConfs, BigInt(msg.seed), 'MMFF94s', 100);
    if (!r.get_success() || !r.get_n_confs()) { self.postMessage({ type: 'done', idx: msg.idx, failed: 'conformers' }); return; }
    const template = r.get_template_sdf();
    const coords = r.get_coordinates();
    const nAtoms = parseInt(sdfLines(template)[3].substring(0, 3));
    let best = null, bestSdf = null;
    for (let j = 0; j < r.get_n_confs(); j++) {
      const csdf = buildSdfFromCoords(Array.from(coords.slice(j * nAtoms * 3, (j + 1) * nAtoms * 3)), template);
      let res = null;
      try { res = JSON.parse(wasm.shape_align_wasm(msg.parentSdf, csdf, '')); } catch (e) { continue; }
      if (res && (best === null || res.tanimoto > best.tanimoto)) { best = res; bestSdf = csdf; }
    }
    // stage D (flex) needs the winning pose: hand back the raw conformer
    // SDF and its align transform when asked (backward compatible without)
    self.postMessage({
      type: 'done', idx: msg.idx, tanimoto: best ? best.tanimoto : null,
      nConfs: r.get_n_confs(),
      bestSdf: msg.returnBest ? bestSdf : undefined,
      transform: msg.returnBest && best ? best.transform : undefined,
    });
  }
}
