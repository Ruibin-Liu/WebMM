// Shape-search farm worker (ES module): generates an entry's conformer
// ensemble and screens every conformer against the query (cheap proxy),
// and full-aligns phase-2 winner candidates. Parameter strings are kept
// VERBATIM identical to the main-thread sequential implementation so the
// results are deterministic across both paths.
//
// Self-contained copies of sdfLines/buildSdfFromCoords (same convention as
// conf.worker.js; the main thread keeps its own). The wasm module URL is
// layout-dependent — the main thread resolves it and hands it over in the
// 'prep' message (first of a search session).
let wasm = null;
let initPromise = null;
// per-search-session query state (set by 'prep', used by 'full')
let qsdf = null;
let qSitesJson = null;

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
  // preserve the M-tail (formal charges! — see conf.worker.js)
  const tail = lines.slice(4 + na + nb);
  return [...header, ...atomLines, ...bondLines, ...(tail.length ? tail : ['M  END'])].join('\n');
}

// async onmessage handlers do NOT queue: a 'full' arriving while an
// 'initq'/'prep' is still awaiting its wasm import would see wasm === null
// and fail. Chain every message through one promise queue (serialized).
let queue = Promise.resolve();
self.onmessage = (e) => {
  const msg = e.data;
  queue = queue.then(() => handle(msg)).catch(err => {
    self.postMessage({ type: 'error', message: String(err && err.message ? err.message : err), idx: msg.idx });
  });
};

async function handle(msg) {
  {
    if (msg.type === 'initq') {
      // session bootstrap (query state + wasm) — sent at worker creation so
      // fully-cached searches (no 'prep' jobs) can still run phase-2 aligns
      if (!initPromise) initPromise = (async () => {
        wasm = await import(msg.wasmUrl);
        await wasm.default();
      })();
      await initPromise;
      qsdf = msg.qsdf;
      qSitesJson = msg.qSitesJson;
      return;
    }
    if (msg.type === 'screen') {
      // cached conformer ensemble, new query: re-screen only (no generation)
      if (!initPromise) initPromise = (async () => {
        wasm = await import(msg.wasmUrl);
        await wasm.default();
      })();
      await initPromise;
      qsdf = msg.qsdf;
      qSitesJson = msg.qSitesJson;
      const proxies = msg.confSdfs.map(sdf =>
        JSON.parse(wasm.shape_align_wasm(qsdf, sdf, '{"screen":true}')).surrogate_overlap);
      self.postMessage({ type: 'prepped', idx: msg.idx, confSdfs: msg.confSdfs, proxies });
      return;
    }
    if (msg.type === 'prep') {
      if (!initPromise) initPromise = (async () => {
        wasm = await import(msg.wasmUrl);
        await wasm.default();
      })();
      await initPromise;
      qsdf = msg.qsdf;
      qSitesJson = msg.qSitesJson;
      const { idx, mb, n, seed, inject } = msg;
      const r = wasm.generate_optimized_conformers_wasm(mb, n, BigInt(seed), 'MMFF94s', 100);
      if (!r.get_success() || !r.get_n_confs()) throw new Error('conformer failed');
      const flat = r.get_coordinates();
      const na = r.get_n_atoms();
      const tpl = r.get_template_sdf();
      const confSdfs = [];
      const proxies = [];
      for (let i = 0; i < r.get_n_confs(); i++) {
        const sdf = i === 0 && inject ? qsdf
          : buildSdfFromCoords(Array.from(flat.slice(i * na * 3, (i + 1) * na * 3)), tpl);
        confSdfs.push(sdf);
        const res = JSON.parse(wasm.shape_align_wasm(qsdf, sdf, '{"screen":true}'));
        proxies.push(res.surrogate_overlap);
      }
      self.postMessage({ type: 'prepped', idx, confSdfs, proxies });
      return;
    }
    if (msg.type === 'full') {
      const { idx, ci, esdf, tSitesJson, useColor } = msg;
      try {
        if (!initPromise) throw new Error('worker not initialized');
        const res = JSON.parse(useColor
          ? wasm.shape_align_color_wasm(qsdf, esdf, qSitesJson, tSitesJson, '{"random_starts": 8}')
          : wasm.shape_align_wasm(qsdf, esdf, '{"random_starts": 8}'));
        self.postMessage({
          type: 'fulled', idx, ci,
          tanimoto: res.tanimoto, color_tanimoto: res.color_tanimoto, transform: res.transform,
        });
      } catch (err) {
        // sentinel: resolves the pending promise on the main thread; the
        // entry simply loses this candidate (skip semantics, like phase 1)
        self.postMessage({ type: 'fulled', idx, ci, failed: true, message: String(err) });
      }
      return;
    }
  }
}
