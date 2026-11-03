// Parallel iMTD-GC driver (v1.12): scatters each round's short-metad runs
// across a small conf.worker pool and gathers their quenched minima into the
// main-thread ImtdSearch coordinator. The coordinator's fold is
// arrival-order-independent, so this produces the same bit-identical ensemble
// as the sequential wasm.run_imtd_conformer_search call.
//
// Cancellation: IMTD.start() returns { promise, cancel(reason) }. cancel()
// terminates the pool immediately; the promise RESOLVES (not rejects) with
// cancelled: true + cancelReason and the last fully-folded round's ensemble —
// an honest partial result.
(function () {
  'use strict';

  // Mirror of ImtdConfig::default() (Rust, src/metad/imtd.rs). The driver
  // always sends COMPLETE option sets to workers, so this table is the single
  // JS-side default source — keep it in sync with the Rust default.
  const DEFAULTS = {
    mmff_variant: 'MMFF94s',
    max_rounds: 4,
    runs_per_round: 3,
    steps_per_run: 2000,
    quench_interval: 250,
    dt_fs: 1.0,
    temperature_k: 400.0,
    friction_per_ps: 5.0,
    seed: 42,
    hill_height: 0.5,
    hill_width: 0.2,
    hill_width_distance: 0.5,
    deposit_interval: 25,
    bias_factor: 10.0,
    quench_max_iterations: 300,
    energy_window: 6.0,
    rmsd_tol: 0.25,
    max_conformers: 250,
    crossing_pairs: 4,
    convergence_new_minima: 0,
  };

  // Complete [setter, value] pairs for a merged option set (defaults +
  // overrides). 'workers' is a driver-only knob, not an ImtdOptions field.
  function optSets(opts) {
    const merged = Object.assign({}, DEFAULTS, opts || {});
    const out = [];
    for (const [k, v] of Object.entries(merged)) {
      if (v === undefined || k === 'workers') continue;
      out.push(['set_' + k, v]);
    }
    return out;
  }

  function optsToWasm(opts) {
    const wasm = window.webmm;
    const o = new wasm.ImtdOptions();
    for (const [setter, v] of optSets(opts)) o[setter](v);
    return o;
  }

  // Start the search; returns { promise, cancel(reason) }.
  function start(sdf, opts, onRound) {
    let cancelFlag = null;   // becomes the cancel reason string on cancel()
    const cancelLatch = [];  // resolvers fired on cancel() (unblocks in-flight rounds)
    let search = null;       // the coordinator (assigned by the run below; the
                             // verify batches below reuse it after completion)
    const promise = (async () => {
      const wasm = window.webmm;
      if (!wasm) throw new Error('wasm not ready');
      search = new wasm.ImtdSearch(sdf, optsToWasm(opts));
      if (!search.success()) throw new Error(search.error() || 'ImtdSearch failed');

      // small dedicated pool (farm discipline: hardware-bound, ≤4)
      const NW = Math.max(1, Math.min((opts && opts.workers) || 4, navigator.hardwareConcurrency || 4, 4));
      const wasmUrl = window.__webmmUrl;
      const workers = [];
      const waiters = Object.create(null);   // id → resolve
      const used = new Set();                // worker indices that executed a run
      let nextId = 0;

      const ensure = (w) => new Promise((res) => {
        const iv = setInterval(() => { if (w.__ready) { clearInterval(iv); res(); } }, 20);
      });
      for (let i = 0; i < NW; i++) {
        const w = new Worker('conf.worker.js', { type: 'module' });
        w.__ready = false;
        w.onmessage = (ev) => {
          const m = ev.data;
          if (m.type === 'ready') { w.__ready = true; return; }
          if (m.type === 'imtdMinima' || m.type === 'imtdError') {
            const res = waiters[m.id];
            if (res) { delete waiters[m.id]; res(m); }
          }
        };
        w.postMessage({ type: 'init', wasmUrl });
        workers.push(w);
      }
      const terminateAll = () => { for (const w of workers) w.terminate(); };
      if (cancelFlag !== null) terminateAll();
      else await Promise.all(workers.map(ensure));

      const runOne = (spec, specJson, i) => new Promise((res, rej) => {
        const id = nextId++;
        const wIdx = i % NW;
        used.add(wIdx);
        waiters[id] = (m) => (m.type === 'imtdMinima' ? res(m) : rej(new Error(m.error)));
        workers[wIdx].postMessage({ type: 'imtdRun', id, runIndex: spec.run_index, sdf, specJson, optSets: optSets(opts) });
      });

      const t0 = performance.now();
      let cancelled = false, cancelReason = null;
      try {
        while (cancelFlag === null && !search.done()) {
          const plan = JSON.parse(search.plan_json());
          if (!plan.runs.length) break;
          // scatter the whole round, then gather (arrival order is
          // irrelevant); a cancel() mid-round unblocks the await — the
          // incomplete round's runs are discarded, the ensemble is the last
          // fully-folded round
          const cancelledSignal = new Promise(res => cancelLatch.push(res));
          const outs = await Promise.race([
            Promise.all(plan.runs.map((spec, i) => runOne(spec, JSON.stringify(spec), i))),
            cancelledSignal.then(() => null),
          ]);
          if (outs === null) break;
          for (const m of outs) {
            if (!search.gather(m.runIndex, m.minimaJson)) {
              throw new Error('gather rejected worker minima: ' + search.error());
            }
          }
          if (onRound) onRound({ round: plan.round, runs: plan.runs.length });
        }
        cancelled = cancelFlag !== null;
        cancelReason = cancelFlag;
      } finally {
        terminateAll();
      }

      const res = search.result();
      const landscapeJson = search.landscape_json();
      // (verify is attached to the handle below — the coordinator stays
      // alive after the search so MD-occupancy batches can extend it)
      return {
        landscapeJson,
        nConfs: res.n_confs(),
        energies: res.energies(),
        coordinates: res.coordinates(),
        roundFound: res.round_found(),
        converged: res.converged(),
        roundsExecuted: res.rounds_executed(),
        newMinimaPerRound: res.new_minima_per_round(),
        poolSize: res.pool_size(),
        nRotatable: res.n_rotatable(),
        workers: workers.length,
        usedWorkers: used.size,
        ms: performance.now() - t0,
        cancelled, cancelReason,
      };
    })();
    // Plain-MD verification batch (hill_height=0 specs from the finished
    // coordinator): a fresh temporary pool executes them with the caller's
    // cadence; frames fold into the MD occupancy ladder. Shares cancelFlag.
    const verify = async (nRuns, stepsOverride) => {
      if (cancelFlag !== null) throw new Error('cancelled');
      if (!search) throw new Error('search not started');
      const plan = JSON.parse(search.plan_verify_json(nRuns || 4));
      if (!plan.runs || !plan.runs.length) {
        throw new Error('no verification runs planned (search still running?)');
      }
      const baseSteps = stepsOverride || (opts && opts.steps_per_run) || 2000;
      const vOpts = Object.assign({}, opts || {}, {
        steps_per_run: baseSteps,
        quench_interval: Math.max(10, Math.round(baseSteps / 8)),
      });
      const NWv = Math.max(1, Math.min(plan.runs.length, navigator.hardwareConcurrency || 4, 4));
      const vWorkers = [];
      const vWaiters = Object.create(null);
      let vNext = 0;
      const ensureV = (w) => new Promise((res) => {
        const iv = setInterval(() => { if (w.__ready) { clearInterval(iv); res(); } }, 20);
      });
      for (let i = 0; i < NWv; i++) {
        const w = new Worker('conf.worker.js', { type: 'module' });
        w.__ready = false;
        w.onmessage = (ev) => {
          const m = ev.data;
          if (m.type === 'ready') { w.__ready = true; return; }
          if (m.type === 'imtdMinima' || m.type === 'imtdError') {
            const r = vWaiters[m.id];
            if (r) { delete vWaiters[m.id]; r(m); }
          }
        };
        w.postMessage({ type: 'init', wasmUrl: window.__webmmUrl });
        vWorkers.push(w);
      }
      await Promise.all(vWorkers.map(ensureV));
      const runOneV = (spec, specJson, i) => new Promise((res, rej) => {
        const id = 'v' + (vNext++);
        vWaiters[id] = (m) => (m.type === 'imtdMinima' ? res(m) : rej(new Error(m.error)));
        vWorkers[i % NWv].postMessage({ type: 'imtdRun', id, runIndex: spec.run_index, sdf, specJson, optSets: optSets(vOpts) });
      });
      let cancelled = false;
      try {
        const cancelSig = new Promise(res => cancelLatch.push(res));
        const outs = await Promise.race([
          Promise.all(plan.runs.map((sp, i) => runOneV(sp, JSON.stringify(sp), i))),
          cancelSig.then(() => null),
        ]);
        if (outs === null) {
          cancelled = true;
        } else {
          for (const m of outs) {
            if (!search.gather_verify(m.runIndex, m.minimaJson)) {
              throw new Error('gather_verify rejected worker minima: ' + search.error());
            }
          }
        }
      } finally {
        for (const w of vWorkers) w.terminate();
      }
      return { landscapeJson: search.landscape_json(), cancelled };
    };

    // Saddle refinement batch: one CI-NEB job per observed transition
    // pair, scattered over a fresh temporary pool; results fold as refined
    // saddle energies (stable basin ids). Resolves with the refreshed
    // landscape JSON. MMFF engines only — refused honestly upstream.
    const refine = async () => {
      if (cancelFlag !== null) throw new Error('cancelled');
      if (!search) throw new Error('search not started');
      const plan = JSON.parse(search.plan_refine_json());
      if (!plan.jobs || !plan.jobs.length) {
        throw new Error('no transition pairs to refine (search running?)');
      }
      const NWn = Math.max(1, Math.min(plan.jobs.length, navigator.hardwareConcurrency || 4, 4));
      const nWorkers = [];
      const nWaiters = Object.create(null);
      let nNext = 0;
      const ensureN = (w) => new Promise((res) => {
        const iv = setInterval(() => { if (w.__ready) { clearInterval(iv); res(); } }, 20);
      });
      for (let i = 0; i < NWn; i++) {
        const w = new Worker('conf.worker.js', { type: 'module' });
        w.__ready = false;
        w.onmessage = (ev) => {
          const m = ev.data;
          if (m.type === 'ready') { w.__ready = true; return; }
          if (m.type === 'imtdNebDone' || m.type === 'imtdNebError') {
            const r = nWaiters[m.id];
            if (r) { delete nWaiters[m.id]; r(m); }
          }
        };
        w.postMessage({ type: 'init', wasmUrl: window.__webmmUrl });
        nWorkers.push(w);
      }
      await Promise.all(nWorkers.map(ensureN));
      const runOneN = (job, i) => new Promise((res, rej) => {
        const id = 'n' + (nNext++);
        nWaiters[id] = (m) => (m.type === 'imtdNebDone' ? res(m) : rej(new Error(m.error)));
        nWorkers[i % NWn].postMessage({ type: 'imtdNeb', id, index: job.index, sdf, jobJson: JSON.stringify(job), optSets: optSets(opts) });
      });
      let cancelled = false;
      try {
        const cancelSig = new Promise(res => cancelLatch.push(res));
        const outs = await Promise.race([
          Promise.all(plan.jobs.map((job, i) => runOneN(job, i))),
          cancelSig.then(() => null),
        ]);
        if (outs === null) {
          cancelled = true;
        } else {
          for (const m of outs) {
            if (!search.gather_refine(m.index, m.resultJson)) {
              throw new Error('gather_refine rejected result: ' + search.error());
            }
          }
        }
      } finally {
        for (const w of nWorkers) w.terminate();
      }
      return { landscapeJson: search.landscape_json(), cancelled, jobs: plan.jobs.length };
    };

    return {
      promise,
      verify,
      refine,
      // pass-through to the live coordinator (thermo is a view derivative,
      // not a batch)
      thermoJson: (t) => search.thermo_json(t),
      cancel(reason) {
        if (cancelFlag !== null) return;
        cancelFlag = reason === undefined ? 'user cancel' : String(reason);
        for (const res of cancelLatch.splice(0)) res();
      },
    };
  }

  // Cancellation-less sugar over start().
  async function runParallel(sdf, opts, onRound) {
    return start(sdf, opts, onRound).promise;
  }

  window.IMTD = { start, runParallel, optsToWasm, DEFAULTS };
})();
