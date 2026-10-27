// app/tasks/pools.js — per-module worker pools under GLOBAL scheduler
// tokens (docs/task-contract.md §8). Generalizes the M3 analog pool:
// + each in-flight job acquires a token (background/interactive budget is
//   enforced across ALL pools, not per pool)
// + jobs-per-worker recycling (idle terminate + respawn; fresh RNG/heap)
// + one warm spare beyond pool size when saturated (hard-cancel recovery)
// Same poison/redispatch/heap logic as the M3 pool. stats() is a superset
// of the historical __analogPool shape (m6 contract).
// Browser-only (Worker); Node tests mock the Worker via opts.workerCtor.
(function (root, factory) {
  if (typeof module === 'object' && module.exports) module.exports = factory();
  else { (root.Tasks = root.Tasks || {}).Pools = factory(); }
})(typeof self !== 'undefined' ? self : this, function () {
  'use strict';

  function createPool(opts) {
    const {
      url, initPayload, scheduler,
      workerCtor,               // injectable for tests; defaults to Worker
      jobsPerRecycle = 200,
      heapLimit = 256 * 1024 * 1024,
      heapProbeEvery = 25,
    } = opts;
    const W = workerCtor || (typeof Worker !== 'undefined' ? Worker : null);
    if (!W) throw new Error('pools: no Worker constructor (Node needs workerCtor)');

    const recs = [];
    const queue = [];            // {id, payload, priority, resolve, tries}
    const inflight = new Map();  // id -> entry
    let seq = 0, generation = 0, poisons = 0, heapRecycles = 0, jobRecycles = 0, initFails = 0;
    let initFailed = false, readyPromise = null, doneSinceProbe = 0;

    function poolSize() {
      const cores = (typeof navigator !== 'undefined' && navigator.hardwareConcurrency) || 4;
      let n = Math.max(1, Math.min(cores - 1, 4));
      const dm = typeof navigator !== 'undefined' && navigator.deviceMemory;
      if (dm && dm <= 2) n = Math.max(1, n >> 1);
      return n;
    }
    const baseSize = opts.size || poolSize();

    async function drain() {
      // dispatch while a ready idle worker exists AND the scheduler grants a
      // token for this job's priority; spawn one warm spare if saturated.
      // acquire() is async — the world can move under the await, so the
      // grant is re-checked against queue/worker state and RELEASED when it
      // can't be used (a leaked token would wedge the whole background
      // budget — the exact bug the poison-chaos test caught).
      for (;;) {
        const job = queue[0];
        if (!job) return;
        const idle = recs.find(r => !r.dead && r.ready && !r.busy);
        if (!idle) {
          if (recs.length < baseSize + 1 && !initFailed) { spawn(); continue; }
          return;
        }
        let release;
        try { release = await scheduler.acquire(job.priority); }
        catch (e) { queue.shift(); job.resolve({ failed: 'token' }); continue; }
        const taken = queue.shift();
        if (!taken) { release(); return; }   // another drain instance dispatched it
        const idle2 = recs.find(r => !r.dead && r.ready && !r.busy);
        if (!idle2) { queue.unshift(taken); release(); return; }   // worker vanished mid-grant
        taken._release = release;
        idle2.busy = true; idle2.job = taken.id; idle2.jobs++;
        inflight.set(taken.id, taken);
        idle2.w.postMessage({ ...taken.payload, __job: taken.id });
      }
    }

    function spawn() {
      generation++;
      const w = new W(url, { type: 'module' });
      const rec = { w, ready: false, busy: false, job: null, heap: 0, dead: false, jobs: 0 };
      w.onmessage = ev => {
        const m = ev.data;
        if (m.type === 'ready') { rec.ready = true; rec.heap = m.wasmHeap || 0; drain(); return; }
        if (m.type === 'heap') {
          rec.heap = m.bytes || 0;
          if (rec.heap > heapLimit && !rec.busy && !rec.dead) { heapRecycles++; recycle(rec); }
          return;
        }
        if (m.type === 'poison') {
          poisons++;
          if (!rec.ready) initFails++;
          if (initFails >= 3) initFailed = true;
          const jid = m.__job != null ? m.__job : rec.job;
          recycle(rec);
          if (jid != null) redispatch(jid);
          return;
        }
        const job = inflight.get(m.__job);
        if (job) {
          inflight.delete(m.__job);
          rec.busy = false; rec.job = null;
          try { job._release && job._release(job.priority); } catch (e) {}
          if (++doneSinceProbe >= heapProbeEvery) { doneSinceProbe = 0; probeHeap(); }
          maybeRetire(rec);          // jobs-per-recycle / spare retirement
          job.resolve(m);
          drain();
        }
      };
      w.onerror = () => {
        poisons++;
        const jid = rec.job;
        recycle(rec);
        if (jid != null) redispatch(jid);
      };
      w.postMessage(initPayload);
      recs.push(rec);
      return rec;
    }

    // idle-time retirement: a worker past its job budget, or a spare when
    // the pool is no longer saturated, terminates and (budget case only)
    // respawns fresh
    function maybeRetire(rec) {
      if (rec.dead || rec.busy) return;
      const saturated = queue.length > 0 && recs.filter(r => !r.dead && r.busy).length >= baseSize;
      if (rec.jobs >= jobsPerRecycle) { jobRecycles++; recycle(rec, true); return; }
      if (!saturated && recs.filter(r => !r.dead).length > baseSize) { retireOne(); }
    }
    function retireOne() {
      const idle = recs.filter(r => !r.dead && !r.busy);
      if (idle.length <= 0) return;
      const spare = idle[idle.length - 1];
      spare.dead = true;
      try { spare.w.terminate(); } catch (e) {}
      const i = recs.indexOf(spare); if (i >= 0) recs.splice(i, 1);
    }

    function recycle(rec, respawn) {
      if (rec.dead) return;
      rec.dead = true;
      try { rec.w.terminate(); } catch (e) {}
      const i = recs.indexOf(rec);
      if (i >= 0) recs.splice(i, 1);
      if (respawn || !initFailed) { if (recs.filter(r => !r.dead).length < baseSize + 1) spawn(); }
      drain();
    }

    function redispatch(jobId) {
      const job = inflight.get(jobId);
      if (!job) return;
      inflight.delete(jobId);
      if (job._release) { try { job._release(job.priority); } catch (e) {} }
      if (job.tries < 1) { job.tries++; queue.push(job); drain(); }
      else job.resolve({ failed: 'poisoned twice' });
    }

    function probeHeap() {
      for (const rec of recs) if (rec.ready && !rec.busy && !rec.dead) { try { rec.w.postMessage({ type: 'heap' }); } catch (e) {} }
    }

    function ready() {
      if (!readyPromise) readyPromise = new Promise((res, rej) => {
        initFailed = false; initFails = 0;
        while (recs.length < baseSize) spawn();
        const t0 = (typeof performance !== 'undefined' ? performance.now() : Date.now());
        const iv = setInterval(() => {
          if (initFailed) {
            clearInterval(iv); readyPromise = null;
            for (const r of [...recs]) { r.dead = true; try { r.w.terminate(); } catch (e) {} }
            recs.length = 0;
            rej(new Error('init')); return;
          }
          if (recs.filter(r => r.ready).length >= baseSize) { clearInterval(iv); res(baseSize); return; }
          if (((typeof performance !== 'undefined' ? performance.now() : Date.now()) - t0) > 30000) {
            clearInterval(iv); readyPromise = null; initFailed = true;
            for (const r of [...recs]) { r.dead = true; try { r.w.terminate(); } catch (e) {} }
            recs.length = 0;
            rej(new Error('timeout')); return;
          }
        }, 50);
      });
      return readyPromise;
    }

    function submit(payload, priority) {
      return new Promise(resolve => {
        queue.push({ id: ++seq, payload, priority: priority === 'interactive' ? 'interactive' : 'background', resolve, tries: 0 });
        drain();
      });
    }

    function stats() {
      return {
        size: recs.length, ready: recs.filter(r => r.ready).length, queued: queue.length,
        generation, poisons, heapRecycles, jobRecycles,
        heaps: recs.map(r => r.heap), jobs: recs.map(r => r.jobs),
      };
    }
    return { ready, submit, stats, baseSize };
  }

  function createPools(scheduler, workerCtor) {
    const pools = new Map();
    function get(name, opts) {
      if (!pools.has(name)) pools.set(name, createPool({ ...opts, scheduler, workerCtor }));
      return pools.get(name);
    }
    function has(name) { return pools.has(name); }
    return { get, has, stats() { const o = {}; for (const [k, p] of pools) o[k] = p.stats(); return o; } };
  }

  return { createPool, createPools };
});
