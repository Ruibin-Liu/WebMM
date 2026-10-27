// app/tasks/funnel.js — generic staged funnel runner (docs/task-contract.md
// §5/§9.1). Barrier semantics: a stage's top-n selection happens ONLY after
// ALL of its candidates completed — streaming feeds the UI, never the
// selection. Deterministic tie-break is contract: (-score, smiles, idx).
// Pure orchestration: executors are injected, Node-testable.
(function (root, factory) {
  if (typeof module === 'object' && module.exports) module.exports = factory();
  else { (root.Tasks = root.Tasks || {}).Funnel = factory(); }
})(typeof self !== 'undefined' ? self : this, function () {
  'use strict';

  function emptyAccount() { return { evaluated: 0, filtered: 0, failed: 0, notEvaluated: 0 }; }

  // stable rank with contract tie-break: score desc, smiles asc, idx asc —
  // identical inputs give identical order at ANY concurrency.
  function rank(entries) {
    return entries.slice().sort((a, b) => {
      if (b.score !== a.score) return b.score - a.score;
      if (a.smiles !== b.smiles) return a.smiles < b.smiles ? -1 : 1;
      return a.idx - b.idx;
    });
  }

  // combine 'max': strict '>' keeps the earlier (B) value on exact ties —
  // deterministic regardless of which ensemble replies first.
  function combineMax(prev, next) { return next > prev ? next : prev; }

  // run one score stage over `list`: dispatches min(list, cap) jobs through
  // opts.submit, barrier-collects ALL replies into a Map(idx -> msg), then
  // returns it. Deterministic-budget cut only: a job that never got a slot
  // is simply absent from the map (counted notEvaluated by the caller).
  async function scoreJobs(list, stage, acc, ctx, opts) {
    const budget = Math.min(list.length, stage.cap);
    let next = 0, done = 0, stop = false;
    const msgs = new Map();

    async function lane() {
      while (!stop) {
        const i = next++;
        if (i >= budget) return;
        ctx.checkpoint(stage.id);
        const e = list[i];
        try {
          const msg = await opts.submit(e, stage);
          done++;
          if (ctx.cancelled()) { stop = true; return; }
          if (msg && msg.score != null) { msgs.set(e.idx, msg); acc[stage.id].evaluated++; }
          else acc[stage.id].failed++;
          if (opts.onProgress) opts.onProgress({ stage: stage.id, done, total: budget, scored: msgs.size });
        } catch (err) {
          done++;
          if (err && err.code === 'TASK_CANCELLED') { stop = true; return; }
          acc[stage.id].failed++;
        }
      }
    }
    const lanes = Math.max(1, Math.min(opts.concurrency() * (stage.inflightPerWorker || 2), budget));
    await Promise.all(Array.from({ length: lanes }, lane));
    if (stop || ctx.cancelled()) { const e = new Error('funnel interrupted'); e.code = 'TASK_INTERRUPTED'; throw e; }
    return msgs;
  }

  // spec: makeSpec() output for analogScan (input.stages = effective snapshot)
  // ctx: scheduler task handle (checkpoint/cancelled)
  // opts: { submit(entry, stage) -> Promise<{score,bestSdf?}>, flex(entry) -> Promise<number|null>,
  //         concurrency() -> int, onProgress, onSettled(stageId, rankedEntries) }
  async function runAnalogFunnel(spec, ctx, opts) {
    const stages = spec.input.stages;
    const acc = { A: emptyAccount(), B: emptyAccount(), C: emptyAccount(), D: emptyAccount() };
    for (const k of Object.keys(acc)) acc[k].stage = k;
    // stage A is the audit snapshot of the exempt (main-thread) enumeration
    acc.A.evaluated = spec.input.candidates.length;
    acc.A.filtered = (spec.input.audit && spec.input.audit.filtered) || 0;
    acc.A.notEvaluated = (spec.input.audit && spec.input.audit.notEvaluated) || 0;

    const entries = spec.input.candidates.map(c => ({ ...c, score: null, flexScore: null, bestSdf: null }));
    const stageB = stages.find(s => s.id === 'B' && s.op === 'score');
    const stageC = stages.find(s => s.id === 'C' && s.op === 'score');
    const stageD = stages.find(s => s.id === 'D' && s.op === 'flex');

    if (stageB) {
      const msgs = await scoreJobs(entries, stageB, acc, ctx, opts);
      for (const e of entries) {
        const m = msgs.get(e.idx);
        if (m) { e.score = m.score; if (m.bestSdf) e.bestSdf = m.bestSdf; }
      }
      for (const e of entries) if (e.score == null) { /* failed or notEvaluated already counted */ }
      if (opts.onSettled) opts.onSettled('B', rank(entries.filter(e => e.score != null)));
    }

    if (stageC) {
      ctx.checkpoint('C');
      // deterministic top-n AFTER the B barrier — never "first n to arrive"
      const subset = rank(entries.filter(e => e.score != null)).slice(0, stageC.cap);
      const msgs = await scoreJobs(subset, stageC, acc, ctx, opts);
      for (const e of subset) {
        const m = msgs.get(e.idx);
        if (m) {
          e.score = e.score != null && stageC.combine === 'max' ? combineMax(e.score, m.score) : m.score;
          if (m.bestSdf) e.bestSdf = m.bestSdf;
        }
      }
      if (opts.onSettled) opts.onSettled('C', rank(entries.filter(e => e.score != null)));
    }

    if (stageD) {
      ctx.checkpoint('D');
      const top = rank(entries.filter(e => e.score != null)).slice(0, stageD.cap);
      let d = 0;
      for (const e of top) {
        ctx.checkpoint('D');
        if (!e.bestSdf || (e.score != null && e.score >= stageD.skipSelfT)) { acc.D.filtered++; continue; }
        d++;
        try { e.flexScore = await opts.flex(e); acc.D.evaluated++; }
        catch (err) {
          if (err && err.code === 'TASK_CANCELLED') { const i = new Error('funnel interrupted'); i.code = 'TASK_INTERRUPTED'; throw i; }
          acc.D.failed++;
        }
        if (opts.onProgress) opts.onProgress({ stage: 'D', done: d, total: top.length, scored: acc.D.evaluated });
      }
      if (opts.onSettled) opts.onSettled('D', rank(entries.filter(e => e.score != null)));
    }

    // notEvaluated (§5): candidates the deterministic budget never reached —
    // counted against the stage's FULL domain (B: all candidates; C: the
    // post-barrier top-n subset), so a cap-cut shows up as notEvaluated
    // rather than vanishing.
    if (stageB) {
      acc.B.notEvaluated = Math.max(0, entries.length - acc.B.evaluated - acc.B.failed);
    }
    if (stageC) {
      const subsetN = Math.min(rank(entries.filter(e => e.score != null)).length, stageC.cap);
      acc.C.notEvaluated = Math.max(0, subsetN - acc.C.evaluated - acc.C.failed);
    }
    if (stageD) {
      acc.D.notEvaluated = 0;   // D walks its whole (small) list deterministically
    }

    const final = rank(entries.filter(e => e.score != null));
    return { entries: final, accounts: acc, interrupted: false };
  }

  return { runAnalogFunnel, rank, combineMax, emptyAccount };
});
