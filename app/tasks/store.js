// app/tasks/store.js — task side table + Done transaction + reload
// reconciliation (docs/task-contract.md §6). Binds to ONE platform storage
// (per-project DB); rebind on project switch. Writes at exactly three
// points: submit (queued), state transition, Done (inside the command+facts
// transaction via applyWithFacts extraOps). Progress never touches disk.
(function (root, factory) {
  if (typeof module === 'object' && module.exports) module.exports = factory();
  else { (root.Tasks = root.Tasks || {}).Store = factory(); }
})(typeof self !== 'undefined' ? self : this, function () {
  'use strict';

  function C() {
    if (typeof module === 'object' && module.exports) return require('./contract.js');
    return (typeof self !== 'undefined' ? self : this).Tasks.Contract;
  }

  const SESSION_ID = 's' + Date.now().toString(36) + '-' + Math.random().toString(36).slice(2, 8);

  function createStore(storage) {
    // records: taskId -> { spec, state, createdAt?, startedAt?, endedAt?,
    //                      durationMs?, failure?, attemptId?, generation? }
    async function get(taskId) { return storage.get('tasks', taskId); }
    async function all() { return storage.all('tasks'); }

    async function submit(spec) {
      const prev = await get(spec.taskId);
      if (prev && prev.state !== 'failed' && prev.state !== 'cancelled' && prev.state !== 'abandoned') {
        throw new Error('task ' + spec.taskId + ' already exists in state ' + prev.state);
      }
      const rec = prev || {};
      rec.spec = spec; rec.state = 'queued';
      rec.submittedAt = Date.now(); rec.attemptId = (rec.attemptId || 0) + 1; rec.generation = rec.generation || 0;
      await storage.put('tasks', spec.taskId, rec);
      return rec;
    }

    async function transition(taskId, to, patch) {
      const rec = await get(taskId);
      if (!rec) throw new Error('transition: unknown task ' + taskId);
      const Contract = C();
      if (!Contract.canTransition(rec.state, to)) {
        throw new Error('transition ' + rec.state + ' -> ' + to + ' forbidden (task ' + taskId + ')');
      }
      rec.state = to;
      Object.assign(rec, patch || {});
      await storage.put('tasks', taskId, rec);
      return rec;
    }

    // Done: the tasks-complete write rides the SAME cross-store transaction
    // as the command + facts (applyWithFacts extraOps). Idempotent: a repeat
    // Done for an already-complete task is a no-op, a stale generation is
    // rejected. Returns 'committed' | 'already' | throws.
    async function commitDone(opts) {
      const { taskId, cmd, factPuts, applyWithFacts, generation } = opts;
      const rec = await get(taskId);
      if (!rec) throw new Error('commitDone: unknown task ' + taskId);
      if (rec.state === 'complete') return 'already';
      if (rec.state !== 'running') throw new Error('commitDone: task in state ' + rec.state + ' (expected running)');
      if (generation != null && rec.generation != null && generation !== rec.generation) {
        throw new Error('commitDone: stale generation ' + generation + ' != ' + rec.generation);
      }
      const doneRec = { ...rec, state: 'complete', endedAt: Date.now() };
      await applyWithFacts(cmd, factPuts, [{ store: 'tasks', key: taskId, value: doneRec }]);
      // applyWithFacts landed the record; mirror locally for immediate reads
      return 'committed';
    }

    // Reload reconciliation (§6): queued/running tasks with their Done fact
    // already present are ghost-completes; the rest are abandoned unless a
    // live session claims them (BroadcastChannel heartbeat, 400ms window).
    // liveIds: taskIds active in THIS session (scheduler-owned).
    async function reconcile(opts) {
      const { factKeyFor, liveIds, claimWindowMs } = opts || {};
      const live = new Set(liveIds || []);
      const rows = await all();
      const orphans = rows.filter(r => (r.value.state === 'queued' || r.value.state === 'running') && !live.has(r.value.spec.taskId));
      const out = { complete: [], abandoned: [] };
      if (!orphans.length) return out;

      // ask other tabs to claim their live tasks
      const claimed = new Set();
      let ch = null;
      if (typeof BroadcastChannel !== 'undefined') {
        ch = new BroadcastChannel('webmm-tasks');
        ch.onmessage = ev => {
          const m = ev.data || {};
          if (m.type === 'claimed' && m.taskId) claimed.add(m.taskId);
        };
        ch.postMessage({ type: 'claim', taskIds: orphans.map(r => r.value.spec.taskId), sessionId: SESSION_ID });
        await new Promise(r => setTimeout(r, claimWindowMs != null ? claimWindowMs : 400));
        try { ch.close(); } catch (e) {}
      }

      for (const r of orphans) {
        const taskId = r.value.spec.taskId;
        if (claimed.has(taskId)) continue;   // another live tab owns it
        const fk = factKeyFor && factKeyFor(r.value.spec);
        const fact = fk ? await storage.get('facts', fk) : null;
        if (fact) {
          // ghost complete: the Done transaction landed, the side-table
          // write raced with the reload
          const rec = { ...r.value, state: 'complete', endedAt: r.value.endedAt || Date.now() };
          await storage.put('tasks', taskId, rec);
          out.complete.push(taskId);
        } else {
          await storage.put('tasks', taskId, { ...r.value, state: 'abandoned', endedAt: Date.now() });
          out.abandoned.push(taskId);
        }
      }
      return out;
    }

    // claim responder: call with this session's live task ids whenever a
    // 'claim' broadcast arrives (wired once by the page shell)
    function attachClaimResponder(getLiveIds) {
      if (typeof BroadcastChannel === 'undefined') return null;
      const ch = new BroadcastChannel('webmm-tasks');
      ch.onmessage = ev => {
        const m = ev.data || {};
        if (m.type === 'claim' && m.sessionId !== SESSION_ID && Array.isArray(m.taskIds)) {
          const mine = (getLiveIds() || []).filter(id => m.taskIds.includes(id));
          for (const id of mine) { try { ch.postMessage({ type: 'claimed', taskId: id, sessionId: SESSION_ID }); } catch (e) {} }
        }
      };
      return ch;
    }

    return { SESSION_ID, get, all, submit, transition, commitDone, reconcile, attachClaimResponder };
  }

  return { createStore, SESSION_ID };
});
