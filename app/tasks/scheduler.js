// app/tasks/scheduler.js — global concurrency tokens + two priorities +
// queue-level cancellation (docs/task-contract.md §8). Pure logic; Worker
// wiring lives in pools.js, persistence in store.js. Node-testable.
(function (root, factory) {
  if (typeof module === 'object' && module.exports) module.exports = factory();
  else { (root.Tasks = root.Tasks || {}).Scheduler = factory(); }
})(typeof self !== 'undefined' ? self : this, function () {
  'use strict';

  // total = min(4, cores); background may hold at most total-1 (one slot
  // reserved for interactive work, which may fill all of them).
  function createScheduler(opts) {
    const o = opts || {};
    const total = Math.max(1, Math.min(4, o.cores || (typeof navigator !== 'undefined' && navigator.hardwareConcurrency) || 4));
    let backgroundHeld = 0, interactiveHeld = 0;
    const waiters = [];            // {priority, resolve, token} FIFO
    const tasks = new Map();       // taskId -> {generation, priority, state, cancelledAt}
    let seq = 0;

    function held(p) { return p === 'interactive' ? interactiveHeld : backgroundHeld; }
    function free(p) {
      const used = backgroundHeld + interactiveHeld;
      if (p === 'interactive') return used < total;
      return used < total && backgroundHeld < total - 1;   // reserved slot
    }
    function pump() {
      for (let i = 0; i < waiters.length;) {
        if (free(waiters[i].priority)) {
          const w = waiters.splice(i, 1)[0];
          if (w.priority === 'interactive') interactiveHeld++; else backgroundHeld++;   // grant = hold
          let released = false;
          w.resolve(() => {   // the acquire() resolution IS the release handle
            if (released) return;   // double-release guard
            released = true;
            release(w.priority);
          });
        }
        else i++;
      }
    }

    // a worker job (or a logical task's coordinator) holds one token while
    // executing. acquire() resolves when a slot is free for that priority.
    function acquire(priority) {
      const token = ++seq;
      return new Promise(resolve => {
        waiters.push({ priority: priority === 'interactive' ? 'interactive' : 'background', resolve, token, held: false });
        pump();
      });
    }
    function release(priority) {
      if (priority === 'interactive') interactiveHeld = Math.max(0, interactiveHeld - 1);
      else backgroundHeld = Math.max(0, backgroundHeld - 1);
      pump();
    }

    // register a logical task (funnel runs register once; their worker jobs
    // each take tokens through acquire). Returns a handle whose generation
    // guards late replies.
    function startTask(taskId, priority) {
      const rec = tasks.get(taskId) || { generation: 0, priority, state: 'queued', cancelled: false };
      rec.priority = priority === 'interactive' ? 'interactive' : 'background';
      rec.state = 'running'; rec.cancelled = false;
      tasks.set(taskId, rec);
      return {
        taskId,
        generation: rec.generation,
        priority: rec.priority,
        cancelled() { return !!rec.cancelled; },
        // cooperative checkpoint between stages; hard cancel is terminate()
        checkpoint(stageId) {
          if (rec.cancelled) { const e = new Error('task cancelled' + (stageId ? ' at ' + stageId : '')); e.code = 'TASK_CANCELLED'; throw e; }
          return true;
        },
        done(state) { rec.state = state || 'complete'; },
      };
    }

    function cancel(taskId) {
      const rec = tasks.get(taskId);
      if (!rec) return 'unknown';
      if (rec.state !== 'running' && rec.state !== 'queued') return rec.state;
      rec.cancelled = true;
      rec.generation++;            // late replies with the old generation are dropped
      return 'cancel-requested';
    }
    function drop(taskId) { tasks.delete(taskId); }
    function gen(taskId) { const r = tasks.get(taskId); return r ? r.generation : -1; }
    function stats() {
      return {
        total, backgroundHeld, interactiveHeld,
        waiting: waiters.length,
        tasks: [...tasks.entries()].map(([id, r]) => ({ id, state: r.state, priority: r.priority, cancelled: r.cancelled, generation: r.generation })),
      };
    }
    return { acquire, release, startTask, cancel, drop, gen, stats };
  }

  return { createScheduler };
});
