// platform/project.js — M1a the store API (M0 spec §4)
// Command log is the source of truth; state is its projection; snapshots
// are derived artifacts taken off the interaction path. Undo appends an
// inverse command — the log is never truncated.
(function (root, factory) {
  if (typeof module === 'object' && module.exports) module.exports = factory(require('./state.js'));
  else { (root.Platform = root.Platform || {}).Project = factory((root.Platform || {}).State); }
})(typeof self !== 'undefined' ? self : this, function (State) {
  'use strict';

  const SNAPSHOT_EVERY = 500;       // hybrid trigger: command count (spec §4)
  const KEEP_SNAPSHOTS = 3;

  async function init(storage, opts) {
    const meta = (await storage.get('meta', 'project')) || { schemaVersion: State.SCHEMA_VERSION, createdAt: Date.now() };
    if (meta.schemaVersion > State.SCHEMA_VERSION) {
      throw new Error('project schema v' + meta.schemaVersion + ' newer than this build (v' + State.SCHEMA_VERSION + ') — update the app');
    }
    // load: latest snapshot (if any) + replay the command tail
    const snaps = (await storage.all('snapshots')).sort((a, b) => b.key - a.key);
    let state = State.initialState();
    let tailStart = 0;
    if (snaps.length) {
      state = deserializeState(snaps[0].value.state);
      tailStart = snaps[0].value.commandId;
    }
    const cmdRecs = (await storage.all('commands'))
      .map(r => ({ key: Number(r.key), value: r.value }))
      .sort((a, b) => a.key - b.key);
    const tail = cmdRecs.filter(r => r.key > tailStart);
    for (const r of tail) {
      const res = State.applyCommand(state, r.value);
      if (res.error) { console.warn('[platform] replay skipped illegal command #' + r.key + ':', res.error); continue; }
      state = res.state;
    }
    const seqRec = await storage.get('meta', 'cmdSeq');
    const lastSeq = seqRec || (cmdRecs.length ? cmdRecs[cmdRecs.length - 1].key : 0);
    const p = {
      storage, meta, state,
      withLock: opts && opts.withLock ? opts.withLock : (async (n, fn) => fn()),
      channel: opts && opts.channel ? opts.channel : null,
      _commandsSinceSnapshot: tail.length,
      _undoStack: [],           // {inverse, describes} of undoable commands
      _redoStack: [],           // undone commands, redoable in order
      _cmdId: lastSeq,
      _hadData: cmdRecs.length > 0 || snaps.length > 0,
    };
    try { await storage.put('meta', 'project', meta); } catch (e) {}
    return api(p);
  }

  function api(p) {
    return {
      get state() { return p.state; },
      get schemaVersion() { return p.state.schemaVersion; },

      // apply a user command: validate -> reduce -> persist -> snapshot.
      // A new user action clears the redo stack (standard undo/redo model).
      async apply(cmd) {
        const r = await applyInternal(p, cmd, { pushUndo: true });
        p._redoStack = [];
        return r;
      },

      // Undo = append the inverse command to the LOG (never truncate), but
      // the undo STACK bookkeeping is separate: the inverse application must
      // not push itself back (that would ping-pong forever). The undone
      // command moves to the redo stack.
      async undo() {
        while (p._undoStack.length) {
          const top = p._undoStack.pop();
          if (State.hasDependents(top.describes, p.state)) {
            p._undoStack.push(top);      // put it back — refusal is not consumption
            return { ok: false, reason: 'has-dependents', describes: top.describes };
          }
          await applyInternal(p, top.inverse, { pushUndo: false });
          p._redoStack.push(top.describes);
          return { ok: true, undone: top.describes };
        }
        return { ok: false, reason: 'empty' };
      },

      async redo() {
        const cmd = p._redoStack.pop();
        if (!cmd) return { ok: false, reason: 'empty' };
        await applyInternal(p, cmd, { pushUndo: true });
        return { ok: true, redone: cmd };
      },

      // checkpoint events (import/reconcile/round complete) also trigger
      async checkpoint() { await snapshot(p); },

      // export = facts only (commands + meta). Golden fixtures + CI use
      // this; the M2 UI export is a wrapper.
      async exportProject() {
        const cmds = (await p.storage.all('commands'))
          .map(r => ({ key: Number(r.key), value: r.value }))
          .sort((a, b) => a.key - b.key);
        return { schemaVersion: State.SCHEMA_VERSION, meta: p.meta, commands: cmds.map(c => c.value) };
      },

      // snapshot: the projection itself, tagged with the command id; keep 3
      // (spec §4). Runs on the caller's context (a worker in M1c); never
      // inline on a hot interaction path beyond M1a's trivial sizes.
      snapshot: () => snapshot(p),
      storagePut: (store, key, val) => p.storage.put(store, key, val),
      storageGet: (store, key) => p.storage.get(store, key),
    };
  }

  // shared application path: prior-value capture for Unpin/Include inverses,
  // reduce, persist to the command log, optional undo-stack push.
  // M2b: the append is SERIALIZED under the Web Lock and takes its id from
  // the STORE's cmdSeq (not the local counter) — two tabs with divergent
  // in-memory states can no longer collide on command ids or clobber the
  // log tail. Read-only tabs learn about the change via the broadcast.
  async function applyInternal(p, cmd, opts) {
    let pre = null;
    if (cmd.type === 'Unpin' && p.state.pins[cmd.molId]) {
      pre = { type: 'Pin', molId: cmd.molId, ...p.state.pins[cmd.molId] };
    } else if (cmd.type === 'Include' && p.state.excludes[cmd.molId]) {
      pre = { type: 'Exclude', molId: cmd.molId, ...p.state.excludes[cmd.molId] };
    }
    const res = State.applyCommand(p.state, cmd);
    if (res.error) throw new Error('command rejected: ' + res.error);
    p.state = res.state;
    const id = await p.withLock('project', async () => {
      const seqRec = await p.storage.get('meta', 'cmdSeq');
      const seq = (seqRec || 0) + 1;
      await p.storage.put('commands', seq, cmd);
      await p.storage.put('meta', 'cmdSeq', seq);
      if (!p._hadData) { p._hadData = true; try { await p.storage.put('meta', 'hadData', true); } catch (e) {} }
      return seq;
    });
    p._cmdId = id;
    try { if (p.channel) p.channel.postMessage({ kind: 'applied', commandId: id }); } catch (e) {}
    if (opts && opts.pushUndo) {
      const inv = pre || State.inverse(cmd, res.state);
      if (inv) p._undoStack.push({ inverse: inv, describes: cmd });
    }
    p._commandsSinceSnapshot += 1;
    if (p._commandsSinceSnapshot >= SNAPSHOT_EVERY) await snapshot(p);
    return { commandId: id };
  }

  async function snapshot(p) {
    const id = p._cmdId;
    await p.storage.put('snapshots', id, { commandId: id, state: serializeState(p.state), at: Date.now() });
    const snaps = (await p.storage.all('snapshots')).sort((a, b) => b.key - a.key);
    for (const s of snaps.slice(KEEP_SNAPSHOTS)) await p.storage.del('snapshots', s.key);
    p._commandsSinceSnapshot = 0;
  }

  // Maps serialize as entries pairs for structured clone / JSON.
  function serializeState(state) {
    return {
      ...state,
      libraries: Object.fromEntries(Object.entries(state.libraries).map(([k, lib]) =>
        [k, { ...lib, entries: [...lib.entries.entries()] }]))
    };
  }
  function deserializeState(raw) {
    const libraries = {};
    for (const [k, lib] of Object.entries(raw.libraries || {})) {
      libraries[k] = { ...lib, entries: new Map(lib.entries) };
    }
    return { ...raw, libraries };
  }
  init._deserializeState = deserializeState;

  // M2b import: replace the whole project with an exported one. Validates
  // the schema version (explicit refusal, never silent), clears all stores
  // (commands/snapshots/facts — a new project, not a merge), then writes
  // the commands sequentially under the lock.
  async function replaceProject(storage, data, opts) {
    if (!data || !Array.isArray(data.commands)) throw new Error('import: no commands array');
    if (data.schemaVersion > State.SCHEMA_VERSION) {
      throw new Error('import: project schema v' + data.schemaVersion + ' is newer than this build (v' + State.SCHEMA_VERSION + ')');
    }
    const withLock = opts && opts.withLock ? opts.withLock : (async (n, fn) => fn());
    return withLock('project', async () => {
      for (const store of ['commands', 'snapshots', 'facts']) {
        const all = await storage.all(store);
        for (const rec of all) await storage.del(store, rec.key);
      }
      await storage.del('meta', 'cmdSeq');
      let seq = 0;
      for (const cmd of data.commands) await storage.put('commands', ++seq, cmd);
      await storage.put('meta', 'cmdSeq', seq);
      await storage.put('meta', 'project', { schemaVersion: State.SCHEMA_VERSION, createdAt: Date.now(), importedAt: Date.now() });
      await storage.put('meta', 'hadData', true);
      return { applied: seq };
    });
  }

  return { init, serializeState, deserializeState, replaceProject, SNAPSHOT_EVERY, KEEP_SNAPSHOTS };
});
