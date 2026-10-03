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

  async function init(storage) {
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
    const p = {
      storage, meta, state,
      _commandsSinceSnapshot: tail.length,
      _undoStack: [],           // {inverse, describes} of undoable commands
      _cmdId: cmdRecs.length ? cmdRecs[cmdRecs.length - 1].key : 0,
      _hadData: cmdRecs.length > 0 || snaps.length > 0,
    };
    try { await storage.put('meta', 'project', meta); } catch (e) {}
    return api(p);
  }

  function api(p) {
    return {
      get state() { return p.state; },
      get schemaVersion() { return p.state.schemaVersion; },

      // apply a command: validate -> reduce -> persist -> snapshot trigger
      async apply(cmd) {
        // prior-value capture for inverses that need it (Unpin needs the old
        // pin, Include the old exclude) — BEFORE the reduce (spec §4)
        let pre = null;
        if (cmd.type === 'Unpin' && p.state.pins[cmd.molId]) {
          pre = { type: 'Pin', molId: cmd.molId, ...p.state.pins[cmd.molId] };
        } else if (cmd.type === 'Include' && p.state.excludes[cmd.molId]) {
          pre = { type: 'Exclude', molId: cmd.molId, ...p.state.excludes[cmd.molId] };
        }
        const res = State.applyCommand(p.state, cmd);
        if (res.error) throw new Error('command rejected: ' + res.error);
        p._cmdId += 1;
        p.state = res.state;
        await p.storage.put('commands', p._cmdId, cmd);
        if (!p._hadData) { p._hadData = true; try { await p.storage.put('meta', 'hadData', true); } catch (e) {} }
        const inv = pre || State.inverse(cmd, res.state);
        if (inv) p._undoStack.push({ inverse: inv, describes: cmd });
        p._commandsSinceSnapshot += 1;
        if (p._commandsSinceSnapshot >= SNAPSHOT_EVERY) await snapshot(p);
        return { commandId: p._cmdId };
      },

      // Undo = append an inverse command (spec §4). Structural commands
      // with dependents are refused with their dependent list.
      async undo() {
        while (p._undoStack.length) {
          const top = p._undoStack.pop();
          if (State.hasDependents(top.describes, p.state)) {
            return { ok: false, reason: 'has-dependents', describes: top.describes };
          }
          await this.apply(top.inverse);
          return { ok: true, undone: top.describes };
        }
        return { ok: false, reason: 'empty' };
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
    };
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

  return { init, serializeState, deserializeState, SNAPSHOT_EVERY, KEEP_SNAPSHOTS };
});
