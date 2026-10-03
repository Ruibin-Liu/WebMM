// platform/state.js — M1a command set + pure reducer (M0 spec §4)
// The command log is the source of truth; this module projects it. Pure:
// no I/O, no DOM, Node-testable. Inverses implement the undo rule
// (undo = append an inverse command, never truncate).
(function (root, factory) {
  if (typeof module === 'object' && module.exports) module.exports = factory();
  else { (root.Platform = root.Platform || {}).State = factory(); }
})(typeof self !== 'undefined' ? self : this, function () {
  'use strict';

  const SCHEMA_VERSION = 1;

  // ---- command types (frozen at M1a; spec addendum: + RemoveLibrary) ----
  // ImportLibrary {libraryId, name, entries:[{ordinal, raw, structureKey?}]}
  // RemoveLibrary {libraryId}
  // IdentityChanged {molId, fromKey, toKey}
  // IdentityMerged {intoMolId, fromMolId}
  // IdentitySplit {fromMolId, toMolIds:[...]}
  // CreateRound {roundId, parentId|null, inputRef, querySpec}
  // SetThreshold {roundId, value}
  // Pin {molId, note?, provenanceRound?}
  // Unpin {molId}
  // Exclude {molId, reason?}
  // Include {molId}          (per-round override arrives with rounds in M1c)
  // ArchiveSubtree {roundId}
  // RenameRound {roundId, name}

  function initialState() {
    return {
      schemaVersion: SCHEMA_VERSION,
      libraries: {},           // libraryId -> {name, importedAt, entries: Map(molId -> {ordinal, raw, structureKey, rowId})}
      rounds: {},              // roundId -> {parentId, inputRef, querySpec, threshold, archived, name}
      pins: {},                // molId -> {note, provenanceRound}
      excludes: {},            // molId -> {reason}
      maxCommandId: 0,
    };
  }

  function err(msg) { return { __error: msg }; }

  // Returns {state, error?} — reducer never throws; errors are returned so
  // the log replay can surface the first illegal transition deterministically.
  function applyCommand(prev, cmd) {
    const s = cloneShallow(prev);
    s.maxCommandId = prev.maxCommandId + 1;
    const id = cmd && cmd.libraryId;
    switch (cmd && cmd.type) {
      case 'ImportLibrary': {
        if (!cmd.libraryId || !Array.isArray(cmd.entries)) return { state: prev, error: 'ImportLibrary: libraryId and entries required' };
        const entries = new Map();
        for (const e of cmd.entries) {
          if (e.molId == null || e.raw == null) return { state: prev, error: 'ImportLibrary: each entry needs molId and raw' };
          entries.set(e.molId, { ordinal: e.ordinal, raw: e.raw, structureKey: e.structureKey || null, rowId: e.rowId });
        }
        s.libraries[id] = { name: cmd.name || id, importedAt: cmd.importedAt || Date.now(), entries };
        return { state: s };
      }
      case 'RemoveLibrary': {
        if (!s.libraries[id]) return { state: prev, error: 'RemoveLibrary: unknown library' };
        delete s.libraries[id];
        return { state: s };
      }
      case 'IdentityChanged': {
        const hit = findEntry(s, cmd.molId);
        if (!hit) return { state: prev, error: 'IdentityChanged: unknown molId' };
        hit.entry.structureKey = cmd.toKey;
        return { state: s };
      }
      case 'IdentityMerged': {
        const from = findEntry(s, cmd.fromMolId);
        const into = findEntry(s, cmd.intoMolId);
        if (!from || !into) return { state: prev, error: 'IdentityMerged: unknown molId' };
        // from-row is tombstoned at the allocator level by the caller; the
        // overlay keeps its provenance history (spec: keep both, surface new)
        if (s.pins[cmd.fromMolId] && !s.pins[cmd.intoMolId]) s.pins[cmd.intoMolId] = s.pins[cmd.fromMolId];
        delete s.pins[cmd.fromMolId];
        delete s.excludes[cmd.fromMolId];
        from.lib.entries.delete(cmd.fromMolId);
        return { state: s };
      }
      case 'IdentitySplit': {
        const src = findEntry(s, cmd.fromMolId);
        if (!src) return { state: prev, error: 'IdentitySplit: unknown molId' };
        if (!Array.isArray(cmd.children) || !cmd.children.length) return { state: prev, error: 'IdentitySplit: children required' };
        // neither child inherits the old row (spec §2): new molIds arrive as
        // entries via their own ImportLibrary/synthetic command; the source
        // row is removed here (allocator tombstones by the caller).
        src.lib.entries.delete(cmd.fromMolId);
        delete s.pins[cmd.fromMolId];
        delete s.excludes[cmd.fromMolId];
        return { state: s };
      }
      case 'CreateRound': {
        if (!cmd.roundId || s.rounds[cmd.roundId]) return { state: prev, error: 'CreateRound: unique roundId required' };
        if (cmd.parentId && !s.rounds[cmd.parentId]) return { state: prev, error: 'CreateRound: unknown parentId' };
        s.rounds[cmd.roundId] = {
          parentId: cmd.parentId || null, inputRef: cmd.inputRef || null,
          querySpec: cmd.querySpec || null, threshold: null, archived: false, name: null,
        };
        return { state: s };
      }
      case 'SetThreshold': {
        const r = s.rounds[cmd.roundId];
        if (!r) return { state: prev, error: 'SetThreshold: unknown round' };
        r.threshold = cmd.value;
        return { state: s };
      }
      case 'Pin': {
        if (!findEntry(s, cmd.molId) && !s.pins[cmd.molId]) return { state: prev, error: 'Pin: unknown molId' };
        s.pins[cmd.molId] = { note: cmd.note || null, provenanceRound: cmd.provenanceRound || null };
        return { state: s };
      }
      case 'Unpin': {
        delete s.pins[cmd.molId];
        return { state: s };
      }
      case 'Exclude': {
        if (!findEntry(s, cmd.molId) && !s.excludes[cmd.molId]) return { state: prev, error: 'Exclude: unknown molId' };
        s.excludes[cmd.molId] = { reason: cmd.reason || null };
        return { state: s };
      }
      case 'Include': {
        delete s.excludes[cmd.molId];
        return { state: s };
      }
      case 'ArchiveSubtree': {
        // archive the node and all descendants (reversible — M0 spec §6)
        const mark = (rid) => {
          const r = s.rounds[rid];
          if (!r) return;
          r.archived = true;
          for (const [k, v] of Object.entries(s.rounds)) if (v.parentId === rid) mark(k);
        };
        if (!s.rounds[cmd.roundId]) return { state: prev, error: 'ArchiveSubtree: unknown round' };
        mark(cmd.roundId);
        return { state: s };
      }
      case 'RenameRound': {
        const r = s.rounds[cmd.roundId];
        if (!r) return { state: prev, error: 'RenameRound: unknown round' };
        r.name = cmd.name;
        return { state: s };
      }
      default:
        return { state: prev, error: 'unknown command type: ' + (cmd && cmd.type) };
    }
  }

  function cloneShallow(prev) {
    return {
      schemaVersion: prev.schemaVersion,
      libraries: cloneLibs(prev.libraries),
      rounds: cloneRounds(prev.rounds),
      pins: { ...prev.pins },
      excludes: { ...prev.excludes },
      maxCommandId: prev.maxCommandId,
    };
  }
  function cloneLibs(libs) {
    const out = {};
    for (const [k, lib] of Object.entries(libs)) {
      out[k] = { name: lib.name, importedAt: lib.importedAt, entries: new Map(lib.entries) };
    }
    return out;
  }
  function cloneRounds(rounds) {
    const out = {};
    for (const [k, r] of Object.entries(rounds)) out[k] = { ...r };
    return out;
  }
  function findEntry(s, molId) {
    for (const lib of Object.values(s.libraries)) {
      const entry = lib.entries.get(molId);
      if (entry) return { lib, entry };
    }
    return null;
  }

  // ---- inverses (undo = append the inverse; M0 spec §4) ----
  // Structural commands with dependents must REFUSE (the caller checks
  // `hasDependents` before offering undo; the inverse builder here encodes
  // the policy for the commands that are always undoable).
  function inverse(cmd, state) {
    switch (cmd.type) {
      case 'Pin': return { type: 'Unpin', molId: cmd.molId };
      case 'Unpin': return null;                       // needs prior value; caller snapshots overlay for undo stack
      case 'Exclude': return { type: 'Include', molId: cmd.molId };
      case 'Include': return null;
      case 'SetThreshold': return { type: 'SetThreshold', roundId: cmd.roundId, value: state.rounds[cmd.roundId] ? state.rounds[cmd.roundId].threshold : null };
      case 'RenameRound': return { type: 'RenameRound', roundId: cmd.roundId, name: state.rounds[cmd.roundId] ? state.rounds[cmd.roundId].name : null };
      case 'ImportLibrary': return null;               // dependents possible (rounds reference it) — policy: refuse
      case 'RemoveLibrary': return null;               // destructive; refuse (archive-style alternative arrives with rounds)
      case 'CreateRound': return null;                 // children possible — refuse per spec
      case 'ArchiveSubtree': return null;              // reversible via Unarchive (M1c); not in frozen set yet
      default: return null;
    }
  }

  // Dependents check for refusal messages (undo-with-children rule).
  function hasDependents(cmd, state) {
    if (cmd.type === 'CreateRound') {
      return Object.values(state.rounds).some(r => r.parentId === cmd.roundId);
    }
    if (cmd.type === 'ImportLibrary') {
      return Object.values(state.rounds).some(r => r.inputRef && String(r.inputRef).startsWith('lib:' + cmd.libraryId));
    }
    return false;
  }

  // ---- projections ----
  // The old UI consumes libraries as [{smiles, name}] (input-only
  // persistence: fingerprints are L2 deterministic caches, recomputed).
  function projectLibraryInputs(state, libraryId) {
    const lib = state.libraries[libraryId];
    if (!lib) return null;
    const out = [];
    for (const [molId, e] of lib.entries) out.push({ molId, ...e.raw });
    out.sort((a, b) => (a.ordinal ?? 0) - (b.ordinal ?? 0));
    return out;
  }
  function activeLibraryIds(state) { return Object.keys(state.libraries); }

  return {
    SCHEMA_VERSION, initialState, applyCommand, inverse, hasDependents,
    projectLibraryInputs, activeLibraryIds, findEntry,
  };
});
