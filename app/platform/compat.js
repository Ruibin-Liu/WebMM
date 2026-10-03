// platform/compat.js — M1a bridge: the old UI's library path writes THROUGH
// the new store from day one (strangler rule, M0 spec §5). localStorage
// becomes a read-only one-time migration source.
(function (root, factory) {
  if (typeof module === 'object' && module.exports) module.exports = factory();
  else { (root.Platform = root.Platform || {}).Compat = factory(); }
})(typeof self !== 'undefined' ? self : this, function () {
  'use strict';
  const LEGACY_KEY = 'wb-searchLibrary';
  const LIB_ID = 'main';     // the single working library of the legacy UI

  let project = null;
  let identity = null;
  let initPromise = null;
  let activeStorage = null;
  let activeLock = null;
  let projectModule = null;
  let needsReinit = false;   // another tab wrote; callers may poll stale()

  // M2b import: replace the whole project with an exported one; the caller
  // then reloads the page (explicit, not silent — the in-memory projection
  // cannot be hot-swapped safely). Uses the storage opened at init.
  async function importProject(data) {
    if (!activeStorage || !projectModule) throw new Error('compat: not initialized');
    return projectModule.replaceProject(activeStorage, data, { withLock: activeLock });
  }

  // deps: { storage, Identity } — browser passes the IDB adapter and the
  // Identity module (with the InChIKey chain already injected).
  function init(deps) {
    if (initPromise) return initPromise;
    initPromise = (async () => {
      identity = deps.Identity;
      const DB = deps.DB || null;
      const storage = deps.storage || (DB ? DB.createIdbStorage() : null);
      if (!storage) throw new Error('compat: no storage');
      if (storage.open) await storage.open();
      const withLock = DB ? DB.withLock : null;
      const channel = DB ? DB.createStoreChannel() : null;
      activeStorage = storage; activeLock = withLock; projectModule = deps.Project;
      // M2b: read-only tabs re-init their projection when another tab writes
      if (channel) channel.onmessage = () => { needsReinit = true; };
      project = await deps.Project.init(storage, { withLock, channel });
      // durability script (persist/estimate/vanish-detect) — surfaced to
      // console + meta in M1a; the UI banner is M1b
      if (DB) { try { await DB.durability(storage); } catch (e) {} }
      return project;
    })();
    return initPromise;
  }

  // Save the parsed library inputs (canonical SMILES + names only — the
  // fingerprints are L2 deterministic caches and never persist).
  async function saveLibraryInputs(entries) {
    if (!project) throw new Error('compat: not initialized');
    const withIds = entries.map((e, i) => {
      const raw = JSON.stringify({ smiles: e.smiles, name: e.name || null });
      const mid = identity.molId(LIB_ID, i, raw);
      let key = null;
      try { key = identity.structureKey(e.smiles); } catch (err) { /* key fn unavailable — entries still load */ }
      return { molId: mid, ordinal: i, raw: { smiles: e.smiles, name: e.name || null }, structureKey: key };
    });
    await project.apply({ type: 'ImportLibrary', libraryId: LIB_ID, name: 'Working library', entries: withIds });
    await project.checkpoint();      // import complete = checkpoint event (spec §4)
    return withIds;                  // enriched: caller attaches molId to in-memory entries
  }

  // ---- M1b triage facade (pin/exclude/notes/undo over the command log) ----
  // Molecule-level overlay + provenance; rounds/propagation arrive in M1c.
  async function pin(molId, note, provenanceRound) {
    await project.apply({ type: 'Pin', molId, note: note || null, provenanceRound: provenanceRound || null });
  }
  async function unpin(molId) { await project.apply({ type: 'Unpin', molId }); }
  async function exclude(molId, reason) {
    await project.apply({ type: 'Exclude', molId, reason: reason || null });
  }
  async function include(molId) { await project.apply({ type: 'Include', molId }); }
  async function undo() { return project.undo(); }
  function getPins() { return project ? project.state.pins : {}; }
  function getExcludes() { return project ? project.state.excludes : {}; }
  function ready() { return !!project; }

  // ---- M2 round scores as L3 facts (per M0 spec §3: stochastic results
  // stored WITH their provenance; deterministic recompute does not apply) ----
  const __scoreCache = new Map();   // roundId -> hits [{molId, score, ...}]
  async function saveRoundScores(roundId, hits) {
    if (!project) return;
    __scoreCache.set(roundId, hits);
    try {
      await project.storagePut('facts', 'scores:' + roundId, {
        engine: 'webmm-' + (typeof WEBMM_VERSION !== 'undefined' ? WEBMM_VERSION : ''),
        hits,
      });
    } catch (e) { console.warn('[platform] round scores not persisted:', e); }
  }
  async function getRoundScores(roundId) {
    if (__scoreCache.has(roundId)) return __scoreCache.get(roundId);
    if (!project) return null;
    const rec = await project.storageGet('facts', 'scores:' + roundId);
    if (rec && Array.isArray(rec.hits)) { __scoreCache.set(roundId, rec.hits); return rec.hits; }
    return null;
  }
  function cachedRoundScores() { return __scoreCache; }
  function projectApi() { return project; }
  function stale() { return needsReinit; }

  async function exportProjectData() {
    if (!project) return null;
    return project.exportProject();
  }

  // Restore on startup: IndexedDB projection first; if never-had, one-time
  // read-only migration from the legacy localStorage key.
  async function restoreInputs(migrateFromLegacy) {
    if (!project) return null;
    const ids = Object.keys(project.state.libraries);
    if (ids.length) {
      const inputs = [];
      for (const id of ids) inputs.push(...libInputs(id));
      return inputs;
    }
    if (typeof migrateFromLegacy === 'function') {
      const legacy = migrateFromLegacy();     // () => [{smiles,name}] | null (read-only)
      if (legacy && legacy.length) {
        await saveLibraryInputs(legacy);
        return legacy;
      }
    }
    return null;
  }
  function libInputs(id) {
    const lib = project.state.libraries[id];
    if (!lib) return [];
    const out = [];
    for (const [molId, e] of lib.entries) out.push({ molId, ordinal: e.ordinal, smiles: e.raw.smiles, name: e.raw.name });
    out.sort((a, b) => a.ordinal - b.ordinal);
    return out;
  }

  async function clearLibrary() {
    if (!project) return;
    if (project.state.libraries[LIB_ID]) {
      await project.apply({ type: 'RemoveLibrary', libraryId: LIB_ID });
      await project.checkpoint();
    }
  }

  return { init, saveLibraryInputs, restoreInputs, clearLibrary, LEGACY_KEY, LIB_ID,
    pin, unpin, exclude, include, undo, getPins, getExcludes, ready, projectApi,
    saveRoundScores, getRoundScores, cachedRoundScores, exportProjectData,
    importProject, stale };
});
