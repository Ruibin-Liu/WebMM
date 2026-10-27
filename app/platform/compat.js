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
  let activeProjectId = 'default';
  let writableGuardRef = null;

  // M2b import: replace the whole project with an exported one; the caller
  // then reloads the page (explicit, not silent — the in-memory projection
  // cannot be hot-swapped safely). Uses the storage opened at init.
  async function importProject(data) {
    if (!activeStorage || !projectModule) throw new Error('compat: not initialized');
    const r = await projectModule.replaceProject(activeStorage, data, { withLock: activeLock });
    if (project && project.invalidate) project.invalidate();   // old projection must never write again
    return r;
  }

  // deps: { storage, Identity } — browser passes the IDB adapter and the
  // Identity module (with the InChIKey chain already injected).
  function init(deps) {
    if (initPromise) return initPromise;
    initPromise = (async () => {
      identity = deps.Identity;
      const DB = deps.DB || null;
      // M4 P1: per-project storage/lock/channel (docs/project-management.md
      // §1.4). deps.storage still wins (tests inject mem storage); absent
      // it, dbName selects the project database. The lock adapter maps the
      // store's internal 'project' lock onto the project-scoped name — the
      // default keeps the legacy 'project' name so old tabs stay exclusive.
      const projectId = deps.projectId || 'default';
      const dbName = deps.dbName || (projectId === 'default' ? undefined : 'webmm-project-' + projectId);
      const storage = deps.storage || (DB ? DB.createIdbStorage(dbName) : null);
      if (!storage) throw new Error('compat: no storage');
      if (storage.open) await storage.open();
      // M4 P2 §3 write gate: reads the registry entry straight from IDB on
      // every write — a cached copy would make the gate worthless
      const registryStorage = deps.registryStorage || null;
      const writableGuard = registryStorage
        ? async () => {
            const rec = await registryStorage.get(projectId);
            if (rec && (rec.state === 'repointing' || rec.state === 'deleting')) {
              throw new Error('project is being repointed — writes refused, reload after it completes');
            }
          }
        : null;
      const lockNameFor = pid => pid === 'default' ? 'project' : 'project:' + pid;
      const withLock = DB
        ? (name, fn) => DB.withLock(name === 'project' ? lockNameFor(projectId) : name, fn)
        : null;
      writableGuardRef = writableGuard;
      const channel = DB ? DB.createStoreChannel(dbName) : null;
      activeProjectId = projectId;
      activeStorage = storage; activeLock = withLock; projectModule = deps.Project;
      // Top-5 (docs/task-contract.md §6): the task side table rides the same
      // per-project storage as commands/facts — Done stays one transaction.
      activeStorage = storage;
      // M2b: read-only tabs re-init their projection when another tab writes
      if (channel) channel.onmessage = () => { needsReinit = true; };
      project = await deps.Project.init(storage, { withLock, channel, writableGuard });
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
  async function pin(molId, note, provenanceRound, lineage) {
    // lineage (M4 P2): {productSmiles, analogRun} makes the pin
    // self-contained — it resolves without the run's facts (docs §4)
    await project.apply({ type: 'Pin', molId, note: note || null, provenanceRound: provenanceRound || null,
      productSmiles: lineage && lineage.productSmiles || null,
      provenanceAnalogRun: lineage && lineage.analogRun || null });
  }

  async function applyWithFacts(cmd, factPuts, extraOps) {
    if (!project || !project.applyWithFacts) throw new Error('compat: applyWithFacts unavailable');
    return project.applyWithFacts(cmd, factPuts, extraOps);   // Top-5: extraOps (task side table) ride the same transaction
  }
  async function getFact(key) {
    if (!project || !project.storageGet) return null;
    return project.storageGet('facts', key);
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
    // R2-1 (r2): L1 facts are never silently dropped — a gate refusal
    // THROWS so the caller can surface "retry after the restore finishes"
    if (writableGuardRef) await writableGuardRef();
    __scoreCache.set(roundId, hits);
    try {
      await project.storagePut('facts', 'scores:' + roundId, {
        engine: 'webmm-' + (typeof WEBMM_VERSION !== 'undefined' ? WEBMM_VERSION : ''),
        layer: 'stochastic',   // L3 result with lineage (docs/project-management.md §2.1)
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

  // M4 P1: all facts records for export; untagged legacy records are
  // upcast to 'stochastic' at read time (docs/project-management.md §2.1)
  async function listFacts(layers) {
    // a (review r1): export filters layer ∈ {fact, stochastic}; L2/derived
    // facts stay out of the envelope. Untagged legacy records upcast.
    const allow = layers || ['fact', 'stochastic'];
    if (!project || !project.storageAll) return [];
    const out = [];
    for (const r of await project.storageAll('facts')) {
      const v = r.value || {};
      const layer = v.layer || 'stochastic';
      if (!allow.includes(layer)) continue;
      out.push({ key: r.key, value: { ...v, layer } });
    }
    return out;
  }

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
    }
    await project.checkpoint();   // ensure a removed-library projection wins replay
  }

  return { init, saveLibraryInputs, restoreInputs, clearLibrary, LEGACY_KEY, LIB_ID,
    pin, unpin, exclude, include, undo, getPins, getExcludes, ready, projectApi, applyWithFacts, getFact, getStorage: () => activeStorage,
    saveRoundScores, getRoundScores, cachedRoundScores, exportProjectData, listFacts,
    importProject, stale };
});
