// platform/registry.js — M4 P1 project registry (docs/project-management.md)
// State machine: creating -> ready; ready -> trashed -> deleting -> (gone);
// repointing arrives with the P2 pointer-switch flow. Reconciliation tags:
// orphan (DB without an entry) / missing (entry without a DB).
// Every mutation runs inside the 'registry' lock; a registry broadcast
// channel notifies other tabs (lock = exclusion, channel = notification,
// same discipline as the store).
(function (root, factory) {
  if (typeof module === 'object' && module.exports) module.exports = factory();
  else { (root.Platform = root.Platform || {}).Registry = factory(); }
})(typeof self !== 'undefined' ? self : this, function () {
  'use strict';

  function uuid() {
    if (typeof crypto !== 'undefined' && crypto.randomUUID) return crypto.randomUUID();
    return 'p' + Date.now().toString(36) + '-' + Math.random().toString(36).slice(2, 10);
  }
  const dbOf = id => id === 'default' ? 'webmm-platform' : 'webmm-project-' + id;
  // lock names are GLOBAL to the origin (not per-DB) — the default project
  // keeps the legacy 'project' name forever so old tabs stay mutually
  // exclusive with new code (docs §1.4)
  const lockOf = id => id === 'default' ? 'project' : 'project:' + id;

  function create(opts) {
    opts = opts || {};
    const DB = opts.DB;
    const mem = opts.mem;                       // tests: mem registry
    const storage = opts.storage || (mem ? DB.createRegistryStorage('mem') : DB.createRegistryStorage());
    let channel = null;
    // guard on browser context: Node 18+ also has BroadcastChannel and an
    // open channel keeps the Node event loop alive (tests would never exit)
    if (typeof window !== 'undefined') {
      try { channel = new BroadcastChannel('webmm-registry'); } catch (e) { /* worker without BC */ }
    }
    const notify = (kind, projectId, extra) => { try { if (channel) channel.postMessage({ kind, projectId, ...extra }); } catch (e) {} };
    const withRegistryLock = fn => DB.withLock('registry', fn);

    async function open() { if (storage.open) await storage.open(); return storage; }

    // zero-copy, idempotent: the legacy singleton becomes the 'default'
    // project by registry entry alone — no data moves (docs §1.3)
    async function ensureMigrated() {
      return withRegistryLock(async () => {
        if (!(await storage.get('default'))) {
          await storage.put('default', {
            dbName: 'webmm-platform', name: '默认项目', state: 'ready',
            createdAt: Date.now(), lastOpened: Date.now(), approxBytes: null,
          });
        }
        return { id: 'default', ...(await storage.get('default')) };
      });
    }

    async function list() {
      const all = await storage.all();
      return all.sort((a, b) => (b.lastOpened || 0) - (a.lastOpened || 0));
    }
    const get = id => storage.get(id);

    async function create(name) {
      const id = uuid();
      return withRegistryLock(async () => {
        const rec = {
          dbName: dbOf(id), name: name || ('项目 ' + id.slice(0, 6)), state: 'creating',
          createdAt: Date.now(), lastOpened: Date.now(), approxBytes: null,
        };
        await storage.put(id, rec);
        return { id, ...rec };
      });
    }
    async function markReady(id) {
      return withRegistryLock(async () => {
        const rec = await storage.get(id);
        if (!rec) throw new Error('registry: unknown project ' + id);
        rec.state = 'ready';
        await storage.put(id, rec);
        notify('ready', id);
        return rec;
      });
    }
    async function rename(id, name) {
      return withRegistryLock(async () => {
        const rec = await storage.get(id);
        if (!rec) throw new Error('registry: unknown project ' + id);
        rec.name = name;
        await storage.put(id, rec);
        notify('renamed', id);
        return rec;
      });
    }
    async function touch(id, bytes) {
      return withRegistryLock(async () => {
        const rec = await storage.get(id);
        if (!rec) return null;
        rec.lastOpened = Date.now();
        if (typeof bytes === 'number') rec.approxBytes = bytes;
        await storage.put(id, rec);
        return rec;
      });
    }
    async function trash(id) {
      return withRegistryLock(async () => {
        const rec = await storage.get(id);
        if (!rec) throw new Error('registry: unknown project ' + id);
        rec.state = 'trashed';
        rec.trashedAt = Date.now();
        await storage.put(id, rec);
        notify('trashed', id);
        return rec;
      });
    }
    async function purge(id) {
      // drop the database and the entry (soft delete's final step)
      return withRegistryLock(async () => {
        const rec = await storage.get(id);
        await storage.del(id);
        notify('purged', id);
        if (rec && typeof indexedDB !== 'undefined' && indexedDB.deleteDatabase) {
          await new Promise(res => {
            const r = indexedDB.deleteDatabase(rec.dbName);
            r.onsuccess = r.onerror = r.onblocked = () => res();
          });
        }
      });
    }
    async function removeMem(id) {   // mem-registry tests: no real DB to drop
      return withRegistryLock(async () => { await storage.del(id); notify('purged', id); });
    }

    // ---- pointer repoint (docs/project-management.md §3) ----
    // restore-into-current writes a NEW database then atomically swaps the
    // entry's dbName. Between prepare and commit the entry is 'repointing'
    // and every write path must refuse (the write gate reads this DIRECTLY
    // from IDB — never a cached copy).
    async function repointPrepare(id) {
      return withRegistryLock(async () => {
        const rec = await storage.get(id);
        if (!rec) throw new Error('registry: unknown project ' + id);
        if (rec.state !== 'ready') throw new Error('registry: repoint requires ready state, got ' + rec.state);
        const oldDb = rec.dbName, newDb = dbOf(uuid());
        rec.state = 'repointing';
        rec.repointingAt = Date.now();   // R2-2: stale-repointing detection
        await storage.put(id, rec);
        return { id, oldDb, newDb };
      });
    }
    async function repointCommit(id, newDb, oldDb) {
      return withRegistryLock(async () => {
        const rec = await storage.get(id);
        if (!rec || rec.state !== 'repointing') throw new Error('registry: repointCommit without prepare');
        // U1 (review r1): the pointer swap must not race in-flight writes.
        // Writers hold the project lock while putting; taking it here (we
        // already hold 'registry' — ordering registry -> project is safe
        // because writers never take 'registry') drains them, and any NEW
        // writer now sees state 'repointing' and refuses.
        // R2-2 (r2): tombstone FIRST, pointer second. A crash between the
        // two leaves "tombstone + repointing" = interpretable interrupted
        // restore; the reverse order would orphan the old DB (reconcile
        // false-positive).
        await storage.put('db-' + uuid(), {
          dbName: oldDb, name: '(旧数据库)', state: 'trashed', dbOnly: true,
          trashedAt: Date.now(), approxBytes: rec.approxBytes || null,
        });
        await DB.withLock(lockOf(id), async () => {
          rec.dbName = newDb;
          rec.state = 'ready';
          rec.repointedAt = Date.now();
          // U4: carry the size estimate onto the new database entry
          await storage.put(id, rec);
        });
        notify('repointed', id, { oldDb, newDb });
        return rec;
      });
    }
    async function repointAbort(id, newDb) {
      return withRegistryLock(async () => {
        const rec = await storage.get(id);
        if (rec && rec.state === 'repointing') {
          rec.state = 'ready';   // old dbName untouched
          await storage.put(id, rec);
        }
        notify('repoint-aborted', id, {});
      });
    }
    // d (review r1): registry writes go through the registry lock —
    // export lineage (hash + time) is no exception
    // R2-3 (r2): structured export lineage — dedup requires BOTH hashes and
    // the same algorithm version. Legacy scalar values upcast on read.
    async function setExportHash(id, commandsHash, factsHash) {
      return withRegistryLock(async () => {
        const rec = await storage.get(id);
        if (!rec) return null;
        rec.lastExport = { commandsHash, factsHash: factsHash || null, algoVersion: 'v1-facts' };
        rec.lastExportAt = Date.now();
        await storage.put(id, rec);
        return rec;
      });
    }
    function exportLineage(rec) {
      // legacy scalar lastExportHash upcast: commands-only, pre-factsHash
      if (rec && rec.lastExport) return rec.lastExport;
      if (rec && rec.lastExportHash) return { commandsHash: rec.lastExportHash, factsHash: null, algoVersion: 'v1' };
      return null;
    }

    // P3: trashed PROJECT entries purge after the retention window
    // (drops the project DB and the registry entry)
    async function purgeTrashedProjects(retentionMs) {
      const now = Date.now();
      const dropped = [];
      for (const rec of await storage.all()) {
        if (rec.state === 'trashed' && rec.trashedAt && now - rec.trashedAt >= (retentionMs || 0)) dropped.push(rec);
      }
      const remaining = [];
      for (const rec of dropped) {
        if (rec.id === 'default') { remaining.push(rec); continue; }   // never purged
        try {
          if (rec.dbName !== 'webmm-platform') await dropDb(rec.dbName);
          await storage.del(rec.id);
        } catch (e) { remaining.push(rec); }   // blocked: retry next boot
      }
      return dropped.length - remaining.length;
    }

    // U5 (review r1): deleteDatabase must NOT resolve on 'blocked' — a held
    // connection means the DB still exists; keep the entry and retry on a
    // later boot instead of pretending it was dropped.
    const dropDb = dbName => new Promise((resolve, reject) => {
      if (typeof indexedDB === 'undefined' || !indexedDB.deleteDatabase) return resolve();
      const r = indexedDB.deleteDatabase(dbName);
      r.onsuccess = () => resolve();
      r.onerror = () => reject(r.error);
      r.onblocked = () => reject(new Error('deleteDatabase blocked (a tab still holds ' + dbName + ')'));
    });

    // startup reconciliation (docs §1.2): advisory only — the registry is
    // the source of truth, databases() can be stale on Safari
    async function reconcile() {
      const out = { orphans: [], missing: [], interruptedRestores: [] };
      // interrupted restores need no databases() — scan first (Node-safe)
      const now = Date.now();
      for (const r of await storage.all()) {
        if (r.state === 'repointing' && r.repointingAt && now - r.repointingAt > 5 * 60 * 1000) {
          out.interruptedRestores.push(r.id);
        }
      }
      if (typeof indexedDB === 'undefined' || !indexedDB.databases) return out;
      let dbs = [];
      try { dbs = await indexedDB.databases(); } catch (e) { return out; }
      const names = new Set(dbs.map(d => d.name));
      const registered = new Set((await storage.all()).map(r => r.dbName));
      for (const n of names) {
        if (/^webmm-project-/.test(n) && !registered.has(n)) out.orphans.push(n);
      }
      for (const r of await storage.all()) {
        if (r.state !== 'trashed' && r.state !== 'deleting' && !names.has(r.dbName)) out.missing.push(r.id);
      }
      return out;
    }

    const onRegistryMessage = fn => { if (channel) channel.onmessage = e => fn(e.data); };
    return {
      open, ensureMigrated, list, get, create, markReady, rename, touch, trash, purge,
      removeMem, reconcile, dbOf, lockOf, onRegistryMessage,
      repointPrepare, repointCommit, repointAbort, purgeTrashedProjects, setExportHash, exportLineage,
      _storage: storage,
    };
  }

  return { create, dbOf, lockOf };
});
