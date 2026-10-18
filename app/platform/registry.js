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
        await storage.put(id, rec);
        return { id, oldDb, newDb };
      });
    }
    async function repointCommit(id, newDb, oldDb) {
      return withRegistryLock(async () => {
        const rec = await storage.get(id);
        if (!rec || rec.state !== 'repointing') throw new Error('registry: repointCommit without prepare');
        rec.dbName = newDb;
        rec.state = 'ready';
        rec.repointedAt = Date.now();
        await storage.put(id, rec);
        // the old database goes to the trash record (retention, then purge)
        const trash = (await storage.get('__trash')) || { dbs: [] };
        trash.dbs.push({ dbName: oldDb, at: Date.now() });
        await storage.put('__trash', trash);
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
    // P3: trashed PROJECT entries purge after the retention window
    // (drops the project DB and the registry entry)
    async function purgeTrashedProjects(retentionMs) {
      const now = Date.now();
      const dropped = [];
      for (const rec of await storage.all()) {
        if (rec.state === 'trashed' && rec.trashedAt && now - rec.trashedAt >= (retentionMs || 0)) dropped.push(rec);
      }
      for (const rec of dropped) {
        if (rec.id === 'default') continue;   // the default project is never purged
        await storage.del(rec.id);
        if (typeof indexedDB !== 'undefined' && indexedDB.deleteDatabase && rec.dbName !== 'webmm-platform') {
          await new Promise(res => { const r = indexedDB.deleteDatabase(rec.dbName); r.onsuccess = r.onerror = r.onblocked = () => res(); });
        }
      }
      if (dropped.length) notify('purged-batch', null, {});
      return dropped.length;
    }

    async function purgeTrashedDbs(retentionMs) {
      // drop trashed databases older than the retention window
      const trash = (await storage.get('__trash')) || { dbs: [] };
      const now = Date.now();
      const keep = [], drop = [];
      for (const t of trash.dbs) (now - t.at >= (retentionMs || 0) ? drop : keep).push(t);
      if (drop.length) {
        await storage.put('__trash', { dbs: keep });
        if (typeof indexedDB !== 'undefined' && indexedDB.deleteDatabase) {
          for (const t of drop) await new Promise(res => {
            const r = indexedDB.deleteDatabase(t.dbName);
            r.onsuccess = r.onerror = r.onblocked = () => res();
          });
        }
      }
      return drop.length;
    }

    // startup reconciliation (docs §1.2): advisory only — the registry is
    // the source of truth, databases() can be stale on Safari
    async function reconcile() {
      const out = { orphans: [], missing: [] };
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
      repointPrepare, repointCommit, repointAbort, purgeTrashedDbs, purgeTrashedProjects,
      _storage: storage,
    };
  }

  return { create, dbOf, lockOf };
});
