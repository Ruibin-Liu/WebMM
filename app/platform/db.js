// platform/db.js — M1a storage adapters (M0 spec §3/§5)
// Two implementations of one narrow interface: real IndexedDB (browser)
// and an in-memory adapter (tests/Node). No IDB shim games.
(function (root, factory) {
  if (typeof module === 'object' && module.exports) module.exports = factory();
  else { (root.Platform = root.Platform || {}).DB = factory(); }
})(typeof self !== 'undefined' ? self : this, function () {
  'use strict';

  const DB_NAME = 'webmm-platform';
  const DB_VERSION = 2;
  const STORES = ['meta', 'commands', 'snapshots', 'facts'];  // facts: L3 results (round scores), keyed 'scores:<roundId>'

  // MIGRATIONS: Map<fromVersion, (txHelpers) => Promise> — the registered
  // hook pattern (M0 spec §5). v1 has none; the registry exists from the
  // first write so v1→v2 is a code path, not a crisis.
  const MIGRATIONS = new Map();

  // ---- in-memory adapter (tests, Node) ----
  function createMemStorage() {
    const data = Object.fromEntries(STORES.map(s => [s, new Map()]));  // mirror IDB schema exactly
    return {
      kind: 'mem',
      async get(store, key) { return data[store].get(key); },
      async put(store, key, val) { data[store].set(key, val); },
      async del(store, key) { data[store].delete(key); },
      async all(store) { return [...data[store].entries()].map(([k, v]) => ({ key: k, value: v })); },
      async count(store) { return data[store].size; },
      async putBatch(ops) { for (const o of ops) data[o.store].set(o.key, o.value); },   // trivially atomic in mem
      _data: data,
    };
  }

  // ---- IndexedDB adapter (browser) ----
  // dbName: each project gets its own database (docs/project-management.md
  // §1); the default (legacy singleton) keeps the historical name.
  function createIdbStorage(dbName) {
    dbName = dbName || DB_NAME;
    let db = null;
    function open() {
      return new Promise((resolve, reject) => {
        const req = indexedDB.open(dbName, DB_VERSION);
        req.onupgradeneeded = (ev) => {
          const d = req.result;
          // run registered migrations sequentially (none at v1)
          for (const store of STORES) if (!d.objectStoreNames.contains(store)) d.createObjectStore(store);
        };
        req.onsuccess = () => {
          db = req.result;
          // let upgrades/deletes from other tabs proceed (§1.2/§3): a held
          // connection would block them forever
          db.onversionchange = () => { try { db.close(); } catch (e) {} };
          resolve(db);
        };
        req.onerror = () => reject(req.error);
      });
    }
    function tx(store, mode, fn) {
      return new Promise((resolve, reject) => {
        const t = db.transaction(store, mode);
        const os = t.objectStore(store);
        const out = fn(os);
        t.oncomplete = () => resolve(out && out.result !== undefined ? out.result : out);
        t.onerror = () => reject(t.error);
      });
    }
    return {
      kind: 'idb',
      open,
      async get(store, key) {
        return new Promise((resolve, reject) => {
          // IDB requests resolve per-request; wrap explicitly
          const t = db.transaction(store, 'readonly');
          const req = t.objectStore(store).get(key);
          req.onsuccess = () => resolve(req.result);
          req.onerror = () => reject(req.error);
        });
      },
      async put(store, key, val) {
        return new Promise((resolve, reject) => {
          const t = db.transaction(store, 'readwrite');
          t.objectStore(store).put(val, key);
          // resolve on TRANSACTION commit (r2: req.onsuccess fires before the
          // tx is durable — the repoint drain relies on commit completion)
          t.oncomplete = () => resolve();
          t.onerror = () => reject(t.error);
          t.onabort = () => reject(t.error || new Error('put aborted'));
        });
      },
      async del(store, key) {
        return new Promise((resolve, reject) => {
          const t = db.transaction(store, 'readwrite');
          const req = t.objectStore(store).delete(key);
          req.onsuccess = () => resolve();
          req.onerror = () => reject(req.error);
        });
      },
      async all(store) {
        return new Promise((resolve, reject) => {
          if (!db) return resolve([]);   // reload can detach before cursor settles
          const t = db.transaction(store, 'readonly');
          const out = [];
          const req = t.objectStore(store).openCursor();
          req.onsuccess = () => {
            const cur = req.result;
            if (cur) { out.push({ key: cur.key, value: cur.value }); cur.continue(); }
            else resolve(out);
          };
          req.onerror = () => reject(req.error);
        });
      },
      async count(store) {
        return new Promise((resolve, reject) => {
          const t = db.transaction(store, 'readonly');
          const req = t.objectStore(store).count();
          req.onsuccess = () => resolve(req.result);
          req.onerror = () => reject(req.error);
        });
      },
      // ONE transaction across stores — the analog-run commit path (docs §4:
      // command + facts land atomically or not at all)
      async putBatch(ops) {
        return new Promise((resolve, reject) => {
          const stores = [...new Set(ops.map(o => o.store))];
          const t = db.transaction(stores, 'readwrite');
          for (const o of ops) t.objectStore(o.store).put(o.value, o.key);
          t.oncomplete = () => resolve();
          t.onerror = () => reject(t.error);
          t.onabort = () => reject(t.error || new Error('putBatch aborted'));
        });
      },
    };
  }

  // ---- durability script (M0 spec §3): persist/estimate/vanish-detect ----
  // M1a surfaces results via console + meta records; the UI banner is M1b.
  async function durability(storage) {
    const report = { persisted: null, estimate: null, vanished: 'never-had' };
    try {
      if (navigator.storage && navigator.storage.persist) {
        report.persisted = await navigator.storage.persist();
      }
      if (navigator.storage && navigator.storage.estimate) {
        report.estimate = await navigator.storage.estimate();
      }
    } catch (e) { /* non-fatal */ }
    // vanished detection: meta.hadData set on first successful write;
    // hadData && stores empty => evicted (not "never had")
    const had = await storage.get('meta', 'hadData');
    const nCmd = await storage.count('commands');
    const nSnap = await storage.count('snapshots');
    report.vanished = had ? (nCmd + nSnap === 0 ? 'evicted' : 'ok') : (nCmd + nSnap > 0 ? 'ok' : 'never-had');
    try { await storage.put('meta', 'durability', { ...report, at: Date.now() }); } catch (e) {}
    if (report.vanished === 'evicted') console.warn('[platform] IndexedDB was evicted (browser storage pressure) — the last session\'s facts are gone; consider installing the app or exporting projects.');
    if (report.persisted === false) console.info('[platform] navigator.storage.persist() was NOT granted — browser may evict this origin\'s data under pressure.');
    return report;
  }

  // ---- M2b multi-tab safety (M0 spec §3): Web Locks serialize writes;
  // BroadcastChannel notifies read-only tabs (a messaging primitive is NOT
  // a lock — the lock does exclusion, the channel does notification) ----
  async function withLock(name, fn) {
    if (typeof navigator !== 'undefined' && navigator.locks && navigator.locks.request) {
      return navigator.locks.request('webmm:' + name, { mode: 'exclusive' }, () => fn());
    }
    return fn();   // Node / old browsers: single-writer assumption
  }
  function createStoreChannel(dbName) {
    // the legacy singleton keeps its historical channel name (multi-tab
    // tests + old tabs keep hearing each other); per-project channels are
    // namespaced by database
    const nm = !dbName || dbName === DB_NAME ? 'webmm-platform-store' : 'webmm-store:' + dbName;
    try { return new BroadcastChannel(nm); }
    catch (e) { return null; }   // Node: no BroadcastChannel
  }

  // ---- registry adapter: one DB ('webmm-projects'), one store, one
  // record per project {id, dbName, name, createdAt, lastOpened, state,
  // approxBytes, schemaVersion} — docs/project-management.md §1 ----
  const REGISTRY_DB = 'webmm-projects';
  function createRegistryStorage(kind) {
    if (kind === 'mem') {
      const m = new Map();
      return {
        kind: 'mem',
        async get(id) { return m.get(id); },
        async put(id, val) { m.set(id, val); },
        async del(id) { m.delete(id); },
        async all() { return [...m.entries()].map(([k, v]) => ({ id: k, ...v })); },
      };
    }
    let db = null;
    function open() {
      return new Promise((resolve, reject) => {
        const req = indexedDB.open(REGISTRY_DB, 1);
        req.onupgradeneeded = () => {
          if (!req.result.objectStoreNames.contains('projects')) req.result.createObjectStore('projects');
        };
        req.onsuccess = () => { db = req.result; resolve(db); };
        req.onerror = () => reject(req.error);
      });
    }
    const req1 = (mode, fn) => new Promise((resolve, reject) => {
      const t = db.transaction('projects', mode);
      const r = fn(t.objectStore('projects'));
      r.onsuccess = () => resolve(r.result);
      r.onerror = () => reject(r.error);
    });
    return {
      kind: 'idb', open,
      async get(id) { return req1('readonly', os => os.get(id)); },
      async put(id, val) { return req1('readwrite', os => os.put(val, id)); },
      async del(id) { return req1('readwrite', os => os.delete(id)); },
      async all() {
        return new Promise((resolve, reject) => {
          const out = [];
          const r = db.transaction('projects', 'readonly').objectStore('projects').openCursor();
          r.onsuccess = () => { const c = r.result; if (c) { out.push({ id: c.key, ...c.value }); c.continue(); } else resolve(out); };
          r.onerror = () => reject(r.error);
        });
      },
    };
  }

  // scratch pad (docs §5): a global cross-project molecule clipboard, its
  // own database — same tiny one-store shape as the registry
  const SCRATCH_DB = 'webmm-scratch';
  function createScratchStorage(kind) {
    if (kind === 'mem') {
      const m = new Map();
      return { kind: 'mem', async get(k) { return m.get(k); }, async put(k, v) { m.set(k, v); },
        async del(k) { m.delete(k); }, async all() { return [...m.entries()].map(([k, v]) => ({ key: k, ...v })); } };
    }
    let db = null;
    function open() {
      return new Promise((resolve, reject) => {
        const req = indexedDB.open(SCRATCH_DB, 1);
        req.onupgradeneeded = () => { if (!req.result.objectStoreNames.contains('molecules')) req.result.createObjectStore('molecules'); };
        req.onsuccess = () => { db = req.result; resolve(db); };
        req.onerror = () => reject(req.error);
      });
    }
    const req1 = (mode, fn) => new Promise((resolve, reject) => {
      const t = db.transaction('molecules', mode);
      const r = fn(t.objectStore('molecules'));
      r.onsuccess = () => resolve(r.result);
      r.onerror = () => reject(r.error);
    });
    return { kind: 'idb', open,
      async get(k) { return req1('readonly', os => os.get(k)); },
      async put(k, v) { return req1('readwrite', os => os.put(v, k)); },
      async del(k) { return req1('readwrite', os => os.delete(k)); },
      async all() {
        return new Promise((resolve, reject) => {
          const out = [];
          const r = db.transaction('molecules', 'readonly').objectStore('molecules').openCursor();
          r.onsuccess = () => { const c = r.result; if (c) { out.push({ key: c.key, ...c.value }); c.continue(); } else resolve(out); };
          r.onerror = () => reject(r.error);
        });
      } };
  }

  return { DB_NAME, DB_VERSION, REGISTRY_DB, SCRATCH_DB, STORES, MIGRATIONS, createMemStorage, createIdbStorage, createRegistryStorage, createScratchStorage, durability, withLock, createStoreChannel };
});
