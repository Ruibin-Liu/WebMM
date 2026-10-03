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
    const data = { meta: new Map(), commands: new Map(), snapshots: new Map() };
    return {
      kind: 'mem',
      async get(store, key) { return data[store].get(key); },
      async put(store, key, val) { data[store].set(key, val); },
      async del(store, key) { data[store].delete(key); },
      async all(store) { return [...data[store].entries()].map(([k, v]) => ({ key: k, value: v })); },
      async count(store) { return data[store].size; },
      _data: data,
    };
  }

  // ---- IndexedDB adapter (browser) ----
  function createIdbStorage() {
    let db = null;
    function open() {
      return new Promise((resolve, reject) => {
        const req = indexedDB.open(DB_NAME, DB_VERSION);
        req.onupgradeneeded = (ev) => {
          const d = req.result;
          // run registered migrations sequentially (none at v1)
          for (const store of STORES) if (!d.objectStoreNames.contains(store)) d.createObjectStore(store);
        };
        req.onsuccess = () => { db = req.result; resolve(db); };
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
          const r = tx(store, 'readonly', os => os.get(key));
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
          const req = t.objectStore(store).put(val, key);
          req.onsuccess = () => resolve();
          req.onerror = () => reject(req.error);
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

  return { DB_NAME, DB_VERSION, STORES, MIGRATIONS, createMemStorage, createIdbStorage, durability };
});
