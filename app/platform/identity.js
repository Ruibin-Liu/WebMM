// platform/identity.js — M1a identity model (M0 spec §2)
// UMD: browser global `Platform.Identity` / Node require.
(function (root, factory) {
  if (typeof module === 'object' && module.exports) module.exports = factory();
  else { (root.Platform = root.Platform || {}).Identity = factory(); }
})(typeof self !== 'undefined' ? self : this, function () {
  'use strict';

  // Deterministic 64-bit-ish hash as "hi:lo" hex (two 32-bit FNV-1a passes —
  // stable across browsers/machines, enough for identity, not cryptography).
  function fnv1a32(str, seed) {
    let h = (seed === undefined ? 0x811c9dc5 : seed) >>> 0;
    for (let i = 0; i < str.length; i++) {
      h ^= str.charCodeAt(i);
      h = Math.imul(h, 0x01000193) >>> 0;
    }
    return h >>> 0;
  }
  function hash64(s) {
    const hi = fnv1a32(s, 0x811c9dc5);
    const lo = fnv1a32(s, 0x9e3779b9);
    return hi.toString(16).padStart(8, '0') + ':' + lo.toString(16).padStart(8, '0');
  }

  // molId = hash(libraryId, record ordinal, raw record) — same-library
  // duplicate rows no longer collide (M0 spec §2).
  function molId(libraryId, ordinal, rawRecord) {
    return hash64(libraryId + '\u0000' + ordinal + '\u0000' + rawRecord);
  }

  // Synthetic-library identity (M3 hook, frozen now to avoid a migration):
  // hash(parentMolId, swapSiteIdx, groupId, groupSetVersion, canonical).
  function syntheticMolId(parentMolId, swapSiteIdx, groupId, groupSetVersion, canonical) {
    return hash64(['syn', parentMolId, swapSiteIdx, groupId, groupSetVersion, canonical].join('\u0000'));
  }

  // structureKey: project-immutable InChIKey level. The key function is
  // INJECTED (browser: the vendored RDKit chain; tests: a stub) — this
  // module stays pure and Node-testable.
  let keyFn = null;
  function setKeyFn(fn) { keyFn = fn; }
  function structureKey(molInput) {
    if (!keyFn) throw new Error('identity: structureKey fn not set');
    return keyFn(molInput);
  }

  // Append-only row-id allocator with tombstones (M0 spec §2): Merge
  // tombstones rows (never physically removes), Split allocates fresh ids,
  // compaction happens only at snapshot time and only for unreferenced rows.
  function createRowAllocator() {
    let next = 1;
    const tombstones = new Set();   // rowIds retired by Merge
    return {
      alloc() { const id = next++; tombstones.delete(id); return id; },
      tombstone(id) { tombstones.add(id); },
      isTombstoned(id) { return tombstones.has(id); },
      // snapshot-time compaction: caller supplies the referenced-id set;
      // returns the renumbering (Map oldId->newId) — tombstoned AND
      // unreferenced rows drop out, survivors renumber densely.
      compact(referencedIds) {
        const map = new Map();
        let n = 1;
        for (let id = 1; id < next; id++) {
          if (tombstones.has(id)) continue;
          if (!referencedIds || referencedIds.has(id)) map.set(id, n++);
        }
        next = n;
        tombstones.clear();
        return map;
      },
      // serialization for snapshots
      dump() { return { next, tombstones: [...tombstones] }; },
      load(d) {
        const a = createRowAllocator();
        if (d) { a._set(d.next, new Set(d.tombstones || [])); }
        return a;
      },
      _set(n, t) { next = n; tombstones.clear(); for (const x of t) tombstones.add(x); },
    };
  }

  return { hash64, molId, syntheticMolId, setKeyFn, structureKey, createRowAllocator };
});
