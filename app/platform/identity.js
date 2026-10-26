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

  // ---- M5 identity layering (docs/assay-foundation.md §1) — pure, Node-testable ----
  const IK_RE = /^[A-Z]{14}-[A-Z]{10}-[A-Z]$/;
  // 'SMI:<canonical>' fallback keys (InChI generation failed) are valid
  // identity but carry no block structure
  function isInchiKey(k) { return typeof k === 'string' && IK_RE.test(k); }
  // {block1, block2, block3} or null for non-InChI keys. Block2's last two
  // chars are version/flag markers — excluded from comparison everywhere.
  function inchiKeyBlocks(k) {
    if (!isInchiKey(k)) return null;
    return { block1: k.slice(0, 14), block2: k.slice(15, 23), block3: k.slice(26, 27), full: k };
  }
  // Block2 comparable hash of a layerless molecule (no stereo / isotope /
  // charge); the full marker is UHFFFAOYSA — the last two chars are the
  // version/flag pair excluded from comparison, so the comparable form is 8
  const NO_LAYER_BLOCK2 = 'UHFFFAOY';
  function isLayerless(b) { return !!b && b.block2 === NO_LAYER_BLOCK2; }
  // Duplicate tier between two keys (either may be null → 'unknown'):
  //   'same'      full key equal (same compound)
  //   'block2'    same skeleton, different non-main layers (stereo/isotope/
  //               charge — disambiguated by the caller with stereoSpecified)
  //   'proton'    same block1+block2, different protonation layer
  //   'different' different skeleton
  function dupeTier(a, b) {
    const A = inchiKeyBlocks(a), B = inchiKeyBlocks(b);
    if (!A || !B) return 'unknown';
    if (A.full === B.full) return 'same';
    if (A.block1 !== B.block1) return 'different';
    if (A.block2 !== B.block2) return 'block2';
    return 'proton';
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

  return { hash64, molId, syntheticMolId, setKeyFn, structureKey, createRowAllocator,
    isInchiKey, inchiKeyBlocks, isLayerless, dupeTier, NO_LAYER_BLOCK2 };
});
