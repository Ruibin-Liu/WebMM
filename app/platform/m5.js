// platform/m5.js — 结构-数据底座 pure logic (docs/assay-foundation.md)
// Everything here is deterministic and Node-testable: CSV parsing (with the
// qualifier/unicode/e-notation traps), the row→molId match ladder, and the
// duplicate-family audit grouping. RDKit-dependent steps (canonical SMILES,
// fragment parent keys) stay in the page and inject results via arguments.
(function (root, factory) {
  if (typeof module === 'object' && module.exports) module.exports = factory();
  else { (root.Platform = root.Platform || {}).M5 = factory(); }
})(typeof self !== 'undefined' ? self : this, function () {
  'use strict';

  // ---- CSV (RFC4180 + BOM + quoted fields + ; or , sniffing) ----
  function parseCsv(text) {
    text = String(text == null ? '' : text);
    if (text.charCodeAt(0) === 0xFEFF) text = text.slice(1);
    const delim = (text.split('\n')[0].match(/;/g) || []).length >
                  (text.split('\n')[0].match(/,/g) || []).length ? ';' : ',';
    const rows = []; let row = []; let cell = ''; let q = false;
    for (let i = 0; i < text.length; i++) {
      const c = text[i];
      if (q) {
        if (c === '"') { if (text[i + 1] === '"') { cell += '"'; i++; } else q = false; }
        else cell += c;
      } else if (c === '"') q = true;
      else if (c === delim) { row.push(cell); cell = ''; }
      else if (c === '\n') { row.push(cell); rows.push(row); row = []; cell = ''; }
      else if (c === '\r') { /* skip */ }
      else cell += c;
    }
    if (cell !== '' || row.length) { row.push(cell); rows.push(row); }
    return rows.filter(r => r.some(c => c.trim() !== ''));
  }

  // ---- cell → {qualifier, value} ----
  // Handles: ">10", "≥ 3.1", "<=2", "~5", "1.2e-6", "1,2" (euro), "10 µM"
  // (trailing unit tolerated only in the numeric-with-unit form), "ND",
  // "inactive", "N/A" → ND. Returns null for unparseable.
  function parseActivityCell(raw) {
    const valueRaw = String(raw == null ? '' : raw).trim();
    if (!valueRaw) return null;
    const up = valueRaw.toUpperCase();
    if (up === 'ND' || up === 'N/A' || up === 'INACTIVE' || up === 'NOT DETERMINED') {
      return { qualifier: 'ND', value: null, valueRaw };
    }
    let s = valueRaw
      .replace(/≥/g, '>=').replace(/≤/g, '<=')
      .replace(/µ/g, 'μ');                       // U+00B5 → U+03BC, one codepoint
    let qualifier = '=';
    const qm = s.match(/^(>=|<=|>|<|~)\s*(.*)$/);
    if (qm) { qualifier = qm[1]; s = qm[2].trim(); }
    // trailing unit ("8.2 µM", "8.2 uM") — one alnum/unit token after the number
    const um = s.match(/^([-+0-9.,eE]+)\s*(?:n|m|µ|μ|u|f|p)?[MmolLrK%]?\b.*$/);
    if (!um) return null;
    let num = um[1];
    // euro decimal comma: "8,2" (single comma, no dot, not a thousand sep)
    if (/^-?\d+,\d+$/.test(num)) num = num.replace(',', '.');
    else num = num.replace(/,/g, '');            // thousands "1,234.5"
    const v = Number(num);
    if (!isFinite(v)) return null;
    return { qualifier, value: v, valueRaw };
  }

  function normalizeName(s) {
    return String(s == null ? '' : s).normalize('NFKC').trim().toLowerCase()
      .replace(/[\s_-]+/g, ' ');
  }

  // ---- match ladder (docs §2 CSV 导入): tier written into the record ----
  // libIndex: { byMolId: Map, byKey: Map(full IK → [molId]), byParentKey: Map,
  //            byBlock1: Map, byName: Map }  (built by the caller)
  // Each csv row: {molId?, name?, smiles?, key?, parentKey?, block1?, value...}
  // Returns {resolved: [{molId, tier, row}], ambiguous: [...], unresolved: [...]}
  // tier: 'molId' | 'exact' | 'parent' | 'block1' | 'name' (the last four are
  // advisory matches the user must confirm in the preview step).
  function matchLadder(rows, libIndex) {
    const resolved = [], ambiguous = [], unresolved = [];
    for (const row of rows) {
      const candidates = (tier) => {
        const key = tier === 'molId' ? row.molId
          : tier === 'exact' ? row.key
          : tier === 'parent' ? row.parentKey
          : tier === 'block1' ? row.block1
          : null;
        if (!key) return [];
        const m = (tier === 'molId') ? (libIndex.byMolId.get(key) ? [key] : [])
          : (libIndex[tier === 'exact' ? 'byKey' : tier === 'parent' ? 'byParentKey' : 'byBlock1'].get(key) || []);
        return m.slice();
      };
      let tier = null, hits = [];
      for (const t of ['molId', 'exact', 'parent', 'block1']) {
        hits = candidates(t);
        if (hits.length === 1) { tier = t; break; }
        if (hits.length > 1) { ambiguous.push({ row, tier: t, hits }); break; }
      }
      if (!tier && hits.length <= 1) {
        // name fallback (normalized)
        const nm = normalizeName(row.name);
        if (nm && libIndex.byName.has(nm)) {
          const h = libIndex.byName.get(nm);
          if (h.length === 1) { tier = 'name'; hits = h; }
          else { ambiguous.push({ row, tier: 'name', hits: h }); }
        }
      }
      if (tier) resolved.push({ molId: hits[0], tier, row });
      else if (!ambiguous.some(a => a.row === row)) unresolved.push(row);
    }
    return { resolved, ambiguous, unresolved };
  }

  // ---- identity audit grouping (docs §1) ----
  // entries: [{molId, smiles, structureKey, parentKey?, stereoSpecified?}]
  // groupBy: 'full' | 'parent'. Pure: builds families where >1 member.
  // QC flags per row computed alongside (salt / stereo unspecified /
  // non-InChI key / parent-uncertain via flags on the entry).
  function identityAudit(entries, groupBy) {
    groupBy = groupBy || 'full';
    const groups = new Map();
    for (const e of entries) {
      let g = null;
      if (groupBy === 'parent' && e.parentKey) g = e.parentKey;
      else if (e.structureKey) g = e.structureKey;
      if (!g) g = 'nok:' + (e.smiles || e.molId);     // degraded keys group apart
      if (!groups.has(g)) groups.set(g, []);
      groups.get(g).push(e);
    }
    const families = [];
    for (const [key, members] of groups) {
      if (members.length < 2) continue;
      // pairwise tier classification inside the family
      const pairs = [];
      for (let i = 0; i < members.length; i++) {
        for (let j = i + 1; j < members.length; j++) {
          const a = members[i], b = members[j];
          let tier;
          if (key.startsWith('nok:')) tier = 'unknown';
          else if (groupBy === 'parent') {
            // grouped by parent: sub-classify by full key relation
            tier = require_dupe(a.structureKey, b.structureKey);
          } else tier = require_dupe(a.structureKey, b.structureKey);
          let note = '';
          if (tier === 'block2') {
            const aS = !!a.stereoSpecified, bS = !!b.stereoSpecified;
            note = aS && bS ? '立体异构对' : (aS !== bS ? '指定/未指定' : '同位素/电荷层');
          } else if (tier === 'proton') note = '质子化层';
          else if (tier === 'same') note = '同物';
          pairs.push({ a: a.molId, b: b.molId, tier, note });
        }
      }
      families.push({ key, members, pairs });
    }
    // QC flag summary
    const flags = { salt: [], stereoUnspecified: [], nonInchiKey: [], parentUncertain: [] };
    for (const e of entries) {
      if (e.salt) flags.salt.push(e.molId);
      if (e.stereoSpecified === false) flags.stereoUnspecified.push(e.molId);
      if (e.structureKey && String(e.structureKey).startsWith('SMI:')) flags.nonInchiKey.push(e.molId);
      if (e.parentUncertain) flags.parentUncertain.push(e.molId);
    }
    return { families, flags };
  }
  // injected to avoid a hard dependency; the page passes Identity.dupeTier
  let _dupe = null;
  function require_dupe(a, b) {
    if (!_dupe) throw new Error('m5: setDupeFn not configured');
    return _dupe(a, b);
  }
  function setDupeFn(fn) { _dupe = fn; }

  return { parseCsv, parseActivityCell, normalizeName, matchLadder, identityAudit, setDupeFn };
});
