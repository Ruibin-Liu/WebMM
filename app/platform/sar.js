// platform/sar.js — Top-4: single-series SAR tools (pure, Node-testable)
// Fixed-core R-group decomposition (user-specified core; MCS is a suggestion
// source, never an automatic truth) and transparent MPO (every criterion
// shows its raw value, its desirability and its weighted contribution —
// never a lone total). All chemistry graph work is local JS; RDKit only
// re-canonicalizes the emitted fragment SMILES on the page side.
(function (root, factory) {
  if (typeof module === 'object' && module.exports) module.exports = factory();
  else { (root.Platform = root.Platform || {}).SAR = factory(); }
})(typeof self !== 'undefined' ? self : this, function () {
  'use strict';

  // ---- V2000 molblock graph (atoms: symbol + heavy neighbors; bonds: order) ----
  // Mirrors the page's existing fixed-column parsing semantics; bonds may be
  // aromatic (4) — kept as 4 and matched tolerantly.
  function parseMolGraph(molblock) {
    const lines = String(molblock || '').replace(/\r/g, '').split('\n');
    const counts = lines[3] || '';
    const na = parseInt(counts.substring(0, 3)) || 0;
    const nb = parseInt(counts.substring(3, 6)) || 0;
    if (lines.length < 4 + na + nb) return null;
    const atoms = [];
    for (let i = 0; i < na; i++) {
      atoms.push({ sym: (lines[4 + i] || '').substring(31, 34).trim(), adj: [] });
    }
    const bonds = [];
    for (let i = 0; i < nb; i++) {
      const l = lines[4 + na + i] || '';
      const a = (parseInt(l.substring(0, 3)) || 0) - 1;
      const b = (parseInt(l.substring(3, 6)) || 0) - 1;
      const o = parseInt(l.substring(6, 9)) || 1;
      if (a >= 0 && b >= 0 && a < na && b < na && a !== b) {
        bonds.push({ a, b, o });
        atoms[a].adj.push(b); atoms[b].adj.push(a);
      }
    }
    return { atoms, bonds };
  }

  // ---- subgraph matching (small cores, backtracking with bond tolerance) ----
  // aromatic bond (4) matches 1 or 2 in the target; aromatic-ring target
  // bonds (perceived by the caller as 4 in the CORE) also accept kekulé
  // 1/2 targets. Symmetric tolerance both directions.
  function bondOk(oc, ot) {
    if (oc === ot) return true;
    if (oc === 4) return ot === 1 || ot === 2;
    if (ot === 4) return oc === 1 || oc === 2;
    return false;
  }
  // Returns an array of matches; each match = Map(coreAtom -> targetAtom).
  // maxMatches caps the search (first match wins when 1 requested).
  function matchCore(core, target, maxMatches) {
    maxMatches = maxMatches || 1;
    const out = [];
    const nc = core.atoms.length, nt = target.atoms.length;
    if (nc === 0 || nc > nt) return out;
    // candidate seeds: same symbol
    const map = new Array(nc).fill(-1);
    const used = new Array(nt).fill(false);
    const cAdj = core.atoms.map(a => a.adj.slice());
    function extend(ci) {
      if (out.length >= maxMatches) return;
      if (ci === nc) {
        out.push(new Map(map.map((t, c) => [c, t])));
        return;
      }
      // pick the next core atom connected to the already-mapped part when
      // possible (connectivity pruning); first atom = seed
      let next = -1, parent = -1;
      if (ci === 0) { next = 0; }
      else {
        let bestDeg = -1;
        for (let c = 0; c < nc; c++) {
          if (map[c] >= 0) continue;
          let p = -1, deg = 0;
          for (const a of cAdj[c]) if (map[a] >= 0) { deg++; if (p < 0) p = a; }
          if (p >= 0 && deg > bestDeg) { bestDeg = deg; next = c; parent = p; }
        }
        if (next < 0) next = map.findIndex(m => m < 0), parent = -1;
      }
      const ca = core.atoms[next];
      for (let t = 0; t < nt; t++) {
        if (used[t]) continue;
        if (target.atoms[t].sym !== ca.sym) continue;
        // bond consistency with the mapped parent (and any other mapped nbr)
        let ok = true;
        if (parent >= 0) {
          const oc = core.bonds.find(bb => (bb.a === next && bb.b === parent) || (bb.b === next && bb.a === parent));
          const ot = target.bonds.find(bb => (bb.a === t && bb.b === map[parent]) || (bb.b === t && bb.a === map[parent]));
          if (!oc || !ot || !bondOk(oc.o, ot.o)) continue;
        }
        for (const a of cAdj[next]) {
          if (map[a] < 0 || a === parent) continue;
          const oc = core.bonds.find(bb => (bb.a === next && bb.b === a) || (bb.b === next && bb.a === a));
          const ot = target.bonds.find(bb => (bb.a === t && bb.b === map[a]) || (bb.b === t && bb.a === map[a]));
          if (!oc || !ot || !bondOk(oc.o, ot.o)) { ok = false; break; }
        }
        if (!ok) continue;
        // degree pruning: a core atom cannot map onto a target atom with
        // FEWER heavy neighbors (substituents are allowed, missing bonds not)
        if (target.atoms[t].adj.length < ca.adj.length) continue;
        map[next] = t; used[t] = true;
        extend(ci + 1);
        map[next] = -1; used[t] = false;
      }
    }
    extend(0);
    return out;
  }

  // ---- R-group extraction ----
  // For each core atom whose matched target atom carries MORE heavy
  // neighbors than the core atom: every extra neighbor starts an external
  // component (the substituent). Returns { rgroups: Map(coreAtomIdx ->
  // { smi, atoms: [targetIdx...] }), match } or null when no match.
  function extractRgroups(core, target) {
    const matches = matchCore(core, target, 1);
    if (!matches.length) return null;
    const m = matches[0];
    const inCore = new Set([...m.values()]);
    const rgroups = new Map();
    for (const [c, t] of m.entries()) {
      const extras = target.atoms[t].adj.filter(n => !inCore.has(n));
      if (!extras.length) continue;
      // collect the whole external component(s) attached here (follow
      // through non-core atoms only)
      const frag = new Set();
      const stack = extras.slice();
      while (stack.length) {
        const x = stack.pop();
        if (inCore.has(x) || frag.has(x)) continue;
        frag.add(x);
        for (const n of target.atoms[x].adj) if (!inCore.has(n) && !frag.has(n)) stack.push(n);
      }
      const smi = fragmentSmiles(target, frag, t);
      rgroups.set(c, { smi, atoms: [...frag], attach: t });
    }
    return { rgroups, match: m };
  }

  // ---- ring-aware kekulé SMILES writer for a fragment ----
  // The fragment's attachment to the core is emitted as the leading dummy
  // [*]. Ring edges get numbered closure labels: whichever endpoint sees the
  // edge after both are visited opens the label, the other closes it.
  // Output is kekulé (explicit = / #; aromatic input bonds land as 1/2) —
  // good enough for display and for RDKit re-canonicalization.
  function fragmentSmiles(g, fragSet, attachTarget) {
    const atoms = [...fragSet];
    if (!atoms.length) return '';
    const localOf = new Map(atoms.map((t, i) => [t, i]));   // target idx -> fragment idx
    const adjF = atoms.map(t => g.atoms[t].adj
      .filter(n => fragSet.has(n) && n !== attachTarget)
      .map(n => localOf.get(n)));
    const bondOf = (t1, t2) => {
      const bb = g.bonds.find(b => (b.a === t1 && b.b === t2) || (b.b === t1 && b.a === t2));
      return bb ? bb.o : 1;
    };
    const n = atoms.length;
    // 1) classify edges: DFS tree edges vs ring (back) edges
    const parent = new Array(n).fill(-1);
    const depth = new Array(n).fill(-1);
    const treeEdge = new Set();          // "i-j" (i<j) local
    const ringEdge = new Map();          // "i-j" -> {order}
    const ek = (i, j) => Math.min(i, j) + '-' + Math.max(i, j);
    (function dfs(u, p) {
      depth[u] = p < 0 ? 0 : depth[p] + 1;
      parent[u] = p;
      for (const v of adjF[u]) {
        if (v === p) continue;
        if (depth[v] < 0) { treeEdge.add(ek(u, v)); dfs(v, u); }
        else if (depth[v] < depth[u]) { ringEdge.set(ek(u, v), { order: bondOf(atoms[u], atoms[v]) }); }
      }
    })(0, -1);
    // 2) emission walk: ring-closure digits are emitted right after the
    //    atom symbol — the ancestor (written first) OPENS the label, the
    //    descendant CLOSES it when written
    const labels = new Map();            // edgeKey -> {label, used}
    let nextRing = 1;
    let out = '';
    const bondChar = o => (o === 2 ? '=' : o === 3 ? '#' : '');
    const symOf = t => {
      const s = g.atoms[t].sym || 'C';
      return /^(C|N|O|S|P|F|Cl|Br|I)$/.test(s) ? s : '[' + s + ']';
    };
    function ringDigits(u) {
      for (const v of adjF[u]) {
        const k = ek(u, v);
        if (!ringEdge.has(k)) continue;
        let L = labels.get(k);
        if (!L) { L = { label: nextRing++, used: false }; labels.set(k, L); out += bondChar(ringEdge.get(k).order) + L.label; }
        else if (!L.used) { L.used = true; out += bondChar(ringEdge.get(k).order) + L.label; }
      }
    }
    function visit(u) {
      out += symOf(atoms[u]);
      ringDigits(u);
      const kids = adjF[u].filter(v => treeEdge.has(ek(u, v)) && parent[v] === u);
      for (let k = 0; k < kids.length; k++) {
        const v = kids[k];
        const o = bondOf(atoms[u], atoms[v]);
        if (k > 0) out += '(';
        out += bondChar(o);
        visit(v);
        if (k > 0) out += ')';
      }
    }
    out = '[*]';
    visit(0);
    for (let i = 1; i < n; i++) if (depth[i] < 0) { out += '.'; (function r(u) { out += symOf(atoms[u]); ringDigits(u); const kids = adjF[u].filter(v => treeEdge.has(ek(u, v)) && parent[v] === u); for (let k = 0; k < kids.length; k++) { if (k > 0) out += '('; out += bondChar(bondOf(atoms[u], atoms[kids[k]])); r(kids[k]); if (k > 0) out += ')'; } })(i); }
    return out;
  }

  // ---- transparent MPO ----
  // criteria: [{key, label, lo, hi, softLo, softHi, weight, dir?}]
  // score(x): 1 inside [lo,hi], linear to 0 at soft bounds (outside soft = 0).
  // dir 'high'/'low' use half-trapezoids: high = ramp up from softLo to lo,
  // 1 above; low = 1 below hi, ramp down to softHi.
  function desirability(x, c) {
    if (x == null || !isFinite(x)) return null;
    if (c.dir === 'high') {
      if (x >= c.lo) return 1;
      if (x <= c.softLo) return 0;
      return (x - c.softLo) / (c.lo - c.softLo);
    }
    if (c.dir === 'low') {
      if (x <= c.hi) return 1;
      if (x >= c.softHi) return 0;
      return (c.softHi - x) / (c.softHi - c.hi);
    }
    if (x >= c.lo && x <= c.hi) return 1;
    if (x < c.lo) return x <= c.softLo ? 0 : (x - c.softLo) / (c.lo - c.softLo);
    return x >= c.softHi ? 0 : (c.softHi - x) / (c.softHi - c.hi);
  }
  // rows: [{molId, name, smiles, props: {MW, cLogP, ...}, activity?}]
  // Returns per row: {raw: {key: value}, scores: {key: s}, contrib: {key: w*s}, total}
  function scoreMpo(rows, criteria) {
    const wSum = criteria.reduce((s, c) => s + (c.weight || 1), 0) || 1;
    return rows.map(r => {
      const raw = {}, scores = {}, contrib = {};
      let total = 0, nScored = 0;
      for (const c of criteria) {
        const x = r.props && r.props[c.key] != null ? r.props[c.key]
          : (c.key === 'activity' ? r.activity : null);
        raw[c.key] = x;
        const s = desirability(x, c);
        scores[c.key] = s;
        if (s != null) { const cn = (c.weight || 1) * s; contrib[c.key] = cn; total += cn; nScored++; }
        else contrib[c.key] = null;
      }
      return { molId: r.molId, name: r.name, raw, scores, contrib,
        total: nScored ? total / wSum : null, nScored };
    });
  }

  const LEAD_LIKE = [
    { key: 'MW', label: 'MW', lo: 200, hi: 450, softLo: 100, softHi: 550, weight: 1 },
    { key: 'cLogP', label: 'cLogP', lo: 1, hi: 3.5, softLo: -1, softHi: 5.5, weight: 1 },
    { key: 'TPSA', label: 'TPSA', lo: 40, hi: 140, softLo: 10, softHi: 180, weight: 1 },
    { key: 'HBD', label: 'HBD', dir: 'low', hi: 2, softHi: 5, weight: 1 },
  ];

  return { parseMolGraph, matchCore, extractRgroups, fragmentSmiles, desirability, scoreMpo, LEAD_LIKE };
});
