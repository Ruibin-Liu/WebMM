// app/agent/aza.js — headless aromatic aza/deaza enumeration extracted from
// platform.html. No DOM; RDKit/WebMM are injected.
(function (root, factory) {
  if (typeof module === 'object' && module.exports) module.exports = factory(require('./rgroup.js'));
  else { (root.Agent = root.Agent || {}).Aza = factory(root.Agent.RGroup); }
})(typeof self !== 'undefined' ? self : this, function (RGroup) {
  'use strict';

  function linesOf(text) { return text.split('\n'); }

  async function enumerateAza(options) {
    const { rdkitModule, webmm, parentSdf, parentMolblock, wantDeaza } = options || {};
    if (!rdkitModule || !webmm) throw new Error('RDKit and WebMM are required');
    const maxCandidates = options.maxCandidates === undefined ? 400 : options.maxCandidates;
    if (!Number.isInteger(maxCandidates) || maxCandidates < 1 || maxCandidates > 400) throw new Error('maxCandidates out of 1..400');
    const ensureSaTable = options.ensureSaTable || (async () => null);
    let parentSmiles = options.parentSmiles || RGroup.canonicalParentSmiles(rdkitModule, parentMolblock);
    let ringCH = [], ringN = [], surgeryBase = null, m0 = null, query = null;
    try {
      query = rdkitModule.get_qmol('[c;H1;r5,r6]');
      m0 = rdkitModule.get_mol(parentSdf);
      if (m0 && m0.is_valid()) {
        const matches = JSON.parse(m0.get_substruct_matches(query));
        if (Array.isArray(matches)) ringCH = matches.map(match => match.atoms[0]);
        // 3D parents may yield explicit-H SMILES. Rebuild the surgery base
        // from an implicit-H canonical structure before using atom indices.
        const surgerySmiles = parentSmiles.replace(/\[H\]/g, '').replace(/\(\)/g, '');
        const cleaned = rdkitModule.get_mol(surgerySmiles);
        if (cleaned && cleaned.is_valid()) {
          parentSmiles = cleaned.get_smiles();
          cleaned.delete();
          // get_molblock preserves input order, which may differ from the
          // canonical SMILES order. Re-parse canonical form before indexing.
          const heavy = rdkitModule.get_mol(parentSmiles);
          if (heavy && heavy.is_valid()) { surgeryBase = linesOf(heavy.get_molblock()); heavy.delete(); }
          else if (heavy) heavy.delete();
        } else if (cleaned) cleaned.delete();
      }
    } finally { if (query) query.delete(); }
    if (wantDeaza && m0) {
      let qn = null;
      try { qn = rdkitModule.get_qmol('[n;r5,r6;!H1]'); } catch { qn = null; }
      if (!qn) { try { qn = rdkitModule.get_qmol('[n;r5,r6]'); } catch { qn = null; } }
      if (qn) {
        try {
          const matches = JSON.parse(m0.get_substruct_matches(qn));
          if (Array.isArray(matches)) ringN = matches.map(match => match.atoms[0]);
        } catch { /* no reverse sites */ }
        qn.delete();
      }
    }
    if (m0) m0.delete();
    if (!ringCH.length && !ringN.length) throw new Error('No aromatic ring carbons (C-H) or pyridine-type nitrogens to swap');
    if (!surgeryBase) throw new Error('Could not build the heavy-atom surgery base');

    const atomCount = parseInt(surgeryBase[3].substring(0, 3));
    if (!Number.isInteger(atomCount)) throw new Error('Invalid surgery base atom count');
    const candidates = [], seen = new Set();
    const trySwap = (indices, element) => {
      const lines = surgeryBase.slice();
      for (const index of indices) {
        const line = lines[4 + index];
        lines[4 + index] = line.substring(0, 31) + element + line.substring(35);
      }
      return lines.join('\n');
    };
    const addCandidate = molblock => {
      if (candidates.length >= maxCandidates) return;
      let mol = null;
      try { mol = rdkitModule.get_mol(molblock); } catch { return; }
      if (!mol || !mol.is_valid()) { if (mol) mol.delete(); return; }
      const canonical = mol.get_smiles();
      const cleaned = mol.get_molblock();
      const descriptors = JSON.parse(mol.get_descriptors());
      mol.delete();
      if (seen.has(canonical)) return;
      seen.add(canonical);
      candidates.push({ smiles: canonical, molblock: cleaned, mw: +descriptors.amw.toFixed(2), rotb: descriptors.NumRotatableBonds ?? 0, sa: null });
    };
    const uiYield = () => new Promise(resolve => setTimeout(resolve, 0));
    let single = 0, double = 0;
    for (const index of ringCH) { single++; addCandidate(trySwap([index], 'N   ')); if (single % 20 === 0) await uiYield(); }
    for (let x = 0; x < ringCH.length; x++) for (let y = x + 1; y < ringCH.length; y++) {
      double++; addCandidate(trySwap([ringCH[x], ringCH[y]], 'N   ')); if (double % 20 === 0) await uiYield();
    }
    let deazaSingle = 0, deazaDouble = 0;
    for (const index of ringN) { deazaSingle++; addCandidate(trySwap([index], 'C   ')); if (deazaSingle % 20 === 0) await uiYield(); }
    for (let x = 0; x < ringN.length; x++) for (let y = x + 1; y < ringN.length; y++) {
      deazaDouble++; addCandidate(trySwap([ringN[x], ringN[y]], 'C   ')); if (deazaDouble % 20 === 0) await uiYield();
    }
    if (await ensureSaTable()) {
      let count = 0;
      for (const candidate of candidates) {
        try { candidate.sa = +webmm.sa_score_wasm(candidate.molblock).toFixed(2); } catch { /* null */ }
        if (++count % 20 === 0) await uiYield();
      }
    }
    return {
      candidates,
      audit: { filtered: 0, notEvaluated: 0 },
      counts: {
        n: candidates.length, single, double,
        deazaSingle: ringN.length ? deazaSingle : undefined,
        deazaDouble: ringN.length ? deazaDouble : undefined,
      },
      parentSmiles,
    };
  }

  return { enumerateAza };
});
