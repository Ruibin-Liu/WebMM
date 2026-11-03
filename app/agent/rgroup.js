// app/agent/rgroup.js — headless stage-A R-group enumeration extracted from
// platform.html. RDKit/WebMM are injected; no DOM or task-store access.
(function (root, factory) {
  if (typeof module === 'object' && module.exports) module.exports = factory();
  else { (root.Agent = root.Agent || {}).RGroup = factory(); }
})(typeof self !== 'undefined' ? self : this, function () {
  'use strict';

  function canonicalParentSmiles(rdkitModule, parentMolblock) {
    let mol = null;
    try {
      mol = rdkitModule.get_mol(parentMolblock);
      if (!mol || !mol.is_valid()) throw new Error('Current molecule is unparseable');
      return mol.get_smiles();
    } finally { if (mol) mol.delete(); }
  }

  async function enumerateRgroup(options) {
    const { rdkitModule, webmm, coreMolblock, parentMolblock, fragments } = options || {};
    if (!rdkitModule || !webmm) throw new Error('RDKit and WebMM are required');
    if (!coreMolblock) throw new Error('Site surgery failed (multi-attachment?)');
    if (!Array.isArray(fragments) || !fragments.length) throw new Error('No fragments in category');
    const maxCandidates = options.maxCandidates === undefined ? 400 : options.maxCandidates;
    if (!Number.isInteger(maxCandidates) || maxCandidates < 1 || maxCandidates > 400) throw new Error('maxCandidates out of 1..400');
    const ensureSaTable = options.ensureSaTable || (async () => null);

    let core = null;
    let coreSmiles;
    try {
      core = rdkitModule.get_mol(coreMolblock);
      if (!core || !core.is_valid()) throw new Error('Core (capped with R#) is unparseable');
      coreSmiles = core.get_smiles();
    } finally { if (core) core.delete(); }
    if (!coreSmiles.includes('*')) throw new Error('Core SMILES lacks the attachment dummy');

    let parent = null;
    let parentHeavy = 0;
    try {
      parent = rdkitModule.get_mol(parentMolblock);
      if (!parent || !parent.is_valid()) throw new Error('Current molecule is unparseable');
      const descriptors = JSON.parse(parent.get_descriptors());
      parentHeavy = descriptors.NumHeavyAtoms || 0;
    } finally { if (parent) parent.delete(); }

    const candidates = [];
    let enumerated = 0, filtered = 0;
    for (const [body] of fragments) {
      if (candidates.length >= maxCandidates) break;
      if (++enumerated % 20 === 0) await new Promise(resolve => setTimeout(resolve, 0));
      const reject = () => { filtered++; };
      const smiles = coreSmiles.replace('*', body);
      let mol = null;
      try { mol = rdkitModule.get_mol(smiles); } catch (error) { reject(); continue; }
      if (!mol || !mol.is_valid()) { if (mol) mol.delete(); reject(); continue; }
      const descriptors = JSON.parse(mol.get_descriptors());
      const rotb = descriptors.NumRotatableBonds ?? 0;
      const heavy = descriptors.NumHeavyAtoms || 0;
      if (rotb > 10 || heavy > parentHeavy + 9 || heavy < Math.max(1, parentHeavy - 6)) { mol.delete(); reject(); continue; }
      let canonical = mol.get_smiles();
      mol.delete();

      // The explicit-H surgical core must be reduced to the same H-free
      // molblock contract used by every conformer pipeline.
      const canonicalClean = canonical.replace(/\[H\]/g, '').replace(/\(\)/g, '');
      let cleaned = null;
      try { cleaned = rdkitModule.get_mol(canonicalClean); } catch (error) { /* filtered below */ }
      if (!cleaned || !cleaned.is_valid()) { if (cleaned) cleaned.delete(); reject(); continue; }
      canonical = cleaned.get_smiles();
      const molblock = cleaned.get_molblock();
      cleaned.delete();
      let sa = null;
      if (await ensureSaTable()) {
        try { sa = +webmm.sa_score_wasm(molblock).toFixed(2); } catch (error) { /* null is honest */ }
      }
      candidates.push({ smiles: canonical, molblock, mw: +descriptors.amw.toFixed(2), rotb, sa });
    }
    return {
      candidates,
      audit: { filtered, notEvaluated: Math.max(0, fragments.length - enumerated) },
      coreSmiles,
      parentSmiles: canonicalParentSmiles(rdkitModule, parentMolblock),
    };
  }

  return { canonicalParentSmiles, enumerateRgroup };
});
