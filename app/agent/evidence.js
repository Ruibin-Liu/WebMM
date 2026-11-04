// app/agent/evidence.js — citation-ready evidence bundles and review checks.
// Numeric evidence is copied from real task artifacts, never inferred here.
(function (root, factory) {
  if (typeof module === 'object' && module.exports) module.exports = factory(require('./capabilities.js'));
  else { (root.Agent = root.Agent || {}).Evidence = factory(root.Agent.Capabilities); }
})(typeof self !== 'undefined' ? self : this, function (Capabilities) {
  'use strict';

  const clone = value => JSON.parse(JSON.stringify(value));
  const nonEmpty = value => typeof value === 'string' && value.trim().length > 0;

  function makeAnalogEvidence(options) {
    const { spec, fact } = options || {};
    if (!spec || spec.kind !== 'analogScan') throw new Error('analogScan spec required');
    if (!fact || !Array.isArray(fact.results)) throw new Error('analog fact with results required');
    const rootId = 'evidence:analog:' + spec.taskId;
    const capability = Capabilities.getCapability('analogScan');
    const candidates = fact.results.map(result => ({
      evidenceId: rootId + ':candidate:' + result.idx,
      candidateRef: 'candidate:' + spec.taskId + ':' + result.idx,
      idx: result.idx,
      smiles: result.smiles,
      scores: {
        shapeTanimoto: result.tanimoto == null ? null : result.tanimoto,
        flexTanimoto: result.flexT == null ? null : result.flexT,
      },
      descriptors: {
        mw: result.mw == null ? null : result.mw,
        rotb: result.rotb == null ? null : result.rotb,
        sa: result.sa == null ? null : result.sa,
      },
    }));
    return {
      evidenceId: rootId,
      kind: 'analogScanEvidence',
      taskId: spec.taskId,
      plan: {
        planHash: spec.provenance && spec.provenance.planHash || null,
        kind: spec.provenance && spec.provenance.planKind || null,
        hypothesis: spec.provenance && spec.provenance.hypothesis || null,
        baseRevision: spec.provenance && spec.provenance.baseRevision || null,
        modelProvider: spec.provenance && spec.provenance.modelProvider || null,
        model: spec.provenance && spec.provenance.model || null,
        providerRevision: spec.provenance && spec.provenance.providerRevision || null,
      },
      policyRev: spec.provenance && spec.provenance.policyRev || null,
      versions: clone(fact.versions || {}),
      accounts: clone(fact.accounts || {}),
      truncated: !!fact.truncated,
      resultCount: candidates.length,
      unverifiedGoals: clone(spec.provenance && spec.provenance.unverifiedGoals || []),
      limitations: capability.limitations.slice(),
      candidates,
    };
  }

  function evidenceIds(bundle) {
    const ids = new Set([bundle.evidenceId]);
    for (const candidate of bundle.candidates || []) ids.add(candidate.evidenceId);
    return ids;
  }

  function validateReview(raw, bundle) {
    const errors = [];
    if (!raw || typeof raw !== 'object' || Array.isArray(raw)) return { ok: false, errors: ['review must be a JSON object'], review: null };
    const allowed = new Set(['claims', 'counterexamples', 'unresolved', 'nextPlanDraft']);
    for (const key of Object.keys(raw)) if (!allowed.has(key)) errors.push('unknown review field: ' + key);
    const ids = evidenceIds(bundle);
    const candidateIds = new Set((bundle.candidates || []).map(c => c.evidenceId));
    const claims = Array.isArray(raw.claims) ? raw.claims : [];
    if (!Array.isArray(raw.claims)) errors.push('claims must be an array');
    const normalizedClaims = [];
    for (const claim of claims) {
      if (!claim || !nonEmpty(claim.text)) { errors.push('every claim needs non-empty text'); continue; }
      for (const key of Object.keys(claim)) if (key !== 'text' && key !== 'evidenceRefs') errors.push('unknown claim field: ' + key);
      if (!Array.isArray(claim.evidenceRefs) || !claim.evidenceRefs.length) { errors.push('every claim needs evidenceRefs'); continue; }
      for (const ref of claim.evidenceRefs) if (!ids.has(ref)) errors.push('unknown evidenceRef: ' + ref);
      normalizedClaims.push({ text: claim.text.trim(), evidenceRefs: claim.evidenceRefs.slice() });
    }
    const counterexamples = raw.counterexamples === undefined ? [] : raw.counterexamples;
    if (!Array.isArray(counterexamples)) errors.push('counterexamples must be an array');
    else for (const ref of counterexamples) if (!candidateIds.has(ref)) errors.push('unknown counterexample evidenceRef: ' + ref);
    const unresolved = raw.unresolved === undefined ? [] : raw.unresolved;
    if (!Array.isArray(unresolved) || unresolved.some(item => !nonEmpty(item))) errors.push('unresolved must be a string array');
    const review = {
      claims: normalizedClaims,
      counterexamples: Array.isArray(counterexamples) ? counterexamples.slice() : [],
      unresolved: Array.isArray(unresolved) ? unresolved.map(item => typeof item === 'string' ? item.trim() : item) : [],
      nextPlanDraft: raw.nextPlanDraft === undefined ? null : raw.nextPlanDraft,
    };
    return { ok: errors.length === 0, errors, review };
  }

  return { makeAnalogEvidence, evidenceIds, validateReview };
});
