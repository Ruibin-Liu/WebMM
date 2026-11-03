// app/agent/capabilities.js — model-visible WebMM capability registry.
// Pure metadata: no DOM, no execution, no model-specific protocol.
(function (root, factory) {
  if (typeof module === 'object' && module.exports) module.exports = factory();
  else { (root.Agent = root.Agent || {}).Capabilities = factory(); }
})(typeof self !== 'undefined' ? self : this, function () {
  'use strict';

  const ANALOG_SCAN = {
    id: 'analogScan',
    revision: 1,
    taskKind: 'analogScan',
    planKind: 'rgroupOptimization',
    description: 'Enumerate local analogs around one parent and evaluate them through the deterministic A/B/C/D funnel.',
    requiredContext: ['parentRef', 'siteRef', 'parentSmiles', 'parentSdf', 'candidates', 'fragmentSet'],
    supportedGoals: ['descriptors', 'shapeTanimoto', 'SA'],
    knownUnverifiedGoals: ['experimentalActivity', 'solubility', 'selectivity'],
    fragmentSets: [{ id: 'builtin-v1', kind: 'rgroup', description: 'Built-in categorized R-group fragments; the host supplies category details.' }],
    budgetLimits: { maxCandidates: 400, refineTop: 30, finalTop: 10 },
    stages: [
      { id: 'A', operation: 'enumerate-and-filter' },
      { id: 'B', operation: 'coarse-conformer-shape-score' },
      { id: 'C', operation: 'refine-top-candidates' },
      { id: 'D', operation: 'flex-top-candidates' },
    ],
    limitations: [
      'Experimental activity is not predicted by this capability.',
      'Shape similarity is not evidence of equal or improved activity.',
      'MMFF energies must not be compared as a cross-molecule activity ranking.',
      'A receptor-free run makes no binding-affinity claim.',
    ],
  };

  const REGISTRY = { analogScan: ANALOG_SCAN };
  const clone = value => JSON.parse(JSON.stringify(value));

  function listCapabilities() { return Object.keys(REGISTRY).sort().map(id => clone(REGISTRY[id])); }
  function getCapability(id) { return REGISTRY[id] ? clone(REGISTRY[id]) : null; }
  function getLibrarySummary(id) {
    const capability = getCapability('analogScan');
    return capability.fragmentSets.find(f => f.id === id) || null;
  }

  return { listCapabilities, getCapability, getLibrarySummary };
});
