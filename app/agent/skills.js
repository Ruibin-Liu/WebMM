// app/agent/skills.js — versioned LBDD method guidance for plan/review calls.
// Scientific rules are duplicated in validators where enforcement is required.
(function (root, factory) {
  if (typeof module === 'object' && module.exports) module.exports = factory();
  else { (root.Agent = root.Agent || {}).Skills = factory(); }
})(typeof self !== 'undefined' ? self : this, function () {
  'use strict';

  const RGROUP_OPTIMIZATION = {
    id: 'rgroup-optimization',
    version: 1,
    planKind: 'rgroupOptimization',
    capability: 'analogScan',
    capabilityRevision: 1,
    whenToUse: 'One parent has a defined attachment site and the user wants a bounded local analog exploration.',
    requiredContext: ['parentRef', 'siteRef', 'baseRevision', 'designGoals', 'availableEvidence'],
    planRules: [
      'Keep the parent and attachment identity explicit.',
      'Apply deterministic sanity and property filters before 3D work.',
      'Use coarse conformer/shape scoring before refinement.',
      'Preserve the approved candidate and refinement budgets.',
      'Declare every requested goal that the capability cannot evaluate as unverified.',
    ],
    reviewRules: [
      'Attach every claim to real evidence identifiers.',
      'Include counterexamples and unresolved goals.',
      'Do not equate shape similarity with experimental activity.',
      'Do not compare MMFF energies across different molecules as an activity ranking.',
      'A next-round plan is a new draft requiring independent approval.',
    ],
  };

  const SKILLS = { 'rgroup-optimization': RGROUP_OPTIMIZATION };
  const clone = value => JSON.parse(JSON.stringify(value));

  function listSkills() { return Object.keys(SKILLS).sort().map(id => clone(SKILLS[id])); }
  function getSkill(id) { return SKILLS[id] ? clone(SKILLS[id]) : null; }

  return { listSkills, getSkill };
});
