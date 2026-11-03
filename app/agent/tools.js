// app/agent/tools.js — allow-listed model tools. There is deliberately no
// execute tool: approval and compilation belong to the UI/execution layer.
(function (root, factory) {
  if (typeof module === 'object' && module.exports) {
    module.exports = factory(require('./capabilities.js'), require('./plans.js'));
  } else {
    (root.Agent = root.Agent || {}).Tools = factory(root.Agent.Capabilities, root.Agent.Plans);
  }
})(typeof self !== 'undefined' ? self : this, function (Capabilities, Plans) {
  'use strict';

  // Only the selected design context may cross the model boundary. Execution
  // payloads (SDF, candidate molblocks, build IDs) stay in the local layer.
  function sanitizeContext(context) {
    const allowed = ['revision', 'parentRef', 'siteRef', 'fragmentSet', 'fragmentCategory', 'parentSmiles', 'projectId', 'roundId', 'designGoals', 'availableEvidence', 'audit'];
    const out = {};
    for (const key of allowed) if (context[key] !== undefined) out[key] = context[key];
    if (Array.isArray(context.candidates)) out.candidateCount = context.candidates.length;
    else if (context.candidateCount !== undefined) out.candidateCount = context.candidateCount;
    return out;
  }

  function createTools(host) {
    if (!host || typeof host.getDesignContext !== 'function') throw new Error('getDesignContext host required');
    if (typeof host.getEvidence !== 'function') throw new Error('getEvidence host required');

    return Object.freeze({
      getDesignContext: async () => sanitizeContext(await host.getDesignContext()),
      listCapabilities: () => Capabilities.listCapabilities(),
      getCapability: id => Capabilities.getCapability(id),
      getLibrarySummary: async id => {
        if (typeof host.getLibrarySummary === 'function') {
          const summary = await host.getLibrarySummary(id);
          if (summary) return summary;
        }
        return Capabilities.getLibrarySummary(id);
      },
      validatePlan: async plan => Plans.validatePlan(plan, sanitizeContext(await host.getDesignContext())),
      getEvidence: taskId => host.getEvidence(taskId),
    });
  }

  return { createTools, sanitizeContext };
});
