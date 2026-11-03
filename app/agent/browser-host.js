// app/agent/browser-host.js — headless browser adapter over the current
// R-group workspace. No DOM reads or writes; all state access is injected.
(function (root, factory) {
  if (typeof module === 'object' && module.exports) {
    module.exports = factory(require('../tasks/contract.js'), require('./rgroup.js'), require('./evidence.js'));
  } else {
    (root.Agent = root.Agent || {}).BrowserHost = factory(root.Tasks.Contract, root.Agent.RGroup, root.Agent.Evidence);
  }
})(typeof self !== 'undefined' ? self : this, function (Contract, RGroup, Evidence) {
  'use strict';

  function createBrowserHost(env) {
    const required = ['rdkitModule', 'webmm', 'getParentMolblock', 'getParentSdf', 'getSite', 'getFragmentCategory', 'getFragments', 'coreMolblock', 'ensureSaTable', 'ensureBuildId', 'getTaskRecord', 'getFact'];
    for (const key of required) if (!env || typeof env[key] !== 'function') throw new Error('browser host requires ' + key);

    function identity(categoryOverride) {
      const parentMolblock = env.getParentMolblock();
      const parentSdf = env.getParentSdf();
      const site = env.getSite();
      const fragmentCategory = categoryOverride || env.getFragmentCategory();
      if (!parentMolblock || !parentSdf) throw new Error('Load the current molecule and detect sites first');
      if (!site) throw new Error('Choose an attachment site first');
      const parentSmiles = RGroup.canonicalParentSmiles(env.rdkitModule(), parentMolblock);
      const parentRef = 'mol:' + Contract.fnv1a(Contract.canon(parentSmiles)).toString(16);
      const siteRef = 'attachment:' + site;
      const fragmentSet = 'builtin-v1';
      const revision = 'rev-' + Contract.fnv1a(Contract.canon({ parentSdf, siteRef, fragmentSet, fragmentCategory })).toString(16);
      return {
        revision, parentRef, siteRef, fragmentSet, fragmentCategory,
        parentSmiles, parentSdf, parentMolblock,
        projectId: env.getProjectId ? env.getProjectId() : null,
        roundId: env.getRoundId ? env.getRoundId() : null,
        availableEvidence: ['descriptors', 'shapeTanimoto', 'SA'],
      };
    }

    async function getDesignContext() { return identity(); }

    async function createLocalPlan(options) {
      const current = identity();
      const goals = ['descriptors', 'shapeTanimoto', 'SA'];
      const unverifiedGoals = [];
      if (options && options.includeActivity) { goals.push('experimentalActivity'); unverifiedGoals.push('experimentalActivity'); }
      return {
        version: 1,
        kind: 'rgroupOptimization',
        capability: 'analogScan',
        capabilityRevision: 1,
        skillId: 'rgroup-optimization',
        skillVersion: 1,
        baseRevision: current.revision,
        parentRef: current.parentRef,
        siteRef: current.siteRef,
        fragmentSet: current.fragmentSet,
        fragmentCategory: current.fragmentCategory,
        hypothesis: options && options.hypothesis || 'Explore a bounded local R-group series and compare only the evidence produced by WebMM.',
        goals,
        unverifiedGoals,
        budget: options && options.budget || { maxCandidates: 40, refineTop: 10, finalTop: 5 },
      };
    }

    async function enumerate(plan) {
      const current = identity();
      if (current.revision !== plan.baseRevision || current.parentRef !== plan.parentRef || current.siteRef !== plan.siteRef || current.fragmentCategory !== plan.fragmentCategory) {
        const error = new Error('The molecule, site, or fragment category changed after approval'); error.code = 'STALE_CONTEXT'; throw error;
      }
      const fragments = env.getFragments().filter(fragment => plan.fragmentCategory === 'all' || fragment[1] === plan.fragmentCategory);
      const core = env.coreMolblock(plan.siteRef.slice('attachment:'.length).split(',').map(Number));
      const result = await RGroup.enumerateRgroup({
        rdkitModule: env.rdkitModule(), webmm: env.webmm(), coreMolblock: core,
        parentMolblock: current.parentMolblock, fragments,
        ensureSaTable: env.ensureSaTable, maxCandidates: plan.budget.maxCandidates,
      });
      const buildId = await env.ensureBuildId();
      return {
        ...result,
        runCtx: { mode: 'swap', site: plan.siteRef.slice('attachment:'.length), fragCat: plan.fragmentCategory, deaza: false, parentSmiles: result.parentSmiles },
        context: { ...current, buildId, candidates: result.candidates, audit: result.audit },
      };
    }

    async function getEvidence(taskId) {
      const record = await env.getTaskRecord(taskId);
      const fact = await env.getFact('analog:' + taskId);
      if (!record || !record.spec || !fact) throw new Error('No complete analog evidence for task ' + taskId);
      return Evidence.makeAnalogEvidence({ spec: record.spec, fact });
    }

    async function getLibrarySummary(id) {
      if (id !== 'builtin-v1') return null;
      const fragments = env.getFragments();
      const categories = {};
      for (const fragment of fragments) categories[fragment[1]] = (categories[fragment[1]] || 0) + 1;
      return { id, kind: 'rgroup', count: fragments.length, categories };
    }

    return Object.freeze({ getDesignContext, createLocalPlan, enumerate, getEvidence, getLibrarySummary });
  }

  return { createBrowserHost };
});
