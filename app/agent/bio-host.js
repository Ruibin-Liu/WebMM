// app/agent/bio-host.js — bioisostere-specific identity/enumeration over the
// shared persistence, evidence, and approval infrastructure.
(function (root, factory) {
  if (typeof module === 'object' && module.exports) {
    module.exports = factory(require('../tasks/contract.js'), require('./rgroup.js'), require('./aza.js'), require('./browser-host.js'));
  } else {
    (root.Agent = root.Agent || {}).BioHost = factory(root.Tasks.Contract, root.Agent.RGroup, root.Agent.Aza, root.Agent.BrowserHost);
  }
})(typeof self !== 'undefined' ? self : this, function (Contract, RGroup, Aza, BrowserHost) {
  'use strict';

  function createBioBrowserHost(env) {
    const shared = BrowserHost.createBrowserHost(env);

    function identity() {
      const parentMolblock = env.getParentMolblock(), parentSdf = env.getParentSdf();
      if (!parentMolblock || !parentSdf) throw new Error('Load the current molecule first');
      const parentSmiles = RGroup.canonicalParentSmiles(env.rdkitModule(), parentMolblock);
      const parentRef = 'mol:' + Contract.fnv1a(Contract.canon(parentSmiles)).toString(16);
      const siteRef = 'aromatic-ring:auto';
      const fragmentSet = 'builtin-aza-v1';
      const fragmentCategory = env.getAzaDeaza() ? 'aza+deaza' : 'aza';
      const revision = 'rev-' + Contract.fnv1a(Contract.canon({ parentSdf, siteRef, fragmentSet, fragmentCategory })).toString(16);
      return {
        revision, parentRef, siteRef, fragmentSet, fragmentCategory,
        capability: 'azaScan', skillId: 'bioisostere-exploration',
        parentSmiles, parentSdf, parentMolblock,
        projectId: env.getProjectId ? env.getProjectId() : null,
        roundId: env.getRoundId ? env.getRoundId() : null,
        availableEvidence: ['descriptors', 'shapeTanimoto', 'SA'],
      };
    }

    async function getDesignContext() { return identity(); }

    async function createLocalPlan(options) {
      const current = identity();
      const goals = ['descriptors', 'shapeTanimoto', 'SA'], unverifiedGoals = [];
      if (options && options.includeActivity) { goals.push('experimentalActivity'); unverifiedGoals.push('experimentalActivity'); }
      return {
        version: 1, kind: 'bioisostereExploration', capability: 'azaScan', capabilityRevision: 1,
        skillId: 'bioisostere-exploration', skillVersion: 1,
        baseRevision: current.revision, parentRef: current.parentRef, siteRef: current.siteRef,
        fragmentSet: current.fragmentSet, fragmentCategory: current.fragmentCategory,
        hypothesis: options && options.hypothesis || 'Explore bounded aromatic aza/deaza substitutions and compare only computed evidence.',
        goals, unverifiedGoals,
        budget: options && options.budget || { maxCandidates: 40, refineTop: 10, finalTop: 5 },
      };
    }

    async function enumerate(plan) {
      const current = identity();
      if (current.revision !== plan.baseRevision || current.parentRef !== plan.parentRef || current.siteRef !== plan.siteRef || current.fragmentCategory !== plan.fragmentCategory) {
        const error = new Error('The molecule or aza direction changed after approval'); error.code = 'STALE_CONTEXT'; throw error;
      }
      const result = await Aza.enumerateAza({
        rdkitModule: env.rdkitModule(), webmm: env.webmm(),
        parentSdf: current.parentSdf, parentMolblock: current.parentMolblock,
        parentSmiles: current.parentSmiles, wantDeaza: plan.fragmentCategory === 'aza+deaza',
        ensureSaTable: env.ensureSaTable, maxCandidates: plan.budget.maxCandidates,
      });
      const buildId = await env.ensureBuildId();
      return {
        ...result,
        runCtx: { mode: 'aza', site: null, fragCat: null, deaza: plan.fragmentCategory === 'aza+deaza', parentSmiles: result.parentSmiles },
        context: { ...current, buildId, candidates: result.candidates, audit: result.audit },
      };
    }

    async function getLibrarySummary(id) {
      if (id !== 'builtin-aza-v1') return null;
      return { id, kind: 'bioisostere', categories: ['aza', 'aza+deaza'], transformations: 'single/double aromatic C-H→N and optional pyridine-type N→C-H' };
    }

    return Object.freeze({
      getDesignContext, createLocalPlan, enumerate, getLibrarySummary,
      getEvidence: shared.getEvidence,
      recordPlanDraft: shared.recordPlanDraft, recordPlanApproval: shared.recordPlanApproval,
      recordPlanComplete: shared.recordPlanComplete, recordReview: shared.recordReview,
      recordDecision: shared.recordDecision, getArtifacts: shared.getArtifacts,
    });
  }

  return { createBioBrowserHost };
});
