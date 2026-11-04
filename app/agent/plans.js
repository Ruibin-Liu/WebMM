// app/agent/plans.js — DesignPlan validation and approval-bound compilation.
// The model may validate drafts; only the UI/execution layer should call
// issueApproval() and compileApprovedPlan().
(function (root, factory) {
  if (typeof module === 'object' && module.exports) {
    module.exports = factory(require('../tasks/contract.js'), require('./capabilities.js'), require('./skills.js'));
  } else {
    (root.Agent = root.Agent || {}).Plans = factory(root.Tasks.Contract, root.Agent.Capabilities, root.Agent.Skills);
  }
})(typeof self !== 'undefined' ? self : this, function (Contract, Capabilities, Skills) {
  'use strict';

  const PLAN_VERSION = 1;
  const FIELDS = new Set(['version', 'kind', 'capability', 'capabilityRevision', 'skillId', 'skillVersion', 'baseRevision', 'parentRef', 'siteRef', 'fragmentSet', 'fragmentCategory', 'hypothesis', 'goals', 'unverifiedGoals', 'budget']);
  const BUDGET_FIELDS = new Set(['maxCandidates', 'refineTop', 'finalTop']);
  const clone = value => JSON.parse(JSON.stringify(value));
  const nonEmpty = value => typeof value === 'string' && value.trim().length > 0;
  const integer = (value, lo, hi) => Number.isInteger(value) && value >= lo && value <= hi;

  function fail(code, message) { const error = new Error(message); error.code = code; throw error; }

  function normalizePlan(raw) {
    const errors = [];
    if (!raw || typeof raw !== 'object' || Array.isArray(raw)) return { plan: null, errors: ['plan must be a JSON object'] };
    for (const key of Object.keys(raw)) if (!FIELDS.has(key)) errors.push('unknown plan field: ' + key);

    const capability = Capabilities.getCapability(raw.capability || 'analogScan');
    const skill = Skills.getSkill(raw.skillId || 'rgroup-optimization');
    if (!capability) errors.push('unknown capability: ' + raw.capability);
    if (!skill) errors.push('unknown skill: ' + raw.skillId);
    if (raw.kind !== 'rgroupOptimization') errors.push('kind must be rgroupOptimization');
    if (raw.version !== undefined && raw.version !== PLAN_VERSION) errors.push('unsupported plan version: ' + raw.version);
    if (capability && raw.capabilityRevision !== capability.revision) errors.push('capabilityRevision must be ' + capability.revision);
    if (skill && raw.skillVersion !== undefined && raw.skillVersion !== skill.version) errors.push('skillVersion must be ' + skill.version);
    if (capability && skill && (skill.capability !== capability.id || skill.capabilityRevision !== capability.revision)) errors.push('skill and capability revisions do not match');

    for (const field of ['baseRevision', 'parentRef', 'siteRef', 'hypothesis']) {
      if (!nonEmpty(raw[field])) errors.push(field + ' must be a non-empty string');
    }
    const fragmentSet = raw.fragmentSet === undefined ? 'builtin-v1' : raw.fragmentSet;
    const fragmentCategory = raw.fragmentCategory === undefined ? 'all' : raw.fragmentCategory;
    const fragmentSetInfo = capability ? capability.fragmentSets.find(f => f.id === fragmentSet) : null;
    if (capability && !fragmentSetInfo) errors.push('unsupported fragmentSet: ' + fragmentSet);
    if (fragmentSetInfo && !fragmentSetInfo.categories.includes(fragmentCategory)) errors.push('unsupported fragmentCategory: ' + fragmentCategory);

    const goals = Array.isArray(raw.goals) ? [...new Set(raw.goals)] : [];
    if (!Array.isArray(raw.goals) || !goals.length) errors.push('goals must be a non-empty string array');
    const unverifiedGoals = raw.unverifiedGoals === undefined ? [] : (Array.isArray(raw.unverifiedGoals) ? [...new Set(raw.unverifiedGoals)] : []);
    if (!Array.isArray(raw.unverifiedGoals) && raw.unverifiedGoals !== undefined) errors.push('unverifiedGoals must be a string array');
    for (const goal of [...goals, ...unverifiedGoals]) if (!nonEmpty(goal)) errors.push('goal names must be non-empty strings');
    if (capability) {
      for (const goal of goals) {
        const supported = capability.supportedGoals.includes(goal);
        const canBeUnverified = capability.knownUnverifiedGoals.includes(goal);
        if (!supported && !canBeUnverified) errors.push('unknown goal: ' + goal);
        if (canBeUnverified && !unverifiedGoals.includes(goal)) errors.push('goal must be declared unverified: ' + goal);
      }
      for (const goal of unverifiedGoals) {
        if (!capability.knownUnverifiedGoals.includes(goal)) errors.push('not an unverified-goal type: ' + goal);
        if (!goals.includes(goal)) errors.push('unverified goal is absent from goals: ' + goal);
      }
    }

    const sourceBudget = raw.budget === undefined ? {} : raw.budget;
    if (!sourceBudget || typeof sourceBudget !== 'object' || Array.isArray(sourceBudget)) errors.push('budget must be an object');
    const budgetInput = sourceBudget && typeof sourceBudget === 'object' && !Array.isArray(sourceBudget) ? sourceBudget : {};
    for (const key of Object.keys(budgetInput)) if (!BUDGET_FIELDS.has(key)) errors.push('unknown budget field: ' + key);
    const limits = capability ? capability.budgetLimits : { maxCandidates: 400, refineTop: 30, finalTop: 10 };
    const maxCandidates = budgetInput.maxCandidates === undefined ? limits.maxCandidates : budgetInput.maxCandidates;
    const refineTop = budgetInput.refineTop === undefined ? Math.min(limits.refineTop, maxCandidates) : budgetInput.refineTop;
    const finalTop = budgetInput.finalTop === undefined ? Math.min(limits.finalTop, refineTop) : budgetInput.finalTop;
    if (!integer(maxCandidates, 1, limits.maxCandidates)) errors.push('maxCandidates out of 1..' + limits.maxCandidates);
    if (!integer(refineTop, 1, limits.refineTop)) errors.push('refineTop out of 1..' + limits.refineTop);
    if (!integer(finalTop, 1, limits.finalTop)) errors.push('finalTop out of 1..' + limits.finalTop);
    if (integer(refineTop, 1, limits.refineTop) && integer(maxCandidates, 1, limits.maxCandidates) && refineTop > maxCandidates) errors.push('refineTop cannot exceed maxCandidates');
    if (integer(finalTop, 1, limits.finalTop) && integer(refineTop, 1, limits.refineTop) && finalTop > refineTop) errors.push('finalTop cannot exceed refineTop');

    const plan = {
      version: PLAN_VERSION,
      kind: raw.kind,
      capability: raw.capability || 'analogScan',
      capabilityRevision: raw.capabilityRevision,
      skillId: raw.skillId || 'rgroup-optimization',
      skillVersion: raw.skillVersion === undefined ? (skill ? skill.version : null) : raw.skillVersion,
      baseRevision: nonEmpty(raw.baseRevision) ? raw.baseRevision.trim() : raw.baseRevision,
      parentRef: nonEmpty(raw.parentRef) ? raw.parentRef.trim() : raw.parentRef,
      siteRef: nonEmpty(raw.siteRef) ? raw.siteRef.trim() : raw.siteRef,
      fragmentSet,
      fragmentCategory,
      hypothesis: nonEmpty(raw.hypothesis) ? raw.hypothesis.trim() : raw.hypothesis,
      goals,
      unverifiedGoals,
      budget: { maxCandidates, refineTop, finalTop },
    };
    return { plan, errors };
  }

  function validatePlan(raw, context) {
    const normalized = normalizePlan(raw);
    const errors = normalized.errors.slice();
    const plan = normalized.plan;
    if (plan && context) {
      if (context.revision !== undefined && context.revision !== plan.baseRevision) errors.push('stale baseRevision: plan=' + plan.baseRevision + ' context=' + context.revision);
      if (context.parentRef !== undefined && context.parentRef !== plan.parentRef) errors.push('parentRef does not match the current context');
      if (context.siteRef !== undefined && context.siteRef !== plan.siteRef) errors.push('siteRef does not match the current context');
      if (context.fragmentSet !== undefined && context.fragmentSet !== plan.fragmentSet) errors.push('fragmentSet does not match the current context');
      if (context.fragmentCategory !== undefined && context.fragmentCategory !== plan.fragmentCategory) errors.push('fragmentCategory does not match the current context');
    }
    let hash = null;
    if (plan) {
      try { hash = 'dp' + Contract.fnv1a(Contract.canon(plan)).toString(16); }
      catch (error) { errors.push('plan is not fingerprintable: ' + error.message); }
    }
    return { ok: errors.length === 0, errors, plan, planHash: hash };
  }

  function issueApproval(plan, options) {
    const validation = validatePlan(plan);
    if (!validation.ok) fail('PLAN_INVALID', validation.errors.join('; '));
    const approvedBy = options && options.approvedBy;
    if (!nonEmpty(approvedBy)) fail('APPROVER_REQUIRED', 'approvedBy must be a non-empty string');
    return {
      planHash: validation.planHash,
      baseRevision: validation.plan.baseRevision,
      approvedBy: approvedBy.trim(),
      modelProvider: options.modelProvider || null,
      model: options.model || null,
      providerRevision: options.providerRevision || null,
    };
  }

  function compileApprovedPlan(options) {
    const { plan, context, approval, taskId } = options || {};
    const validation = validatePlan(plan, context);
    if (!validation.ok) fail('PLAN_INVALID', validation.errors.join('; '));
    if (!approval || approval.planHash !== validation.planHash) fail('APPROVAL_MISMATCH', 'approval planHash does not match the normalized plan');
    if (approval.baseRevision !== validation.plan.baseRevision) fail('APPROVAL_STALE', 'approval baseRevision does not match the plan');
    if (!nonEmpty(approval.approvedBy)) fail('APPROVER_REQUIRED', 'approvedBy must be a non-empty string');
    if (!nonEmpty(taskId)) fail('TASK_ID_REQUIRED', 'taskId must be a non-empty string');
    if (!context || !nonEmpty(context.parentSmiles) || !nonEmpty(context.parentSdf)) fail('CONTEXT_INVALID', 'parentSmiles and parentSdf are required');
    for (const field of ['revision', 'parentRef', 'siteRef', 'fragmentSet', 'fragmentCategory', 'buildId']) {
      if (!nonEmpty(context[field])) fail('CONTEXT_INVALID', field + ' is required');
    }
    if (!Array.isArray(context.candidates) || !context.candidates.length) fail('CANDIDATES_REQUIRED', 'at least one enumerated candidate is required');
    const budget = validation.plan.budget;
    if (context.candidates.length > budget.maxCandidates) fail('BUDGET_EXCEEDED', 'candidate count exceeds maxCandidates');
    const ids = new Set();
    for (const candidate of context.candidates) {
      if (!candidate || !Number.isInteger(candidate.idx) || candidate.idx < 0 || ids.has(candidate.idx)) fail('CANDIDATE_INVALID', 'candidate idx values must be unique non-negative integers');
      if (!nonEmpty(candidate.smiles) || !nonEmpty(candidate.molblock)) fail('CANDIDATE_INVALID', 'every candidate needs smiles and molblock');
      for (const field of ['mw', 'rotb', 'sa']) {
        if (candidate[field] != null && (typeof candidate[field] !== 'number' || !Number.isFinite(candidate[field]))) {
          fail('CANDIDATE_INVALID', field + ' must be finite when present');
        }
      }
      ids.add(candidate.idx);
    }
    const seed = options.seed === undefined ? 7 : options.seed;
    if (!Number.isInteger(seed)) fail('SEED_INVALID', 'seed must be an integer');
    const p = validation.plan;
    const spec = Contract.makeSpec({
      taskId,
      kind: 'analogScan',
      seed,
      input: {
        parentSmiles: context.parentSmiles,
        parentSdf: context.parentSdf,
        candidates: context.candidates.map(candidate => ({
          idx: candidate.idx,
          smiles: candidate.smiles,
          molblock: candidate.molblock,
          mw: candidate.mw == null ? null : candidate.mw,
          rotb: candidate.rotb == null ? null : candidate.rotb,
          sa: candidate.sa == null ? null : candidate.sa,
        })),
        audit: clone(context.audit || { filtered: 0, notEvaluated: 0 }),
        inputRevision: { fragmentSet: p.fragmentSet },
        stages: [
          { id: 'A', op: 'enumerate', caps: { products: budget.maxCandidates } },
          { id: 'B', op: 'score', cap: budget.maxCandidates },
          { id: 'C', op: 'score', cap: budget.refineTop },
          { id: 'D', op: 'flex', cap: budget.finalTop },
        ],
      },
      provenance: {
        projectId: context.projectId || null,
        roundId: context.roundId || null,
        molId: p.parentRef,
        buildId: context.buildId,
        planHash: validation.planHash,
        planKind: p.kind,
        capabilityRevision: p.capabilityRevision,
        baseRevision: p.baseRevision,
        parentRef: p.parentRef,
        siteRef: p.siteRef,
        fragmentCategory: p.fragmentCategory,
        skillId: p.skillId,
        skillVersion: p.skillVersion,
        hypothesis: p.hypothesis,
        unverifiedGoals: p.unverifiedGoals.slice(),
        approvedBy: approval.approvedBy,
        modelProvider: approval.modelProvider || null,
        model: approval.model || null,
        providerRevision: approval.providerRevision || null,
      },
    });
    spec.provenance.policyRev = Contract.policyRev(spec.input.stages);
    return spec;
  }

  return { PLAN_VERSION, validatePlan, issueApproval, compileApprovedPlan };
});
