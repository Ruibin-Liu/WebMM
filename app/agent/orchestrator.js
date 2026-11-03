// app/agent/orchestrator.js — two bounded model calls, no autonomous loop.
// Provider interface: plan(request) -> JSON, review(request) -> JSON.
(function (root, factory) {
  if (typeof module === 'object' && module.exports) {
    module.exports = factory(require('./capabilities.js'), require('./skills.js'), require('./plans.js'), require('./evidence.js'));
  } else {
    (root.Agent = root.Agent || {}).Orchestrator = factory(root.Agent.Capabilities, root.Agent.Skills, root.Agent.Plans, root.Agent.Evidence);
  }
})(typeof self !== 'undefined' ? self : this, function (Capabilities, Skills, Plans, Evidence) {
  'use strict';

  function createOrchestrator(options) {
    const { provider, tools } = options || {};
    if (!provider || typeof provider.plan !== 'function' || typeof provider.review !== 'function') {
      throw new Error('provider.plan and provider.review are required');
    }
    if (!tools) throw new Error('model tools required');

    async function proposePlan(options) {
      const context = options.context || await tools.getDesignContext();
      const capability = Capabilities.getCapability('analogScan');
      const skill = Skills.getSkill(options.skillId || 'rgroup-optimization');
      if (!skill) throw new Error('unknown skill: ' + options.skillId);
      const raw = await provider.plan({ intent: options.intent, context, capability, skill });
      const validation = Plans.validatePlan(raw, context);
      if (!validation.ok) return { status: 'rejected', errors: validation.errors, plan: validation.plan, planHash: validation.planHash };
      return {
        status: 'awaiting_approval',
        approvalRequired: true,
        plan: validation.plan,
        planHash: validation.planHash,
        unverifiedGoals: validation.plan.unverifiedGoals.slice(),
      };
    }

    async function reviewResults(options) {
      const evidence = await tools.getEvidence(options.taskId);
      const capability = Capabilities.getCapability('analogScan');
      const skill = Skills.getSkill(options.skillId || 'rgroup-optimization');
      if (!skill) throw new Error('unknown skill: ' + options.skillId);
      if (options.plan) {
        const planValidation = Plans.validatePlan(options.plan);
        if (!planValidation.ok || planValidation.planHash !== evidence.plan.planHash) {
          return { status: 'rejected', errors: ['review plan does not match the evidence planHash'] };
        }
      }
      const raw = await provider.review({ evidence, capability, skill });
      const validation = Evidence.validateReview(raw, evidence);
      if (!validation.ok) return { status: 'rejected', errors: validation.errors };
      const result = { status: 'reviewed', review: validation.review };
      if (validation.review.nextPlanDraft != null) {
        const context = options.context || await tools.getDesignContext();
        const next = Plans.validatePlan(validation.review.nextPlanDraft, context);
        result.nextPlan = next.ok
          ? { status: 'awaiting_approval', approvalRequired: true, plan: next.plan, planHash: next.planHash }
          : { status: 'rejected', errors: next.errors };
      }
      return result;
    }

    return Object.freeze({ proposePlan, reviewResults });
  }

  return { createOrchestrator };
});
