// LBDD AI collaboration contract tests — Node, no browser/model/network.
'use strict';
const path = require('path');
const Capabilities = require(path.join(__dirname, '../../app/agent/capabilities.js'));
const Skills = require(path.join(__dirname, '../../app/agent/skills.js'));
const Plans = require(path.join(__dirname, '../../app/agent/plans.js'));
const Evidence = require(path.join(__dirname, '../../app/agent/evidence.js'));
const Tools = require(path.join(__dirname, '../../app/agent/tools.js'));
const Orchestrator = require(path.join(__dirname, '../../app/agent/orchestrator.js'));
const Contract = require(path.join(__dirname, '../../app/tasks/contract.js'));

let pass = 0, fail = 0;
function check(name, condition, detail) {
  if (condition) { pass++; console.log('  ✓ ' + name); }
  else { fail++; console.log('  ✗ ' + name + (detail ? ' — ' + JSON.stringify(detail).slice(0, 180) : '')); }
}
function rejects(name, fn, code) {
  try { fn(); check(name, false, 'did not throw'); }
  catch (error) { check(name, !code || error.code === code, { code: error.code, message: error.message }); }
}
const clone = value => JSON.parse(JSON.stringify(value));

function basePlan() {
  return {
    version: 1,
    kind: 'rgroupOptimization',
    capability: 'analogScan',
    capabilityRevision: 1,
    skillId: 'rgroup-optimization',
    skillVersion: 1,
    baseRevision: 'rev-8',
    parentRef: 'mol:lead-01',
    siteRef: 'attachment:2',
    fragmentSet: 'builtin-v1',
    hypothesis: 'Bounded local substitutions preserve the parent scaffold while improving the evidence profile.',
    goals: ['descriptors', 'shapeTanimoto'],
    budget: { maxCandidates: 20, refineTop: 5, finalTop: 3 },
  };
}
function baseContext() {
  return {
    revision: 'rev-8', parentRef: 'mol:lead-01', siteRef: 'attachment:2', fragmentSet: 'builtin-v1',
    parentSmiles: 'CCO', parentSdf: 'PARENT-SDF', buildId: 'build-a', projectId: 'project-a', roundId: 'round-8',
    audit: { filtered: 2, notEvaluated: 1 },
    candidates: [
      { idx: 0, smiles: 'CC', molblock: 'MB0', mw: 30, rotb: 0, sa: 1 },
      { idx: 1, smiles: 'CCC', molblock: 'MB1', mw: 44, rotb: 1, sa: 1.2 },
      { idx: 2, smiles: 'CCO', molblock: 'MB2', mw: 46, rotb: 0, sa: 1.1 },
    ],
  };
}
function approvedSpec() {
  const plan = basePlan();
  plan.goals.push('experimentalActivity');
  plan.unverifiedGoals = ['experimentalActivity'];
  const approval = Plans.issueApproval(plan, { approvedBy: 'user:test' });
  return Plans.compileApprovedPlan({ plan, context: baseContext(), approval, taskId: 'task-agent-1' });
}

async function main() {
  console.log('capability registry + skill:');
  const capability = Capabilities.getCapability('analogScan');
  check('analogScan exposes revision, task kind, limits and scientific boundaries',
    capability.revision === 1 && capability.taskKind === 'analogScan' && capability.budgetLimits.maxCandidates === 400 &&
    capability.limitations.some(x => /Experimental activity/.test(x)));
  check('registry returns defensive copies',
    (() => { const a = Capabilities.getCapability('analogScan'); a.budgetLimits.maxCandidates = 1; return Capabilities.getCapability('analogScan').budgetLimits.maxCandidates === 400; })());
  const skill = Skills.getSkill('rgroup-optimization');
  check('R-group skill binds the same capability revision and review rules',
    skill.capability === 'analogScan' && skill.capabilityRevision === 1 && skill.reviewRules.some(x => /evidence identifiers/.test(x)));
  check('library summary is metadata, not a full fragment dump', Capabilities.getLibrarySummary('builtin-v1').kind === 'rgroup');

  console.log('DesignPlan normalization + scientific honesty:');
  const valid = Plans.validatePlan(basePlan(), baseContext());
  check('valid R-group plan normalizes with deterministic planHash', valid.ok && valid.plan.skillVersion === 1 && valid.planHash.startsWith('dp'), valid.errors);
  check('planHash is stable across input key order and sensitive to budget',
    valid.planHash === Plans.validatePlan(Object.assign({}, ...Object.entries(basePlan()).reverse().map(([k, v]) => ({ [k]: v })))).planHash &&
    valid.planHash !== Plans.validatePlan({ ...basePlan(), budget: { maxCandidates: 19, refineTop: 5, finalTop: 3 } }).planHash);
  const defaults = Plans.validatePlan({ ...basePlan(), budget: { maxCandidates: 8 } });
  check('small candidate budget shrinks default refinement budgets', defaults.ok && defaults.plan.budget.refineTop === 8 && defaults.plan.budget.finalTop === 8);
  const activity = basePlan(); activity.goals.push('experimentalActivity');
  check('experimentalActivity is rejected when presented as computable', !Plans.validatePlan(activity).ok);
  activity.unverifiedGoals = ['experimentalActivity'];
  const activityValid = Plans.validatePlan(activity);
  check('experimentalActivity is allowed only as an explicit unverified goal', activityValid.ok && activityValid.plan.unverifiedGoals.includes('experimentalActivity'));
  const unknownGoal = basePlan(); unknownGoal.goals.push('bindingAffinity');
  check('unknown goals are rejected rather than silently relabelled', !Plans.validatePlan(unknownGoal).ok);
  check('unknown plan fields are rejected', !Plans.validatePlan({ ...basePlan(), autoExecute: true }).ok);
  check('budget hard limits match the existing analog funnel',
    !Plans.validatePlan({ ...basePlan(), budget: { maxCandidates: 401 } }).ok &&
    !Plans.validatePlan({ ...basePlan(), budget: { refineTop: 31 } }).ok &&
    !Plans.validatePlan({ ...basePlan(), budget: { finalTop: 11 } }).ok);
  check('inconsistent refinement budget is rejected', !Plans.validatePlan({ ...basePlan(), budget: { maxCandidates: 5, refineTop: 6, finalTop: 2 } }).ok);

  console.log('approval-bound PlanCompiler:');
  const spec = approvedSpec();
  check('compiler emits the frozen analogScan input v1 with expanded stage defaults',
    spec.kind === 'analogScan' && spec.inputVersion === 1 && spec.input.stages.length === 4 &&
    spec.input.stages[1].nConfs === 10 && spec.input.stages[2].nConfs === 30 && spec.input.stages[2].seedOffset === 14200);
  check('compiler maps plan budgets to A/B/C/D without arbitrary stages',
    spec.input.stages[0].caps.products === 20 && spec.input.stages[1].cap === 20 && spec.input.stages[2].cap === 5 && spec.input.stages[3].cap === 3);
  check('compiler records plan, approval, unresolved goals and policy provenance',
    spec.provenance.planHash.startsWith('dp') && spec.provenance.approvedBy === 'user:test' &&
    spec.provenance.unverifiedGoals.includes('experimentalActivity') && spec.provenance.policyRev === Contract.policyRev(spec.input.stages));
  check('compiled spec retains analog inputRevision semantics', Contract.inputRevisionFor('analogScan', spec.input).policy === spec.provenance.policyRev);
  rejects('compile requires an approval object', () => Plans.compileApprovedPlan({ plan: basePlan(), context: baseContext(), taskId: 't' }), 'APPROVAL_MISMATCH');
  rejects('approval cannot be reused after a plan edit', () => {
    const plan = basePlan(); const approval = Plans.issueApproval(plan, { approvedBy: 'u' }); plan.budget.maxCandidates = 19;
    Plans.compileApprovedPlan({ plan, context: baseContext(), approval, taskId: 't' });
  }, 'APPROVAL_MISMATCH');
  rejects('stale context revision is rejected before execution', () => {
    const plan = basePlan(); const approval = Plans.issueApproval(plan, { approvedBy: 'u' }); const context = { ...baseContext(), revision: 'rev-9' };
    Plans.compileApprovedPlan({ plan, context, approval, taskId: 't' });
  }, 'PLAN_INVALID');
  rejects('candidate count cannot exceed the approved budget', () => {
    const plan = { ...basePlan(), budget: { maxCandidates: 2, refineTop: 2, finalTop: 1 } };
    const approval = Plans.issueApproval(plan, { approvedBy: 'u' });
    Plans.compileApprovedPlan({ plan, context: baseContext(), approval, taskId: 't' });
  }, 'BUDGET_EXCEEDED');
  rejects('compile requires complete current context identity', () => {
    const plan = basePlan(); const approval = Plans.issueApproval(plan, { approvedBy: 'u' }); const context = baseContext(); delete context.revision;
    Plans.compileApprovedPlan({ plan, context, approval, taskId: 't' });
  }, 'CONTEXT_INVALID');

  console.log('evidence bundle + cited review:');
  const evidence = Evidence.makeAnalogEvidence({ spec, fact: {
    versions: { buildId: 'build-a' }, truncated: true, accounts: { B: { evaluated: 3, filtered: 0, failed: 0, notEvaluated: 0 } },
    results: [{ idx: 0, smiles: 'CC', tanimoto: 0.8, flexT: 0.85, mw: 30, rotb: 0, sa: 1 }],
  } });
  check('evidence bundle preserves real results, accounts, truncation and unresolved goals',
    evidence.resultCount === 1 && evidence.truncated && evidence.accounts.B.evaluated === 3 && evidence.unverifiedGoals.includes('experimentalActivity'));
  const review = Evidence.validateReview({
    claims: [{ text: 'One candidate has computed shape evidence; activity remains unverified.', evidenceRefs: [evidence.evidenceId, evidence.candidates[0].evidenceId] }],
    counterexamples: [evidence.candidates[0].evidenceId], unresolved: ['experimentalActivity'],
  }, evidence);
  check('review accepts claims backed by real evidence identifiers', review.ok, review.errors);
  check('review rejects fabricated evidence references', !Evidence.validateReview({ claims: [{ text: 'invented', evidenceRefs: ['evidence:fake'] }] }, evidence).ok);
  check('review rejects the summary as a candidate counterexample', !Evidence.validateReview({ claims: [], counterexamples: [evidence.evidenceId] }, evidence).ok);
  check('malformed unresolved data returns errors without throwing', !Evidence.validateReview({ claims: [], unresolved: 'activity' }, evidence).ok);
  check('malformed counterexample data returns errors without throwing', !Evidence.validateReview({ claims: [], counterexamples: {} }, evidence).ok);

  console.log('allow-listed tools + two-call orchestrator:');
  const host = {
    getDesignContext: async () => baseContext(),
    getEvidence: async taskId => { if (taskId !== evidence.taskId) throw new Error('unknown task'); return evidence; },
  };
  const tools = Tools.createTools(host);
  check('model tools expose no execution or approval capability',
    !Object.keys(tools).some(k => /execute|approve|compile/i.test(k)));
  const modelContext = await tools.getDesignContext();
  check('model context tool drops execution payloads', modelContext.candidateCount === 3 && modelContext.candidates === undefined && modelContext.parentSdf === undefined && modelContext.buildId === undefined);
  check('validatePlan tool uses current headless context', (await tools.validatePlan(basePlan())).ok);
  check('unknown library summary is honest', (await tools.getLibrarySummary('missing')) === null);

  let planCalls = 0, reviewCalls = 0;
  const provider = {
    plan: async request => { planCalls++; check('planning provider receives intent, context, capability and skill', request.intent === 'explore local analogs' && request.context.revision === 'rev-8' && request.capability.id === 'analogScan' && request.skill.id === 'rgroup-optimization'); check('planning context is whitelist-sanitized', request.context.candidateCount === 3 && request.context.candidates === undefined && request.context.parentSdf === undefined && request.context.buildId === undefined); return basePlan(); },
    review: async request => { reviewCalls++; check('review provider receives the real evidence bundle and method skill', request.evidence.evidenceId === evidence.evidenceId && request.skill.id === 'rgroup-optimization'); return {
      claims: [{ text: 'Computed evidence is available; activity remains unresolved.', evidenceRefs: [evidence.evidenceId] }],
      counterexamples: [], unresolved: ['experimentalActivity'], nextPlanDraft: basePlan(),
    }; },
  };
  const orchestrator = Orchestrator.createOrchestrator({ provider, tools });
  const proposed = await orchestrator.proposePlan({ intent: 'explore local analogs' });
  check('planning call stops at awaiting_approval', proposed.status === 'awaiting_approval' && proposed.approvalRequired && planCalls === 1);
  const reviewed = await orchestrator.reviewResults({ taskId: evidence.taskId });
  check('review call validates citations and produces only a new approval-bound draft',
    reviewed.status === 'reviewed' && reviewed.nextPlan.status === 'awaiting_approval' && reviewed.nextPlan.approvalRequired && reviewCalls === 1);
  const badProvider = { plan: async () => { const p = basePlan(); p.goals.push('experimentalActivity'); return p; }, review: provider.review };
  const rejected = await Orchestrator.createOrchestrator({ provider: badProvider, tools }).proposePlan({ intent: 'overclaim activity' });
  check('orchestrator rejects unsupported model claims without execution', rejected.status === 'rejected' && rejected.errors.some(e => /unverified/.test(e)));
  const badReviewProvider = { plan: provider.plan, review: async () => ({ claims: [{ text: 'fake', evidenceRefs: ['evidence:fake'] }] }) };
  const badReview = await Orchestrator.createOrchestrator({ provider: badReviewProvider, tools }).reviewResults({ taskId: evidence.taskId });
  check('orchestrator rejects fabricated review citations', badReview.status === 'rejected');
  check('test flow used exactly one planning call and one review call per provider', planCalls === 1 && reviewCalls === 1);
}

main().then(() => {
  console.log('RESULT: ' + pass + ' passed, ' + fail + ' failed');
  process.exitCode = fail ? 1 : 0;
}).catch(error => { console.error(error); console.log('RESULT: ' + pass + ' passed, ' + (fail + 1) + ' failed'); process.exitCode = 1; });
