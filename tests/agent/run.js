// LBDD AI collaboration contract tests — Node, no browser/model/network.
'use strict';
const path = require('path');
const Capabilities = require(path.join(__dirname, '../../app/agent/capabilities.js'));
const Skills = require(path.join(__dirname, '../../app/agent/skills.js'));
const Plans = require(path.join(__dirname, '../../app/agent/plans.js'));
const Evidence = require(path.join(__dirname, '../../app/agent/evidence.js'));
const Tools = require(path.join(__dirname, '../../app/agent/tools.js'));
const Orchestrator = require(path.join(__dirname, '../../app/agent/orchestrator.js'));
const MagpieProvider = require(path.join(__dirname, '../../app/agent/magpie-provider.js'));
const Contract = require(path.join(__dirname, '../../app/tasks/contract.js'));
const State = require(path.join(__dirname, '../../app/platform/state.js'));
const Project = require(path.join(__dirname, '../../app/platform/project.js'));

let pass = 0, fail = 0;
function check(name, condition, detail) {
  if (condition) { pass++; console.log('  ✓ ' + name); }
  else { fail++; console.log('  ✗ ' + name + (detail ? ' — ' + JSON.stringify(detail).slice(0, 180) : '')); }
}
function rejects(name, fn, code) {
  try { fn(); check(name, false, 'did not throw'); }
  catch (error) { check(name, !code || error.code === code, { code: error.code, message: error.message }); }
}
async function rejectsAsync(name, fn, code) {
  try { await fn(); check(name, false, 'did not throw'); }
  catch (error) { check(name, !code || error.code === code, { code: error.code, message: error.message }); }
}
function fakeResponse(data, ok = true, status = 200) {
  return { ok, status, text: async () => JSON.stringify(data) };
}
function fakeChat(content) { return fakeResponse({ choices: [{ message: { content } }] }); }
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
    fragmentCategory: 'alkyl',
    hypothesis: 'Bounded local substitutions preserve the parent scaffold while improving the evidence profile.',
    goals: ['descriptors', 'shapeTanimoto'],
    budget: { maxCandidates: 20, refineTop: 5, finalTop: 3 },
  };
}
function baseContext() {
  return {
    revision: 'rev-8', parentRef: 'mol:lead-01', siteRef: 'attachment:2', fragmentSet: 'builtin-v1', fragmentCategory: 'alkyl',
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
  const approval = Plans.issueApproval(plan, { approvedBy: 'user:test', modelProvider: 'pi/magpie', model: 'zcode/GLM-5.3-Flash', providerRevision: 'magpie-provider-v1' });
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
  check('fragmentCategory defaults to all and rejects unknown categories',
    Plans.validatePlan((() => { const p = basePlan(); delete p.fragmentCategory; return p; })()).plan.fragmentCategory === 'all' &&
    !Plans.validatePlan({ ...basePlan(), fragmentCategory: 'unknown' }).ok);
  check('fragmentCategory participates in planHash', valid.planHash !== Plans.validatePlan({ ...basePlan(), fragmentCategory: 'aryl' }).planHash);
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
    spec.provenance.unverifiedGoals.includes('experimentalActivity') && spec.provenance.policyRev === Contract.policyRev(spec.input.stages) &&
    spec.provenance.model === 'zcode/GLM-5.3-Flash' && spec.provenance.modelProvider === 'pi/magpie');
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

  console.log('project Plan/Review/Decision persistence:');
  const persistedPlan = basePlan(), persistedHash = Plans.validatePlan(persistedPlan).planHash;
  let projectState = State.initialState();
  const reduce = cmd => { const result = State.applyCommand(projectState, cmd); if (result.error) throw new Error(result.error); projectState = result.state; return result; };
  reduce({ type: 'RecordDesignPlan', planHash: persistedHash, plan: persistedPlan, modelProvider: 'pi/magpie', model: 'zcode/GLM-5.3-Flash', createdAt: 1 });
  check('RecordDesignPlan projects a model-tagged draft', projectState.designPlans[persistedHash].status === 'draft' && projectState.designPlans[persistedHash].model.includes('GLM'));
  check('duplicate plan and approval-before-draft invariants hold',
    !!State.applyCommand(projectState, { type: 'RecordDesignPlan', planHash: persistedHash, plan: persistedPlan }).error &&
    !!State.applyCommand(State.initialState(), { type: 'ApproveDesignPlan', planHash: persistedHash, approval: {} }).error);
  const persistedApproval = Plans.issueApproval(persistedPlan, { approvedBy: 'user:test', modelProvider: 'pi/magpie', model: 'zcode/GLM-5.3-Flash' });
  reduce({ type: 'ApproveDesignPlan', planHash: persistedHash, approval: persistedApproval });
  check('approval moves draft to approved without deleting plan', projectState.designPlans[persistedHash].status === 'approved' && projectState.designPlans[persistedHash].approval.approvedBy === 'user:test');
  check('review cannot attach before task completion', !!State.applyCommand(projectState, { type: 'RecordDesignReview', reviewId: 'r1', taskId: 'task-1', planHash: persistedHash, review: {} }).error);
  reduce({ type: 'CompleteDesignPlan', planHash: persistedHash, taskId: 'task-1' });
  check('completion links the real task and blocks approval undo as dependent',
    projectState.designPlans[persistedHash].status === 'complete' && projectState.designPlans[persistedHash].taskId === 'task-1' &&
    State.hasDependents({ type: 'ApproveDesignPlan', planHash: persistedHash }, projectState));
  reduce({ type: 'RecordDesignReview', reviewId: 'r1', taskId: 'task-1', planHash: persistedHash, review: { claims: [] }, model: 'zcode/GLM-5.3-Flash', createdAt: 2 });
  const nextDraft = { ...persistedPlan, hypothesis: 'A distinct next-round hypothesis.' };
  reduce({ type: 'RecordDesignDecision', decisionId: 'd1', taskId: 'task-1', planHash: persistedHash, action: 'next', actor: 'user:local', nextPlanDraft: nextDraft, createdAt: 3 });
  check('review and next-round decision persist with the old task identity',
    projectState.designReviews.r1.taskId === 'task-1' && projectState.designDecisions.d1.nextPlanDraft.hypothesis.includes('distinct'));
  check('AI cannot record a human Decision', !!State.applyCommand(projectState, { type: 'RecordDesignDecision', decisionId: 'd2', taskId: 'task-1', planHash: persistedHash, action: 'accept', actor: 'ai:glm' }).error);
  const decisionInverse = State.inverse({ type: 'RecordDesignDecision', decisionId: 'd1' }, projectState);
  reduce(decisionInverse);
  check('decision undo hides but never deletes the fact', projectState.hiddenDesignDecisions.d1 === true && !!projectState.designDecisions.d1);
  reduce(State.inverse(decisionInverse, projectState));
  check('decision redo unhides the same fact', !projectState.hiddenDesignDecisions.d1);
  const reviewInverse = State.inverse({ type: 'RecordDesignReview', reviewId: 'r1' }, projectState);
  reduce(reviewInverse);
  check('review undo also uses the hidden overlay', projectState.hiddenDesignReviews.r1 === true && !!projectState.designReviews.r1);
  const restoredState = Project.deserializeState(Project.serializeState(projectState));
  check('snapshot round-trip preserves Plan/Review/Decision projections',
    restoredState.designPlans[persistedHash].status === 'complete' && !!restoredState.designReviews.r1 && !!restoredState.designDecisions.d1);
  const legacyState = Project.deserializeState({ libraries: {}, activities: {} });
  check('old snapshots deserialize without schema migration', legacyState.designPlans === undefined && State.initialState().schemaVersion === 1);

  console.log('evidence bundle + cited review:');
  const evidence = Evidence.makeAnalogEvidence({ spec, fact: {
    versions: { buildId: 'build-a' }, truncated: true, accounts: { B: { evaluated: 3, filtered: 0, failed: 0, notEvaluated: 0 } },
    results: [{ idx: 0, smiles: 'CC', tanimoto: 0.8, flexT: 0.85, mw: 30, rotb: 0, sa: 1 }],
  } });
  check('evidence bundle preserves real results, accounts, truncation and unresolved goals',
    evidence.resultCount === 1 && evidence.truncated && evidence.accounts.B.evaluated === 3 && evidence.unverifiedGoals.includes('experimentalActivity') &&
    evidence.plan.model === 'zcode/GLM-5.3-Flash' && evidence.plan.providerRevision === 'magpie-provider-v1');
  const review = Evidence.validateReview({
    claims: [{ text: 'One candidate has computed shape evidence; activity remains unverified.', evidenceRefs: [evidence.evidenceId, evidence.candidates[0].evidenceId] }],
    counterexamples: [evidence.candidates[0].evidenceId], unresolved: ['experimentalActivity'],
  }, evidence);
  check('review accepts claims backed by real evidence identifiers', review.ok, review.errors);
  check('review rejects fabricated evidence references', !Evidence.validateReview({ claims: [{ text: 'invented', evidenceRefs: ['evidence:fake'] }] }, evidence).ok);
  check('review rejects the summary as a candidate counterexample', !Evidence.validateReview({ claims: [], counterexamples: [evidence.evidenceId] }, evidence).ok);
  check('malformed unresolved data returns errors without throwing', !Evidence.validateReview({ claims: [], unresolved: 'activity' }, evidence).ok);
  check('malformed counterexample data returns errors without throwing', !Evidence.validateReview({ claims: [], counterexamples: {} }, evidence).ok);

  console.log('pi/magpie provider + corrective retry:');
  check('provider exposes only status/plan/review', (() => { const p = MagpieProvider.createMagpieProvider({ fetchImpl: async () => fakeResponse({}) }); return ['status', 'plan', 'review'].every(k => typeof p[k] === 'function') && !Object.keys(p).some(k => /execute|approve/i.test(k)); })());
  check('JSON parser recovers fenced output and embedded object',
    MagpieProvider.parseJsonText('```json\n{"ok":true}\n```').value.ok === true && MagpieProvider.parseJsonText('note {"x":1} tail').value.x === 1);
  const statusProvider = MagpieProvider.createMagpieProvider({ fetchImpl: async url => fakeResponse({ ok: true, model: url.endsWith('/status') ? 'zcode/GLM-5.3-Flash' : null }) });
  check('provider status uses the local pi/magpie route', (await statusProvider.status()).model === 'zcode/GLM-5.3-Flash');
  let firstPlanCalls = 0;
  const firstPlanProvider = MagpieProvider.createMagpieProvider({ fetchImpl: async () => { firstPlanCalls++; return fakeChat(JSON.stringify(basePlan())); } });
  const firstPlan = await firstPlanProvider.plan({ intent: 'test', context: baseContext(), capability, skill });
  check('valid GLM-shaped plan passes without retry', firstPlan.fragmentCategory === 'alkyl' && firstPlanCalls === 1);
  let retryPlanCalls = 0;
  const retryPlanProvider = MagpieProvider.createMagpieProvider({ fetchImpl: async (url, init) => {
    retryPlanCalls++;
    if (retryPlanCalls === 1) { const invalid = basePlan(); invalid.goals.push('experimentalActivity'); return fakeChat(JSON.stringify(invalid)); }
    const messages = JSON.parse(init.body).messages;
    check('corrective retry includes validation errors and prior model output', messages.length === 4 && /unverified/.test(messages[3].content));
    return fakeChat(JSON.stringify(basePlan()));
  } });
  const retriedPlan = await retryPlanProvider.plan({ intent: 'test', context: baseContext(), capability, skill });
  check('invalid plan is corrected once and normalized', retriedPlan.budget.maxCandidates === 20 && retryPlanCalls === 2);
  let badPlanCalls = 0;
  const badPlanProvider = MagpieProvider.createMagpieProvider({ fetchImpl: async () => { badPlanCalls++; const invalid = basePlan(); invalid.goals.push('experimentalActivity'); return fakeChat(JSON.stringify(invalid)); } });
  await rejectsAsync('two invalid plan responses fail closed', () => badPlanProvider.plan({ intent: 'test', context: baseContext(), capability, skill }), 'PROVIDER_SCHEMA');
  check('two invalid plan responses use exactly one corrective retry', badPlanCalls === 2);
  let retryReviewCalls = 0;
  const retryReviewProvider = MagpieProvider.createMagpieProvider({ fetchImpl: async () => {
    retryReviewCalls++;
    if (retryReviewCalls === 1) return fakeChat(JSON.stringify({ version: 1, claims: [] }));
    return fakeChat(JSON.stringify({ claims: [{ text: 'Evidence only.', evidenceRefs: [evidence.evidenceId] }], counterexamples: [], unresolved: ['experimentalActivity'], nextPlanDraft: null }));
  } });
  const retriedReview = await retryReviewProvider.review({ evidence, capability, skill });
  check('invalid review is corrected and retains real evidence refs', retriedReview.claims[0].evidenceRefs[0] === evidence.evidenceId && retryReviewCalls === 2);

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
