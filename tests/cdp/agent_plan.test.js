// agent_plan.test.js — AI collaboration browser vertical slice.
// Uses a schema-valid mock provider; no network/model credentials.
const { chromium } = require('/opt/homebrew/lib/node_modules/playwright');
const EXE = '/Users/rliu/Library/Caches/ms-playwright/chromium_headless_shell-1243/chrome-headless-shell-mac-arm64/chrome-headless-shell';
const URL = process.env.WEBMM_AGENT_URL || 'http://localhost:8901/app/platform.html';

let pass = 0, fail = 0;
function check(name, cond, detail) {
  if (cond) { pass++; console.log('  ✓ ' + name); }
  else { fail++; console.log('  ✗ ' + name + (detail ? ' — ' + detail : '')); }
}

let browser;
(async () => {
  browser = await chromium.launch({ executablePath: EXE });
  const ctx = await browser.newContext({ viewport: { width: 1440, height: 900 } });
  const page = await ctx.newPage();
  const errors = [];
  page.on('pageerror', e => errors.push('PAGEERROR: ' + e.message.slice(0, 220)));
  page.on('console', m => { if (m.type() === 'error') errors.push('CONSOLE: ' + m.text().slice(0, 180)); });

  await page.goto(URL, { waitUntil: 'load' });
  await page.waitForFunction(() => window.webmm !== undefined && window.Agent && window.WebMMAgentTools, null, { timeout: 120000 });
  await page.evaluate(() => { indexedDB.deleteDatabase('webmm-platform'); indexedDB.deleteDatabase('webmm-projects'); });
  await page.waitForTimeout(300);
  await page.reload({ waitUntil: 'load' });
  await page.waitForFunction(() => window.webmm !== undefined && window.Agent && window.WebMMAgentTools, null, { timeout: 120000 });
  await page.evaluate(() => { document.getElementById('input').value = 'CC(=O)Oc1ccccc1O'; process(true); });
  await page.waitForFunction(() => !!document.querySelector('#svg svg'), null, { timeout: 60000 });
  await page.evaluate(() => embed3D());
  await page.waitForFunction(() => !!document.querySelector('#viewer3d canvas'), null, { timeout: 60000 });
  await page.evaluate(() => { switchMode('analog'); analogDetectSites(); });
  await page.waitForFunction(() => [...document.getElementById('analogSite').options].some(o => o.value), null, { timeout: 30000 });
  await page.evaluate(() => {
    document.getElementById('analogSite').selectedIndex = 1;
    document.getElementById('analogFragCat').value = 'alkyl';
    document.getElementById('agentMaxCandidates').value = '4';
    document.getElementById('agentRefineTop').value = '2';
    document.getElementById('agentFinalTop').value = '1';
  });

  const publicContext = await page.evaluate(async () => {
    const context = await WebMMAgentTools.getDesignContext();
    return {
      revision: context.revision, category: context.fragmentCategory,
      hasSdf: 'parentSdf' in context, hasCandidates: 'candidates' in context,
      hasBuildId: 'buildId' in context,
    };
  });
  check('headless context has current category and no execution payloads',
    publicContext.revision.startsWith('rev-') && publicContext.category === 'alkyl' && !publicContext.hasSdf && !publicContext.hasCandidates && !publicContext.hasBuildId,
    JSON.stringify(publicContext));

  await page.evaluate(() => {
    window.__providerPlanCalls = 0; window.__providerReviewCalls = 0;
    window.WebMMAgentProvider = {
      plan: async request => {
        window.__providerPlanCalls++;
        if (request.context.parentSdf || request.context.candidates || request.context.buildId) throw new Error('unsanitized model context');
        const c = request.context;
        return {
          version: 1, kind: 'rgroupOptimization', capability: 'analogScan', capabilityRevision: 1,
          skillId: 'rgroup-optimization', skillVersion: 1,
          baseRevision: c.revision, parentRef: c.parentRef, siteRef: c.siteRef,
          fragmentSet: c.fragmentSet, fragmentCategory: c.fragmentCategory,
          hypothesis: 'Mock GLM plan: bounded local analog exploration with explicit uncertainty.',
          goals: ['descriptors', 'shapeTanimoto', 'SA', 'experimentalActivity'],
          unverifiedGoals: ['experimentalActivity'],
          budget: { maxCandidates: 4, refineTop: 2, finalTop: 1 },
        };
      },
      review: async request => {
        window.__providerReviewCalls++;
        return {
          claims: [{ text: 'Mock GLM review: computed evidence is available; activity is not claimed.', evidenceRefs: [request.evidence.evidenceId] }],
          counterexamples: [], unresolved: ['experimentalActivity'], nextPlanDraft: null,
        };
      },
    };
  });
  await page.evaluate(() => agentPlanDraft());
  await page.waitForFunction(() => /等待本地批准/.test(document.getElementById('agentPlanStatus').textContent), null, { timeout: 30000 });
  const draft = await page.evaluate(() => ({
    category: agentDraftPlan.fragmentCategory,
    budget: agentDraftPlan.budget,
    calls: window.__providerPlanCalls,
  }));
  check('provider draft binds category and approved budgets', draft.category === 'alkyl' && draft.budget.maxCandidates === 4 && draft.budget.refineTop === 2 && draft.budget.finalTop === 1 && draft.calls === 1, JSON.stringify(draft));

  const staleCategory = await page.evaluate(async () => {
    document.getElementById('analogFragCat').value = 'aryl';
    try { await WebMMAgentHost.enumerate(agentDraftPlan); return 'allowed'; }
    catch (error) { return error.code || error.message; }
    finally { document.getElementById('analogFragCat').value = 'alkyl'; }
  });
  check('changing fragment category after approval is stale and cannot enumerate', staleCategory === 'STALE_CONTEXT', staleCategory);

  await page.evaluate(() => {
    document.getElementById('agentHypothesis').value += ' · edited-before-approval';
    agentPlanApproveAndRun();
  });
  await page.waitForFunction(() => /Done: \d+ scored/.test(document.getElementById('analogStatus').textContent), null, { timeout: 180000 });
  await page.waitForFunction(() => !document.getElementById('agentReviewBtn').disabled, null, { timeout: 30000 });
  await page.waitForTimeout(800);
  const run = await page.evaluate(async () => {
    const storage = window.__platformCompat.getStorage();
    const records = (await storage.all('tasks')).map(r => r.value).filter(r => r.state === 'complete');
    const record = records[records.length - 1];
    const fact = record ? await storage.get('facts', 'analog:' + record.spec.taskId) : null;
    return record ? { record, fact } : null;
  });
  check('approved plan completes through the existing task executor', !!run && !!run.fact, run ? '' : 'no completed run');
  if (run) {
    const spec = run.record.spec;
    check('spec preserves plan/approval/category provenance and edited hypothesis',
      spec.provenance.planHash.startsWith('dp') && spec.provenance.approvedBy === 'user:local' && spec.provenance.fragmentCategory === 'alkyl' && /edited-before-approval/.test(spec.provenance.hypothesis));
    check('approved budgets cap candidates and C/D refinement',
      spec.input.candidates.length <= 4 && spec.input.stages[1].cap === 4 && spec.input.stages[2].cap === 2 && spec.input.stages[3].cap === 1,
      JSON.stringify({ candidates: spec.input.candidates.length, stages: spec.input.stages }));
    check('real evidence points back to the approved plan', run.fact.results.length > 0 && run.fact.accounts && spec.provenance.unverifiedGoals.includes('experimentalActivity'));
  }

  await page.evaluate(() => agentPlanReview());
  await page.waitForFunction(() => /Mock GLM review/.test(document.getElementById('agentReviewOut').textContent), null, { timeout: 30000 });
  const review = await page.evaluate(() => ({ text: document.getElementById('agentReviewOut').textContent, calls: window.__providerReviewCalls }));
  check('review uses real evidence references and leaves activity unresolved', review.text.includes('evidence:analog:') && review.text.includes('experimentalActivity') && review.calls === 1, review.text.slice(0, 220));

  // Legacy explorer must use the same extracted enumeration but no plan approval.
  await page.evaluate(() => { analogStatus(''); runAnalogExplorer(); });
  await page.waitForFunction(() => /Done: \d+ scored/.test(document.getElementById('analogStatus').textContent), null, { timeout: 180000 });
  await page.waitForTimeout(800);
  const legacy = await page.evaluate(async () => {
    const records = (await window.__platformCompat.getStorage().all('tasks')).map(r => r.value).filter(r => r.state === 'complete');
    return records[records.length - 1].spec;
  });
  check('legacy explorer still works without plan provenance', legacy.kind === 'analogScan' && !legacy.provenance.planHash && legacy.input.candidates.length > 0);
  check('zero page errors across provider, approved run, review and legacy run', errors.length === 0, errors.slice(0, 3).join(' | '));

  await browser.close();
  console.log('RESULT: ' + pass + ' passed, ' + fail + ' failed');
  process.exitCode = fail ? 1 : 0;
})().catch(async error => { console.error(error); if (browser) await browser.close(); console.log('RESULT: ' + pass + ' passed, ' + (fail + 1) + ' failed'); process.exitCode = 1; });
