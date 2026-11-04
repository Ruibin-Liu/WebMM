// bio_agent_plan.test.js — second LBDD Skill: aza/deaza bioisostere.
const { chromium } = require('/opt/homebrew/lib/node_modules/playwright');
const EXE = '/Users/rliu/Library/Caches/ms-playwright/chromium_headless_shell-1243/chrome-headless-shell-mac-arm64/chrome-headless-shell';
const URL = process.env.WEBMM_AGENT_URL || 'http://localhost:8901/app/platform.html';
let pass = 0, fail = 0, browser;
function check(name, cond, detail) {
  if (cond) { pass++; console.log('  ✓ ' + name); }
  else { fail++; console.log('  ✗ ' + name + (detail ? ' — ' + detail : '')); }
}
async function processAndEmbed(page, smiles) {
  await page.evaluate(s => { document.getElementById('input').value = s; process(true); }, smiles);
  await page.waitForFunction(() => !!document.querySelector('#svg svg'), null, { timeout: 60000 });
  await page.evaluate(() => embed3D());
  await page.waitForFunction(() => !!document.querySelector('#viewer3d canvas'), null, { timeout: 60000 });
}
(async () => {
  browser = await chromium.launch({ executablePath: EXE });
  const page = await (await browser.newContext({ viewport: { width: 1440, height: 900 } })).newPage();
  const errors = [];
  page.on('pageerror', e => errors.push(e.message.slice(0, 220)));
  page.on('console', m => { if (m.type() === 'warning' || m.type() === 'error') console.log('BROWSER', m.type(), m.text().slice(0, 300)); });
  await page.goto(URL, { waitUntil: 'load' });
  await page.waitForFunction(() => window.webmm && window.Agent && window.WebMMAgentTools, null, { timeout: 120000 });
  await processAndEmbed(page, 'CC(=O)Oc1ccccc1O');
  await page.evaluate(() => { switchMode('analog'); analogDetectSites(); agentSetSkill('bioisostere-exploration'); });
  await page.waitForFunction(() => [...document.getElementById('analogSite').options].some(o => o.value), null, { timeout: 30000 });
  await page.evaluate(() => {
    document.getElementById('azaDeaza').checked = false;
    window.WebMMAgentProvider = {
      plan: async request => {
        if (request.capability.id !== 'azaScan' || request.skill.id !== 'bioisostere-exploration') throw new Error('wrong skill route');
        const c = request.context;
        return {
          version: 1, kind: 'bioisostereExploration', capability: 'azaScan', capabilityRevision: 1,
          skillId: 'bioisostere-exploration', skillVersion: 1,
          baseRevision: c.revision, parentRef: c.parentRef, siteRef: c.siteRef,
          fragmentSet: c.fragmentSet, fragmentCategory: c.fragmentCategory,
          hypothesis: 'Mock bio plan: bounded aromatic aza substitutions.',
          goals: ['descriptors', 'shapeTanimoto', 'SA', 'experimentalActivity'], unverifiedGoals: ['experimentalActivity'],
          budget: { maxCandidates: 4, refineTop: 2, finalTop: 1 },
        };
      },
      review: async request => ({
        claims: [{ text: 'Mock bio review: aza evidence only; no activity or selectivity claim.', evidenceRefs: [request.evidence.evidenceId] }],
        counterexamples: [], unresolved: ['experimentalActivity', 'selectivity'], nextPlanDraft: null,
      }),
    };
    agentPlanDraft();
  });
  await page.waitForFunction(() => /草稿已保存，等待本地批准/.test(document.getElementById('agentPlanStatus').textContent), null, { timeout: 30000 });
  const draft = await page.evaluate(() => ({ kind: agentDraftPlan.kind, capability: agentDraftPlan.capability, category: agentDraftPlan.fragmentCategory }));
  check('bio draft uses independent kind/capability/direction', draft.kind === 'bioisostereExploration' && draft.capability === 'azaScan' && draft.category === 'aza', JSON.stringify(draft));

  await page.evaluate(() => agentPlanApproveAndRun());
  await page.waitForFunction(() => /Done: \d+ scored/.test(document.getElementById('analogStatus').textContent), null, { timeout: 180000 });
  await page.waitForTimeout(800);
  const run = await page.evaluate(async () => {
    const records = (await window.__platformCompat.getStorage().all('tasks')).map(r => r.value).filter(r => r.state === 'complete').sort((a, b) => (a.submittedAt || 0) - (b.submittedAt || 0));
    const record = records[records.length - 1];
    return { record, fact: record ? await window.__platformCompat.getStorage().get('facts', 'analog:' + record.spec.taskId) : null, mode: window.__analogRunCtx && window.__analogRunCtx.mode, counts: window.__azaCounts };
  });
  check('approved bio plan completes through the shared analog executor', !!run.record && !!run.fact, run.record ? '' : 'no run');
  if (run.record) {
    const spec = run.record.spec;
    check('bio spec retains capability and approved budget',
      spec.kind === 'analogScan' && spec.provenance.capability === 'azaScan' && spec.provenance.planKind === 'bioisostereExploration' && spec.provenance.fragmentCategory === 'aza' && spec.input.candidates.length <= 4 && spec.input.stages[2].cap === 2 && spec.input.stages[3].cap === 1);
    check('aza run has real generated candidates and counts', run.fact.results.length > 0 && run.fact.results.every(r => r.smiles.includes('n')) && run.mode === 'aza' && run.counts.n > 0,
      JSON.stringify({ mode: run.mode, counts: run.counts, results: run.fact.results.map(r => r.smiles) }));
  }
  await page.evaluate(() => agentPlanReview());
  await page.waitForFunction(() => /Mock bio review/.test(document.getElementById('agentReviewOut').textContent), null, { timeout: 30000 });
  const reviewArtifacts = await page.evaluate(() => Object.values(window.__platformCompat.projectApi().state.designReviews));
  check('bio review persists with unresolved activity/selectivity', reviewArtifacts.length === 1 && reviewArtifacts[0].review.unresolved.includes('selectivity'));

  await page.evaluate(() => { document.getElementById('agentDecisionReason').value = 'Next aza round with explicit experimental uncertainty.'; agentDecision('next'); });
  await page.waitForFunction(() => /草稿已保存，等待本地批准/.test(document.getElementById('agentPlanStatus').textContent), null, { timeout: 30000 });
  const next = await page.evaluate(() => ({ plan: agentDraftPlan, decisions: Object.values(window.__platformCompat.projectApi().state.designDecisions) }));
  check('bio next Decision creates another bio draft without execution', next.plan.kind === 'bioisostereExploration' && next.plan.capability === 'azaScan' && next.decisions.length === 1 && next.decisions[0].nextPlanDraft.kind === 'bioisostereExploration');

  // Legacy deaza proves the extracted headless module preserves reverse N→CH.
  await processAndEmbed(page, 'Cc1ccncc1');
  await page.evaluate(() => { document.getElementById('azaDeaza').checked = true; runAzaScan(); });
  await page.waitForFunction(() => /Done: \d+ scored/.test(document.getElementById('analogStatus').textContent), null, { timeout: 180000 });
  await page.waitForFunction(async () => {
    const records = (await window.__platformCompat.getStorage().all('tasks')).map(r => r.value);
    return records.some(r => r.state === 'complete' && !r.spec.provenance.planHash);
  }, null, { timeout: 30000 });
  const legacy = await page.evaluate(async () => {
    const all = (await window.__platformCompat.getStorage().all('tasks')).map(r => r.value);
    const records = all.filter(r => r.state === 'complete').sort((a, b) => (a.submittedAt || 0) - (b.submittedAt || 0));
    return { spec: records[records.length - 1].spec, rows: [...document.querySelectorAll('#analogRows .sm')].map(e => e.title || e.textContent), counts: window.__azaCounts, all: all.map(r => ({ id: r.spec.taskId, state: r.state, submittedAt: r.submittedAt, planHash: r.spec.provenance.planHash || null })) };
  });
  check('legacy deaza still produces N-free and aza products with no plan provenance',
    legacy.rows.some(x => /^Cc1ccccc1$/.test(x)) && legacy.counts.deazaSingle === 1 && !legacy.spec.provenance.planHash,
    JSON.stringify({ rows: legacy.rows.slice(0, 4), counts: legacy.counts, provenance: legacy.spec.provenance, all: legacy.all }));
  check('zero page errors across bio Agent and legacy deaza', errors.length === 0, errors.slice(0, 3).join(' | '));

  await browser.close();
  console.log('RESULT: ' + pass + ' passed, ' + fail + ' failed');
  process.exitCode = fail ? 1 : 0;
})().catch(async error => { console.error(error); if (browser) await browser.close(); console.log('RESULT: ' + pass + ' passed, ' + (fail + 1) + ' failed'); process.exitCode = 1; });
