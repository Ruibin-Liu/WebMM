// agent_glm_live.test.js — MANUAL live acceptance through pi/magpie.
// Requires: npm run agent:dev (default :8902) and a ready local magpie.
// Not registered in run.sh: it consumes a real model call and may be slow.
const { chromium } = require('/opt/homebrew/lib/node_modules/playwright');
const EXE = '/Users/rliu/Library/Caches/ms-playwright/chromium_headless_shell-1243/chrome-headless-shell-mac-arm64/chrome-headless-shell';
const URL = process.env.WEBMM_AGENT_URL || 'http://127.0.0.1:8902/app/platform.html';
const SKILL = process.env.WEBMM_AGENT_SKILL || 'rgroup-optimization';
const EXPECTED_CAPABILITY = SKILL === 'bioisostere-exploration' ? 'azaScan' : 'analogScan';
let pass = 0, fail = 0, browser;
function check(name, cond, detail) {
  if (cond) { pass++; console.log('  ✓ ' + name); }
  else { fail++; console.log('  ✗ ' + name + (detail ? ' — ' + detail : '')); }
}
(async () => {
  browser = await chromium.launch({ executablePath: EXE });
  const page = await (await browser.newContext({ viewport: { width: 1440, height: 900 } })).newPage();
  const errors = [];
  page.on('pageerror', error => errors.push(error.message.slice(0, 200)));
  await page.goto(URL, { waitUntil: 'load' });
  await page.waitForFunction(() => window.webmm && window.Agent && window.WebMMAgentTools, null, { timeout: 120000 });
  await page.evaluate(() => { document.getElementById('input').value = 'CC(=O)Oc1ccccc1O'; process(true); });
  await page.waitForFunction(() => !!document.querySelector('#svg svg'), null, { timeout: 60000 });
  await page.evaluate(() => embed3D());
  await page.waitForFunction(() => !!document.querySelector('#viewer3d canvas'), null, { timeout: 60000 });
  await page.evaluate(skill => {
    switchMode('analog'); analogDetectSites(); agentSetSkill(skill);
    if (skill === 'bioisostere-exploration') document.getElementById('agentIntent').value = '探索芳环 aza 等排替换；实验活性与选择性仅作未验证目标';
  }, SKILL);
  await page.waitForFunction(() => [...document.getElementById('analogSite').options].some(o => o.value), null, { timeout: 30000 });
  await page.evaluate(() => {
    document.getElementById('analogSite').selectedIndex = 1;
    document.getElementById('analogFragCat').value = 'alkyl';
    agentPlanConnect();
  });
  await page.waitForFunction(() => /已连接 .*GLM/.test(document.getElementById('agentPlanStatus').textContent), null, { timeout: 30000 });
  check('Plan Card connects to the configured pi/magpie GLM', await page.evaluate(() => window.WebMMAgentModel.includes('GLM')));

  await page.evaluate(() => agentPlanDraft());
  await page.waitForFunction(() => /等待本地批准|失败|被拒/.test(document.getElementById('agentPlanStatus').textContent), null, { timeout: 240000 });
  const draftStatus = await page.evaluate(() => document.getElementById('agentPlanStatus').textContent);
  check('real GLM produces a schema-valid plan (one corrective retry allowed)', draftStatus.includes('等待本地批准'), draftStatus);
  if (!draftStatus.includes('等待本地批准')) throw new Error(draftStatus);
  await page.evaluate(() => {
    document.getElementById('agentMaxCandidates').value = '4';
    document.getElementById('agentRefineTop').value = '2';
    document.getElementById('agentFinalTop').value = '1';
    agentPlanApproveAndRun();
  });
  await page.waitForFunction(() => /Done: \d+ scored|失败|取消/.test(document.getElementById('analogStatus').textContent), null, { timeout: 300000 });
  await page.waitForTimeout(1000);
  const run = await page.evaluate(async () => {
    const records = (await window.__platformCompat.getStorage().all('tasks')).map(r => r.value).filter(r => r.state === 'complete');
    const record = records[records.length - 1];
    const fact = record ? await window.__platformCompat.getStorage().get('facts', 'analog:' + record.spec.taskId) : null;
    return record ? { record, fact } : null;
  });
  check('real GLM plan completes after local approval', !!run && !!run.fact, run ? '' : 'no complete fact');
  if (run) {
    check('task and evidence retain pi/magpie model provenance',
      run.record.spec.provenance.modelProvider === 'pi/magpie' && run.record.spec.provenance.model.includes('GLM') && run.record.spec.provenance.providerRevision === 'magpie-provider-v1' && run.record.spec.provenance.planHash.startsWith('dp') && run.record.spec.provenance.capability === EXPECTED_CAPABILITY);
    check('live run uses the edited small budget', run.record.spec.input.candidates.length <= 4 && run.record.spec.input.stages[2].cap === 2 && run.record.spec.input.stages[3].cap === 1);
  }
  await page.evaluate(() => agentPlanReview());
  await page.waitForFunction(() => document.getElementById('agentReviewOut').textContent.length > 20, null, { timeout: 240000 });
  const review = await page.evaluate(() => document.getElementById('agentReviewOut').textContent);
  check('real GLM review cites actual evidence without execution', review.includes('evidence:analog:') && !/复盘被拒|复盘失败/.test(review), review.slice(0, 400));
  const persisted = await page.evaluate(hash => ({
    plan: window.__platformCompat.projectApi().state.designPlans[hash],
    reviews: Object.values(window.__platformCompat.projectApi().state.designReviews),
  }), run.record.spec.provenance.planHash);
  check('real GLM Plan/Review persist with task and model lineage', persisted.plan && persisted.plan.status === 'complete' && persisted.reviews.length === 1 && persisted.reviews[0].model === run.record.spec.provenance.model);
  check('zero page errors in the live path', errors.length === 0, errors.join(' | '));
  await browser.close();
  console.log('RESULT: ' + pass + ' passed, ' + fail + ' failed');
  process.exitCode = fail ? 1 : 0;
})().catch(async error => { console.error(error); if (browser) await browser.close(); console.log('RESULT: ' + pass + ' passed, ' + (fail + 1) + ' failed'); process.exitCode = 1; });
