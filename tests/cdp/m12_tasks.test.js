// m12_tasks.test.js — Top-5 unified task layer acceptance
// (docs/task-contract.md). Covers: task lifecycle through the analogScan
// migration (spec at submit, Done = command+facts+task-complete in one
// transaction, 4-way stage accounting), cancellation (no partial facts),
// reload reconciliation (abandoned / ghost-complete), the minimal task tray,
// and the determinism red line: the same spec at pool size 1 vs 4 gives
// BITWISE-identical results.
const { chromium } = require('/opt/homebrew/lib/node_modules/playwright');
const EXE = '/Users/rliu/Library/Caches/ms-playwright/chromium_headless_shell-1243/chrome-headless-shell-mac-arm64/chrome-headless-shell';
const URL = 'http://localhost:8901/app/platform.html';

let pass = 0, fail = 0;
function check(name, cond, detail) {
  if (cond) { pass++; console.log('  ✓ ' + name); }
  else { fail++; console.log('  ✗ ' + name + (detail ? ' — ' + detail : '')); }
}

async function setup(page, poolSize) {
  await page.goto(URL, { waitUntil: 'load' });
  await page.waitForFunction(() => window.webmm !== undefined, null, { timeout: 120000 });
  await page.evaluate(sz => {
    indexedDB.deleteDatabase('webmm-platform'); indexedDB.deleteDatabase('webmm-projects');
    if (sz) window.__analogPoolSize = sz;
  }, poolSize);
  await page.waitForTimeout(300);   // let the delete requests register before unload
  await page.reload({ waitUntil: 'load' });
  await page.waitForFunction(() => window.webmm !== undefined, null, { timeout: 120000 });
  await page.evaluate(() => { document.getElementById('input').value = 'CC(=O)Oc1ccccc1O'; process(true); });
  await page.waitForFunction(() => !!document.querySelector('#svg svg'), null, { timeout: 60000 });
  await page.evaluate(() => embed3D());
  await page.waitForFunction(() => !!document.querySelector('#viewer3d canvas'), null, { timeout: 60000 });
  await page.evaluate(() => { switchMode('analog'); analogDetectSites(); });
  await page.waitForFunction(() => [...document.getElementById('analogSite').options].some(o => o.value), null, { timeout: 30000 });
}
const kickRun = page => page.evaluate(() => {
  document.getElementById('analogSite').selectedIndex = 1;
  document.getElementById('analogFragCat').value = 'alkyl';
  runAnalogExplorer();
});
const awaitDone = (page, t) => page.waitForFunction(() => /Done: \d+ scored/.test(document.getElementById('analogStatus').textContent), null, { timeout: t || 150000 });
const taskRows = page => page.evaluate(async () => {
  const st = window.__platformCompat.getStorage();
  return (await st.all('tasks')).map(r => r.value);
});

(async () => {
  const browser = await chromium.launch({ executablePath: EXE });

  // ============ 1. lifecycle + accounting (pool size 1 for determinism run A) ============
  {
    const ctx = await browser.newContext({ viewport: { width: 1440, height: 900 } });
    const page = await ctx.newPage();
    const errors = [];
    page.on('pageerror', e => errors.push('PAGEERROR: ' + e.message.slice(0, 200)));
    await setup(page, 1);
    kickRun(page);
    // tray shows the live task
    await page.waitForFunction(() => document.getElementById('taskTray').style.display !== 'none' && +document.getElementById('taskTrayCount').textContent >= 1, null, { timeout: 30000 })
      .then(() => check('task tray shows the running task', true))
      .catch(() => check('task tray shows the running task', false));
    await awaitDone(page);
    await page.waitForTimeout(1500);
    const rows = await taskRows(page);
    const rec = rows.find(r => r.state === 'complete');
    check('run completes; task record complete (spec persisted at submit)', !!rec, JSON.stringify(rows.map(r => r.state)));
    if (rec) {
      check('spec: effective stages snapshot + provenance.buildId + policy-relevant input',
        rec.spec.input.stages.length === 4 && rec.spec.input.stages[1].nConfs === 10 &&
        /^b/.test(rec.spec.provenance.buildId || '') && rec.spec.inputVersion === 1,
        JSON.stringify({ stages: rec.spec.input.stages.length, buildId: rec.spec.provenance.buildId }));
      const fact = await page.evaluate(async id => (await window.__platformCompat.getStorage().get('facts', 'analog:' + id)),
        rec.spec.taskId);
      const acc = fact && fact.accounts;
      check('Done fact: results + 4-way accounts + effective stages embedded',
        !!fact && Array.isArray(fact.results) && fact.results.length > 0 &&
        acc && acc.A.evaluated === rec.spec.input.candidates.length &&
        acc.B.evaluated + acc.B.failed + acc.B.notEvaluated === acc.A.evaluated &&
        acc.B.evaluated === fact.results.length &&
        (acc.D.evaluated + acc.D.failed + acc.D.filtered) <= 10 &&
        fact.stages.length === 4, JSON.stringify(acc));
      check('versions carry buildId', fact && fact.versions && /b/.test(fact.versions.buildId || ''), fact && JSON.stringify(fact.versions).slice(0, 90));
      // keep results for the determinism comparison
      await page.evaluate(res => { window.__runA = res; }, fact.results.map(r => [r.smiles, r.tanimoto, r.flexT]));
    }
    check('zero page errors (lifecycle)', errors.length === 0, errors.slice(0, 2).join(' | '));
    await ctx.close();
  }

  // ============ 2. determinism red line: same spec, pool size 1 vs 4 ============
  {
    // The two contexts above used DIFFERENT IndexedDBs (both default DB name,
    // but contexts are isolated). Re-do the comparison inside ONE context:
    const ctx = await browser.newContext({ viewport: { width: 1440, height: 900 } });
    const page = await ctx.newPage();
    const errors = [];
    page.on('pageerror', e => errors.push('PAGEERROR: ' + e.message.slice(0, 200)));
    // size-1 page
    await setup(page, 1);
    kickRun(page); await awaitDone(page); await page.waitForTimeout(1200);
    const runA = await page.evaluate(async () => {
      const st = window.__platformCompat.getStorage();
      const rec = (await st.all('tasks')).map(r => r.value).find(r => r.state === 'complete');
      const fact = rec ? await st.get('facts', 'analog:' + rec.spec.taskId) : null;
      return fact ? fact.results.map(r => [r.smiles, r.tanimoto, r.flexT]) : null;
    });
    // size-4 page (fresh DB via reload+delete)
    await setup(page, 4);
    kickRun(page); await awaitDone(page); await page.waitForTimeout(1200);
    const runB = await page.evaluate(async () => {
      const st = window.__platformCompat.getStorage();
      const rec = (await st.all('tasks')).map(r => r.value).find(r => r.state === 'complete');
      const fact = rec ? await st.get('facts', 'analog:' + rec.spec.taskId) : null;
      return fact ? fact.results.map(r => [r.smiles, r.tanimoto, r.flexT]) : null;
    });
    check('DETERMINISM RED LINE: pool 1 vs pool 4 — bitwise-identical results',
      !!runA && !!runB && JSON.stringify(runA) === JSON.stringify(runB),
      runA && runB ? ('first diff: ' + JSON.stringify(runA.map((r, i) => r.join('='))).slice(0, 60) + ' vs ' + JSON.stringify(runB.map((r, i) => r.join('='))).slice(0, 60)) : 'missing runs');
    check('zero page errors (determinism pair)', errors.length === 0, errors.slice(0, 2).join(' | '));
    await ctx.close();
  }

  // ============ 3. cancel mid-run: no partial facts ============
  {
    const ctx = await browser.newContext({ viewport: { width: 1440, height: 900 } });
    const page = await ctx.newPage();
    const errors = [];
    page.on('pageerror', e => errors.push('PAGEERROR: ' + e.message.slice(0, 200)));
    await setup(page, 0);
    kickRun(page);
    await page.waitForFunction(() => /B: \d+\/\d+/.test(document.getElementById('analogStatus').textContent), null, { timeout: 60000 });
    await page.waitForTimeout(1500);   // let a few score
    await page.evaluate(() => analogCancelRun());
    await page.waitForFunction(() => /取消|中断/.test(document.getElementById('analogStatus').textContent), null, { timeout: 30000 });
    await page.waitForTimeout(800);
    const rows = await taskRows(page);
    const cancelled = rows.find(r => r.state === 'cancelled');
    const facts = await page.evaluate(async () => (await window.__platformCompat.getStorage().all('facts')).filter(f => f.key.startsWith('analog:')));
    check('cancelled run: task record state=cancelled', !!cancelled, JSON.stringify(rows.map(r => r.state)));
    check('cancelled run: ZERO analog facts (nothing partial landed)', facts.length === 0, JSON.stringify(facts.map(f => f.key)));
    check('zero page errors (cancel)', errors.length === 0, errors.slice(0, 2).join(' | '));
    await ctx.close();
  }

  // ============ 4. reload reconciliation: abandoned + ghost-complete ============
  let reconWarn = [];
  {
    const ctx = await browser.newContext({ viewport: { width: 1440, height: 900 } });
    const page = await ctx.newPage();
    const errors = [];
    page.on('pageerror', e => errors.push('PAGEERROR: ' + e.message.slice(0, 200)));
    await setup(page, 0);
    // 5a: kick a run, reload mid-flight -> abandoned
    kickRun(page);
    await page.waitForFunction(() => document.getElementById('taskTrayCount').textContent === '1', null, { timeout: 60000 });
    await page.waitForTimeout(1200);
    await page.reload({ waitUntil: 'load' });
    await page.waitForFunction(() => window.webmm !== undefined, null, { timeout: 120000 });
    await page.waitForTimeout(1500);   // reconcile runs after platform store init
    await page.waitForFunction(() => window.__platformCompat && window.__platformCompat.getStorage(), null, { timeout: 20000 });
    const after = await page.evaluate(async () => {
      const st = window.__platformCompat.getStorage();
      return (await st.all('tasks')).map(r => ({ id: r.value.spec.taskId.slice(0, 6), state: r.value.state }));
    });
    check('reload mid-run: task reconciled to abandoned (never silently lost)', after.length === 1 && after[0].state === 'abandoned', JSON.stringify(after));
    // 5b: ghost-complete — inject a running task whose Done fact already exists
    await page.evaluate(async () => {
      const c = window.__platformCompat;
      const st = c.getStorage();
      const spec = window.Tasks.Contract.makeSpec({ taskId: 'aghost1', kind: 'analogScan', seed: 7,
        input: { parentSmiles: 'C', parentSdf: 'S', candidates: [{ idx: 0, smiles: 'CC', molblock: 'M' }], audit: { filtered: 0, notEvaluated: 0 }, inputRevision: { fragmentSet: 'builtin-v1' } },
        provenance: { buildId: 'btest' } });
      const ts = window.Tasks.Store.createStore(st);
      await ts.submit(spec);
      await ts.transition('aghost1', 'running');
      await st.put('facts', 'analog:aghost1', { layer: 'stochastic', results: [{ smiles: 'CC', tanimoto: 0.5 }] });
    });
    await page.reload({ waitUntil: 'load' });
    await page.waitForFunction(() => window.webmm !== undefined, null, { timeout: 120000 });
    await page.waitForTimeout(1500);
    await page.waitForFunction(() => window.__platformCompat && window.__platformCompat.getStorage(), null, { timeout: 20000 });
    const ghost = await page.evaluate(async () => {
      const st = window.__platformCompat.getStorage();
      return (await st.all('tasks')).map(r => ({ id: r.value.spec.taskId, state: r.value.state }));
    });
    check('ghost-complete: fact present reconciles to complete (no phantom abandonment)',
      ghost.some(r => r.id === 'aghost1' && r.state === 'complete'), JSON.stringify(ghost));
    check('zero page errors (reconcile)', errors.length === 0, errors.slice(0, 2).join(' | '));
    await ctx.close();
  }

  await browser.close();
  console.log('RESULT: ' + pass + ' passed, ' + fail + ' failed');
  process.exit(fail ? 1 : 0);
})().catch(e => { console.error('ERR', e); process.exit(1); });
