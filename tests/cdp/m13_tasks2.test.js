// m13_tasks2.test.js — Top-5 stage 3: optimize + shapeSearch kinds.
// optimize: typed task with the application guard (structure hash; stale on
// mid-run change, result NEVER silently applied). shapeSearch: orchestration
// shell (spec + tray + W global tokens for the farm; results still land via
// CreateRound). rerun reads the SNAPSHOT spec from the side table.
const { chromium } = require('/opt/homebrew/lib/node_modules/playwright');
const EXE = '/Users/rliu/Library/Caches/ms-playwright/chromium_headless_shell-1243/chrome-headless-shell-mac-arm64/chrome-headless-shell';
const URL = 'http://localhost:8901/app/platform.html';

let pass = 0, fail = 0;
function check(name, cond, detail) {
  if (cond) { pass++; console.log('  ✓ ' + name); }
  else { fail++; console.log('  ✗ ' + name + (detail ? ' — ' + detail : '')); }
}

(async () => {
  const browser = await chromium.launch({ executablePath: EXE });
  const ctx = await browser.newContext({ viewport: { width: 1440, height: 900 } });
  const page = await ctx.newPage();
  const errors = [];
  const dialogs = [];
  page.on('pageerror', e => errors.push('PAGEERROR: ' + e.message.slice(0, 200)));
  page.on('dialog', d => { dialogs.push(d.message().slice(0, 120)); d.accept().catch(() => {}); });

  await page.goto(URL, { waitUntil: 'load' });
  await page.waitForFunction(() => window.webmm !== undefined, null, { timeout: 60000 });
  await page.evaluate(() => { indexedDB.deleteDatabase('webmm-platform'); indexedDB.deleteDatabase('webmm-projects'); });
  await page.waitForTimeout(300);
  await page.reload({ waitUntil: 'load' });
  await page.waitForFunction(() => window.webmm !== undefined, null, { timeout: 60000 });

  const taskRows = () => page.evaluate(async () => {
    await new Promise(r => { (function w() { window.__platformCompat && window.__platformCompat.getStorage() ? r() : setTimeout(w, 100); })(); });
    const st = window.__platformCompat.getStorage();
    return (await st.all('tasks')).map(r => ({ id: r.value.spec.taskId, kind: r.value.spec.kind, state: r.value.state,
      E: r.value.E, converged: r.value.converged, durationMs: r.value.durationMs, input: r.value.spec.input }));
  });

  // ================= 1. optimize: normal lifecycle =================
  await page.evaluate(() => { document.getElementById('input').value = 'CC(=O)Oc1ccccc1O paracetamol'; process(true); });
  await page.waitForFunction(() => !!document.querySelector('#svg svg'), null, { timeout: 60000 });
  await page.evaluate(() => embed3D());
  await page.waitForFunction(() => !!document.querySelector('#viewer3d canvas'), null, { timeout: 60000 });
  await page.evaluate(() => optimize3D());
  await page.waitForFunction(() => document.getElementById('status3d').textContent.length > 0, null, { timeout: 90000 });
  await page.waitForFunction(() => !/Optimizing/.test(document.getElementById('status3d').textContent), null, { timeout: 90000 });
  await page.waitForTimeout(2500);   // let the side-table transition land
  let rows = await taskRows();
  let opt = rows.find(r => r.kind === 'optimize' && r.state === 'complete');
  check('optimize task completes with E/converged/duration recorded', !!opt && opt.converged === true && isFinite(opt.E) && opt.durationMs > 0, JSON.stringify(rows.map(r => [r.kind, r.state])));

  // ================= 2. optimize: stale (structure changed mid-run) =================
  // Fully deterministic: intercept the opt worker and HOLD the 'optimize'
  // message until after the molecule swap — the reply then arrives with
  // molEpoch bumped and the structure hash changed, so the guard MUST fire.
  await page.evaluate(() => { document.getElementById('input').value = 'CC(=O)Oc1ccccc1O paracetamol'; process(true); });
  await page.waitForFunction(() => !!document.querySelector('#svg svg'), null, { timeout: 60000 });
  await page.evaluate(() => embed3D());
  await page.waitForFunction(() => !!document.querySelector('#viewer3d canvas'), null, { timeout: 60000 });
  await page.evaluate(() => {
    const orig = getOptWorker;
    window.__optHeld = [];
    window.__optHold = true;
    getOptWorker = function () {
      const w = orig();
      const realPost = w.postMessage.bind(w);
      w.postMessage = m => { if (window.__optHold && m && m.type === 'optimize') window.__optHeld.push([realPost, m]); else realPost(m); };
      return w;
    };
    optimize3D();
  });
  await page.waitForFunction(() => taskScheduler.stats().tasks.some(t => t.id.startsWith('o')), null, { timeout: 20000 });
  await page.evaluate(() => { document.getElementById('input').value = 'CCO swapped'; process(true); });   // molEpoch++ + structure swapped
  await page.evaluate(() => { window.__optHold = false; for (const [post, m] of window.__optHeld) post(m); window.__optHeld = []; });   // release the job NOW
  await page.waitForTimeout(6000);   // reply arrives; guard must reject application
  rows = await taskRows();
  const stale = rows.find(r => r.kind === 'optimize' && r.state === 'stale');
  check('mid-run structure change → task marked stale (no silent application)', !!stale, JSON.stringify(rows.map(r => [r.kind, r.state])));
  const staleText = await page.evaluate(() => document.body.textContent.slice(0, 20000));
  check('stale surfaced to the user with rerun guidance (stale + 未套用)',
    /stale/.test(staleText) && (/未套用|切回原结构/.test(staleText)), 'no visible marker');
  // restore the real worker getter for later sections
  await page.evaluate(() => { delete window.__optHold; delete window.__optHeld; });

  // ================= 3. rerun: hash guard (structure still swapped) =================
  dialogs.length = 0;
  await page.evaluate(id => rerunTask(id), stale ? stale.id : 'x');
  await page.waitForTimeout(800);
  check('rerun on optimize refuses to run against a changed structure (hash guard alert)',
    dialogs.some(d => /哈希|变化/.test(d)), JSON.stringify(dialogs));

  // ================= 4. shapeSearch: orchestration-shell task =================
  await page.evaluate(async () => {
    switchMode('search');
    await new Promise(r => loadSearchLibraryFrom([
      { smiles: 'CC(=O)Oc1ccccc1O', name: 'paracetamol' },
      { smiles: 'CC(=O)Oc1ccccc1C(=O)O', name: 'aspirin' },
      { smiles: 'Oc1ccccc1O', name: 'resorcinol' },
      { smiles: 'CCO', name: 'ethanol' },
    ]).persisted.then(r));
  });
  await page.waitForTimeout(500);
  // single-conformer path (nConfs=1)
  await page.evaluate(() => {
    document.getElementById('searchQuery').value = 'CC(=O)Oc1ccccc1O';
    document.getElementById('searchMode').value = 'shape';
    document.getElementById('shapeConfs').value = '1';
    runSearch();
  });
  await page.waitForFunction(() => document.querySelectorAll('#searchRows tr').length >= 1, null, { timeout: 120000 });
  await page.waitForTimeout(1500);
  rows = await taskRows();
  const ss1 = rows.find(r => r.kind === 'shapeSearch' && r.state === 'complete');
  check('shapeSearch (single-conf) completes as a task with revision + snapshot input',
    !!ss1 && ss1.input.threshold > 0 && ss1.input.nConfs === 1 && /^h/.test(ss1.input.qsdfHash) && /^h/.test(ss1.input.libraryFingerprint),
    JSON.stringify(rows.filter(r => r.kind === 'shapeSearch').map(r => [r.state, r.input.nConfs])));

  // ensemble path (nConfs=4) — the farm must not wedge the global tokens
  await page.evaluate(() => { document.getElementById('shapeConfs').value = '4'; runSearch(); });
  await page.waitForFunction(() => document.querySelectorAll('#searchRows tr').length >= 1, null, { timeout: 180000 });
  await page.waitForTimeout(2500);
  rows = await taskRows();
  const ss2 = rows.filter(r => r.kind === 'shapeSearch' && r.state === 'complete');
  const schedOk = await page.evaluate(() => { const s = taskScheduler.stats(); return s.interactiveHeld === 0 && s.backgroundHeld === 0; });
  check('shapeSearch (ensemble farm) completes and RELEASES all global tokens', ss2.length >= 2 && schedOk,
    JSON.stringify({ n: ss2.length, sched: schedOk }));

  // ================= 5. rerun: shapeSearch from the snapshot =================
  dialogs.length = 0;
  const before = rows.filter(r => r.kind === 'shapeSearch' && r.state === 'complete').length;
  await page.evaluate(id => rerunTask(id), ss2[ss2.length - 1].id);
  await page.waitForFunction(() => document.querySelectorAll('#searchRows tr').length >= 1, null, { timeout: 120000 });
  await page.waitForTimeout(1500);
  rows = await taskRows();
  const after = rows.filter(r => r.kind === 'shapeSearch' && r.state === 'complete').length;
  check('rerun re-executes from the snapshot spec (new complete record)', after === before + 1, JSON.stringify({ before, after }));

  check('zero page errors', errors.length === 0, errors.slice(0, 3).join(' | '));
  await browser.close();
  console.log('RESULT: ' + pass + ' passed, ' + fail + ' failed');
  process.exit(fail ? 1 : 0);
})().catch(e => { console.error('ERR', e); process.exit(1); });
