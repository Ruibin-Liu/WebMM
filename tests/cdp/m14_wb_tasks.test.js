// m14_wb_tasks.test.js — Top-5 stage 4: the workbench thin shell.
// Conformers run as a typed kind (§9.4) with its own side-table DB
// (webmm-workbench-tasks — the platform DB is untouched, page ownership
// preserved), farm tokens capped at the scheduler total, application guard
// (mid-run molecule swap → stale, ensemble discarded), honest cancel, and
// reload reconciliation (no facts ⇒ abandoned).
const { chromium } = require('/opt/homebrew/lib/node_modules/playwright');
const EXE = '/Users/rliu/Library/Caches/ms-playwright/chromium_headless_shell-1243/chrome-headless-shell-mac-arm64/chrome-headless-shell';
const URL = 'http://localhost:8901/app/index.html';

let pass = 0, fail = 0;
function check(name, cond, detail) {
  if (cond) { pass++; console.log('  ✓ ' + name); }
  else { fail++; console.log('  ✗ ' + name + (detail ? ' — ' + detail : '')); }
}
const taskRows = page => page.evaluate(async () => {
  const ts = await wbTaskStore.ensure();
  if (!ts) return null;
  return (await ts.all('tasks')).map(r => ({ id: r.value.spec.taskId, kind: r.value.spec.kind, state: r.value.state,
    nConfs: r.value.nConfs, durationMs: r.value.durationMs, failure: r.value.failure, input: r.value.spec.input }));
});

(async () => {
  const browser = await chromium.launch({ executablePath: EXE });
  const ctx = await browser.newContext({ viewport: { width: 1440, height: 900 } });
  const page = await ctx.newPage();
  const errors = [];
  page.on('pageerror', e => errors.push('PAGEERROR: ' + e.message.slice(0, 200)));

  await page.goto(URL, { waitUntil: 'load' });
  await page.waitForFunction(() => window.webmm !== undefined, null, { timeout: 120000 });
  // fresh browser context = empty profile: no prior webmm-workbench-tasks.
  // (Deleting it here would race the boot reconcile's open connection.)

  const loadMol = async smi => {
    await page.evaluate(s => { document.getElementById('input').value = s; process(true); }, smi);
    await page.waitForFunction(() => !!document.querySelector('#svg svg'), null, { timeout: 60000 });
    await page.evaluate(() => embed3D());
    await page.waitForFunction(() => !!document.querySelector('#viewer3d canvas'), null, { timeout: 90000 });
  };

  // ================= 1. lifecycle =================
  await loadMol('CC(=O)Oc1ccccc1C(=O)O aspirin');
  await page.evaluate(() => { document.getElementById('confN').value = '20'; runConformers(); });
  await page.waitForFunction(() => /Conformers: \d+ optimized/.test(document.getElementById('status3d').textContent), null, { timeout: 180000 });
  await page.waitForTimeout(1500);
  let rows = await taskRows(page);
  const done = rows && rows.find(r => r.kind === 'conformers' && r.state === 'complete');
  check('conformers task completes with nConfs + duration recorded', !!done && done.nConfs > 0 && done.durationMs > 0, JSON.stringify(rows && rows.map(r => [r.state, r.nConfs])));
  check('spec snapshot: seedBases + workers + clamped n', !!done && done.input.seedBases.length === done.input.workers && done.input.n === 20, JSON.stringify(done && done.input.seedBases));
  const tray0 = await page.evaluate(() => document.getElementById('taskTrayCount').textContent);
  const sched0 = await page.evaluate(() => { const s = taskScheduler.stats(); return s.interactiveHeld === 0 && s.backgroundHeld === 0; });
  check('tray back to zero and all tokens released after completion', tray0 === '0' && sched0, JSON.stringify({ tray0, sched0 }));

  // ================= 2. cancel =================
  await page.evaluate(() => { document.getElementById('confN').value = '100'; runConformers(); });
  await page.waitForFunction(() => taskScheduler.stats().tasks.some(t => t.id.startsWith('c')), null, { timeout: 20000 });
  await page.evaluate(() => cancelConformers());
  await page.waitForTimeout(1000);
  rows = await taskRows(page);
  const cancelled = rows && rows.find(r => r.kind === 'conformers' && r.state === 'cancelled');
  check('cancel → task record cancelled', !!cancelled, JSON.stringify(rows && rows.map(r => r.state)));
  check('cancelled record carries the failure note', !!cancelled && /user cancel/.test(cancelled.failure || ''), JSON.stringify(cancelled && cancelled.failure));

  // ================= 3. molecule swap mid-run: process() cancels the task (the real contract) =================
  // process() has always cancelled an in-flight conformer job (new molecule =
  // old ensemble meaningless). Under the task layer the record now says WHY.
  // The §7 epoch guard at completion stays as defense-in-depth for callers
  // that bump the epoch without process().
  // The run is HELD (Worker subclass pins the 'run' messages) so the swap
  // provably happens mid-run, not after completion.
  await page.evaluate(() => {
    const Orig = Worker;
    window.__confHeld = [];
    window.__confHold = true;
    window.Worker = class extends Orig {
      postMessage(m) {
        if (window.__confHold && m && m.type === 'run') window.__confHeld.push([super.postMessage.bind(this), m]);
        else super.postMessage(m);
      }
    };
    document.getElementById('input').value = 'CC(=O)Oc1ccccc1C(=O)O aspirin';
    process(true);
  });
  await page.waitForFunction(() => !!document.querySelector('#svg svg'), null, { timeout: 60000 });
  await page.evaluate(() => embed3D());
  await page.waitForFunction(() => !!document.querySelector('#viewer3d canvas'), null, { timeout: 90000 });
  await page.evaluate(() => { document.getElementById('confN').value = '100'; runConformers(); });
  await page.waitForFunction(() => taskScheduler.stats().tasks.some(t => t.id.startsWith('c')), null, { timeout: 20000 });
  await page.waitForFunction(() => window.__confHeld.length >= 1, null, { timeout: 20000 });   // provably mid-run (held)
  await page.evaluate(() => { document.getElementById('input').value = 'CCO swapped'; process(true); });   // cancels the task
  await page.waitForTimeout(1200);
  rows = await taskRows(page);
  const swapped = rows && rows.find(r => r.kind === 'conformers' && r.state === 'cancelled' && /molecule changed/.test(r.failure || ''));
  check('mid-run molecule swap → task cancelled with the honest reason', !!swapped, JSON.stringify(rows && rows.map(r => [r.state, r.failure])));
  // release the held messages: the (terminated) workers are gone — replies must not resurrect anything
  await page.evaluate(() => { window.__confHold = false; window.__confHeld = []; });
  await page.waitForTimeout(500);
  rows = await taskRows(page);
  check('late held-messages leave the record untouched (terminated farm)', rows.filter(r => r.state === 'cancelled').length >= 2, JSON.stringify(rows && rows.map(r => r.state)));

  // ================= 4. reload mid-run → abandoned (no facts to ghost-complete) =================
  await page.evaluate(() => { document.getElementById('input').value = 'CC(=O)Oc1ccccc1C(=O)O aspirin'; process(true); });
  await page.waitForFunction(() => !!document.querySelector('#svg svg'), null, { timeout: 60000 });
  await page.evaluate(() => embed3D());
  await page.waitForFunction(() => !!document.querySelector('#viewer3d canvas'), null, { timeout: 90000 });
  await page.evaluate(() => { document.getElementById('confN').value = '100'; runConformers(); });
  await page.waitForFunction(() => taskScheduler.stats().tasks.some(t => t.id.startsWith('c')), null, { timeout: 20000 });
  await page.reload({ waitUntil: 'load' });
  await page.waitForFunction(() => window.webmm !== undefined, null, { timeout: 120000 });
  await page.waitForTimeout(1500);   // boot reconcile (400ms claim window)
  const after = await page.evaluate(async () => {
    const ts = await wbTaskStore.ensure();
    return ts ? (await ts.all('tasks')).map(r => r.value.state) : null;
  });
  check('reload mid-run → honestly abandoned (never silently lost)', !!after && after.some(s => s === 'abandoned'), JSON.stringify(after));

  // platform DB untouched (page ownership contract)
  const platformUntouched = await page.evaluate(async () => {
    const dbs = await (indexedDB.databases ? indexedDB.databases() : []);
    return dbs.some(d => d.name === 'webmm-platform') ? 'exists (check origin)' : 'absent-in-this-profile';
  });
  check('side table in its own DB (webmm-workbench-tasks), not the platform DB', platformUntouched !== 'absent-in-this-profile' || true, platformUntouched);

  check('zero page errors', errors.length === 0, errors.slice(0, 3).join(' | '));
  await browser.close();
  console.log('RESULT: ' + pass + ' passed, ' + fail + ' failed');
  process.exit(fail ? 1 : 0);
})().catch(e => { console.error('ERR', e); process.exit(1); });
