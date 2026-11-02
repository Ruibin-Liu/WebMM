// m16_imtd_ui.test.js — iMTD-GC workbench integration (contract §9.5).
// The iMTD button runs a typed task: parallel scatter-gather through
// window.IMTD, ensemble feeds the shared conformer UI, honest cancellation
// on molecule swap (Worker-subclass hold makes the swap provably mid-run),
// GFN-FF engine smoke, spec snapshot + token discipline.
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
    nConfs: r.value.nConfs, roundsExecuted: r.value.roundsExecuted, converged: r.value.converged,
    durationMs: r.value.durationMs, failure: r.value.failure, input: r.value.spec.input }));
});

(async () => {
  const browser = await chromium.launch({ executablePath: EXE });
  const ctx = await browser.newContext({ viewport: { width: 1440, height: 900 } });
  const page = await ctx.newPage();
  const errors = [];
  page.on('pageerror', e => errors.push('PAGEERROR: ' + e.message.slice(0, 200)));

  await page.goto(URL, { waitUntil: 'load' });
  await page.waitForFunction(() => window.webmm !== undefined, null, { timeout: 120000 });

  const loadMol = async smi => {
    await page.evaluate(s => { document.getElementById('input').value = s; process(true); }, smi);
    await page.waitForFunction(() => !!document.querySelector('#svg svg'), null, { timeout: 60000 });
    await page.evaluate(() => embed3D());
    await page.waitForFunction(() => !!document.querySelector('#viewer3d canvas'), null, { timeout: 90000 });
  };

  // ================= 1. lifecycle: typed task + shared conformer UI =================
  await loadMol('CCO ethanol');
  check('iMTD button enabled after embed', await page.evaluate(() => !document.getElementById('btnImtd').disabled));
  await page.evaluate(() => {
    document.getElementById('imtdRounds').value = '2';
    document.getElementById('imtdSteps').value = '600';
    runImtd();
  });
  await page.waitForFunction(() => /unique minima/.test(document.getElementById('status3d').textContent), null, { timeout: 180000 });
  await page.waitForTimeout(1200);   // side-table transition settles
  const st = await page.evaluate(() => ({
    status: document.getElementById('status3d').textContent,
    n: confEnsemble.length,
    Es: confEnsemble.map(k => k.E),
    panel: document.getElementById('confPanel').style.display,
    button: document.getElementById('btnImtd').textContent,
    hasSdf: confEnsemble.every(k => k.sdf && k.sdf.includes('V2000')),
  }));
  check('ensemble found (ethanol anti/gauche ≥2)', st.n >= 2, `${st.n} conformers`);
  check('energies sorted ascending', st.Es.every((e, i) => i === 0 || st.Es[i - 1] <= e), JSON.stringify(st.Es.map(e => e.toFixed(3))));
  check('conformer panel visible with per-conformer SDFs', st.panel === 'block' && st.hasSdf);
  check('button back to idle label', st.button === 'iMTD', st.button);
  check('status reports converged/budget + rounds + workers', /round/.test(st.status) && /worker/.test(st.status), st.status.slice(0, 120));

  let rows = await taskRows(page);
  const done = rows && rows.find(r => r.kind === 'imtd' && r.state === 'complete');
  check('imtd task completes with nConfs/rounds/duration', !!done && done.nConfs >= 2 && done.roundsExecuted >= 1 && done.durationMs > 0,
        JSON.stringify(rows && rows.map(r => [r.state, r.nConfs, r.roundsExecuted])));
  check('spec snapshot: rounds clamped + workers ≤ 4 + engine', !!done && done.input.rounds === 2 && done.input.workers >= 1 && done.input.workers <= 4 && done.input.engine === 'mmff',
        JSON.stringify(done && done.input));
  const tray0 = await page.evaluate(() => document.getElementById('taskTrayCount').textContent);
  const sched0 = await page.evaluate(() => { const s = taskScheduler.stats(); return s.interactiveHeld === 0 && s.backgroundHeld === 0; });
  check('tray zero + all tokens released', tray0 === '0' && sched0, JSON.stringify({ tray0, sched0 }));

  // ================= 2. GFN-FF engine smoke =================
  await page.selectOption('#engineSel', 'GFNFF');
  await page.evaluate(() => runImtd());
  await page.waitForFunction(() => /unique minima/.test(document.getElementById('status3d').textContent) && /GFNFF/.test(document.getElementById('status3d').textContent),
    null, { timeout: 180000 });
  const g = await page.evaluate(() => ({ n: confEnsemble.length, E0: confEnsemble[0].E }));
  check('GFN-FF engine path finds minima (finite energies)', g.n >= 1 && isFinite(g.E0), `n=${g.n} E0=${g.E0}`);
  await page.selectOption('#engineSel', 'MMFF94s');

  // ================= 3. mid-run molecule swap cancels honestly =================
  // Hold the workers' imtdRun messages so the swap provably happens mid-round.
  await page.evaluate(() => {
    const Orig = Worker;
    window.__imtdHeld = [];
    window.__imtdHold = true;
    window.Worker = class extends Orig {
      postMessage(m) {
        if (window.__imtdHold && m && m.type === 'imtdRun') window.__imtdHeld.push([super.postMessage.bind(this), m]);
        else super.postMessage(m);
      }
    };
    document.getElementById('imtdRounds').value = '4';
    document.getElementById('imtdSteps').value = '2000';
    runImtd();
  });
  await page.waitForFunction(() => taskScheduler.stats().tasks.some(t => t.id.startsWith('i')), null, { timeout: 20000 });
  await page.waitForFunction(() => window.__imtdHeld.length >= 1, null, { timeout: 20000 });   // provably mid-run (held)
  await page.evaluate(() => { document.getElementById('input').value = 'CCO swapped-molecule'; process(true); });   // reset3DState → cancelImtd('molecule changed')
  await page.waitForTimeout(1500);
  rows = await taskRows(page);
  const swapped = rows && rows.find(r => r.kind === 'imtd' && r.state === 'cancelled' && /molecule changed/.test(r.failure || ''));
  check('mid-run molecule swap → cancelled with the honest reason', !!swapped, JSON.stringify(rows && rows.map(r => [r.state, r.failure])));
  const btnAfter = await page.evaluate(() => document.getElementById('btnImtd').textContent);
  check('button reset after cancellation', btnAfter === 'iMTD', btnAfter);
  const sched1 = await page.evaluate(() => { const s = taskScheduler.stats(); return s.interactiveHeld === 0; });
  check('tokens released after cancellation', sched1);
  // release the held messages: workers were terminated — nothing may resurrect
  await page.evaluate(() => { window.__imtdHold = false; window.__imtdHeld = []; });
  await page.waitForTimeout(500);
  rows = await taskRows(page);
  const stillCancelled = rows.filter(r => r.kind === 'imtd' && r.state === 'cancelled').length;
  check('late held-messages leave the record untouched', stillCancelled >= 1, JSON.stringify(rows && rows.map(r => r.state)));

  check('zero page errors', errors.length === 0, errors.join(' | '));

  console.log(`\nRESULT: ${pass} passed, ${fail} failed`);
  await browser.close();
  process.exit(fail ? 1 : 0);
})().catch(e => { console.error('FATAL:', e); process.exit(1); });
