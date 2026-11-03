// m23_hardening.test.js — freeze-round hardening: conformers rerun path,
// reload abandonment of phase records, mid-run verify cancellation with the
// message-hold trick, malformed gather rejection through the UI path, and a
// bigger-molecule end-to-end smoke (performance envelope).
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
  return (await ts.all('tasks')).map(r => ({ id: r.value.spec.taskId, kind: r.value.spec.kind,
    phase: r.value.spec.input.phase || null, state: r.value.state, failure: r.value.failure,
    mdFrames: r.value.mdFrames, nConfs: r.value.nConfs, input: r.value.spec.input }));
});

(async () => {
  const browser = await chromium.launch({ executablePath: EXE });

  // ================= 1. conformers rerun path =================
  {
    const ctx = await browser.newContext({ viewport: { width: 1440, height: 900 } });
    const page = await ctx.newPage();
    const errors = [];
    page.on('pageerror', e => errors.push(e.message.slice(0, 150)));
    await page.goto(URL, { waitUntil: 'load' });
    await page.waitForFunction(() => window.webmm !== undefined, null, { timeout: 120000 });
    await page.evaluate(() => { document.getElementById('input').value = 'CCO ethanol'; process(true); });
    await page.waitForFunction(() => !!document.querySelector('#svg svg'), null, { timeout: 60000 });
    await page.evaluate(() => embed3D());
    await page.waitForFunction(() => !!document.querySelector('#viewer3d canvas'), null, { timeout: 90000 });
    await page.evaluate(() => { document.getElementById('confN').value = '8'; runConformers(); });
    await page.waitForFunction(() => /Conformers: \d+ optimized/.test(document.getElementById('status3d').textContent), null, { timeout: 180000 });
    const reran = await page.evaluate(async () => {
      const ts = await wbTaskStore.ensure();
      const rows = (await ts.all('tasks')).map(r => r.value);
      const rec = rows.filter(v => v.spec.kind === 'conformers' && v.state === 'complete').pop();
      if (!rec) return { ok: false, why: 'no record' };
      document.getElementById('confN').value = '50';   // perturb: rerun must override
      await window.rerunWbTask(rec.spec.taskId);
      return { ok: true, n: rec.spec.input.n };
    });
    await page.waitForFunction(() => /Conformers: \d+ optimized/.test(document.getElementById('status3d').textContent), null, { timeout: 180000 });
    const restored = await page.evaluate(() => ({
      n: document.getElementById('confN').value,
      kept: confEnsemble.length > 0 && confEnsemble.length <= 10,
    }));
    check('conformers rerun replays snapshot n (perturbed 50 → 8)',
      reran.ok && reran.n === 8 && restored.n === '8', JSON.stringify(reran) + restored.n);
    check('conformers rerun keeps the ensemble bounded by snapshot n', restored.kept,
          `kept=${restored.kept}`);
    check('conformers rerun: zero page errors', errors.length === 0, errors.join('|'));
    await ctx.close();
  }

  // ================= 2. reload abandons in-flight phase records =================
  {
    const ctx = await browser.newContext({ viewport: { width: 1440, height: 900 } });
    const page = await ctx.newPage();
    await page.goto(URL, { waitUntil: 'load' });
    await page.waitForFunction(() => window.webmm !== undefined, null, { timeout: 120000 });
    await page.evaluate(() => { document.getElementById('input').value = 'CCO ethanol'; process(true); });
    await page.waitForFunction(() => !!document.querySelector('#svg svg'), null, { timeout: 60000 });
    await page.evaluate(() => embed3D());
    await page.waitForFunction(() => !!document.querySelector('#viewer3d canvas'), null, { timeout: 90000 });
    await page.evaluate(() => { document.getElementById('imtdRounds').value = '2'; document.getElementById('imtdSteps').value = '800'; runImtd(); });
    await page.waitForFunction(() => /unique minima/.test(document.getElementById('status3d').textContent), null, { timeout: 180000 });
    // hold verify workers, start verify, reload mid-flight
    await page.evaluate(() => {
      const Orig = Worker;
      window.__held = [];
      window.Worker = class extends Orig {
        postMessage(m) {
          if (m && m.type === 'imtdRun') window.__held.push([super.postMessage.bind(this), m]);
          else super.postMessage(m);
        }
      };
      window.runMdVerify();
    });
    await page.waitForFunction(() => window.__held && window.__held.length >= 1, null, { timeout: 20000 });
    await page.reload({ waitUntil: 'load' });
    await page.waitForFunction(() => window.webmm !== undefined, null, { timeout: 120000 });
    await page.waitForTimeout(1500);   // boot reconcile settles
    const rows = await taskRows(page);
    const abandoned = rows && rows.filter(r => r.phase === 'verify' && r.state === 'abandoned');
    check('reload mid-verify → phase record abandoned (no facts to resurrect)',
      abandoned && abandoned.length >= 1, JSON.stringify(rows && rows.map(r => [r.phase, r.state])));
    await ctx.close();
  }

  // ================= 3. mid-run verify cancel → failed record with reason =================
  {
    const ctx = await browser.newContext({ viewport: { width: 1440, height: 900 } });
    const page = await ctx.newPage();
    const errors = [];
    page.on('pageerror', e => errors.push(e.message.slice(0, 150)));
    await page.goto(URL, { waitUntil: 'load' });
    await page.waitForFunction(() => window.webmm !== undefined, null, { timeout: 120000 });
    await page.evaluate(() => { document.getElementById('input').value = 'CCO ethanol'; process(true); });
    await page.waitForFunction(() => !!document.querySelector('#svg svg'), null, { timeout: 60000 });
    await page.evaluate(() => embed3D());
    await page.waitForFunction(() => !!document.querySelector('#viewer3d canvas'), null, { timeout: 90000 });
    await page.evaluate(() => { document.getElementById('imtdRounds').value = '2'; document.getElementById('imtdSteps').value = '800'; runImtd(); });
    await page.waitForFunction(() => /unique minima/.test(document.getElementById('status3d').textContent), null, { timeout: 180000 });
    await page.evaluate(() => {
      const Orig = Worker;
      window.__held2 = [];
      window.Worker = class extends Orig {
        postMessage(m) {
          if (m && m.type === 'imtdRun') window.__held2.push([super.postMessage.bind(this), m]);
          else super.postMessage(m);
        }
      };
      window.runMdVerify();
    });
    await page.waitForFunction(() => taskScheduler.stats().tasks.some(t => t.id.startsWith('v')), null, { timeout: 20000 });
    await page.waitForFunction(() => window.__held2 && window.__held2.length >= 1, null, { timeout: 20000 });
    // structure change mid-verify: the handle dies; the record must fail honestly
    await page.evaluate(() => { document.getElementById('input').value = 'CCC propane'; process(true); });
    await page.waitForTimeout(1500);
    const rows = await taskRows(page);
    const vrec = rows && rows.find(r => r.phase === 'verify' && (r.state === 'failed' || r.state === 'cancelled'));
    check('structure change mid-verify → phase record terminal with reason',
      !!vrec && (vrec.failure || '').length > 0, JSON.stringify(vrec && [vrec.state, vrec.failure]));
    const tokensOk = await page.evaluate(() => taskScheduler.stats().interactiveHeld === 0);
    check('verify tokens released after the failure', tokensOk);
    check('mid-run cancel: zero page errors', errors.length === 0, errors.join('|'));
    await ctx.close();
  }

  // ================= 4. bigger molecule end-to-end (perf envelope smoke) =================
  {
    const ctx = await browser.newContext({ viewport: { width: 1440, height: 900 } });
    const page = await ctx.newPage();
    const errors = [];
    page.on('pageerror', e => errors.push(e.message.slice(0, 150)));
    await page.goto(URL, { waitUntil: 'load' });
    await page.waitForFunction(() => window.webmm !== undefined, null, { timeout: 120000 });
    await page.evaluate(() => { document.getElementById('input').value = 'CC(C)Cc1ccc(cc1)C(C)C(=O)O ibuprofen'; process(true); });
    await page.waitForFunction(() => !!document.querySelector('#svg svg'), null, { timeout: 60000 });
    await page.evaluate(() => embed3D());
    await page.waitForFunction(() => !!document.querySelector('#viewer3d canvas'), null, { timeout: 120000 });
    const t0 = Date.now();
    await page.evaluate(() => { document.getElementById('imtdRounds').value = '2'; document.getElementById('imtdSteps').value = '600'; runImtd(); });
    await page.waitForFunction(() => /unique minima/.test(document.getElementById('status3d').textContent), null, { timeout: 300000 });
    const dt = ((Date.now() - t0) / 1000).toFixed(1);
    const ibu = await page.evaluate(() => ({
      n: confEnsemble.length,
      rounds: lastLandscape ? (lastLandscape.frames || 0) : 0,
      g: lastLandscape ? lastLandscape.basins.map(b => b.g) : [],
    }));
    check(`ibuprofen iMTD end-to-end (${dt}s)`, ibu.n >= 1 && ibu.rounds > 0,
          `n=${ibu.n} frames=${ibu.rounds}`);
    const dtVerify = await page.evaluate(async () => {
      const t = performance.now();
      await window.runMdVerify();
      return ((performance.now() - t) / 1000).toFixed(1);
    });
    await page.waitForFunction(() => /MD 占据验证(完成|失败)/.test(document.getElementById('status3d').textContent), null, { timeout: 300000 });
    const mdOk = await page.evaluate(() => {
      const md = lastLandscape && lastLandscape.md;
      return md && md.frames > 0;
    });
    check(`ibuprofen MD verify folds (${dtVerify}s)`, mdOk);
    check('ibuprofen E2E: zero page errors', errors.length === 0, errors.join('|'));
    await ctx.close();
  }

  console.log(`\nRESULT: ${pass} passed, ${fail} failed`);
  await browser.close();
  process.exit(fail ? 1 : 0);
})().catch(e => { console.error('FATAL:', e); process.exit(1); });
