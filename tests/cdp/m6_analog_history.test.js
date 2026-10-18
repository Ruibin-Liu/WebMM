// m6_analog_history.test.js — M4 P2: explorer runs persist into the project.
// A completed funnel auto-commits (command + facts, one transaction); the
// history list renders in analog mode; reload restores it; restore brings
// the results table back (rebuilding molblocks for SDF export); hide =
// undo semantics (facts never deleted). Docs: project-management.md §4.
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
  const page = await (await browser.newContext()).newPage();
  const errors = [];
  page.on('pageerror', e => errors.push('PAGEERROR: ' + e.message.slice(0, 200)));
  page.on('console', m => { if (m.type() === 'error') errors.push(m.text().slice(0, 150)); });

  const boot = async () => {
    await page.goto(URL, { waitUntil: 'load' });
    await page.waitForFunction(() => window.webmm !== undefined, null, { timeout: 60000 });
    await page.evaluate(() => { try { switchMode('search'); } catch (e) {} });
    await page.waitForFunction(() => window.__platformCompat, null, { timeout: 60000 });
  };
  await boot();

  // ---- run an explorer pass on paracetamol (alkyl, small + fast) ----
  await page.evaluate(() => { document.getElementById('input').value = 'CC(=O)Oc1ccccc1O'; process(true); });
  await page.waitForFunction(() => !!document.querySelector('#svg svg'), null, { timeout: 60000 });
  await page.evaluate(() => embed3D());
  await page.waitForFunction(() => !!document.querySelector('#viewer3d canvas'), null, { timeout: 60000 });
  await page.evaluate(() => { switchMode('analog'); analogDetectSites(); });
  await page.waitForFunction(() => [...document.getElementById('analogSite').options].some(o => o.value), null, { timeout: 30000 });
  await page.evaluate(() => {
    document.getElementById('analogSite').selectedIndex = 1;
    document.getElementById('analogFragCat').value = 'alkyl';
    runAnalogExplorer();
  });
  await page.waitForFunction(() => /Done: \d+ scored/.test(document.getElementById('analogStatus').textContent), null, { timeout: 180000 });

  // ---- auto-commit: status carries the save marker + history renders ----
  const committed = await page.evaluate(() => ({
    status: document.getElementById('analogStatus').textContent,
    histWrap: document.getElementById('analogHistoryWrap').style.display !== 'none',
    histItems: document.getElementById('analogHistory').children.length,
    runs: Object.keys((window.__platformCompat.projectApi().state.analogRuns) || {}).length,
  }));
  check('completed run auto-commits (status marker + history entry)', /已保存/.test(committed.status) && committed.histWrap && committed.histItems === 1 && committed.runs === 1, JSON.stringify(committed).slice(0, 160));

  // ---- facts landed with lineage + minimal fields ----
  const fact = await page.evaluate(async () => {
    const c = window.__platformCompat;
    const id = Object.keys(c.projectApi().state.analogRuns)[0];
    const f = await c.getFact('analog:' + id);
    return { id, layer: f && f.layer, versions: !!(f && f.versions && f.versions.rdkit), n: f ? f.results.length : 0,
      fields: f && f.results[0] ? Object.keys(f.results[0]).sort().join(',') : '' };
  });
  check('facts: stochastic layer + version lineage + minimal fields',
    fact.layer === 'stochastic' && fact.versions && fact.n > 0 && /^flexT,idx,mw,rotb,sa,smiles,tanimoto$/.test(fact.fields), JSON.stringify(fact));

  // ---- reload: history persists, restore rebuilds the table ----
  await boot();
  await page.evaluate(() => switchMode('analog'));
  await page.waitForFunction(() => document.getElementById('analogHistory').children.length === 1, null, { timeout: 20000 });
  check('history survives reload', true);
  await page.evaluate(() => analogRestoreRun(Object.keys(window.__platformCompat.projectApi().state.analogRuns)[0]));
  await page.waitForFunction(() => document.querySelectorAll('#analogRows tr').length > 0, null, { timeout: 20000 });
  const restored = await page.evaluate(() => ({
    rows: document.querySelectorAll('#analogRows tr').length,
    status: document.getElementById('analogStatus').textContent.slice(0, 50),
    status2: document.getElementById('analogStatus').textContent,
  }));
  check('restore rebuilds the results table', restored.rows > 0 && /已恢复|旧版本/.test(restored.status2), JSON.stringify(restored));
  // restored rows support export (molblock rebuilt from SMILES)
  const dl = page.waitForEvent('download', { timeout: 15000 });
  await page.evaluate(() => exportAnalogSDF());
  await (await dl).saveAs('/tmp/analog-restore.sdf');
  const sdf = require('fs').readFileSync('/tmp/analog-restore.sdf', 'utf8');
  check('restored rows export as SDF (molblocks rebuilt)', (sdf.match(/V2000/g) || []).length === restored.rows, 'v2k=' + (sdf.match(/V2000/g) || []).length);

  // ---- hide = undo semantics: list empties, facts survive ----
  await page.evaluate(() => analogHideRun(Object.keys(window.__platformCompat.projectApi().state.analogRuns)[0]));
  await page.waitForFunction(() => document.getElementById('analogHistory').children.length === 0, null, { timeout: 20000 });
  const afterHide = await page.evaluate(async () => {
    const c = window.__platformCompat;
    const id = Object.keys(c.projectApi().state.analogRuns)[0];
    return { hidden: !!c.projectApi().state.hiddenRuns[id], fact: !!(await c.getFact('analog:' + id)) };
  });
  check('hide removes from list; run + facts survive (undoable)', afterHide.hidden && afterHide.fact, JSON.stringify(afterHide));

  check('no page errors', errors.length === 0, errors.slice(0, 2).join(' | '));
  await browser.close();
  console.log('RESULT: ' + pass + ' passed, ' + fail + ' failed');
  process.exit(fail ? 1 : 0);
})().catch(e => { console.error(e); process.exit(1); });
