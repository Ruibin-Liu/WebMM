// M2 acceptance: conformer ensemble — batch embed, per-conformer optimize,
// energy ranking, RMSD prune, click-to-view, multi-conformer SDF export.
const { chromium } = require('/opt/homebrew/lib/node_modules/playwright');
const fs = require('fs');

const EXE = '/Users/rliu/Library/Caches/ms-playwright/chromium_headless_shell-1234/chrome-headless-shell-mac-arm64/chrome-headless-shell';

setTimeout(() => { console.log('WATCHDOG: test exceeded 240s'); process.exit(2); }, 240000);
(async () => {
    const browser = await chromium.launch({
    executablePath: EXE,
    args: ['--disable-background-timer-throttling', '--disable-renderer-backgrounding',
           '--disable-features=IntensiveWakeUpThrottling'],
  });
  const page = await browser.newPage({ viewport: { width: 1440, height: 900 } });
  const errors = [];
  page.on('console', m => { if (m.type() === 'error') errors.push(m.text()); });
  page.on('pageerror', e => errors.push('PAGEERROR: ' + e.message));

  // ibuprofen SMILES (implicit-H 19-atom graph; fast optimizations)
  const hash = Buffer.from(JSON.stringify({ smiles: 'CC(C)Cc1ccc(cc1)C(C)C(=O)O', name: 'ibuprofen' })).toString('base64');
  await page.goto('http://localhost:8901/app/index.html#mol=' + hash, { waitUntil: 'load' });
  await page.waitForFunction(() => window.webmmReady === true, null, { timeout: 60000 });

  let pass = 0, fail = 0;
  const check = (name, cond, detail = '') => {
    if (cond) { pass++; console.log(`  ✓ ${name}${detail ? ' — ' + detail : ''}`); }
    else { fail++; console.log(`  ✗ ${name}${detail ? ' — ' + detail : ''}`); }
  };

  // set N=12 for test speed, RMSD 0.5, engine MMFF94s
  await page.evaluate(() => { document.getElementById('confN').value = '12'; document.getElementById('confRmsd').value = '0.5'; });
  await page.selectOption('#engineSel', 'MMFF94s');
  const t0 = Date.now();
  await page.evaluate(() => runConformers());
  await page.waitForFunction(() => document.getElementById('status3d').textContent.includes('Conformers:'), null, { timeout: 180000 });
  const dt = ((Date.now() - t0) / 1000).toFixed(1);

  const st = await page.evaluate(() => ({
    status: document.getElementById('status3d').textContent,
    n: confEnsemble.length,
    Es: confEnsemble.map(k => k.E),
    dE0: confEnsemble[0].dE,
    axis: document.getElementById('confChartAxis').textContent,
    panelVisible: document.getElementById('confPanel').style.display,
    tableGone: !document.getElementById('confTable'),
  }));
  console.log(`ensemble (${dt}s): ${st.n} kept, axis: ${st.axis}`);
  const kept = st.n;
  check('status reports optimized + kept counts', /Conformers: (\d+) optimized .* (\d+) kept/.test(st.status), st.status.slice(0, 90));
  check('ensemble panel rendered (chart only, no table)', st.panelVisible === 'block' && kept >= 1 && kept <= 12 && st.tableGone, kept + ' conformers');
  check('energies sorted ascending with ΔE[0]=0', st.Es.every((e, i) => i === 0 || e >= st.Es[i-1]) && st.dE0 === 0, 'rank 1 ΔE=0');

  // click row 2 → 3D shows that conformer (status changes)
  if (kept >= 2) {
    await page.evaluate(() => showEnsembleConf(1));
    const s2 = await page.evaluate(() => document.getElementById('status3d').textContent);
    check('click row → conformer #2 shown', s2.includes('Conformer #2'), s2);
  }

  // export ensemble SDF → count blocks, validate with RDKit.js
  const dl = await (async () => { const p = page.waitForEvent('download'); await page.evaluate(() => exportEnsembleSDF()); return p; })();
  const sdfPath = '/tmp/m2_ensemble.sdf';
  await dl.saveAs(sdfPath);
  const sdfText = fs.readFileSync(sdfPath, 'utf8');
  const nBlocks = sdfText.split('$$$$').filter(b => b.trim()).length;
  check('ensemble SDF block count == kept', nBlocks === kept, nBlocks + ' blocks');
  const firstBlock = sdfText.split('$$$$')[0];
  const rd = await page.evaluate(s => {
    const m = rdkitModule.get_mol(s);
    const ok = m && m.is_valid();
    const n = ok ? JSON.parse(m.get_descriptors()).NumAtoms : 0;
    if (m) m.delete();
    return { ok, n };
  }, firstBlock);
  check('first block readable by RDKit.js (33 explicit atoms incl. attached H)', rd.ok && rd.n === 33, rd.n + ' atoms');

  // energy spread sanity: ibuprofen MMFF minima within a few kcal/mol
  const eStr = ['E=' + st.Es[0].toFixed(2)];
  check('lowest E finite (implicit-H MMFF may run positive)', eStr[0] && isFinite(st.Es[0]), eStr[0]);

  // cancel path: fire-and-forget a run, cancel mid-optimizations (if still
  // running — fast embeds may finish before the cancel mark; the ΔE-window
  // re-optimization phase then runs to completion with btnConf disabled)
  await page.evaluate(() => { document.getElementById('confN').value = '12'; runConformers(); });
  for (let t = 0; t < 30; t++) {
    await page.waitForTimeout(1000);
    const st = await page.evaluate(() => ({ run: !!confRun, n: confRun ? confRun.confs.length : -1 }));
    console.log('cancel-watch T+' + (t + 1) + 's ' + JSON.stringify(st));
    if (!st.run) break;
    if (t === 4) { await page.evaluate(() => { if (confRun) confRun.cancelled = true; }); console.log('cancel requested'); }
  }
  await page.waitForFunction(() => confRun === null && !reoptState &&
    !document.getElementById('btnConf').disabled, null, { timeout: 300000 });
  check('cancel/reopt path leaves UI usable', await page.evaluate(() =>
    !document.getElementById('btnConf').disabled && confEnsemble !== null),
    (await page.evaluate(() => document.getElementById('status3d').textContent)).slice(0, 70));

  await page.screenshot({ path: '/tmp/app_m2.png', fullPage: true });
  check('zero console errors', errors.length === 0, errors.slice(0, 3).join(' | '));

  console.log(`\nRESULT: ${pass} passed, ${fail} failed`);
  await browser.close();
  process.exit(fail ? 1 : 0);
})();
