// M4 regression: batch mode — descriptors, drug rules, CSV export, 3D + SDF
// provenance export. Run: node m4_batch.test.js (needs the :8901 test server)
const { chromium } = require('/opt/homebrew/lib/node_modules/playwright');
const fs = require('fs');
const EXE = '/Users/rliu/Library/Caches/ms-playwright/chromium_headless_shell-1234/chrome-headless-shell-mac-arm64/chrome-headless-shell';

const INPUT = [
  'CN1C=NC2=C1C(=O)N(C)C(=O)N2 caffeine',
  'CC(=O)Oc1ccccc1C(=O)O aspirin',
  'Cn1cnc2c1c(=O)n(C)c(=O)n2C theophylline',
].join('\n');

(async () => {
  const browser = await chromium.launch({ executablePath: EXE });
  const page = await browser.newPage({ viewport: { width: 1440, height: 900 } });
  const errors = [];
  page.on('pageerror', e => errors.push('PAGEERROR: ' + e.message.slice(0, 120)));

  let pass = 0, fail = 0;
  const check = (name, cond, detail = '') => {
    if (cond) { pass++; console.log(`  ✓ ${name}${detail ? ' — ' + detail : ''}`); }
    else { fail++; console.log(`  ✗ ${name}${detail ? ' — ' + detail : ''}`); }
  };

  await page.goto('http://localhost:8901/app/index.html', { waitUntil: 'load' });
  await page.waitForFunction(() => window.webmmReady === true, null, { timeout: 60000 });

  // 1. switch to batch tab, enter 3 molecules
  const tabs = await page.evaluate(() => {
    const tabs = [...document.querySelectorAll('.mode-tab')];
    const batch = tabs.find(t => /batch/i.test(t.textContent));
    if (batch) batch.click();
    return tabs.map(t => t.textContent.trim());
  });
  await page.evaluate(s => { document.getElementById('input').value = s; }, INPUT);
  await page.waitForTimeout(200);

  // 2. run batch (descriptors only)
  await page.evaluate(() => runBatch());
  await page.waitForFunction(() => {
    const rows = [...document.querySelectorAll('#batchRows tr, #batchTable tr')];
    return rows.length >= 3 && !batchRun;
  }, null, { timeout: 120000 });
  const rows1 = await page.evaluate(() => [...document.querySelectorAll('#batchRows tr, #batchTable tr')]
    .map(tr => tr.textContent.trim()).slice(0, 6));
  check('batch renders 3 rows', rows1.length >= 3, rows1.length + ' rows');
  check('rows carry names + descriptors', rows1.some(r => r.includes('caffeine')) && rows1.some(r => r.includes('aspirin')) && rows1.some(r => r.includes('MW') || r.includes('194') || r.includes('180')), rows1[1]?.slice(0, 50));

  // 3. CSV export
  const dl = await (async () => { const p = page.waitForEvent('download'); await page.evaluate(() => exportBatchCSV()); return p; })();
  await dl.saveAs('/tmp/batch_test.csv');
  const csv = fs.readFileSync('/tmp/batch_test.csv', 'utf8');
  check('CSV export has header + 3 rows', csv.split('\n').filter(l => l.trim()).length >= 4 && /name|smiles/i.test(csv.split('\n')[0]), csv.split('\n')[0].slice(0, 60));
  check('CSV includes caffeine row', /caffeine/.test(csv));

  // 4. 3D batch run (2 molecules: drop theophylline for speed) + SDF provenance
  await page.evaluate((s) => { document.getElementById('input').value = s.split('\n').slice(0, 2).join('\n'); }, INPUT);
  await page.evaluate(() => { document.getElementById('batch3d').checked = true; });
  await page.evaluate(() => runBatch());
  await page.waitForFunction(() => !batchRun, null, { timeout: 300000 });
  const sdfReady = await page.evaluate(() => batchResults.filter(r => r.sdf).length);
  check('3D batch produces SDFs', sdfReady >= 2, sdfReady + ' sdf rows');
  const dl2 = await (async () => { const p = page.waitForEvent('download'); await page.evaluate(() => exportBatchSDF()); return p; })();
  await dl2.saveAs('/tmp/batch_test.sdf');
  const sdf = fs.readFileSync('/tmp/batch_test.sdf', 'utf8');
  check('batch SDF carries provenance', sdf.includes('WEBMM_PROGRAM') && sdf.includes('WEBMM_ENGINE') && sdf.includes('$$$$'), 'WEBMM_PROGRAM/ENGINE present');

  // 5. zero page errors
  check('zero page errors', errors.length === 0, errors.slice(0, 2).join(' | '));

  console.log(`\nRESULT: ${pass} passed, ${fail} failed`);
  await browser.close();
  process.exit(fail ? 1 : 0);
})();
