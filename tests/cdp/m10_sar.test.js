// m10_sar.test.js — Top-4: fixed-core R-group table + transparent MPO.
// Seed library with a toluene series; core = benzene; assert the R table
// (positions, fragment SMILES re-canonicalized, unmatched NOT silently
// dropped), MPO transparency (raw + per-criterion score + total), the
// activity criterion from the M5 assay store, and CSV export.
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
  page.on('pageerror', e => errors.push('PAGEERROR: ' + e.message.slice(0, 200)));

  await page.goto(URL + '?m=search', { waitUntil: 'load' });
  await page.waitForFunction(() => document.getElementById('rdkitVersion').textContent !== 'Loading...', null, { timeout: 60000 });
  await page.evaluate(() => { indexedDB.deleteDatabase('webmm-platform'); indexedDB.deleteDatabase('webmm-projects'); });
  await page.reload({ waitUntil: 'load' });
  await page.waitForFunction(() => document.getElementById('rdkitVersion').textContent !== 'Loading...', null, { timeout: 60000 });
  await page.evaluate(() => switchMode('search'));

  await page.evaluate(async () => {
    await new Promise(r => loadSearchLibraryFrom([
      { smiles: 'Cc1ccccc1', name: 'toluene' },
      { smiles: 'Cc1ccc(C)cc1', name: 'p-xylene' },
      { smiles: 'COc1ccccc1', name: 'anisole' },
      { smiles: 'Oc1ccccc1', name: 'phenol' },
      { smiles: 'CCO', name: 'ethanol' },          // does not match the core
    ]).persisted.then(r));
  });
  await page.waitForTimeout(600);
  // one real search first — searchPanel (and its cards) render on results
  await page.evaluate(() => { document.getElementById('searchQuery').value = 'Cc1ccccc1'; runSearch(); });
  await page.waitForFunction(() => document.querySelectorAll('#searchRows tr').length >= 1, null, { timeout: 60000 });
  await page.waitForTimeout(300);

  // ---- 1. R-group table ----
  await page.fill('#sarCore', 'c1ccccc1');
  await page.evaluate(() => runSarAnalysis());
  await page.waitForTimeout(400);
  const sar1 = await page.evaluate(() => ({
    status: document.getElementById('sarStatus').textContent,
    unmatched: document.getElementById('sarUnmatched').textContent,
    header: [...document.querySelectorAll('#sarTable thead th')].map(th => th.textContent),
    rows: [...document.querySelectorAll('#sarTable tbody tr')].map(tr => [...tr.querySelectorAll('td')].map(td => td.textContent.trim())),
    tips: [...document.querySelectorAll('#sarTable tbody tr td[title]')].map(td => td.title),
  }));
  check('status: 4/5 matched, 1 R position + p-xylene second', /匹配 4\/5/.test(sar1.status) && /R 位 2/.test(sar1.status), sar1.status);
  check('unmatched listed (ethanol NOT silently dropped)', /未匹配 1/.test(sar1.unmatched) && /ethanol/.test(sar1.unmatched), sar1.unmatched);
  check('R header: R1 + R2 columns', sar1.header.includes('R1') && sar1.header.includes('R2'), JSON.stringify(sar1.header));
  const tolueneRow = sar1.rows.find(r => r.some(c => /toluene/.test(c)));
  const phenolRow = sar1.rows.find(r => r.some(c => /phenol/.test(c)));
  const anisoleRow = sar1.rows.find(r => r.some(c => /anisole/.test(c)));
  check('toluene R1 = methyl (canonicalized fragment)', tolueneRow && /C/.test(tolueneRow[2] || ''), JSON.stringify(tolueneRow));
  check('phenol R1 = O', phenolRow && /O/.test(phenolRow[2] || ''), JSON.stringify(phenolRow));
  check('anisole R1 = O(C) fragment', anisoleRow && /OC|CO/.test(anisoleRow[2] || ''), JSON.stringify(anisoleRow));
  check('raw [*] fragment SMILES preserved in tooltips', sar1.tips.some(t => t.includes('[*]')), sar1.tips.slice(0, 2).join(' | '));

  // ---- 2. MPO transparency ----
  const mpo1 = await page.evaluate(() => ({
    header: [...document.querySelectorAll('#mpoTable thead th')].map(th => th.textContent),
    firstRow: (document.querySelector('#mpoTable tbody tr') || { textContent: '' }).textContent.trim(),
    nRows: document.querySelectorAll('#mpoTable tbody tr').length,
  }));
  check('MPO columns: MW/cLogP/TPSA/HBD + total', ['MW', 'cLogP', 'TPSA', 'HBD', 'MPO'].every(h => mpo1.header.some(x => x.includes(h))), JSON.stringify(mpo1.header));
  check('MPO transparency: raw + (score) + total + n/4 in every row', /\(/.test(mpo1.firstRow) && /\d\.\d{3}/.test(mpo1.firstRow) && /\/4/.test(mpo1.firstRow), mpo1.firstRow.slice(0, 90));
  check('MPO rows = matched entries only', mpo1.nRows === 4, mpo1.nRows);

  // ---- 3. activity criterion from the M5 assay store ----
  await page.evaluate(() => activityImportUI());
  await page.fill('#actImpName', 'IC50');
  await page.fill('#actImpUnit', 'µM');
  await page.fill('#actImpCsv', 'name,value\nphenol,0.8\nanisole,12\n');
  await page.evaluate(() => activityImportPreview());
  await page.evaluate(() => activityImportCommit());
  await page.waitForTimeout(500);
  await page.check('#sarUseActivity');
  await page.evaluate(() => runSarAnalysis());
  await page.waitForTimeout(300);
  const mpo2 = await page.evaluate(() => ({
    header: [...document.querySelectorAll('#mpoTable thead th')].map(th => th.textContent),
    phenol: (() => { const tr = [...document.querySelectorAll('#mpoTable tbody tr')].find(t => /phenol/.test(t.textContent)); return tr ? tr.textContent : null; })(),
  }));
  check('activity criterion appears when assay selected', mpo2.header.some(h => /活性/.test(h)), JSON.stringify(mpo2.header));
  check('phenol activity raw = 0.8 with score 1.0 (best in set)', !!mpo2.phenol && /0\.80/.test(mpo2.phenol), mpo2.phenol);

  // ---- 4. export ----
  const dl = page.waitForEvent('download', { timeout: 15000 }).catch(() => null);
  await page.evaluate(() => sarExportCsv());
  const d = await dl;
  check('CSV export downloads (sar-mpo.csv)', d && /sar-mpo\.csv/.test(d.suggestedFilename()), d && d.suggestedFilename());

  // ---- 5. MCS guess fills the core input (advisory, not truth) ----
  await page.evaluate(() => sarGuessMcs());
  const mcs = await page.evaluate(() => ({ v: document.getElementById('sarCore').value, s: document.getElementById('sarStatus').textContent }));
  check('MCS guess fills the input with an advisory note', mcs.v.length > 3 && /建议/.test(mcs.s), JSON.stringify(mcs).slice(0, 100));

  check('zero page errors', errors.length === 0, errors.slice(0, 3).join(' | '));
  await browser.close();
  console.log('RESULT: ' + pass + ' passed, ' + fail + ' failed');
  process.exit(fail ? 1 : 0);
})().catch(e => { console.error('ERR', e); process.exit(1); });
