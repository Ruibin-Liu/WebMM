// m6_scratch.test.js — M4 P3 (docs §5): the molecule clipboard moves to a
// global webmm-scratch IDB with one-time legacy migration, survives reload,
// and gains the scratch→project send gesture. Plus: trashed-project purge.
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
  page.on('dialog', d => d.accept());

  // seed the legacy localStorage history BEFORE first load (migration input)
  await page.addInitScript(() => {
    localStorage.setItem('wb-compoundHistory', JSON.stringify([
      { smiles: 'CCO', name: '乙醇(草稿)', timestamp: '2026/10/18 09:00:00' },
      { smiles: 'c1ccccc1', name: '苯(草稿)', timestamp: '2026/10/18 09:05:00' },
      { smiles: 'CCO', name: '乙醇重复', timestamp: '2026/10/18 09:10:00' },   // dedup input
    ]));
  });
  await page.goto(URL, { waitUntil: 'load' });
  await page.waitForFunction(() => window.webmm !== undefined, null, { timeout: 60000 });

  // ---- migration: history modal shows the deduped legacy entries ----
  await page.evaluate(() => showHistory());
  await page.waitForFunction(() => document.getElementById('historyList').children.length === 2, null, { timeout: 20000 });
  const hist = await page.evaluate(() => [...document.querySelectorAll('#historyList strong')].map(e => e.textContent));
  check('legacy localStorage migrated (deduped by SMILES)', hist.includes('乙醇(草稿)') && hist.includes('苯(草稿)'), JSON.stringify(hist));

  // ---- new entries persist in the scratch DB, not localStorage ----
  await page.evaluate(() => { document.getElementById('input').value = 'CCCO'; process(true); });   // NOT in the demo library
  await page.waitForFunction(() => !!document.querySelector('#svg svg'), null, { timeout: 60000 });
  await page.evaluate(() => showHistory());
  await page.waitForFunction(() => document.querySelectorAll('#historyList .hist-item').length === 3, null, { timeout: 20000 });
  const ls = await page.evaluate(() => (localStorage.getItem('wb-compoundHistory') || '').includes('CCCO'));
  check('new entries live in webmm-scratch (localStorage untouched)', !ls);

  // ---- reload: scratch survives (IDB), migration does not re-import ----
  await page.reload({ waitUntil: 'load' });
  await page.waitForFunction(() => window.webmm !== undefined, null, { timeout: 60000 });
  await page.evaluate(() => showHistory());
  await page.waitForFunction(() => document.querySelectorAll('#historyList .hist-item').length === 3, null, { timeout: 20000 });
  check('scratch survives reload (IDB-backed)', true);

  // ---- send-to-library: scratch → current project ----
  await page.evaluate(() => { switchMode('search'); loadDemoLibrary(); });
  await page.waitForFunction(() => (window.__search.getState() || []).length === 55, null, { timeout: 60000 });
  await page.evaluate(() => showHistory());
  await page.evaluate(() => scratchToLibrary(0));   // newest entry (CCCO, not in demo lib)
  await page.waitForFunction(() => (window.__search.getState() || []).length === 56, null, { timeout: 30000 });
  const libHas = await page.evaluate(() => (window.__search.getState() || []).some(e => e.smiles === 'CCCO'));
  check('send-to-library appends into the project (55 → 56, persisted)', libHas);

  // ---- clear wipes the scratch DB ----
  await page.evaluate(() => showHistory());
  await page.evaluate(() => clearHistory());
  await page.waitForFunction(() => document.querySelectorAll('#historyList .hist-item').length === 0, null, { timeout: 20000 });
  check('clear wipes webmm-scratch', true);

  check('no page errors', errors.length === 0, errors.slice(0, 2).join(' | '));
  await browser.close();
  console.log('RESULT: ' + pass + ' passed, ' + fail + ' failed');
  process.exit(fail ? 1 : 0);
})().catch(e => { console.error(e); process.exit(1); });
