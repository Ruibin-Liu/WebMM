// m6_csv_export.test.js — CSV export hardening (formula-injection neutralization
// + RFC 4180 quoting). Both pages share the csvField() helper; unit checks run
// against each page, plus one real download flow on platform.html with
// hostile-named molecules. (index.html export paths are already covered by
// m4/m5 content assertions; its call sites are byte-identical to the unit-
// tested helper.)
const { chromium } = require('/opt/homebrew/lib/node_modules/playwright');
const fs = require('fs');
const EXE = '/Users/rliu/Library/Caches/ms-playwright/chromium_headless_shell-1243/chrome-headless-shell-mac-arm64/chrome-headless-shell';
const PAGES = ['http://localhost:8901/app/platform.html', 'http://localhost:8901/app/index.html'];

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

  const CASES = [
    // [input, expectation lambda, label]
    ['=SUM(1+1)*cmd|\' /C calc\'!A0', v => v.startsWith("'="), 'formula = neutralized'],
    ['+CONCAT(1,2)', v => v.replace(/^"|"$/g, '').startsWith("'+"), 'formula + neutralized'],
    ['-2+cmd|', v => v.startsWith("'-"), 'formula - neutralized'],
    ['@SUM(A1)', v => v.startsWith("'@"), 'formula @ neutralized'],
    ['\tTAB', v => v.startsWith("'\t") || v.startsWith("'\t"), 'tab-prefixed neutralized'],
    ['a,b', v => v === '"a,b"', 'comma quoted'],
    ['say "hi"', v => v === '"say ""hi"""', 'quotes escaped'],
    ['line1\nline2', v => v === '"line1\nline2"', 'newline stays one field'],
    ['line1\r\nline2', v => v.startsWith('"line1\r\n'), 'CRLF stays one field'],
    ['plain benzene', v => v === 'plain benzene', 'plain name untouched'],
    [null, v => v === '', 'null -> empty'],
    [123, v => v === '123', 'number as string'],
  ];

  for (const url of PAGES) {
    const tag = url.split('/').pop();
    await page.goto(url, { waitUntil: 'load' });
    await page.waitForFunction(() => typeof csvField === 'function', null, { timeout: 60000 });
    for (const [input, expect, label] of CASES) {
      const got = await page.evaluate(v => csvField(v), input);
      check(`${tag}: ${label}`, expect(got), JSON.stringify(got).slice(0, 60));
    }
  }

  // real flow: hostile names in the library -> sim search -> results CSV download
  await page.goto(PAGES[0], { waitUntil: 'load' });
  await page.waitForFunction(() => window.webmm !== undefined && typeof loadSearchLibraryFrom === 'function', null, { timeout: 60000 });
  await page.evaluate(() => {
    switchMode('search');
    loadSearchLibraryFrom([
      { smiles: 'CCO', name: "=cmd|' /C calc'!A0" },
      { smiles: 'CCN', name: 'plain ethylamine' },
      { smiles: 'CCCO', name: 'comma,name "quoted"' },
    ]);
    document.getElementById('searchMode').value = 'sim'; onSearchModeChange();
    document.getElementById('searchQuery').value = 'CCO';
    runSearch();
  });
  await page.waitForFunction(() => /at T|sorted|Done/i.test(document.getElementById('searchResultStatus').textContent), null, { timeout: 120000 });
  console.log('   status:', await page.evaluate(() => document.getElementById('searchResultStatus').textContent));
  const dl = page.waitForEvent('download', { timeout: 15000 });
  await page.evaluate(() => exportResultsCSV());
  await (await dl).saveAs('/tmp/csv-hardening.csv');
  const csv = fs.readFileSync('/tmp/csv-hardening.csv', 'utf8');
  const lines = csv.split('\n');
  check('download: hostile name neutralized (leading apostrophe)', csv.includes("'=cmd|"), csv.slice(0, 80));
  check('download: NO line begins with a live formula', !lines.some(l => /^[=+@]/.test(l)), lines.slice(0, 3).join(' | ').slice(0, 100));
  check('download: comma+quote name intact as one field', csv.includes('"comma,name ""quoted"""'));
  check('download: plain name row present', /plain ethylamine/.test(csv));
  check('download: no page errors during flow', errors.length === 0, errors.join(' | '));

  await browser.close();
  console.log('RESULT: ' + pass + ' passed, ' + fail + ' failed');
  process.exit(fail ? 1 : 0);
})().catch(e => { console.error(e); process.exit(1); });
