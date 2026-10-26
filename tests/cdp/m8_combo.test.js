// m8_combo.test.js — Top-2 color deepening (docs/assay-foundation.md followup):
// the shape↔color combo weight slider (fixed-pose re-rank vs re-alignment),
// per-feature pharmacophore contribution tooltips, ensemble-fairness surface.
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
  await page.waitForFunction(() => window.webmm !== undefined, null, { timeout: 60000 });
  await page.evaluate(() => { indexedDB.deleteDatabase('webmm-platform'); indexedDB.deleteDatabase('webmm-projects'); });
  await page.reload({ waitUntil: 'load' });
  await page.waitForFunction(() => document.getElementById('rdkitVersion').textContent !== 'Loading...', null, { timeout: 60000 });
  await page.evaluate(() => switchMode('search'));

  await page.evaluate(async () => {
    await new Promise(r => loadSearchLibraryFrom([
      { smiles: 'CC(O)c1ccccc1', name: 'phenol-adj' },       // donor+acceptor+ring (query neighbor)
      { smiles: 'CC(=O)Nc1ccccc1', name: 'acetanilide' },    // donor+acceptor+amides+ring
      { smiles: 'CC(=O)Oc1ccccc1', name: 'phenyl-acetate' }, // acceptor+ring (no donor)
      { smiles: 'c1ccc2ccccc2c1', name: 'naphthalene' },     // ring-only, no polarity
    ]).persisted.then(r));
  });
  await page.waitForTimeout(500);

  // ---- 1. slider present, default 50/50, realign button hidden at default ----
  const ui0 = await page.evaluate(() => ({
    hasSlider: !!document.getElementById('shapeColorW'),
    readout: document.getElementById('shapeColorWVal').textContent,
    realignHidden: document.getElementById('shapeRealignBtn').style.display === 'none',
  }));
  check('slider present, default 50/50, realign hidden at default',
    ui0.hasSlider && /50\/50/.test(ui0.readout) && ui0.realignHidden, JSON.stringify(ui0));

  // column helper: map a visible header label to its row-cell ordinal
  const colOf = (label) => page.evaluate((lb) => {
    const ths = [...document.querySelectorAll('#searchTable thead th')].filter(th => th.style.display !== 'none');
    const i = ths.findIndex(th => th.textContent.includes(lb));
    return i;
  }, label);
  // ---- 2. shape search (ensemble path nConfs=2 for conf surface) ----
  await page.evaluate(() => segModeSet('shape'));
  await page.evaluate(() => { document.getElementById('searchQuery').value = 'CC(O)c1ccccc1'; document.getElementById('shapeConfs').value = '2'; });
  await page.evaluate(() => runSearch());
  await page.waitForFunction(() => document.querySelectorAll('#searchRows tr').length >= 2, null, { timeout: 90000 });
  await page.waitForTimeout(400);

  const iShape = await colOf('Shape T'), iColor = await colOf('Color T'), iCombo = await colOf('Combo'), iPharm = await colOf('Pharm');
  const rows0 = await page.evaluate(([iS, iC, iB]) => [...document.querySelectorAll('#searchRows tr')].map(tr => {
    const tds = tr.querySelectorAll('td');
    return { shape: tds[iS] && tds[iS].textContent, color: tds[iC] && tds[iC].textContent, combo: tds[iB] && tds[iB].textContent };
  }), [iShape, iColor, iCombo]);
  check('shape+color+combo columns render', rows0.every(r => r.shape && r.color && r.combo), JSON.stringify(rows0));

  // ---- 3. w=1: combo column equals shape values (fixed-pose re-rank) ----
  await page.evaluate(() => { document.getElementById('shapeColorW').value = '1'; onComboWChange(); });
  await page.waitForTimeout(300);
  const rows1 = await page.evaluate(([iS, iB]) => ({
    label: document.getElementById('searchResultStatus').textContent,
    combo: [...document.querySelectorAll('#searchRows tr')].map(tr => tr.querySelectorAll('td')[iB].textContent),
    shape: [...document.querySelectorAll('#searchRows tr')].map(tr => tr.querySelectorAll('td')[iS].textContent),
    realignVisible: document.getElementById('shapeRealignBtn').style.display !== 'none',
  }), [iShape, iCombo]);
  check('w=1: status line marks the re-rank as fixed-pose (not re-aligned)',
    /重排/.test(rows1.label) && /固定姿势/.test(rows1.label) && /Shape 100\/Color 0/.test(rows1.label), rows1.label.slice(0, 120));
  check('w=1: combo column == shape column (pure-shape blend)',
    rows1.combo.length === rows1.shape.length && rows1.combo.every((c, i) => c === rows1.shape[i]),
    JSON.stringify({ combo: rows1.combo, shape: rows1.shape }));
  check('w≠0.5 exposes the re-alignment search button', rows1.realignVisible);

  // ---- 4. w=0.25: combo = 0.25*shape + 0.75*color on the displayed percents ----
  await page.evaluate(() => { document.getElementById('shapeColorW').value = '0.25'; onComboWChange(); });
  await page.waitForTimeout(300);
  const blend = await page.evaluate(([iS, iC, iB]) => [...document.querySelectorAll('#searchRows tr')].slice(0, 3).map(tr => {
    const tds = tr.querySelectorAll('td');
    const pct = s => parseFloat(String(s).replace('%', ''));
    const sh = pct(tds[iS].textContent), co = tds[iC].textContent === '—' ? null : pct(tds[iC].textContent), cb = pct(tds[iB].textContent);
    return { sh, co, cb, expect: co == null ? sh : 0.25 * sh + 0.75 * co };
  }), [iShape, iColor, iCombo]);
  check('w=0.25: combo = 25% shape + 75% color (missing color falls back to shape)',
    blend.every(b => Math.abs(b.cb - b.expect) < 0.15), JSON.stringify(blend));

  // ---- 5. per-feature contribution tooltip (donor/acceptor ✓/✗ with ratios) ----
  const pharmTip = await page.evaluate((iP) => {
    const tr = document.querySelector('#searchRows tr');
    const p = tr && tr.querySelectorAll('td')[iP];
    return p && p.title ? p.title : '';
  }, iPharm);
  check('pharm column tooltip carries per-feature ✓/✗ with match ratios',
    /donor/.test(pharmTip) && /✓|✗/.test(pharmTip) && /\d\.\d\d/.test(pharmTip), pharmTip.slice(0, 100));

  // ---- 6. ensemble-fairness surface: conf tooltip names effective conformers ----
  const confTip = await page.evaluate(() => {
    const tr = document.querySelector('#searchRows tr');
    const confTd = tr && tr.querySelectorAll('td')[4];
    return confTd && confTd.title ? confTd.title : '';
  });
  check('conf tooltip: winning conformer + effective ensemble size + fairness note',
    /获胜构象/.test(confTip) && /有效构象/.test(confTip) && /采样优势/.test(confTip), confTip);

  check('zero page errors', errors.length === 0, errors.slice(0, 3).join(' | '));
  await browser.close();
  console.log('RESULT: ' + pass + ' passed, ' + fail + ' failed');
  process.exit(fail ? 1 : 0);
})().catch(e => { console.error('ERR', e); process.exit(1); });
