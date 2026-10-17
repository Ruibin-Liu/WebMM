// m6_analog_ui.test.js — explorer UI round: thumbnail rendering quality
// (RDKit renders at the display scale), process-label tooltips, and
// CSV/SDF export of the results table.
const { chromium } = require('/opt/homebrew/lib/node_modules/playwright');
const fs = require('fs');
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

  await page.goto(URL, { waitUntil: 'load' });
  await page.waitForFunction(() => window.webmm !== undefined, null, { timeout: 60000 });

  // ---- tooltips (static, before any run) ----
  const tips = await page.evaluate(() => {
    const pick = el => el && el.title ? el.title.length : 0;
    return {
      swap: pick(document.querySelector('#analogPanel .path-head .card-sub')),
      aza: pick([...document.querySelectorAll('#analogPanel .path-head .card-sub')][1]),
      seg: [...document.querySelectorAll('#segMode .seg-btn')].map(b => b.title.length),
    };
  });
  check('swap-card tooltip is detailed (>=60 chars, mentions the funnel)', tips.swap >= 60, JSON.stringify(tips));
  check('aza-card tooltip is detailed (>=60 chars)', tips.aza >= 60);
  check('all three search-mode seg buttons have tooltips', tips.seg.length === 3 && tips.seg.every(n => n >= 40), JSON.stringify(tips.seg));

  // ---- run a small explorer for thumbnails + exports ----
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

  // row thumbnails fill asynchronously (queueAnalogThumbs, 40ms/tick budget)
  await page.waitForFunction(() => !!document.querySelector('#analogRows .rowthumb svg'), null, { timeout: 20000 }).catch(() => {});
  const rowSvg = await page.evaluate(() => {
    const el = document.querySelector('#analogRows .rowthumb svg');
    return el ? { w: el.viewBox.baseVal.width, box: document.querySelector('#analogRows .rowthumb').getBoundingClientRect().width } : null;
  });
  check('row thumbnail renders at 128 wide (2x retina of 64px box)', rowSvg && Math.abs(rowSvg.w - 128) < 1 && Math.abs(rowSvg.box - 64) < 1, JSON.stringify(rowSvg));

  // hover preview: 400x272 render shown at exactly 200x136
  await page.hover('#analogRows tr td.cell-smiles');
  await page.waitForTimeout(300);
  const pop = await page.evaluate(() => {
    const el = document.querySelector('.thumb-pop svg');
    if (!el) return null;
    const box = el.getBoundingClientRect();
    return { w: el.viewBox.baseVal.width, h: el.viewBox.baseVal.height, bw: box.width, bh: box.height };
  });
  check('hover preview renders at 400x272', pop && Math.abs(pop.w - 400) < 1 && Math.abs(pop.h - 272) < 1, JSON.stringify(pop));
  check('hover preview displays at exactly 200x136 (aspect-true)', pop && Math.abs(pop.bw - 200) < 1 && Math.abs(pop.bh - 136) < 1, JSON.stringify(pop));
  await page.mouse.move(5, 5);

  // ---- exports ----
  const rowCount = await page.evaluate(() => document.querySelectorAll('#analogRows tr').length);
  let dl = page.waitForEvent('download', { timeout: 15000 });
  await page.evaluate(() => exportAnalogCSV());
  await (await dl).saveAs('/tmp/analog_test.csv');
  const csv = fs.readFileSync('/tmp/analog_test.csv', 'utf8');
  const csvRows = csv.trim().split('\n');
  check('CSV: header + one row per result', csvRows.length === rowCount + 1 && /^SMILES,ShapeT,FlexT,MW,RotB,SA/.test(csvRows[0]), csvRows[0]);
  check('CSV: T column is numeric', rowCount > 0 && /^"?[CCN][^,]*"?,\d\.\d{4}/.test(csvRows[1]), csvRows[1].slice(0, 50));

  dl = page.waitForEvent('download', { timeout: 15000 });
  await page.evaluate(() => exportAnalogSDF());
  await (await dl).saveAs('/tmp/analog_test.sdf');
  const sdf = fs.readFileSync('/tmp/analog_test.sdf', 'utf8');
  const v2k = (sdf.match(/V2000/g) || []).length;
  const dollars = (sdf.match(/\$\$\$\$/g) || []).length;
  check('SDF: one V2000 record per result, $$$$ separated', v2k === rowCount && dollars === rowCount, 'v2k=' + v2k + ' $=' + dollars + ' rows=' + rowCount);
  check('SDF: metric data fields present', /<ShapeT>/.test(sdf) && /<SA>/.test(sdf));
  check('no page errors', errors.length === 0, errors.slice(0, 2).join(' | '));

  // ---- bidirectional aza-scan (deaza N->CH) ----
  const setupPicoline = async () => {
    await page.evaluate(() => { document.getElementById('input').value = 'Cc1ccncc1'; process(true); });
    await page.waitForFunction(() => !!document.querySelector('#svg svg'), null, { timeout: 60000 });
    await page.evaluate(() => embed3D());
    await page.waitForFunction(() => !!document.querySelector('#viewer3d canvas'), null, { timeout: 60000 });
  };
  const runAza = () => page.evaluate(() => runAzaScan());
  const awaitAzaDone = () => page.waitForFunction(() => /Done: \d+ scored/.test(document.getElementById('analogStatus').textContent), null, { timeout: 180000 });

  // default (aza only): every candidate keeps at least one ring N
  await setupPicoline();
  await runAza();
  await awaitAzaDone();
  const azaOnly = await page.evaluate(() => ({
    rows: [...document.querySelectorAll('#analogRows .sm')].map(e => e.title || e.textContent),
    counts: window.__azaCounts,
  }));
  check('aza-only run: every candidate contains a ring N', azaOnly.rows.length > 0 && azaOnly.rows.every(x => x.includes('n')), azaOnly.rows.slice(0, 3).join(' '));
  check('aza-only run: no deaza counts reported', azaOnly.counts.deazaSingle == null, JSON.stringify(azaOnly.counts));

  // deaza checked: toluene (no N) appears, counts reported
  await page.evaluate(() => { document.getElementById('azaDeaza').checked = true; });
  await runAza();
  await awaitAzaDone();
  const withDeaza = await page.evaluate(() => ({
    rows: [...document.querySelectorAll('#analogRows .sm')].map(e => e.title || e.textContent),
    counts: window.__azaCounts,
    status: document.getElementById('analogStatus').textContent,
  }));
  check('deaza run: N-free candidate (toluene) present', withDeaza.rows.some(x => /^Cc1ccccc1$/.test(x)), withDeaza.rows.slice(0, 4).join(' '));
  check('deaza run: counts reported (4-picoline has exactly 1 ring N)', withDeaza.counts.deazaSingle === 1 && withDeaza.counts.deazaDouble === 0, JSON.stringify(withDeaza.counts));
  check('deaza run: Done line mentions deaza', /deaza 1 \+ 0/.test(withDeaza.status), withDeaza.status.slice(-90));
  await page.evaluate(() => { document.getElementById('azaDeaza').checked = false; });
  check('no page errors (deaza round)', errors.length === 0, errors.slice(0, 2).join(' | '));

  await browser.close();
  console.log('RESULT: ' + pass + ' passed, ' + fail + ' failed');
  process.exit(fail ? 1 : 0);
})().catch(e => { console.error(e); process.exit(1); });
