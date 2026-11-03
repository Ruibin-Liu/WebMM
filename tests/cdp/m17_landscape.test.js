// m17_landscape.test.js — energy landscape view (iMTD runs only): the
// disconnectivity graph + basin statistics render inside the conformer panel,
// leaves/table link to the conformer selection, the honest-bias banner is
// present, and the ETKDG conformers path hides the section.
const { chromium } = require('/opt/homebrew/lib/node_modules/playwright');
const EXE = '/Users/rliu/Library/Caches/ms-playwright/chromium_headless_shell-1243/chrome-headless-shell-mac-arm64/chrome-headless-shell';
const URL = 'http://localhost:8901/app/index.html';

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

  await page.goto(URL, { waitUntil: 'load' });
  await page.waitForFunction(() => window.webmm !== undefined, null, { timeout: 120000 });

  const loadMol = async smi => {
    await page.evaluate(s => { document.getElementById('input').value = s; process(true); }, smi);
    await page.waitForFunction(() => !!document.querySelector('#svg svg'), null, { timeout: 60000 });
    await page.evaluate(() => embed3D());
    await page.waitForFunction(() => !!document.querySelector('#viewer3d canvas'), null, { timeout: 90000 });
  };

  // ================= 1. iMTD run produces the landscape section =================
  await loadMol('CCO ethanol');
  await page.evaluate(() => {
    document.getElementById('imtdRounds').value = '2';
    document.getElementById('imtdSteps').value = '800';
    runImtd();
  });
  await page.waitForFunction(() => /unique minima/.test(document.getElementById('status3d').textContent), null, { timeout: 180000 });
  const st = await page.evaluate(() => {
    const wrap = document.getElementById('landscapeWrap');
    const svg = document.getElementById('dgraphSvg');
    const L = lastLandscape;
    const occ = (L && L.frames) ? L.basins.reduce((a, b) => a + b.visits, 0) / L.frames : 0;
    return {
      visible: wrap && wrap.style.display !== 'none',
      leaves: svg ? svg.querySelectorAll('circle[data-basin]').length : 0,
      svgLeavesAttr: svg ? svg.dataset.leaves : '',
      n: confEnsemble.length,
      frames: L ? L.frames : 0,
      pruned: L ? L.pruned : 0,
      occSum: occ,
      tableRows: document.querySelectorAll('#basinTableWrap tbody tr').length,
      banner: wrap ? wrap.textContent : '',
      meta: document.getElementById('landscapeMeta').textContent,
      hasEdges: svg ? svg.querySelectorAll('path').length : 0,
    };
  });
  check('landscape section visible after iMTD', st.visible);
  check('leaf count == conformer count', st.leaves === st.n && +st.svgLeavesAttr === st.n, `leaves=${st.leaves} attr=${st.svgLeavesAttr} n=${st.n}`);
  check('basin table rows == conformer count', st.tableRows === st.n, `${st.tableRows} vs ${st.n}`);
  check('frames recorded (≥ 15), occupancy sums ≤ 100%', st.frames >= 15 && st.occSum <= 1.0 + 1e-9, `frames=${st.frames} occ=${(st.occSum * 100).toFixed(0)}%`);
  check('transition edges drawn', st.hasEdges >= 1, `paths=${st.hasEdges} · ${st.meta}`);
  check('honest banner: metad bias + qualitative + not saddle', /metad/.test(st.banner) && /定性/.test(st.banner) && /鞍点|上界/.test(st.banner));

  // ================= 2. leaf click → conformer selection linkage =================
  const linked = await page.evaluate(() => {
    const leaf = document.querySelector('#dgraphSvg circle[data-basin="1"]');
    if (!leaf) return { ok: false, why: 'no leaf #1' };
    leaf.dispatchEvent(new MouseEvent('click', { bubbles: true }));
    return { ok: true, sel: selConf, status: document.getElementById('status3d').textContent.slice(0, 60) };
  });
  check('leaf click selects conformer #2 (selConf=1)', linked.ok && linked.sel === 1, JSON.stringify(linked));
  const highlighted = await page.evaluate(() => {
    const sel = document.querySelector('#dgraphSvg circle[data-basin="1"]');
    return sel && sel.getAttribute('r') === '6' && sel.getAttribute('fill') === '#b45309';
  });
  check('clicked leaf highlighted (selection ring)', !!highlighted);

  // ================= 3. ETKDG conformers path hides the landscape =================
  await page.evaluate(() => { document.getElementById('confN').value = '10'; runConformers(); });
  await page.waitForFunction(() => /Conformers: \d+ optimized/.test(document.getElementById('status3d').textContent), null, { timeout: 180000 });
  const hidden = await page.evaluate(() => document.getElementById('landscapeWrap').style.display === 'none' && lastLandscape === null);
  check('ETKDG conformers path hides the landscape (no trace)', hidden);
  const imtdBack = await page.evaluate(() => {
    runImtd();
    return true;
  });
  await page.waitForFunction(() => /unique minima/.test(document.getElementById('status3d').textContent), null, { timeout: 180000 });
  const backVis = await page.evaluate(() => document.getElementById('landscapeWrap').style.display !== 'none');
  check('re-running iMTD brings the landscape back', backVis && imtdBack);

  check('zero page errors', errors.length === 0, errors.join(' | '));

  console.log(`\nRESULT: ${pass} passed, ${fail} failed`);
  await browser.close();
  process.exit(fail ? 1 : 0);
})().catch(e => { console.error('FATAL:', e); process.exit(1); });
