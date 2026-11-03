// m20_degeneracy_thermo.test.js — degeneracy recovery (mirror twins fold
// with g=2 and lift the Boltzmann share) and the harmonic free-energy tier
// (G column, Boltzmann(G), imaginary-mode flags, temperature refresh).
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

  await loadMol('CCO ethanol');
  await page.evaluate(() => {
    document.getElementById('imtdRounds').value = '2';
    document.getElementById('imtdSteps').value = '800';
    runImtd();
  });
  await page.waitForFunction(() => /unique minima/.test(document.getElementById('status3d').textContent), null, { timeout: 180000 });

  // ================= 1. degeneracy + Boltzmann(E) =================
  const st = await page.evaluate(() => {
    const th = [...document.querySelectorAll('#basinTableWrap th')].map(t => t.textContent);
    const rows = [...document.querySelectorAll('#basinTableWrap tbody tr')].map(r =>
      [...r.children].map(c => c.textContent));
    const L = lastLandscape;
    return { th, rows, g: (L.basins || []).map(b => b.g), n: confEnsemble.length, Es: confEnsemble.map(c => c.E) };
  });
  const gi = st.th.findIndex(t => /^g$/.test(t.trim()));
  const bzi = st.th.findIndex(t => /Boltz\.E/.test(t));
  check('g + Boltz.E + Boltz.G columns present', gi >= 0 && bzi >= 0 && st.th.some(t => /Boltz\.G/.test(t)),
        st.th.join('|'));
  check('g column valid (≥1, sane upper bound on observed wells)', st.g.length === st.n && st.g.every(x => x >= 1 && x <= 24), JSON.stringify(st.g));
  const boltzOk = await page.evaluate(() => {
    const L = lastLandscape;
    const gs = L.basins.map(b => b.g || 1);
    const expect = JSON.parse(window.webmm.imtd_boltzmann_populations_g_json(confEnsemble.map(c => c.E), gs, 298));
    const shown = [...document.querySelectorAll('#basinTableWrap tbody tr')].map(r => parseFloat(r.children[6].textContent));
    return shown.every((v, i) => Math.abs(v - 100 * expect[i]) < 1.0);
  });
  check('Boltz.E column matches g·exp weights from the wasm fn', boltzOk);

  // ================= 2. harmonic tier =================
  const thermo = await page.evaluate(() => {
    const rows = [...document.querySelectorAll('#basinTableWrap tbody tr')].map(r =>
      [...r.children].map(c => c.textContent));
    return { rows, n: confEnsemble.length, E0: confEnsemble[0].E };
  });
  const gCol = thermo.rows.map(r => r[2]);
  check('G column finite and above E (ZPE-dominated)', thermo.rows.length === thermo.n &&
    gCol.every((v, i) => isFinite(parseFloat(v)) && parseFloat(v) > thermo.E0 + 10), gCol.join('|'));
  const bzg = thermo.rows.map(r => parseFloat(r[7]));
  const bzgSum = bzg.reduce((a, b) => a + (isFinite(b) ? b : 0), 0);
  check('Boltz.G populations sum to 100%', Math.abs(bzgSum - 100) < 2, bzg.join('|'));
  const noImag = await page.evaluate(() =>
    ![...document.querySelectorAll('#basinTableWrap span')].some(s => /⚠/.test(s.textContent)));
  check('no imaginary-frequency flags on ethanol minima', noImag);

  // T change refreshes G and Boltzmann
  const hot = await page.evaluate(() => {
    document.getElementById('imtdTemp').value = '500';
    document.getElementById('imtdTemp').onchange();
    const rows = [...document.querySelectorAll('#basinTableWrap tbody tr')].map(r =>
      [...r.children].map(c => c.textContent));
    return rows;
  });
  const hotG = hot.map(r => parseFloat(r[2]));
  check('T=500 K refreshes the harmonic tier (G shifts)', hotG.every(v => isFinite(v)) &&
    hotG.some((v, i) => Math.abs(v - parseFloat(gCol[i])) > 0.05), hotG.join('|'));

  // ================= 3. degeneracy lifts a share when a twin exists =================
  const twinLift = await page.evaluate(() => {
    const L = lastLandscape;
    const gs = L.basins.map(b => b.g || 1);
    if (!gs.some(x => x > 1)) return null;   // no twin folded in this run
    const withG = JSON.parse(window.webmm.imtd_boltzmann_populations_g_json(confEnsemble.map(c => c.E), gs, 298));
    const noG = JSON.parse(window.webmm.imtd_boltzmann_populations_json(confEnsemble.map(c => c.E), 298));
    return { gs, withG, noG };
  });
  check('mirror twin (when folded) lifts its Boltzmann share',
    twinLift === null || twinLift.withG.some((v, i) => v > twinLift.noG[i]),
    twinLift ? JSON.stringify(twinLift.gs) : 'no twin in this deterministic run (allowed)');

  check('zero page errors', errors.length === 0, errors.join(' | '));

  console.log(`\nRESULT: ${pass} passed, ${fail} failed`);
  await browser.close();
  process.exit(fail ? 1 : 0);
})().catch(e => { console.error('FATAL:', e); process.exit(1); });
