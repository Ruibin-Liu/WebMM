// m19_saddle_boltzmann.test.js — saddle refinement (CI-NEB band peaks) over
// the observed transition pairs, and the Boltzmann population reconstruction
// (temperature input, ideal populations over the found minima).
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

  // ================= 1. Boltzmann column =================
  const bz = await page.evaluate(() => {
    const th = [...document.querySelectorAll('#basinTableWrap th')];
    const vals = [...document.querySelectorAll('#basinTableWrap tbody tr')]
      .map(r => parseFloat(r.children[6].textContent));
    return { col: th.some(t => /Boltz\.E/.test(t.textContent)), vals, n: confEnsemble.length };
  });
  const bzSum = bz.vals.reduce((a, b) => a + (isFinite(b) ? b : 0), 0);
  check('Boltzmann column present, populations sum to 100%', bz.col && bz.vals.length === bz.n && Math.abs(bzSum - 100) < 2,
        `${bz.vals} (sum ${bzSum.toFixed(1)})`);
  // temperature effect: hotter T flattens — basin #2 population must rise
  const hot = await page.evaluate(() => {
    document.getElementById('imtdTemp').value = '3000';
    drawLandscape();
    return parseFloat([...document.querySelectorAll('#basinTableWrap tbody tr')][1].children[6].textContent);
  });
  const cold = await page.evaluate(() => {
    document.getElementById('imtdTemp').value = '298';
    drawLandscape();
    return parseFloat([...document.querySelectorAll('#basinTableWrap tbody tr')][1].children[6].textContent);
  });
  check('higher T flattens populations (basin #2 rises)', hot >= cold, `T=3000: ${hot}% vs T=298: ${cold}%`);

  // ================= 2. saddle refinement =================
  const before = await page.evaluate(() => ({
    refined: (lastLandscape.transitions || []).filter(t => t.refined != null).length,
    tableHidden: document.getElementById('transTableWrap').style.display === 'none',
  }));
  check('no refined saddles before the run (table hidden)', before.refined === 0 && before.tableHidden);
  await page.evaluate(() => window.runNebRefine());
  await page.waitForFunction(() => /势垒精修完成|势垒精修失败/.test(document.getElementById('status3d').textContent), null, { timeout: 300000 });
  const after = await page.evaluate(() => {
    const rows = [...document.querySelectorAll('#transTableWrap tbody tr')].map(r => r.textContent);
    const refined = (lastLandscape.transitions || []).filter(t => t.refined != null);
    const okSaddle = refined.every(t => {
      const ea = confEnsemble[t.a] ? confEnsemble[t.a].E : -Infinity;
      const eb = confEnsemble[t.b] ? confEnsemble[t.b].E : -Infinity;
      return t.refined >= Math.max(ea, eb) - 1e-6;
    });
    return {
      n: refined.length,
      tableShown: document.getElementById('transTableWrap').style.display !== 'none',
      rows, okSaddle,
      tokens: taskScheduler.stats().interactiveHeld === 0,
      btn: document.getElementById('btnNebRefine').textContent,
    };
  });
  check('NEB refined at least one edge', after.n >= 1, `${after.n} refined`);
  check('transition table shown with observed + refined', after.tableShown && after.rows.length >= 1, after.rows.join(' | ').slice(0, 120));
  check('refined saddles bound both endpoint energies', after.okSaddle);
  check('tokens released, button restored', after.tokens && after.btn === '势垒精修(NEB)');

  // ================= 3. GFN-FF honest refusal =================
  await page.selectOption('#engineSel', 'GFNFF');
  await page.evaluate(() => runImtd());
  await page.waitForFunction(() => /unique minima/.test(document.getElementById('status3d').textContent), null, { timeout: 180000 });
  const gf = await page.evaluate(() => {
    window.runNebRefine();
    return document.getElementById('status3d').textContent;
  });
  check('GFN-FF engine refuses refinement with an honest reason', /仅支持 MMFF/.test(gf), gf.slice(0, 80));
  await page.selectOption('#engineSel', 'MMFF94s');

  // ================= 4. structure change clears the handle =================
  await page.evaluate(() => { document.getElementById('input').value = 'CCC propane'; process(true); });
  await page.waitForTimeout(800);
  const gone = await page.evaluate(() => imtdVerifyHandle === null && lastLandscape === null);
  check('molecule swap drops the landscape + handles', gone);

  check('zero page errors', errors.length === 0, errors.join(' | '));

  console.log(`\nRESULT: ${pass} passed, ${fail} failed`);
  await browser.close();
  process.exit(fail ? 1 : 0);
})().catch(e => { console.error('FATAL:', e); process.exit(1); });
