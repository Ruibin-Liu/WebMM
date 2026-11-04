// m24 — design loop: neighborhood-explorer results push to the workbench
// (analogToBench: input.value + switchMode('single') + process(true))
// 用户点题功能;纯前端 DOM 接线,m22 golden 不涉 app DOM。两页(index 工作台壳 +
// platform 平台壳)各有内联副本,均验证。
const { chromium } = require('/opt/homebrew/lib/node_modules/playwright');
const EXE = '/Users/rliu/Library/Caches/ms-playwright/chromium_headless_shell-1243/chrome-headless-shell-mac-arm64/chrome-headless-shell';

let pass = 0, fail = 0;
function check(name, cond, detail) {
  if (cond) { pass++; console.log('  ✓ ' + name); }
  else { fail++; console.log('  ✗ ' + name + (detail ? ' — ' + detail : '')); }
}

async function openPage(browser, url) {
  const page = await (await browser.newContext()).newPage();
  const errors = [];
  page.on('pageerror', e => errors.push('PAGEERROR: ' + e.message.slice(0, 200)));
  await page.goto(url, { waitUntil: 'load' });
  await page.waitForFunction(() => window.webmm !== undefined, null, { timeout: 60000 });
  return { page, errors };
}

(async () => {
  const browser = await chromium.launch({ executablePath: EXE });

  for (const [label, url] of [
    ['workbench (index.html)', 'http://localhost:8901/app/index.html'],
    ['platform (platform.html)', 'http://localhost:8901/app/platform.html'],
  ]) {
    const { page, errors } = await openPage(browser, url);
    await page.waitForFunction(() => typeof analogToBench === 'function', null, { timeout: 60000 });

    // ---- wiring: mocked results, direct render, click the push button ----
    const res = await page.evaluate(async () => {
      const fake = [
        { idx: 'a1', smiles: 'CCO', mw: 46.1, tanimoto: 0.9, rotb: 0, sa: 1.0, flexT: null },
        { idx: 'a2', smiles: 'CCCO', mw: 60.1, tanimoto: 0.85, rotb: 1, sa: 1.0, flexT: null },
      ];
      renderAnalogRows(fake, analogRunToken);
      document.getElementById('analogResultsWrap').style.display = '';
      const btns = document.querySelectorAll('#analogRows span[onclick^="analogToBench"]');
      if (btns.length !== 2) return { ok: false, why: 'push buttons', n: btns.length };
      btns[1].click();
      await new Promise(r => setTimeout(r, 300));
      return {
        ok: true,
        input: document.getElementById('input').value,
        singleActive: document.getElementById('tabSingle').classList.contains('active'),
        status: document.getElementById('analogStatus').textContent,
      };
    });
    check(`[${label}] mocked render → ⇗ fills #input with the row SMILES`, res.ok && res.input === 'CCCO', JSON.stringify(res));
    check(`[${label}] mode switches to single`, res.singleActive);
    check(`[${label}] feedback in #analogStatus`, /analog #a2 → 工作台/.test(res.status || ''), res.status);
    check(`[${label}] zero page errors`, errors.length === 0, errors.join(' | '));
    await page.context().close();
  }

  // ---- E2E on the workbench: real exploration → push → single view renders ----
  const { page, errors } = await openPage(browser, 'http://localhost:8901/app/index.html');
  await page.evaluate(() => { document.getElementById('input').value = 'CCO ethanol'; process(true); });
  await page.waitForFunction(() => !!document.querySelector('#svg svg'), null, { timeout: 60000 });
  await page.evaluate(() => embed3D());
  await page.waitForFunction(() => !!document.querySelector('#viewer3d canvas'), null, { timeout: 60000 });
  await page.evaluate(() => { switchMode('analog'); analogDetectSites(); });
  await page.waitForFunction(() => [...document.getElementById('analogSite').options].some(o => o.value), null, { timeout: 30000 });
  await page.evaluate(() => {
    document.getElementById('analogSite').selectedIndex = 1;
    document.getElementById('analogFragCat').value = 'halo';
    runAnalogExplorer();
  });
  await page.waitForFunction(
    () => document.querySelectorAll('#analogRows span[onclick^="analogToBench"]').length > 0,
    null, { timeout: 90000 });
  const out = await page.evaluate(async () => {
    const btn = document.querySelector('#analogRows span[onclick^="analogToBench"]');
    const rowSmiles = btn.closest('tr').querySelector('td:nth-child(2)').textContent;
    btn.click();
    await new Promise(r => setTimeout(r, 2000));
    return {
      rowSmiles,
      input: document.getElementById('input').value,
      singleActive: document.getElementById('tabSingle').classList.contains('active'),
      svg: !!document.querySelector('#svg svg'),
    };
  });
  check('[E2E] real explorer row → #input matches the row SMILES', out.input === out.rowSmiles, JSON.stringify(out));
  check('[E2E] single view active and molecule rendered', out.singleActive && out.svg, JSON.stringify(out));
  check('[E2E] zero page errors', errors.length === 0, errors.join(' | '));

  console.log(`\nRESULT: ${pass} passed, ${fail} failed`);
  await browser.close();
  process.exit(fail ? 1 : 0);
})().catch(e => { console.error('FATAL', e.message); process.exit(1); });
