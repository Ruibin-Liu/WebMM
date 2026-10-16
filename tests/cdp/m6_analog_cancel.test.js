// m6_analog_cancel.test.js — analog worker lifecycle. Two checks:
// 1) engine-fetch failure inside the workers must fail fast with a status
//    message (old code: ensure() polled forever, run hung, intervals leaked);
// 2) re-running while a run is active must not break the second run (stale
//    token now terminates the OLD farm — the observable contract is: second
//    run completes with Done and the page stays error-free).
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
  const page = await (await browser.newContext()).newPage();
  const errors = [];
  page.on('pageerror', e => errors.push('PAGEERROR: ' + e.message.slice(0, 200)));
  page.on('console', m => { if (m.type() === 'error') errors.push(m.text().slice(0, 150)); });

  await page.goto(URL, { waitUntil: 'load' });
  await page.waitForFunction(() => window.webmm !== undefined, null, { timeout: 60000 });

  // ---- 1) worker-side engine fetch failure -> fast, visible error ----
  await page.route('**/pkg/webmm.js', r => r.abort());   // main thread already loaded; workers' import fails
  await page.route('**/pkg/webmm_bg.wasm', r => r.abort());
  await page.evaluate(() => { document.getElementById('input').value = 'CC(=O)Oc1ccccc1O'; process(true); });
  await page.waitForFunction(() => !!document.querySelector('#svg svg'), null, { timeout: 60000 });
  await page.evaluate(() => embed3D());   // site detection needs the parent 3D conformer
  await page.waitForFunction(() => !!document.querySelector('#viewer3d canvas'), null, { timeout: 60000 });
  await page.evaluate(() => { switchMode('analog'); analogDetectSites(); });   // the 用当前分子 button handler
  await page.waitForFunction(() => [...document.getElementById('analogSite').options].some(o => o.value), null, { timeout: 30000 });
  await page.evaluate(() => {
    document.getElementById('analogSite').selectedIndex = 1;
    document.getElementById('analogFragCat').value = 'halo';
    runAnalogExplorer();
  });
  let failedFast = false;
  try {
    await page.waitForFunction(() => /Worker init failed/.test(document.getElementById('analogStatus').textContent), null, { timeout: 15000 });
    failedFast = true;
  } catch (e) { /* keep false */ }
  check('engine-fetch failure fails fast with a status message', failedFast,
    await page.evaluate(() => document.getElementById('analogStatus').textContent.slice(0, 60)));
  await page.unroute('**/pkg/webmm.js'); await page.unroute('**/pkg/webmm_bg.wasm');

  // ---- 2) re-run while active: second run completes, page healthy ----
  await page.evaluate(() => {
    document.getElementById('analogSite').selectedIndex = 1;
    runAnalogExplorer();          // run #1
    runAnalogExplorer();          // immediately superseded -> stale farm
  });
  await page.waitForFunction(() => /Done: \d+ scored/.test(document.getElementById('analogStatus').textContent), null, { timeout: 180000 });
  check('second run completes after superseding the first', true);
  check('no page errors across both scenarios', errors.length === 0, errors.slice(0, 2).join(' | '));

  await browser.close();
  console.log('RESULT: ' + pass + ' passed, ' + fail + ' failed');
  process.exit(fail ? 1 : 0);
})().catch(e => { console.error(e); process.exit(1); });
