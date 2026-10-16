// m6_sa_warmup.test.js — the SA-table idle warmup must survive the engine WASM
// not being ready yet. Old behavior: idle fires while window.webmm is still
// undefined → fetch resolves, window.webmm.sa_load_table_wasm(...) throws,
// warning logged, warmup wasted (first analog run then eats the ~1.4s sync
// parse). Fix: warmup polls until webmm exists. Deterministic: the engine
// script is held so webmm stays undefined across the idle window.
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
  const warnings = [], errors = [];
  page.on('console', m => {
    if (/\[analog\] SA table fetch failed/.test(m.text())) warnings.push(m.text());
    if (m.type() === 'error') errors.push(m.text().slice(0, 150));
  });
  page.on('pageerror', e => errors.push('PAGEERROR: ' + e.message.slice(0, 200)));

  // hold the engine script(s) so window.webmm stays undefined past idle
  let release; const held = new Promise(r => (release = r)); let continued;
  const hold = route => {
    continued = held.then(() => route.continue()).catch(e => { if (!/already handled/.test(e.message)) throw e; });
    return continued;
  };
  await page.route('**/pkg/webmm.js', hold);
  await page.route('**/webmm.js', hold);

  await page.goto(URL, { waitUntil: 'load' });
  await page.waitForTimeout(1800);          // let idle fire while webmm is undefined
  const webmmDuringHold = await page.evaluate(() => window.webmm === undefined);
  check('race forced: webmm undefined during idle window', webmmDuringHold);
  check('no wasted-warmup warning while engine loads', warnings.length === 0, warnings[0]);

  release(); await continued;
  await page.unroute('**/pkg/webmm.js'); await page.unroute('**/webmm.js');
  await page.waitForFunction(() => window.webmm !== undefined, null, { timeout: 60000 });
  await page.waitForFunction(() => document.getElementById('rdkitVersion').textContent !== 'Loading...', null, { timeout: 60000 });
  // the warmup retry (500ms cadence) should load the table on its own — prove
  // it by scoring directly; sa_score_wasm requires the loaded table.
  let sa = null, saErr = null;
  for (let i = 0; i < 20 && sa == null; i++) {
    await page.waitForTimeout(500);
    ({ sa, saErr } = await page.evaluate(() => {
      try {
        const mb = rdkitModule.get_mol('CC(=O)Oc1ccccc1C(=O)O').get_molblock();
        const v = window.webmm.sa_score_wasm(mb);
        return { sa: isFinite(v) ? v : null };
      } catch (e) { return { sa: null, saErr: String(e).slice(0, 120) }; }
    }));
  }
  check('warmup loaded the SA table by itself (direct score works)', sa !== null && sa > 0, 'sa=' + sa + (saErr ? ' err=' + saErr : ''));
  check('no page errors', errors.length === 0, errors.join(' | '));

  await browser.close();
  console.log('RESULT: ' + pass + ' passed, ' + fail + ' failed');
  process.exit(fail ? 1 : 0);
})().catch(e => { console.error(e); process.exit(1); });
