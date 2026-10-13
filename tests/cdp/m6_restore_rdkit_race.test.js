// m6_restore_rdkit_race.test.js — regression: restoring the search library
// while RDKit's WASM is still loading must not leave the library permanently
// empty. Deterministic: holds the RDKit_minimal.wasm response so restore is
// guaranteed to run first, then releases it and expects the library back.
// Fixture: /tmp/m5-project.json, produced by m6_platform.test.js. Run m6 first.
const { chromium } = require('/opt/homebrew/lib/node_modules/playwright');
const fs = require('fs');
const EXE = '/Users/rliu/Library/Caches/ms-playwright/chromium_headless_shell-1243/chrome-headless-shell-mac-arm64/chrome-headless-shell';
const URL = 'http://localhost:8901/app/platform.html';
const FIXTURE = '/tmp/m5-project.json';
const HOLD_MS = 3000;

let pass = 0, fail = 0;
function check(name, cond, detail) {
  if (cond) { pass++; console.log('  \u2713 ' + name); }
  else { fail++; console.log('  \u2717 ' + name + (detail ? ' \u2014 ' + detail : '')); }
}
const rdkitReady = page => page.evaluate(() => document.getElementById('rdkitVersion').textContent !== 'Loading...');
const libLen = page => page.evaluate(() => (window.__search?.getState() || []).length);

(async () => {
  if (!fs.existsSync(FIXTURE)) { console.log('missing ' + FIXTURE + ' \u2014 run tests/cdp/m6_platform.test.js first'); process.exit(2); }
  const browser = await chromium.launch({ executablePath: EXE });
  const page = await (await browser.newContext({ viewport: { width: 1440, height: 900 } })).newPage();
  const errors = [];
  page.on('pageerror', e => errors.push('PAGEERROR: ' + e.message.slice(0, 300)));
  page.on('console', m => { if (m.type() === 'error') errors.push(m.text().slice(0, 200)); });
  page.on('dialog', d => d.accept());

  // baseline: import the project normally
  await page.goto(URL, { waitUntil: 'load' });
  await page.waitForFunction(() => document.getElementById('rdkitVersion').textContent !== 'Loading...', null, { timeout: 60000 });
  await page.evaluate(() => switchMode('search'));
  await page.waitForFunction(() => window.__platformCompat, null, { timeout: 30000 });
  const nav = page.waitForNavigation({ waitUntil: 'load', timeout: 60000 });
  await page.setInputFiles('#projectImportFile', FIXTURE);
  await nav;
  await page.waitForFunction(() => document.getElementById('rdkitVersion').textContent !== 'Loading...', null, { timeout: 60000 });
  await page.evaluate(() => switchMode('search'));
  await page.waitForFunction(() => (window.__search?.getState() || []).length > 0, null, { timeout: 30000 }).catch(() => {});
  const expected = await libLen(page);
  check('baseline: imported library is non-empty', expected > 0, 'lib=' + expected);

  // race: hold the WASM (routing also bypasses the HTTP cache), reload into search
  let release; const held = new Promise(r => (release = r));
  let wasmHeld = false; let continued;
  await page.route('**/RDKit_minimal.wasm', route => {
    wasmHeld = true;
    continued = held.then(() => route.continue()).catch(e => { if (!/already handled/.test(e.message)) throw e; });
    return continued;
  });
  await page.reload({ waitUntil: 'domcontentloaded' });
  await page.waitForFunction(() => typeof switchMode === 'function', null, { timeout: 30000 });
  await page.evaluate(() => switchMode('search'));
  await page.waitForFunction(() => window.__platformCompat, null, { timeout: 30000 });
  await page.waitForTimeout(HOLD_MS);

  // guard against vacuous passes: the race must actually have been forced
  const duringLib = await libLen(page);
  check('race forced: wasm request intercepted', wasmHeld);
  check('race forced: RDKit not ready during hold', !(await rdkitReady(page)));
  check('race forced: library empty during hold', duringLib === 0, 'lib=' + duringLib);

  release();
  await continued;                                   // let the held request go through before unrouting
  await page.unroute('**/RDKit_minimal.wasm');
  await page.waitForFunction(() => document.getElementById('rdkitVersion').textContent !== 'Loading...', null, { timeout: 60000 });
  await page.waitForFunction(n => (window.__search?.getState() || []).length === n, expected, { timeout: 15000 }).catch(() => {});
  const got = await libLen(page);
  check('library restored after RDKit becomes ready', got === expected, 'lib=' + got + ' expected=' + expected);
  check('no page errors', errors.length === 0, errors.join(' | '));

  await browser.close();
  console.log('RESULT: ' + pass + ' passed, ' + fail + ' failed');
  process.exit(fail ? 1 : 0);
})().catch(e => { console.error(e); process.exit(1); });
