// m6_import_loop.test.js — stress the import→reload race (snapshot guard +
// projection invalidation). Repeats import + immediate reload N times and
// checks the reloaded state equals the imported project every time.
// Fixture: /tmp/m5-project.json, produced by m6_platform.test.js. Run m6 first.
const { chromium } = require('/opt/homebrew/lib/node_modules/playwright');
const fs = require('fs');
const EXE = '/Users/rliu/Library/Caches/ms-playwright/chromium_headless_shell-1243/chrome-headless-shell-mac-arm64/chrome-headless-shell';
const URL = 'http://localhost:8901/app/platform.html';
const FIXTURE = '/tmp/m5-project.json';
const N = Number(process.env.LOOP_N || 10);

let pass = 0, fail = 0;
function check(name, cond, detail) {
  if (cond) { pass++; console.log('  ✓ ' + name); }
  else { fail++; console.log('  ✗ ' + name + (detail ? ' — ' + detail : '')); }
}

async function waitReady(page) {
  await page.waitForFunction(() => document.getElementById('rdkitVersion').textContent !== 'Loading...', null, { timeout: 60000 });
  await page.waitForFunction(() => window.webmm !== undefined, null, { timeout: 60000 });
  await page.evaluate(() => switchMode('search'));             // __platformCompat is created on entering search
  await page.waitForFunction(() => window.__platformCompat, null, { timeout: 30000 });
}

(async () => {
  if (!fs.existsSync(FIXTURE)) { console.log('missing ' + FIXTURE + ' — run tests/cdp/m6_platform.test.js first'); process.exit(2); }
  const browser = await chromium.launch({ executablePath: EXE });
  const ctx = await browser.newContext({ viewport: { width: 1440, height: 900 } });
  const page = await ctx.newPage();
  const errors = []; let skipped = 0;
  const logs = [];
  page.on('console', m => {
    const t = m.text(); logs.push(m.type() + ': ' + t.slice(0, 200));
    if (/snapshot skipped/.test(t)) skipped++;
    if (m.type() === 'error') errors.push(t.slice(0, 200));
  });
  page.on('pageerror', e => errors.push('PAGEERROR: ' + e.message.slice(0, 300)));
  page.on('dialog', d => d.accept());

  await page.goto(URL, { waitUntil: 'load' });
  await waitReady(page);

  console.log('import/reload loop x' + N + ':');
  for (let i = 1; i <= N; i++) {
    const errBefore = errors.length; logs.length = 0;
    let st;
    try {
      await page.evaluate(() => switchMode('search'));
      await page.evaluate(() => clearSearchLibrary());           // sentinel: lib -> 0
      await page.waitForTimeout(200);
      const nav = page.waitForNavigation({ waitUntil: 'load', timeout: 60000 });
      await page.setInputFiles('#projectImportFile', FIXTURE);
      await nav;
      await page.reload({ waitUntil: 'load' });                  // race: no lock flush
      await waitReady(page);
      await page.evaluate(() => switchMode('search'));
      // lib restores first, the pin-count render (renderWorkSet) lands a beat
      // later — wait for BOTH before sampling (r8 flake)
      await page.waitForFunction(() => (window.__search?.getState() || []).length > 0
        && /已钉/.test((document.getElementById('workSetCounts') || { textContent: '' }).textContent), null, { timeout: 30000 }).catch(() => {});
      st = await page.evaluate(() => ({
        lib: (window.__search.getState() || []).length,
        hist: document.querySelectorAll('#queryHistoryPanel div[onclick]').length,
        pins: document.getElementById('workSetCounts').textContent,
      }));
    } catch (e) { st = { threw: String(e).slice(0, 200) }; }
    const newErr = errors.slice(errBefore);
    if (!(st.lib === 55) && process.env.STOP_ON_FAIL) {
      await page.waitForTimeout(15000);
      const diag = await page.evaluate(() => {
        const P = window.__platformCompat; const api = P?.projectApi?.(); const S = api?.state || {};
        return {
          libAfter15s: (window.__search?.getState() || []).length,
          ready: !!P?.ready?.(), restored: !!window.__search?.restored,
          stateKeys: Object.keys(S),
          rounds: Object.keys(S.rounds || {}).length,
          libraries: Object.fromEntries(Object.entries(S.libraries || {}).map(([k, v]) => [k, Array.isArray(v) ? v.length : (v && typeof v === 'object' ? Object.keys(v).length + ' keys: ' + Object.keys(v).slice(0, 6).join(',') : v)])),
          pins: Array.isArray(S.pins) ? S.pins.length : Object.keys(S.pins || {}).length,
          maxCommandId: S.maxCommandId,
          mainEntries: S.libraries?.main?.entries?.size ?? null,
        };
      }).catch(e => ({ diagFailed: String(e) }));
      console.log('DIAG #' + i + ' ' + JSON.stringify(diag));
      console.log('LOGS #' + i + '\n  ' + logs.join('\n  '));
      check('#' + i + ' (diag)', false, JSON.stringify(st)); break;
    }
    check('#' + i + ' reloaded state == imported project, no errors',
      st.lib === 55 && st.hist >= 1 && /1 已钉/.test(st.pins || '') && newErr.length === 0,
      JSON.stringify(st) + (newErr.length ? ' errors=' + JSON.stringify(newErr) : ''));
  }
  console.log('snapshot skipped (guard fired, allowed): ' + skipped);
  console.log('RESULT: ' + pass + ' passed, ' + fail + ' failed');
  await browser.close();
  process.exit(fail ? 1 : 0);
})().catch(e => { console.log('ERR ' + (e.stack || e)); process.exit(1); });
