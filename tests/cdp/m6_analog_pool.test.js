// m6_analog_pool.test.js — the M3 persistent worker pool. Verifies the spec
// behaviors that matter: workers are REUSED across runs (cancel/finish no
// longer throw WASM heaps away), a poisoned worker is recycled and its job
// re-dispatched, and the pool stays healthy after poison chaos.
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

  // molecule + embed (site detection needs the parent 3D conformer) + sites
  await page.evaluate(() => { document.getElementById('input').value = 'CC(=O)Oc1ccccc1O'; process(true); });
  await page.waitForFunction(() => !!document.querySelector('#svg svg'), null, { timeout: 60000 });
  await page.evaluate(() => embed3D());
  await page.waitForFunction(() => !!document.querySelector('#viewer3d canvas'), null, { timeout: 60000 });
  await page.evaluate(() => { switchMode('analog'); analogDetectSites(); });
  await page.waitForFunction(() => [...document.getElementById('analogSite').options].some(o => o.value), null, { timeout: 30000 });

  const runOnce = () => page.evaluate(() => {
    document.getElementById('analogSite').selectedIndex = 1;
    document.getElementById('analogFragCat').value = 'alkyl';
    runAnalogExplorer();
  });
  const awaitDone = t => page.waitForFunction(() => /Done: \d+ scored/.test(document.getElementById('analogStatus').textContent), null, { timeout: t });

  // ---- 1/2: two consecutive runs — the second must reuse the pool ----
  await runOnce();
  await awaitDone(180000);
  const s1 = await page.evaluate(() => window.__analogPool.stats());
  check('run 1 completes with a warm pool', s1.size > 0 && s1.ready === s1.size, JSON.stringify(s1));
  await runOnce();
  await awaitDone(180000);
  const s2 = await page.evaluate(() => window.__analogPool.stats());
  check('run 2 reuses the pool (zero respawns)', s2.generation === s1.generation && s2.ready === s2.size,
    JSON.stringify({ gen1: s1.generation, gen2: s2.generation }));

  // ---- 3: poison — a bad-typed payload throws inside the worker ----
  const genBefore = s2.generation;
  const poisonRes = await page.evaluate(() =>
    window.__analogPool.submit({ type: 'score', idx: 'probe', molblock: 12345, nConfs: 10, seed: 1, parentSdf: 'x' }));
  const s3 = await page.evaluate(() => window.__analogPool.stats());
  check('poison job settles as failed (one re-dispatch)', !!poisonRes.failed, JSON.stringify(poisonRes).slice(0, 80));
  check('poisoned workers were recycled + respawned', s3.poisons >= 2 && s3.generation >= genBefore + 2, JSON.stringify(s3));
  // respawned workers need ~1-2s to re-import the engine — wait for readiness
  await page.waitForFunction(() => { const t = window.__analogPool.stats(); return t.ready === t.size; }, null, { timeout: 20000 })
    .catch(() => {});
  const s3b = await page.evaluate(() => window.__analogPool.stats());
  check('pool back at full size after poison chaos', s3b.ready === s3b.size, JSON.stringify(s3b));

  // ---- 4: the pool still works — a real run completes ----
  await runOnce();
  await awaitDone(180000);
  check('real run completes after poison recycling', true);
  const s4 = await page.evaluate(() => window.__analogPool.stats());
  check('heap reports are finite numbers', s4.heaps.length > 0 && s4.heaps.every(h => isFinite(h) && h >= 0), JSON.stringify(s4.heaps));
  // NOTE: heaps read 0 until the engine exports its wasm memory (the 'heap'
  // protocol + threshold recycling are wired and activate when it does).
  check('no page errors', errors.length === 0, errors.slice(0, 2).join(' | '));

  await browser.close();
  console.log('RESULT: ' + pass + ' passed, ' + fail + ' failed');
  process.exit(fail ? 1 : 0);
})().catch(e => { console.error(e); process.exit(1); });
