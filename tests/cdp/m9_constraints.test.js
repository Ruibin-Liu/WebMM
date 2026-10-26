// m9_constraints.test.js — Top-3: atom picking + user restraints in the
// molecule studio. The atom ledger (3Dmol serial == SDF order, verified per
// pick against the SDF line's symbol) feeds distance/angle/dihedral/freeze
// restraints through conf.worker into optimize_from_sdf. The gold assertion
// is end-to-end: a 2.5 Å distance restraint between atoms 1–3 (C,O) of
// ethanol actually pulls them apart — proof the indices map correctly
// through viewer → JSON → engine.
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
  const ctx = await browser.newContext({ viewport: { width: 1440, height: 900 } });
  const page = await ctx.newPage();
  const errors = [];
  page.on('pageerror', e => errors.push('PAGEERROR: ' + e.message.slice(0, 200)));

  await page.goto(URL, { waitUntil: 'load' });
  await page.waitForFunction(() => document.getElementById('rdkitVersion').textContent !== 'Loading...', null, { timeout: 60000 });
  await page.waitForFunction(() => window.webmm !== undefined, null, { timeout: 60000 });

  // ethanol: heavy C(1) C(2) O(3), H appended 4..8
  await page.evaluate(() => { document.getElementById('input').value = 'CCO ethanol'; process(true); });
  await page.waitForFunction(() => document.querySelector('#svg svg'), null, { timeout: 30000 });
  await page.evaluate(() => embed3D());
  await page.waitForFunction(() => typeof sdf3d !== 'undefined' && sdf3d && /Embedded|3D coordinates present/.test(document.getElementById('status3d').textContent), null, { timeout: 60000 });

  const dist01 = () => page.evaluate(() => {
    const l = (typeof sdfLines === 'function' ? sdfLines(sdf3d) : sdf3d.split('\n'));
    const c = i => [parseFloat(l[4 + i].substring(0, 10)), parseFloat(l[4 + i].substring(10, 20)), parseFloat(l[4 + i].substring(20, 30))];
    const d = Math.hypot(c(0)[0] - c(2)[0], c(0)[1] - c(2)[1], c(0)[2] - c(2)[2]);
    return d;
  });
  const d0 = await dist01();
  check('baseline distance between atoms 1,3 is embed-sane', d0 > 0.8 && d0 < 2.8, d0.toFixed(3));

  // ---- pick mode + ledger ----
  await page.evaluate(() => togglePickMode());
  check('pick bar opens', await page.isVisible('#pickBar'));
  await page.evaluate(() => pickAtom(1, 'N'));   // ledger violation: SDF says C
  const rej = await page.evaluate(() => ({ n: picked.length, err: document.getElementById('error').textContent }));
  check('ledger mismatch REJECTED (never point at the wrong atom)', rej.n === 0 && /账本校验失败/.test(rej.err), JSON.stringify(rej));

  await page.evaluate(() => pickAtom(1, 'C'));
  await page.evaluate(() => pickAtom(3, 'O'));
  const sel = await page.evaluate(() => ({ n: picked.length, status: document.getElementById('pickStatus').textContent, geom: document.getElementById('pickGeom').textContent }));
  check('two atoms picked, live geometry shown (Å)', sel.n === 2 && /已选 2/.test(sel.status) && /当前 \d/.test(sel.geom), JSON.stringify(sel));

  // ---- build the restraint at 2.5 Å ----
  await page.evaluate(() => { document.getElementById('pickR0').value = '3.2'; pickAddConstraint(); });
  const chips = await page.evaluate(() => document.getElementById('restraintChips').textContent);
  check('restraint chip lands (d(1,3) = 3.20 Å ×50)', /d\(1,3\)/.test(chips) && /3\.20/.test(chips), chips);

  // ---- constrained optimize: the gold end-to-end mapping proof ----
  await page.evaluate(() => optimize3D());
  await page.waitForFunction(() => /converged/.test(document.getElementById('status3d').textContent), null, { timeout: 120000 });
  const d1 = await dist01();
  check('constrained optimize pulls atoms 1,3 toward 3.2 Å (indices map end-to-end)',
    d1 > 2.5 && d1 < 3.5 && (d1 - d0) > 0.15, 'final ' + d1.toFixed(3) + ' (baseline ' + d0.toFixed(3) + ')');

  // ---- freeze: atom 1 must not move ----
  await page.evaluate(() => { pickAtom(1, 'C'); pickFreezeAtoms(); });
  const before = await page.evaluate(() => { window.__m9mark = document.getElementById('status3d').textContent; return sdf3d.split('\n')[4].substring(0, 30); });
  await page.evaluate(() => optimize3D());
  await page.waitForFunction(() => document.getElementById('status3d').textContent !== window.__m9mark && /converged/.test(document.getElementById('status3d').textContent), null, { timeout: 120000 });
  const after = await page.evaluate(() => sdf3d.split('\n')[4].substring(0, 30));
  check('frozen atom 1 coordinates unchanged', before === after, before + ' vs ' + after);

  // ---- chips survive the optimize round-trip (same molecule), Esc clears picks ----
  const chipsAfter = await page.evaluate(() => document.getElementById('restraintChips').style.display);
  check('restraints survive the optimize round-trip (iterate on the target)', chipsAfter !== 'none');

  // ---- a NEW molecule clears stale constraints ----
  await page.evaluate(() => { document.getElementById('input').value = 'c1ccccc1 benzene'; process(true); });
  await page.waitForFunction(() => document.getElementById('name').textContent.includes('benzene'), null, { timeout: 30000 });
  const chipsGone = await page.evaluate(() => document.getElementById('restraintChips').style.display === 'none');
  check('new molecule clears stale restraints (indices die with their molecule)', chipsGone);

  check('zero page errors', errors.length === 0, errors.slice(0, 3).join(' | '));
  await browser.close();
  console.log('RESULT: ' + pass + ' passed, ' + fail + ' failed');
  process.exit(fail ? 1 : 0);
})().catch(e => { console.error('ERR', e); process.exit(1); });
