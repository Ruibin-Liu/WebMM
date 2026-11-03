// m18_landscape_v2.test.js — landscape v2: plain-MD occupancy verification
// (unbiased ladder via live-coordinator verify batches) and the force-directed
// network view, side by side with the disconnectivity graph.
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

  // ================= 1. network view toggle =================
  const before = await page.evaluate(() => ({
    view: typeof landscapeView !== 'undefined' ? landscapeView : '?',
    leaves: document.querySelectorAll('#dgraphSvg circle[data-basin]').length,
    n: confEnsemble.length,
  }));
  check('default view is the disconnectivity graph', before.view === 'dg', before.view);
  await page.evaluate(() => setLandscapeView('net'));
  const net = await page.evaluate(() => ({
    nodes: document.querySelectorAll('#dgraphSvg circle[data-basin]').length,
    edges: document.querySelectorAll('#dgraphSvg line[data-edge]').length,
    legend: (document.querySelector('#dgraphSvg text[fill="#64748b"]') || {}).textContent || '',
  }));
  check('network view renders one node per basin', net.nodes === before.n, `nodes=${net.nodes} n=${before.n}`);
  check('network edges = transition pairs', net.edges >= 0 && /metad/.test(net.legend), `edges=${net.edges} · ${net.legend}`);
  const netLink = await page.evaluate(() => {
    const node = document.querySelector('#dgraphSvg circle[data-basin="1"]');
    if (!node) return { ok: false };
    node.dispatchEvent(new MouseEvent('click', { bubbles: true }));
    return { ok: true, sel: selConf };
  });
  check('network node click selects conformer #2', netLink.ok && netLink.sel === 1, JSON.stringify(netLink));
  await page.evaluate(() => setLandscapeView('dg'));
  const backDg = await page.evaluate(() => landscapeView === 'dg' &&
    document.querySelectorAll('#dgraphSvg circle[data-basin]').length === confEnsemble.length);
  check('toggle back to disconnectivity graph', backDg);

  // ================= 2. MD occupancy verification =================
  const noMdBefore = await page.evaluate(() => {
    const th = [...document.querySelectorAll('#basinTableWrap th')];
    return { mdCol: th.some(t => /MD 占据/.test(t.textContent)), md: lastLandscape && lastLandscape.md };
  });
  check('no MD column before verification', !noMdBefore.mdCol && !(noMdBefore.md && noMdBefore.md.frames));
  await page.evaluate(() => window.runMdVerify());
  await page.waitForFunction(() => /MD 占据验证完成/.test(document.getElementById('status3d').textContent), null, { timeout: 180000 });
  const after = await page.evaluate(() => {
    const th = [...document.querySelectorAll('#basinTableWrap th')];
    const rows = [...document.querySelectorAll('#basinTableWrap tbody tr')];
    const md = lastLandscape && lastLandscape.md;
    const mdOcc = md ? md.basins.reduce((a, b) => a + b.visits, 0) / Math.max(1, md.frames) : 0;
    return {
      mdCol: th.some(t => /MD 占据/.test(t.textContent)),
      mdFrames: md ? md.frames : 0,
      mdOccSum: mdOcc,
      meta: document.getElementById('landscapeMeta').textContent,
      tokens: (() => { const s = taskScheduler.stats(); return s.interactiveHeld === 0; })(),
      btnLabel: document.getElementById('btnMdVerify').textContent,
    };
  });
  check('MD column appears in the basin table', after.mdCol);
  check('MD ladder folded (frames > 0, occupancy ≤ 100%)', after.mdFrames > 0 && after.mdOccSum <= 1.0 + 1e-9,
        `frames=${after.mdFrames} occ=${(after.mdOccSum * 100).toFixed(0)}%`);
  check('meta line reports the MD frame count', /MD \d+ 帧无偏占据/.test(after.meta), after.meta);
  check('verify tokens fully released, button restored', after.tokens && after.btnLabel === 'MD 占据验证');

  // network now prefers the MD ladder
  await page.evaluate(() => setLandscapeView('net'));
  const netMd = await page.evaluate(() => {
    const t = (document.querySelector('#dgraphSvg text[fill="#64748b"]') || {}).textContent || '';
    return t;
  });
  check('network view switched to the MD ladder', /MD 无偏/.test(netMd), netMd);
  await page.evaluate(() => setLandscapeView('dg'));

  // a second batch is plannable (batches repeat)
  await page.evaluate(() => window.runMdVerify());
  await page.waitForFunction(() => /MD 占据验证完成/.test(document.getElementById('status3d').textContent), null, { timeout: 180000 });
  const md2 = await page.evaluate(() => lastLandscape && lastLandscape.md ? lastLandscape.md.frames : 0);
  check('second verification batch folds more frames', md2 > after.mdFrames, `${md2} > ${after.mdFrames}`);

  // ================= 3. structure change invalidates the handle =================
  await page.evaluate(() => { document.getElementById('input').value = 'CCC propane'; process(true); });
  await page.waitForTimeout(800);
  const gone = await page.evaluate(() => imtdVerifyHandle === null && lastLandscape === null);
  check('molecule swap drops the verify handle + landscape', gone);

  check('zero page errors', errors.length === 0, errors.join(' | '));

  console.log(`\nRESULT: ${pass} passed, ${fail} failed`);
  await browser.close();
  process.exit(fail ? 1 : 0);
})().catch(e => { console.error('FATAL:', e); process.exit(1); });
