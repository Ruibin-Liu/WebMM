// M3 regression: URL deep-links, history modal, axe a11y, full-chain smoke,
// conformer chart keyboard navigation. Run: node m3_site_nav.test.js
const { chromium } = require('/opt/homebrew/lib/node_modules/playwright');
const AXE = '/opt/homebrew/lib/node_modules/axe-core/axe.min.js';
const EXE = '/Users/rliu/Library/Caches/ms-playwright/chromium_headless_shell-1234/chrome-headless-shell-mac-arm64/chrome-headless-shell';

(async () => {
  const browser = await chromium.launch({ executablePath: EXE });
  const page = await browser.newPage({ viewport: { width: 1440, height: 900 } });
  const errors = [];
  page.on('pageerror', e => errors.push('PAGEERROR: ' + e.message.slice(0, 120)));
  page.on('console', m => { if (m.type() === 'error') errors.push(m.text().slice(0, 120)); });

  let pass = 0, fail = 0;
  const check = (name, cond, detail = '') => {
    if (cond) { pass++; console.log(`  ✓ ${name}${detail ? ' — ' + detail : ''}`); }
    else { fail++; console.log(`  ✗ ${name}${detail ? ' — ' + detail : ''}`); }
  };

  const hash = Buffer.from(JSON.stringify({ smiles: 'CN1C=NC2=C1C(=O)N(C)C(=O)N2', name: 'caffeine', engine: 'GFNFF', confN: '7', confRmsd: '0.6' })).toString('base64');
  await page.goto('http://localhost:8901/app/index.html#mol=' + hash, { waitUntil: 'load' });
  await page.waitForFunction(() => window.webmmReady === true, null, { timeout: 60000 });
  await page.waitForTimeout(500);

  // 1. URL restores molecule + engine + conformer params
  const restored = await page.evaluate(() => ({
    smiles: document.getElementById('input').value,
    engine: document.getElementById('engineSel').value,
    confN: document.getElementById('confN').value,
    confRmsd: document.getElementById('confRmsd').value,
  }));
  check('URL restores molecule', restored.smiles.includes('CN1C=NC2') || restored.smiles.includes('Cn1c'), restored.smiles.slice(0, 40));
  check('URL restores engine/confN/confRmsd', restored.engine === 'GFNFF' && restored.confN === '7' && restored.confRmsd === '0.6',
        JSON.stringify({ engine: restored.engine, confN: restored.confN, rmsd: restored.confRmsd }));

  // 2. setURL round-trip (share link copies the current state)
  const urlState = await page.evaluate(() => {
    setURL(document.getElementById('input').value.split(/\s+/)[0], 'caffeine');
    return JSON.parse(decodeURIComponent(escape(atob(location.hash.slice(5)))));
  });
  check('setURL round-trip includes engine/confN', urlState.engine === 'GFNFF' && urlState.confN === '7', JSON.stringify(urlState).slice(0, 80));

  // 3. history modal opens / Escape closes
  await page.evaluate(() => showHistory());
  await page.waitForTimeout(200);
  const modalOpen = await page.evaluate(() => {
    for (const id of ['historyModal', 'modal', 'history']) {
      const m = document.getElementById(id);
      if (m && m.style.display === 'block') return id;
    }
    return [...document.querySelectorAll('.modal')].find(m => m.style.display === 'block')?.id || null;
  });
  check('history modal opens', !!modalOpen, modalOpen || 'not found');
  await page.keyboard.press('Escape');
  await page.waitForTimeout(200);
  const modalClosed = await page.evaluate(() =>
    ![...document.querySelectorAll('.modal')].some(m => m.style.display === 'block'));
  check('Escape closes modal', modalClosed);

  // 4. axe: no serious/critical violations
  await page.addScriptTag({ path: AXE });
  const axe = await page.evaluate(() => axe.run(document, { resultTypes: ['violations'] })).then(r =>
    r.violations.filter(v => ['serious', 'critical'].includes(v.impact)));
  check('axe: no serious/critical violations', axe.length === 0, axe.map(v => v.id).join(',') || 'clean');

  // 5. full-chain smoke: embed -> optimize -> conformers
  await page.evaluate(() => embed3D());
  await page.waitForFunction(() => (sdf3d || '').includes('M  END'), null, { timeout: 120000 });
  await page.selectOption('#engineSel', 'MMFF94s');
  await page.evaluate(() => optimize3D());
  await page.waitForFunction(() => document.getElementById('status3d').textContent.includes('E(MMFF94s)'), null, { timeout: 180000 });
  await page.evaluate(() => { document.getElementById('confN').value = '6'; runConformers(); });
  await page.waitForFunction(() => !confRun && !reoptState && !document.getElementById('btnConf').disabled, null, { timeout: 300000 });
  const chain = await page.evaluate(() => document.getElementById('status3d').textContent);
  check('full-chain smoke (embed→opt→conformers→reopt)', chain.includes('kept') || chain.includes('Re-optimization done'), chain.slice(0, 80));

  // 6. conformer chart keyboard navigation
  const kb = await page.evaluate(() => {
    const cv = document.getElementById('confChart');
    cv.focus();
    cv.dispatchEvent(new KeyboardEvent('keydown', { key: 'ArrowRight', bubbles: true }));
    return { focused: document.activeElement === cv, status: document.getElementById('status3d').textContent };
  });
  check('conformer chart keyboard navigation (→)', kb.focused && kb.status.includes('Conformer #2'), kb.status.slice(0, 50));

  // 7. topnav present with the unified links
  const nav = await page.evaluate(() => [...document.querySelectorAll('.topnav .links a')].map(a => a.textContent.trim()));
  check('unified topnav (Workbench/Demo/Playground/GitHub)',
        ['Workbench', 'Demo', 'Playground', 'GitHub'].every(t => nav.includes(t)), nav.join(','));

  // 8. zero page errors
  check('zero page errors', errors.length === 0, errors.slice(0, 2).join(' | '));

  console.log(`\nRESULT: ${pass} passed, ${fail} failed`);
  await browser.close();
  process.exit(fail ? 1 : 0);
})();
