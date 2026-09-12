// M0 acceptance: app/ is behaviorally identical to MC for the RDKit chain.
// Drives BOTH pages with the same inputs and asserts per-field equality.
const { chromium } = require('/opt/homebrew/lib/node_modules/playwright');

const EXE = '/Users/rliu/Library/Caches/ms-playwright/chromium_headless_shell-1234/chrome-headless-shell-mac-arm64/chrome-headless-shell';
const INPUTS = [
  { label: 'caffeine (kekule)', smiles: 'CN1C=NC2=C1C(=O)N(C)C(=O)N2' },
  { label: 'ibuprofen', smiles: 'CC(C)Cc1ccc(cc1)C(C)C(=O)O' },
  { label: 'ascorbic acid', smiles: 'OCC(O)C1OC(=O)C(O)=C1O' },
  { label: 'caffeine (aromatic, known get_mol NULL — parity of failure)', smiles: 'Cn1cnc2c1c(=O)n(C)c(=O)n2' },
];

async function collect(page, smiles) {
  return page.evaluate((smi) => {
    document.getElementById('input').value = smi;
    try { process(true); } catch (e) {}
    const grab = (id) => {
      const el = document.getElementById(id);
      return el ? [...el.querySelectorAll('tr')].map(tr => [...tr.querySelectorAll('td')].map(td => td.textContent.trim())) : null;
    };
    return {
      smiles: document.getElementById('smiles')?.textContent ?? '',
      inchi: document.getElementById('inchi')?.textContent ?? '',
      inchikey: document.getElementById('inchikey')?.textContent ?? '',
      props: grab('props'), rules: grab('rules'),
      moreProps: grab('moreProps'), moreRules: grab('moreRules'),
      error: document.getElementById('error').textContent,
      svgAtoms: document.querySelectorAll('#svg svg [class*="atom"]').length,
      outputVisible: document.getElementById('output').style.display,
    };
  }, smiles);
}

(async () => {
  const browser = await chromium.launch({ executablePath: EXE });
  const mc = await browser.newPage({ viewport: { width: 1440, height: 900 } });
  const app = await browser.newPage({ viewport: { width: 1440, height: 900 } });
  const appErrors = [];
  app.on('console', m => { if (m.type() === 'error') appErrors.push(m.text()); });
  app.on('pageerror', e => appErrors.push('PAGEERROR: ' + e.message));

  await mc.goto('https://ruibin-liu.github.io/molecule-clipboard/', { waitUntil: 'load' });
  await app.goto('http://localhost:8901/app/index.html', { waitUntil: 'load' });
  await mc.waitForFunction(() => document.getElementById('rdkitVersion').textContent !== 'Loading...', null, { timeout: 30000 });
  await app.waitForFunction(() => document.getElementById('rdkitVersion').textContent !== 'Loading...', null, { timeout: 30000 });
  console.log('RDKit versions — MC:', await mc.evaluate(() => document.getElementById('rdkitVersion').textContent),
              '| app:', await app.evaluate(() => document.getElementById('rdkitVersion').textContent));

  let pass = 0, fail = 0;
  const check = (name, cond, detail = '') => {
    if (cond) { pass++; console.log(`  ✓ ${name}${detail ? ' — ' + detail : ''}`); }
    else { fail++; console.log(`  ✗ ${name}${detail ? ' — ' + detail : ''}`); }
  };

  for (const inp of INPUTS) {
    const a = await collect(mc, inp.smiles);
    const b = await collect(app, inp.smiles);
    console.log(`\n${inp.label}:`);
    check('identical parse outcome', (a.smiles !== '') === (b.smiles !== ''),
          a.smiles ? `MC=${a.smiles} app=${b.smiles}` : `both fail (MC err: ${a.error})`);
    check('props identical', JSON.stringify(a.props) === JSON.stringify(b.props), JSON.stringify(b.props));
    check('rules identical', JSON.stringify(a.rules) === JSON.stringify(b.rules), b.rules.map(r => r.join('=')).join(' '));
    check('moreProps identical', JSON.stringify(a.moreProps) === JSON.stringify(b.moreProps), JSON.stringify(b.moreProps));
    check('moreRules identical', JSON.stringify(a.moreRules) === JSON.stringify(b.moreRules));
    check('InChIKey identical', a.inchikey === b.inchikey, b.inchikey);
    if (a.smiles) check('SVG rendered', b.svgAtoms > 0, `${b.svgAtoms} atoms`);
  }

  // ---- app-only UI checks ----
  console.log('\nUI checks (app):');
  const ui = await app.evaluate(() => {
    const strongs = [...document.querySelectorAll('#output strong')].map(s => s.textContent.trim());
    const btnTexts = [...document.querySelectorAll('#output button')].map(b => b.textContent.trim());
    return {
      strongs,
      has2D: strongs.includes('2D Structure'),
      has3D: strongs.includes('3D Structure'),
      hasProps: strongs.includes('Properties'),
      btn2D: ['SVG', 'PNG', 'Copy', 'Link'].every(t => btnTexts.includes(t)),
      spatial: [...document.querySelectorAll('.action-btn.spatial')].map(b => b.textContent.trim()),
      export3d: [...document.querySelectorAll('.action-btn.export')].map(b => b.textContent.trim()),
      ready3d: (() => {
        // Embed/Conformers unlock on parse; Optimize/exports unlock after a
        // 3D structure exists (reset3DState + enable3dActions/showSdfIn3D).
        const byLabel = (t) => [...document.querySelectorAll('.panel .action-btn')]
          .find(b => b.textContent.trim().startsWith(t));
        const embed = byLabel('Embed 3D'), conf = byLabel('Conformers');
        return embed && conf && !embed.disabled && !conf.disabled;
      })(),
    };
  });
  check('panel titles: 2D Structure / 3D Structure / Properties', ui.has2D && ui.has3D && ui.hasProps, ui.strongs.join(' / '));
  check('2D control buttons', ui.btn2D);
  const expSet = new Set(ui.export3d);
  const hasExports = ['Export SDF', 'XYZ', 'PNG', 'Export CSV'].every(t => expSet.has(t));
  check('3D buttons present + Embed/Conformers enabled on parse', ui.spatial.length === 3 && hasExports && ui.ready3d,
        [...ui.spatial, ...ui.export3d].join(','));

  // URL state + error parity
  await app.evaluate(() => process(true));
  check('URL state encoded', (await app.evaluate(() => location.hash)).startsWith('#mol='));
  const mcErr = await collect(mc, 'not_a_molecule!!');
  const appErr = await collect(app, 'not_a_molecule!!');
  check('bad input error identical', mcErr.error === appErr.error && appErr.error.length > 0, appErr.error);

  // JSME vendored
  await app.evaluate(() => openJSME());
  await app.waitForFunction(() => typeof JSApplet !== 'undefined', null, { timeout: 20000 }).catch(() => {});
  await app.waitForTimeout(1200);
  check('JSME loads (vendored)', await app.evaluate(() => typeof JSApplet !== 'undefined'));

  // History (3 valid molecules processed above)
  const hist = await app.evaluate(() => { loadHistory(); return compoundHistory.length; });
  check('history persisted', hist >= 3, hist + ' entries');

  // Screenshots at two breakpoints
  await app.evaluate(s => { document.getElementById('input').value = s; process(false); }, 'CN1C=NC2=C1C(=O)N(C)C(=O)N2');
  await app.waitForTimeout(300);
  await app.screenshot({ path: '/tmp/app_1440.png', fullPage: true });
  await app.setViewportSize({ width: 820, height: 900 });
  await app.waitForTimeout(300);
  await app.screenshot({ path: '/tmp/app_820.png', fullPage: true });
  check('screenshots taken', true);

  check('zero console errors (app)', appErrors.length === 0, appErrors.slice(0, 3).join(' | '));

  console.log(`\nRESULT: ${pass} passed, ${fail} failed`);
  await browser.close();
  process.exit(fail ? 1 : 0);
})();
