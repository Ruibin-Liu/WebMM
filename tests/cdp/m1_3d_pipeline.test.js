// M1 acceptance v3.
//  (1) SMILES flow (13 heavy atoms): embed → MMFF94s UI optimize; deterministic
//      path identity vs optimize_from_sdf(2D molblock), checked at full precision
//      in-page. GFN-FF UI optimize + panel rendering.
//  (2) Demo-preset 24-atom caffeine SDF: MMFF94s must reproduce demo E = -123.49;
//      GFN-FF optimized geometry → xtb single-point cross-check (|ΔE| < 0.02 Eh).
//  (3) Exports: SDF readback via RDKit.js, XYZ, PNG. Zero console errors.
const { chromium } = require('/opt/homebrew/lib/node_modules/playwright');
const { execSync } = require('child_process');
const fs = require('fs');

const EXE = '/Users/rliu/Library/Caches/ms-playwright/chromium_headless_shell-1234/chrome-headless-shell-mac-arm64/chrome-headless-shell';
const XTB = process.env.HOME + '/.local/xtb-gxtb/bin/xtb';
const CAFF_SDF3D = fs.readFileSync('/tmp/caff24.sdf', 'utf8');

async function optimizeWith(page, engine, sdf) {
  return page.evaluate(([eng, sdfText]) => {
    const o = new webmm.OptimizationOptions();
    o.engine = eng;
    o.set_max_iterations(500);
    const res = webmm.optimize_from_sdf(sdfText, o);
    const terms = JSON.parse(res.get_energy_terms_json());
    const sum = Object.values(terms).reduce((a, b) => a + b, 0);
    return { E: res.final_energy, conv: res.get_converged(), iters: res.get_iterations(),
             terms, sum, engine: res.get_engine() };
  }, [engine, sdf]);
}

(async () => {
  const browser = await chromium.launch({ executablePath: EXE });
  const page = await browser.newPage({ viewport: { width: 1440, height: 900 } });
  const errors = [];
  page.on('console', m => { if (m.type() === 'error') errors.push(m.text()); });
  page.on('pageerror', e => errors.push('PAGEERROR: ' + e.message));

  const hash = Buffer.from(JSON.stringify({ smiles: 'CN1C=NC2=C1C(=O)N(C)C(=O)N2', name: 'caffeine' })).toString('base64');
  await page.goto('http://localhost:8901/app/index.html#mol=' + hash, { waitUntil: 'load' });
  await page.waitForFunction(() => window.webmmReady === true, null, { timeout: 60000 });

  let pass = 0, fail = 0;
  const check = (name, cond, detail = '') => {
    if (cond) { pass++; console.log(`  ✓ ${name}${detail ? ' — ' + detail : ''}`); }
    else { fail++; console.log(`  ✗ ${name}${detail ? ' — ' + detail : ''}`); }
  };

  // ---- (1) SMILES flow, heavy-atom caffeine ----
  await page.evaluate(() => embed3D());
  await page.waitForFunction(() => (sdf3d || '').includes('M  END'), null, { timeout: 30000 });
  const emb = await page.evaluate(() => ({
    n: parseInt(sdf3d.split('\n')[3].substring(0, 3)),
    hasCanvas: !!document.querySelector('#viewer3d canvas'),
  }));
  check('embed: 21 atoms (13 heavy + 8 H attached in 3D)', emb.n === 21, String(emb.n));
  check('embed: 3Dmol canvas created', emb.hasCanvas);

  await page.selectOption('#engineSel', 'MMFF94s');
  await page.evaluate(() => optimize3D());
  await page.waitForFunction(() => document.getElementById('status3d').textContent.includes('E(MMFF94s)'), null, { timeout: 60000 });
  // deterministic-path identity at full precision, in-page
  const ident = await page.evaluate(() => {
    const mb = sourceMolblock;
    // replicate the app path: embed(seed 42) → optimize(that geometry)
    const res = webmm.generate_initial_coordinates_wasm(mb);
    const coords = res.get_coordinates();
    const na = parseInt(mb.trim().split('\n')[3].substring(0, 3));
    const nb = parseInt(mb.trim().split('\n')[3].substring(3, 6));
    const lines = mb.trim().split('\n');
    const atomLines = [];
    for (let i = 0; i < na; i++)
      atomLines.push(coords[i*3].toFixed(4).padStart(10) + coords[i*3+1].toFixed(4).padStart(10) + coords[i*3+2].toFixed(4).padStart(10) + lines[4+i].substring(30));
    const sdfEmbed = [...lines.slice(0, 4), ...atomLines, ...lines.slice(4+na, 4+na+nb), 'M  END'].join('\n');
    // NOTE: OptimizationOptions is consumed by value on each call — fresh object per call
    const o1 = new webmm.OptimizationOptions();
    o1.engine = 'MMFF94s'; o1.set_max_iterations(500);
    const eApp = webmm.optimize_from_sdf(sdfEmbed, o1).final_energy;
    const o2 = new webmm.OptimizationOptions();
    o2.engine = 'MMFF94s'; o2.set_max_iterations(500);
    const eDirect = webmm.optimize_from_sdf(mb, o2).final_energy;  // internal ETKDG(42) path
    return { eApp, eDirect, d: Math.abs(eApp - eDirect) };
  });
  check('MMFF94s: app path == optimize_from_sdf(2D) identity (|ΔE| < 1e-5)',
        ident.d < 1e-5, `Δ = ${ident.d.toExponential(2)} (E = ${ident.eApp.toFixed(6)})`);

  await page.selectOption('#engineSel', 'GFNFF');
  await page.evaluate(() => optimize3D());
  await page.waitForFunction(() => document.getElementById('status3d').textContent.includes('E(GFNFF)'), null, { timeout: 120000 });
  const gfn13 = await optimizeWith(page, 'GFNFF', await page.evaluate(() => sdf3d));
  check('GFN-FF (13 atoms): converged + 9 terms + sum == total',
        gfn13.conv && Object.keys(gfn13.terms).length === 9 && Math.abs(gfn13.sum - gfn13.E) < 1e-4,
        `E=${gfn13.E.toFixed(3)} kcal/mol`);
  check('energy panel renders 9 terms + total', await page.evaluate(() =>
        document.getElementById('energyPanel').style.display === 'block' &&
        document.getElementById('eterms').rows.length === 10));

  // ---- (2) 24-atom demo caffeine ----
  await page.evaluate(s => { document.getElementById('input').value = s; process(true); }, CAFF_SDF3D);
  await page.evaluate(() => embed3D());
  await page.waitForFunction(() => document.getElementById('status3d').textContent.includes('3D coordinates present'), null, { timeout: 10000 });
  await page.selectOption('#engineSel', 'MMFF94s');
  await page.evaluate(() => optimize3D());
  await page.waitForFunction(() => document.getElementById('status3d').textContent.includes('E(MMFF94s)'), null, { timeout: 60000 });
  const s24 = await page.evaluate(() => document.getElementById('status3d').textContent);
  const e24 = parseFloat(s24.match(/E\(MMFF94s\) = (-?[\d.]+)/)[1]);
  check('24-atom caffeine MMFF94s reproduces demo E = -123.49', Math.abs(e24 - (-123.49)) < 0.5, `E = ${e24}`);

  // GFN-FF on the 24-atom molecule → XYZ → xtb cross-check
  await page.selectOption('#engineSel', 'GFNFF');
  await page.evaluate(() => optimize3D());
  await page.waitForFunction(() => document.getElementById('status3d').textContent.includes('E(GFNFF)'), null, { timeout: 120000 });
  const gfn24 = await optimizeWith(page, 'GFNFF', await page.evaluate(() => sdf3d));
  const dl = await (async () => { const p = page.waitForEvent('download'); await page.evaluate(() => exportXYZ()); return p; })();
  const xyzPath = '/tmp/m1_caffeine.xyz';
  await dl.saveAs(xyzPath);
  const out = execSync(`rm -rf /tmp/m1_xtb_run && mkdir -p /tmp/m1_xtb_run && cp ${xyzPath} /tmp/m1_xtb_run/in.xyz && cd /tmp/m1_xtb_run && ${XTB} in.xyz --gfnff 2>&1`, { timeout: 60000 }).toString();
  const xtbE = parseFloat(out.match(/TOTAL ENERGY\s+(-?\d+\.\d+)/)[1]);
  const dE = Math.abs(xtbE - gfn24.E / 627.5094740631);
  check('GFN-FF vs xtb @ same geometry (24-atom)', dE < 0.02,
        `webmm ${(gfn24.E / 627.5094740631).toFixed(6)} Eh vs xtb ${xtbE.toFixed(6)} Eh, Δ=${dE.toExponential(2)}`);

  // ---- (3) exports ----
  const dl2 = await (async () => { const p = page.waitForEvent('download'); await page.evaluate(() => exportSDF3D()); return p; })();
  const sdfPath = '/tmp/m1_caffeine.sdf';
  await dl2.saveAs(sdfPath);
  const sdfText = fs.readFileSync(sdfPath, 'utf8');
  const rd = await page.evaluate(s => {
    const m = rdkitModule.get_mol(s);
    const ok = m && m.is_valid();
    const n = ok ? JSON.parse(m.get_descriptors()).NumAtoms : 0;
    if (m) m.delete();
    return { ok, n };
  }, sdfText);
  check('exported SDF readable by RDKit.js (24 atoms)', rd.ok && rd.n === 24, rd.n + ' atoms');

  const dl3 = await (async () => { const p = page.waitForEvent('download'); await page.evaluate(() => exportPNG3D()); return p; })();
  const pngPath = '/tmp/m1_caffeine_3d.png';
  await dl3.saveAs(pngPath);
  check('PNG export saved', fs.existsSync(pngPath) && fs.statSync(pngPath).size > 10000, fs.statSync(pngPath).size + ' bytes');

  await page.screenshot({ path: '/tmp/app_m1.png', fullPage: false });
  check('zero console errors', errors.length === 0, errors.slice(0, 3).join(' | '));

  console.log(`\nRESULT: ${pass} passed, ${fail} failed`);
  await browser.close();
  process.exit(fail ? 1 : 0);
})();
