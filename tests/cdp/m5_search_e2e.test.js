// M5 acceptance: Search tab end-to-end against the frozen references.
// Covers: demo-library load + localStorage restore + clear; similarity
// parity (2 queries × 5 fingerprints × 55 entries, exact ===); the UI sim
// path (threshold/sort/self-top); substructure parity (4 queries, hit-name
// sets); Shape (3D) search (self 100% top, color/pharm/combo columns, the
// ≥60% pharm filter, aligned-pose row reload with auto-embed); RGD against
// rgd_refs.json (hit names + fragment SMILES, exact); Auto core; scaffold
// frequency (counts sum = cyclic count, row reload); and the pharmacophore
// query editor (default feature set, ±1.5 Å/N=1 hits, Conf column, row
// reload). Needs the usual :8901 repo-root server (see README).
const { chromium } = require('/opt/homebrew/lib/node_modules/playwright');
const fs = require('fs');
const path = require('path');
const EXE = '/Users/rliu/Library/Caches/ms-playwright/chromium_headless_shell-1234/chrome-headless-shell-mac-arm64/chrome-headless-shell';
const REFS = JSON.parse(fs.readFileSync(path.join(__dirname, '../fixtures/lbdd/search_refs.json'), 'utf8'));
const RGD = JSON.parse(fs.readFileSync(path.join(__dirname, '../fixtures/lbdd/rgd_refs.json'), 'utf8'));
const ASPIRIN = 'CC(=O)Oc1ccccc1C(=O)O';

let pass = 0, fail = 0;
const failures = [];
const check = (name, cond, detail = '') => {
  if (cond) { pass++; console.log(`  ✓ ${name}${detail ? ' — ' + detail : ''}`); }
  else { fail++; failures.push(name); console.log(`  ✗ ${name}${detail ? ' — ' + detail : ''}`); }
};
const eqSet = (a, b) => a.length === b.length && [...a].sort().join('|') === [...b].sort().join('|');

(async () => {
  const browser = await chromium.launch({ executablePath: EXE });
  const page = await browser.newPage({ viewport: { width: 1440, height: 900 } });
  const errors = [];
  page.on('console', m => { if (m.type() === 'error') errors.push(m.text().slice(0, 200)); });
  page.on('pageerror', e => errors.push('PAGEERROR: ' + e.message.slice(0, 300)));

  await page.goto('http://localhost:8901/app/index.html', { waitUntil: 'load' });
  await page.waitForFunction(() => document.getElementById('rdkitVersion').textContent !== 'Loading...', null, { timeout: 60000 });
  await page.waitForFunction(() => window.webmm !== undefined, null, { timeout: 60000 });

  // ---- library load + persistence ----
  console.log('library load / persistence:');
  await page.evaluate(() => switchMode('search'));
  await page.evaluate(() => loadDemoLibrary());
  await page.waitForFunction(() => document.getElementById('searchStatus').textContent.includes('55'), null, { timeout: 60000 });
  check('demo library loads 55', true);
  await page.reload({ waitUntil: 'load' });
  await page.waitForFunction(() => document.getElementById('rdkitVersion').textContent !== 'Loading...', null, { timeout: 60000 });
  await page.waitForFunction(() => window.webmm !== undefined, null, { timeout: 60000 });
  await page.evaluate(() => switchMode('search'));
  await page.waitForFunction(() => document.getElementById('searchStatus').textContent.includes('restored'), null, { timeout: 30000 });
  const restored = await page.evaluate(() => (window.__search.getState() || []).length);
  check('library auto-restored after reload (localStorage)', restored === 55, 'n=' + restored);

  // ---- similarity parity (exact doubles equality by construction) ----
  console.log('similarity parity vs search_refs.json:');
  const simParity = await page.evaluate((refs) => {
    const fps = { morgan: 'get_morgan_fp', rdkit: 'get_rdkit_fp', maccs: 'get_maccs_fp', atompair: 'get_atom_pair_fp', topologicaltorsion: 'get_topological_torsion_fp' };
    const lib = window.__search.getState();
    const byName = {}; lib.forEach(e => byName[e.name] = e);
    let mism = [], checked = 0;
    for (const [qname, q] of Object.entries(refs.queries)) {
      if (q.kind !== 'sim') continue;
      const qm = rdkitModule.get_mol(q.smiles);
      for (const [fpk, fn] of Object.entries(fps)) {
        const qbits = qm[fn]();
        refs.library.forEach((name, i) => {
          const e = byName[name];
          if (!e) { mism.push(qname + '/' + fpk + '/' + name + ' MISSING'); return; }
          const t = window.__search.tanimotoBits(qbits, e.fps[fpk]);
          checked++;
          if (t !== q.tanimoto[fpk][i]) mism.push(`${qname}/${fpk}/${name}: page=${t} ref=${q.tanimoto[fpk][i]}`);
        });
      }
      qm.delete();
    }
    return { mism, checked };
  }, REFS);
  check(`all ${simParity.checked} Tanimoto values exactly equal the refs`, simParity.mism.length === 0, simParity.mism.slice(0, 3).join(' ; ') || 'exact');

  // ---- UI similarity path: threshold, descending sort, self-match top ----
  console.log('UI similarity path:');
  await page.evaluate((q) => {
    document.getElementById('searchMode').value = 'sim'; onSearchModeChange();
    document.getElementById('searchQuery').value = q; runSearch();
  }, ASPIRIN);
  await page.waitForFunction(() => document.getElementById('searchRows').children.length > 0, null, { timeout: 15000 });
  const simUI = await page.evaluate(() => ({
    top: document.querySelector('#searchRows tr').children[1].textContent,
    topScore: document.querySelector('#searchRows tr').children[4].textContent,
    scores: [...document.querySelectorAll('#searchRows tr')].map(tr => parseFloat(tr.children[4].textContent)),
  }));
  check('self-match tops the table at 100.0%', simUI.top === 'aspirin' && simUI.topScore === '100.0%', simUI.top + ' ' + simUI.topScore);
  check('rows sorted by descending Tanimoto', simUI.scores.every((v, i, a) => i === 0 || a[i - 1] >= v), 'n=' + simUI.scores.length);
  check('T ≥ 0.30 threshold respected', simUI.scores.every(v => v >= 30), 'min=' + Math.min(...simUI.scores));

  // ---- substructure parity ----
  console.log('substructure parity vs search_refs.json:');
  for (const [qn, q] of Object.entries(REFS.queries)) {
    if (q.kind !== 'sub') continue;
    await page.evaluate((qq) => {
      document.getElementById('searchMode').value = 'sub'; onSearchModeChange();
      document.getElementById('searchQuery').value = qq; runSearch();
    }, q.smiles);
    await page.waitForFunction(() => document.getElementById('searchResultStatus').textContent.includes('matching'), null, { timeout: 15000 });
    const names = await page.evaluate(() => [...document.querySelectorAll('#searchRows tr')].map(tr => tr.children[1].textContent));
    check(`${qn}: hit set exactly matches ref (${q.hits.length})`, eqSet(names, q.hits), eqSet(names, q.hits) ? '' : 'page=' + JSON.stringify(names));
  }

  // ---- Shape (3D) search ----
  console.log('shape (3D) search:');
  await page.evaluate((q) => {
    document.getElementById('searchMode').value = 'shape'; onSearchModeChange();
    document.getElementById('searchQuery').value = q; runSearch();
  }, ASPIRIN);
  await page.waitForFunction(() => /sorted by combo/.test(document.getElementById('searchResultStatus').textContent), null, { timeout: 180000 });
  const shape = await page.evaluate(() => ({
    top: document.querySelector('#searchRows tr').children[1].textContent,
    topScore: document.querySelector('#searchRows tr').children[4].textContent,
    colorShown: document.getElementById('colorCol').style.display !== 'none',
    pharmShown: document.getElementById('pharmCol').style.display !== 'none',
    comboShown: document.getElementById('comboCol').style.display !== 'none',
    status: document.getElementById('searchResultStatus').textContent,
  }));
  check('self-match tops the table at 100.0%', shape.top === 'aspirin' && shape.topScore === '100.0%', shape.top + ' ' + shape.topScore);
  check('Color T / Pharm / Combo columns shown', shape.colorShown && shape.pharmShown && shape.comboShown, shape.status.slice(0, 90));
  await page.evaluate(() => { document.getElementById('pharmFilter').value = '0.6'; runSearch(); });
  await page.waitForFunction(() => /pharm filter/.test(document.getElementById('searchResultStatus').textContent), null, { timeout: 180000 });
  const pfNames = await page.evaluate(() => [...document.querySelectorAll('#searchRows tr')].map(tr => tr.children[1].textContent));
  check('pharm ≥60% filter keeps salicylic, drops glucose', pfNames.includes('salicylic acid') && !pfNames.includes('glucose'), JSON.stringify(pfNames));
  await page.evaluate(() => { document.getElementById('pharmFilter').value = '0'; });

  // aligned-pose row reload: auto-embed + feature spheres on
  await page.evaluate(() => document.querySelector('#searchRows tr').dispatchEvent(new MouseEvent('click', { bubbles: true })));
  await page.waitForFunction(() => !!document.querySelector('#viewer3d canvas'), null, { timeout: 60000 });
  const rowLoad = await page.evaluate(() => ({
    single: document.getElementById('tabSingle').classList.contains('active'),
    output: document.getElementById('output').style.display,
    feat: document.getElementById('featSpheres').checked,
  }));
  check('row click: single view + aligned pose auto-embedded + Features on', rowLoad.single && rowLoad.output === 'block' && rowLoad.feat, JSON.stringify(rowLoad));

  // ---- RGD ----
  console.log('R-group decomposition:');
  await page.evaluate(() => switchMode('search'));
  await page.evaluate(() => { fillRgdExample(); runRGD(); });
  await page.waitForFunction(() => document.getElementById('rgdRows').children.length > 0, null, { timeout: 15000 });
  const rgd = await page.evaluate(() => ({
    status: document.getElementById('rgdStatus').textContent,
    rows: [...document.querySelectorAll('#rgdRows tr')].map(tr => [...tr.querySelectorAll('td')].map(td => td.textContent.trim()).slice(1)),
  }));
  const rgdMap = Object.fromEntries(rgd.rows.map(r => [r[0], r.slice(2)]));
  const rgdNames = Object.keys(RGD.refs).filter(k => RGD.refs[k] !== null);
  check('hit set == rgd_refs (6/55)', eqSet(Object.keys(rgdMap), rgdNames), rgd.status);
  let rgdFragOK = true, rgdFragBad = '';
  for (const n of rgdNames) {
    if (!eqSet(rgdMap[n] || [], RGD.refs[n])) { rgdFragOK = false; rgdFragBad = n + ': page=' + JSON.stringify(rgdMap[n]) + ' ref=' + JSON.stringify(RGD.refs[n]); break; }
  }
  check('fragment SMILES == rgd_refs', rgdFragOK, rgdFragBad || 'all exact');
  const autoCore = await page.evaluate(() => { document.getElementById('rgdCore').value = ''; autoRgdCore(); return document.getElementById('rgdCore').value; });
  check('Auto core derives the labeled para-benzene core', autoCore === '[*:1]c1ccc([*:2])cc1', autoCore);

  // ---- scaffold frequency ----
  console.log('scaffold frequency:');
  await page.evaluate(() => analyzeScaffolds());
  await page.waitForFunction(() => document.getElementById('scaffoldRows').children.length > 0, null, { timeout: 15000 });
  const scaf = await page.evaluate(() => ({
    status: document.getElementById('scaffoldStatus').textContent,
    counts: [...document.querySelectorAll('#scaffoldRows tr')].map(tr => parseInt(tr.children[2].textContent)),
  }));
  check('counts sum to the cyclic-molecule count (34)', scaf.counts.reduce((a, b) => a + b, 0) === 34, scaf.status.slice(0, 70));
  await page.evaluate(() => document.querySelector('#scaffoldRows tr').dispatchEvent(new MouseEvent('click', { bubbles: true })));
  await page.waitForTimeout(500);
  const scafLoad = await page.evaluate(() => ({
    single: document.getElementById('tabSingle').classList.contains('active'),
    output: document.getElementById('output').style.display,
    input: document.getElementById('input').value,
  }));
  check('row click loads the scaffold into the single view', scafLoad.single && scafLoad.output === 'block' && scafLoad.input.length > 2, scafLoad.input.slice(0, 40));

  // ---- pharmacophore query editor ----
  console.log('pharmacophore query:');
  await page.evaluate((s) => {
    switchMode('single');
    document.getElementById('input').value = s; process(true);
  }, ASPIRIN);
  await page.waitForTimeout(400);
  await page.evaluate(() => embed3D());
  await page.waitForFunction(() => !!document.querySelector('#viewer3d canvas'), null, { timeout: 60000 });
  await page.evaluate(() => switchMode('search'));
  await page.evaluate(() => buildPharmQuery());
  await page.waitForTimeout(300);
  const pq = await page.evaluate(() => ({
    feats: document.querySelectorAll('#pharmQFeatures input[type=checkbox]').length,
    checked: document.querySelectorAll('#pharmQFeatures input[type=checkbox]:checked').length,
    dist: document.getElementById('pharmQDist').textContent.length,
  }));
  check('query built: 8 features, default compact set of 4 checked', pq.feats === 8 && pq.checked === 4, JSON.stringify(pq));
  check('distance-matrix preview rendered', pq.dist > 10);
  await page.evaluate(() => { document.getElementById('pharmTol').value = '1.5'; document.getElementById('pharmConfs').value = '1'; screenPharmQuery(); });
  await page.waitForFunction(() => /match the pharmacophore query/.test(document.getElementById('pharmQStatus').textContent), null, { timeout: 180000 });
  const pharmHits = await page.evaluate(() => [...document.querySelectorAll('#pharmQRows tr')].map(tr => ({ name: tr.children[1].textContent, conf: tr.children[4].textContent })));
  check('±1.5 Å / N=1 hits exactly {aspirin, salicylic acid} (glucose negative control implied)', eqSet(pharmHits.map(h => h.name), ['aspirin', 'salicylic acid']), JSON.stringify(pharmHits));
  check('Conf column is 1 for the single-conformer screen', pharmHits.every(h => h.conf === '1'));
  await page.evaluate(() => document.querySelector('#pharmQRows tr').dispatchEvent(new MouseEvent('click', { bubbles: true })));
  await page.waitForTimeout(600);
  check('row click loads the hit into the single view', await page.evaluate(() => document.getElementById('tabSingle').classList.contains('active') && document.getElementById('output').style.display === 'block'));

  // ---- clear library (also leaves clean storage for the next run) ----
  await page.evaluate(() => switchMode('search'));
  await page.evaluate(() => clearSearchLibrary());
  const cleared = await page.evaluate(() => ({ status: document.getElementById('searchStatus').textContent, stored: localStorage.getItem('wb-searchLibrary') }));
  check('Clear wipes the library and its persisted copy', /cleared/.test(cleared.status) && cleared.stored === null, cleared.status);

  console.log('\npage errors:', JSON.stringify(errors, null, 1));
  check('zero page errors', errors.length === 0, errors.slice(0, 3).join(' ; '));
  console.log(`\nRESULT: ${pass} passed, ${fail} failed`);
  if (failures.length) console.log('failed: ' + failures.join(' | '));
  await browser.close();
  process.exit(fail ? 1 : 0);
})();
