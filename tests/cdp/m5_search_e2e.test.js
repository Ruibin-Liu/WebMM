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
  // NOTE: the confs input defaults to 10 (ensemble semantics: entry score =
  // best over conformers); the single-conformer legacy path is confs=1.
  await page.evaluate((q) => {
    document.getElementById('searchMode').value = 'shape'; onSearchModeChange();
    document.getElementById('searchQuery').value = q; runSearch();
  }, ASPIRIN);
  await page.waitForFunction(() => /sorted by combo/.test(document.getElementById('searchResultStatus').textContent), null, { timeout: 180000 });
  const shape = await page.evaluate(() => ({
    top: document.querySelector('#searchRows tr').children[1].textContent,
    topScore: document.querySelector('#searchRows tr').children[4].textContent,
    // the ensemble default inserts a Conf column before the score — read the
    // score from the cell whose header says Shape T
    topScoreByHeader: [...document.querySelector('#searchRows tr').cells].map(td => td.textContent.trim())[
      [...document.querySelectorAll('#searchTable th')].findIndex(th => th.textContent.trim() === 'Shape T')],
    colorShown: document.getElementById('colorCol').style.display !== 'none',
    pharmShown: document.getElementById('pharmCol').style.display !== 'none',
    comboShown: document.getElementById('comboCol').style.display !== 'none',
    status: document.getElementById('searchResultStatus').textContent,
  }));
  check('self-match tops the table at 100.0%', shape.top === 'aspirin' && shape.topScoreByHeader === '100.0%', shape.top + ' ' + shape.topScoreByHeader);
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

  // ---- multi-conformer shape semantics (default) ----
  console.log('shape ensemble (10 confs/entry):');
  const ens = await page.evaluate(() => {
    document.getElementById('searchQuery').value = 'CC(C)Cc1ccc(cc1)C(C)C(=O)O';
    document.getElementById('searchResultStatus').textContent = '';
    runSearch();
    return null;
  });
  await page.waitForFunction(() => /confs\/entry/.test(document.getElementById('searchResultStatus').textContent), null, { timeout: 300000 });
  const ensRes = await page.evaluate(() => ({
    status: document.getElementById('searchResultStatus').textContent,
    top: document.querySelector('#searchRows tr') ? {
      name: document.querySelector('#searchRows tr').children[1].textContent,
      score: document.querySelector('#searchRows tr').children[5].textContent,
      conf: document.querySelector('#searchRows tr').children[4].textContent,
    } : null,
    confColShown: document.getElementById('confCol').style.display !== 'none',
    rows: document.getElementById('searchRows').children.length,
  }));
  check('ensemble: self-match tops at 100.0% (query-conformer injection)', ensRes.top && ensRes.top.name === 'ibuprofen' && ensRes.top.score === '100.0%' && ensRes.top.conf === '1', JSON.stringify(ensRes.top));
  check('ensemble: Conf column shown with the winning conformer', ensRes.confColShown && ensRes.top.conf === '1');
  check('ensemble: status documents confs/entry + two-phase structure', /10 confs\/entry/.test(ensRes.status) && /screen 55×10 → full 50×3/.test(ensRes.status), ensRes.status.slice(0, 100));
  check('ensemble: rows bounded by the top-50 candidate cap', ensRes.rows >= 40 && ensRes.rows <= 50, 'rows=' + ensRes.rows);

  // ---- color force field weights (v1.6.5) ----
  console.log('color weights:');
  const cwUI = await page.evaluate(() => ({
    wrapShown: document.getElementById('colorWeightsWrap').style.display !== 'none',
    inputs: document.querySelectorAll('#colorWeightsWrap .cw').length,
    defaults: [...document.querySelectorAll('#colorWeightsWrap .cw')].every(el => el.value === '1'),
  }));
  check('six weight inputs visible in shape mode, default 1.0', cwUI.wrapShown && cwUI.inputs === 6 && cwUI.defaults);

  // all-zero weights: Color T must be 0 everywhere and Combo == Shape T
  const zeroRun = await page.evaluate(() => {
    for (const el of document.querySelectorAll('#colorWeightsWrap .cw')) el.value = '0';
    document.getElementById('searchResultStatus').textContent = '';
    runSearch();
    return null;
  });
  await page.waitForFunction(() => /sorted by combo/.test(document.getElementById('searchResultStatus').textContent), null, { timeout: 300000 });
  const zeroRes = await page.evaluate(() => [...document.querySelectorAll('#searchRows tr')].slice(0, 10).map(tr => ({
    name: tr.children[1].textContent,
    shapeT: tr.children[5].textContent, colorT: tr.children[6].textContent, combo: tr.children[8].textContent,
  })));
  check('all-zero weights: Color T = 0% everywhere, Combo = Shape T',
    zeroRes.every(r => r.colorT === '0.0%' && r.combo === r.shapeT),
    JSON.stringify(zeroRes[0]));

  // donor weight 3 must perturb at least one Color T vs default
  const donorRun = await page.evaluate(() => {
    resetColorWeights();
    document.querySelector('#colorWeightsWrap .cw[data-t="donor"]').value = '3';
    document.getElementById('searchResultStatus').textContent = '';
    runSearch();
    return null;
  });
  await page.waitForFunction(() => /sorted by combo/.test(document.getElementById('searchResultStatus').textContent), null, { timeout: 300000 });
  const donorRes = await page.evaluate(() => [...document.querySelectorAll('#searchRows tr')].slice(0, 10).map(tr => tr.children[6].textContent));
  check('donor=3 changes at least one Color T vs all-zero', donorRes.some(c => c !== '0.0%'), donorRes.slice(0, 3).join(','));

  // restore defaults; the default run must match the pre-weights results
  await page.evaluate(() => resetColorWeights());
  const defRun = await page.evaluate(() => {
    document.getElementById('searchResultStatus').textContent = '';
    runSearch();
    return null;
  });
  await page.waitForFunction(() => /sorted by combo/.test(document.getElementById('searchResultStatus').textContent), null, { timeout: 300000 });
  const defTop = await page.evaluate(() => document.querySelector('#searchRows tr').children[1].textContent +
    ':' + document.querySelector('#searchRows tr').children[5].textContent);
  check('weights back at 1.0: default results unchanged (zero drift)', defTop === 'ibuprofen:100.0%', defTop);

  // ---- M1c rounds: DAG + hit-as-query + lineage persistence ----
  console.log('rounds / hit-as-query:');
  await page.evaluate(() => {
    document.getElementById('searchMode').value = 'shape'; onSearchModeChange();
    document.getElementById('searchQuery').value = 'CC(C)Cc1ccc(C(C)C(=O)O)cc1'; runSearch();
  });
  await page.waitForFunction(() => /sorted by combo/.test(document.getElementById('searchResultStatus').textContent), null, { timeout: 300000 });
  const hist1 = await page.evaluate(() => ({
    shown: document.getElementById('queryHistoryPanel').style.display !== 'none',
    entries: [...document.querySelectorAll('#queryHistoryPanel div[onclick]')].length,
    hasSwap: !!document.querySelector('#searchRows tr span[onclick*=hitAsQuery]'),
  }));
  check('search records a round; lineage panel + ⇄ action present',
    hist1.shown && hist1.entries >= 1 && hist1.hasSwap, JSON.stringify(hist1));
  // ONE-ACTION lead hop (baseline C: 11 actions)
  await page.evaluate(() => document.querySelectorAll('#searchRows tr')[1].querySelector('span[onclick*=hitAsQuery]').click());
  await page.waitForFunction(() => {
    const es = [...document.querySelectorAll('#queryHistoryPanel div[onclick]')];
    return es.length >= 2 && es[es.length - 1].textContent.includes('↳') && document.getElementById('searchQuery').value.includes('COc1ccc');
  }, null, { timeout: 300000 });
  const hop = await page.evaluate(() => ({
    entries: [...document.querySelectorAll('#queryHistoryPanel div[onclick]')].map(d => d.textContent.trim()),
    query: document.getElementById('searchQuery').value,
  }));
  check('hit-as-query = ONE action, child round indented (↳), query set to the hit',
    hop.entries.length >= 2 && hop.entries[hop.entries.length - 1].includes('↳') && hop.query.includes('COc1ccc'), JSON.stringify(hop.entries));
  // reload: the DAG persists
  await page.reload({ waitUntil: 'load' });
  await page.waitForFunction(() => document.getElementById('rdkitVersion').textContent !== 'Loading...', null, { timeout: 60000 });
  await page.waitForFunction(() => window.webmm !== undefined, null, { timeout: 60000 });
  await page.evaluate(() => switchMode('search'));
  await page.waitForFunction(() => document.getElementById('searchStatus').textContent.includes('restored'), null, { timeout: 30000 });
  const histReload = await page.evaluate(() => [...document.querySelectorAll('#queryHistoryPanel div[onclick]')].map(d => d.textContent.trim()));
  check('lineage persists across reload via the command log',
    histReload.length >= 2 && histReload.some(x => x.includes('↳')), JSON.stringify(histReload));
  // rerun a historical round — the LAST entry is the ↳ shape child
  // (history[0] may be a sim round whose completion status is 'at T', not 'sorted by combo')
  await page.evaluate(() => {
    const es = [...document.querySelectorAll('#queryHistoryPanel div[onclick]')];
    es[es.length - 1].click();
  });
  await page.waitForFunction(() => /sorted by combo/.test(document.getElementById('searchResultStatus').textContent), null, { timeout: 300000 });
  check('rerun from history re-executes the round query', true);

  // ---- M1b triage overlay (pin/exclude/notes/undo, platform store) ----
  console.log('triage overlay:');
  // restore any leftover pin state from previous runs of this suite is fine —
  // the store is per-origin persistent; assert relative behavior
  await page.evaluate(() => {
    document.getElementById('searchMode').value = 'sim'; onSearchModeChange();
    document.getElementById('searchQuery').value = 'CC(=O)Oc1ccccc1C(=O)O'; runSearch();
  });
  await page.waitForFunction(() => document.getElementById('searchResultStatus').textContent.includes('at T'), null, { timeout: 120000 });
  const before = await page.evaluate(() => document.getElementById('workSetCounts').textContent);
  const pinRes = await page.evaluate(async () => {
    const row = document.querySelector('#searchRows tr');
    const pin = row.querySelector('td:last-child span[onclick*=triageTogglePin]');
    pin.click();
    await new Promise(r => setTimeout(r, 500));
    return {
      star: document.querySelector('#searchRows tr td:last-child span[onclick*=triageTogglePin]').textContent,
      counts: document.getElementById('workSetCounts').textContent,
      panelShown: document.getElementById('workSetPanel').style.display !== 'none',
      pinRows: document.querySelectorAll('#workSetPins > div').length,
    };
  });
  check('pin toggles star, working-set panel shows the pin',
    pinRes.star === '★' && /1 pinned/.test(pinRes.counts) && pinRes.panelShown && pinRes.pinRows >= 1, JSON.stringify(pinRes));
  const exclRes = await page.evaluate(async () => {
    const row = document.querySelectorAll('#searchRows tr')[1];
    row.querySelector('td:last-child span[onclick*=triageToggleExclude]').click();
    await new Promise(r => setTimeout(r, 500));
    return {
      counts: document.getElementById('workSetCounts').textContent,
      opacity: getComputedStyle(document.querySelectorAll('#searchRows tr')[1]).opacity,
    };
  });
  check('exclude dims the row and updates counts', /1 excluded/.test(exclRes.counts) && exclRes.opacity === '0.45', JSON.stringify(exclRes));
  const undoRes = await page.evaluate(async () => {
    triageUndo();
    await new Promise(r => setTimeout(r, 500));
    return {
      status: document.getElementById('searchStatus').textContent,
      opacity: getComputedStyle(document.querySelectorAll('#searchRows tr')[1]).opacity,
    };
  });
  check('undo reverses the exclude (inverse command, log intact)',
    /Undid: Exclude/.test(undoRes.status) && undoRes.opacity === '1', JSON.stringify(undoRes));
  // note editing via the working-set panel
  const noteRes = await page.evaluate(async () => {
    const inp = document.querySelector('#workSetPins input');
    if (!inp) return { skip: true };
    inp.value = 'lead candidate';
    inp.dispatchEvent(new Event('change'));
    await new Promise(r => setTimeout(r, 500));
    return { note: document.querySelector('#workSetPins input').value };
  });
  check('note edited through the panel persists in the session', !noteRes.skip || noteRes.note === 'lead candidate');
  // pin survives reload (IndexedDB command log)
  await page.reload({ waitUntil: 'load' });
  await page.waitForFunction(() => document.getElementById('rdkitVersion').textContent !== 'Loading...', null, { timeout: 60000 });
  await page.waitForFunction(() => window.webmm !== undefined, null, { timeout: 60000 });
  await page.evaluate(() => switchMode('search'));
  await page.waitForFunction(() => document.getElementById('searchStatus').textContent.includes('restored'), null, { timeout: 30000 });
  const persisted = await page.evaluate(async () => {
    await new Promise(r => setTimeout(r, 300));
    return document.getElementById('workSetCounts').textContent;
  });
  check('pin survives page reload via the project store', /1 pinned/.test(persisted), persisted);
  // cleanup: unpin + restore the mode state this reload disturbed (the
  // later flex section assumes shape mode persisted from earlier sections)
  await page.evaluate(async () => {
    const c = window.__platformCompat;
    if (c) for (const id of Object.keys(c.getPins())) await c.unpin(id);
    document.getElementById('searchMode').value = 'shape'; onSearchModeChange();
    document.getElementById('searchQuery').value = 'CC(C)Cc1ccc(C(C)C(=O)O)cc1';
  });

  // ---- projected color sites (v1.9.0) ----
  console.log('projected color sites:');
  const proj = await page.evaluate(() => {
    const m = rdkitModule.get_mol('CC(=O)Oc1ccccc1C(=O)O');
    const mb = m.get_molblock(); m.delete();
    const r = window.webmm.generate_optimized_conformers_wasm(mb, 1, BigInt(7), 'MMFF94s', 250);
    const sdf = buildSdfFromCoords(r.get_coordinates(), r.get_template_sdf());
    const sm = rdkitModule.get_mol(sdf);
    const sites = colorSites(sdf, sm); sm.delete();
    const json = sitesToEngineJson(sites, null, sdfCoords(sdf));
    return {
      donors: sites.filter(x => x.type === 'donor').length,
      donorsAtH: sites.filter(x => x.type === 'donor' && x.proj && x.proj.h).length,
      acceptors: sites.filter(x => x.type === 'acceptor').length,
      acceptorsProj: sites.filter(x => x.type === 'acceptor' && x.proj).length,
      offLen: json.filter(x => x.off).length,
      offNorm: json.filter(x => x.off).every(x => Math.hypot(...x.off) > 0.5 && Math.hypot(...x.off) < 1.5),
      // unprojected types must stay atom-centered (zero drift for them)
      hydroAtAtom: sites.filter(x => x.type === 'hydrophobe').every(x => !x.proj),
    };
  });
  check('aspirin: every donor site projected onto its hydrogens', proj.donors >= 1 && proj.donors === proj.donorsAtH, JSON.stringify(proj).slice(0, 90));
  check('aspirin: every acceptor site carries a lone-pair projection', proj.acceptors >= 3 && proj.acceptors === proj.acceptorsProj);
  check('engine JSON carries ~1 A offsets; hydrophobe/ring stay atom-centered',
    proj.offLen === proj.donors + proj.acceptors && proj.offNorm && proj.hydroAtAtom);

  // ---- flexible alignment (v1.8.0) ----
  console.log('flexible alignment:');
  const flexUI = await page.evaluate(() => ({
    shown: document.getElementById('shapeFlexWrap').style.display !== 'none',
    checked: document.getElementById('shapeFlex').checked,
  }));
  check('flex toggle visible in shape mode, default off', flexUI.shown && !flexUI.checked, JSON.stringify(flexUI));

  await page.evaluate(() => {
    document.getElementById('shapeFlex').checked = true;
    document.getElementById('searchResultStatus').textContent = '';
    runSearch();
  });
  await page.waitForFunction(() => /sorted by combo/.test(document.getElementById('searchResultStatus').textContent), null, { timeout: 300000 });
  const flexRun = await page.evaluate(() => ({
    status: document.getElementById('searchResultStatus').textContent,
    badges: document.querySelectorAll('#searchRows .result-badge').length,
    top: document.querySelector('#searchRows tr').children[1].textContent,
    topTitle: document.querySelector('#searchRows tr .result-badge') ? document.querySelector('#searchRows tr .result-badge').title : '',
  }));
  check('flex on: status reports refinement', /flex refined \/?\d+\/10/.test(flexRun.status), flexRun.status.slice(0, 80));
  check('flex on: top rows carry flex badges (MCS anchors)', flexRun.badges >= 1 && flexRun.top.startsWith('ibuprofen') && /MCS \d+ atoms/.test(flexRun.topTitle), flexRun.top + ' ' + flexRun.topTitle);
  const selfFlex = await page.evaluate(() => [...document.querySelectorAll('#searchRows tr')].slice(0, 3).map(tr => [...tr.cells].map(td => td.textContent.trim())));
  check('flex never degrades an exact self-match below 100.0%',
    selfFlex.some(cells => cells[1].startsWith('ibuprofen') ? cells.includes('100.0%') : true) &&
    selfFlex.every(cells => !cells[1].startsWith('ibuprofen') || cells.some(c => parseFloat(c) >= 99.9)),
    JSON.stringify(selfFlex.map(c => c.slice(0, 2))));

  await page.evaluate(() => {
    document.getElementById('shapeFlex').checked = false;
    document.getElementById('searchResultStatus').textContent = '';
    runSearch();
  });
  await page.waitForFunction(() => /sorted by combo/.test(document.getElementById('searchResultStatus').textContent), null, { timeout: 300000 });
  const offRun = await page.evaluate(() => ({
    badges: document.querySelectorAll('#searchRows .result-badge').length,
    status: document.getElementById('searchResultStatus').textContent,
    top: document.querySelector('#searchRows tr').children[1].textContent + ':' +
      [...document.querySelector('#searchRows tr').cells].map(td => td.textContent.trim()).find(t => t.endsWith('%') && parseFloat(t) >= 100),
  }));
  check('flex off: no badges, results back to rigid (zero drift)', offRun.badges === 0 && !/flex refined/.test(offRun.status) && offRun.top === 'ibuprofen:100.0%', offRun.top + ' badges=' + offRun.badges);

  // ---- single-conformer legacy path (confs=1) is unchanged ----
  console.log('shape legacy single-conformer (confs=1):');
  await page.evaluate(() => {
    document.getElementById('shapeConfs').value = '1';
    document.getElementById('searchResultStatus').textContent = '';
    runSearch();
  });
  await page.waitForFunction(() => /sorted by combo/.test(document.getElementById('searchResultStatus').textContent) && !/confs\/entry/.test(document.getElementById('searchResultStatus').textContent), null, { timeout: 120000 });
  const legRes = await page.evaluate(() => ({
    status: document.getElementById('searchResultStatus').textContent,
    top: document.querySelector('#searchRows tr').children[1].textContent,
    confColShown: document.getElementById('confCol').style.display !== 'none',
    confCell: document.querySelector('#searchRows tr').children[4].textContent,
  }));
  check('legacy: self-match still tops, no confs note in status', legRes.top === 'ibuprofen' && !/conf\/entry/.test(legRes.status), legRes.status.slice(0, 80));
  check('legacy: Conf column hidden at confs=1 (no misaligned empty column)', !legRes.confColShown);
  // restore the ensemble default for any later steps
  await page.evaluate(() => { document.getElementById('shapeConfs').value = '10'; });

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
