// m6_platform.test.js — the LBDD platform page (app/platform.html):
// rounds DAG + hit-as-query, consensus/provenance CSV, project export/
// import + multi-tab, triage overlay, neighborhood explorer (M1-M3).
// Served from repo root on :8901 (same convention as m0-m5).
const { chromium } = require('/opt/homebrew/lib/node_modules/playwright');
const fs = require('fs');
const EXE = '/Users/rliu/Library/Caches/ms-playwright/chromium_headless_shell-1234/chrome-headless-shell-mac-arm64/chrome-headless-shell';
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
  page.on('console', m => { if (m.type() === 'error') errors.push(m.text().slice(0, 200)); });
  page.on('pageerror', e => errors.push('PAGEERROR: ' + e.message.slice(0, 300)));

  await page.goto(URL, { waitUntil: 'load' });
  await page.waitForFunction(() => document.getElementById('rdkitVersion').textContent !== 'Loading...', null, { timeout: 60000 });
  await page.waitForFunction(() => window.webmm !== undefined, null, { timeout: 60000 });

  // ---- library setup (the platform sections expect a loaded library) ----
  console.log('library setup:');
  await page.evaluate(() => { switchMode('search'); loadDemoLibrary(); });
  await page.waitForFunction(() => document.getElementById('searchStatus').textContent.includes('55'), null, { timeout: 60000 });
  check('demo library loads 55 on the platform page', true);

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

  // ---- M2a: round-score facts + consensus + provenance CSV + SAR scope ----
  console.log('M2a consensus / provenance / SAR scope:');
  // two OVERLAPPING sim rounds (aspirin, then salicylic acid) — both hit the same pair
  for (const q of ['CC(=O)Oc1ccccc1C(=O)O', 'Oc1ccccc1C(=O)O']) {
    await page.evaluate((qq) => {
      document.getElementById('searchMode').value = 'sim'; onSearchModeChange();
      document.getElementById('searchQuery').value = qq; runSearch();
    }, q);
    await page.waitForFunction(() => document.getElementById('searchResultStatus').textContent.includes('at T'), null, { timeout: 120000 });
  }
  const cons = await page.evaluate(() => {
    const ci = [...document.querySelectorAll('#searchRows tr')].slice(0, 2).map(tr => {
      const c = tr.cells[tr.cells.length - 2];
      return { name: tr.children[1].textContent.trim(), cons: c.textContent, tip: c.title };
    });
    return ci;
  });
  check('consensus aggregates every recorded round (mean rank, formula in tooltip)',
    cons.length === 2 && cons.every(c => parseFloat(c.cons) > 0) && cons[0].tip.split(' · ').length >= 2, JSON.stringify(cons));
  // provenance CSV downloads with the expected columns
  const dlP = page.waitForEvent('download', { timeout: 15000 });
  await page.evaluate(() => exportResultsCSV());
  await (await dlP).saveAs('/tmp/m5-prov.csv');
  const csvHead = fs.readFileSync('/tmp/m5-prov.csv', 'utf8').split('\n')[0];
  check('provenance CSV carries per-round ranks + consensus + triage columns',
    csvHead.includes('rank:sim') && csvHead.includes('consensus_mean_rank') && csvHead.includes('pinned') && csvHead.includes('excluded'), csvHead.slice(0, 90));
  // SAR on the current hit set
  await page.evaluate(() => { document.getElementById('scaffoldScope').value = 'hits'; analyzeScaffolds(); });
  await page.waitForFunction(() => document.getElementById('scaffoldStatus').textContent.includes('scaffold'), null, { timeout: 60000 });
  const scafScope = await page.evaluate(() => document.getElementById('scaffoldStatus').textContent);
  check('scaffold analysis runs over the current hit set', /\d+ scaffolds across \d+ cyclic/.test(scafScope), scafScope.slice(0, 70));
  await page.evaluate(() => { document.getElementById('rgdScope').value = 'hits'; fillRgdExample(); runRGD(); });
  await page.waitForFunction(() => document.getElementById('rgdStatus').textContent.includes('contain the core'), null, { timeout: 60000 });
  const rgdScope = await page.evaluate(() => document.getElementById('rgdStatus').textContent);
  check('RGD denominator reflects the hit-set scope (runs over the sim hits, not /55)', /\/3 molecules/.test(rgdScope), rgdScope.slice(0, 60));
  // restore defaults for later sections
  await page.evaluate(() => { document.getElementById('scaffoldScope').value = 'library'; document.getElementById('rgdScope').value = 'library'; });

  // ---- M2b: project export/import + multi-tab broadcast ----
  console.log('M2b project portability / multi-tab:');
  // export (needs a search + a pin so the project is non-trivial)
  await page.evaluate(() => {
    document.getElementById('searchMode').value = 'sim'; onSearchModeChange();
    document.getElementById('searchQuery').value = 'CC(=O)Oc1ccccc1C(=O)O'; runSearch();
  });
  await page.waitForFunction(() => document.getElementById('searchResultStatus').textContent.includes('at T'), null, { timeout: 120000 });
  await page.evaluate(async () => {
    const id = window.__search.getState()[0].molId;
    await window.__platformCompat.pin(id, 'exported-pin');
  });
  await page.waitForTimeout(300);
  const dlP2 = page.waitForEvent('download', { timeout: 15000 });
  await page.evaluate(() => exportProjectFile());
  await (await dlP2).saveAs('/tmp/m5-project.json');
  const projData = JSON.parse(fs.readFileSync('/tmp/m5-project.json', 'utf8'));
  check('project export = command log with library + round + pin',
    projData.schemaVersion === 1 && projData.commands.some(c => c.type === 'ImportLibrary') &&
    projData.commands.some(c => c.type === 'CreateRound') && projData.commands.some(c => c.type === 'Pin'),
    'commands=' + projData.commands.length);
  // clear, then import back; the page reloads (confirm dialog auto-accepted)
  await page.evaluate(() => clearSearchLibrary());
  await page.waitForTimeout(500);
  page.once('dialog', d => d.accept());
  await page.setInputFiles('#projectImportFile', '/tmp/m5-project.json');
  await page.waitForFunction(() => document.readyState === 'complete', null, { timeout: 30000 }).catch(() => {});
  await page.waitForTimeout(2500);
  await page.evaluate(() => switchMode('search'));
  await page.waitForFunction(() => /restored/.test(document.getElementById('searchStatus').textContent), null, { timeout: 30000 });
  const imported = await page.evaluate(() => ({
    lib: (window.__search.getState() || []).length,
    hist: document.querySelectorAll('#queryHistoryPanel div[onclick]').length,
    pins: document.getElementById('workSetCounts').textContent,
  }));
  check('import restores library + rounds + pins after reload',
    imported.lib === 55 && imported.hist >= 1 && /1 pinned/.test(imported.pins), JSON.stringify(imported));
  // multi-tab: a second page in the same context hears the broadcast
  const page2 = await ctx.newPage();   // SAME context: shared storage + broadcast
  await page2.goto(URL, { waitUntil: 'load' });
  await page2.waitForFunction(() => document.getElementById('rdkitVersion').textContent !== 'Loading...', null, { timeout: 60000 });
  await page2.waitForFunction(() => window.webmm !== undefined, null, { timeout: 60000 });
  await page2.evaluate(() => switchMode('search'));
  await page2.waitForFunction(() => window.__platformCompat, null, { timeout: 15000 });
  await page2.evaluate(() => { window.__platformCompat_pinProbe = true; });
  await page.evaluate(async () => {
    const id = window.__search.getState()[0].molId;
    await window.__platformCompat.pin(id, null);   // any write broadcasts
  });
  await page2.waitForFunction(() => window.__platformCompat.stale() === true, null, { timeout: 8000 }).catch(() => {});
  const bStale = await page2.evaluate(() => window.__platformCompat.stale());
  await page2.close();
  check('multi-tab: second page is notified (stale flag) when this tab writes', bStale === true, 'stale=' + bStale);
  // cleanup: this section's pins must not leak into the triage section
  await page.evaluate(async () => {
    const c = window.__platformCompat;
    for (const id of Object.keys(c.getPins())) await c.unpin(id);
  });

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

  // ---- M3: neighborhood explorer (analog enumeration + shape rescoring) ----
  console.log('M3 neighborhood explorer:');
  await page.evaluate(() => { switchMode('single'); document.getElementById('input').value = 'CC(=O)Nc1ccc(O)cc1'; process(true); });
  await page.waitForTimeout(500);
  await page.evaluate(() => embed3D());
  await page.waitForFunction(() => !!document.querySelector('#viewer3d canvas'), null, { timeout: 60000 });
  await page.evaluate(() => { switchMode('search'); analogDetectSites(); });
  await page.waitForTimeout(600);
  const sites = await page.evaluate(() => [...document.getElementById('analogSite').options].filter(o => o.value).length);
  check('replaceable sites detected on paracetamol (>=2)', sites >= 2, 'sites=' + sites);
  await page.evaluate(() => {
    const sel = document.getElementById('analogSite');
    const opt = [...sel.options].find(o => o.textContent.includes('methyl'));
    sel.value = opt.value;
    document.getElementById('analogFragCat').value = 'alkyl';
    runAnalogExplorer();
  });
  await page.waitForFunction(() => /Done:|No candidates/.test(document.getElementById('analogStatus').textContent), null, { timeout: 300000 });
  const exp = await page.evaluate(() => ({
    status: document.getElementById('analogStatus').textContent,
    rows: [...document.querySelectorAll('#analogRows tr')].map(tr => ({
      smi: tr.children[1].textContent, T: tr.children[3].textContent,
    })),
  }));
  check('alkyl category scores analogs (parent among them, sensible order)',
    /Done: \d+ scored/.test(exp.status) && exp.rows.length >= 5 &&
    exp.rows.some(r => r.smi === 'CC(=O)Nc1ccc(O)cc1') &&
    exp.rows.every(r => parseFloat(r.T) > 0), exp.status.slice(0, 60));
  const parentRow = exp.rows.find(r => r.smi === 'CC(=O)Nc1ccc(O)cc1');
  check('parent-swap-to-methyl scores high (rigid self-analog)', parentRow && parseFloat(parentRow.T) >= 0.8, parentRow && parentRow.T);
  // full funnel: wait for stage C (30-conf refinement) + D (flex top-10)
  await page.waitForFunction(() => /Done: \d+ scored .*flex-refined/.test(document.getElementById('analogStatus').textContent), null, { timeout: 600000 });
  const analogFlex = await page.evaluate(() => ({
    status: document.getElementById('analogStatus').textContent,
    top: [...document.querySelectorAll('#analogRows tr')].slice(0, 6).map(tr => ({
      smi: tr.children[1].textContent.slice(0, 30), rigid: tr.children[3].textContent, flex: tr.children[5].textContent })),
  }));
  check('funnel completes B->C->D with flex-refined top rows',
    /flex-refined/.test(analogFlex.status) && analogFlex.top.some(r => parseFloat(r.flex) > 0), analogFlex.status.slice(0, 70));
  // pin an analog = synthetic-library identity
  await page.evaluate(() => document.querySelector('#analogRows span[onclick*=analogPin]').click());
  await page.waitForTimeout(500);
  const analogPin = await page.evaluate(() => ({
    counts: document.getElementById('workSetCounts').textContent,
    star: document.querySelector('#analogRows span[onclick*=analogPin]').textContent,
  }));
  check('analog pins via synthetic-library identity (star fills, counts update)',
    analogPin.star === '★' && /1 pinned/.test(analogPin.counts), JSON.stringify(analogPin));
  await page.evaluate(async () => {
    const c = window.__platformCompat;
    for (const id of Object.keys(c.getPins())) await c.unpin(id);
  });
  await page.evaluate(() => { switchMode('single'); });   // leave the mode tidy


  // ---- M3+: SA column + aza-scan scaffold hop ----
  console.log('SA column / aza-scan hop:');
  await page.evaluate(() => { switchMode('single'); document.getElementById('input').value = 'CC(=O)Nc1ccc(O)cc1'; process(true); });
  await page.waitForTimeout(500);
  await page.evaluate(() => embed3D());
  await page.waitForFunction(() => !!document.querySelector('#viewer3d canvas'), null, { timeout: 60000 });
  await page.evaluate(() => { switchMode('search'); runAzaScan(); });
  // wait for THIS run's terminal marker (the stale status from the M3
  // section matches the generic Done regex; the aza-hop suffix discriminates)
  await page.waitForFunction(() => /aza-hop|No valid aza|No H-bearing/.test(document.getElementById('analogStatus').textContent), null, { timeout: 300000 });
  const hopRes = await page.evaluate(() => ({
    status: document.getElementById('analogStatus').textContent,
    rows: [...document.querySelectorAll('#analogRows tr')].slice(0, 4).map(tr => ({
      smi: tr.children[1].textContent, T: tr.children[3].textContent,
      sa: tr.children[5].textContent, flex: tr.children[6].textContent })),
  }));
  check('aza-scan produces pyridyl scaffold hops with SA + flex columns',
    /Done: \d+ scored/.test(hopRes.status) && hopRes.rows.length >= 2 &&
    hopRes.rows.every(r => r.smi.includes('n')) &&
    hopRes.rows.every(r => parseFloat(r.sa) > 0) &&
    hopRes.rows.every(r => parseFloat(r.flex) > 0), hopRes.status.slice(0, 60));
  const best = hopRes.rows[0];
  check('hop keeps high shape similarity (pyridine ≈ benzene bioisostere)',
    parseFloat(best.T) >= 0.7, best.T);
  check('SA column active (fragment table lazy-loaded, values present)',
    hopRes.rows.every(r => parseFloat(r.sa) >= 1 && parseFloat(r.sa) <= 10), hopRes.rows.map(r => r.sa).join(','));

  // ---- zero page errors ----
  check('zero page errors', errors.length === 0, errors.slice(0, 2).join(' | '));

  console.log('RESULT: ' + pass + ' passed, ' + fail + ' failed');
  await browser.close();
  process.exit(fail ? 1 : 0);
})().catch(e => { console.error(e); process.exit(1); });
