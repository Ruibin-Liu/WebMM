// M0 baseline sessions — recorded on the CURRENT app (pre-platform UI).
// Re-runnable: node scripts/baseline/m0_baseline.js  (needs :8901 serving repo root)
// Purpose: falsifiable before/after metrics for the platform redesign (v0.5 §8).
// Each session counts USER GESTURES (clicks + input edits + select changes)
// and wall time; "time-to-shortlist" = from first search action to a
// filtered, ordered, exported-ready list of <=10 candidates.
const { chromium } = require('/opt/homebrew/lib/node_modules/playwright');
const EXE = '/Users/rliu/Library/Caches/ms-playwright/chromium_headless_shell-1243/chrome-headless-shell-mac-arm64/chrome-headless-shell';
const URL = 'http://localhost:8901/app/index.html';

const r = (fn) => ({ t: Date.now(), n: 1 }); // helper marker (unused)

async function newPage(b) {
  const pg = await b.newPage();
  await pg.goto(URL, { waitUntil: 'load' });
  await pg.waitForFunction(() => document.getElementById('rdkitVersion').textContent !== 'Loading...', null, { timeout: 60000 });
  await pg.waitForFunction(() => window.webmm !== undefined, null, { timeout: 60000 });
  return pg;
}

(async () => {
  const b = await chromium.launch({ executablePath: EXE });
  const out = {};

  // ---------- Session A: hit discovery (query -> escalate -> shortlist) ----------
  {
    const pg = await newPage(b);
    let actions = 0; const t0 = Date.now(); let tShort = null;
    const act = () => actions++;
    await pg.evaluate(() => { switchMode('search'); }); act();                       // 1 tab switch
    await pg.evaluate(() => loadDemoLibrary()); act();                               // 2 load library
    await pg.waitForFunction(() => document.getElementById('searchStatus').textContent.includes('55'), null, { timeout: 60000 });
    // similarity pass
    await pg.evaluate(() => {
      document.getElementById('searchMode').value = 'sim'; onSearchModeChange();
      document.getElementById('searchQuery').value = 'CC(=O)Oc1ccccc1C(=O)O'; runSearch();
    }); act(); act(); act();                                                        // 3-5 mode+query+search
    await pg.waitForFunction(() => document.getElementById('searchResultStatus').textContent.includes('at T'), null, { timeout: 120000 });
    // tighten threshold
    await pg.evaluate(() => { document.getElementById('searchThreshold').value = '0.5'; runSearch(); }); act(); act();  // 6-7
    await pg.waitForFunction(() => document.getElementById('searchResultStatus').textContent.includes('at T'), null, { timeout: 120000 });
    // escalate to shape 3D
    await pg.evaluate(() => {
      document.getElementById('searchMode').value = 'shape'; onSearchModeChange(); runSearch();
    }); act(); act();                                                                // 8-9
    await pg.waitForFunction(() => /sorted by combo/.test(document.getElementById('searchResultStatus').textContent), null, { timeout: 300000 });
    // pharmacophore filter -> shortlist
    await pg.evaluate(() => { document.getElementById('pharmFilter').value = '0.6'; runSearch(); }); act(); act();      // 10-11
    await pg.waitForFunction(() => /sorted by combo/.test(document.getElementById('searchResultStatus').textContent), null, { timeout: 300000 });
    tShort = Date.now() - t0;
    const shortlist = await pg.evaluate(() => [...document.querySelectorAll('#searchRows tr')].map(tr => tr.children[1].textContent.trim()).slice(0, 10));
    out.A_hit_discovery = { actions, wall_s: +((Date.now() - t0) / 1000).toFixed(1), time_to_shortlist_s: +(tShort / 1000).toFixed(1), shortlist, shortlist_n: shortlist.length };
    await pg.close();
  }

  // ---------- Session B: batch triage + SAR ----------
  {
    const pg = await newPage(b);
    let actions = 0; const t0 = Date.now(); let tShort = null;
    const act = () => actions++;
    const BATCH = 'Cn1c(=O)c2c(ncn2C)n(C)c1=O caffeine\nCC(=O)Oc1ccccc1C(=O)O aspirin\nCC(C)Cc1ccc(C(C)C(=O)O)cc1 ibuprofen\nCOc1ccc2cc(C(C)C(=O)O)ccc2c1 naproxen\nCC(=O)Nc1ccc(O)cc1 paracetamol\nCN(C)C(=N)N=C(N)N metformin\nCC(C)CCCC(C)C1CCC2(C)C1CCC1C2CC=C2CC(O)CCC21C cholesterol\nCCCCCCCCCC decane\nOc1ccccc1C(=O)O salicylic\nCCOC(=O)c1ccc(N)cc1 benzocaine';
    await pg.evaluate((s) => { switchMode('batch'); document.getElementById('input').value = s; runBatch(); }, BATCH); act(); act(); act(); // 1-3
    await pg.waitForFunction(() => /processed|Cancelled/.test(document.getElementById('batchStatus').textContent), null, { timeout: 300000 });
    // triage: sort by QED desc
    await pg.evaluate(() => { const h = [...document.querySelectorAll('#batchTable th')].find(th => th.textContent.trim() === 'QED'); h.click(); }); act(); // 4
    // filter to drug-like
    await pg.evaluate(() => { document.getElementById('batchFilter').value = 'lipinski'; }); act();  // 5
    const rows = await pg.evaluate(() => [...document.querySelectorAll('#batchRows tr')].map(tr => tr.children[1].textContent.trim()).slice(0, 10));
    tShort = Date.now() - t0;
    // inspect the top candidate
    await pg.evaluate(() => document.querySelector('#batchRows tr').dispatchEvent(new MouseEvent('click', { bubbles: true }))); act(); // 6
    await pg.waitForTimeout(400);
    // SAR: scaffold frequencies over the demo library
    await pg.evaluate(() => { switchMode('search'); loadDemoLibrary(); }); act(); act(); // 7-8
    await pg.waitForFunction(() => document.getElementById('searchStatus').textContent.includes('55'), null, { timeout: 60000 });
    const tSar = Date.now();
    out.B_triage_sar = {
      actions, wall_s: +((Date.now() - t0) / 1000).toFixed(1), time_to_shortlist_s: +(tShort / 1000).toFixed(1),
      shortlist_n: rows.length, top_rows: rows.slice(0, 5),
    };
    await pg.close();
  }

  // ---------- Session C: lead hopping (search -> hit as new query -> compare) ----------
  {
    const pg = await newPage(b);
    let actions = 0; const t0 = Date.now(); let tShort = null;
    const act = () => actions++;
    await pg.evaluate(() => { switchMode('search'); loadDemoLibrary(); }); act(); act();
    await pg.waitForFunction(() => document.getElementById('searchStatus').textContent.includes('55'), null, { timeout: 60000 });
    // round 1: shape search with ibuprofen
    await pg.evaluate(() => {
      document.getElementById('searchMode').value = 'shape'; onSearchModeChange();
      document.getElementById('searchQuery').value = 'CC(C)Cc1ccc(C(C)C(=O)O)cc1'; runSearch();
    }); act(); act(); act();
    await pg.waitForFunction(() => /sorted by combo/.test(document.getElementById('searchResultStatus').textContent), null, { timeout: 300000 });
    const r1 = await pg.evaluate(() => [...document.querySelectorAll('#searchRows tr')].slice(0, 3).map(tr => tr.children[1].textContent.trim()));
    // pick hit #2 (not the self-match) -> load into single view -> use as new query
    await pg.evaluate(() => document.querySelectorAll('#searchRows tr')[1].dispatchEvent(new MouseEvent('click', { bubbles: true }))); act();
    await pg.waitForTimeout(600);
    await pg.evaluate(() => { switchMode('search'); document.getElementById('searchMode').value = 'shape'; onSearchModeChange(); }); act(); act();
    await pg.evaluate(() => document.getElementById('searchQuery').value = ''); act();
    await pg.evaluate(() => useCurrentAsQuery ? useCurrentAsQuery() : null).catch(() => {});
    // fall back to typing the SMILES (as a real user would if no shortcut existed)
    const hit2 = r1[1] === 'ibuprofen' ? (r1[2] || r1[1]) : r1[1];
    const smiOf = { 'naproxen': 'COc1ccc2cc(C(C)C(=O)O)ccc2c1', 'ibuprofen': 'CC(C)Cc1ccc(C(C)C(=O)O)cc1', 'ketoprofen': 'CC(C)c1cccc(c1)C(=O)O', 'flurbiprofen': 'CC(C)c1ccc(c(c1)c1ccccc1)C(=O)O', 'aspirin': 'CC(=O)Oc1ccccc1C(=O)O' };
    const smi = smiOf[hit2] || 'COc1ccc2cc(C(C)C(=O)O)ccc2c1';
    await pg.evaluate((s) => { document.getElementById('searchQuery').value = s; runSearch(); }, smi); act(); act();
    await pg.waitForFunction(() => /sorted by combo/.test(document.getElementById('searchResultStatus').textContent), null, { timeout: 300000 });
    const r2 = await pg.evaluate(() => [...document.querySelectorAll('#searchRows tr')].slice(0, 5).map(tr => tr.children[1].textContent.trim()));
    tShort = Date.now() - t0;
    out.C_lead_hopping = {
      actions, wall_s: +((Date.now() - t0) / 1000).toFixed(1), time_to_shortlist_s: +(tShort / 1000).toFixed(1),
      round1_top3: r1, round2_top5: r2, overlap_top5: r2.filter(x => r1.includes(x)).length,
      hop_via: hit2,
    };
    await pg.close();
  }

  // ---------- Session C': lead hopping AFTER M1c (hit-as-query) ----------
  {
    const pg = await newPage(b);
    let actions = 0; const t0 = Date.now(); let tShort = null;
    const act = () => actions++;
    await pg.evaluate(() => { switchMode('search'); loadDemoLibrary(); }); act(); act();
    await pg.waitForFunction(() => document.getElementById('searchStatus').textContent.includes('55'), null, { timeout: 60000 });
    await pg.evaluate(() => {
      document.getElementById('searchMode').value = 'shape'; onSearchModeChange();
      document.getElementById('searchQuery').value = 'CC(C)Cc1ccc(C(C)C(=O)O)cc1'; runSearch();
    }); act(); act(); act();
    await pg.waitForFunction(() => /sorted by combo/.test(document.getElementById('searchResultStatus').textContent), null, { timeout: 300000 });
    const r1 = await pg.evaluate(() => [...document.querySelectorAll('#searchRows tr')].slice(0, 3).map(tr => tr.children[1].textContent.trim()));
    // ONE-ACTION lead hop via ⇄
    await pg.evaluate(() => document.querySelectorAll('#searchRows tr')[1].querySelector('span[onclick*=hitAsQuery]').click()); act();
    await pg.waitForFunction(() => {
      const es = [...document.querySelectorAll('#queryHistoryPanel div[onclick]')];
      return es.length >= 2 && es[es.length - 1].textContent.includes('↳');
    }, null, { timeout: 300000 });
    const r2 = await pg.evaluate(() => [...document.querySelectorAll('#searchRows tr')].slice(0, 5).map(tr => tr.children[1].textContent.trim()));
    tShort = Date.now() - t0;
    out.C_prime_hit_as_query = {
      actions, wall_s: +((Date.now() - t0) / 1000).toFixed(1), time_to_shortlist_s: +(tShort / 1000).toFixed(1),
      round1_top3: r1, round2_top5: r2, overlap_top5: r2.filter(x => r1.includes(x)).length,
      hop_actions: 1, hop_via: r1[1],
    };
    await pg.close();
  }

  await b.close();
  console.log(JSON.stringify(out, null, 1));
})().catch(e => { console.error('BASELINE ERR', e.message); process.exit(1); });
