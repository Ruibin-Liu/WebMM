// m7_activity.test.js — M5 milestone (docs/assay-foundation.md; suite number
// m7 = next in the CDP sequence): measured-activity CSV import (two-step
// preview with the match ladder), the activity result column, identity audit
// duplicate families (full vs parent grouping), undo, persistence.
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

  await page.goto(URL + '?m=search', { waitUntil: 'load' });
  await page.waitForFunction(() => document.getElementById('rdkitVersion').textContent !== 'Loading...', null, { timeout: 60000 });
  await page.evaluate(() => { indexedDB.deleteDatabase('webmm-platform'); indexedDB.deleteDatabase('webmm-projects'); });
  await page.reload({ waitUntil: 'load' });
  await page.waitForFunction(() => document.getElementById('rdkitVersion').textContent !== 'Loading...', null, { timeout: 60000 });

  // seed: salt family (propanol vs propanol·HCl) + stereo-specified/unspecified pair
  await page.evaluate(async () => {
    await new Promise(r => loadSearchLibraryFrom([
      { smiles: 'CCO', name: 'ethanol' },
      { smiles: 'CCCO.Cl', name: 'propanol HCl' },
      { smiles: 'CCCO', name: 'propanol' },
      { smiles: 'C/C=C/C', name: 'trans-but2ene' },
    ]).persisted.then(r));
  });
  await page.waitForTimeout(500);

  // ---- 1. import preview: dry-run report, ghost row excluded, tiers recorded ----
  await page.evaluate(() => activityImportUI());
  await page.fill('#actImpName', 'IC50');
  await page.fill('#actImpUnit', 'µM');
  await page.fill('#actImpSource', 'm7-test');
  await page.fill('#actImpCsv', 'name,value,note\nethanol,8.2,\npropanol,>10,batch2\nghost,5.0,\n');
  await page.evaluate(() => activityImportPreview());
  const report = await page.textContent('#actImpReport');
  check('preview: 3 parsed, 2 matched, 1 unresolved (ghost NOT silently dropped)',
    /解析 3 行/.test(report) && /匹配 2/.test(report) && /未解析 1/.test(report) && /ghost/.test(report), report.replace(/\s+/g, ' ').slice(0, 120));
  check('preview: loose matches warned (name tier)', /松匹配 2/.test(report));
  await page.evaluate(() => activityImportCommit());
  await page.waitForTimeout(600);

  const panel = await page.evaluate(() => ({
    meta: document.getElementById('activityMeta').textContent,
    actColShown: document.getElementById('actCol').style.display !== 'none',
  }));
  check('panel: assay meta lands (name/unit/source/2 molecules)',
    /IC50/.test(panel.meta) && /µM/.test(panel.meta) && /2 分子/.test(panel.meta), panel.meta);
  check('result table: activity column visible', panel.actColShown);

  // ---- 2. activity column renders qualifiers; missing rows show — ----
  await page.evaluate(() => { document.getElementById('searchQuery').value = 'CCO'; runSearch(); });
  await page.waitForFunction(() => document.querySelectorAll('#searchRows tr').length >= 3, null, { timeout: 30000 });
  await page.waitForTimeout(300);
  const cells = await page.evaluate(() => [...document.querySelectorAll('#searchRows tr')].map(tr => tr.querySelectorAll('td')[4] && tr.querySelectorAll('td')[4].textContent.trim()));
  check('column: qualifier preserved (8.2 / >10), missing rows dash',
    cells.includes('8.2') && cells.includes('>10') && cells.includes('—'), JSON.stringify(cells));
  const tip = await page.evaluate(() => {
    const tr = [...document.querySelectorAll('#searchRows tr')].find(t => { const td = t.querySelectorAll('td')[4]; return td && td.textContent.trim() === '>10'; });
    return tr ? tr.querySelectorAll('td')[4].title : '';
  });
  check('column: matchTier in tooltip (loose match is auditable)', /match:name/.test(tip), tip);

  // ---- 3. identity audit: full vs parent grouping ----
  await page.evaluate(() => renderIdentityAudit());
  const fullAudit = await page.evaluate(() => ({
    flags: document.getElementById('auditFlags').textContent,
    noFamilies: /无重复家族/.test(document.getElementById('auditFamilies').textContent),
  }));
  check('audit(full): salt 1, stereo-unspecified 3, no full-key dupes',
    /盐\/多组分 1/.test(fullAudit.flags) && /立体未指定 3/.test(fullAudit.flags) && fullAudit.noFamilies, fullAudit.flags);
  await page.selectOption('#auditGroupBy', 'parent');
  const parentAudit = await page.evaluate(() => ({
    families: document.querySelectorAll('#auditFamilies > div').length,
    txt: document.getElementById('auditFamilies').textContent,
  }));
  check('audit(parent): salt family surfaces (propanol ⇄ propanol·HCl)',
    parentAudit.families >= 1 && /质子化|同物|成员/.test(parentAudit.txt), parentAudit.txt.replace(/\s+/g, ' ').slice(0, 80));
  check('audit: InChI tautomer scope disclaimer present',
    /互变/.test(await page.textContent('#identityAuditPanel')));

  // ---- 4. undo removes the assay (search rounds push SetThreshold etc. — undo past them) ----
  for (let i = 0; i < 4; i++) {
    const gone = await page.evaluate(() => document.getElementById('activityAssaySel').options.length === 1 && !document.getElementById('activityAssaySel').value);
    if (gone) break;
    await page.evaluate(() => triageUndo());
    await page.waitForTimeout(300);
  }
  check('undo removes the imported assay', await page.evaluate(() => document.getElementById('activityAssaySel').options.length === 1 && !document.getElementById('activityAssaySel').value));
  await page.evaluate(() => activityImportUI());
  await page.fill('#actImpName', 'IC50');
  await page.fill('#actImpCsv', 'name,value\nethanol,9.9\n');
  await page.evaluate(() => activityImportPreview());
  await page.evaluate(() => activityImportCommit());
  await page.waitForTimeout(500);
  const afterReplace = await page.evaluate(() => projStateM5().activities && Object.entries(projStateM5().activities).map(([k, a]) => a.meta.name + ':' + a.records.size).join('|'));
  check('same-name re-import replaces (one assay, new value)', /^IC50:1$/.test(afterReplace), afterReplace);

  // ---- 5. persistence: reload keeps the assay + audit ----
  await page.reload({ waitUntil: 'load' });
  await page.waitForFunction(() => document.getElementById('rdkitVersion').textContent !== 'Loading...', null, { timeout: 60000 });
  await page.evaluate(() => switchMode('search'));   // compat boots on entering search (established m6 pattern)
  await page.waitForFunction(() => typeof projStateM5 === 'function' && projStateM5() && Object.keys(projStateM5().activities).length > 0, null, { timeout: 30000 });
  const persisted = await page.evaluate(() => {
    const acts = projStateM5().activities;
    const a = Object.values(acts)[0];
    return { n: Object.keys(acts).length, name: a.meta.name, recs: a.records.size, val: [...a.records.values()][0][0].value };
  });
  check('reload: assay replays from the command log', persisted.n === 1 && persisted.name === 'IC50' && persisted.recs === 1 && persisted.val === 9.9, JSON.stringify(persisted));

  // ---- 6. export CSV carries the molId column (re-import ladder priority 1) ----
  const dl = page.waitForEvent('download', { timeout: 15000 }).catch(() => null);
  await page.evaluate(() => activityExportCsv());
  const d = await dl;
  check('export: activity-IC50.csv downloads', d && /activity-IC50\.csv/.test(d.suggestedFilename()), d && d.suggestedFilename());

  check('zero page errors', errors.length === 0, errors.slice(0, 3).join(' | '));
  await browser.close();
  console.log('RESULT: ' + pass + ' passed, ' + fail + ' failed');
  process.exit(fail ? 1 : 0);
})().catch(e => { console.error('ERR', e); process.exit(1); });
