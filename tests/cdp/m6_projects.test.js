// m6_projects.test.js — M4 P1 project management: registry + per-project
// DBs + zero-copy migration + URL-scoped active project + export v1
// envelope + import-as-new-project. Docs: project-management.md.
const { chromium } = require('/opt/homebrew/lib/node_modules/playwright');
const fs = require('fs');
const EXE = '/Users/rliu/Library/Caches/ms-playwright/chromium_headless_shell-1243/chrome-headless-shell-mac-arm64/chrome-headless-shell';
const URL = 'http://localhost:8901/app/platform.html';

let pass = 0, fail = 0;
function check(name, cond, detail) {
  if (cond) { pass++; console.log('  ✓ ' + name); }
  else { fail++; console.log('  ✗ ' + name + (detail ? ' — ' + detail : '')); }
}
// the store initializes lazily on entering search mode (same trigger the
// other suites use) — entering it first, then waiting for both signals
const ready = async page => {
  await page.waitForFunction(() => window.webmm !== undefined, null, { timeout: 60000 });
  await page.evaluate(() => { try { switchMode('search'); } catch (e) {} });
  await page.waitForFunction(() => window.__platformCompat, null, { timeout: 60000 });
};

(async () => {
  const browser = await chromium.launch({ executablePath: EXE });
  const ctx = await browser.newContext({ viewport: { width: 1440, height: 900 } });
  const page = await ctx.newPage();
  const errors = [];
  page.on('pageerror', e => errors.push('PAGEERROR: ' + e.message.slice(0, 200)));
  page.on('dialog', d => d.accept(d.type() === 'prompt' ? '测试项目B' : undefined));

  await page.goto(URL, { waitUntil: 'load' });
  await page.waitForFunction(() => typeof switchMode === 'function', null, { timeout: 60000 });

  // ---- 0. fresh page straight into the overview boots the store ----
  // (regression: the lazy init was only tripped by search mode; the project
  // overview used to show "注册库不可用" and 新建 alerted)
  await page.evaluate(() => switchMode('project'));
  await page.waitForFunction(() => (document.getElementById('projectList') || { textContent: '' }).textContent.includes('默认项目'), null, { timeout: 30000 });
  const freshBoot = await page.evaluate(() => ({ compat: !!window.__platformCompat, reg: !!window.__projectRegistry }));
  check('fresh page → overview directly: store + registry boot', freshBoot.compat && freshBoot.reg, JSON.stringify(freshBoot));

  // ---- 1. default boot: zero-copy migration + list rendering ----
  await ready(page);
  await page.evaluate(() => { switchMode('search'); loadDemoLibrary(); });
  await page.waitForFunction(() => (window.__search.getState() || []).length === 55, null, { timeout: 60000 });
  await page.evaluate(() => {    // a round so export carries facts
    document.getElementById('searchMode').value = 'sim'; onSearchModeChange();
    document.getElementById('searchQuery').value = 'CCO';
    runSearch();
  });
  await page.waitForFunction(() => /at T/.test(document.getElementById('searchResultStatus').textContent), null, { timeout: 120000 });
  await page.evaluate(() => switchMode('project'));
  await page.waitForFunction(() => (document.getElementById('projectList') || {}).children !== undefined && document.getElementById('projectList').children.length > 0, null, { timeout: 15000 });
  const list1 = await page.evaluate(() => [...document.getElementById('projectList').children].map(d => d.textContent));
  check('migration: default project listed as current', list1.length === 1 && /默认项目/.test(list1[0]) && /当前/.test(list1[0]), JSON.stringify(list1));

  // ---- 2. create a project: per-tab URL + storage isolation ----
  const nav = page.waitForNavigation({ waitUntil: 'load', timeout: 60000 });
  await page.evaluate(() => createProjectUI());
  await nav;
  await ready(page);
  check('create switches the URL to ?project=<id>', /project=/.test(page.url()), page.url());
  const iso = await page.evaluate(async () => {
    switchMode('search');
    await new Promise(r => setTimeout(r, 1500));   // restore attempt on an empty store
    return {
      lib: (window.__search.getState() || []).length,
      proj: (window.__activeProject || {}).name,
      id: (window.__activeProject || {}).id,
    };
  });
  check('new project starts empty (storage isolation)', iso.lib === 0 && iso.proj === '测试项目B', JSON.stringify(iso));

  // ---- 3. switch back to default: data intact ----
  const nav2 = page.waitForNavigation({ waitUntil: 'load', timeout: 60000 });
  await page.goto(URL, { waitUntil: 'load' });
  await nav2.catch(() => {});
  await ready(page);
  await page.evaluate(() => switchMode('search'));
  await page.waitForFunction(() => (window.__search.getState() || []).length === 55, null, { timeout: 30000 });
  check('default project restores its 55-molecule library', true);

  // ---- 4. export: v1 envelope ----
  const dl = page.waitForEvent('download', { timeout: 15000 });
  await page.evaluate(() => exportProjectFile());
  await (await dl).saveAs('/tmp/m4-export.json');
  const env = JSON.parse(fs.readFileSync('/tmp/m4-export.json', 'utf8'));
  check('export: v1 envelope (format/version/logHash)',
    env.format === 'webmm-project' && env.version === 1 && /^[0-9a-f]{16}$/.test(env.logHash) && env.sourceProjectId === 'default',
    'format=' + env.format + ' v' + env.version + ' hash=' + env.logHash);
  check('export: commands + layer-tagged facts included',
    env.commands.some(c => c.type === 'ImportLibrary') && env.facts.length >= 1 && env.facts.every(f => f.value.layer === 'stochastic'),
    'cmds=' + env.commands.length + ' facts=' + env.facts.length);

  // ---- 5. import the file: dedup prompt (same source+hash) then a NEW project ----
  await page.setInputFiles('#projectImportFile', '/tmp/m4-export.json');
  await page.waitForNavigation({ waitUntil: 'load', timeout: 60000 }).catch(() => {});
  await ready(page);
  const imp = await page.evaluate(async () => {
    switchMode('search');
    await new Promise(r => setTimeout(r, 1500));
    return { url: location.search, lib: (window.__search.getState() || []).length, proj: (window.__activeProject || {}).name, id: (window.__activeProject || {}).id };
  });
  check('import lands in a NEW project (copy, after dedup prompt)',
    /project=/.test(imp.url) && imp.id && imp.id !== 'default' && imp.proj === '默认项目' /* named after the file */,
    JSON.stringify(imp));
  check('imported project restores the library', imp.lib === 55, 'lib=' + imp.lib);
  const rounds = await page.evaluate(() => {
    const c = window.__platformCompat;
    return Object.keys(c.projectApi().state.rounds).length;
  });
  check('imported project restores the round (command log)', rounds >= 1, 'rounds=' + rounds);

  // ---- 6. registry consistency + no errors ----
  const regList = await page.evaluate(async () => (await window.__projectRegistry.list()).map(r => r.name + '/' + r.state));
  check('registry: three ready projects', regList.length === 3 && regList.every(x => x.endsWith('/ready')), JSON.stringify(regList));
  check('no page errors', errors.length === 0, errors.slice(0, 2).join(' | '));

  await browser.close();
  console.log('RESULT: ' + pass + ' passed, ' + fail + ' failed');
  process.exit(fail ? 1 : 0);
})().catch(e => { console.error(e); process.exit(1); });
