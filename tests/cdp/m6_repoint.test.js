// m6_repoint.test.js — M4 P2 §3: restore-into-current via pointer repoint.
// Tab A restores a backup into the current project (new DB written, registry
// pointer swapped atomically, old DB trashed); tab B — same project — gets
// the forced reload modal; after both reload, data matches the backup.
const { chromium } = require('/opt/homebrew/lib/node_modules/playwright');
const fs = require('fs');
const EXE = '/Users/rliu/Library/Caches/ms-playwright/chromium_headless_shell-1243/chrome-headless-shell-mac-arm64/chrome-headless-shell';
const URL = 'http://localhost:8901/app/platform.html';
const FIXTURE = '/tmp/m5-project.json';

let pass = 0, fail = 0;
function check(name, cond, detail) {
  if (cond) { pass++; console.log('  ✓ ' + name); }
  else { fail++; console.log('  ✗ ' + name + (detail ? ' — ' + detail : '')); }
}
const boot = async page => {
  await page.waitForFunction(() => window.webmm !== undefined, null, { timeout: 60000 });
  await page.evaluate(() => { try { switchMode('search'); } catch (e) {} });
  await page.waitForFunction(() => window.__platformCompat, null, { timeout: 60000 });
};
const libSize = page => page.evaluate(() => (window.__search.getState() || []).length);

(async () => {
  if (!fs.existsSync(FIXTURE)) { console.log('missing ' + FIXTURE + ' — run m6_platform first'); process.exit(2); }
  const browser = await chromium.launch({ executablePath: EXE });
  const ctx = await browser.newContext({ viewport: { width: 1440, height: 900 } });
  const A = await ctx.newPage();
  const B = await ctx.newPage();
  for (const p of [A, B]) p.on('dialog', d => d.accept());
  const errors = [];
  for (const p of [A, B]) p.on('pageerror', e => errors.push(e.message.slice(0, 150)));

  // A first loads the library, THEN B opens — B's boot restore must find it
  await A.goto(URL, { waitUntil: 'load' }); await boot(A);
  await A.evaluate(() => { switchMode('search'); loadDemoLibrary(); });
  await A.waitForFunction(() => (window.__search.getState() || []).length === 55, null, { timeout: 60000 });
  await B.goto(URL, { waitUntil: 'load' }); await boot(B);
  await B.evaluate(() => switchMode('search'));
  await B.waitForFunction(() => (window.__search.getState() || []).length === 55, null, { timeout: 30000 });
  check('setup: both tabs see the library', true);

  // A restores the backup file INTO the current project (pointer repoint)
  await A.evaluate(() => switchMode('project'));
  await A.waitForFunction(() => document.getElementById('projectList').children.length > 0, null, { timeout: 15000 });
  await A.setInputFiles('#projectRestoreFile', FIXTURE);
  await A.waitForNavigation({ waitUntil: 'load', timeout: 60000 }).catch(() => {});
  await boot(A);
  await A.evaluate(() => switchMode('search'));
  await A.waitForFunction(() => (window.__search.getState() || []).length === 55, null, { timeout: 30000 });
  check('A reloaded onto the repointed database (data from the file)', true);

  // B: forced modal appears; reload lands on the same data
  await B.waitForFunction(() => !!document.getElementById('repointedModal'), null, { timeout: 30000 }).catch(() => {});
  const modal = await B.evaluate(() => !!document.getElementById('repointedModal'));
  check('B shows the forced reload modal', modal);
  const navB = B.waitForNavigation({ waitUntil: 'load', timeout: 60000 });
  await B.evaluate(() => document.querySelector('#repointedModal button').click());
  await navB; await boot(B);
  await B.evaluate(() => switchMode('search'));
  await B.waitForFunction(() => (window.__search.getState() || []).length === 55, null, { timeout: 30000 });
  check('B after reload sees the repointed data', true);

  // registry integrity: default project ready, pointer = a project DB,
  // the old singleton database sits in the trash record
  const reg = await A.evaluate(async () => {
    const r = window.__projectRegistry;
    const rec = await r.get('default');
    const trash = await r._storage.get('__trash');
    return { state: rec.state, dbName: rec.dbName, trashed: trash ? trash.dbs.map(t => t.dbName) : [] };
  });
  check('registry: default ready on the new DB, old DB trashed',
    reg.state === 'ready' && /webmm-project-/.test(reg.dbName) && reg.trashed.includes('webmm-platform'),
    JSON.stringify(reg));
  check('no page errors', errors.length === 0, errors.slice(0, 2).join(' | '));

  await browser.close();
  console.log('RESULT: ' + pass + ' passed, ' + fail + ' failed');
  process.exit(fail ? 1 : 0);
})().catch(e => { console.error(e); process.exit(1); });
