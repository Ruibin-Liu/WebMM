// m21_final_round.test.js — closing round: verify/refine as typed phase
// records (§9.5 v1.1), task-tray rerun from snapshots (structure-hash guard),
// quasi-harmonic G column, and the C3/mirror degeneracy in a real run.
const { chromium } = require('/opt/homebrew/lib/node_modules/playwright');
const EXE = '/Users/rliu/Library/Caches/ms-playwright/chromium_headless_shell-1243/chrome-headless-shell-mac-arm64/chrome-headless-shell';
const URL = 'http://localhost:8901/app/index.html';

let pass = 0, fail = 0;
function check(name, cond, detail) {
  if (cond) { pass++; console.log('  ✓ ' + name); }
  else { fail++; console.log('  ✗ ' + name + (detail ? ' — ' + detail : '')); }
}
const taskRows = page => page.evaluate(async () => {
  const ts = await wbTaskStore.ensure();
  if (!ts) return null;
  return (await ts.all('tasks')).map(r => ({ id: r.value.spec.taskId, kind: r.value.spec.kind,
    phase: r.value.spec.input.phase || null, state: r.value.state,
    mdFrames: r.value.mdFrames, refinedCount: r.value.refinedCount, jobs: r.value.jobs,
    input: r.value.spec.input }));
});

(async () => {
  const browser = await chromium.launch({ executablePath: EXE });
  const ctx = await browser.newContext({ viewport: { width: 1440, height: 900 } });
  const page = await ctx.newPage();
  const errors = [];
  page.on('pageerror', e => errors.push('PAGEERROR: ' + e.message.slice(0, 200)));

  await page.goto(URL, { waitUntil: 'load' });
  await page.waitForFunction(() => window.webmm !== undefined, null, { timeout: 120000 });

  const loadMol = async smi => {
    await page.evaluate(s => { document.getElementById('input').value = s; process(true); }, smi);
    await page.waitForFunction(() => !!document.querySelector('#svg svg'), null, { timeout: 60000 });
    await page.evaluate(() => embed3D());
    await page.waitForFunction(() => !!document.querySelector('#viewer3d canvas'), null, { timeout: 90000 });
  };

  await loadMol('CCO ethanol');
  await page.evaluate(() => {
    document.getElementById('imtdRounds').value = '2';
    document.getElementById('imtdSteps').value = '800';
    runImtd();
  });
  await page.waitForFunction(() => /unique minima/.test(document.getElementById('status3d').textContent), null, { timeout: 180000 });

  // ================= 1. verify + refine leave typed phase records =================
  await page.evaluate(() => window.runMdVerify());
  await page.waitForFunction(() => /MD 占据验证完成/.test(document.getElementById('status3d').textContent), null, { timeout: 180000 });
  await page.evaluate(() => window.runNebRefine());
  await page.waitForFunction(() => /势垒精修完成/.test(document.getElementById('status3d').textContent), null, { timeout: 300000 });
  let rows = await taskRows(page);
  const vrec = rows && rows.find(r => r.phase === 'verify' && r.state === 'complete');
  const rrec = rows && rows.find(r => r.phase === 'refine' && r.state === 'complete');
  check('verify phase record complete with mdFrames', !!vrec && vrec.mdFrames > 0, JSON.stringify(vrec && [vrec.state, vrec.mdFrames]));
  check('refine phase record complete with refinedCount/jobs', !!rrec && rrec.refinedCount >= 1 && rrec.jobs >= 1,
        JSON.stringify(rrec && [rrec.state, rrec.refinedCount, rrec.jobs]));

  // ================= 2. rerun from snapshot =================
  const trayHasRerun = await page.evaluate(() => {
    taskTrayRender();
    return [...document.querySelectorAll('#taskTrayList button')].some(b => /重跑/.test(b.textContent));
  });
  check('terminal tray rows expose 重跑', trayHasRerun);
  // rerun the search-phase imtd record: fields replay from the snapshot
  const reran = await page.evaluate(async () => {
    const ts = await wbTaskStore.ensure();
    const rows = (await ts.all('tasks')).map(r => r.value);
    const search = rows.filter(v => v.spec.kind === 'imtd' && !v.spec.input.phase && v.state === 'complete').pop();
    if (!search) return { ok: false, why: 'no search record' };
    // perturb the fields first — the rerun must override them from the snapshot
    document.getElementById('imtdRounds').value = '1';
    document.getElementById('imtdSteps').value = '200';
    await window.rerunWbTask(search.spec.taskId);
    return { ok: true };
  });
  check('rerun(search record) starts', reran.ok, JSON.stringify(reran));
  await page.waitForFunction(() => /unique minima/.test(document.getElementById('status3d').textContent), null, { timeout: 180000 });
  const fieldsRestored = await page.evaluate(() => ({
    rounds: document.getElementById('imtdRounds').value,
    steps: document.getElementById('imtdSteps').value,
  }));
  check('rerun restores the SNAPSHOT fields (not the perturbed ones)',
    fieldsRestored.rounds === '2' && fieldsRestored.steps === '800', JSON.stringify(fieldsRestored));
  rows = await taskRows(page);
  const rerunRec = rows && rows.filter(r => r.kind === 'imtd' && !r.phase).pop();
  check('rerun produces a fresh complete search record', !!rerunRec && rerunRec.state === 'complete',
        JSON.stringify(rerunRec && rerunRec.state));

  // phase records refuse standalone rerun
  const phaseRefuse = await page.evaluate(async () => {
    const ts = await wbTaskStore.ensure();
    const rows = (await ts.all('tasks')).map(r => r.value);
    const v = rows.filter(x => x.spec.input.phase === 'verify').pop();
    await window.rerunWbTask(v.spec.taskId);
    return document.getElementById('status3d').textContent;
  });
  check('phase record rerun refused with guidance', /阶段记录/.test(phaseRefuse), phaseRefuse.slice(0, 60));

  // structure change → hash guard refuses
  await page.evaluate(() => { document.getElementById('input').value = 'CCC propane'; process(true); });
  await page.waitForTimeout(900);
  const guardRefuse = await page.evaluate(async () => {
    const ts = await wbTaskStore.ensure();
    const rows = (await ts.all('tasks')).map(r => r.value);
    const s = rows.filter(v => v.spec.kind === 'conformers' || (v.spec.kind === 'imtd' && !v.spec.input.phase)).pop();
    await window.rerunWbTask(s.spec.taskId);
    return document.getElementById('status3d').textContent;
  });
  check('structure mismatch refuses rerun with both hashes', /重跑被拒/.test(guardRefuse) && /快照/.test(guardRefuse),
        guardRefuse.slice(0, 70));

  // ================= 3. QH + degeneracy visible in a fresh run =================
  await loadMol('CO methanol');
  await page.evaluate(() => {
    document.getElementById('imtdRounds').value = '2';
    document.getElementById('imtdSteps').value = '800';
    runImtd();
  });
  await page.waitForFunction(() => /unique minima/.test(document.getElementById('status3d').textContent), null, { timeout: 180000 });
  const meth = await page.evaluate(() => {
    const rows = [...document.querySelectorAll('#basinTableWrap tbody tr')].map(r => [...r.children].map(c => c.textContent));
    const title = document.querySelector('#basinTableWrap tbody span');
    return { rows, g: lastLandscape.basins.map(b => b.g), qhTitle: title ? title.title : '' };
  });
  check('methanol run: degeneracy recorded (rotor wells observed)',
    meth.g.length > 0 && meth.g.every(x => x >= 1) && Math.max(...meth.g) >= 2, JSON.stringify(meth.g));
  check('G column carries the quasi-harmonic tooltip with the plain value', /纯谐振/.test(meth.qhTitle), meth.qhTitle.slice(0, 60));

  check('zero page errors', errors.length === 0, errors.join(' | '));

  console.log(`\nRESULT: ${pass} passed, ${fail} failed`);
  await browser.close();
  process.exit(fail ? 1 : 0);
})().catch(e => { console.error('FATAL:', e); process.exit(1); });
