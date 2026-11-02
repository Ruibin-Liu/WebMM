// M15 acceptance: parallel multi-metad — the IMTD scatter-gather driver runs
// each round's short metadynamics runs across a conf.worker pool and produces
// the bit-identical ensemble of the sequential one-shot wasm call.
const { chromium } = require('/opt/homebrew/lib/node_modules/playwright');

const EXE = '/Users/rliu/Library/Caches/ms-playwright/chromium_headless_shell-1243/chrome-headless-shell-mac-arm64/chrome-headless-shell';

setTimeout(() => { console.log('WATCHDOG: test exceeded 420s'); process.exit(2); }, 420000);
(async () => {
    const browser = await chromium.launch({
    executablePath: EXE,
    args: ['--disable-background-timer-throttling', '--disable-renderer-backgrounding',
           '--disable-features=IntensiveWakeUpThrottling'],
  });
    const page = await browser.newPage({ viewport: { width: 1440, height: 900 } });
    const errors = [];
    page.on('console', m => { if (m.type() === 'error') errors.push(m.text()); });
    page.on('pageerror', e => errors.push('PAGEERROR: ' + e.message));

    await page.goto('http://localhost:8901/app/index.html', { waitUntil: 'load' });
    await page.waitForFunction(() => window.webmmReady === true, null, { timeout: 60000 });

    let pass = 0, fail = 0;
    const check = (name, cond, detail = '') => {
        if (cond) { pass++; console.log(`  ✓ ${name}${detail ? ' — ' + detail : ''}`); }
        else { fail++; console.log(`  ✗ ${name}${detail ? ' — ' + detail : ''}`); }
    };

    // ethanol (explicit H) — 2 rotatable torsions, fast quenches
    const out = await page.evaluate(async () => {
        const sdf = await (await fetch('/scripts/val_set/ethanol.sdf')).text();
        const opts = {
            max_rounds: 2, runs_per_round: 3, steps_per_run: 600,
            quench_interval: 200, quench_max_iterations: 200,
            crossing_pairs: 2, seed: 42, workers: 3,
        };
        const rounds = [];
        const par = await window.IMTD.runParallel(sdf, opts, r => rounds.push(r));
        const seq = window.webmm.run_imtd_conformer_search(sdf, window.IMTD.optsToWasm(opts));
        return {
            exposed: !!window.IMTD && typeof window.IMTD.runParallel === 'function',
            par: {
                nConfs: par.nConfs, energies: par.energies, converged: par.converged,
                roundsExecuted: par.roundsExecuted, newMinima: par.newMinimaPerRound,
                poolSize: par.poolSize, workers: par.workers, usedWorkers: par.usedWorkers,
                roundReports: rounds,
            },
            seq: { nConfs: seq.n_confs(), energies: seq.energies(), roundsExecuted: seq.rounds_executed(), converged: seq.converged() },
        };
    });

    check('window.IMTD exposed (driver + optsToWasm)', out.exposed);
    check('parallel search found H rotamers (anti/gauche)', out.par.nConfs >= 2, `${out.par.nConfs} conformers`);
    check('energies sorted ascending', out.par.energies.every((e, i) => i === 0 || out.par.energies[i - 1] <= e));
    check('rounds reported through onRound', out.par.roundReports.length === out.par.roundsExecuted,
          `${JSON.stringify(out.par.roundReports)} vs ${out.par.roundsExecuted} rounds`);
    check('3 runs per round scattered', out.par.roundReports.every(r => r.runs === 3), JSON.stringify(out.par.roundReports));
    check('pool reports 3 workers', out.par.workers === 3, `${out.par.workers}`);
    check('all 3 workers actually executed runs (real parallelism)', out.par.usedWorkers === 3, `${out.par.usedWorkers}/3 used`);
    check('pool size = 2 rotatable torsions (ethanol)', out.par.poolSize === 2, `${out.par.poolSize}`);

    // parity: parallel scatter-gather == sequential one-shot, bit-for-bit
    const sameN = out.par.nConfs === out.seq.nConfs;
    const sameE = sameN && out.par.energies.every((e, i) => e === out.seq.energies[i]);
    check('parallel == sequential (bit-identical energies)', sameE,
          `${out.par.nConfs} vs ${out.seq.nConfs} conformers`);
    check('rounds executed match', out.par.roundsExecuted === out.seq.roundsExecuted,
          `${out.par.roundsExecuted} vs ${out.seq.roundsExecuted}`);
    check('convergence flag match', out.par.converged === out.seq.converged);

    // error path: bad SDF rejects with a surfaced message (no hang)
    const bad = await page.evaluate(async () => {
        try { await window.IMTD.runParallel('not an sdf', { workers: 2 }); return 'no-throw'; }
        catch (e) { return String(e.message || e); }
    });
    check('bad SDF rejects with error message', bad !== 'no-throw' && /parse|error/i.test(bad), bad);

    check('zero page errors', errors.length === 0, errors.join(' | '));

    console.log(`\nRESULT: ${pass} passed, ${fail} failed`);
    await browser.close();
    process.exit(fail ? 1 : 0);
})().catch(e => { console.error('FATAL:', e); process.exit(1); });
