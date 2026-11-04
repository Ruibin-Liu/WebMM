#![allow(clippy::unwrap_used, reason = "diagnostic CLI on constant inputs")]
//! Climbing-image validation across multi-CV molecules: full native
//! iMTD → landscape → refine pipeline, per-edge NEB with climbing on.
use webmm::metad::imtd::{execute_run, Quencher};
use webmm::metad::imtd::{ImtdConfig, ImtdSearch, SearchEngine};
use webmm::neb::{neb, NebConfig};

fn run_molecule(name: &str, sdf: &str, cfg: &ImtdConfig) {
    let heavy = webmm::molecule::parser::parse_sdf(sdf).unwrap();
    let mol = webmm::molecule::hydrogens::add_hydrogens(&heavy);
    let engine = SearchEngine::MMFF94s;
    let quencher = Quencher::new(&mol, engine);
    let mut search = ImtdSearch::new(&mol, engine, cfg);
    let t0 = std::time::Instant::now();
    while !search.done() {
        let specs = search.plan().to_vec();
        if specs.is_empty() {
            break;
        }
        for spec in &specs {
            let minima = execute_run(&mol, engine, &quencher, cfg, spec);
            search.gather(spec.run_index, minima);
        }
    }
    let t_search = t0.elapsed();
    // MD verify ladder (feeds transitions into the landscape)
    let t1 = std::time::Instant::now();
    let specs = search.plan_verify(2);
    for spec in &specs {
        let minima = execute_run(&mol, engine, &quencher, cfg, spec);
        search.gather_verify(spec.run_index, minima);
    }
    let t_verify = t1.elapsed();
    let jobs = search.plan_refine();
    let n_jobs = jobs.len();
    let t2 = std::time::Instant::now();
    let ff = webmm::mmff::MMFFForceField::new(&mol, webmm::MMFFVariant::MMFF94s);
    let mut n_conv = 0usize;
    let mut n_climb = 0usize;
    let mut n_fallback = 0usize;
    let mut phase2_share = 0.0f64;
    let mut rows: Vec<(f64, f64, bool, usize, usize)> = Vec::new(); // (Ea, saddle, conv, iters, climb_from)
    let climbing_toggle = std::env::var_os("WEBMM_NO_CLIMB").is_none();
    for job in &jobs {
        let r = neb(
            &ff,
            &job.xyz_a,
            &job.xyz_b,
            &NebConfig {
                climbing: climbing_toggle,
                ..Default::default()
            },
        );
        if r.converged {
            n_conv += 1;
        }
        if r.climb_started > 0 {
            n_climb += 1;
        } else {
            n_fallback += 1;
        }
        if r.climb_started > 0 {
            phase2_share += (r.iterations - r.climb_started) as f64 / r.iterations as f64;
        }
        rows.push((
            r.energies[0],
            r.saddle_energy,
            r.converged,
            r.iterations,
            r.climb_started,
        ));
        search.gather_refine(job.index, Some(r.saddle_energy));
    }
    let t_refine = t2.elapsed();
    let n_above = rows.iter().filter(|r| r.1 > r.0).count();
    println!(
        "== {name}: {} atoms, {} minima, {} refine edges",
        mol.atoms.len(),
        search.result().conformers.len(),
        n_jobs
    );
    println!(
        "   search {:?} | verify {:?} | refine {:?} ({} edges)",
        t_search, t_verify, t_refine, n_jobs
    );
    println!(
        "   converged {n_conv}/{} | climb {n_climb} | no-climb(fallback/guard) {n_fallback} | phase2 share {:.0}% | saddle>Ea {n_above}/{}",
        n_jobs,
        if n_climb > 0 { phase2_share / n_climb as f64 * 100.0 } else { 0.0 },
        n_jobs
    );
    for (i, (ea, s, c, it, cf)) in rows.iter().enumerate() {
        println!(
        "   edge {i:>2}: Ea={ea:7.3}  saddle={s:7.3}  Δ={:7.3}  conv={c} iters={it:>5} climb_from={cf}",
        s - ea
    );
    }
}

fn main() {
    let cfg = ImtdConfig {
        max_rounds: 2,
        runs_per_round: 3,
        steps_per_run: 800,
        ..Default::default()
    };
    for (name, path) in [
        ("ibuprofen", "tests/fixtures/conformers/ibuprofen.sdf"),
        ("threonine", "tests/fixtures/conformers/threonine.sdf"),
        ("aspirin", "tests/fixtures/conformers/aspirin.sdf"),
        ("naphthalene", "tests/fixtures/conformers/naphthalene.sdf"),
    ] {
        let sdf = std::fs::read_to_string(path).unwrap();
        run_molecule(name, &sdf, &cfg);
    }
}
