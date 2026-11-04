// Diagnostic perf probe (freeze round). Inputs are compile-time constants.
#![allow(
    clippy::unwrap_used,
    reason = "diagnostic CLI on constant inputs — house examples convention"
)]
// Measures the universal kernel
// (MMFF gradient eval), the optimizer path, NEB per-edge, iMTD search and
// the Hessian tier. Run with `cargo run --release --example perf_probe`.
use std::time::Instant;
use webmm::forces::ForceField;

fn timeit(label: &str, n: usize, mut f: impl FnMut()) {
    f(); // warmup
    let t0 = Instant::now();
    for _ in 0..n {
        f();
    }
    let us = t0.elapsed().as_micros() as f64 / n as f64;
    println!("{label:<42} {us:>10.1} µs/call");
}

const ETHANOL: &str = "ethanol\n     RDKit          3D\n\n  9  8  0  0  0  0  0  0  0  0999 V2000\n    0.8867    0.1097    0.1398 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.4649   -0.1632   -0.4871 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.4908    0.4341    0.2924 O   0  0  0  0  0  0  0  0  0  0  0  0\n    1.0600    1.1872    0.2278 H   0  0  0  0  0  0  0  0  0  0  0  0\n    1.6897   -0.3301   -0.4581 H   0  0  0  0  0  0  0  0  0  0  0  0\n    0.9343   -0.3042    1.1524 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.5119    0.2599   -1.4947 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.6562   -1.2384   -0.5550 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.4469    0.0448    1.1825 H   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0\n  2  3  1  0\n  1  4  1  0\n  1  5  1  0\n  1  6  1  0\n  2  7  1  0\n  2  8  1  0\n  3  9  1  0\nM  END\n";

fn main() {
    // ── 1. the universal kernel: MMFF energy_and_gradient ──
    let mol = webmm::molecule::parser::parse_sdf(ETHANOL).unwrap();
    let ff = webmm::mmff::MMFFForceField::new(&mol, webmm::MMFFVariant::MMFF94s);
    let coords: Vec<[f64; 3]> = mol.atoms.iter().map(|a| a.position).collect();
    let mut grad = vec![[0.0f64; 3]; coords.len()];
    timeit("MMFF energy_and_gradient (ethanol 9at)", 20000, || {
        let _ = ff.energy_and_gradient(&coords, &mut grad);
    });

    // ── 2. one quench from a displaced geometry (300-iter cap) ──
    let mut disp = coords.clone();
    for (i, p) in disp.iter_mut().enumerate() {
        p[0] += 0.02 * (i as f64 % 3.0) - 0.02;
        p[1] += 0.015 * (i as f64 % 2.0) - 0.015;
    }
    let conv = webmm::ConvergenceOptions {
        max_iterations: 300,
        ..Default::default()
    };
    timeit("L-BFGS quench 300-cap (displaced)", 300, || {
        let _ = webmm::optimizer::optimize(&ff, &disp, &conv);
    });

    // ── 3. NEB single edge: butane anti↔gauche, default config ──
    let bsdf = "Butane\n     RDKit          3D\n\n  4  3  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.5500    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    3.0500    0.1000    0.9000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    4.5500    0.1000    0.9000 C   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0\n  2  3  1  0\n  3  4  1  0\nM  END\n";
    let heavy = webmm::molecule::parser::parse_sdf(bsdf).unwrap();
    let bmol = webmm::molecule::hydrogens::add_hydrogens(&heavy);
    let bff = webmm::mmff::MMFFForceField::new(&bmol, webmm::MMFFVariant::MMFF94s);
    let base: Vec<[f64; 3]> = bmol.atoms.iter().map(|a| a.position).collect();
    let qconv = webmm::ConvergenceOptions {
        max_iterations: 500,
        ..Default::default()
    };
    let ra = webmm::optimizer::optimize(&bff, &base, &qconv);
    let (j, k) = (1usize, 2usize);
    let a0 = ra.optimized_coords;
    let axis = [
        a0[k][0] - a0[j][0],
        a0[k][1] - a0[j][1],
        a0[k][2] - a0[j][2],
    ];
    let nn = (axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2]).sqrt();
    let u = [axis[0] / nn, axis[1] / nn, axis[2] / nn];
    let th = 2.0 * std::f64::consts::PI / 3.0;
    let (sn, cs) = (th.sin(), th.cos());
    let mut seen = vec![false; bmol.atoms.len()];
    seen[j] = true;
    seen[k] = true;
    let mut side = vec![k];
    let mut queue = std::collections::VecDeque::new();
    queue.push_back(k);
    while let Some(z) = queue.pop_front() {
        for &w in &bmol.adjacency[z] {
            if !seen[w] {
                seen[w] = true;
                side.push(w);
                queue.push_back(w);
            }
        }
    }
    let mut b0 = a0.clone();
    for &i in &side {
        let d = [
            a0[i][0] - a0[k][0],
            a0[i][1] - a0[k][1],
            a0[i][2] - a0[k][2],
        ];
        let cr = [
            u[1] * d[2] - u[2] * d[1],
            u[2] * d[0] - u[0] * d[2],
            u[0] * d[1] - u[1] * d[0],
        ];
        let dt = u[0] * d[0] + u[1] * d[1] + u[2] * d[2];
        for x in 0..3 {
            b0[i][x] = a0[k][x] + cs * d[x] + sn * cr[x] + u[x] * dt * (1.0 - cs);
        }
    }
    let rb = webmm::optimizer::optimize(&bff, &b0, &qconv);
    let neb_cfg = webmm::neb::NebConfig::default();
    let r = webmm::neb::neb(&bff, &a0, &rb.optimized_coords, &neb_cfg); // warmup
    let t0 = Instant::now();
    let r2 = webmm::neb::neb(&bff, &a0, &rb.optimized_coords, &neb_cfg);
    let neb_ms = t0.elapsed().as_micros() as f64 / 1000.0;
    println!(
        "{:<42} {:>10.1} ms  (iters={}, conv={})",
        "NEB edge butane (default cfg)", neb_ms, r2.iterations, r2.converged
    );
    let _ = r;

    // ── 4. iMTD ethanol search (2×3×800) ──
    let cfg = webmm::metad::imtd::ImtdConfig {
        max_rounds: 2,
        runs_per_round: 3,
        steps_per_run: 800,
        quench_interval: 200,
        quench_max_iterations: 200,
        crossing_pairs: 2,
        seed: 42,
        ..Default::default()
    };
    let t0 = Instant::now();
    let res = webmm::metad::imtd::imtd_conformer_search(
        &mol,
        webmm::metad::imtd::SearchEngine::MMFF94s,
        &cfg,
    );
    println!(
        "{:<42} {:>10.1} ms  (n={})",
        "iMTD search ethanol 2×3×800",
        t0.elapsed().as_micros() as f64 / 1000.0,
        res.conformers.len()
    );

    // ── 5. Hessian tier (ethanol, 27 coords → 54 grad evals) ──
    let masses: Vec<f64> = mol
        .atoms
        .iter()
        .map(|a| webmm::md::atomic_mass(a.atomic_number))
        .collect();
    let rc = webmm::optimizer::optimize(&ff, &coords, &qconv);
    let qc = rc.optimized_coords;
    timeit("vib::frequencies (ethanol)", 50, || {
        let _ = webmm::vib::frequencies(&ff, &qc, &masses);
    });
}
