//! Benchmark WebMM MMFF94s energy+gradient throughput over an SDF directory.
//! Run: cargo run --release --example bench_mmff -- <SDF_DIR> [--min-ms 30]
//! Emits JSON: {name: {n_atoms, build_us, eg_us, eg_ops, us_per_op, us_per_atom}}
//!
//! Protocol (mirrored by the RDKit side in scripts/benchmark_mmff.py):
//!   * setup: parse SDF + MMFFForceField::new (timed separately, excluded from
//!     steady-state numbers)
//!   * warmup: 20 ops
//!   * steady state: calculate_energy_and_gradient until >= min-ms elapsed
use std::hint::black_box;
use std::io::Write;
use std::time::{Duration, Instant};

use webmm::mmff::MMFFForceField;
use webmm::molecule::parser::parse_sdf;
use webmm::MMFFVariant;

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let mut dir = "scripts/val_set".to_string();
    let mut min_ms = 30u64;
    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "--min-ms" => {
                i += 1;
                min_ms = args.get(i).and_then(|s| s.parse().ok()).unwrap_or(30);
            }
            s if !s.starts_with("--") => dir = s.to_string(),
            _ => {}
        }
        i += 1;
    }
    let min_dur = Duration::from_millis(min_ms);

    let mut files: Vec<String> = std::fs::read_dir(&dir)
        .expect("SDF dir")
        .flatten()
        .filter(|e| e.path().extension().and_then(|s| s.to_str()) == Some("sdf"))
        .map(|e| e.path().file_stem().unwrap().to_string_lossy().into_owned())
        .collect();
    files.sort();

    let mut out = String::from("{");
    let mut first = true;
    for name in &files {
        if !first {
            out.push(',');
        }
        first = false;
        let path = format!("{dir}/{name}.sdf");
        let mol = match std::fs::read_to_string(&path)
            .ok()
            .and_then(|t| parse_sdf(&t).ok())
        {
            Some(m) => m,
            None => {
                out.push_str(&format!("\"{name}\":{{\"error\":\"parse\"}}"));
                continue;
            }
        };
        let coords: Vec<[f64; 3]> = mol.atoms.iter().map(|a| a.position).collect();

        let t0 = Instant::now();
        let ff = MMFFForceField::new(&mol, MMFFVariant::MMFF94s);
        let build_us = t0.elapsed().as_micros() as f64;

        // Warmup
        for _ in 0..20 {
            let (e, _g) = ff.calculate_energy_and_gradient(&coords);
            black_box(e);
        }

        // Steady state: run until min-ms of measured time
        let t_start = Instant::now();
        let mut ops = 0u64;
        let mut esum = 0.0;
        while t_start.elapsed() < min_dur {
            let (e, _g) = ff.calculate_energy_and_gradient(&coords);
            esum += e;
            ops += 1;
        }
        let eg_us = t_start.elapsed().as_micros() as f64;
        let us_per_op = eg_us / ops as f64;
        let us_per_atom = us_per_op / mol.atoms.len() as f64;

        // Suppress dead-code elimination; energy sum doubles as a sanity check.
        let _ = black_box(esum);
        out.push_str(&format!(
            "\"{name}\":{{\"n_atoms\":{},\"build_us\":{:.1},\"eg_us\":{:.1},\
             \"eg_ops\":{},\"us_per_op\":{:.3},\"us_per_atom\":{:.4}}}",
            mol.atoms.len(),
            build_us,
            eg_us,
            ops,
            us_per_op,
            us_per_atom
        ));
    }
    out.push('}');
    let _ = write!(std::io::stdout(), "{out}");
}
