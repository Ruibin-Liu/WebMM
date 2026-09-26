//! MD / metadynamics performance audit.
//! - MD throughput: steps/ms x molecule x engine (NVE, fixed seed)
//! - Metad degradation: per-step cost vs deposited hill count
use std::rc::Rc;
use webmm::forces::ForceField;
use webmm::md::{MDConfig, MDRunner};
use webmm::metad::{DihedralCV, MetaDConfig, MetaDynamics};
use webmm::mmff::{MMFFForceField, MMFFVariant};
use webmm::molecule::parser::parse_sdf;

fn mol_of(name: &str) -> webmm::molecule::Molecule {
    let dir = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/fixtures/conformers/");
    parse_sdf(&std::fs::read_to_string(format!("{dir}{name}.sdf")).unwrap()).unwrap()
}

fn md_bench() {
    for name in ["ethanol", "aspirin", "ibuprofen"] {
        let mol = mol_of(name);
        let coords: Vec<[f64; 3]> = mol.atoms.iter().map(|a| a.position).collect();
        for engine in ["MMFF94s", "GFN-FF"] {
            let ff: Rc<dyn ForceField> = if engine == "MMFF94s" {
                Rc::new(MMFFForceField::new(&mol, MMFFVariant::MMFF94s))
            } else {
                let at: Vec<usize> = mol.atoms.iter().map(|a| a.atomic_number as usize).collect();
                let q = mol.atoms.iter().map(|a| a.charge).sum::<f64>().round();
                Rc::new(webmm::gfnff::GfnffForceField::new(&at, &coords, q))
            };
            let cfg = MDConfig {
                dt_fs: 0.5,
                temperature_k: 300.0,
                friction_per_ps: 0.0,
                ..Default::default()
            };
            let mut runner = MDRunner::from_molecule(ff, &mol, cfg);
            runner.rescale_temperature(300.0);
            for _ in 0..200 {
                runner.step();
            }
            let n = 3000;
            let t = std::time::Instant::now();
            for _ in 0..n {
                runner.step();
            }
            let us_per = t.elapsed().as_secs_f64() * 1e6 / n as f64;
            println!(
                "{name:>9} {engine:>8}: {us_per:7.1} us/step ({:6.0} steps/ms)",
                1000.0 / us_per
            );
        }
    }
}

fn metad_bench() {
    let mol = mol_of("ibuprofen");

    // dihedral CV over a side-chain torsion (any 4 heavy atoms)
    let n = mol.atoms.len();
    let cv = DihedralCV::new(0, 1, 2, 9.min(n - 1));
    let mc = MetaDConfig {
        hill_height: 0.3,
        hill_width: 0.2,
        deposit_interval: 10, // build hills fast for the audit
        bias_factor: 0.0,
        temperature_k: 300.0,
    };
    let metad = MetaDynamics::new(
        Box::new(MMFFForceField::new(&mol, MMFFVariant::MMFF94s)),
        Box::new(cv),
        mc,
    );
    let cfg = MDConfig {
        dt_fs: 0.5,
        temperature_k: 300.0,
        friction_per_ps: 5.0,
        ..Default::default()
    };
    let mut runner = MDRunner::from_molecule(Rc::new(metad), &mol, cfg);
    runner.rescale_temperature(300.0);
    // measure per-step time as hills accumulate: windows of 2000 steps
    for window in 0..10 {
        let t = std::time::Instant::now();
        for _ in 0..2000 {
            runner.step();
        }
        let us = t.elapsed().as_secs_f64() * 1e6 / 2000.0;
        let hills = (window + 1) * 2000 / mc.deposit_interval;
        println!("metad window {window:>2} (~{hills:>5} hills): {us:7.1} us/step");
    }
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if !args.contains(&"--metad".to_string()) {
        md_bench();
    }
    if !args.contains(&"--md".to_string()) {
        metad_bench();
    }
}
