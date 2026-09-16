/// Conformer-ensemble end-to-end parity vs RDKit (ETKDGv3 + MMFF94s
/// optimize). Same seeds, N=30 per molecule; compare ensemble statistics.
use webmm::molecule::parser::parse_sdf;
use webmm::MMFFVariant;

fn main() {
    let refj: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string("/tmp/conf_parity/rdkit_ref.json").unwrap())
            .unwrap();
    let names = [
        "ethanol",
        "n_butane",
        "naphthalene",
        "aspirin",
        "threonine",
        "ibuprofen",
    ];
    for name in names {
        let sdf = std::fs::read_to_string(format!("/tmp/conf_parity/{name}.sdf")).unwrap();
        let mol = parse_sdf(&sdf).unwrap();
        let mut energies: Vec<f64> = Vec::new();
        let mut best: Option<(f64, Vec<[f64; 3]>)> = None;
        for seed in 42u64..42 + 30 {
            let config = webmm::etkdg::ETKDGConfig {
                random_seed: seed as i64,
                ..Default::default()
            };
            let mut coords = webmm::etkdg::generate_initial_coords_with_config(&mol, &config);
            if coords.is_empty() {
                continue;
            }
            // MMFF94s optimize
            let ff = webmm::mmff::MMFFForceField::new(&mol, MMFFVariant::MMFF94s);
            let conv = webmm::ConvergenceOptions {
                max_iterations: 2000,
                ..Default::default()
            };
            let r = webmm::optimizer::optimize(&ff, &coords, &conv);
            if !r.converged {
                continue;
            }
            let flat: Vec<[f64; 3]> = r.optimized_coords.clone();
            let e = ff.calculate_energy(&flat);
            energies.push(e);
            if best.as_ref().map(|(be, _)| e < *be).unwrap_or(true) {
                best = Some((e, flat.clone()));
            }
            coords = flat;
        }
        energies.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let n = energies.len();
        if n == 0 {
            println!("{name:10} NO CONFORMERS");
            continue;
        }
        let stats = &refj[name];
        println!("{name:10} n={:2}  min {:+8.3} (RDKit {:+8.3})  med {:+8.3} ({:+8.3})  p90 {:+8.3} ({:+8.3})",
            n, energies[0], stats["min"].as_f64().unwrap(),
            energies[n/2], stats["median"].as_f64().unwrap(),
            energies[(9*(n-1))/10], stats["p90"].as_f64().unwrap());
    }
}
