//! Dump WebMM ETKDG-embedded-conformer MMFF energy for each validation SDF -> JSON.
//!
//! Mirrors gen_etkdg_ref.py: EMBEDS each molecule with WebMM ETKDG at K seeds
//! (env WEBMM_SEEDS, default "42,43,100,7") and reports the MMFF94s single-point
//! energy of each embedded conformer plus the mean/stddev across seeds. The mean
//! is the strain/quality measure compared by validate_etkdg.py.
//!
//! Run: cargo run --release --example dump_etkdg_geom [sdf_dir] > scripts/etkdg_val/webmm_ref.json
//!   WEBMM_SEEDS=42,43,100,7 cargo run --release --example dump_etkdg_geom
use std::io::Write;
use webmm::etkdg::{generate_initial_coords_with_config, ETKDGConfig};
use webmm::mmff::MMFFForceField;
use webmm::molecule::parser::parse_sdf;
use webmm::MMFFVariant;

fn parse_seeds() -> Vec<i64> {
    match std::env::var("WEBMM_SEEDS") {
        Ok(s) => s
            .split(',')
            .filter_map(|t| t.trim().parse::<i64>().ok())
            .collect(),
        Err(_) => vec![42, 43, 100, 7],
    }
}

fn main() {
    let dir = std::env::args()
        .nth(1)
        .unwrap_or_else(|| "scripts/val_set".to_string());
    let seeds = parse_seeds();
    eprintln!("# seeds {:?}", seeds);
    // M0+: toggle the RDKit-faithful embedding path for A/B measurement.
    let rdkit_faithful = std::env::var("WEBMM_RDKIT_FAITHFUL")
        .map(|v| !v.is_empty() && v != "0")
        .unwrap_or(false);
    eprintln!("# rdkit_faithful={rdkit_faithful}");
    let mut files: Vec<String> = std::fs::read_dir(&dir)
        .expect("sdf dir")
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
        let path = format!("{}/{}.sdf", dir, name);
        let txt = match std::fs::read_to_string(&path) {
            Ok(t) => t,
            Err(_) => {
                out.push_str(&format!("\"{}\":{{\"error\":\"read\"}}", name));
                continue;
            }
        };
        let mol = match parse_sdf(&txt) {
            Ok(m) => m,
            Err(e) => {
                out.push_str(&format!(
                    "\"{}\":{{\"error\":\"parse:{}\"}}",
                    name,
                    e.replace('"', "'")
                ));
                continue;
            }
        };
        let n = mol.atoms.len();
        let ff = MMFFForceField::new(&mol, MMFFVariant::MMFF94s);
        // Embed fresh with WebMM ETKDG at each seed. WebMM re-embeds from
        // connectivity, ignoring whatever coords the SDF carries.
        let mut seed_energies: Vec<(i64, f64)> = Vec::new();
        for &seed in &seeds {
            let config = ETKDGConfig {
                random_seed: seed,
                use_macrocycle_torsions: true,
                use_macrocycle_14config: true,
                force_trans_amides: true,
                et_version: 2,
                max_attempts: 10,
                max_iterations: 2000,
                rdkit_faithful,
                ..Default::default()
            };
            let coords = generate_initial_coords_with_config(&mol, &config);
            let embed_ok = coords.len() == n
                && coords
                    .iter()
                    .all(|c| c[0].is_finite() && c[1].is_finite() && c[2].is_finite());
            if embed_ok {
                seed_energies.push((seed, ff.calculate_energy(&coords)));
            }
        }
        let n_emb = seed_energies.len();
        if n_emb == 0 {
            out.push_str(&format!(
                "\"{}\":{{\"mean\":null,\"stddev\":null,\"seeds\":{{}},\"n_embedded\":0,\"n_atoms\":{},\"embed_ok\":false}}",
                name, n
            ));
            continue;
        }
        let mean = seed_energies.iter().map(|&(_, e)| e).sum::<f64>() / n_emb as f64;
        let var = seed_energies
            .iter()
            .map(|&(_, e)| (e - mean).powi(2))
            .sum::<f64>()
            / n_emb as f64;
        let stddev = var.sqrt();
        let mut seeds_json = String::new();
        for (i, &(s, e)) in seed_energies.iter().enumerate() {
            if i > 0 {
                seeds_json.push(',');
            }
            seeds_json.push_str(&format!("\"{}\":{:.4}", s, e));
        }
        out.push_str(&format!(
            "\"{}\":{{\"mean\":{:.4},\"stddev\":{:.4},\"seeds\":{{{}}},\"n_embedded\":{},\"n_atoms\":{},\"embed_ok\":{}}}",
            name,
            mean,
            stddev,
            seeds_json,
            n_emb,
            n,
            n_emb == seeds.len()
        ));
    }
    out.push('}');
    let _ = write!(std::io::stdout(), "{}", out);
}
