//! Dump WebMM ETKDG-embedded conformer MMFF energy for each validation SDF -> JSON.
//! Mirrors dump_types_energy.rs, but EMBEDS each molecule with WebMM ETKDG (rather
//! than reading fixed coords) and reports the single-point MMFF energy of the
//! embedded conformer — a strain/quality measure for the embedding pipeline.
//!
//! Run: cargo run --release --example dump_etkdg_geom [sdf_dir] > scripts/etkdg_val/webmm_ref.json
use std::io::Write;
use webmm::etkdg::{generate_initial_coords_with_config, ETKDGConfig};
use webmm::mmff::MMFFForceField;
use webmm::molecule::parser::parse_sdf;
use webmm::MMFFVariant;

fn main() {
    let dir = std::env::args()
        .nth(1)
        .unwrap_or_else(|| "scripts/val_set".to_string());
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
        // Embed fresh with WebMM ETKDG (seed 42, matching the RDKit ref). WebMM
        // re-embeds from connectivity, ignoring whatever coords the SDF carries.
        let config = ETKDGConfig {
            random_seed: 42,
            use_macrocycle_torsions: true,
            use_macrocycle_14config: true,
            force_trans_amides: true,
            et_version: 2,
            max_attempts: 10,
            max_iterations: 2000,
            ..Default::default()
        };
        let coords = generate_initial_coords_with_config(&mol, &config);
        let embed_ok = coords.len() == n
            && coords
                .iter()
                .all(|c| c[0].is_finite() && c[1].is_finite() && c[2].is_finite());
        if embed_ok {
            let ff = MMFFForceField::new(&mol, MMFFVariant::MMFF94s);
            let energy = ff.calculate_energy(&coords);
            out.push_str(&format!(
                "\"{}\":{{\"energy\":{:.4},\"n_atoms\":{},\"embed_ok\":true}}",
                name, energy, n
            ));
        } else {
            out.push_str(&format!(
                "\"{}\":{{\"energy\":null,\"n_atoms\":{},\"embed_ok\":false}}",
                name, n
            ));
        }
    }
    out.push('}');
    let _ = write!(std::io::stdout(), "{}", out);
}
