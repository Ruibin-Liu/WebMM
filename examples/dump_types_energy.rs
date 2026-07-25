//! Dump WebMM atom types + energy (at SDF coords) for all validation SDFs -> JSON.
//! Run: cargo run --release --example dump_types_energy > scripts/val_set/webmm_ref.json
use std::io::Write;
use webmm::mmff::MMFFForceField;
use webmm::molecule::parser::parse_sdf;
use webmm::MMFFVariant;

fn main() {
    let dir = "scripts/val_set";
    let mut files: Vec<String> = std::fs::read_dir(dir)
        .expect("val_set dir")
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
        let coords: Vec<[f64; 3]> = mol.atoms.iter().map(|a| a.position).collect();
        let ff = MMFFForceField::new(&mol, MMFFVariant::MMFF94s);
        let energy = ff.calculate_energy(&coords);
        let types: Vec<u32> = ff.type_ids.iter().map(|&t| t as u32).collect();
        let charges: Vec<f64> = ff.charges.iter().map(|c| (c * 1e4).round() / 1e4).collect();
        out.push_str(&format!(
            "\"{}\":{{\"types\":{:?},\"charges\":{:?},\"energy\":{:.4},\"n_atoms\":{}}}",
            name, types, charges, energy, mol.atoms.len()
        ));
    }
    out.push('}');
    let _ = write!(std::io::stdout(), "{}", out);
}
