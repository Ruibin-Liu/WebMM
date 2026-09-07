//! Print WebMM MMFF94s atom types for all test molecules -> JSON.
//! Run: cargo run --release --example dump_atom_types > scripts/webmm_atom_types.json

use std::io::Write;
use webmm::mmff::{MMFFAtomType, MMFFForceField};
use webmm::molecule::parser::parse_sdf;
use webmm::MMFFVariant;

fn type_id(t: &MMFFAtomType) -> u32 {
    // mirror mmff_type_id but that's pub; use debug-name map via the variant
    // simplest: use the existing public field ff.type_ids
    let _ = t;
    0
}

fn main() {
    let dirs = ["test_molecules", "scripts/test_molecules"];
    let mut seen = std::collections::BTreeSet::new();
    let mut files: Vec<String> = Vec::new();
    for d in dirs {
        if let Ok(entries) = std::fs::read_dir(d) {
            for e in entries.flatten() {
                let p = e.path();
                if p.extension().and_then(|s| s.to_str()) == Some("sdf") {
                    if let Some(name) = p.file_stem().and_then(|s| s.to_str()) {
                        let key = name.to_string();
                        if seen.insert(key.clone()) {
                            files.push(p.to_string_lossy().into_owned());
                        }
                    }
                }
            }
        }
    }
    files.sort();

    let mut out = String::from("{");
    let mut first_file = true;
    for f in &files {
        let stem = std::path::Path::new(f)
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or("?");
        if !first_file {
            out.push(',');
        }
        first_file = false;
        out.push_str(&format!("\"{}\":[", stem));
        let txt = match std::fs::read_to_string(f) {
            Ok(t) => t,
            Err(_) => {
                out.push_str("{\"error\":\"read failed\"}]");
                continue;
            }
        };
        // split multi-mol SDF on "$$$$"
        let mols: Vec<&str> = txt.split("$$$$").filter(|s| !s.trim().is_empty()).collect();
        let mut first_mol = true;
        for m in &mols {
            if !first_mol {
                out.push(',');
            }
            first_mol = false;
            let mol = match parse_sdf(m) {
                Ok(m) => m,
                Err(e) => {
                    out.push_str(&format!("{{\"error\":\"parse: {}\"}}", e.replace('"', "'")));
                    continue;
                }
            };
            let ff = MMFFForceField::new(&mol, MMFFVariant::MMFF94s);
            let syms: Vec<&str> = mol.atoms.iter().map(|a| a.symbol.as_str()).collect();
            let charges: Vec<i32> = mol.atoms.iter().map(|a| a.charge as i32).collect();
            let type_ids: Vec<u32> = ff.type_ids.iter().map(|&t| t as u32).collect();
            let bonds: Vec<[u32; 3]> = mol
                .bonds
                .iter()
                .map(|b| {
                    let order = match b.bond_type {
                        webmm::molecule::BondType::Single => 1,
                        webmm::molecule::BondType::Double => 2,
                        webmm::molecule::BondType::Triple => 3,
                        webmm::molecule::BondType::Aromatic => 4,
                    };
                    [b.atom1 as u32, b.atom2 as u32, order]
                })
                .collect();
            let aromatic: Vec<u32> = mol
                .atoms
                .iter()
                .enumerate()
                .map(|(i, _)| {
                    if ff.rings.iter().any(|r| r.contains(&i)) {
                        // crude: use graph aromaticity
                    }
                    0
                })
                .collect();
            let _ = (
                &syms,
                &charges,
                &bonds,
                &aromatic,
                type_id(&MMFFAtomType::C_3),
            );
            out.push_str(&format!(
                "{{\"types\":{:?},\"syms\":{:?},\"charges\":{:?},\"bonds\":{:?}}}",
                type_ids, syms, charges, bonds
            ));
        }
        out.push(']');
    }
    out.push('}');
    let stdout = std::io::stdout();
    let mut h = stdout.lock();
    let _ = write!(h, "{}", out);
}
