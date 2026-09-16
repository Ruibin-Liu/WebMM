fn main() {
    let refs: std::collections::HashMap<String, f64> =
        serde_json::from_str(&std::fs::read_to_string("/tmp/charged/refs.json").unwrap()).unwrap();
    for (name, r) in &refs {
        let sdf = std::fs::read_to_string(format!("/tmp/charged/{name}.sdf")).unwrap();
        let mol = webmm::molecule::parser::parse_sdf(&sdf).unwrap();
        let coords: Vec<[f64; 3]> = mol.atoms.iter().map(|a| a.position).collect();
        if let Err(e) = webmm::mmff::MMFFForceField::check_mmff_support(&mol) {
            println!("{name:24} REFUSED: {e}");
            continue;
        }
        let ff = webmm::mmff::MMFFForceField::new(&mol, webmm::MMFFVariant::MMFF94s);
        let e = ff.calculate_energy(&coords);
        println!("{name:24} ours {e:+.4}  RDKit {r:+.4}  d {:+.4}", e - r);
    }
}
