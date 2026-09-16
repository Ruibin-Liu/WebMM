fn main() {
    for name in ["nitrate", "dihydrogen_phosphate"] {
        let sdf = std::fs::read_to_string(format!("/tmp/charged/{name}.sdf")).unwrap();
        let mol = webmm::molecule::parser::parse_sdf(&sdf).unwrap();
        let coords: Vec<[f64; 3]> = mol.atoms.iter().map(|a| a.position).collect();
        let ff = webmm::mmff::MMFFForceField::new(&mol, webmm::MMFFVariant::MMFF94s);
        let e = ff.calculate_energy(&coords);
        println!("== {name}  E = {e:+.4}");
        for (i, t) in ff.type_ids.iter().enumerate() {
            println!(
                "  {:2} {:2} deg {} chg {:+} -> type {}",
                i + 1,
                mol.atoms[i].symbol,
                mol.adjacency[i].len(),
                mol.atoms[i].charge as i32,
                t
            );
        }
    }
}
