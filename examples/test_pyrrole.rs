use webmm::molecule::parser::parse_sdf;

fn main() {
    let sdf = r#"Pyrrole
     RDKit          3D
 10 10  0  0  0  0  0  0  0  0999 V2000
    1.0759    0.3167    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.1438    1.3335    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0759   -0.3167    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.1438   -1.3335    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    2.1483    0.2560    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.2948    2.4051    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1483   -0.2560    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2948   -2.4051    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    1.0100 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0
  2  3  1  0
  3  4  2  0
  4  5  1  0
  5  1  1  0
  1  6  1  0
  2  7  1  0
  3  8  1  0
  4  9  1  0
  5 10  1  0
M  END"#;
    let mol = parse_sdf(sdf).unwrap();
    let aromatic = webmm::molecule::graph::get_aromatic_atoms(&mol);
    println!("Aromatic atoms: {:?}", aromatic);
    for i in 0..mol.atoms.len() {
        let hyb = webmm::molecule::graph::determine_hybridization(i, &mol);
        println!("Atom {} ({}): {:?}", i, mol.atoms[i].symbol, hyb);
    }
}
