use webmm::molecule::parser::parse_sdf;

fn main() {
    let sdf = r#"Thiophene
     RDKit          3D

  9  9  0  0  0  0  0  0  0  0999 V2000
    1.1149    0.2097    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.2650    1.2632    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.1149   -0.2097    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2650   -1.2632    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 S   0  0  0  0  0  0  0  0  0  0  0  0
    2.1778    0.3308    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.5666    2.2962    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1778   -0.3308    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5666   -2.2962    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0
  2  3  1  0
  3  4  2  0
  4  5  1  0
  5  1  1  0
  1  6  1  0
  2  7  1  0
  3  8  1  0
  4  9  1  0
M  END"#;
    let mol = parse_sdf(sdf).unwrap();
    let aromatic = webmm::molecule::graph::get_aromatic_atoms(&mol);
    println!("Aromatic atoms: {:?}", aromatic);
    println!("Atoms: {:?}", mol.atoms.iter().map(|a| &a.symbol).collect::<Vec<_>>());
}
