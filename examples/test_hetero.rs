use webmm::molecule::parser::parse_sdf;
use webmm::etkdg::generate_initial_coords;

fn max_planar_deviation(coords: &[[f64; 3]], atoms: &[usize]) -> f64 {
    let n = atoms.len() as f64;
    let cx = atoms.iter().map(|&i| coords[i][0]).sum::<f64>() / n;
    let cy = atoms.iter().map(|&i| coords[i][1]).sum::<f64>() / n;
    let cz = atoms.iter().map(|&i| coords[i][2]).sum::<f64>() / n;
    let mut cov = [[0.0f64; 3]; 3];
    for &idx in atoms {
        let dx = coords[idx][0] - cx; let dy = coords[idx][1] - cy; let dz = coords[idx][2] - cz;
        cov[0][0] += dx * dx; cov[0][1] += dx * dy; cov[0][2] += dx * dz;
        cov[1][1] += dy * dy; cov[1][2] += dy * dz; cov[2][2] += dz * dz;
    }
    cov[1][0] = cov[0][1]; cov[2][0] = cov[0][2]; cov[2][1] = cov[1][2];
    let normal = webmm::etkdg::eigenvector_smallest_eigenvalue_3x3(&cov);
    atoms.iter().map(|&idx| {
        let dx = coords[idx][0] - cx; let dy = coords[idx][1] - cy; let dz = coords[idx][2] - cz;
        (dx * normal[0] + dy * normal[1] + dz * normal[2]).abs()
    }).fold(0.0, f64::max)
}

fn test_molecule(name: &str, sdf: &str, all_atoms: &[usize], ring_atoms: &[usize]) {
    let mol = parse_sdf(sdf).unwrap();
    let mut max_ring = 0.0; let mut max_all = 0.0;
    for _ in 0..10 {
        let coords = generate_initial_coords(&mol);
        let rd = max_planar_deviation(&coords, ring_atoms);
        let ad = max_planar_deviation(&coords, all_atoms);
        if rd > max_ring { max_ring = rd; }
        if ad > max_all { max_all = ad; }
    }
    println!("{:15} ring_max={:.6}  all_max={:.6}", name, max_ring, max_all);
}

fn main() {
    let thiophene = r#"Thiophene
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

    let furan = r#"Furan
     RDKit          3D
  9  9  0  0  0  0  0  0  0  0999 V2000
    1.0874    0.4433    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.2774    1.5678    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0874   -0.4433    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2774   -1.5678    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.1630    0.4506    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.5847    2.5974    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1630   -0.4506    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5847   -2.5974    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
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

    let pyrrole = r#"Pyrrole
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

    println!("Planarity stress test (10 runs each):");
    test_molecule("Thiophene", thiophene, &[0,1,2,3,4,5,6,7,8], &[0,1,2,3,4]);
    test_molecule("Furan", furan, &[0,1,2,3,4,5,6,7,8], &[0,1,2,3,4]);
    test_molecule("Pyrrole", pyrrole, &[0,1,2,3,4,5,6,7,8,9], &[0,1,2,3,4]);
}
