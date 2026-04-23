use webmm::molecule::parser::parse_sdf;
use webmm::etkdg::generate_initial_coords;

fn max_planar_deviation(coords: &[[f64; 3]], atoms: &[usize]) -> f64 {
    let n = atoms.len() as f64;
    let cx = atoms.iter().map(|&i| coords[i][0]).sum::<f64>() / n;
    let cy = atoms.iter().map(|&i| coords[i][1]).sum::<f64>() / n;
    let cz = atoms.iter().map(|&i| coords[i][2]).sum::<f64>() / n;
    let mut cov = [[0.0f64; 3]; 3];
    for &idx in atoms {
        let dx = coords[idx][0] - cx;
        let dy = coords[idx][1] - cy;
        let dz = coords[idx][2] - cz;
        cov[0][0] += dx * dx; cov[0][1] += dx * dy; cov[0][2] += dx * dz;
        cov[1][1] += dy * dy; cov[1][2] += dy * dz; cov[2][2] += dz * dz;
    }
    cov[1][0] = cov[0][1]; cov[2][0] = cov[0][2]; cov[2][1] = cov[1][2];
    let normal = webmm::etkdg::eigenvector_smallest_eigenvalue_3x3(&cov);
    atoms.iter().map(|&idx| {
        let dx = coords[idx][0] - cx;
        let dy = coords[idx][1] - cy;
        let dz = coords[idx][2] - cz;
        (dx * normal[0] + dy * normal[1] + dz * normal[2]).abs()
    }).fold(0.0, f64::max)
}

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
    let all_atoms: Vec<usize> = (0..9).collect();
    let ring_atoms: Vec<usize> = vec![0, 1, 2, 3, 4];
    
    for run in 0..5 {
        let coords = generate_initial_coords(&mol);
        let all_dev = max_planar_deviation(&coords, &all_atoms);
        let ring_dev = max_planar_deviation(&coords, &ring_atoms);
        println!("Run {}: ring_dev={:.6} Å  all_dev={:.6} Å", run, ring_dev, all_dev);
    }
}
