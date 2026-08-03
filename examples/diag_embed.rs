//! General ETKDG embed diagnostic: parse SDF, embed with WebMM, print MMFF energy
//! breakdown + all bond lengths + all bond angles (to localize which geometry is
//! blown up for an outlier molecule). Usage: diag_embed <sdf>
use webmm::etkdg::{generate_initial_coords_with_config, ETKDGConfig};
use webmm::mmff::MMFFForceField;
use webmm::molecule::parser::parse_sdf;
use webmm::MMFFVariant;

fn dist(c: &[[f64; 3]], i: usize, j: usize) -> f64 {
    let d = [c[i][0] - c[j][0], c[i][1] - c[j][1], c[i][2] - c[j][2]];
    (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt()
}
fn angle(c: &[[f64; 3]], i: usize, j: usize, k: usize) -> f64 {
    let v1 = [c[i][0] - c[j][0], c[i][1] - c[j][1], c[i][2] - c[j][2]];
    let v2 = [c[k][0] - c[j][0], c[k][1] - c[j][1], c[k][2] - c[j][2]];
    let m = (v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2]).sqrt();
    let n = (v2[0] * v2[0] + v2[1] * v2[1] + v2[2] * v2[2]).sqrt();
    if m < 1e-12 || n < 1e-12 {
        return f64::NAN;
    }
    let dot = (v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2]) / (m * n);
    dot.clamp(-1.0, 1.0).acos().to_degrees()
}

fn main() {
    let path = std::env::args()
        .nth(1)
        .unwrap_or_else(|| "scripts/val_set/cyclopropane.sdf".to_string());
    let sdf = std::fs::read_to_string(&path).expect("read");
    let mol = parse_sdf(&sdf).expect("parse");
    let cfg = ETKDGConfig {
        random_seed: 42,
        use_macrocycle_torsions: true,
        use_macrocycle_14config: true,
        force_trans_amides: true,
        et_version: 2,
        max_attempts: 10,
        max_iterations: 2000,
        ..Default::default()
    };
    let coords = generate_initial_coords_with_config(&mol, &cfg);
    let ff = MMFFForceField::new(&mol, MMFFVariant::MMFF94s);
    let e = ff.calculate_energy(&coords);
    let bd = ff.calculate_energy_breakdown(&coords);
    println!(
        "total={e:.2} bond={:.1} angle={:.1} sb={:.1} tor={:.1} oop={:.1} vdw={:.1} elec={:.1}",
        bd.bond, bd.angle, bd.stretch_bend, bd.torsion, bd.oop, bd.vdw, bd.electrostatic
    );
    println!("\nbond lengths (atom i-j sym sym : len):");
    for b in &mol.bonds {
        println!(
            "  {}{}-{}{} : {:.3}",
            mol.atoms[b.atom1].symbol,
            b.atom1,
            mol.atoms[b.atom2].symbol,
            b.atom2,
            dist(&coords, b.atom1, b.atom2)
        );
    }
    println!("\nangles (i-j-k sym : deg):");
    for j in 0..mol.atoms.len() {
        let nbrs: Vec<usize> = mol.adjacency[j].to_vec();
        for xi in 0..nbrs.len() {
            for xk in (xi + 1)..nbrs.len() {
                let (i, k) = (nbrs[xi], nbrs[xk]);
                println!(
                    "  {}{}-{}{}-{}{} : {:.1}",
                    mol.atoms[i].symbol,
                    i,
                    mol.atoms[j].symbol,
                    j,
                    mol.atoms[k].symbol,
                    k,
                    angle(&coords, i, j, k)
                );
            }
        }
    }
}
