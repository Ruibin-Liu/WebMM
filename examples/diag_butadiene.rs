//! Diagnose WebMM ETKDG embedding of a small conjugated molecule (butadiene default):
//! prints MMFF energy breakdown + key bond lengths + the conjugated dihedral.
use webmm::etkdg::{generate_initial_coords_with_config, ETKDGConfig};
use webmm::mmff::MMFFForceField;
use webmm::molecule::parser::parse_sdf;
use webmm::MMFFVariant;

fn dist(c: &[[f64; 3]], i: usize, j: usize) -> f64 {
    let d = [c[i][0] - c[j][0], c[i][1] - c[j][1], c[i][2] - c[j][2]];
    (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt()
}
fn dihedral(c: &[[f64; 3]], i: usize, j: usize, k: usize, l: usize) -> f64 {
    let b1 = [c[j][0] - c[i][0], c[j][1] - c[i][1], c[j][2] - c[i][2]];
    let b2 = [c[k][0] - c[j][0], c[k][1] - c[j][1], c[k][2] - c[j][2]];
    let b3 = [c[l][0] - c[k][0], c[l][1] - c[k][1], c[l][2] - c[k][2]];
    let n1 = [
        b1[1] * b2[2] - b1[2] * b2[1],
        b1[2] * b2[0] - b1[0] * b2[2],
        b1[0] * b2[1] - b1[1] * b2[0],
    ];
    let n2 = [
        b2[1] * b3[2] - b2[2] * b3[1],
        b2[2] * b3[0] - b2[0] * b3[2],
        b2[0] * b3[1] - b2[1] * b3[0],
    ];
    let m = (n1[0] * n1[0] + n1[1] * n1[1] + n1[2] * n1[2]).sqrt();
    let nn = (n2[0] * n2[0] + n2[1] * n2[1] + n2[2] * n2[2]).sqrt();
    if m < 1e-12 || nn < 1e-12 {
        return f64::NAN;
    }
    let dot = (n1[0] * n2[0] + n1[1] * n2[1] + n1[2] * n2[2]) / (m * nn);
    dot.clamp(-1.0, 1.0).acos().to_degrees()
}

fn main() {
    let path = std::env::args()
        .nth(1)
        .unwrap_or_else(|| "scripts/val_set/butadiene.sdf".to_string());
    let sdf = std::fs::read_to_string(&path).expect("read");
    let mol = parse_sdf(&sdf).expect("parse");
    for seed in [42i64, 1, 7, 100] {
        let cfg = ETKDGConfig {
            random_seed: seed,
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
        if seed == 42 {
            println!("(types: {:?})", ff.atom_types);
        }
        let e = ff.calculate_energy(&coords);
        let bd = ff.calculate_energy_breakdown(&coords);
        println!("\n=== seed {seed}: total energy = {e:.2} ===");
        println!(
            "  breakdown: bond={:.1} angle={:.1} sb={:.1} tor={:.1} oop={:.1} vdw={:.1} elec={:.1}",
            bd.bond, bd.angle, bd.stretch_bend, bd.torsion, bd.oop, bd.vdw, bd.electrostatic
        );
        println!(
            "  C0=C1={:.3} (exp 1.34)  C1-C2={:.3} (exp 1.47)  C2=C3={:.3} (exp 1.34)",
            dist(&coords, 0, 1),
            dist(&coords, 1, 2),
            dist(&coords, 2, 3)
        );
        println!(
            "  dihedral C0-C1-C2-C3 = {:.1}° (s-trans ~180 / s-cis ~0 / twisted ~90)",
            dihedral(&coords, 0, 1, 2, 3)
        );
    }
}
