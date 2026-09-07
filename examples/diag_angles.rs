//! Angle strain diagnosis for any molecule: worst angles + 1-4 bound check.
use std::collections::HashMap;
use webmm::etkdg::{generate_initial_coords_with_config, ETKDGConfig};
use webmm::mmff::MMFFForceField;
use webmm::molecule::parser::parse_sdf;
use webmm::MMFFVariant;

fn bond_angle_deg(c: &[[f64; 3]], i: usize, j: usize, k: usize) -> f64 {
    let v1 = [c[i][0] - c[j][0], c[i][1] - c[j][1], c[i][2] - c[j][2]];
    let v2 = [c[k][0] - c[j][0], c[k][1] - c[j][1], c[k][2] - c[j][2]];
    let n1 = (v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2]).sqrt();
    let n2 = (v2[0] * v2[0] + v2[1] * v2[1] + v2[2] * v2[2]).sqrt();
    ((v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2]) / (n1 * n2))
        .clamp(-1.0, 1.0)
        .acos()
        * 180.0
        / std::f64::consts::PI
}

fn main() {
    let name = std::env::args()
        .nth(1)
        .unwrap_or_else(|| "caffeine".to_string());
    let sdf = std::fs::read_to_string(format!("scripts/val_set/{}.sdf", name)).unwrap();
    let mol = parse_sdf(&sdf).unwrap();
    let ff = MMFFForceField::new(&mol, MMFFVariant::MMFF94s);
    let cfg = ETKDGConfig {
        random_seed: 42,
        ..Default::default()
    };
    let coords = generate_initial_coords_with_config(&mol, &cfg);
    let sdf_c: Vec<[f64; 3]> = mol.atoms.iter().map(|a| a.position).collect();

    let (e_emb, _) = ff.calculate_energy_and_gradient(&coords);
    let (e_sdf, _) = ff.calculate_energy_and_gradient(&sdf_c);
    println!(
        "{}: embedded {:.1}  SDF {:.1}  diff {:.1}",
        name,
        e_emb,
        e_sdf,
        e_emb - e_sdf
    );

    // Worst angles: embedded vs equilibrium
    let theta0 = ff.per_angle_theta0();
    let mut t0map: HashMap<(usize, usize, usize), f64> = HashMap::new();
    for &(i, j, k, t0) in &theta0 {
        t0map.insert((i, j, k), t0); // theta0 already in degrees
    }
    let mut strained: Vec<(f64, usize, usize, usize, f64, f64, f64)> = Vec::new();
    for &(i, j, k, _t0) in &theta0 {
        let t0d = t0map[&(i, j, k)];
        let emb = bond_angle_deg(&coords, i, j, k);
        let sdf = bond_angle_deg(&sdf_c, i, j, k);
        strained.push((
            (emb - t0d).abs() + (emb - sdf).abs(),
            i,
            j,
            k,
            emb,
            sdf,
            t0d,
        ));
    }
    strained.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap());
    let sym = |idx: usize| mol.atoms[idx].symbol.clone();
    println!("\n  worst angles:                emb    SDF   theta0  |emb-eq|");
    for &(_, i, j, k, emb, sdf, t0d) in strained.iter().take(8) {
        println!(
            "    {}({})-{}({})-{}({}):  {:5.1} {:5.1} {:6.1}   {:.1}",
            sym(i),
            i,
            sym(j),
            j,
            sym(k),
            k,
            emb,
            sdf,
            t0d,
            (emb - t0d).abs()
        );
    }
}
