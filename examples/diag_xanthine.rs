//! Diagnose xanthine ETKDG embedding: ring torsions + MMFF breakdown.
use webmm::etkdg::{generate_initial_coords_with_config, ETKDGConfig};
use webmm::mmff::MMFFForceField;
use webmm::molecule::parser::parse_sdf;
use webmm::molecule::graph::{find_rings, get_aromatic_atoms};
use webmm::MMFFVariant;

fn dihedral(c: &[[f64; 3]], i: usize, j: usize, k: usize, l: usize) -> f64 {
    let b1 = [c[j][0] - c[i][0], c[j][1] - c[i][1], c[j][2] - c[i][2]];
    let b2 = [c[k][0] - c[j][0], c[k][1] - c[j][1], c[k][2] - c[j][2]];
    let b3 = [c[l][0] - c[k][0], c[l][1] - c[k][1], c[l][2] - c[k][2]];
    let b2n = (b2[0] * b2[0] + b2[1] * b2[1] + b2[2] * b2[2]).sqrt();
    let b2u = [b2[0] / b2n, b2[1] / b2n, b2[2] / b2n];
    let v = [
        b3[0] - b2u[0] * (b3[0] * b2u[0] + b3[1] * b2u[1] + b3[2] * b2u[2]),
        b3[1] - b2u[1] * (b3[0] * b2u[0] + b3[1] * b2u[1] + b3[2] * b2u[2]),
        b3[2] - b2u[2] * (b3[0] * b2u[0] + b3[1] * b2u[1] + b3[2] * b2u[2]),
    ];
    let w = [
        b1[0] - b2u[0] * (b1[0] * b2u[0] + b1[1] * b2u[1] + b1[2] * b2u[2]),
        b1[1] - b2u[1] * (b1[0] * b2u[0] + b1[1] * b2u[1] + b1[2] * b2u[2]),
        b1[2] - b2u[2] * (b1[0] * b2u[0] + b1[1] * b2u[1] + b1[2] * b2u[2]),
    ];
    let x = v[0] * w[0] + v[1] * w[1] + v[2] * w[2];
    let cx = b2u[1] * v[2] - b2u[2] * v[1];
    let cy = b2u[2] * v[0] - b2u[0] * v[2];
    let cz = b2u[0] * v[1] - b2u[1] * v[0];
    let y = cx * w[0] + cy * w[1] + cz * w[2];
    y.atan2(x)
}

fn main() {
    let sdf = std::fs::read_to_string("scripts/val_set/xanthine.sdf").unwrap();
    let mol = parse_sdf(&sdf).unwrap();
    let rings = find_rings(&mol);
    let arom = get_aromatic_atoms(&mol);

    // Embed
    let cfg = ETKDGConfig { random_seed: 42, ..Default::default() };
    let coords = generate_initial_coords_with_config(&mol, &cfg);

    let ff = MMFFForceField::new(&mol, MMFFVariant::MMFF94s);
    let (e, _) = ff.calculate_energy_and_gradient(&coords);
    println!("embedded MMFF energy: {:.3} (RDKit ref ~-142.8)", e);

    // Ring torsions
    for (ri, ring) in rings.iter().enumerate() {
        let is_arom = ring.iter().all(|a| arom.contains(a));
        if ring.len() < 4 { continue; }
        let tors: Vec<f64> = (0..ring.len()).map(|s| {
            dihedral(&coords, ring[s], ring[(s+1)%ring.len()], ring[(s+2)%ring.len()], ring[(s+3)%ring.len()]) * 180.0 / std::f64::consts::PI
        }).collect();
        println!("ring {} (size {}, arom={}): torsions {:?}", ri, ring.len(), is_arom,
            tors.iter().map(|t| format!("{:.0}", t)).collect::<Vec<_>>());
    }

    // Breakdown: embedded vs SDF (optimized) coords
    let b_emb = ff.calculate_energy_breakdown(&coords);
    let sdf_coords: Vec<[f64; 3]> = mol.atoms.iter().map(|a| a.position).collect();
    let b_sdf = ff.calculate_energy_breakdown(&sdf_coords);
    println!("\nbreakdown        embedded    SDF-ref    diff");
    for (name, e, s) in [
        ("bond", b_emb.bond, b_sdf.bond),
        ("angle", b_emb.angle, b_sdf.angle),
        ("strbnd", b_emb.stretch_bend, b_sdf.stretch_bend),
        ("torsion", b_emb.torsion, b_sdf.torsion),
        ("oop", b_emb.oop, b_sdf.oop),
        ("vdw", b_emb.vdw, b_sdf.vdw),
        ("elec", b_emb.electrostatic, b_sdf.electrostatic),
    ] {
        println!("  {:<10} {:9.2} {:9.2} {:9.2}", name, e, s, e - s);
    }
}
