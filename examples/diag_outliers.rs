//! Diagnostic: embed ONE molecule with WebMM ETKDG and dump geometry + MMFF energy.
//!
//! Usage: cargo run --release --example diag_outliers -- <mol_name> [seed]
//! Prints per-atom (idx symbol hybridization x y z), bond lengths, key 1-3 angles,
//! and the MMFF94s single-point energy, as TSV/JSON-ish lines for comparison
//! against RDKit embeddings (see scripts/diag_compare.py).
use std::env;

use webmm::etkdg::{generate_initial_coords_with_config, ETKDGConfig};
use webmm::mmff::MMFFForceField;
use webmm::molecule::graph::determine_hybridization;
use webmm::molecule::parser::parse_sdf;
use webmm::MMFFVariant;

fn main() {
    let name = env::args()
        .nth(1)
        .expect("usage: diag_outliers <mol_name> [seed]");
    // Optional 3rd arg: xyz file with "idx x y z" lines to score with webmm MMFF
    // (for cross-checking webmm-MMFF on RDKit geometry).
    if let Some(xyz) = env::args().nth(3) {
        let txt = std::fs::read_to_string(format!("scripts/val_set/{}.sdf", name)).unwrap();
        let mol = parse_sdf(&txt).unwrap();
        let mut coords = vec![[0.0f64; 3]; mol.atoms.len()];
        for line in std::fs::read_to_string(&xyz).unwrap().lines() {
            let p: Vec<&str> = line.split_whitespace().collect();
            if p.len() >= 5 && (p[0] == "C" || p[0] == "A") {
                let i: usize = p[1].parse().unwrap();
                // formats: "C <idx> <sym> <x> <y> <z>" or "C <idx> <x> <y> <z>"
                let (x, y, z) = if p[2].parse::<f64>().is_ok() {
                    (p[2], p[3], p[4])
                } else {
                    (p[3], p[4], p[5])
                };
                coords[i] = [x.parse().unwrap(), y.parse().unwrap(), z.parse().unwrap()];
            }
        }
        let ff = MMFFForceField::new(&mol, MMFFVariant::MMFF94s);
        let bd = ff.calculate_energy_breakdown(&coords);
        println!(
            "# webmm-MMFF on external coords: {:.4} | bond={:.3} angle={:.3} sb={:.3} tor={:.3} oop={:.3} vdw={:.3} ele={:.3}",
            bd.total(), bd.bond, bd.angle, bd.stretch_bend, bd.torsion, bd.oop, bd.vdw, bd.electrostatic
        );
        for (i, a) in mol.atoms.iter().enumerate() {
            println!("CHG\t{}\t{}\t{:.6}", i, a.symbol, ff.charges[i]);
        }
        return;
    }
    let seed: i64 = env::args()
        .nth(2)
        .unwrap_or_else(|| "42".into())
        .parse()
        .unwrap();
    let txt = std::fs::read_to_string(format!("scripts/val_set/{}.sdf", name))
        .unwrap_or_else(|_| panic!("no scripts/val_set/{}.sdf", name));
    let mol = parse_sdf(&txt).unwrap();
    let config = ETKDGConfig {
        random_seed: seed,
        use_macrocycle_torsions: true,
        use_macrocycle_14config: true,
        force_trans_amides: true,
        et_version: 2,
        max_attempts: 10,
        max_iterations: 2000,
        ..Default::default()
    };
    let coords = generate_initial_coords_with_config(&mol, &config);
    let ff = MMFFForceField::new(&mol, MMFFVariant::MMFF94s);
    let e = ff.calculate_energy(&coords);
    println!(
        "# name={name} seed={seed} n_atoms={} mmff={e:.4}",
        mol.atoms.len()
    );

    // Atoms with WebMM hybridization
    println!("# atom idx sym hyb");
    for (i, a) in mol.atoms.iter().enumerate() {
        let hyb = determine_hybridization(i, &mol);
        let c = coords[i];
        println!(
            "A\t{}\t{}\t{:?}\t{:.4}\t{:.4}\t{:.4}",
            i, a.symbol, hyb, c[0], c[1], c[2]
        );
    }

    // Bonds
    println!("# bond i j order length");
    for b in &mol.bonds {
        let d = ((coords[b.atom1][0] - coords[b.atom2][0]).powi(2)
            + (coords[b.atom1][1] - coords[b.atom2][1]).powi(2)
            + (coords[b.atom1][2] - coords[b.atom2][2]).powi(2))
        .sqrt();
        println!("B\t{}\t{}\t{:?}\t{:.4}", b.atom1, b.atom2, b.bond_type, d);
    }

    // All 1-3 angles (i-j-k)
    let bonds: Vec<(usize, usize)> = mol.bonds.iter().map(|b| (b.atom1, b.atom2)).collect();
    let mut adj: Vec<Vec<usize>> = vec![Vec::new(); mol.atoms.len()];
    for &(a, b) in &bonds {
        adj[a].push(b);
        adj[b].push(a);
    }
    for i in 0..mol.atoms.len() {
        for a in &adj[i] {
            for b in &adj[i] {
                if a < b {
                    let v1 = [
                        coords[*a][0] - coords[i][0],
                        coords[*a][1] - coords[i][1],
                        coords[*a][2] - coords[i][2],
                    ];
                    let v2 = [
                        coords[*b][0] - coords[i][0],
                        coords[*b][1] - coords[i][1],
                        coords[*b][2] - coords[i][2],
                    ];
                    let n1 = (v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2]).sqrt();
                    let n2 = (v2[0] * v2[0] + v2[1] * v2[1] + v2[2] * v2[2]).sqrt();
                    let cos =
                        (v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2]) / (n1 * n2).max(1e-12);
                    let ang = cos.clamp(-1.0, 1.0).acos().to_degrees();
                    println!("ANG\t{}\t{}\t{}\t{:.2}", *a, i, *b, ang);
                }
            }
        }
    }
}
