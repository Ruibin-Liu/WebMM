//! Run WebMM ETKDG on the RDKit #9143 reporter macrocycle and report amide
//! bond planarity. Amide is planar when the O=C-N-X dihedral ~ 0°/180°;
//! the #9143 bug produces ~90° (twisted).
use webmm::etkdg::{generate_initial_coords_with_config, ETKDGConfig};
use webmm::molecule::graph::find_rings;
use webmm::molecule::parser::parse_sdf;
use webmm::molecule::BondType;

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
    let path = std::env::args().nth(1).unwrap_or_else(|| "/tmp/macrocycle.sdf".to_string());
    let sdf = std::fs::read_to_string(&path).expect("read sdf");
    let mol = parse_sdf(&sdf).expect("parse");
    eprintln!("parsed: {} atoms, {} bonds", mol.atoms.len(), mol.bonds.len());

    let cfg = ETKDGConfig {
        random_seed: 0xf00d, // matches the RDKit #9143 regression-test seed
        use_macrocycle_torsions: true,
        use_macrocycle_14config: true,
        force_trans_amides: true,
        max_attempts: 50,
        max_iterations: 2000,
        et_version: 2,
        ..Default::default()
    };
    let coords = generate_initial_coords_with_config(&mol, &cfg);
    eprintln!("coords: {}", coords.len());
    if coords.len() != mol.atoms.len() {
        eprintln!("WARN: coord/atom count mismatch — embedding may have failed");
        return;
    }
    let finite = coords
        .iter()
        .all(|c| c[0].is_finite() && c[1].is_finite() && c[2].is_finite());
    eprintln!("all coords finite: {finite}");

    // find amide C-N single bonds where C has a =O neighbor
    let mut amides = Vec::new();
    for b in &mol.bonds {
        if b.bond_type != BondType::Single {
            continue;
        }
        let (a, c) = (b.atom1, b.atom2);
        let (sa, sc) = (&mol.atoms[a].symbol, &mol.atoms[c].symbol);
        let (cc, nidx) = if sa == "C" && sc == "N" {
            (a, c)
        } else if sa == "N" && sc == "C" {
            (c, a)
        } else {
            continue;
        };
        let mut ox = None;
        for &nb in &mol.adjacency[cc] {
            if nb == nidx {
                continue;
            }
            for bb in &mol.bonds {
                if ((bb.atom1 == cc && bb.atom2 == nb) || (bb.atom2 == cc && bb.atom1 == nb))
                    && bb.bond_type == BondType::Double
                    && mol.atoms[nb].symbol == "O"
                {
                    ox = Some(nb);
                }
            }
        }
        if let Some(o) = ox {
            amides.push((cc, nidx, o));
        }
    }
    eprintln!("amide C-N bonds found (C,N,O): {amides:?}");

    // ---- DIAGNOSE: does SSSR find the macrocycle? does ring_size_of_bond >= 9 for amides? ----
    let rings = find_rings(&mol);
    let mut sizes: Vec<usize> = rings.iter().map(|r| r.len()).collect();
    sizes.sort();
    eprintln!("\n=== RING DETECTION ===");
    eprintln!("find_rings: {} rings, sizes = {:?}", rings.len(), sizes);
    eprintln!("any ring >= 9 (macrocycle)? {}", rings.iter().any(|r| r.len() >= 9));
    for &(cc, nn, _o) in &amides {
        let containing: Vec<usize> = rings
            .iter()
            .filter(|r| r.contains(&cc) && r.contains(&nn))
            .map(|r| r.len())
            .collect();
        let first = rings
            .iter()
            .find(|r| r.contains(&cc) && r.contains(&nn))
            .map(|r| r.len());
        eprintln!(
            "  amide C{cc}-N{nn}: rings containing both = {:?}; ring_size_of_bond (first) = {:?} -> macrocycle branch fires? {}",
            containing, first, first.is_some_and(|s| s >= 9)
        );
    }
    eprintln!();

    // try a few seeds in case one embedding is unlucky
    println!("\n=== seed 0xf00d: amide planarity (O=C-N-X; ~0/180=planar, ~90=twisted) ===");
    report(&mol, &coords, &amides);

    let mut worst_overall = 0.0f64;
    for seed in [1i64, 7, 42, 100, 12345] {
        let cfg2 = ETKDGConfig {
            random_seed: seed,
            use_macrocycle_torsions: true,
            use_macrocycle_14config: true,
            force_trans_amides: true,
            max_attempts: 50,
            max_iterations: 2000,
            et_version: 2,
            ..Default::default()
        };
        let coords2 = generate_initial_coords_with_config(&mol, &cfg2);
        let mut worst = 0.0f64;
        for &(cc, nn, ox) in &amides {
            let x = n_sub(&mol, nn, cc);
            if let Some(xx) = x {
                let d = dihedral(&coords2, ox, cc, nn, xx);
                if d.is_finite() {
                    worst = worst.max((90.0 - (d - 90.0).abs()).abs()); // distance from 0/180
                    worst = worst.max(if d > 45.0 && d < 135.0 { 1.0 } else { 0.0 });
                }
            }
        }
        let any_twisted = amides.iter().any(|&(cc, nn, ox)| {
            n_sub(&mol, nn, cc).is_some_and(|xx| {
                let d = dihedral(&coords2, ox, cc, nn, xx);
                d.is_finite() && d > 60.0 && d < 120.0
            })
        });
        println!("seed {seed:>6}: any amide twisted(60-120°)? {any_twisted}");
        worst_overall = worst_overall.max(if any_twisted { 1.0 } else { 0.0 });
    }
    println!(
        "\nSUMMARY: {}",
        if worst_overall > 0.0 {
            "at least one seed produced a twisted amide"
        } else {
            "all seeds produced planar amides"
        }
    );
}

fn n_sub(mol: &webmm::molecule::Molecule, nn: usize, cc: usize) -> Option<usize> {
    let mut heavy = None;
    let mut hyd = None;
    for &nb in &mol.adjacency[nn] {
        if nb == cc {
            continue;
        }
        if mol.atoms[nb].symbol == "H" {
            hyd = Some(nb);
        } else {
            heavy = Some(nb);
        }
    }
    heavy.or(hyd)
}

fn report(mol: &webmm::molecule::Molecule, coords: &[[f64; 3]], amides: &[(usize, usize, usize)]) {
    for &(cc, nn, ox) in amides {
        let x = n_sub(mol, nn, cc);
        if let Some(xx) = x {
            let d = dihedral(coords, ox, cc, nn, xx);
            let verdict = if !(20.0..=160.0).contains(&d) {
                "PLANAR"
            } else if d > 70.0 && d < 110.0 {
                "*** TWISTED ~90 ***"
            } else {
                "nonplanar"
            };
            println!(
                "  amide C{cc}-N{nn}: O=C-N-{} dihedral = {d:6.1}°  [{verdict}]",
                mol.atoms[xx].symbol
            );
        }
    }
}
