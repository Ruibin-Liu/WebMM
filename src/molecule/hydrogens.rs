//! Explicit-hydrogen addition for heavy-atom SDF inputs (Workbench v1.1).
//!
//! H counts come from standard organic valences (aromatic bonds count 1.5);
//! the graph is extended with single bonds to the new H atoms. Coordinates
//! are placed away from the neighbor centroid — good enough as an L-BFGS
//! start, and irrelevant for the ETKDG path which regenerates all coordinates
//! anyway. Molecules that already contain explicit H are returned unchanged
//! (no double addition).
//!
//! Formal charges round-trip via `M  CHG` records (the parser's authoritative
//! source), with charge-dependent standard valences (O⁻ → 1, N⁺ → 4, …) so
//! carboxylates/nitro groups get no spurious H.

use super::{Atom, Bond, BondStereo, BondType, Molecule};

/// Standard valence for the organic subset, charge-adjusted.
fn std_valence(z: usize, charge: f64) -> f64 {
    let q = charge as i32;
    match z {
        1 => 1.0,
        5 => 3.0, // B
        6 => {
            if q != 0 {
                3.0
            } else {
                4.0
            }
        } // C±
        7 => match q {
            1 => 4.0,
            -1 => 2.0,
            -2 => 1.0,
            _ => 3.0,
        },
        8 => match q {
            1 => 3.0,
            -1 => 1.0,
            _ => 2.0,
        },
        9 | 17 | 35 | 53 => 1.0, // F Cl Br I
        14 => 4.0,               // Si
        15 => 3.0,               // P
        16 => match q {
            1 => 3.0,
            -1 => 1.0,
            _ => 2.0,
        }, // S
        _ => 0.0,                // unknown → no H
    }
}

fn bond_order(bt: &BondType) -> f64 {
    match bt {
        BondType::Single => 1.0,
        BondType::Double => 2.0,
        BondType::Triple => 3.0,
        BondType::Aromatic => 1.5,
    }
}

fn bond_code(bt: &BondType) -> u8 {
    match bt {
        BondType::Single => 1,
        BondType::Double => 2,
        BondType::Triple => 3,
        BondType::Aromatic => 4,
    }
}

/// Add explicit hydrogens to `mol`. No-op when any explicit H is present.
pub fn add_hydrogens(mol: &Molecule) -> Molecule {
    if mol.atoms.iter().any(|a| a.atomic_number == 1) {
        return mol.clone();
    }
    let n = mol.atoms.len();

    // explicit valence per atom (aromatic counts 1.5)
    let mut explicit = vec![0.0f64; n];
    for b in &mol.bonds {
        let o = bond_order(&b.bond_type);
        explicit[b.atom1] += o;
        explicit[b.atom2] += o;
    }

    // H count per heavy atom
    let n_h_per: Vec<usize> = mol
        .atoms
        .iter()
        .enumerate()
        .map(|(i, a)| {
            let need = (std_valence(a.atomic_number as usize, a.charge) - explicit[i]).round();
            if need > 0.0 {
                need as usize
            } else {
                0
            }
        })
        .collect();
    let total_h: usize = n_h_per.iter().sum();
    if total_h == 0 {
        return mol.clone();
    }

    // build the extended molecule
    let mut atoms = Vec::with_capacity(n + total_h);
    atoms.extend(mol.atoms.iter().cloned());
    let mut bonds = mol.bonds.clone();

    // H placement: for 3D inputs, standard geometric attachment (unused-valence
    // directions; tetrahedral spread for one heavy neighbor); the L-BFGS
    // optimization relaxes them, and the ETKDG path regenerates all
    // coordinates anyway (heavy-only embedding keeps it fast).
    for (i, &n_hi) in n_h_per.iter().enumerate() {
        if n_hi == 0 {
            continue;
        }
        let p = mol.atoms[i].position;
        // unit vectors from this atom to its heavy neighbors + bond length
        let dirs: Vec<[f64; 3]> = mol.adjacency[i]
            .iter()
            .filter(|&&j| mol.atoms[j].atomic_number != 1)
            .map(|&j| {
                let v = [
                    mol.atoms[j].position[0] - p[0],
                    mol.atoms[j].position[1] - p[1],
                    mol.atoms[j].position[2] - p[2],
                ];
                let len = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt().max(1e-6);
                [v[0] / len, v[1] / len, v[2] / len]
            })
            .collect();
        let l = match mol.atoms[i].atomic_number as usize {
            8 => 0.97,  // O-H
            7 => 1.01,  // N-H
            16 => 1.34, // S-H
            _ => 1.09,  // C-H default
        };
        // k=1: three/few H in a tetrahedral fan opposite the single bond;
        // k=2: one on the reversed bisector, one perpendicular;
        // k>=3: one on the reversed vector sum (the only unused direction).
        let h_dirs: Vec<[f64; 3]> = match dirs.len() {
            1 => {
                let d = dirs[0];
                let refv = if d[2].abs() < 0.9 {
                    [0.0, 0.0, 1.0]
                } else {
                    [1.0, 0.0, 0.0]
                };
                let u = normalize3(cross3(d, refv));
                let w = cross3(d, u);
                // tetrahedral: 109.5° from the bond axis, cos = -1/3
                let n_fan = n_hi.min(3);
                (0..n_fan)
                    .map(|k| {
                        let phi = k as f64 * 2.0 * std::f64::consts::PI / n_fan as f64;
                        let cu = phi.cos() * (2.0f64.sqrt() * 2.0 / 3.0);
                        let sw = phi.sin() * (2.0f64.sqrt() * 2.0 / 3.0);
                        normalize3([
                            -d[0] / 3.0 + cu * u[0] + sw * w[0],
                            -d[1] / 3.0 + cu * u[1] + sw * w[1],
                            -d[2] / 3.0 + cu * u[2] + sw * w[2],
                        ])
                    })
                    .chain(std::iter::repeat(normalize3([-d[0], -d[1], -d[2]])))
                    .take(n_hi)
                    .collect()
            }
            2 => {
                let b = normalize3([
                    dirs[0][0] + dirs[1][0],
                    dirs[0][1] + dirs[1][1],
                    dirs[0][2] + dirs[1][2],
                ]);
                let perp = normalize3(cross3(dirs[0], dirs[1]));
                let mut v = vec![b];
                v.push(perp);
                v
            }
            _ => {
                let mut sum = [0.0f64; 3];
                for d in &dirs {
                    sum[0] -= d[0];
                    sum[1] -= d[1];
                    sum[2] -= d[2];
                }
                vec![normalize3(sum)]
            }
        };
        for k in 0..n_hi {
            let hd = h_dirs[k.min(h_dirs.len() - 1)];
            atoms.push(Atom {
                symbol: "H".to_string(),
                atomic_number: 1,
                mass: 1.008,
                charge: 0.0,
                position: [p[0] + l * hd[0], p[1] + l * hd[1], p[2] + l * hd[2]],
                index: atoms.len(),
                stereo_parity: 0,
            });
            bonds.push(Bond {
                atom1: i,
                atom2: atoms.len() - 1,
                bond_type: BondType::Single,
                stereo: BondStereo::None,
            });
        }
    }

    // adjacency rebuild (cached lists of NEIGHBOR ATOM indices, same
    // convention as build_adjacency_list_from_bonds)
    let mut adjacency = vec![Vec::new(); atoms.len()];
    for b in &bonds {
        adjacency[b.atom1].push(b.atom2);
        adjacency[b.atom2].push(b.atom1);
    }

    Molecule {
        atoms,
        bonds,
        name: mol.name.clone(),
        adjacency,
    }
}

fn normalize3(v: [f64; 3]) -> [f64; 3] {
    let len = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt().max(1e-12);
    [v[0] / len, v[1] / len, v[2] / len]
}

fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// Serialize a Molecule back to a V2000 molblock (title + program + comment +
/// counts + atoms + bonds + M  CHG records + M  END). Formal charges from
/// `atom.charge` are emitted as M  CHG (the parser's authoritative form).
pub fn to_molblock(mol: &Molecule) -> String {
    let mut out = String::new();
    out.push_str(&mol.name);
    out.push('\n');
    out.push_str("     WebMM          3D\n");
    out.push('\n');
    let na = mol.atoms.len();
    let nb = mol.bonds.len();
    out.push_str(&format!(
        "{:>3}{:>3}  0  0  0  0  0  0  0  0999 V2000\n",
        na, nb
    ));
    for a in &mol.atoms {
        out.push_str(&format!(
            "{:>10.4}{:>10.4}{:>10.4} {:<3} 0  0  0  0  0  0  0  0  0  0  0  0\n",
            a.position[0], a.position[1], a.position[2], a.symbol
        ));
    }
    for b in &mol.bonds {
        out.push_str(&format!(
            "{:>3}{:>3}{:>3}  0\n",
            b.atom1 + 1,
            b.atom2 + 1,
            bond_code(&b.bond_type)
        ));
    }
    // formal charges
    let charged: Vec<(usize, i64)> = mol
        .atoms
        .iter()
        .enumerate()
        .filter(|(_, a)| a.charge.abs() > 0.25)
        .map(|(i, a)| (i + 1, a.charge.round() as i64))
        .collect();
    if !charged.is_empty() {
        out.push_str(&format!("M  CHG{:>3}", charged.len()));
        for (idx, q) in &charged {
            out.push_str(&format!("{:>4}{:>4}", idx, q));
        }
        out.push('\n');
    }
    out.push_str("M  END\n");
    out
}

/// Convenience: parse an SDF, add hydrogens, return the molblock text.
pub fn molblock_with_h(sdf: &str) -> Result<String, String> {
    let mol = super::parser::parse_sdf(sdf)?;
    let with_h = add_hydrogens(&mol);
    Ok(to_molblock(&with_h))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn parse_mb(mb: &str) -> Molecule {
        super::super::parser::parse_sdf(mb).unwrap()
    }

    const ETHANOL_HEAVY: &str = "ethanol\n  test\n\n  3  2  0  0  0  0  0  0  0  0999 V2000\n    1.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    2.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    3.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0\n  2  3  1  0\nM  END";

    #[test]
    fn ethanol_gets_five_hydrogens() {
        let with_h = molblock_with_h(ETHANOL_HEAVY).unwrap();
        eprintln!("MB:\n{}", with_h);
        let m = parse_mb(&with_h);
        eprintln!("parsed {} atoms, {} bonds", m.atoms.len(), m.bonds.len());
        assert_eq!(m.atoms.len(), 9, "CCO + 5H");
        let h_on = |a: usize| {
            m.adjacency[a]
                .iter()
                .filter(|&&j| m.atoms[j].atomic_number == 1)
                .count()
        };
        eprintln!(
            "atoms={}, H_on(C0,C1,O2)=({},{},{})",
            m.atoms.len(),
            h_on(0),
            h_on(1),
            h_on(2)
        );
        assert_eq!(h_on(0), 3, "CH3");
        assert_eq!(h_on(1), 2, "CH2");
        assert_eq!(h_on(2), 1, "OH");
        // O gets exactly one H
        assert_eq!(h_on(2), 1, "OH");
    }

    #[test]
    fn benzene_aromatic_ch_gets_one_h_each() {
        // heavy-only benzene, aromatic bond type 4
        let mut mb = String::from("benzene\n  test\n\n  6  6  0  0  0  0  0  0  0  0999 V2000\n");
        for k in 0..6 {
            let a = k as f64 * std::f64::consts::PI / 3.0;
            mb.push_str(&format!(
                "{:>10.4}{:>10.4}{:>10.4} C   0  0  0  0  0  0  0  0  0  0  0  0\n",
                1.39 * a.cos(),
                1.39 * a.sin(),
                0.0
            ));
        }
        for k in 0..6 {
            mb.push_str(&format!("{:>3}{:>3}  4  0\n", k + 1, (k + 1) % 6 + 1));
        }
        mb.push_str("M  END\n");
        let with_h = molblock_with_h(&mb).unwrap();
        let m = parse_mb(&with_h);
        assert_eq!(m.atoms.len(), 12, "C6 + 6H");
    }

    #[test]
    fn carboxylate_anion_o_gets_no_h() {
        // acetate C(=O)[O-]: M CHG on O3 = -1
        let mb = "acetate\n  test\n\n  5  4  0  0  0  0  0  0  0  0999 V2000\n    1.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    2.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    3.0000    1.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n    3.2000   -1.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n    0.5000   -1.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0\n  2  3  2  0\n  2  4  1  0\n  1  5  1  0\nM  CHG  1   4  -1\nM  END";
        let with_h = molblock_with_h(mb).unwrap();
        let m = parse_mb(&with_h);
        // C1: 2 heavy bonds → 2H; C2 (carboxyl): 4 explicit → 0H; O3 (=O): 0H;
        // O4(-1): 1 bond, valence(O-) = 1 → 0H; C5 methyl: 1 heavy bond → 3H
        assert_eq!(m.atoms.len(), 5 + 5, "acetate: 2(C1) + 3(C5) = 5 H");
        // total charge preserved: -1 on O4
        let q: f64 = m.atoms.iter().map(|a| a.charge).sum();
        assert!((q - (-1.0)).abs() < 1e-6, "charge {}", q);
    }

    #[test]
    fn hydrogens_3d_placement_geometry() {
        // all-atom ethanol (9 atoms) with a 3D conformer: H bond lengths sane,
        // no two H overlapping, every H within bonding distance of its parent
        let text = std::fs::read_to_string("/tmp/caff_heavy.sdf");
        let mol = match text {
            Ok(t) => {
                let m = super::super::parser::parse_sdf(&t).unwrap();
                add_hydrogens(&m)
            }
            Err(_) => {
                let heavy = parse_mb(ETHANOL_HEAVY);
                // give it a fake 3D conformer
                let mut heavy = heavy;
                for (i, a) in heavy.atoms.iter_mut().enumerate() {
                    a.position = [i as f64 * 1.4, 0.0, (i as f64).sin()];
                }
                add_hydrogens(&heavy)
            }
        };
        let hs: Vec<usize> = (0..mol.atoms.len())
            .filter(|&i| mol.atoms[i].atomic_number == 1)
            .collect();
        // bond lengths
        for &h in &hs {
            let &p = mol.adjacency[h].first().unwrap();
            let d = ((mol.atoms[h].position[0] - mol.atoms[p].position[0]).powi(2)
                + (mol.atoms[h].position[1] - mol.atoms[p].position[1]).powi(2)
                + (mol.atoms[h].position[2] - mol.atoms[p].position[2]).powi(2))
            .sqrt();
            assert!((0.9..1.45).contains(&d), "H{}-parent bond length {}", h, d);
        }
        // H-H separation
        for (a, ha) in hs.iter().enumerate() {
            for hb in hs.iter().skip(a + 1) {
                let d = ((mol.atoms[*ha].position[0] - mol.atoms[*hb].position[0]).powi(2)
                    + (mol.atoms[*ha].position[1] - mol.atoms[*hb].position[1]).powi(2)
                    + (mol.atoms[*ha].position[2] - mol.atoms[*hb].position[2]).powi(2))
                .sqrt();
                assert!(d > 0.5, "H{} and H{} too close: {}", ha, hb, d);
            }
        }
    }

    #[test]
    fn explicit_h_input_is_unchanged() {
        let with_h = molblock_with_h(ETHANOL_HEAVY).unwrap();
        let again = molblock_with_h(&with_h).unwrap();
        let m = parse_mb(&again);
        assert_eq!(m.atoms.len(), 9, "no double addition");
    }
}
