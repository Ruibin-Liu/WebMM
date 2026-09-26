//! Gradient audit: analytic vs central-finite-difference per atom across
//! the MMFF validation corpus (scripts/val_set*) and the GFN-FF fixtures,
//! at multiple perturbed geometries. Correctness infrastructure — run with:
//!   cargo run --release --example grad_audit [-- --mmff | --gfnff]
use std::path::PathBuf;

fn corpus_sdfs() -> Vec<PathBuf> {
    let mut out = Vec::new();
    for set in [
        "val_set",
        "val_set_new",
        "val_set_new2",
        "val_set_new3",
        "val_set_new4",
        "val_set_new5",
        "val_set_bulk",
        "val_set_new6",
    ] {
        let dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("scripts")
            .join(set);
        if let Ok(rd) = std::fs::read_dir(&dir) {
            for e in rd.flatten() {
                let p = e.path();
                if p.extension().map(|x| x == "sdf").unwrap_or(false) {
                    out.push(p);
                }
            }
        }
    }
    out.sort();
    out
}

/// Small deterministic LCG so audits are reproducible.
struct Rng(u64);
impl Rng {
    fn next_f64(&mut self) -> f64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (self.0 >> 11) as f64 / (1u64 << 53) as f64
    }
}

fn perturb(coords: &[[f64; 3]], rng: &mut Rng, amp: f64) -> Vec<[f64; 3]> {
    coords
        .iter()
        .map(|c| {
            [
                c[0] + (rng.next_f64() - 0.5) * amp,
                c[1] + (rng.next_f64() - 0.5) * amp,
                c[2] + (rng.next_f64() - 0.5) * amp,
            ]
        })
        .collect()
}

fn audit_mmff() {
    use webmm::forces::ForceField;
    use webmm::mmff::{MMFFForceField, MMFFVariant};
    use webmm::molecule::parser::parse_sdf;

    let eps = 1e-6; // central FD
    let mut worst_overall = 0.0f64;
    let mut worst_mol = String::new();
    let mut n_mols = 0usize;
    let mut failures = Vec::new();

    for sdf_path in corpus_sdfs() {
        let name = sdf_path.file_stem().unwrap().to_string_lossy().to_string();
        let Ok(sdf) = std::fs::read_to_string(&sdf_path) else {
            continue;
        };
        let Ok(mol) = parse_sdf(&sdf) else {
            continue;
        };
        let ff = MMFFForceField::new(&mol, MMFFVariant::MMFF94s);
        let base: Vec<[f64; 3]> = mol.atoms.iter().map(|a| a.position).collect();
        n_mols += 1;

        for trial in 0..3 {
            let mut rng = Rng(0x9E3779B97F4A7C15 ^ ((trial as u64) * 31 + name.len() as u64));
            let coords = perturb(&base, &mut rng, 0.35);
            let n = coords.len();
            let mut g = vec![[0.0f64; 3]; n];
            let _ = ff.energy_and_gradient(&coords, &mut g);

            let mut worst_here = 0.0f64;
            'atoms: for a in 0..n {
                for d in 0..3 {
                    let mut cp = coords.clone();
                    cp[a][d] += eps;
                    let ep = ff.energy(&cp);
                    cp[a][d] -= 2.0 * eps;
                    let em = ff.energy(&cp);
                    let fd = (ep - em) / (2.0 * eps);
                    let an = g[a][d];
                    let diff = (an - fd).abs();
                    // relative to the larger magnitude, with an absolute floor
                    let scale = an.abs().max(fd.abs()).max(1e-3);
                    let rel = diff / scale;
                    if rel > worst_here {
                        worst_here = rel;
                    }
                    // FD truncation O(eps^2 * f''') ~ 1e-8-ish on smooth
                    // surfaces; clamp regions and degenerate guards are
                    // legitimately discontinuous — flag above 1e-4.
                    if rel > 1e-4 && diff > 1e-6 {
                        failures.push(format!(
                            "{name} trial {trial} atom {a} dim {d}: analytic {an:+.6} fd {fd:+.6} (rel {rel:.2e})"
                        ));
                        break 'atoms;
                    }
                }
            }
            if worst_here > worst_overall {
                worst_overall = worst_here;
                worst_mol = name.clone();
            }
        }
    }
    println!(
        "MMFF: {n_mols} molecules x 3 geometries | worst rel err {worst_overall:.3e} ({worst_mol}) | failures {}",
        failures.len()
    );
    for f in failures.iter().take(20) {
        println!("  FAIL {f}");
    }
}

fn audit_gfnff() {
    use webmm::forces::ForceField;
    use webmm::gfnff::GfnffForceField;
    use webmm::molecule::parser::parse_sdf;

    let mut fixtures: Vec<PathBuf> = vec![
        PathBuf::from("tests/fixtures/gfnff/n_methylformamide.mol"),
        PathBuf::from("tests/fixtures/gfnff/thiophene.mol"),
        PathBuf::from("tests/fixtures/gfnff/pyrrole.mol"),
    ];
    for m in [
        "cobalt_ammine",
        "ferricyanide",
        "ferrocene",
        "nicarbonyl",
        "zinc_ammine",
    ] {
        fixtures.push(PathBuf::from(format!(
            "tests/fixtures/gfnff/metals/{m}.mol"
        )));
    }
    let eps = 1e-6;
    let mut failures = Vec::new();
    let mut worst = 0.0f64;
    for f in &fixtures {
        let name = f.file_stem().unwrap().to_string_lossy().to_string();
        let path = if f.exists() {
            f.clone()
        } else {
            PathBuf::from(env!("CARGO_MANIFEST_DIR")).join(f)
        };
        let Ok(text) = std::fs::read_to_string(&path) else {
            println!("  (skip {name}: not found)");
            continue;
        };
        let Ok(mol) = parse_sdf(&text) else {
            println!("  (skip {name}: parse failed)");
            continue;
        };
        let at: Vec<usize> = mol.atoms.iter().map(|a| a.atomic_number as usize).collect();
        let base: Vec<[f64; 3]> = mol.atoms.iter().map(|a| a.position).collect();
        let charge = mol.atoms.iter().map(|a| a.charge).sum::<f64>().round();
        let ff = GfnffForceField::new(&at, &base, charge);
        for trial in 0..2 {
            let mut rng = Rng(0xDEADBEEF ^ ((trial as u64) * 97 + name.len() as u64));
            let coords = perturb(&base, &mut rng, 0.25);
            let n = coords.len();
            let mut g = vec![[0.0f64; 3]; n];
            let _ = ff.energy_and_gradient(&coords, &mut g);
            'atoms: for a in 0..n {
                for d in 0..3 {
                    let mut cp = coords.clone();
                    cp[a][d] += eps;
                    let ep = ff.energy(&cp);
                    cp[a][d] -= 2.0 * eps;
                    let em = ff.energy(&cp);
                    let fd = (ep - em) / (2.0 * eps);
                    let an = g[a][d];
                    let diff = (an - fd).abs();
                    let scale = an.abs().max(fd.abs()).max(1e-3);
                    let rel = diff / scale;
                    if rel > worst {
                        worst = rel;
                    }
                    if rel > 1e-4 && diff > 1e-6 {
                        failures.push(format!(
                            "{name} trial {trial} atom {a} dim {d}: analytic {an:+.6} fd {fd:+.6} (rel {rel:.2e})"
                        ));
                        break 'atoms;
                    }
                }
            }
        }
    }
    println!(
        "GFN-FF: fixtures x 2 geometries | worst rel err {worst:.3e} | failures {}",
        failures.len()
    );
    for f in failures.iter().take(20) {
        println!("  FAIL {f}");
    }
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let do_mmff = !args.contains(&"--gfnff".to_string());
    let do_gfnff = !args.contains(&"--mmff".to_string());
    if do_mmff {
        audit_mmff();
    }
    if do_gfnff {
        audit_gfnff();
    }
}
