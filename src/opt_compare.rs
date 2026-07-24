#[cfg(test)]
mod opt_compare {
    use crate::mmff::{MMFFForceField, MMFFVariant};
    use crate::molecule::{Atom, Bond, BondStereo, BondType, Molecule};
    use crate::optimizer;
    use crate::ConvergenceOptions;

    fn dist(a: &[f64; 3], b: &[f64; 3]) -> f64 {
        ((a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2) + (a[2] - b[2]).powi(2)).sqrt()
    }

    fn bd(name: &str, mol: &Molecule, coords: &[[f64; 3]]) {
        let ff = MMFFForceField::new(mol, MMFFVariant::MMFF94s);
        let b = ff.calculate_energy_breakdown(coords);
        println!("\n  BD {}: bond={:.3} angle={:.3} sb={:.3} tor={:.3} oop={:.3} vdw={:.3} elec={:.3} total={:.3}",
            name, b.bond, b.angle, b.stretch_bend, b.torsion, b.oop, b.vdw, b.electrostatic, b.total());
    }

    fn make_mol(
        name: &str,
        atom_syms: &[&str],
        bonds: &[(usize, usize, u8)],
        coords: &[[f64; 3]],
    ) -> Molecule {
        let anum = |s: &str| -> u8 {
            match s {
                "H" => 1,
                "C" => 6,
                "N" => 7,
                "O" => 8,
                "F" => 9,
                "P" => 15,
                "S" => 16,
                "Cl" => 17,
                "Br" => 35,
                "I" => 53,
                _ => 0,
            }
        };
        let atoms: Vec<Atom> = atom_syms
            .iter()
            .enumerate()
            .map(|(i, s)| Atom {
                symbol: s.to_string(),
                atomic_number: anum(s),
                mass: 0.0,
                charge: 0.0,
                position: coords[i],
                index: i,
                stereo_parity: 0,
            })
            .collect();
        let bonds: Vec<Bond> = bonds
            .iter()
            .map(|(a, b, bo)| Bond {
                atom1: *a,
                atom2: *b,
                bond_type: match bo {
                    20 => BondType::Double,
                    30 => BondType::Triple,
                    15 => BondType::Aromatic,
                    _ => BondType::Single,
                },
                stereo: BondStereo::None,
            })
            .collect();
        let n = atoms.len();
        let mut adj = vec![vec![]; n];
        for b in &bonds {
            adj[b.atom1].push(b.atom2);
            adj[b.atom2].push(b.atom1);
        }
        Molecule {
            atoms,
            bonds,
            name: name.to_string(),
            adjacency: adj,
        }
    }

    #[test]
    fn opt_water() {
        let syms: &[&str] = &["O", "H", "H"];
        let bonds: &[(usize, usize, u8)] = &[(0, 1, 10), (0, 2, 10)];
        let coords: &[[f64; 3]] = &[
            [0.007544, 0.397743, 0.000000],
            [-0.767103, -0.184393, 0.000000],
            [0.759559, -0.213350, 0.000000],
        ];
        let mol = make_mol("water", syms, bonds, coords);
        let ff = MMFFForceField::new(&mol, MMFFVariant::MMFF94s);
        let rdkit_coords: Vec<[f64; 3]> = coords.to_vec();
        let opts = ConvergenceOptions {
            max_iterations: 2000,
            ..Default::default()
        };
        let result = optimizer::optimize(&ff, &rdkit_coords, &opts);
        let rdkit_e = 0.000000_f64;
        let webmm_e = result.final_energy;
        println!("\n=== WATER ===");
        println!(
            "RDKit E: {:.6}  WebMM E: {:.6}  Delta: {:.6}",
            rdkit_e,
            webmm_e,
            webmm_e - rdkit_e
        );
        println!(
            "Converged: {} ({} iters)",
            result.converged, result.iterations
        );
        {
            let dr = dist(&rdkit_coords[0], &rdkit_coords[1]);
            let dw = dist(&result.optimized_coords[0], &result.optimized_coords[1]);
            let err = (dr - dw).abs();
            if err > 0.01 {
                println!("  Bond 0-1: RDKit={:.4} WebMM={:.4} err={:.4}", dr, dw, err);
            }
        }
        {
            let dr = dist(&rdkit_coords[0], &rdkit_coords[2]);
            let dw = dist(&result.optimized_coords[0], &result.optimized_coords[2]);
            let err = (dr - dw).abs();
            if err > 0.01 {
                println!("  Bond 0-2: RDKit={:.4} WebMM={:.4} err={:.4}", dr, dw, err);
            }
        }
    }

    #[test]
    fn opt_methane() {
        let syms: &[&str] = &["C", "H", "H", "H", "H"];
        let bonds: &[(usize, usize, u8)] = &[(0, 1, 10), (0, 2, 10), (0, 3, 10), (0, 4, 10)];
        let coords: &[[f64; 3]] = &[
            [-0.000000, -0.000000, -0.000000],
            [-0.675343, 0.854123, -0.085361],
            [-0.394946, -0.835938, -0.581485],
            [0.084754, -0.293741, 1.048538],
            [0.985535, 0.275556, -0.381692],
        ];
        let mol = make_mol("methane", syms, bonds, coords);
        let ff = MMFFForceField::new(&mol, MMFFVariant::MMFF94s);
        let rdkit_coords: Vec<[f64; 3]> = coords.to_vec();
        let opts = ConvergenceOptions {
            max_iterations: 2000,
            ..Default::default()
        };
        let result = optimizer::optimize(&ff, &rdkit_coords, &opts);
        let rdkit_e = 0.026383_f64;
        let webmm_e = result.final_energy;
        println!("\n=== METHANE ===");
        println!(
            "RDKit E: {:.6}  WebMM E: {:.6}  Delta: {:.6}",
            rdkit_e,
            webmm_e,
            webmm_e - rdkit_e
        );
        println!(
            "Converged: {} ({} iters)",
            result.converged, result.iterations
        );
        {
            let dr = dist(&rdkit_coords[0], &rdkit_coords[1]);
            let dw = dist(&result.optimized_coords[0], &result.optimized_coords[1]);
            let err = (dr - dw).abs();
            if err > 0.01 {
                println!("  Bond 0-1: RDKit={:.4} WebMM={:.4} err={:.4}", dr, dw, err);
            }
        }
        {
            let dr = dist(&rdkit_coords[0], &rdkit_coords[2]);
            let dw = dist(&result.optimized_coords[0], &result.optimized_coords[2]);
            let err = (dr - dw).abs();
            if err > 0.01 {
                println!("  Bond 0-2: RDKit={:.4} WebMM={:.4} err={:.4}", dr, dw, err);
            }
        }
        {
            let dr = dist(&rdkit_coords[0], &rdkit_coords[3]);
            let dw = dist(&result.optimized_coords[0], &result.optimized_coords[3]);
            let err = (dr - dw).abs();
            if err > 0.01 {
                println!("  Bond 0-3: RDKit={:.4} WebMM={:.4} err={:.4}", dr, dw, err);
            }
        }
        {
            let dr = dist(&rdkit_coords[0], &rdkit_coords[4]);
            let dw = dist(&result.optimized_coords[0], &result.optimized_coords[4]);
            let err = (dr - dw).abs();
            if err > 0.01 {
                println!("  Bond 0-4: RDKit={:.4} WebMM={:.4} err={:.4}", dr, dw, err);
            }
        }
    }

    #[test]
    fn opt_ammonia() {
        let syms: &[&str] = &["N", "H", "H", "H"];
        let bonds: &[(usize, usize, u8)] = &[(0, 1, 10), (0, 2, 10), (0, 3, 10)];
        let coords: &[[f64; 3]] = &[
            [0.004320, 0.004533, 0.295529],
            [0.917090, -0.199551, -0.108898],
            [-0.632036, -0.697914, -0.078610],
            [-0.289374, 0.892932, -0.108021],
        ];
        let mol = make_mol("ammonia", syms, bonds, coords);
        let ff = MMFFForceField::new(&mol, MMFFVariant::MMFF94s);
        let rdkit_coords: Vec<[f64; 3]> = coords.to_vec();
        let opts = ConvergenceOptions {
            max_iterations: 2000,
            ..Default::default()
        };
        let result = optimizer::optimize(&ff, &rdkit_coords, &opts);
        let rdkit_e = 0.000000_f64;
        let webmm_e = result.final_energy;
        println!("\n=== AMMONIA ===");
        println!(
            "RDKit E: {:.6}  WebMM E: {:.6}  Delta: {:.6}",
            rdkit_e,
            webmm_e,
            webmm_e - rdkit_e
        );
        println!(
            "Converged: {} ({} iters)",
            result.converged, result.iterations
        );
        {
            let dr = dist(&rdkit_coords[0], &rdkit_coords[1]);
            let dw = dist(&result.optimized_coords[0], &result.optimized_coords[1]);
            let err = (dr - dw).abs();
            if err > 0.01 {
                println!("  Bond 0-1: RDKit={:.4} WebMM={:.4} err={:.4}", dr, dw, err);
            }
        }
        {
            let dr = dist(&rdkit_coords[0], &rdkit_coords[2]);
            let dw = dist(&result.optimized_coords[0], &result.optimized_coords[2]);
            let err = (dr - dw).abs();
            if err > 0.01 {
                println!("  Bond 0-2: RDKit={:.4} WebMM={:.4} err={:.4}", dr, dw, err);
            }
        }
        {
            let dr = dist(&rdkit_coords[0], &rdkit_coords[3]);
            let dw = dist(&result.optimized_coords[0], &result.optimized_coords[3]);
            let err = (dr - dw).abs();
            if err > 0.01 {
                println!("  Bond 0-3: RDKit={:.4} WebMM={:.4} err={:.4}", dr, dw, err);
            }
        }
    }

    #[test]
    fn opt_ethane() {
        let syms: &[&str] = &["C", "C", "H", "H", "H", "H", "H", "H"];
        let bonds: &[(usize, usize, u8)] = &[
            (0, 1, 10),
            (0, 2, 10),
            (0, 3, 10),
            (0, 4, 10),
            (1, 5, 10),
            (1, 6, 10),
            (1, 7, 10),
        ];
        let coords: &[[f64; 3]] = &[
            [-0.755815, 0.007100, -0.016504],
            [0.755816, -0.007100, 0.016504],
            [-1.163351, -0.100381, 0.993141],
            [-1.122263, 0.948089, -0.437543],
            [-1.134621, -0.815581, -0.630281],
            [1.134621, 0.815583, 0.630279],
            [1.163351, 0.100379, -0.993141],
            [1.122263, -0.948088, 0.437545],
        ];
        let mol = make_mol("ethane", syms, bonds, coords);
        let ff = MMFFForceField::new(&mol, MMFFVariant::MMFF94s);
        let rdkit_coords: Vec<[f64; 3]> = coords.to_vec();
        let opts = ConvergenceOptions {
            max_iterations: 2000,
            ..Default::default()
        };
        let result = optimizer::optimize(&ff, &rdkit_coords, &opts);
        let rdkit_e = -4.734365_f64;
        let webmm_e = result.final_energy;
        println!("\n=== ETHANE ===");
        println!(
            "RDKit E: {:.6}  WebMM E: {:.6}  Delta: {:.6}",
            rdkit_e,
            webmm_e,
            webmm_e - rdkit_e
        );
        println!(
            "Converged: {} ({} iters)",
            result.converged, result.iterations
        );
        {
            let dr = dist(&rdkit_coords[0], &rdkit_coords[1]);
            let dw = dist(&result.optimized_coords[0], &result.optimized_coords[1]);
            let err = (dr - dw).abs();
            if err > 0.01 {
                println!("  Bond 0-1: RDKit={:.4} WebMM={:.4} err={:.4}", dr, dw, err);
            }
        }
        {
            let dr = dist(&rdkit_coords[0], &rdkit_coords[2]);
            let dw = dist(&result.optimized_coords[0], &result.optimized_coords[2]);
            let err = (dr - dw).abs();
            if err > 0.01 {
                println!("  Bond 0-2: RDKit={:.4} WebMM={:.4} err={:.4}", dr, dw, err);
            }
        }
        {
            let dr = dist(&rdkit_coords[0], &rdkit_coords[3]);
            let dw = dist(&result.optimized_coords[0], &result.optimized_coords[3]);
            let err = (dr - dw).abs();
            if err > 0.01 {
                println!("  Bond 0-3: RDKit={:.4} WebMM={:.4} err={:.4}", dr, dw, err);
            }
        }
        {
            let dr = dist(&rdkit_coords[0], &rdkit_coords[4]);
            let dw = dist(&result.optimized_coords[0], &result.optimized_coords[4]);
            let err = (dr - dw).abs();
            if err > 0.01 {
                println!("  Bond 0-4: RDKit={:.4} WebMM={:.4} err={:.4}", dr, dw, err);
            }
        }
        {
            let dr = dist(&rdkit_coords[1], &rdkit_coords[5]);
            let dw = dist(&result.optimized_coords[1], &result.optimized_coords[5]);
            let err = (dr - dw).abs();
            if err > 0.01 {
                println!("  Bond 1-5: RDKit={:.4} WebMM={:.4} err={:.4}", dr, dw, err);
            }
        }
        {
            let dr = dist(&rdkit_coords[1], &rdkit_coords[6]);
            let dw = dist(&result.optimized_coords[1], &result.optimized_coords[6]);
            let err = (dr - dw).abs();
            if err > 0.01 {
                println!("  Bond 1-6: RDKit={:.4} WebMM={:.4} err={:.4}", dr, dw, err);
            }
        }
        {
            let dr = dist(&rdkit_coords[1], &rdkit_coords[7]);
            let dw = dist(&result.optimized_coords[1], &result.optimized_coords[7]);
            let err = (dr - dw).abs();
            if err > 0.01 {
                println!("  Bond 1-7: RDKit={:.4} WebMM={:.4} err={:.4}", dr, dw, err);
            }
        }
    }

    #[test]
    fn opt_ethanol() {
        let syms: &[&str] = &["C", "C", "O", "H", "H", "H", "H", "H", "H"];
        let bonds: &[(usize, usize, u8)] = &[
            (0, 1, 10),
            (1, 2, 10),
            (0, 3, 10),
            (0, 4, 10),
            (0, 5, 10),
            (1, 6, 10),
            (1, 7, 10),
            (2, 8, 10),
        ];
        let coords: &[[f64; 3]] = &[
            [-0.888311, 0.167003, -0.027316],
            [0.465753, -0.511559, -0.036795],
            [1.431075, 0.322916, 0.586670],
            [-0.848741, 1.117480, -0.569524],
            [-1.647121, -0.470443, -0.489637],
            [-1.196397, 0.397845, 0.997723],
            [0.791997, -0.722428, -1.059726],
            [0.424604, -1.455862, 0.513791],
            [1.467142, 1.155048, 0.084814],
        ];
        let mol = make_mol("ethanol", syms, bonds, coords);
        let ff = MMFFForceField::new(&mol, MMFFVariant::MMFF94s);
        let rdkit_coords: Vec<[f64; 3]> = coords.to_vec();
        let opts = ConvergenceOptions {
            max_iterations: 2000,
            ..Default::default()
        };
        let result = optimizer::optimize(&ff, &rdkit_coords, &opts);
        let rdkit_e = -1.336857_f64;
        let webmm_e = result.final_energy;
        println!("\n=== ETHANOL ===");
        println!(
            "RDKit E: {:.6}  WebMM E: {:.6}  Delta: {:.6}",
            rdkit_e,
            webmm_e,
            webmm_e - rdkit_e
        );
        println!(
            "Converged: {} ({} iters)",
            result.converged, result.iterations
        );
        {
            let dr = dist(&rdkit_coords[0], &rdkit_coords[1]);
            let dw = dist(&result.optimized_coords[0], &result.optimized_coords[1]);
            let err = (dr - dw).abs();
            if err > 0.01 {
                println!("  Bond 0-1: RDKit={:.4} WebMM={:.4} err={:.4}", dr, dw, err);
            }
        }
        {
            let dr = dist(&rdkit_coords[1], &rdkit_coords[2]);
            let dw = dist(&result.optimized_coords[1], &result.optimized_coords[2]);
            let err = (dr - dw).abs();
            if err > 0.01 {
                println!("  Bond 1-2: RDKit={:.4} WebMM={:.4} err={:.4}", dr, dw, err);
            }
        }
        {
            let dr = dist(&rdkit_coords[0], &rdkit_coords[3]);
            let dw = dist(&result.optimized_coords[0], &result.optimized_coords[3]);
            let err = (dr - dw).abs();
            if err > 0.01 {
                println!("  Bond 0-3: RDKit={:.4} WebMM={:.4} err={:.4}", dr, dw, err);
            }
        }
        {
            let dr = dist(&rdkit_coords[0], &rdkit_coords[4]);
            let dw = dist(&result.optimized_coords[0], &result.optimized_coords[4]);
            let err = (dr - dw).abs();
            if err > 0.01 {
                println!("  Bond 0-4: RDKit={:.4} WebMM={:.4} err={:.4}", dr, dw, err);
            }
        }
        {
            let dr = dist(&rdkit_coords[0], &rdkit_coords[5]);
            let dw = dist(&result.optimized_coords[0], &result.optimized_coords[5]);
            let err = (dr - dw).abs();
            if err > 0.01 {
                println!("  Bond 0-5: RDKit={:.4} WebMM={:.4} err={:.4}", dr, dw, err);
            }
        }
        {
            let dr = dist(&rdkit_coords[1], &rdkit_coords[6]);
            let dw = dist(&result.optimized_coords[1], &result.optimized_coords[6]);
            let err = (dr - dw).abs();
            if err > 0.01 {
                println!("  Bond 1-6: RDKit={:.4} WebMM={:.4} err={:.4}", dr, dw, err);
            }
        }
        {
            let dr = dist(&rdkit_coords[1], &rdkit_coords[7]);
            let dw = dist(&result.optimized_coords[1], &result.optimized_coords[7]);
            let err = (dr - dw).abs();
            if err > 0.01 {
                println!("  Bond 1-7: RDKit={:.4} WebMM={:.4} err={:.4}", dr, dw, err);
            }
        }
        {
            let dr = dist(&rdkit_coords[2], &rdkit_coords[8]);
            let dw = dist(&result.optimized_coords[2], &result.optimized_coords[8]);
            let err = (dr - dw).abs();
            if err > 0.01 {
                println!("  Bond 2-8: RDKit={:.4} WebMM={:.4} err={:.4}", dr, dw, err);
            }
        }
    }

    #[test]
    fn opt_formaldehyde() {
        let syms: &[&str] = &["C", "O", "H", "H"];
        let bonds: &[(usize, usize, u8)] = &[(0, 1, 20), (0, 2, 10), (0, 3, 10)];
        let coords: &[[f64; 3]] = &[
            [-0.012188, 0.001876, 0.000183],
            [1.198062, -0.184373, -0.018005],
            [-0.451296, 1.012282, 0.002316],
            [-0.734578, -0.829785, 0.015506],
        ];
        let mol = make_mol("formaldehyde", syms, bonds, coords);
        let ff = MMFFForceField::new(&mol, MMFFVariant::MMFF94s);
        let rdkit_coords: Vec<[f64; 3]> = coords.to_vec();
        let opts = ConvergenceOptions {
            max_iterations: 2000,
            ..Default::default()
        };
        let result = optimizer::optimize(&ff, &rdkit_coords, &opts);
        let rdkit_e = 0.054161_f64;
        let webmm_e = result.final_energy;
        println!("\n=== FORMALDEHYDE ===");
        println!(
            "RDKit E: {:.6}  WebMM E: {:.6}  Delta: {:.6}",
            rdkit_e,
            webmm_e,
            webmm_e - rdkit_e
        );
        println!(
            "Converged: {} ({} iters)",
            result.converged, result.iterations
        );
        {
            let dr = dist(&rdkit_coords[0], &rdkit_coords[1]);
            let dw = dist(&result.optimized_coords[0], &result.optimized_coords[1]);
            let err = (dr - dw).abs();
            if err > 0.01 {
                println!("  Bond 0-1: RDKit={:.4} WebMM={:.4} err={:.4}", dr, dw, err);
            }
        }
        {
            let dr = dist(&rdkit_coords[0], &rdkit_coords[2]);
            let dw = dist(&result.optimized_coords[0], &result.optimized_coords[2]);
            let err = (dr - dw).abs();
            if err > 0.01 {
                println!("  Bond 0-2: RDKit={:.4} WebMM={:.4} err={:.4}", dr, dw, err);
            }
        }
        {
            let dr = dist(&rdkit_coords[0], &rdkit_coords[3]);
            let dw = dist(&result.optimized_coords[0], &result.optimized_coords[3]);
            let err = (dr - dw).abs();
            if err > 0.01 {
                println!("  Bond 0-3: RDKit={:.4} WebMM={:.4} err={:.4}", dr, dw, err);
            }
        }
    }

    #[test]
    fn opt_acetic_acid() {
        let syms: &[&str] = &["C", "C", "O", "O", "H", "H", "H", "H"];
        let bonds: &[(usize, usize, u8)] = &[
            (0, 1, 10),
            (1, 2, 20),
            (1, 3, 10),
            (0, 4, 10),
            (0, 5, 10),
            (0, 6, 10),
            (3, 7, 10),
        ];
        let coords: &[[f64; 3]] = &[
            [-0.933543, -0.060063, -0.230376],
            [0.493608, 0.278933, 0.046914],
            [1.032525, 1.356629, -0.136136],
            [1.181401, -0.764470, 0.546212],
            [-1.442687, 0.820299, -0.632706],
            [-1.430457, -0.357583, 0.696271],
            [-0.986731, -0.862547, -0.970180],
            [2.085884, -0.411198, 0.680001],
        ];
        let mol = make_mol("acetic_acid", syms, bonds, coords);
        let ff = MMFFForceField::new(&mol, MMFFVariant::MMFF94s);
        let rdkit_coords: Vec<[f64; 3]> = coords.to_vec();
        let opts = ConvergenceOptions {
            max_iterations: 2000,
            ..Default::default()
        };
        let result = optimizer::optimize(&ff, &rdkit_coords, &opts);
        let rdkit_e = -26.405933_f64;
        let webmm_e = result.final_energy;
        println!("\n=== ACETIC_ACID ===");
        println!(
            "RDKit E: {:.6}  WebMM E: {:.6}  Delta: {:.6}",
            rdkit_e,
            webmm_e,
            webmm_e - rdkit_e
        );
        println!(
            "Converged: {} ({} iters)",
            result.converged, result.iterations
        );
        {
            let dr = dist(&rdkit_coords[0], &rdkit_coords[1]);
            let dw = dist(&result.optimized_coords[0], &result.optimized_coords[1]);
            let err = (dr - dw).abs();
            if err > 0.01 {
                println!("  Bond 0-1: RDKit={:.4} WebMM={:.4} err={:.4}", dr, dw, err);
            }
        }
        {
            let dr = dist(&rdkit_coords[1], &rdkit_coords[2]);
            let dw = dist(&result.optimized_coords[1], &result.optimized_coords[2]);
            let err = (dr - dw).abs();
            if err > 0.01 {
                println!("  Bond 1-2: RDKit={:.4} WebMM={:.4} err={:.4}", dr, dw, err);
            }
        }
        {
            let dr = dist(&rdkit_coords[1], &rdkit_coords[3]);
            let dw = dist(&result.optimized_coords[1], &result.optimized_coords[3]);
            let err = (dr - dw).abs();
            if err > 0.01 {
                println!("  Bond 1-3: RDKit={:.4} WebMM={:.4} err={:.4}", dr, dw, err);
            }
        }
        {
            let dr = dist(&rdkit_coords[0], &rdkit_coords[4]);
            let dw = dist(&result.optimized_coords[0], &result.optimized_coords[4]);
            let err = (dr - dw).abs();
            if err > 0.01 {
                println!("  Bond 0-4: RDKit={:.4} WebMM={:.4} err={:.4}", dr, dw, err);
            }
        }
        {
            let dr = dist(&rdkit_coords[0], &rdkit_coords[5]);
            let dw = dist(&result.optimized_coords[0], &result.optimized_coords[5]);
            let err = (dr - dw).abs();
            if err > 0.01 {
                println!("  Bond 0-5: RDKit={:.4} WebMM={:.4} err={:.4}", dr, dw, err);
            }
        }
        {
            let dr = dist(&rdkit_coords[0], &rdkit_coords[6]);
            let dw = dist(&result.optimized_coords[0], &result.optimized_coords[6]);
            let err = (dr - dw).abs();
            if err > 0.01 {
                println!("  Bond 0-6: RDKit={:.4} WebMM={:.4} err={:.4}", dr, dw, err);
            }
        }
        {
            let dr = dist(&rdkit_coords[3], &rdkit_coords[7]);
            let dw = dist(&result.optimized_coords[3], &result.optimized_coords[7]);
            let err = (dr - dw).abs();
            if err > 0.01 {
                println!("  Bond 3-7: RDKit={:.4} WebMM={:.4} err={:.4}", dr, dw, err);
            }
        }
    }

    #[test]
    fn opt_acetone() {
        let syms: &[&str] = &["C", "C", "O", "C", "H", "H", "H", "H", "H", "H"];
        let bonds: &[(usize, usize, u8)] = &[
            (0, 1, 10),
            (1, 2, 20),
            (1, 3, 10),
            (0, 4, 10),
            (0, 5, 10),
            (0, 6, 10),
            (3, 7, 10),
            (3, 8, 10),
            (3, 9, 10),
        ];
        let coords: &[[f64; 3]] = &[
            [-1.265992, 0.241506, 0.092371],
            [0.068495, 0.155258, -0.597576],
            [0.204122, 0.462684, -1.780855],
            [1.228925, -0.325526, 0.231017],
            [-2.022220, 0.601195, -0.610859],
            [-1.200942, 0.940754, 0.929593],
            [-1.558302, -0.748971, 0.449621],
            [1.026390, -1.336406, 0.593257],
            [2.135769, -0.343814, -0.379795],
            [1.383755, 0.353319, 1.073226],
        ];
        let mol = make_mol("acetone", syms, bonds, coords);
        let ff = MMFFForceField::new(&mol, MMFFVariant::MMFF94s);
        let rdkit_coords: Vec<[f64; 3]> = coords.to_vec();
        let opts = ConvergenceOptions {
            max_iterations: 2000,
            ..Default::default()
        };
        let result = optimizer::optimize(&ff, &rdkit_coords, &opts);
        let rdkit_e = 1.652258_f64;
        let webmm_e = result.final_energy;
        println!("\n=== ACETONE ===");
        println!(
            "RDKit E: {:.6}  WebMM E: {:.6}  Delta: {:.6}",
            rdkit_e,
            webmm_e,
            webmm_e - rdkit_e
        );
        println!(
            "Converged: {} ({} iters)",
            result.converged, result.iterations
        );
        {
            let dr = dist(&rdkit_coords[0], &rdkit_coords[1]);
            let dw = dist(&result.optimized_coords[0], &result.optimized_coords[1]);
            let err = (dr - dw).abs();
            if err > 0.01 {
                println!("  Bond 0-1: RDKit={:.4} WebMM={:.4} err={:.4}", dr, dw, err);
            }
        }
        {
            let dr = dist(&rdkit_coords[1], &rdkit_coords[2]);
            let dw = dist(&result.optimized_coords[1], &result.optimized_coords[2]);
            let err = (dr - dw).abs();
            if err > 0.01 {
                println!("  Bond 1-2: RDKit={:.4} WebMM={:.4} err={:.4}", dr, dw, err);
            }
        }
        {
            let dr = dist(&rdkit_coords[1], &rdkit_coords[3]);
            let dw = dist(&result.optimized_coords[1], &result.optimized_coords[3]);
            let err = (dr - dw).abs();
            if err > 0.01 {
                println!("  Bond 1-3: RDKit={:.4} WebMM={:.4} err={:.4}", dr, dw, err);
            }
        }
        {
            let dr = dist(&rdkit_coords[0], &rdkit_coords[4]);
            let dw = dist(&result.optimized_coords[0], &result.optimized_coords[4]);
            let err = (dr - dw).abs();
            if err > 0.01 {
                println!("  Bond 0-4: RDKit={:.4} WebMM={:.4} err={:.4}", dr, dw, err);
            }
        }
        {
            let dr = dist(&rdkit_coords[0], &rdkit_coords[5]);
            let dw = dist(&result.optimized_coords[0], &result.optimized_coords[5]);
            let err = (dr - dw).abs();
            if err > 0.01 {
                println!("  Bond 0-5: RDKit={:.4} WebMM={:.4} err={:.4}", dr, dw, err);
            }
        }
        {
            let dr = dist(&rdkit_coords[0], &rdkit_coords[6]);
            let dw = dist(&result.optimized_coords[0], &result.optimized_coords[6]);
            let err = (dr - dw).abs();
            if err > 0.01 {
                println!("  Bond 0-6: RDKit={:.4} WebMM={:.4} err={:.4}", dr, dw, err);
            }
        }
        {
            let dr = dist(&rdkit_coords[3], &rdkit_coords[7]);
            let dw = dist(&result.optimized_coords[3], &result.optimized_coords[7]);
            let err = (dr - dw).abs();
            if err > 0.01 {
                println!("  Bond 3-7: RDKit={:.4} WebMM={:.4} err={:.4}", dr, dw, err);
            }
        }
        {
            let dr = dist(&rdkit_coords[3], &rdkit_coords[8]);
            let dw = dist(&result.optimized_coords[3], &result.optimized_coords[8]);
            let err = (dr - dw).abs();
            if err > 0.01 {
                println!("  Bond 3-8: RDKit={:.4} WebMM={:.4} err={:.4}", dr, dw, err);
            }
        }
        {
            let dr = dist(&rdkit_coords[3], &rdkit_coords[9]);
            let dw = dist(&result.optimized_coords[3], &result.optimized_coords[9]);
            let err = (dr - dw).abs();
            if err > 0.01 {
                println!("  Bond 3-9: RDKit={:.4} WebMM={:.4} err={:.4}", dr, dw, err);
            }
        }
    }

    #[test]
    fn opt_methylamine() {
        let syms: &[&str] = &["C", "N", "H", "H", "H", "H", "H"];
        let bonds: &[(usize, usize, u8)] = &[
            (0, 1, 10),
            (0, 2, 10),
            (0, 3, 10),
            (0, 4, 10),
            (1, 5, 10),
            (1, 6, 10),
        ];
        let coords: &[[f64; 3]] = &[
            [-0.574373, 0.003792, 0.015807],
            [0.840454, 0.165146, 0.299415],
            [-1.052213, -0.587049, 0.802165],
            [-0.719360, -0.504165, -0.941862],
            [-1.064694, 0.980272, -0.027267],
            [1.290827, -0.749083, 0.306947],
            [1.279359, 0.691087, -0.455206],
        ];
        let mol = make_mol("methylamine", syms, bonds, coords);
        let ff = MMFFForceField::new(&mol, MMFFVariant::MMFF94s);
        let rdkit_coords: Vec<[f64; 3]> = coords.to_vec();
        let opts = ConvergenceOptions {
            max_iterations: 2000,
            ..Default::default()
        };
        let result = optimizer::optimize(&ff, &rdkit_coords, &opts);
        let rdkit_e = -1.694530_f64;
        let webmm_e = result.final_energy;
        println!("\n=== METHYLAMINE ===");
        println!(
            "RDKit E: {:.6}  WebMM E: {:.6}  Delta: {:.6}",
            rdkit_e,
            webmm_e,
            webmm_e - rdkit_e
        );
        println!(
            "Converged: {} ({} iters)",
            result.converged, result.iterations
        );
        {
            let dr = dist(&rdkit_coords[0], &rdkit_coords[1]);
            let dw = dist(&result.optimized_coords[0], &result.optimized_coords[1]);
            let err = (dr - dw).abs();
            if err > 0.01 {
                println!("  Bond 0-1: RDKit={:.4} WebMM={:.4} err={:.4}", dr, dw, err);
            }
        }
        {
            let dr = dist(&rdkit_coords[0], &rdkit_coords[2]);
            let dw = dist(&result.optimized_coords[0], &result.optimized_coords[2]);
            let err = (dr - dw).abs();
            if err > 0.01 {
                println!("  Bond 0-2: RDKit={:.4} WebMM={:.4} err={:.4}", dr, dw, err);
            }
        }
        {
            let dr = dist(&rdkit_coords[0], &rdkit_coords[3]);
            let dw = dist(&result.optimized_coords[0], &result.optimized_coords[3]);
            let err = (dr - dw).abs();
            if err > 0.01 {
                println!("  Bond 0-3: RDKit={:.4} WebMM={:.4} err={:.4}", dr, dw, err);
            }
        }
        {
            let dr = dist(&rdkit_coords[0], &rdkit_coords[4]);
            let dw = dist(&result.optimized_coords[0], &result.optimized_coords[4]);
            let err = (dr - dw).abs();
            if err > 0.01 {
                println!("  Bond 0-4: RDKit={:.4} WebMM={:.4} err={:.4}", dr, dw, err);
            }
        }
        {
            let dr = dist(&rdkit_coords[1], &rdkit_coords[5]);
            let dw = dist(&result.optimized_coords[1], &result.optimized_coords[5]);
            let err = (dr - dw).abs();
            if err > 0.01 {
                println!("  Bond 1-5: RDKit={:.4} WebMM={:.4} err={:.4}", dr, dw, err);
            }
        }
        {
            let dr = dist(&rdkit_coords[1], &rdkit_coords[6]);
            let dw = dist(&result.optimized_coords[1], &result.optimized_coords[6]);
            let err = (dr - dw).abs();
            if err > 0.01 {
                println!("  Bond 1-6: RDKit={:.4} WebMM={:.4} err={:.4}", dr, dw, err);
            }
        }
    }

    #[test]
    fn opt_dimethyl_ether() {
        let syms: &[&str] = &["C", "O", "C", "H", "H", "H", "H", "H", "H"];
        let bonds: &[(usize, usize, u8)] = &[
            (0, 1, 10),
            (1, 2, 10),
            (0, 3, 10),
            (0, 4, 10),
            (0, 5, 10),
            (2, 6, 10),
            (2, 7, 10),
            (2, 8, 10),
        ];
        let coords: &[[f64; 3]] = &[
            [-1.168191, -0.002702, -0.404053],
            [-0.026714, 0.839290, -0.315816],
            [1.165550, 0.085665, -0.141863],
            [-2.051827, 0.627311, -0.536561],
            [-1.084258, -0.672845, -1.265092],
            [-1.287694, -0.583654, 0.515588],
            [1.121349, -0.492440, 0.786234],
            [2.007003, 0.780997, -0.080557],
            [1.324783, -0.581623, -0.994445],
        ];
        let mol = make_mol("dimethyl_ether", syms, bonds, coords);
        let ff = MMFFForceField::new(&mol, MMFFVariant::MMFF94s);
        let rdkit_coords: Vec<[f64; 3]> = coords.to_vec();
        let opts = ConvergenceOptions {
            max_iterations: 2000,
            ..Default::default()
        };
        let result = optimizer::optimize(&ff, &rdkit_coords, &opts);
        let rdkit_e = 6.089114_f64;
        let webmm_e = result.final_energy;
        println!("\n=== DIMETHYL_ETHER ===");
        println!(
            "RDKit E: {:.6}  WebMM E: {:.6}  Delta: {:.6}",
            rdkit_e,
            webmm_e,
            webmm_e - rdkit_e
        );
        println!(
            "Converged: {} ({} iters)",
            result.converged, result.iterations
        );
        {
            let dr = dist(&rdkit_coords[0], &rdkit_coords[1]);
            let dw = dist(&result.optimized_coords[0], &result.optimized_coords[1]);
            let err = (dr - dw).abs();
            if err > 0.01 {
                println!("  Bond 0-1: RDKit={:.4} WebMM={:.4} err={:.4}", dr, dw, err);
            }
        }
        {
            let dr = dist(&rdkit_coords[1], &rdkit_coords[2]);
            let dw = dist(&result.optimized_coords[1], &result.optimized_coords[2]);
            let err = (dr - dw).abs();
            if err > 0.01 {
                println!("  Bond 1-2: RDKit={:.4} WebMM={:.4} err={:.4}", dr, dw, err);
            }
        }
        {
            let dr = dist(&rdkit_coords[0], &rdkit_coords[3]);
            let dw = dist(&result.optimized_coords[0], &result.optimized_coords[3]);
            let err = (dr - dw).abs();
            if err > 0.01 {
                println!("  Bond 0-3: RDKit={:.4} WebMM={:.4} err={:.4}", dr, dw, err);
            }
        }
        {
            let dr = dist(&rdkit_coords[0], &rdkit_coords[4]);
            let dw = dist(&result.optimized_coords[0], &result.optimized_coords[4]);
            let err = (dr - dw).abs();
            if err > 0.01 {
                println!("  Bond 0-4: RDKit={:.4} WebMM={:.4} err={:.4}", dr, dw, err);
            }
        }
        {
            let dr = dist(&rdkit_coords[0], &rdkit_coords[5]);
            let dw = dist(&result.optimized_coords[0], &result.optimized_coords[5]);
            let err = (dr - dw).abs();
            if err > 0.01 {
                println!("  Bond 0-5: RDKit={:.4} WebMM={:.4} err={:.4}", dr, dw, err);
            }
        }
        {
            let dr = dist(&rdkit_coords[2], &rdkit_coords[6]);
            let dw = dist(&result.optimized_coords[2], &result.optimized_coords[6]);
            let err = (dr - dw).abs();
            if err > 0.01 {
                println!("  Bond 2-6: RDKit={:.4} WebMM={:.4} err={:.4}", dr, dw, err);
            }
        }
        {
            let dr = dist(&rdkit_coords[2], &rdkit_coords[7]);
            let dw = dist(&result.optimized_coords[2], &result.optimized_coords[7]);
            let err = (dr - dw).abs();
            if err > 0.01 {
                println!("  Bond 2-7: RDKit={:.4} WebMM={:.4} err={:.4}", dr, dw, err);
            }
        }
        {
            let dr = dist(&rdkit_coords[2], &rdkit_coords[8]);
            let dw = dist(&result.optimized_coords[2], &result.optimized_coords[8]);
            let err = (dr - dw).abs();
            if err > 0.01 {
                println!("  Bond 2-8: RDKit={:.4} WebMM={:.4} err={:.4}", dr, dw, err);
            }
        }
    }

    #[test]
    fn opt_propane() {
        let syms: &[&str] = &["C", "C", "C", "H", "H", "H", "H", "H", "H", "H", "H"];
        let bonds: &[(usize, usize, u8)] = &[
            (0, 1, 10),
            (1, 2, 10),
            (0, 3, 10),
            (0, 4, 10),
            (0, 5, 10),
            (1, 6, 10),
            (1, 7, 10),
            (2, 8, 10),
            (2, 9, 10),
            (2, 10, 10),
        ];
        let coords: &[[f64; 3]] = &[
            [1.210836, -0.174866, -0.371642],
            [0.033860, -0.189097, 0.589031],
            [-1.236390, 0.317571, -0.072876],
            [1.411933, 0.839536, -0.730748],
            [2.112380, -0.543514, 0.127572],
            [1.017132, -0.814150, -1.238939],
            [-0.126974, -1.209651, 0.953738],
            [0.265678, 0.435035, 1.459162],
            [-2.068373, 0.297751, 0.637974],
            [-1.507443, -0.306150, -0.930730],
            [-1.112641, 1.347537, -0.422541],
        ];
        let mol = make_mol("propane", syms, bonds, coords);
        let ff = MMFFForceField::new(&mol, MMFFVariant::MMFF94s);
        let rdkit_coords: Vec<[f64; 3]> = coords.to_vec();
        let opts = ConvergenceOptions {
            max_iterations: 2000,
            ..Default::default()
        };
        let result = optimizer::optimize(&ff, &rdkit_coords, &opts);
        let rdkit_e = -4.897293_f64;
        let webmm_e = result.final_energy;
        println!("\n=== PROPANE ===");
        println!(
            "RDKit E: {:.6}  WebMM E: {:.6}  Delta: {:.6}",
            rdkit_e,
            webmm_e,
            webmm_e - rdkit_e
        );
        println!(
            "Converged: {} ({} iters)",
            result.converged, result.iterations
        );
        {
            let dr = dist(&rdkit_coords[0], &rdkit_coords[1]);
            let dw = dist(&result.optimized_coords[0], &result.optimized_coords[1]);
            let err = (dr - dw).abs();
            if err > 0.01 {
                println!("  Bond 0-1: RDKit={:.4} WebMM={:.4} err={:.4}", dr, dw, err);
            }
        }
        {
            let dr = dist(&rdkit_coords[1], &rdkit_coords[2]);
            let dw = dist(&result.optimized_coords[1], &result.optimized_coords[2]);
            let err = (dr - dw).abs();
            if err > 0.01 {
                println!("  Bond 1-2: RDKit={:.4} WebMM={:.4} err={:.4}", dr, dw, err);
            }
        }
        {
            let dr = dist(&rdkit_coords[0], &rdkit_coords[3]);
            let dw = dist(&result.optimized_coords[0], &result.optimized_coords[3]);
            let err = (dr - dw).abs();
            if err > 0.01 {
                println!("  Bond 0-3: RDKit={:.4} WebMM={:.4} err={:.4}", dr, dw, err);
            }
        }
        {
            let dr = dist(&rdkit_coords[0], &rdkit_coords[4]);
            let dw = dist(&result.optimized_coords[0], &result.optimized_coords[4]);
            let err = (dr - dw).abs();
            if err > 0.01 {
                println!("  Bond 0-4: RDKit={:.4} WebMM={:.4} err={:.4}", dr, dw, err);
            }
        }
        {
            let dr = dist(&rdkit_coords[0], &rdkit_coords[5]);
            let dw = dist(&result.optimized_coords[0], &result.optimized_coords[5]);
            let err = (dr - dw).abs();
            if err > 0.01 {
                println!("  Bond 0-5: RDKit={:.4} WebMM={:.4} err={:.4}", dr, dw, err);
            }
        }
        {
            let dr = dist(&rdkit_coords[1], &rdkit_coords[6]);
            let dw = dist(&result.optimized_coords[1], &result.optimized_coords[6]);
            let err = (dr - dw).abs();
            if err > 0.01 {
                println!("  Bond 1-6: RDKit={:.4} WebMM={:.4} err={:.4}", dr, dw, err);
            }
        }
        {
            let dr = dist(&rdkit_coords[1], &rdkit_coords[7]);
            let dw = dist(&result.optimized_coords[1], &result.optimized_coords[7]);
            let err = (dr - dw).abs();
            if err > 0.01 {
                println!("  Bond 1-7: RDKit={:.4} WebMM={:.4} err={:.4}", dr, dw, err);
            }
        }
        {
            let dr = dist(&rdkit_coords[2], &rdkit_coords[8]);
            let dw = dist(&result.optimized_coords[2], &result.optimized_coords[8]);
            let err = (dr - dw).abs();
            if err > 0.01 {
                println!("  Bond 2-8: RDKit={:.4} WebMM={:.4} err={:.4}", dr, dw, err);
            }
        }
        {
            let dr = dist(&rdkit_coords[2], &rdkit_coords[9]);
            let dw = dist(&result.optimized_coords[2], &result.optimized_coords[9]);
            let err = (dr - dw).abs();
            if err > 0.01 {
                println!("  Bond 2-9: RDKit={:.4} WebMM={:.4} err={:.4}", dr, dw, err);
            }
        }
        {
            let dr = dist(&rdkit_coords[2], &rdkit_coords[10]);
            let dw = dist(&result.optimized_coords[2], &result.optimized_coords[10]);
            let err = (dr - dw).abs();
            if err > 0.01 {
                println!(
                    "  Bond 2-10: RDKit={:.4} WebMM={:.4} err={:.4}",
                    dr, dw, err
                );
            }
        }
    }

    #[test]
    fn opt_benzene() {
        let syms: &[&str] = &["C", "C", "C", "C", "C", "C", "H", "H", "H", "H", "H", "H"];
        let bonds: &[(usize, usize, u8)] = &[
            (0, 1, 15),
            (1, 2, 15),
            (2, 3, 15),
            (3, 4, 15),
            (4, 5, 15),
            (5, 0, 15),
            (0, 6, 10),
            (1, 7, 10),
            (2, 8, 10),
            (3, 9, 10),
            (4, 10, 10),
            (5, 11, 10),
        ];
        let coords: &[[f64; 3]] = &[
            [0.803512, -1.140105, -0.008249],
            [1.388856, 0.125801, -0.028129],
            [0.585344, 1.265905, -0.019873],
            [-0.803512, 1.140105, 0.008259],
            [-1.388857, -0.125800, 0.028123],
            [-0.585344, -1.265905, 0.019872],
            [1.429541, -2.028378, -0.014665],
            [2.470935, 0.223813, -0.050057],
            [1.041395, 2.252191, -0.035361],
            [-1.429541, 2.028378, 0.014703],
            [-2.470936, -0.223813, 0.050023],
            [-1.041394, -2.252192, 0.035354],
        ];
        let mol = make_mol("benzene", syms, bonds, coords);
        let ff = MMFFForceField::new(&mol, MMFFVariant::MMFF94s);
        let rdkit_coords: Vec<[f64; 3]> = coords.to_vec();
        let opts = ConvergenceOptions {
            max_iterations: 2000,
            ..Default::default()
        };
        let result = optimizer::optimize(&ff, &rdkit_coords, &opts);
        let rdkit_e = 16.226967_f64;
        let webmm_e = result.final_energy;
        println!("\n=== BENZENE ===");
        println!(
            "RDKit E: {:.6}  WebMM E: {:.6}  Delta: {:.6}",
            rdkit_e,
            webmm_e,
            webmm_e - rdkit_e
        );
        println!(
            "Converged: {} ({} iters)",
            result.converged, result.iterations
        );
        {
            let dr = dist(&rdkit_coords[0], &rdkit_coords[1]);
            let dw = dist(&result.optimized_coords[0], &result.optimized_coords[1]);
            let err = (dr - dw).abs();
            if err > 0.01 {
                println!("  Bond 0-1: RDKit={:.4} WebMM={:.4} err={:.4}", dr, dw, err);
            }
        }
        {
            let dr = dist(&rdkit_coords[1], &rdkit_coords[2]);
            let dw = dist(&result.optimized_coords[1], &result.optimized_coords[2]);
            let err = (dr - dw).abs();
            if err > 0.01 {
                println!("  Bond 1-2: RDKit={:.4} WebMM={:.4} err={:.4}", dr, dw, err);
            }
        }
        {
            let dr = dist(&rdkit_coords[2], &rdkit_coords[3]);
            let dw = dist(&result.optimized_coords[2], &result.optimized_coords[3]);
            let err = (dr - dw).abs();
            if err > 0.01 {
                println!("  Bond 2-3: RDKit={:.4} WebMM={:.4} err={:.4}", dr, dw, err);
            }
        }
        {
            let dr = dist(&rdkit_coords[3], &rdkit_coords[4]);
            let dw = dist(&result.optimized_coords[3], &result.optimized_coords[4]);
            let err = (dr - dw).abs();
            if err > 0.01 {
                println!("  Bond 3-4: RDKit={:.4} WebMM={:.4} err={:.4}", dr, dw, err);
            }
        }
        {
            let dr = dist(&rdkit_coords[4], &rdkit_coords[5]);
            let dw = dist(&result.optimized_coords[4], &result.optimized_coords[5]);
            let err = (dr - dw).abs();
            if err > 0.01 {
                println!("  Bond 4-5: RDKit={:.4} WebMM={:.4} err={:.4}", dr, dw, err);
            }
        }
        {
            let dr = dist(&rdkit_coords[5], &rdkit_coords[0]);
            let dw = dist(&result.optimized_coords[5], &result.optimized_coords[0]);
            let err = (dr - dw).abs();
            if err > 0.01 {
                println!("  Bond 5-0: RDKit={:.4} WebMM={:.4} err={:.4}", dr, dw, err);
            }
        }
        {
            let dr = dist(&rdkit_coords[0], &rdkit_coords[6]);
            let dw = dist(&result.optimized_coords[0], &result.optimized_coords[6]);
            let err = (dr - dw).abs();
            if err > 0.01 {
                println!("  Bond 0-6: RDKit={:.4} WebMM={:.4} err={:.4}", dr, dw, err);
            }
        }
        {
            let dr = dist(&rdkit_coords[1], &rdkit_coords[7]);
            let dw = dist(&result.optimized_coords[1], &result.optimized_coords[7]);
            let err = (dr - dw).abs();
            if err > 0.01 {
                println!("  Bond 1-7: RDKit={:.4} WebMM={:.4} err={:.4}", dr, dw, err);
            }
        }
        {
            let dr = dist(&rdkit_coords[2], &rdkit_coords[8]);
            let dw = dist(&result.optimized_coords[2], &result.optimized_coords[8]);
            let err = (dr - dw).abs();
            if err > 0.01 {
                println!("  Bond 2-8: RDKit={:.4} WebMM={:.4} err={:.4}", dr, dw, err);
            }
        }
        {
            let dr = dist(&rdkit_coords[3], &rdkit_coords[9]);
            let dw = dist(&result.optimized_coords[3], &result.optimized_coords[9]);
            let err = (dr - dw).abs();
            if err > 0.01 {
                println!("  Bond 3-9: RDKit={:.4} WebMM={:.4} err={:.4}", dr, dw, err);
            }
        }
        {
            let dr = dist(&rdkit_coords[4], &rdkit_coords[10]);
            let dw = dist(&result.optimized_coords[4], &result.optimized_coords[10]);
            let err = (dr - dw).abs();
            if err > 0.01 {
                println!(
                    "  Bond 4-10: RDKit={:.4} WebMM={:.4} err={:.4}",
                    dr, dw, err
                );
            }
        }
        {
            let dr = dist(&rdkit_coords[5], &rdkit_coords[11]);
            let dw = dist(&result.optimized_coords[5], &result.optimized_coords[11]);
            let err = (dr - dw).abs();
            if err > 0.01 {
                println!(
                    "  Bond 5-11: RDKit={:.4} WebMM={:.4} err={:.4}",
                    dr, dw, err
                );
            }
        }
    }

    #[test]
    fn opt_phenol() {
        let syms: &[&str] = &[
            "C", "C", "C", "C", "O", "C", "C", "H", "H", "H", "H", "H", "H",
        ];
        let bonds: &[(usize, usize, u8)] = &[
            (0, 1, 15),
            (1, 2, 15),
            (2, 3, 15),
            (3, 4, 10),
            (3, 5, 15),
            (5, 6, 15),
            (6, 0, 15),
            (0, 7, 10),
            (1, 8, 10),
            (2, 9, 10),
            (4, 10, 10),
            (5, 11, 10),
            (6, 12, 10),
        ];
        let coords: &[[f64; 3]] = &[
            [-1.583016, 0.472765, 0.034105],
            [-0.531677, 1.391140, 0.036561],
            [0.789937, 0.942578, 0.006990],
            [1.046563, -0.425283, -0.024844],
            [2.321647, -0.908022, -0.054391],
            [0.005052, -1.346417, -0.027508],
            [-1.314635, -0.896178, 0.002067],
            [-2.611602, 0.823582, 0.057154],
            [-0.743749, 2.456871, 0.061552],
            [1.597221, 1.667606, 0.009279],
            [2.935484, -0.155754, -0.048557],
            [0.220683, -2.410578, -0.052520],
            [-2.131908, -1.612309, 0.000111],
        ];
        let mol = make_mol("phenol", syms, bonds, coords);
        let ff = MMFFForceField::new(&mol, MMFFVariant::MMFF94s);
        let rdkit_coords: Vec<[f64; 3]> = coords.to_vec();
        let opts = ConvergenceOptions {
            max_iterations: 2000,
            ..Default::default()
        };
        let result = optimizer::optimize(&ff, &rdkit_coords, &opts);
        let rdkit_e = 5.713261_f64;
        let webmm_e = result.final_energy;
        println!("\n=== PHENOL ===");
        println!(
            "RDKit E: {:.6}  WebMM E: {:.6}  Delta: {:.6}",
            rdkit_e,
            webmm_e,
            webmm_e - rdkit_e
        );
        println!(
            "Converged: {} ({} iters)",
            result.converged, result.iterations
        );
        {
            let dr = dist(&rdkit_coords[0], &rdkit_coords[1]);
            let dw = dist(&result.optimized_coords[0], &result.optimized_coords[1]);
            let err = (dr - dw).abs();
            if err > 0.01 {
                println!("  Bond 0-1: RDKit={:.4} WebMM={:.4} err={:.4}", dr, dw, err);
            }
        }
        {
            let dr = dist(&rdkit_coords[1], &rdkit_coords[2]);
            let dw = dist(&result.optimized_coords[1], &result.optimized_coords[2]);
            let err = (dr - dw).abs();
            if err > 0.01 {
                println!("  Bond 1-2: RDKit={:.4} WebMM={:.4} err={:.4}", dr, dw, err);
            }
        }
        {
            let dr = dist(&rdkit_coords[2], &rdkit_coords[3]);
            let dw = dist(&result.optimized_coords[2], &result.optimized_coords[3]);
            let err = (dr - dw).abs();
            if err > 0.01 {
                println!("  Bond 2-3: RDKit={:.4} WebMM={:.4} err={:.4}", dr, dw, err);
            }
        }
        {
            let dr = dist(&rdkit_coords[3], &rdkit_coords[4]);
            let dw = dist(&result.optimized_coords[3], &result.optimized_coords[4]);
            let err = (dr - dw).abs();
            if err > 0.01 {
                println!("  Bond 3-4: RDKit={:.4} WebMM={:.4} err={:.4}", dr, dw, err);
            }
        }
        {
            let dr = dist(&rdkit_coords[3], &rdkit_coords[5]);
            let dw = dist(&result.optimized_coords[3], &result.optimized_coords[5]);
            let err = (dr - dw).abs();
            if err > 0.01 {
                println!("  Bond 3-5: RDKit={:.4} WebMM={:.4} err={:.4}", dr, dw, err);
            }
        }
        {
            let dr = dist(&rdkit_coords[5], &rdkit_coords[6]);
            let dw = dist(&result.optimized_coords[5], &result.optimized_coords[6]);
            let err = (dr - dw).abs();
            if err > 0.01 {
                println!("  Bond 5-6: RDKit={:.4} WebMM={:.4} err={:.4}", dr, dw, err);
            }
        }
        {
            let dr = dist(&rdkit_coords[6], &rdkit_coords[0]);
            let dw = dist(&result.optimized_coords[6], &result.optimized_coords[0]);
            let err = (dr - dw).abs();
            if err > 0.01 {
                println!("  Bond 6-0: RDKit={:.4} WebMM={:.4} err={:.4}", dr, dw, err);
            }
        }
        {
            let dr = dist(&rdkit_coords[0], &rdkit_coords[7]);
            let dw = dist(&result.optimized_coords[0], &result.optimized_coords[7]);
            let err = (dr - dw).abs();
            if err > 0.01 {
                println!("  Bond 0-7: RDKit={:.4} WebMM={:.4} err={:.4}", dr, dw, err);
            }
        }
        {
            let dr = dist(&rdkit_coords[1], &rdkit_coords[8]);
            let dw = dist(&result.optimized_coords[1], &result.optimized_coords[8]);
            let err = (dr - dw).abs();
            if err > 0.01 {
                println!("  Bond 1-8: RDKit={:.4} WebMM={:.4} err={:.4}", dr, dw, err);
            }
        }
        {
            let dr = dist(&rdkit_coords[2], &rdkit_coords[9]);
            let dw = dist(&result.optimized_coords[2], &result.optimized_coords[9]);
            let err = (dr - dw).abs();
            if err > 0.01 {
                println!("  Bond 2-9: RDKit={:.4} WebMM={:.4} err={:.4}", dr, dw, err);
            }
        }
        {
            let dr = dist(&rdkit_coords[4], &rdkit_coords[10]);
            let dw = dist(&result.optimized_coords[4], &result.optimized_coords[10]);
            let err = (dr - dw).abs();
            if err > 0.01 {
                println!(
                    "  Bond 4-10: RDKit={:.4} WebMM={:.4} err={:.4}",
                    dr, dw, err
                );
            }
        }
        {
            let dr = dist(&rdkit_coords[5], &rdkit_coords[11]);
            let dw = dist(&result.optimized_coords[5], &result.optimized_coords[11]);
            let err = (dr - dw).abs();
            if err > 0.01 {
                println!(
                    "  Bond 5-11: RDKit={:.4} WebMM={:.4} err={:.4}",
                    dr, dw, err
                );
            }
        }
        {
            let dr = dist(&rdkit_coords[6], &rdkit_coords[12]);
            let dw = dist(&result.optimized_coords[6], &result.optimized_coords[12]);
            let err = (dr - dw).abs();
            if err > 0.01 {
                println!(
                    "  Bond 6-12: RDKit={:.4} WebMM={:.4} err={:.4}",
                    dr, dw, err
                );
            }
        }
    }

    #[test]
    fn opt_acetaldehyde() {
        let syms: &[&str] = &["C", "C", "O", "H", "H", "H", "H"];
        let bonds: &[(usize, usize, u8)] = &[
            (0, 1, 10),
            (1, 2, 20),
            (0, 3, 10),
            (0, 4, 10),
            (0, 5, 10),
            (1, 6, 10),
        ];
        let coords: &[[f64; 3]] = &[
            [-0.649314, -0.041826, -0.011756],
            [0.844468, 0.060709, 0.016859],
            [1.591518, -0.883737, -0.216234],
            [-1.040888, 0.650177, -0.760908],
            [-0.950404, -1.061159, -0.265659],
            [-1.044533, 0.219461, 0.972603],
            [1.249152, 1.056375, 0.265094],
        ];
        let mol = make_mol("acetaldehyde", syms, bonds, coords);
        let ff = MMFFForceField::new(&mol, MMFFVariant::MMFF94s);
        let rdkit_coords: Vec<[f64; 3]> = coords.to_vec();
        let opts = ConvergenceOptions {
            max_iterations: 2000,
            ..Default::default()
        };
        let result = optimizer::optimize(&ff, &rdkit_coords, &opts);
        let rdkit_e = -0.137887_f64;
        let webmm_e = result.final_energy;
        println!("\n=== ACETALDEHYDE ===");
        println!(
            "RDKit E: {:.6}  WebMM E: {:.6}  Delta: {:.6}",
            rdkit_e,
            webmm_e,
            webmm_e - rdkit_e
        );
        println!(
            "Converged: {} ({} iters)",
            result.converged, result.iterations
        );
        {
            let dr = dist(&rdkit_coords[0], &rdkit_coords[1]);
            let dw = dist(&result.optimized_coords[0], &result.optimized_coords[1]);
            let err = (dr - dw).abs();
            if err > 0.01 {
                println!("  Bond 0-1: RDKit={:.4} WebMM={:.4} err={:.4}", dr, dw, err);
            }
        }
        {
            let dr = dist(&rdkit_coords[1], &rdkit_coords[2]);
            let dw = dist(&result.optimized_coords[1], &result.optimized_coords[2]);
            let err = (dr - dw).abs();
            if err > 0.01 {
                println!("  Bond 1-2: RDKit={:.4} WebMM={:.4} err={:.4}", dr, dw, err);
            }
        }
        {
            let dr = dist(&rdkit_coords[0], &rdkit_coords[3]);
            let dw = dist(&result.optimized_coords[0], &result.optimized_coords[3]);
            let err = (dr - dw).abs();
            if err > 0.01 {
                println!("  Bond 0-3: RDKit={:.4} WebMM={:.4} err={:.4}", dr, dw, err);
            }
        }
        {
            let dr = dist(&rdkit_coords[0], &rdkit_coords[4]);
            let dw = dist(&result.optimized_coords[0], &result.optimized_coords[4]);
            let err = (dr - dw).abs();
            if err > 0.01 {
                println!("  Bond 0-4: RDKit={:.4} WebMM={:.4} err={:.4}", dr, dw, err);
            }
        }
        {
            let dr = dist(&rdkit_coords[0], &rdkit_coords[5]);
            let dw = dist(&result.optimized_coords[0], &result.optimized_coords[5]);
            let err = (dr - dw).abs();
            if err > 0.01 {
                println!("  Bond 0-5: RDKit={:.4} WebMM={:.4} err={:.4}", dr, dw, err);
            }
        }
        {
            let dr = dist(&rdkit_coords[1], &rdkit_coords[6]);
            let dw = dist(&result.optimized_coords[1], &result.optimized_coords[6]);
            let err = (dr - dw).abs();
            if err > 0.01 {
                println!("  Bond 1-6: RDKit={:.4} WebMM={:.4} err={:.4}", dr, dw, err);
            }
        }
    }

    #[test]
    fn opt_glycine() {
        let syms: &[&str] = &["N", "C", "C", "O", "O", "H", "H", "H", "H", "H"];
        let bonds: &[(usize, usize, u8)] = &[
            (0, 1, 10),
            (1, 2, 10),
            (2, 3, 20),
            (2, 4, 10),
            (0, 5, 10),
            (0, 6, 10),
            (1, 7, 10),
            (1, 8, 10),
            (4, 9, 10),
        ];
        let coords: &[[f64; 3]] = &[
            [-1.619902, 0.408956, 0.152410],
            [-0.534743, -0.560379, -0.090505],
            [0.828079, 0.112804, -0.012738],
            [1.048344, 1.292172, 0.216620],
            [1.843398, -0.745836, -0.230233],
            [-1.450537, 0.853652, 1.057399],
            [-1.518263, 1.172983, -0.519436],
            [-0.659660, -0.994458, -1.086560],
            [-0.584386, -1.349233, 0.665405],
            [2.647670, -0.190660, -0.152362],
        ];
        let mol = make_mol("glycine", syms, bonds, coords);
        let ff = MMFFForceField::new(&mol, MMFFVariant::MMFF94s);
        let rdkit_coords: Vec<[f64; 3]> = coords.to_vec();
        let opts = ConvergenceOptions {
            max_iterations: 2000,
            ..Default::default()
        };
        let result = optimizer::optimize(&ff, &rdkit_coords, &opts);
        let rdkit_e = 20.929656_f64;
        let webmm_e = result.final_energy;
        println!("\n=== GLYCINE ===");
        println!(
            "RDKit E: {:.6}  WebMM E: {:.6}  Delta: {:.6}",
            rdkit_e,
            webmm_e,
            webmm_e - rdkit_e
        );
        println!(
            "Converged: {} ({} iters)",
            result.converged, result.iterations
        );
        {
            let dr = dist(&rdkit_coords[0], &rdkit_coords[1]);
            let dw = dist(&result.optimized_coords[0], &result.optimized_coords[1]);
            let err = (dr - dw).abs();
            if err > 0.01 {
                println!("  Bond 0-1: RDKit={:.4} WebMM={:.4} err={:.4}", dr, dw, err);
            }
        }
        {
            let dr = dist(&rdkit_coords[1], &rdkit_coords[2]);
            let dw = dist(&result.optimized_coords[1], &result.optimized_coords[2]);
            let err = (dr - dw).abs();
            if err > 0.01 {
                println!("  Bond 1-2: RDKit={:.4} WebMM={:.4} err={:.4}", dr, dw, err);
            }
        }
        {
            let dr = dist(&rdkit_coords[2], &rdkit_coords[3]);
            let dw = dist(&result.optimized_coords[2], &result.optimized_coords[3]);
            let err = (dr - dw).abs();
            if err > 0.01 {
                println!("  Bond 2-3: RDKit={:.4} WebMM={:.4} err={:.4}", dr, dw, err);
            }
        }
        {
            let dr = dist(&rdkit_coords[2], &rdkit_coords[4]);
            let dw = dist(&result.optimized_coords[2], &result.optimized_coords[4]);
            let err = (dr - dw).abs();
            if err > 0.01 {
                println!("  Bond 2-4: RDKit={:.4} WebMM={:.4} err={:.4}", dr, dw, err);
            }
        }
        {
            let dr = dist(&rdkit_coords[0], &rdkit_coords[5]);
            let dw = dist(&result.optimized_coords[0], &result.optimized_coords[5]);
            let err = (dr - dw).abs();
            if err > 0.01 {
                println!("  Bond 0-5: RDKit={:.4} WebMM={:.4} err={:.4}", dr, dw, err);
            }
        }
        {
            let dr = dist(&rdkit_coords[0], &rdkit_coords[6]);
            let dw = dist(&result.optimized_coords[0], &result.optimized_coords[6]);
            let err = (dr - dw).abs();
            if err > 0.01 {
                println!("  Bond 0-6: RDKit={:.4} WebMM={:.4} err={:.4}", dr, dw, err);
            }
        }
        {
            let dr = dist(&rdkit_coords[1], &rdkit_coords[7]);
            let dw = dist(&result.optimized_coords[1], &result.optimized_coords[7]);
            let err = (dr - dw).abs();
            if err > 0.01 {
                println!("  Bond 1-7: RDKit={:.4} WebMM={:.4} err={:.4}", dr, dw, err);
            }
        }
        {
            let dr = dist(&rdkit_coords[1], &rdkit_coords[8]);
            let dw = dist(&result.optimized_coords[1], &result.optimized_coords[8]);
            let err = (dr - dw).abs();
            if err > 0.01 {
                println!("  Bond 1-8: RDKit={:.4} WebMM={:.4} err={:.4}", dr, dw, err);
            }
        }
        {
            let dr = dist(&rdkit_coords[4], &rdkit_coords[9]);
            let dw = dist(&result.optimized_coords[4], &result.optimized_coords[9]);
            let err = (dr - dw).abs();
            if err > 0.01 {
                println!("  Bond 4-9: RDKit={:.4} WebMM={:.4} err={:.4}", dr, dw, err);
            }
        }
    }

    #[test]
    fn opt_cyclohexane() {
        let syms: &[&str] = &[
            "C", "C", "C", "C", "C", "C", "H", "H", "H", "H", "H", "H", "H", "H", "H", "H", "H",
            "H",
        ];
        let bonds: &[(usize, usize, u8)] = &[
            (0, 1, 10),
            (1, 2, 10),
            (2, 3, 10),
            (3, 4, 10),
            (4, 5, 10),
            (5, 0, 10),
            (0, 6, 10),
            (0, 7, 10),
            (1, 8, 10),
            (1, 9, 10),
            (2, 10, 10),
            (2, 11, 10),
            (3, 12, 10),
            (3, 13, 10),
            (4, 14, 10),
            (4, 15, 10),
            (5, 16, 10),
            (5, 17, 10),
        ];
        let coords: &[[f64; 3]] = &[
            [0.330467, -1.376449, -0.294026],
            [1.456688, -0.345627, -0.191969],
            [0.965662, 1.068917, 0.123219],
            [-0.438260, 1.314553, -0.412492],
            [-1.456688, 0.345627, 0.191967],
            [-0.857868, -1.007020, 0.583302],
            [-0.002500, -1.448099, -1.337242],
            [0.712385, -2.366536, -0.020562],
            [2.169780, -0.650197, 0.583746],
            [2.016114, -0.342985, -1.135386],
            [1.665642, 1.802304, -0.292818],
            [0.966144, 1.219452, 1.210251],
            [-0.747665, 2.346277, -0.210679],
            [-0.430957, 1.199203, -1.503795],
            [-2.270014, 0.197194, -0.528773],
            [-1.915882, 0.795988, 1.080401],
            [-0.532688, -0.970558, 1.630789],
            [-1.630361, -1.782045, 0.524066],
        ];
        let mol = make_mol("cyclohexane", syms, bonds, coords);
        let ff = MMFFForceField::new(&mol, MMFFVariant::MMFF94s);
        let rdkit_coords: Vec<[f64; 3]> = coords.to_vec();
        let opts = ConvergenceOptions {
            max_iterations: 2000,
            ..Default::default()
        };
        let result = optimizer::optimize(&ff, &rdkit_coords, &opts);
        let rdkit_e = 2.368812_f64;
        let webmm_e = result.final_energy;
        println!("\n=== CYCLOHEXANE ===");
        println!(
            "RDKit E: {:.6}  WebMM E: {:.6}  Delta: {:.6}",
            rdkit_e,
            webmm_e,
            webmm_e - rdkit_e
        );
        println!(
            "Converged: {} ({} iters)",
            result.converged, result.iterations
        );
        {
            let dr = dist(&rdkit_coords[0], &rdkit_coords[1]);
            let dw = dist(&result.optimized_coords[0], &result.optimized_coords[1]);
            let err = (dr - dw).abs();
            if err > 0.01 {
                println!("  Bond 0-1: RDKit={:.4} WebMM={:.4} err={:.4}", dr, dw, err);
            }
        }
        {
            let dr = dist(&rdkit_coords[1], &rdkit_coords[2]);
            let dw = dist(&result.optimized_coords[1], &result.optimized_coords[2]);
            let err = (dr - dw).abs();
            if err > 0.01 {
                println!("  Bond 1-2: RDKit={:.4} WebMM={:.4} err={:.4}", dr, dw, err);
            }
        }
        {
            let dr = dist(&rdkit_coords[2], &rdkit_coords[3]);
            let dw = dist(&result.optimized_coords[2], &result.optimized_coords[3]);
            let err = (dr - dw).abs();
            if err > 0.01 {
                println!("  Bond 2-3: RDKit={:.4} WebMM={:.4} err={:.4}", dr, dw, err);
            }
        }
        {
            let dr = dist(&rdkit_coords[3], &rdkit_coords[4]);
            let dw = dist(&result.optimized_coords[3], &result.optimized_coords[4]);
            let err = (dr - dw).abs();
            if err > 0.01 {
                println!("  Bond 3-4: RDKit={:.4} WebMM={:.4} err={:.4}", dr, dw, err);
            }
        }
        {
            let dr = dist(&rdkit_coords[4], &rdkit_coords[5]);
            let dw = dist(&result.optimized_coords[4], &result.optimized_coords[5]);
            let err = (dr - dw).abs();
            if err > 0.01 {
                println!("  Bond 4-5: RDKit={:.4} WebMM={:.4} err={:.4}", dr, dw, err);
            }
        }
        {
            let dr = dist(&rdkit_coords[5], &rdkit_coords[0]);
            let dw = dist(&result.optimized_coords[5], &result.optimized_coords[0]);
            let err = (dr - dw).abs();
            if err > 0.01 {
                println!("  Bond 5-0: RDKit={:.4} WebMM={:.4} err={:.4}", dr, dw, err);
            }
        }
        {
            let dr = dist(&rdkit_coords[0], &rdkit_coords[6]);
            let dw = dist(&result.optimized_coords[0], &result.optimized_coords[6]);
            let err = (dr - dw).abs();
            if err > 0.01 {
                println!("  Bond 0-6: RDKit={:.4} WebMM={:.4} err={:.4}", dr, dw, err);
            }
        }
        {
            let dr = dist(&rdkit_coords[0], &rdkit_coords[7]);
            let dw = dist(&result.optimized_coords[0], &result.optimized_coords[7]);
            let err = (dr - dw).abs();
            if err > 0.01 {
                println!("  Bond 0-7: RDKit={:.4} WebMM={:.4} err={:.4}", dr, dw, err);
            }
        }
        {
            let dr = dist(&rdkit_coords[1], &rdkit_coords[8]);
            let dw = dist(&result.optimized_coords[1], &result.optimized_coords[8]);
            let err = (dr - dw).abs();
            if err > 0.01 {
                println!("  Bond 1-8: RDKit={:.4} WebMM={:.4} err={:.4}", dr, dw, err);
            }
        }
        {
            let dr = dist(&rdkit_coords[1], &rdkit_coords[9]);
            let dw = dist(&result.optimized_coords[1], &result.optimized_coords[9]);
            let err = (dr - dw).abs();
            if err > 0.01 {
                println!("  Bond 1-9: RDKit={:.4} WebMM={:.4} err={:.4}", dr, dw, err);
            }
        }
        {
            let dr = dist(&rdkit_coords[2], &rdkit_coords[10]);
            let dw = dist(&result.optimized_coords[2], &result.optimized_coords[10]);
            let err = (dr - dw).abs();
            if err > 0.01 {
                println!(
                    "  Bond 2-10: RDKit={:.4} WebMM={:.4} err={:.4}",
                    dr, dw, err
                );
            }
        }
        {
            let dr = dist(&rdkit_coords[2], &rdkit_coords[11]);
            let dw = dist(&result.optimized_coords[2], &result.optimized_coords[11]);
            let err = (dr - dw).abs();
            if err > 0.01 {
                println!(
                    "  Bond 2-11: RDKit={:.4} WebMM={:.4} err={:.4}",
                    dr, dw, err
                );
            }
        }
        {
            let dr = dist(&rdkit_coords[3], &rdkit_coords[12]);
            let dw = dist(&result.optimized_coords[3], &result.optimized_coords[12]);
            let err = (dr - dw).abs();
            if err > 0.01 {
                println!(
                    "  Bond 3-12: RDKit={:.4} WebMM={:.4} err={:.4}",
                    dr, dw, err
                );
            }
        }
        {
            let dr = dist(&rdkit_coords[3], &rdkit_coords[13]);
            let dw = dist(&result.optimized_coords[3], &result.optimized_coords[13]);
            let err = (dr - dw).abs();
            if err > 0.01 {
                println!(
                    "  Bond 3-13: RDKit={:.4} WebMM={:.4} err={:.4}",
                    dr, dw, err
                );
            }
        }
        {
            let dr = dist(&rdkit_coords[4], &rdkit_coords[14]);
            let dw = dist(&result.optimized_coords[4], &result.optimized_coords[14]);
            let err = (dr - dw).abs();
            if err > 0.01 {
                println!(
                    "  Bond 4-14: RDKit={:.4} WebMM={:.4} err={:.4}",
                    dr, dw, err
                );
            }
        }
        {
            let dr = dist(&rdkit_coords[4], &rdkit_coords[15]);
            let dw = dist(&result.optimized_coords[4], &result.optimized_coords[15]);
            let err = (dr - dw).abs();
            if err > 0.01 {
                println!(
                    "  Bond 4-15: RDKit={:.4} WebMM={:.4} err={:.4}",
                    dr, dw, err
                );
            }
        }
        {
            let dr = dist(&rdkit_coords[5], &rdkit_coords[16]);
            let dw = dist(&result.optimized_coords[5], &result.optimized_coords[16]);
            let err = (dr - dw).abs();
            if err > 0.01 {
                println!(
                    "  Bond 5-16: RDKit={:.4} WebMM={:.4} err={:.4}",
                    dr, dw, err
                );
            }
        }
        {
            let dr = dist(&rdkit_coords[5], &rdkit_coords[17]);
            let dw = dist(&result.optimized_coords[5], &result.optimized_coords[17]);
            let err = (dr - dw).abs();
            if err > 0.01 {
                println!(
                    "  Bond 5-17: RDKit={:.4} WebMM={:.4} err={:.4}",
                    dr, dw, err
                );
            }
        }
    }

    #[test]
    fn debug_breakdowns() {
        // Cyclohexane
        {
            let syms: &[&str] = &[
                "C", "C", "C", "C", "C", "C", "H", "H", "H", "H", "H", "H", "H", "H", "H", "H",
                "H", "H",
            ];
            let bonds: &[(usize, usize, u8)] = &[
                (0, 1, 10),
                (1, 2, 10),
                (2, 3, 10),
                (3, 4, 10),
                (4, 5, 10),
                (5, 0, 10),
                (0, 6, 10),
                (0, 7, 10),
                (1, 8, 10),
                (1, 9, 10),
                (2, 10, 10),
                (2, 11, 10),
                (3, 12, 10),
                (3, 13, 10),
                (4, 14, 10),
                (4, 15, 10),
                (5, 16, 10),
                (5, 17, 10),
            ];
            let coords: &[[f64; 3]] = &[
                [0.330467, -1.376449, -0.294026],
                [1.456688, -0.345627, -0.191969],
                [0.965662, 1.068917, 0.123219],
                [-0.438260, 1.314553, -0.412492],
                [-1.456688, 0.345627, 0.191967],
                [-0.857868, -1.007020, 0.583302],
                [-0.002500, -1.448099, -1.337242],
                [0.712385, -2.366536, -0.020562],
                [2.169780, -0.650197, 0.583746],
                [2.016114, -0.342985, -1.135386],
                [1.665642, 1.802304, -0.292818],
                [0.966144, 1.219452, 1.210251],
                [-0.747665, 2.346277, -0.210679],
                [-0.430957, 1.199203, -1.503795],
                [-2.270014, 0.197194, -0.528773],
                [-1.915882, 0.795988, 1.080401],
                [-0.532688, -0.970558, 1.630789],
                [-1.630361, -1.782045, 0.524066],
            ];
            let mol = make_mol("cyclohexane", syms, bonds, coords);
            bd("cyclohexane@rdkit_geom", &mol, coords);
            let ff = MMFFForceField::new(&mol, MMFFVariant::MMFF94s);
            let opts = ConvergenceOptions {
                max_iterations: 2000,
                ..Default::default()
            };
            let result = optimizer::optimize(&ff, &coords.to_vec(), &opts);
            bd("cyclohexane@webmm_opt", &mol, &result.optimized_coords);
            println!(
                "RDKit E: 2.368812  WebMM E: {:.6}  Delta: {:.6}",
                result.final_energy,
                result.final_energy - 2.368812
            );
        }
        // Phenol
        {
            let syms: &[&str] = &[
                "C", "C", "C", "C", "O", "C", "C", "H", "H", "H", "H", "H", "H",
            ];
            let bonds: &[(usize, usize, u8)] = &[
                (0, 1, 15),
                (1, 2, 15),
                (2, 3, 15),
                (3, 4, 10),
                (3, 5, 15),
                (5, 6, 15),
                (6, 0, 15),
                (0, 7, 10),
                (1, 8, 10),
                (2, 9, 10),
                (4, 10, 10),
                (5, 11, 10),
                (6, 12, 10),
            ];
            let coords: &[[f64; 3]] = &[
                [-1.583016, 0.472765, 0.034105],
                [-0.531677, 1.391140, 0.036561],
                [0.789937, 0.942578, 0.006990],
                [1.046563, -0.425283, -0.024844],
                [2.321647, -0.908022, -0.054391],
                [0.005052, -1.346417, -0.027508],
                [-1.314635, -0.896178, 0.002067],
                [-2.611602, 0.823582, 0.057154],
                [-0.743749, 2.456871, 0.061552],
                [1.597221, 1.667606, 0.009279],
                [2.935484, -0.155754, -0.048557],
                [0.220683, -2.410578, -0.052520],
                [-2.131908, -1.612309, 0.000111],
            ];
            let mol = make_mol("phenol", syms, bonds, coords);
            bd("phenol@rdkit_geom", &mol, coords);
            let ff = MMFFForceField::new(&mol, MMFFVariant::MMFF94s);
            let opts = ConvergenceOptions {
                max_iterations: 2000,
                ..Default::default()
            };
            let result = optimizer::optimize(&ff, &coords.to_vec(), &opts);
            bd("phenol@webmm_opt", &mol, &result.optimized_coords);
            println!(
                "RDKit E: 5.713261  WebMM E: {:.6}  Delta: {:.6}",
                result.final_energy,
                result.final_energy - 5.713261
            );
        }
        // Benzene
        {
            let syms: &[&str] = &["C", "C", "C", "C", "C", "C", "H", "H", "H", "H", "H", "H"];
            let bonds: &[(usize, usize, u8)] = &[
                (0, 1, 15),
                (1, 2, 15),
                (2, 3, 15),
                (3, 4, 15),
                (4, 5, 15),
                (5, 0, 15),
                (0, 6, 10),
                (1, 7, 10),
                (2, 8, 10),
                (3, 9, 10),
                (4, 10, 10),
                (5, 11, 10),
            ];
            let coords: &[[f64; 3]] = &[
                [0.803512, -1.140105, -0.008249],
                [1.388856, 0.125801, -0.028129],
                [0.585344, 1.265905, -0.019873],
                [-0.803512, 1.140105, 0.008259],
                [-1.388857, -0.125800, 0.028123],
                [-0.585344, -1.265905, 0.019872],
                [1.429541, -2.028378, -0.014665],
                [2.470935, 0.223813, -0.050057],
                [1.041395, 2.252191, -0.035361],
                [-1.429541, 2.028378, 0.014703],
                [-2.470936, -0.223813, 0.050023],
                [-1.041394, -2.252192, 0.035354],
            ];
            let mol = make_mol("benzene", syms, bonds, coords);
            bd("benzene@rdkit_geom", &mol, coords);
            let ff = MMFFForceField::new(&mol, MMFFVariant::MMFF94s);
            let opts = ConvergenceOptions {
                max_iterations: 2000,
                ..Default::default()
            };
            let result = optimizer::optimize(&ff, &coords.to_vec(), &opts);
            bd("benzene@webmm_opt", &mol, &result.optimized_coords);
            println!(
                "RDKit E: 16.227397  WebMM E: {:.6}  Delta: {:.6}",
                result.final_energy,
                result.final_energy - 16.227397
            );
        }
        // Glycine
        {
            let syms: &[&str] = &["N", "C", "C", "O", "O", "H", "H", "H", "H", "H"];
            let bonds: &[(usize, usize, u8)] = &[
                (0, 1, 10),
                (1, 2, 10),
                (2, 3, 20),
                (2, 4, 10),
                (0, 5, 10),
                (0, 6, 10),
                (1, 7, 10),
                (1, 8, 10),
                (4, 9, 10),
            ];
            let coords: &[[f64; 3]] = &[
                [-1.619902, 0.408956, 0.152410],
                [-0.534743, -0.560379, -0.090505],
                [0.828079, 0.112804, -0.012738],
                [1.048344, 1.292172, 0.216620],
                [1.843398, -0.745836, -0.230233],
                [-1.450537, 0.853652, 1.057399],
                [-1.518263, 1.172983, -0.519436],
                [-0.659660, -0.994458, -1.086560],
                [-0.584386, -1.349233, 0.665405],
                [2.647670, -0.190660, -0.152362],
            ];
            let mol = make_mol("glycine", syms, bonds, coords);
            bd("glycine@rdkit_geom", &mol, coords);
            let ff = MMFFForceField::new(&mol, MMFFVariant::MMFF94s);
            let opts = ConvergenceOptions {
                max_iterations: 2000,
                ..Default::default()
            };
            let result = optimizer::optimize(&ff, &coords.to_vec(), &opts);
            bd("glycine@webmm_opt", &mol, &result.optimized_coords);
            println!(
                "RDKit E: 20.929656  WebMM E: {:.6}  Delta: {:.6}",
                result.final_energy,
                result.final_energy - 20.929656
            );
        }
    }
}
