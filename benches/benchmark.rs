//! Benchmarks for MolGopt
//!
//! Run with: `cargo bench`

use std::time::Instant;

fn main() {
    println!("MolGopt Benchmark Suite");
    println!("======================\n");

    // Load test molecules
    let water_sdf = r#"Water
     RDKit          3D

  3  2  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0
    0.9580    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0
   -0.2390    0.9270    0.0000 H   0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  1  3  1  0  0  0  0
M  END"#;

    let benzene_sdf = r#"Benzene
     RDKit          3D

  6  6  0  0  0  0  0  0  0  0999 V2000
    0.0000    1.4010    0.0000 C   0  0  0  0  0  0  0  0  0
    1.2115    0.7035    0.0000 C   0  0  0  0  0  0  0  0  0
   -0.6060    1.0493    0.0000 C   0  0  0  0  0  0  0  0  0
   -1.2115    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0
   -0.6060   -1.0493    0.0000 C   0  0  0  0  0  0  0  0  0
    1.2115   -0.7035    0.0000 C   0  0  0  0  0  0  0  0  0
  1  2  4  0  0  0  0
  2  3  4  0  0  0  0
  3  4  4  0  0  0  0
  4  5  4  0  0  0  0
  5  6  4  0  0  0  0
  6  1  4  0  0  0  0
M  END"#;

    let ethanol_sdf = r#"Ethanol
     RDKit          3D

  9  8  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0
    1.5260    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0
    2.1450    1.0190    0.0000 O   0  0  0  0  0  0  0  0  0
   -0.5430    1.0200    0.0000 H   0  0  0  0  0  0  0  0  0
   -0.5430   -1.0200    0.0000 H   0  0  0  0  0  0  0  0  0
    1.8860   -1.0200    0.0000 H   0  0  0  0  0  0  0  0  0
    2.0190   -0.5430    0.0000 H   0  0  0  0  0  0  0  0  0
    1.5270    1.5640    0.9170 H   0  0  0  0  0  0  0  0  0
    3.0290    1.3430    0.3590 H   0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  1  0  0  0  0
  1  4  1  0  0  0  0
  1  5  1  0  0  0  0
  2  6  1  0  0  0  0
  2  7  1  0  0  0  0
  3  8  1  0  0  0  0
  3  9  1  0  0  0  0
M  END"#;

    // Run benchmarks
    benchmark_parsing(water_sdf, benzene_sdf, ethanol_sdf);
    benchmark_etkdg(water_sdf, benzene_sdf, ethanol_sdf);
    benchmark_energy(water_sdf, benzene_sdf, ethanol_sdf);
    benchmark_optimization(water_sdf, benzene_sdf, ethanol_sdf);
}

fn benchmark_parsing(water: &str, benzene: &str, ethanol: &str) {
    println!("SDF Parsing Benchmarks");
    println!("----------------------");

    let iterations = 1000;

    let start = Instant::now();
    for _ in 0..iterations {
        let _ = molgopt::molecule::parser::parse_sdf(water).unwrap();
    }
    let water_time = start.elapsed().as_micros() as f64 / iterations as f64;
    println!("Water (3 atoms):       {:>8.2} µs/op", water_time);

    let start = Instant::now();
    for _ in 0..iterations {
        let _ = molgopt::molecule::parser::parse_sdf(benzene).unwrap();
    }
    let benzene_time = start.elapsed().as_micros() as f64 / iterations as f64;
    println!("Benzene (6 atoms):     {:>8.2} µs/op", benzene_time);

    let start = Instant::now();
    for _ in 0..iterations {
        let _ = molgopt::molecule::parser::parse_sdf(ethanol).unwrap();
    }
    let ethanol_time = start.elapsed().as_micros() as f64 / iterations as f64;
    println!("Ethanol (9 atoms):     {:>8.2} µs/op", ethanol_time);

    println!();
}

fn benchmark_etkdg(water: &str, benzene: &str, ethanol: &str) {
    println!("ETKDG Embedding Benchmarks");
    println!("---------------------------");

    let water_mol = molgopt::molecule::parser::parse_sdf(water).unwrap();
    let benzene_mol = molgopt::molecule::parser::parse_sdf(benzene).unwrap();
    let ethanol_mol = molgopt::molecule::parser::parse_sdf(ethanol).unwrap();

    let iterations = 100;

    let start = Instant::now();
    for _ in 0..iterations {
        let _ = molgopt::etkdg::generate_initial_coords(&water_mol);
    }
    let water_time = start.elapsed().as_millis() as f64 / iterations as f64;
    println!("Water (3 atoms):       {:>8.2} ms/op", water_time);

    let start = Instant::now();
    for _ in 0..iterations {
        let _ = molgopt::etkdg::generate_initial_coords(&benzene_mol);
    }
    let benzene_time = start.elapsed().as_millis() as f64 / iterations as f64;
    println!("Benzene (6 atoms):     {:>8.2} ms/op", benzene_time);

    let start = Instant::now();
    for _ in 0..iterations {
        let _ = molgopt::etkdg::generate_initial_coords(&ethanol_mol);
    }
    let ethanol_time = start.elapsed().as_millis() as f64 / iterations as f64;
    println!("Ethanol (9 atoms):     {:>8.2} ms/op", ethanol_time);

    println!();
}

fn benchmark_energy(water: &str, benzene: &str, ethanol: &str) {
    println!("MMFF Energy Calculation Benchmarks");
    println!("-----------------------------------");

    let water_mol = molgopt::molecule::parser::parse_sdf(water).unwrap();
    let benzene_mol = molgopt::molecule::parser::parse_sdf(benzene).unwrap();
    let ethanol_mol = molgopt::molecule::parser::parse_sdf(ethanol).unwrap();

    let water_coords = molgopt::etkdg::generate_initial_coords(&water_mol);
    let benzene_coords = molgopt::etkdg::generate_initial_coords(&benzene_mol);
    let ethanol_coords = molgopt::etkdg::generate_initial_coords(&ethanol_mol);

    let water_ff = molgopt::mmff::MMFFForceField::new(&water_mol, molgopt::MMFFVariant::MMFF94s);
    let benzene_ff =
        molgopt::mmff::MMFFForceField::new(&benzene_mol, molgopt::MMFFVariant::MMFF94s);
    let ethanol_ff =
        molgopt::mmff::MMFFForceField::new(&ethanol_mol, molgopt::MMFFVariant::MMFF94s);

    let iterations = 1000;

    let start = Instant::now();
    for _ in 0..iterations {
        let _ = water_ff.calculate_energy(&water_coords);
    }
    let water_time = start.elapsed().as_micros() as f64 / iterations as f64;
    println!("Water (3 atoms):       {:>8.2} µs/op", water_time);

    let start = Instant::now();
    for _ in 0..iterations {
        let _ = benzene_ff.calculate_energy(&benzene_coords);
    }
    let benzene_time = start.elapsed().as_micros() as f64 / iterations as f64;
    println!("Benzene (6 atoms):     {:>8.2} µs/op", benzene_time);

    let start = Instant::now();
    for _ in 0..iterations {
        let _ = ethanol_ff.calculate_energy(&ethanol_coords);
    }
    let ethanol_time = start.elapsed().as_micros() as f64 / iterations as f64;
    println!("Ethanol (9 atoms):     {:>8.2} µs/op", ethanol_time);

    println!();
}

fn benchmark_optimization(water: &str, benzene: &str, ethanol: &str) {
    println!("L-BFGS Optimization Benchmarks");
    println!("-------------------------------");

    let water_mol = molgopt::molecule::parser::parse_sdf(water).unwrap();
    let benzene_mol = molgopt::molecule::parser::parse_sdf(benzene).unwrap();
    let ethanol_mol = molgopt::molecule::parser::parse_sdf(ethanol).unwrap();

    let water_coords = molgopt::etkdg::generate_initial_coords(&water_mol);
    let benzene_coords = molgopt::etkdg::generate_initial_coords(&benzene_mol);
    let ethanol_coords = molgopt::etkdg::generate_initial_coords(&ethanol_mol);

    let water_ff = molgopt::mmff::MMFFForceField::new(&water_mol, molgopt::MMFFVariant::MMFF94s);
    let benzene_ff =
        molgopt::mmff::MMFFForceField::new(&benzene_mol, molgopt::MMFFVariant::MMFF94s);
    let ethanol_ff =
        molgopt::mmff::MMFFForceField::new(&ethanol_mol, molgopt::MMFFVariant::MMFF94s);

    let conv = molgopt::ConvergenceOptions::default();

    let iterations = 10;

    let start = Instant::now();
    for _ in 0..iterations {
        let _ = molgopt::optimizer::optimize(&water_ff, &water_coords, &conv);
    }
    let water_time = start.elapsed().as_millis() as f64 / iterations as f64;
    println!("Water (3 atoms):       {:>8.2} ms/op", water_time);

    let start = Instant::now();
    for _ in 0..iterations {
        let _ = molgopt::optimizer::optimize(&benzene_ff, &benzene_coords, &conv);
    }
    let benzene_time = start.elapsed().as_millis() as f64 / iterations as f64;
    println!("Benzene (6 atoms):     {:>8.2} ms/op", benzene_time);

    let start = Instant::now();
    for _ in 0..iterations {
        let _ = molgopt::optimizer::optimize(&ethanol_ff, &ethanol_coords, &conv);
    }
    let ethanol_time = start.elapsed().as_millis() as f64 / iterations as f64;
    println!("Ethanol (9 atoms):     {:>8.2} ms/op", ethanol_time);

    println!();
}
