//! Fixed-pose Gaussian shape Tanimoto between two SDF inputs (validation
//! tooling for the RDKit-grid correlation script — no alignment involved).
//!
//! Usage: cargo run --release --example shape_fixed_pose -- a.sdf b.sdf

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() != 3 {
        eprintln!("usage: shape_fixed_pose <a.sdf> <b.sdf>");
        std::process::exit(2);
    }
    let a =
        webmm::molecule::parser::parse_sdf(&std::fs::read_to_string(&args[1]).unwrap()).unwrap();
    let b =
        webmm::molecule::parser::parse_sdf(&std::fs::read_to_string(&args[2]).unwrap()).unwrap();
    let sa = webmm::shape::shape_atoms(&a);
    let sb = webmm::shape::shape_atoms(&b);
    println!("{:.9}", webmm::shape::shape_tanimoto(&sa, &sb));
}
