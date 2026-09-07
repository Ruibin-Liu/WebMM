use webmm::optimizer;
use webmm::ConvergenceOptions;
use webmm::forces::ForceField;
use webmm::molecule::parser::parse_sdf;

fn main() {
    // stage1 MMFF, stage2 GFNFF (native pipeline, same as the browser)
    let sdf = std::fs::read_to_string("/tmp/caff24.sdf").unwrap();
    let mol = parse_sdf(&sdf).unwrap();
    let coords: Vec<[f64; 3]> = mol.atoms.iter().map(|a| a.position).collect();
    let ff_mmff = webmm::mmff::MMFFForceField::new(&mol, webmm::MMFFVariant::MMFF94s);
    let r1 = optimizer::optimize(&ff_mmff, &coords, &ConvergenceOptions { max_iterations: 500, ..Default::default() });
    let mut mol2 = mol.clone();
    for (i, p) in r1.optimized_coords.iter().enumerate() { mol2.atoms[i].position = *p; }
    let sdf2 = webmm::molecule::hydrogens::to_molblock(&mol2);
    let mol3 = parse_sdf(&sdf2).unwrap();
    let at: Vec<usize> = mol3.atoms.iter().map(|a| a.atomic_number as usize).collect();
    let charge: f64 = mol3.atoms.iter().map(|a| a.charge).sum();
    let c3: Vec<[f64; 3]> = mol3.atoms.iter().map(|a| a.position).collect();
    let ff_gfn = webmm::gfnff::GfnffForceField::new(&at, &c3, charge);
    let r2 = optimizer::optimize(&ff_gfn, &c3, &ConvergenceOptions { max_iterations: 500, ..Default::default() });
    println!("our GFNFF 500it: E={:.6} Eh converged={} minD={:.4}",
        r2.final_energy / 627.5094740631, r2.converged,
        { let mut m = f64::MAX; for i in 0..r2.optimized_coords.len() { for j in (i+1)..r2.optimized_coords.len() {
            let d = ((r2.optimized_coords[i][0]-r2.optimized_coords[j][0]).powi(2) + (r2.optimized_coords[i][1]-r2.optimized_coords[j][1]).powi(2) + (r2.optimized_coords[i][2]-r2.optimized_coords[j][2]).powi(2)).sqrt();
            m = m.min(d); } } m });
    let mut out = String::from("24\ncaffeine gfnff\n");
    for (a, p) in mol3.atoms.iter().zip(&r2.optimized_coords) {
        out.push_str(&format!("{} {:.6} {:.6} {:.6}\n", a.symbol, p[0], p[1], p[2]));
    }
    std::fs::write("/tmp/gfnff_native.xyz", out).unwrap();
}
