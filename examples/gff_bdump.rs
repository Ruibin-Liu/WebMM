use webmm::gfnff::Gfnff;
fn main() {
    let name = std::env::args().nth(1).unwrap_or("zinc_ammine".into());
    let sdf = std::fs::read_to_string(format!("/tmp/metals/{name}/{name}.mol")).unwrap();
    let mol = webmm::molecule::parser::parse_sdf(&sdf).unwrap();
    let at: Vec<usize> = mol.atoms.iter().map(|a| a.atomic_number as usize).collect();
    let xyz: Vec<[f64; 3]> = mol.atoms.iter().map(|a| a.position).collect();
    let g = Gfnff::new(&at, &xyz, mol_total_charge(&mol));
    let _ = g.energy(&xyz);
    let xyzb: Vec<[f64; 3]> = xyz
        .iter()
        .map(|r| [r[0] / 0.52917726, r[1] / 0.52917726, r[2] / 0.52917726])
        .collect();
    let dist = |i: usize, j: usize| -> f64 {
        let dx = xyzb[i][0] - xyzb[j][0];
        let dy = xyzb[i][1] - xyzb[j][1];
        let dz = xyzb[i][2] - xyzb[j][2];
        (dx * dx + dy * dy + dz * dz).sqrt()
    };
    let xyzb: Vec<[f64; 3]> = xyz
        .iter()
        .map(|r| [r[0] / 0.52917726, r[1] / 0.52917726, r[2] / 0.52917726])
        .collect();
    let cn = webmm::gfnff::erf_cn(&g.p, &at, &xyzb);
    let mut sum = 0.0f64;
    println!("ours: bond  R(A) rab0(A) dr(A)  kb     E(Eh)");
    for b in &g.bonds {
        let r = dist(b.i, b.j);
        let rab0 = g.p.gfnffrab(g.at[b.i], g.at[b.j], cn[b.i], cn[b.j], b.r0);
        let dr = r - rab0;
        let e = b.kb * (-b.alp * dr * dr).exp();
        sum += e;
        println!(
            "{:3} {:3} {:2}-{:2} {:.3} {:.3} {:+.4} {:+.5} {:+.6}",
            b.i + 1,
            b.j + 1,
            g.at[b.i],
            g.at[b.j],
            r * 0.52917726,
            rab0 * 0.52917726,
            dr * 0.52917726,
            b.kb,
            e
        );
    }
    println!("sum bond = {sum:.6} Eh");
}
fn mol_total_charge(mol: &webmm::molecule::Molecule) -> f64 {
    mol.atoms.iter().map(|a| a.charge).sum()
}
