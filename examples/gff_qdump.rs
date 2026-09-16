use webmm::gfnff::Gfnff;
fn main() {
    let name = std::env::args().nth(1).unwrap_or("zinc_ammine".into());
    let sdf = std::fs::read_to_string(format!("/tmp/metals/{name}/{name}.mol")).unwrap();
    let mol = webmm::molecule::parser::parse_sdf(&sdf).unwrap();
    let at: Vec<usize> = mol.atoms.iter().map(|a| a.atomic_number as usize).collect();
    let xyz: Vec<[f64; 3]> = mol.atoms.iter().map(|a| a.position).collect();
    let q: f64 = mol.atoms.iter().map(|a| a.charge).sum();
    let g = Gfnff::new(&at, &xyz, q);
    let ec = g.energy(&xyz);
    println!("E = {:.6} Eh, es = {:.6}", ec.total(), ec.es);
    let t = &g.topo;
    for (i, _z) in at.iter().enumerate() {
        println!(
            "{:3} {:2} qa {:+.5} chi {:+.5} gam {:+.5} alp {:.5} imetal {} mchar {:.3} hyb {}",
            i + 1,
            at[i],
            t.qa[i],
            t.chieeq[i],
            t.gameeq[i],
            t.alpeeq[i],
            t.imetal[i],
            t.mchar[i],
            t.hyb[i]
        );
    }
}
