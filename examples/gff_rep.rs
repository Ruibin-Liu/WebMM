use webmm::gfnff::Gfnff;
fn main() {
    let name = std::env::args().nth(1).unwrap_or("ferrocene".into());
    let sdf = std::fs::read_to_string(format!("/tmp/metals/{name}/{name}.mol")).unwrap();
    let mol = webmm::molecule::parser::parse_sdf(&sdf).unwrap();
    let at: Vec<usize> = mol.atoms.iter().map(|a| a.atomic_number as usize).collect();
    let xyz: Vec<[f64; 3]> = mol.atoms.iter().map(|a| a.position).collect();
    let g = Gfnff::new(&at, &xyz, mol.atoms.iter().map(|a| a.charge).sum());
    let xyzb: Vec<[f64; 3]> = xyz
        .iter()
        .map(|r| [r[0] / 0.52917726, r[1] / 0.52917726, r[2] / 0.52917726])
        .collect();
    let mut sum = 0.0f64;
    let mut groups: std::collections::BTreeMap<String, (f64, usize)> = Default::default();
    for i in 0..at.len() {
        for j in 0..i {
            if g.topo.nb[i].contains(&j) {
                continue;
            }
            let r = ((xyzb[i][0] - xyzb[j][0]).powi(2)
                + (xyzb[i][1] - xyzb[j][1]).powi(2)
                + (xyzb[i][2] - xyzb[j][2]).powi(2))
            .sqrt();
            if r > 20.0 {
                continue;
            }
            let (zi, zj) = (g.at[i], g.at[j]);
            let fni = 1.0 + g.p.nrepscal / (1.0 + (g.topo.nb[i].len() as f64).powi(2));
            let fnj = 1.0 + g.p.nrepscal / (1.0 + (g.topo.nb[j].len() as f64).powi(2));
            let di = g.p.repan[zi - 1] * (1.0 + g.topo.qa[i] * g.p.qrepscal) * fni;
            let dj = g.p.repan[zj - 1] * (1.0 + g.topo.qa[j] * g.p.qrepscal) * fnj;
            let mut ff = 1.0;
            if zi == 1 && zj == 1 {
                ff = g.p.hhfac;
                let bp = g
                    .topo
                    .bpair
                    .get(i)
                    .and_then(|r| r.get(j))
                    .copied()
                    .unwrap_or(5);
                if bp == 2 {
                    ff *= g.p.hh13rep;
                }
                if bp == 3 {
                    ff *= g.p.hh14rep;
                }
            }
            if (zi == 1 && g.p.metal[zj - 1] > 0) || (zj == 1 && g.p.metal[zi - 1] > 0) {
                ff = 0.85;
            }
            if (zi == 1 && zj == 6) || (zj == 1 && zi == 6) {
                ff = 0.91;
            }
            if (zi == 1 && zj == 8) || (zj == 1 && zi == 8) {
                ff = 1.04;
            }
            let alpha = (di * dj).sqrt() * ff;
            let t16 = r.powf(1.5);
            let e = (-alpha * t16).exp() * g.p.repz[zi - 1] * g.p.repz[zj - 1] * g.p.repscaln / r;
            sum += e;
            let key = format!("{}-{}", zi.min(zj), zi.max(zj));
            let ent = groups.entry(key).or_insert((0.0, 0));
            ent.0 += e;
            ent.1 += 1;
        }
    }
    println!("rep sum = {sum:.6} Eh");
    for (k, (e, c)) in &groups {
        println!("  pair {k:6}: {e:+.6} Eh over {c} pairs");
    }
}
