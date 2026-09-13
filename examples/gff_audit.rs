fn main() {
    let refs: std::collections::HashMap<String, serde_json::Value> =
        serde_json::from_str(&std::fs::read_to_string("/tmp/gff_audit/refs.json").unwrap())
            .unwrap();
    let eh = 627.5094740631;
    let mut names: Vec<&String> = refs.keys().collect();
    names.sort();
    let mut worst = 0.0f64;
    let mut fails = 0;
    for k in names {
        let sdf = std::fs::read_to_string(format!("/tmp/frag_{k}.mol")).unwrap();
        let r = refs[k.as_str()].clone();
        let xt = r["total energy"].as_f64().unwrap();
        match webmm::energy_terms_wasm(&sdf, "GFNFF".to_string()) {
            Ok(d) => {
                let d: serde_json::Value = serde_json::from_str(&d).unwrap();
                let e = d["E"].as_f64().unwrap() / eh;
                let dl = (e - xt).abs();
                worst = worst.max(dl);
                if dl > 1e-4 {
                    fails += 1;
                    println!("{k:26} FAIL {e:.6} vs {xt:.6} |d|={dl:.5}");
                }
            }
            Err(e2) => {
                fails += 1;
                println!("{k:26} ERR {e2:?}");
            }
        }
    }
    println!(
        "GFN-FF parity: worst |d| = {worst:.6} Eh, fails = {fails}/{}",
        refs.len()
    );
}
