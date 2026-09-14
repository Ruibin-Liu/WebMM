fn main() {
    let refs: std::collections::HashMap<String, Option<f64>> =
        serde_json::from_str(&std::fs::read_to_string("/tmp/frag_ref_2026.json").unwrap()).unwrap();
    let mut worst = 0.0f64;
    let mut fails = 0;
    let mut n = 0;
    for (name, r) in &refs {
        let Some(ref r) = r else { continue };
        let sdf = std::fs::read_to_string(format!("/tmp/frag_{name}.mol")).unwrap();
        match webmm::energy_terms_wasm(&sdf, "MMFF94s".to_string()) {
            Ok(d) => {
                let d: serde_json::Value = serde_json::from_str(&d).unwrap();
                let e = d["E"].as_f64().unwrap();
                let dl = (e - r).abs();
                n += 1;
                if dl > 0.01 {
                    fails += 1;
                    println!("{name:26} FAIL {e:.5} vs {r:.5} |d|={dl:.5}");
                }
                worst = worst.max(dl);
            }
            Err(e) => {
                fails += 1;
                println!("{name:26} ERR {e:?}");
            }
        }
    }
    println!("MMFF parity vs RDKit 2026.03.6: {n} molecules, worst |d| = {worst:.5} kcal, fails = {fails}");
}
