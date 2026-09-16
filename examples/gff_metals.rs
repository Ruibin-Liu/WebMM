fn main() {
    let names = [
        "nicarbonyl",
        "zinc_ammine",
        "cobalt_ammine",
        "ferricyanide",
        "ferrocene",
    ];
    let eh = 627.5094740631;
    let rmap = [
        ("bond", "bond"),
        ("angle", "angle"),
        ("torsion", "torsion"),
        ("repulsion", "rep"),
        ("electrostat", "es"),
        ("dispersion", "disp"),
        ("HB", "hb"),
        ("XB", "xb"),
        ("bonded atm", "batm"),
    ];
    for name in names {
        let sdf = std::fs::read_to_string(format!("/tmp/metals/{name}/{name}.mol")).unwrap();
        let refj: std::collections::HashMap<String, f64> = serde_json::from_str(
            &std::fs::read_to_string(format!("/tmp/metals/{name}/ref.json")).unwrap(),
        )
        .unwrap();
        match webmm::energy_terms_wasm(&sdf, "GFNFF".to_string()) {
            Ok(d) => {
                let d: serde_json::Value = serde_json::from_str(&d).unwrap();
                let t = &d["terms"];
                let e = d["E"].as_f64().unwrap() / eh;
                let xt = refj["total energy"];
                let mut line = format!("{name:14} E {e:+.6} vs {xt:+.6}  d {:+.6}", e - xt);
                for (rk, ok) in rmap {
                    let ours = t[ok].as_f64().unwrap() / eh; // kcal -> Eh
                    let theirs = refj[rk];
                    let dl = ours - theirs;
                    if dl.abs() > 1e-4 {
                        line.push_str(&format!(" |{ok} {dl:+.5}"));
                    }
                }
                println!("{line}");
            }
            Err(e) => println!("{name:14} ERR {e:?}"),
        }
    }
}
