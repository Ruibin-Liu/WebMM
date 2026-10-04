// Native golden test: SA score + bit-exact Morgan vs Python-generated
// references (tests/fixtures/sascore/*.json). Run: cargo test --test sascore_golden
use std::collections::HashMap;

#[path = "../src/sascore.rs"]
mod sascore;

fn build(mb: &str) -> sascore::SaGraph {
    // minimal V2000 parse (aromatic type-4, M CHG)
    let raw: Vec<&str> = mb.lines().collect();
    let ci = raw
        .iter()
        .position(|l| l.contains("V2000"))
        .expect("no V2000");
    let start = ci.saturating_sub(3);
    let lines: Vec<&str> = raw[start..].to_vec();
    let counts = lines[3];
    let na: usize = counts[0..3].trim().parse().unwrap();
    let nb: usize = counts[3..6].trim().parse().unwrap();
    let mut zs = Vec::new();
    let mut chg = vec![0i32; na];
    let mut dm = vec![0i32; na];
    let mut par = vec![0u8; na];
    let mut bonds = Vec::new();
    let sym_z = |s: &str| -> u32 {
        match s {
            "H" => 1,
            "B" => 5,
            "C" => 6,
            "N" => 7,
            "O" => 8,
            "F" => 9,
            "Si" => 14,
            "P" => 15,
            "S" => 16,
            "Cl" => 17,
            "Se" => 34,
            "Br" => 35,
            "I" => 53,
            _ => 0,
        }
    };
    for i in 0..na {
        let l = lines[4 + i];
        let parts: Vec<&str> = l.split_whitespace().collect();
        zs.push(sym_z(parts[3]));
        let b = l.as_bytes();
        if b.len() >= 36 {
            let t: String = String::from_utf8_lossy(&b[34..36]).trim().to_string();
            dm[i] = t.parse().unwrap_or(0);
        }
        if b.len() >= 42 {
            let t: String = String::from_utf8_lossy(&b[39..42]).trim().to_string();
            if let Ok(v) = t.parse::<u8>() {
                par[i] = v;
            }
        }
    }
    for i in 0..nb {
        let l = lines[4 + na + i];
        let a: usize = l[0..3].trim().parse().unwrap();
        let bidx: usize = l[3..6].trim().parse().unwrap();
        let o: u8 = l[6..9].trim().parse().unwrap();
        bonds.push((a - 1, bidx - 1, o));
    }
    // M CHG
    for l in lines.iter().skip(4 + na + nb) {
        let t = l.trim();
        if t.starts_with("M  CHG") {
            let parts: Vec<&str> = t.split_whitespace().collect();
            let cnt: usize = parts[2].parse().unwrap();
            let mut k = 3;
            for _ in 0..cnt {
                let ai: usize = parts[k].parse().unwrap();
                let c: i32 = parts[k + 1].parse().unwrap();
                chg[ai - 1] = c;
                k += 2;
            }
        }
    }
    sascore::build_graph(&zs, &chg, &dm, &par, &bonds)
}

#[test]
fn sa_golden_bits_and_scores() {
    let tbl = std::fs::read("app/fpscores.bin").expect("app/fpscores.bin");
    let n = sascore::sa_load_table(&tbl).expect("table load");
    assert_eq!(n, 705292);
    for fname in [
        "tests/fixtures/sascore/golden.json",
        "tests/fixtures/sascore/probes.json",
    ] {
        let data = std::fs::read_to_string(fname).unwrap();
        let items: Vec<serde_json::Value> = serde_json::from_str(&data).unwrap();
        let mut bit_ok = 0;
        let mut bit_fail: Vec<&str> = Vec::new();
        let mut sa_err_max = 0.0f64;
        let mut sa_worst = "";
        for it in &items {
            let g = build(it["mb"].as_str().unwrap());
            let counts = sascore::morgan_sparse_counts(&g, 2);
            let mut want: HashMap<u32, u32> = HashMap::new();
            for (k, v) in it["fp"].as_object().unwrap() {
                want.insert(k.parse().unwrap(), v.as_u64().unwrap() as u32);
            }
            if counts == want {
                bit_ok += 1;
            } else {
                bit_fail.push(it["smiles"].as_str().unwrap());
            }
            let sa = sascore::sa_score_from_graph(&g).unwrap();
            let ref_sa = it["sa"].as_f64().unwrap();
            let err = (sa - ref_sa).abs();
            if err > sa_err_max {
                sa_err_max = err;
                sa_worst = it["smiles"].as_str().unwrap();
            }
        }
        println!("{}: bits {}/{} ok", fname, bit_ok, items.len());
        assert!(
            bit_fail.is_empty(),
            "bit mismatches: {:?}",
            &bit_fail[..bit_fail.len().min(5)]
        );
        println!(
            "{}: max |SA err| = {:.4} (worst: {})",
            fname, sa_err_max, sa_worst
        );
        // tolerance covers the documented stereo/bridgehead penalty
        // approximation (nonzero only for the steroid in this corpus)
        assert!(
            sa_err_max < 1.5,
            "SA error too large: {} {}",
            sa_err_max,
            sa_worst
        );
    }
}
