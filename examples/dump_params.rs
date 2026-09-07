//! Dump webmm's resolved per-bond / per-angle / per-type-vdW params for an SDF dir.
//! Run: cargo run --release --example dump_params -- <SDF_DIR> > params.json
//! Used by the parameter audit (scripts/audit_params.py) to compare against RDKit.
use std::io::Write;
use webmm::mmff::{
    get_angle_params_with_bond_info, get_bond_params, get_mmff_bond_type, get_vdw_params,
    MMFFForceField,
};
use webmm::molecule::parser::parse_sdf;
use webmm::molecule::BondType;
use webmm::MMFFVariant;

fn bo(b: BondType) -> f64 {
    match b {
        BondType::Single => 1.0,
        BondType::Double => 2.0,
        BondType::Triple => 3.0,
        BondType::Aromatic => 1.5,
    }
}

fn main() {
    let dir = std::env::args()
        .nth(1)
        .unwrap_or_else(|| "scripts/val_set".to_string());
    let mut files: Vec<String> = std::fs::read_dir(&dir)
        .expect("dir")
        .flatten()
        .filter(|e| e.path().extension().and_then(|s| s.to_str()) == Some("sdf"))
        .map(|e| e.path().file_stem().unwrap().to_string_lossy().into_owned())
        .collect();
    files.sort();

    let mut out = String::from("{");
    let mut first = true;
    for name in &files {
        if !first {
            out.push(',');
        }
        first = false;
        let path = format!("{dir}/{name}.sdf");
        let mol = match std::fs::read_to_string(&path)
            .ok()
            .and_then(|t| parse_sdf(&t).ok())
        {
            Some(m) => m,
            None => {
                out.push_str(&format!("\"{name}\":{{\"error\":\"parse\"}}"));
                continue;
            }
        };
        let ff = MMFFForceField::new(&mol, MMFFVariant::MMFF94s);
        let tids = &ff.type_ids;
        let at = &ff.atom_types;
        // bonds
        let mut bonds = String::from("[");
        for (bi, b) in mol.bonds.iter().enumerate() {
            if bi > 0 {
                bonds.push(',');
            }
            let p = get_bond_params(at[b.atom1], at[b.atom2], b.bond_type);
            let (kb, r0) = p.map(|q| (q.k_bond, q.r0)).unwrap_or((f64::NAN, f64::NAN));
            bonds.push_str(&format!(
                "{{\"i\":{},\"j\":{},\"ti\":{},\"tj\":{},\"bo\":{},\"kb\":{:?},\"r0\":{:?}}}",
                b.atom1,
                b.atom2,
                tids[b.atom1],
                tids[b.atom2],
                bo(b.bond_type),
                kb,
                r0
            ));
        }
        bonds.push(']');
        // angles (reproduce new()'s resolution via the public API)
        let mut angles = String::from("[");
        for (ai, a) in ff.angles.iter().enumerate() {
            if ai > 0 {
                angles.push(',');
            }
            let (i, j, k) = (a.atom1, a.atom2, a.atom3);
            let bij = (i.min(j), i.max(j));
            let bkj = (k.min(j), k.max(j));
            let bt_ij = ff
                .bond_map
                .get(&bij)
                .map(|b| get_mmff_bond_type(b.bond_type, tids[i], tids[j]))
                .unwrap_or(0);
            let bt_jk = ff
                .bond_map
                .get(&bkj)
                .map(|b| get_mmff_bond_type(b.bond_type, tids[j], tids[k]))
                .unwrap_or(0);
            let ring = ff.angle_ring_size(i, j, k);
            let r0_ij = ff
                .bond_map
                .get(&bij)
                .and_then(|b| get_bond_params(at[i], at[j], b.bond_type))
                .map(|p| p.r0)
                .unwrap_or(1.5);
            let r0_jk = ff
                .bond_map
                .get(&bkj)
                .and_then(|b| get_bond_params(at[k], at[j], b.bond_type))
                .map(|p| p.r0)
                .unwrap_or(1.5);
            let p = get_angle_params_with_bond_info(
                at[i], at[j], at[k], bt_ij, bt_jk, ring, r0_ij, r0_jk,
            );
            let (ka, t0) = p
                .map(|q| (q.k_theta, q.theta0))
                .unwrap_or((f64::NAN, f64::NAN));
            angles.push_str(&format!(
                "{{\"i\":{},\"j\":{},\"k\":{},\"ti\":{},\"tj\":{},\"tk\":{},\"ka\":{:?},\"t0\":{:?}}}",
                i, j, k, tids[i], tids[j], tids[k], ka, t0
            ));
        }
        angles.push(']');
        // vdW per unique type
        let mut vdw = String::from("{");
        let mut vf = true;
        let mut seen = Vec::new();
        for idx in 0..mol.atoms.len() {
            let t = tids[idx];
            if seen.contains(&t) {
                continue;
            }
            seen.push(t);
            let p = get_vdw_params(at[idx]);
            if !vf {
                vdw.push(',');
            }
            vf = false;
            vdw.push_str(&format!(
                "\"{}\":[{:?},{:?},{:?},{:?},{}]",
                t, p.r_star, p.alpha_i, p.n_i, p.g_i, p.da
            ));
        }
        vdw.push('}');
        out.push_str(&format!(
            "\"{name}\":{{\"bonds\":{bonds},\"angles\":{angles},\"vdw\":{vdw}}}"
        ));
    }
    out.push('}');
    let _ = write!(std::io::stdout(), "{out}");
}
