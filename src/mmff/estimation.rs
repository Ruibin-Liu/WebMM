use super::atom_types::get_atom_type_props;
use super::MMFFAtomType;
use crate::molecule::BondType;

#[rustfmt::skip]
/// MMFFCovRadPauEle (RDKit defaultMMFFCovRadPauEle): (Z, r0 [A], Pauling chi).
const COVRAD_PAUELE: &[(u8, f64, f64); 18] = &[
    (1, 0.33, 2.20),
    (3, 1.34, 0.97),
    (6, 0.77, 2.50),
    (7, 0.73, 3.07),
    (8, 0.72, 3.50),
    (9, 0.74, 4.12),
    (11, 1.54, 1.01),
    (12, 1.30, 1.23),
    (14, 1.15, 1.74),
    (15, 1.09, 2.06),
    (16, 1.03, 2.44),
    (17, 1.01, 2.83),
    (19, 1.96, 0.91),
    (20, 1.74, 1.04),
    (29, 1.38, 1.75),
    (30, 1.31, 1.66),
    (35, 1.15, 2.74),
    (53, 1.33, 2.21),
];

#[rustfmt::skip]
/// MMFFBndk (RDKit defaultMMFFBndk): (Z_min, Z_max, r0, kb) — reference
/// single-bond length/force used to scale kb in the empirical rule.
const BNDK: &[(u8, u8, f64, f64); 58] = &[
    (1, 6, 1.084, 5.15),
    (1, 7, 1.001, 7.35),
    (1, 8, 0.947, 9.10),
    (1, 9, 0.92, 10.6),
    (1, 14, 1.48, 2.3),
    (1, 15, 1.415, 2.95),
    (1, 16, 1.326, 4.30),
    (1, 17, 1.28, 4.3),
    (1, 35, 1.41, 4.2),
    (1, 53, 1.60, 2.7),
    (6, 6, 1.512, 3.80),
    (6, 7, 1.439, 4.55),
    (6, 8, 1.393, 5.40),
    (6, 9, 1.353, 6.20),
    (6, 14, 1.86, 2.6),
    (6, 15, 1.84, 2.7),
    (6, 16, 1.812, 2.85),
    (6, 17, 1.781, 2.75),
    (6, 35, 1.94, 2.6),
    (6, 53, 2.16, 1.4),
    (7, 7, 1.283, 6.00),
    (7, 8, 1.333, 5.90),
    (7, 9, 1.36, 5.9),
    (7, 14, 1.74, 3.7),
    (7, 15, 1.65, 4.8),
    (7, 16, 1.674, 3.75),
    (7, 17, 1.75, 3.5),
    (7, 35, 1.90, 2.9),
    (7, 53, 2.10, 1.6),
    (8, 8, 1.48, 3.6),
    (8, 9, 1.42, 4.6),
    (8, 14, 1.63, 5.2),
    (8, 15, 1.66, 4.7),
    (8, 16, 1.470, 9.90),
    (8, 17, 1.70, 4.1),
    (8, 35, 1.85, 3.4),
    (8, 53, 2.05, 1.6),
    (9, 14, 1.57, 6.4),
    (9, 15, 1.54, 7.1),
    (9, 16, 1.55, 6.9),
    (14, 14, 2.32, 1.3),
    (14, 15, 2.25, 1.5),
    (14, 16, 2.15, 2.0),
    (14, 17, 2.02, 3.1),
    (14, 35, 2.19, 2.1),
    (14, 53, 2.44, 1.5),
    (15, 15, 2.21, 1.7),
    (15, 16, 2.10, 2.4),
    (15, 17, 2.03, 3.0),
    (15, 35, 2.21, 2.0),
    (15, 53, 2.47, 1.4),
    (16, 16, 2.052, 2.50),
    (16, 17, 2.04, 2.9),
    (16, 35, 2.24, 1.9),
    (16, 53, 2.40, 1.7),
    (17, 17, 1.99, 3.5),
    (35, 35, 2.28, 2.4),
    (53, 53, 2.67, 1.6),
];

#[rustfmt::skip]
/// Herschbach-Laurie Badger-rule params (RDKit defaultMMFFHerschbachLaurie):
/// (row_i, row_j, a_ij, d_ij) keyed by RDKit getPeriodicTableRowHL.
const HERSCHBACH_LAURIE: &[(u8, u8, f64, f64); 25] = &[
    (0, 0, 1.26, 0.025),
    (0, 1, 1.66, 0.30),
    (0, 2, 1.84, 0.38),
    (0, 3, 1.98, 0.49),
    (0, 4, 2.03, 0.51),
    (0, 5, 2.03, 0.25),
    (0, 30, 1.85, 0.15),
    (0, 40, 1.84, 0.61),
    (0, 50, 1.78, 0.97),
    (1, 1, 1.91, 0.68),
    (1, 2, 2.28, 0.74),
    (1, 3, 2.35, 0.85),
    (1, 4, 2.33, 0.68),
    (1, 5, 2.50, 0.97),
    (1, 30, 2.08, 1.14),
    (1, 40, 2.34, 1.17),
    (2, 2, 2.41, 1.18),
    (2, 3, 2.52, 1.02),
    (2, 4, 2.61, 1.28),
    (2, 5, 2.60, 0.84),
    (3, 3, 2.58, 1.41),
    (3, 4, 2.66, 0.86),
    (3, 5, 2.75, 1.14),
    (4, 4, 2.85, 1.62),
    (4, 5, 2.76, 1.25),
];

/// RDKit `getPeriodicTableRowHL`: row for the Herschbach-Laurie rule.
/// H = 0, He = 1, then 2..=5; transition metals get their row × 10.
fn hl_row(z: i32) -> u8 {
    let mut row: u8 = if z == 2 {
        1
    } else if (3..=10).contains(&z) {
        2
    } else if (11..=18).contains(&z) {
        3
    } else if (19..=36).contains(&z) {
        4
    } else if (37..=54).contains(&z) {
        5
    } else {
        0 // hydrogen
    };
    if (21..=30).contains(&z) || (39..=48).contains(&z) {
        row *= 10;
    }
    row
}

/// RDKit `getMMFFBondStretchEmpiricalRuleParams` (AtomTyper.cpp) — the exact
/// empirical rule MMFF94 uses for bond pairs missing from the release table
/// (e.g. S(17)-O(6) in sulfite esters, S(17)-S(15) in thiosulfinates):
///   r0 = r0_i + r0_j - c * |chi_i - chi_j|^1.4        (eq. 18, MMFF.V p.625)
///   kb = kb_bndk * (r0_bndk / r0)^6                   (eq. 19)
/// or, when no Bndk row exists for the element pair, the Herschbach-Laurie
/// version of Badger's rule:
///   kb = 10^(-(r0 - a_ij) / d_ij)
/// Note: bond order does NOT enter (RDKit compiles the BO branch out).
pub fn estimate_bond_params_rdkit(z1: i32, z2: i32) -> Option<(f64, f64)> {
    let cr1 = COVRAD_PAUELE.iter().find(|t| t.0 as i32 == z1)?;
    let cr2 = COVRAD_PAUELE.iter().find(|t| t.0 as i32 == z2)?;
    let c = if z1 == 1 || z2 == 1 { 0.050 } else { 0.085 };
    let r0 = cr1.1 + cr2.1 - c * (cr1.2 - cr2.2).abs().powf(1.4);
    let kb = if let Some(b) = BNDK
        .iter()
        .find(|t| t.0 as i32 == z1.min(z2) && t.1 as i32 == z1.max(z2))
    {
        let coeff = b.2 / r0;
        b.3 * coeff.powi(6)
    } else {
        let (r1, r2) = (hl_row(z1), hl_row(z2));
        let hl = HERSCHBACH_LAURIE
            .iter()
            .find(|t| t.0 == r1.min(r2) && t.1 == r1.max(r2))?;
        10f64.powf(-(r0 - hl.2) / hl.3)
    };
    Some((kb, r0))
}

pub fn estimate_bond_params(
    type1: MMFFAtomType,
    type2: MMFFAtomType,
    bond_type: BondType,
) -> Option<(f64, f64)> {
    let props1 = get_atom_type_props(type1)?;
    let props2 = get_atom_type_props(type2)?;

    let bc1 = props1.bond_class as f64;
    let bc2 = props2.bond_class as f64;

    let k_bond = (2.0 * bc1 * bc2) / (bc1 + bc2);

    let mut r0 = props1.crd + props2.crd - 0.01 * (props1.crd - props2.crd).powi(2);

    match bond_type {
        BondType::Double => r0 *= 0.94,
        BondType::Triple => r0 *= 0.90,
        BondType::Aromatic => r0 *= 0.97,
        BondType::Single => {}
    }

    Some((k_bond, r0))
}

#[cfg(test)]
mod tests {

    #[test]
    fn rdkit_bond_empirical_rule_matches_reference() {
        // RDKit getMMFFBondStretchEmpiricalRuleParams replicas (Bndk-scaled);
        // values cross-checked against the formula with RDKit's tables.
        // S-S (thiosulfinate core): r0 = 2.060, kb ~ 2.4423
        let (kb, r0) = estimate_bond_params_rdkit(16, 16).unwrap();
        assert!((r0 - 2.0600).abs() < 1e-3, "S-S r0 = {r0}");
        assert!((kb - 2.4423).abs() < 1e-3, "S-S kb = {kb}");
        // C-C: r0 = 1.540, kb ~ 3.4038 (both from Bndk)
        let (kb, r0) = estimate_bond_params_rdkit(6, 6).unwrap();
        assert!((r0 - 1.5400).abs() < 1e-3, "C-C r0 = {r0}");
        assert!((kb - 3.4038).abs() < 1e-3, "C-C kb = {kb}");
        // O-H goes through the c = 0.05 branch
        let (_, r0) = estimate_bond_params_rdkit(1, 8).unwrap();
        assert!((r0 - 0.978).abs() < 5e-3, "O-H r0 = {r0}");
    }

    use super::*;

    #[test]
    fn test_c3_c3_single_bond() {
        let (k_bond, r0) =
            estimate_bond_params(MMFFAtomType::C_3, MMFFAtomType::C_3, BondType::Single).unwrap();
        assert!((r0 - 1.54).abs() < 0.05, "r0 = {r0}, expected ~1.54");
        // Force constant must be in mdyn/Å (real MMFF bond kb ~1.5–17). The
        // harmonic-mean-of-bond-class estimate gives ~2.0 for C_3-C_3; the
        // true table value is 4.258. It must NOT include the kcal conversion
        // factor (bond_energy applies c1=143.9325 separately).
        assert!(
            (1.0..=5.0).contains(&k_bond),
            "k_bond = {k_bond}, expected ~2 mdyn/Å (table value 4.258)"
        );
    }

    #[test]
    fn test_c3_c3_double_bond_shorter() {
        let (_k_s, r0_s) =
            estimate_bond_params(MMFFAtomType::C_3, MMFFAtomType::C_3, BondType::Single).unwrap();
        let (_k_d, r0_d) =
            estimate_bond_params(MMFFAtomType::C_3, MMFFAtomType::C_3, BondType::Double).unwrap();
        assert!(
            r0_d < r0_s,
            "Double bond ({:.3}) should be shorter than single ({:.3})",
            r0_d,
            r0_s
        );
    }

    #[test]
    fn test_c3_c3_triple_bond_shortest() {
        let (_k_s, r0_s) =
            estimate_bond_params(MMFFAtomType::C_3, MMFFAtomType::C_3, BondType::Single).unwrap();
        let (_k_d, r0_d) =
            estimate_bond_params(MMFFAtomType::C_3, MMFFAtomType::C_3, BondType::Double).unwrap();
        let (_k_t, r0_t) =
            estimate_bond_params(MMFFAtomType::C_3, MMFFAtomType::C_3, BondType::Triple).unwrap();
        assert!(r0_t < r0_d, "Triple ({:.3}) < double ({:.3})", r0_t, r0_d);
        assert!(r0_d < r0_s, "Double ({:.3}) < single ({:.3})", r0_d, r0_s);
    }

    #[test]
    fn test_aromatic_bond_intermediate() {
        let (_k_s, r0_s) =
            estimate_bond_params(MMFFAtomType::C_3, MMFFAtomType::C_3, BondType::Single).unwrap();
        let (_k_a, r0_a) =
            estimate_bond_params(MMFFAtomType::C_3, MMFFAtomType::C_3, BondType::Aromatic).unwrap();
        assert!(
            r0_a < r0_s,
            "Aromatic ({:.3}) should be shorter than single ({:.3})",
            r0_a,
            r0_s
        );
    }

    #[test]
    fn test_heteroatom_bond_estimation() {
        // C-N bond should have finite parameters
        let result = estimate_bond_params(MMFFAtomType::C_3, MMFFAtomType::N_3, BondType::Single);
        assert!(result.is_some(), "C_3-N_3 estimation should succeed");
        let (_k, r0) = result.unwrap();
        assert!(
            r0 > 1.0 && r0 < 2.0,
            "C-N bond length should be reasonable, got {r0}"
        );
    }

    #[test]
    fn test_symmetric_bond_types() {
        let (k1, r0_1) =
            estimate_bond_params(MMFFAtomType::C_3, MMFFAtomType::N_3, BondType::Single).unwrap();
        let (k2, r0_2) =
            estimate_bond_params(MMFFAtomType::N_3, MMFFAtomType::C_3, BondType::Single).unwrap();
        assert!(
            (k1 - k2).abs() < 1e-10,
            "Bond force constant should be symmetric"
        );
        assert!(
            (r0_1 - r0_2).abs() < 1e-10,
            "Bond length should be symmetric"
        );
    }
}
