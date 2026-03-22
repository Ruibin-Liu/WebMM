use super::atom_types::get_atom_type_props;
use super::MMFFAtomType;
use crate::molecule::BondType;

pub fn estimate_bond_params(
    type1: MMFFAtomType,
    type2: MMFFAtomType,
    bond_type: BondType,
) -> Option<(f64, f64)> {
    let props1 = get_atom_type_props(type1)?;
    let props2 = get_atom_type_props(type2)?;

    let bc1 = props1.bond_class as f64;
    let bc2 = props2.bond_class as f64;

    let k_bond = 71.9662 * (2.0 * bc1 * bc2) / (bc1 + bc2);

    let mut r0 = props1.crd + props2.crd - 0.01 * (props1.crd - props2.crd).powi(2);

    match bond_type {
        BondType::Double => r0 *= 0.94,
        BondType::Triple => r0 *= 0.90,
        BondType::Aromatic => r0 *= 0.97,
        BondType::Single => {}
    }

    Some((k_bond, r0))
}

pub fn estimate_angle_params(
    type1: MMFFAtomType,
    type2: MMFFAtomType,
    type3: MMFFAtomType,
    num_neighbors_central: usize,
) -> Option<(f64, f64)> {
    let props1 = get_atom_type_props(type1)?;
    let props2 = get_atom_type_props(type2)?;
    let props3 = get_atom_type_props(type3)?;

    let ac1 = props1.angle_class as f64;
    let ac2 = props2.angle_class as f64;
    let ac3 = props3.angle_class as f64;

    let k_theta = 0.001 * (ac1 * ac2 * ac3).sqrt().max(0.001);

    let theta0 = match num_neighbors_central {
        2 => 180.0,
        3 => 120.0,
        _ => 109.47,
    };

    Some((k_theta, theta0))
}

pub fn default_torsion_params() -> (f64, f64, f64) {
    (0.0, 0.0, 0.0)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c3_c3_single_bond() {
        let (k_bond, r0) =
            estimate_bond_params(MMFFAtomType::C_3, MMFFAtomType::C_3, BondType::Single).unwrap();
        assert!((r0 - 1.54).abs() < 0.05, "r0 = {r0}, expected ~1.54");
        assert!(k_bond > 100.0, "k_bond = {k_bond}, expected > 100");
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
    fn test_tetrahedral_angle() {
        let (k_theta, theta0) =
            estimate_angle_params(MMFFAtomType::C_3, MMFFAtomType::C_3, MMFFAtomType::C_3, 4)
                .unwrap();
        assert!((theta0 - 109.47).abs() < 0.01, "theta0 = {theta0}");
        assert!(k_theta > 0.0);
    }

    #[test]
    fn test_linear_angle() {
        let (_k, theta0) =
            estimate_angle_params(MMFFAtomType::C_2, MMFFAtomType::C_2, MMFFAtomType::C_2, 2)
                .unwrap();
        assert!(
            (theta0 - 180.0).abs() < 0.01,
            "Linear angle theta0 = {theta0}, expected 180"
        );
    }

    #[test]
    fn test_trigonal_planar_angle() {
        let (_k, theta0) =
            estimate_angle_params(MMFFAtomType::C_2, MMFFAtomType::C_2, MMFFAtomType::C_2, 3)
                .unwrap();
        assert!(
            (theta0 - 120.0).abs() < 0.01,
            "Trigonal angle theta0 = {theta0}, expected 120"
        );
    }

    #[test]
    fn test_default_torsion_params() {
        let (v1, v2, v3) = default_torsion_params();
        assert!((v1, v2, v3) == (0.0, 0.0, 0.0));
    }

    #[test]
    fn test_ion_types_return_none() {
        let result = estimate_bond_params(MMFFAtomType::Na, MMFFAtomType::C_3, BondType::Single);
        assert!(result.is_none());

        let result =
            estimate_angle_params(MMFFAtomType::Na, MMFFAtomType::C_3, MMFFAtomType::C_3, 4);
        assert!(result.is_none());
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
