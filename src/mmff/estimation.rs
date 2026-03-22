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
    fn test_tetrahedral_angle() {
        let (k_theta, theta0) =
            estimate_angle_params(MMFFAtomType::C_3, MMFFAtomType::C_3, MMFFAtomType::C_3, 4)
                .unwrap();
        assert!((theta0 - 109.47).abs() < 0.01, "theta0 = {theta0}");
        assert!(k_theta > 0.0);
    }

    #[test]
    fn test_ion_types_return_none() {
        let result = estimate_bond_params(MMFFAtomType::Na, MMFFAtomType::C_3, BondType::Single);
        assert!(result.is_none());

        let result =
            estimate_angle_params(MMFFAtomType::Na, MMFFAtomType::C_3, MMFFAtomType::C_3, 4);
        assert!(result.is_none());
    }
}
