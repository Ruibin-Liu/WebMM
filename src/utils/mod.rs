//! Utility functions

use serde_json::Value;

pub fn load_mmff_params() -> Value {
    let json_str = include_str!("../../data/mmff94_sample_parameters.json");
    serde_json::from_str(json_str).expect("Failed to parse embedded MMFF parameters")
}

pub fn get_bond_params_from_json(
    type1_name: &str,
    type2_name: &str,
    bond_type_name: &str,
) -> Option<crate::mmff::bond::BondParams> {
    let params = load_mmff_params();
    for (_key, val) in params.get("bond_parameters")?.as_object()? {
        if val.get("type1")?.as_str()? == type1_name
            && val.get("type2")?.as_str()? == type2_name
            && val.get("bond_type")?.as_str()? == bond_type_name
        {
            return Some(crate::mmff::bond::BondParams {
                k_bond: val.get("kb")?.as_f64()?,
                r0: val.get("r0")?.as_f64()?,
            });
        }
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_load_mmff_params() {
        let params = load_mmff_params();
        assert!(params.is_object());
        assert!(params.get("bond_parameters").is_some());
        assert!(params.get("angle_parameters").is_some());
        assert!(params.get("vdw_parameters").is_some());
    }

    #[test]
    fn test_get_bond_params_from_json() {
        let result = get_bond_params_from_json("C_3", "C_3", "Single");
        assert!(result.is_some());
        let params = result.unwrap();
        assert!((params.k_bond - 4.7).abs() < 1e-6);
        assert!((params.r0 - 1.526).abs() < 1e-6);
    }
}
