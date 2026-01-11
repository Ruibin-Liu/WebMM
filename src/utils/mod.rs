//! Utility functions

use serde_json::Value;

/// Load MMFF94 parameters from JSON file
pub fn load_mmff_params() -> Value {
    // For now, return empty placeholder
    // TODO: Load from data/mmff94_sample_parameters.json
    serde_json::json!({})
}
