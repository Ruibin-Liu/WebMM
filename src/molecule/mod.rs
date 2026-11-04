//! Molecule module shim: implementation migrated to the `avidity-engine`
//! crate (github.com/Ruibin-Liu/Avidity). Re-exported here so existing
//! `crate::molecule::…` paths keep working unchanged.
pub use avidity_engine::molecule::*;
pub use avidity_engine::molecule::{graph, hydrogens, parser};
