//! Force-source abstraction shared by the optimizer and the MD engine.
//!
/// Anything that can supply a potential energy and its gradient for a set of
/// atomic coordinates. `MMFFForceField` implements this; future additive terms
/// (implicit solvation, a metadynamics bias) will implement it too, so the MD
/// integrator can sum them without knowing their internals.
///
/// Units: energy in kcal/mol, gradient in kcal/(mol·Å) (i.e. dE/dx). The MD
/// integrator treats forces as -gradient.
pub trait ForceField {
    /// Compute the potential energy and write the gradient (dE/dx) into `grad`.
    /// `grad` is zeroed first and must have length == number of atoms.
    /// Returns the energy (kcal/mol).
    fn energy_and_gradient(&self, coords: &[[f64; 3]], grad: &mut [[f64; 3]]) -> f64;
}
