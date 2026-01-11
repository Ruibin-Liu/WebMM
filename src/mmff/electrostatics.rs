//! Electrostatics term for MMFF94

use crate::molecule::Molecule;

/// Electrostatic energy between two charged atoms
pub fn electrostatic_energy_and_gradient(
    coords: &[[f64; 3]],
    charges: &[f64],
    i: usize,
    j: usize,
    dielectric: f64,
) -> (f64, [f64; 3], [f64; 3]) {
    let r_vec = [
        coords[j][0] - coords[i][0],
        coords[j][1] - coords[i][1],
        coords[j][2] - coords[i][2],
    ];
    let r = (r_vec[0].powi(2) + r_vec[1].powi(2) + r_vec[2].powi(2)).sqrt();

    if r < 1e-10 {
        return (0.0, [0.0; 3], [0.0; 3]);
    }

    let q_prod = charges[i] * charges[j];

    // Energy: 332.1 * q_i * q_j / (dielectric * r)  (kcal/mol)
    let energy = 332.1 * q_prod / (dielectric * r);

    // Gradient: dE/dx = -332.1 * q_i * q_j / (dielectric * r^2) * (dx/dr)
    let dE_dr = -332.1 * q_prod / (dielectric * r * r);
    let factor = dE_dr / r;

    // Force on atom i is negative gradient (attractive if charges opposite)
    // Force on atom j is positive gradient
    let grad_i = [factor * r_vec[0], factor * r_vec[1], factor * r_vec[2]];
    let grad_j = [-factor * r_vec[0], -factor * r_vec[1], -factor * r_vec[2]];

    (energy, grad_i, grad_j)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_electrostatic_energy() {
        let coords = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let charges = vec![-1.0, 1.0]; // Opposite charges

        let (energy, grad_i, grad_j) =
            electrostatic_energy_and_gradient(&coords, &charges, 0, 1, 1.0);

        assert!(energy < 0.0); // Attractive energy
                               // Forces should point toward each other (attraction)
                               // grad_i is force on atom i (at index 0, q=-1), should point to j (positive x)
                               // grad_j is force on atom j (at index 1, q=+1), should point to i (negative x)
        assert!(
            grad_i[0] > 0.0,
            "Force on negative charge should point toward positive charge"
        );
        assert!(
            grad_j[0] < 0.0,
            "Force on positive charge should point toward negative charge"
        );
    }
}
