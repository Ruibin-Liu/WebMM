//! Electrostatics term for MMFF94

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

    // Gradient: dE/dr = -332.1 * q_i * q_j / (dielectric * r^2)
    let d_e_dr = -332.1 * q_prod / (dielectric * r * r);

    // grad_i = dE/dx_i = dE/dr * dr/dx_i = dE/dr * (-r_vec_x / r)
    let grad_i = [
        -d_e_dr * r_vec[0] / r,
        -d_e_dr * r_vec[1] / r,
        -d_e_dr * r_vec[2] / r,
    ];
    let grad_j = [
        d_e_dr * r_vec[0] / r,
        d_e_dr * r_vec[1] / r,
        d_e_dr * r_vec[2] / r,
    ];

    (energy, grad_i, grad_j)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_electrostatic_energy() {
        let coords = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let charges = vec![-1.0, 1.0];

        let (energy, grad_i, grad_j) =
            electrostatic_energy_and_gradient(&coords, &charges, 0, 1, 1.0);

        assert!(
            energy < 0.0,
            "Opposite charges should have attractive energy"
        );
        // Gradient convention: dE/dx. For opposite charges, moving atom i toward j decreases energy.
        // So grad_i (negative direction in gradient space) = negative x (descent is toward j).
        assert!(
            grad_i[0] < 0.0,
            "grad_i should be negative (descent direction pulls i toward j)"
        );
        assert!(
            grad_j[0] > 0.0,
            "grad_j should be positive (descent direction pulls j toward i)"
        );
    }
}
