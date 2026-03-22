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
        assert!(
            grad_i[0] < 0.0,
            "grad_i should be negative (descent direction pulls i toward j)"
        );
        assert!(
            grad_j[0] > 0.0,
            "grad_j should be positive (descent direction pulls j toward i)"
        );
    }

    #[test]
    fn test_electrostatic_like_charges() {
        let coords = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let charges = vec![1.0, 1.0];

        let (energy, grad_i, grad_j) =
            electrostatic_energy_and_gradient(&coords, &charges, 0, 1, 1.0);

        assert!(
            energy > 0.0,
            "Like charges should have repulsive energy, got {}",
            energy
        );
        assert!(
            grad_i[0] > 0.0,
            "grad_i should push atom i away from j for like charges"
        );
        assert!(
            grad_j[0] < 0.0,
            "grad_j should push atom j away from i for like charges"
        );
    }

    #[test]
    fn test_electrostatic_energy_magnitude() {
        let coords = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let charges = vec![-1.0, 1.0];

        let (energy, _, _) = electrostatic_energy_and_gradient(&coords, &charges, 0, 1, 1.0);

        let expected = 332.1 * (-1.0) * 1.0 / (1.0 * 1.0);
        assert!(
            (energy - expected).abs() < 0.01,
            "Energy = {}, expected {}",
            energy,
            expected
        );
    }

    #[test]
    fn test_electrostatic_dielectric() {
        let coords = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let charges = vec![-1.0, 1.0];

        let (e_d1, _, _) = electrostatic_energy_and_gradient(&coords, &charges, 0, 1, 1.0);
        let (e_d4, _, _) = electrostatic_energy_and_gradient(&coords, &charges, 0, 1, 4.0);

        assert!(
            (e_d1 / e_d4 - 4.0).abs() < 0.01,
            "Doubling dielectric should halve energy: e_d1={}, e_d4={}",
            e_d1,
            e_d4
        );
    }

    #[test]
    fn test_electrostatic_3d_geometry() {
        let coords = vec![[0.0, 0.0, 0.0], [1.0, 2.0, 3.0]];
        let charges = vec![-1.0, 1.0];

        let (energy, grad_i, grad_j) =
            electrostatic_energy_and_gradient(&coords, &charges, 0, 1, 1.0);

        assert!(energy < 0.0);
        assert!(energy.is_finite());

        // grad_i and grad_j should be equal and opposite
        for dim in 0..3 {
            assert!(
                (grad_i[dim] + grad_j[dim]).abs() < 1e-10,
                "grad_i[{}] + grad_j[{}] should be zero (Newton's 3rd law), got {} + {} = {}",
                dim,
                dim,
                grad_i[dim],
                grad_j[dim],
                grad_i[dim] + grad_j[dim]
            );
        }
    }

    #[test]
    fn test_electrostatic_zero_charge() {
        let coords = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let charges = vec![0.0, 1.0];

        let (energy, grad_i, grad_j) =
            electrostatic_energy_and_gradient(&coords, &charges, 0, 1, 1.0);

        assert!(
            energy.abs() < 1e-10,
            "Energy should be zero when one charge is zero"
        );
        for dim in 0..3 {
            assert!(
                grad_i[dim].abs() < 1e-10,
                "grad_i[{}] should be zero when charge is zero",
                dim
            );
            assert!(
                grad_j[dim].abs() < 1e-10,
                "grad_j[{}] should be zero when other charge is zero",
                dim
            );
        }
    }

    #[test]
    fn test_electrostatic_coincident_atoms() {
        let coords = vec![[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]];
        let charges = vec![-1.0, 1.0];

        let (energy, grad_i, grad_j) =
            electrostatic_energy_and_gradient(&coords, &charges, 0, 1, 1.0);

        assert!(
            energy.abs() < 1e-10,
            "Energy should be zero for coincident atoms"
        );
        for dim in 0..3 {
            assert!(grad_i[dim].abs() < 1e-10);
            assert!(grad_j[dim].abs() < 1e-10);
        }
    }

    #[test]
    fn test_electrostatic_gradient_numerical() {
        let coords = vec![[0.0, 0.0, 0.0], [1.0, 2.0, 3.0]];
        let charges = vec![-0.5, 0.3];

        let (_, grad_i, grad_j) = electrostatic_energy_and_gradient(&coords, &charges, 0, 1, 1.5);

        let eps = 1e-7;
        for (atom_idx, grad) in [(0, grad_i), (1, grad_j)] {
            for dim in 0..3 {
                let mut coords_p = coords.clone();
                coords_p[atom_idx][dim] += eps;
                let (e_plus, _, _) =
                    electrostatic_energy_and_gradient(&coords_p, &charges, 0, 1, 1.5);
                let (e0, _, _) = electrostatic_energy_and_gradient(&coords, &charges, 0, 1, 1.5);
                let num_grad = (e_plus - e0) / eps;
                assert!(
                    (grad[dim] - num_grad).abs() < 1e-3,
                    "Analytical grad[{}] = {} but numerical = {} for atom {}",
                    dim,
                    grad[dim],
                    num_grad,
                    atom_idx
                );
            }
        }
    }
}
