//! Electrostatics term for MMFF94
//!
//! RDKit uses buffered distance: corr_dist = r + 0.05
//! E = 332.0716 * q_i * q_j / corr_dist
//! 1-4 interactions are scaled by 0.75

/// Coulomb conversion constant ( kcal·Å/(mol·e²) )
const COULOMB_CONST: f64 = 332.0716;

/// Electrostatic energy between two charged atoms
/// RDKit-compatible buffered distance: corr_dist = r + 0.05
/// 1-4 interactions are scaled by 0.75
pub fn electrostatic_energy_and_gradient(
    coords: &[[f64; 3]],
    charges: &[f64],
    i: usize,
    j: usize,
    _dielectric: f64,
    is_14: bool,
) -> (f64, [f64; 3], [f64; 3]) {
    let r_vec = [
        coords[j][0] - coords[i][0],
        coords[j][1] - coords[i][1],
        coords[j][2] - coords[i][2],
    ];
    let r_sq = r_vec[0] * r_vec[0] + r_vec[1] * r_vec[1] + r_vec[2] * r_vec[2];
    let r = r_sq.sqrt();

    if r < 1e-10 {
        return (0.0, [0.0; 3], [0.0; 3]);
    }

    let q_prod = charges[i] * charges[j];

    // RDKit buffered distance: corr_dist = r + 0.05
    let buffer = 0.05;
    let corr_dist = r + buffer;

    // 1-4 scaling
    let scale = if is_14 { 0.75 } else { 1.0 };

    // E = 332.0716 * q_i * q_j / (r + 0.05) * scale
    let energy = scale * COULOMB_CONST * q_prod / corr_dist;

    // dE/dr = -332.0716 * q_i * q_j / (r + 0.05)² * scale
    let d_e_dr = -scale * COULOMB_CONST * q_prod / (corr_dist * corr_dist);

    // grad_i = dE/dx_i = dE/dr * (-r_vec / r)
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
    fn test_electrostatic_opposite_charges() {
        let coords = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let charges = vec![-1.0, 1.0];

        let (energy, grad_i, grad_j) =
            electrostatic_energy_and_gradient(&coords, &charges, 0, 1, 1.0, false);

        assert!(
            energy < 0.0,
            "Opposite charges should have attractive energy, got {}",
            energy
        );
        // With buffered distance, opposite charges attract so grad_i pulls i toward j
        assert!(
            grad_i[0] < 0.0,
            "grad_i should pull atom i toward j, got {}",
            grad_i[0]
        );
        assert!(
            grad_j[0] > 0.0,
            "grad_j should pull atom j toward i, got {}",
            grad_j[0]
        );
    }

    #[test]
    fn test_electrostatic_like_charges() {
        let coords = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let charges = vec![1.0, 1.0];

        let (energy, grad_i, grad_j) =
            electrostatic_energy_and_gradient(&coords, &charges, 0, 1, 1.0, false);

        assert!(
            energy > 0.0,
            "Like charges should have repulsive energy, got {}",
            energy
        );
        assert!(
            grad_i[0] > 0.0,
            "grad_i should push atom i away from j, got {}",
            grad_i[0]
        );
        assert!(
            grad_j[0] < 0.0,
            "grad_j should push atom j away from i, got {}",
            grad_j[0]
        );
    }

    #[test]
    fn test_electrostatic_energy_magnitude() {
        let coords = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let charges = vec![-1.0, 1.0];

        let (energy, _, _) = electrostatic_energy_and_gradient(&coords, &charges, 0, 1, 1.0, false);

        // E = 332.0716 * (-1) * 1 / (1.0 + 0.05) = -316.2587
        let expected = COULOMB_CONST * (-1.0) * 1.0 / (1.0 + 0.05);
        assert!(
            (energy - expected).abs() < 0.01,
            "Energy = {}, expected {}",
            energy,
            expected
        );
    }

    #[test]
    fn test_electrostatic_14_scaling() {
        let coords = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let charges = vec![-1.0, 1.0];

        let (e_full, _, _) = electrostatic_energy_and_gradient(&coords, &charges, 0, 1, 1.0, false);
        let (e_14, _, _) = electrostatic_energy_and_gradient(&coords, &charges, 0, 1, 1.0, true);

        assert!(
            (e_14 - 0.75 * e_full).abs() < 1e-10,
            "1-4 electrostatic should be 0.75 * full: got {}, expected {}",
            e_14,
            0.75 * e_full
        );
    }

    #[test]
    fn test_electrostatic_zero_charge() {
        let coords = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let charges = vec![0.0, 1.0];

        let (energy, grad_i, grad_j) =
            electrostatic_energy_and_gradient(&coords, &charges, 0, 1, 1.0, false);

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
            electrostatic_energy_and_gradient(&coords, &charges, 0, 1, 1.0, false);

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
    fn test_electrostatic_3d_geometry() {
        let coords = vec![[0.0, 0.0, 0.0], [1.0, 2.0, 3.0]];
        let charges = vec![-1.0, 1.0];

        let (energy, grad_i, grad_j) =
            electrostatic_energy_and_gradient(&coords, &charges, 0, 1, 1.0, false);

        assert!(energy < 0.0);
        assert!(energy.is_finite());

        // Newton's 3rd law: grad_i + grad_j = 0
        for dim in 0..3 {
            assert!(
                (grad_i[dim] + grad_j[dim]).abs() < 1e-10,
                "grad_i[{}] + grad_j[{}] = {} (should be zero)",
                dim,
                dim,
                grad_i[dim] + grad_j[dim]
            );
        }
    }

    #[test]
    fn test_electrostatic_gradient_numerical() {
        let coords = vec![[0.0, 0.0, 0.0], [1.0, 2.0, 3.0]];
        let charges = vec![-0.5, 0.3];

        let (_, grad_i, grad_j) =
            electrostatic_energy_and_gradient(&coords, &charges, 0, 1, 1.5, false);

        let eps = 1e-7;
        for (atom_idx, grad) in [(0, grad_i), (1, grad_j)] {
            for dim in 0..3 {
                let mut coords_p = coords.clone();
                coords_p[atom_idx][dim] += eps;
                let (e_plus, _, _) =
                    electrostatic_energy_and_gradient(&coords_p, &charges, 0, 1, 1.5, false);
                let (e0, _, _) =
                    electrostatic_energy_and_gradient(&coords, &charges, 0, 1, 1.5, false);
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
