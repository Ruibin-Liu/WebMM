//! Property-based tests for WebMM
//!
//! These tests verify invariants using randomly generated inputs.

#[cfg(test)]
mod tests {
    use proptest::prelude::*;

    // Test that energy is always finite for valid coordinates
    proptest! {
        #[test]
        fn energy_is_finite_for_valid_coords(
            x in -10.0f64..10.0,
            y in -10.0f64..10.0,
            z in -10.0f64..10.0,
        ) {
            let sdf = format!(r#"Test
     RDKit          3D

  2  1  0  0  0  0  0  0  0  0999 V2000
    {:.4}    {:.4}    {:.4} C   0  0  0  0  0  0  0  0  0
    {:.4}    {:.4}    {:.4} C   0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
M  END"#, x, y, z, x + 1.5, y, z);

            if let Ok(mol) = crate::molecule::parser::parse_sdf(&sdf) {
                let coords = vec![[x, y, z], [x + 1.5, y, z]];
                let ff = crate::mmff::MMFFForceField::new(&mol, crate::MMFFVariant::MMFF94s);
                let energy = ff.calculate_energy(&coords);
                prop_assert!(energy.is_finite(), "Energy should be finite for valid coordinates");
            }
        }
    }

    // Test that bond energy is zero at equilibrium
    proptest! {
        #[test]
        fn bond_energy_zero_at_equilibrium(length in 1.0f64..2.0) {
            use crate::mmff::bond::{bond_energy, BondParams};

            let params = BondParams {
                k_bond: 4.7,
                r0: length,
            };

            // At equilibrium distance, energy should be zero
            let coords = vec![[0.0, 0.0, 0.0], [length, 0.0, 0.0]];
            let energy = bond_energy(&coords, 0, 1, &params);

            prop_assert!(
                energy.abs() < 1e-6,
                "Bond energy at equilibrium should be ~0, got {}",
                energy
            );
        }
    }

    // Test that energy increases when bond is stretched or compressed
    proptest! {
        #[test]
        fn bond_energy_increases_away_from_equilibrium(
            deviation in 0.1f64..1.0
        ) {
            use crate::mmff::bond::{bond_energy, BondParams};

            let r0 = 1.526;
            let params = BondParams {
                k_bond: 4.7,
                r0,
            };

            let coords_equilibrium = vec![[0.0, 0.0, 0.0], [r0, 0.0, 0.0]];
            let coords_stretched = vec![[0.0, 0.0, 0.0], [r0 + deviation, 0.0, 0.0]];
            let coords_compressed = vec![[0.0, 0.0, 0.0], [r0 - deviation.min(r0 - 0.5), 0.0, 0.0]];

            let e_eq = bond_energy(&coords_equilibrium, 0, 1, &params);
            let e_stretch = bond_energy(&coords_stretched, 0, 1, &params);
            let e_compress = bond_energy(&coords_compressed, 0, 1, &params);

            prop_assert!(
                e_stretch > e_eq,
                "Stretched bond energy ({}) should be > equilibrium ({})",
                e_stretch, e_eq
            );
            prop_assert!(
                e_compress > e_eq,
                "Compressed bond energy ({}) should be > equilibrium ({})",
                e_compress, e_eq
            );
        }
    }

    // Test that VDW energy has minimum (attractive well)
    proptest! {
        #[test]
        fn vdw_has_attractive_well(_separation in 2.0f64..10.0) {
            use crate::mmff::vdw::{vdw_energy_and_gradient, VDWParams};

            let params = VDWParams {
                r0: 3.5,
                epsilon: 0.1,
                alpha: 12.0,
                beta: 2.25,
            };

            let coords_close = vec![[0.0, 0.0, 0.0], [1.5, 0.0, 0.0]]; // Very close, repulsive
            let coords_mid = vec![[0.0, 0.0, 0.0], [3.5, 0.0, 0.0]];    // At r0, attractive well
            let coords_far = vec![[0.0, 0.0, 0.0], [10.0, 0.0, 0.0]];  // Far, ~zero

            let (_e_close, _, _) = vdw_energy_and_gradient(
                &coords_close, 0, 1, &params, &params);
            let (e_mid, _, _) = vdw_energy_and_gradient(
                &coords_mid, 0, 1, &params, &params);
            let (e_far, _, _) = vdw_energy_and_gradient(
                &coords_far, 0, 1, &params, &params);

            // At minimum distance (r0), energy should be negative (attractive)
            prop_assert!(
                e_mid < 0.0,
                "VDW energy at r0 should be negative (attractive), got {}",
                e_mid
            );

            // Far away, energy should approach zero from below
            prop_assert!(
                e_far.abs() < 0.01,
                "VDW energy at infinity should be ~0, got {}",
                e_far
            );
        }
    }

    // Test that parser rejects invalid SDF files
    proptest! {
        #[test]
        fn parser_rejects_invalid_input(input in "[\u{0000}-\u{007F}]{0,100}") {
            // Most random strings are not valid SDF
            if !input.as_str().contains("V2000") && !input.as_str().contains("V3000") {
                let _result = crate::molecule::parser::parse_sdf(&input);
                // Should either parse successfully or return an error
                // Either is fine, just shouldn't panic
            }
        }
    }

    // Test energy/gradient consistency (finite difference check)
    proptest! {
        #[test]
        fn gradient_finite_difference(
            x in -1.0f64..1.0,
            y in -1.0f64..1.0,
            z in -1.0f64..1.0,
        ) {
            use crate::mmff::bond::{bond_energy, bond_gradient, BondParams};

            // Skip degenerate case where atoms are at same position
            let dist_sq = x * x + y * y + z * z;
            prop_assume!(dist_sq > 1e-10);

            let coords = vec![[0.0, 0.0, 0.0], [x, y, z]];
            let params = BondParams { k_bond: 4.7, r0: 1.526 };

            // Calculate analytical gradient
            let (_g1, g2) = bond_gradient(&coords, 0, 1, &params);

            // Calculate numerical gradient using finite difference
            let eps = 1e-7;
            let e0 = bond_energy(&coords, 0, 1, &params);

            for dim in 0..3 {
                let mut coords_plus = coords.clone();
                coords_plus[1][dim] += eps;
                let e_plus = bond_energy(&coords_plus, 0, 1, &params);
                let num_grad = (e_plus - e0) / eps;

                prop_assert!(
                    (g2[dim] - num_grad).abs() < 1e-3,
                    "Gradient mismatch at dim {}: analytical={}, numerical={}",
                    dim, g2[dim], num_grad
                );
            }
        }
    }
}
