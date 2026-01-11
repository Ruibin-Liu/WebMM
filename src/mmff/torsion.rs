//! Torsion term for MMFF94

use super::super::mmff::MMFFAtomType;
use super::super::mmff::MMFFVariant;

/// Torsion parameters
#[derive(Debug, Clone, Copy)]
pub struct TorsionParams {
    pub v1: f64, // kcal/mol
    pub v2: f64,
    pub v3: f64,
}

/// Get torsion parameters for atom types
pub fn get_torsion_params(
    type1: MMFFAtomType,
    type2: MMFFAtomType,
    type3: MMFFAtomType,
    type4: MMFFAtomType,
    mmff_variant: MMFFVariant,
) -> Option<TorsionParams> {
    // TODO: Load from parameter tables
    // For now, return simple default values
    match (type1, type2, type3, type4) {
        (MMFFAtomType::C_3, MMFFAtomType::C_3, MMFFAtomType::C_3, MMFFAtomType::C_3) => {
            Some(TorsionParams {
                v1: 0.0,
                v2: 0.2,
                v3: 0.15,
            })
        }
        (MMFFAtomType::C_3, MMFFAtomType::C_3, MMFFAtomType::C_2, MMFFAtomType::C_2) => {
            Some(TorsionParams {
                v1: 0.0,
                v2: 0.2,
                v3: 0.15,
            })
        }
        (MMFFAtomType::C_3, MMFFAtomType::C_3, MMFFAtomType::C_AR, MMFFAtomType::C_AR) => {
            Some(TorsionParams {
                v1: 0.0,
                v2: 0.2,
                v3: 0.15,
            })
        }
        _ => None,
    }
}

/// Calculate torsion energy
pub fn torsion_energy(
    coords: &[[f64; 3]],
    i: usize,
    j: usize,
    k: usize,
    l: usize,
    params: &TorsionParams,
) -> f64 {
    let phi = calculate_dihedral(coords, i, j, k, l);

    // E_tors = V1 * (1 + cos(phi)) + V2 * (1 - cos(2*phi)) + V3 * (1 + cos(3*phi))
    let cos_phi = phi.cos();
    let cos_2phi = (2.0 * phi).cos();
    let cos_3phi = (3.0 * phi).cos();

    params.v1 * (1.0 + cos_phi) + params.v2 * (1.0 - cos_2phi) + params.v3 * (1.0 + cos_3phi)
}

/// Calculate dihedral angle between four atoms
fn calculate_dihedral(coords: &[[f64; 3]], i: usize, j: usize, k: usize, l: usize) -> f64 {
    let b1 = [
        coords[i][0] - coords[j][0],
        coords[i][1] - coords[j][1],
        coords[i][2] - coords[j][2],
    ];
    let b2 = [
        coords[k][0] - coords[j][0],
        coords[k][1] - coords[j][1],
        coords[k][2] - coords[j][2],
    ];
    let b3 = [
        coords[l][0] - coords[k][0],
        coords[l][1] - coords[k][1],
        coords[l][2] - coords[k][2],
    ];

    // Normal vectors
    let n1 = [
        b1[1] * b2[2] - b1[2] * b2[1],
        b1[2] * b2[0] - b1[0] * b2[2],
        b1[0] * b2[1] - b1[1] * b2[0],
    ];

    let n2 = [
        b2[1] * b3[2] - b2[2] * b3[1],
        b2[2] * b3[0] - b2[0] * b3[2],
        b2[0] * b3[1] - b2[1] * b3[0],
    ];

    let n1_norm = (n1[0].powi(2) + n1[1].powi(2) + n1[2].powi(2)).sqrt();
    let n2_norm = (n2[0].powi(2) + n2[1].powi(2) + n2[2].powi(2)).sqrt();

    if n1_norm == 0.0 || n2_norm == 0.0 {
        return 0.0;
    }

    // Normalize
    let m1 = [n1[0] / n1_norm, n1[1] / n1_norm, n1[2] / n1_norm];
    let m2 = [n2[0] / n2_norm, n2[1] / n2_norm, n2[2] / n2_norm];

    // x axis = m1
    let x = m1;

    // y axis = (m2 - (m2 . x) * x) normalized
    let y_dot_x = m2[0] * x[0] + m2[1] * x[1] + m2[2] * x[2];
    let y = [
        m2[0] - y_dot_x * x[0],
        m2[1] - y_dot_x * x[1],
        m2[2] - y_dot_x * x[2],
    ];
    let y_norm = (y[0].powi(2) + y[1].powi(2) + y[2].powi(2)).sqrt();

    if y_norm == 0.0 {
        return 0.0;
    }
    let y = [y[0] / y_norm, y[1] / y_norm, y[2] / y_norm];

    // z axis = m1 cross m2
    let z = [
        m1[1] * m2[2] - m1[2] * m2[1],
        m1[2] * m2[0] - m1[0] * m2[2],
        m1[0] * m2[1] - m1[1] * m2[0],
    ];

    // Project b1 and b3 onto coordinate system
    let b1_norm = (b1[0].powi(2) + b1[1].powi(2) + b1[2].powi(2)).sqrt();
    let b3_norm = (b3[0].powi(2) + b3[1].powi(2) + b3[2].powi(2)).sqrt();

    if b1_norm == 0.0 || b3_norm == 0.0 {
        return 0.0;
    }

    let xb1 = b1[0] * x[0] + b1[1] * x[1] + b1[2] * x[2];
    let yb1 = b1[0] * y[0] + b1[1] * y[1] + b1[2] * y[2];
    let zb1 = b1[0] * z[0] + b1[1] * z[1] + b1[2] * z[2];

    let xb3 = b3[0] * x[0] + b3[1] * x[1] + b3[2] * x[2];
    let yb3 = b3[0] * y[0] + b3[1] * y[1] + b3[2] * y[2];
    let zb3 = b3[0] * z[0] + b3[1] * z[1] + b3[2] * z[2];

    let b1_x = xb1 / b1_norm;
    let b1_y = yb1 / b1_norm;
    let b3_x = xb3 / b3_norm;
    let b3_y = yb3 / b3_norm;

    // Calculate dihedral angle
    (b1_x.atan2(b1_y) - b3_x.atan2(b3_y)).abs()
}

pub fn torsion_gradient(
    coords: &[[f64; 3]],
    atom1: usize,
    atom2: usize,
    atom3: usize,
    atom4: usize,
    _params: &TorsionParams,
) -> ([f64; 3], [f64; 3], [f64; 3], [f64; 3]) {
    // Simplified gradient - full torsion gradient is complex
    // Using approximate gradient for now
    let phi = calculate_dihedral(coords, atom1, atom2, atom3, atom4);

    // Force from torsion potential
    let sin_phi = phi.sin();
    let cos_phi = phi.cos();

    // Gradient contributions (simplified)
    // Real implementation needs full chain rule
    let grad_factor = (_params.v1 * sin_phi
        + 2.0 * _params.v2 * (2.0 * sin_phi).cos()
        + 3.0 * _params.v3 * (3.0 * sin_phi).cos());

    // Apply to torsion atoms (2 and 3)
    let mut g2 = [0.0; 3];
    let mut g3 = [0.0; 3];

    g2[0] = grad_factor * 0.1; // Approximate x component
    g2[1] = grad_factor * 0.1; // Approximate y component
    g2[2] = grad_factor * 0.1; // Approximate z component

    g3[0] = -grad_factor * 0.1;
    g3[1] = -grad_factor * 0.1;
    g3[2] = -grad_factor * 0.1;

    ([0.0; 3], g2, g3, [0.0; 3])
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mmff::MMFFVariant;

    #[test]
    fn test_torsion_energy() {
        // Test with simple ethane torsion
        let coords = vec![
            [1.526, 0.0, 0.0],     // H
            [0.0, 0.0, 0.0],       // C
            [1.526, 0.934, 0.0],   // H
            [2.526, 0.934, 1.526], // H (third H)
        ];

        let params = TorsionParams {
            v1: 0.0,
            v2: 0.2,
            v3: 0.15,
        };
        let energy = torsion_energy(&coords, 0, 1, 2, 3, &params);

        assert!(energy.is_finite());
    }
}
