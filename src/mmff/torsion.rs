//! Torsion term for MMFF94

use super::MMFFAtomType;
use super::MMFFVariant;

/// Torsion parameters
#[derive(Debug, Clone, Copy)]
pub struct TorsionParams {
    pub v1: f64,
    pub v2: f64,
    pub v3: f64,
}

/// Get torsion parameters for atom types
pub fn get_torsion_params(
    type1: MMFFAtomType,
    type2: MMFFAtomType,
    type3: MMFFAtomType,
    type4: MMFFAtomType,
    _mmff_variant: MMFFVariant,
) -> Option<TorsionParams> {
    match (type1, type2, type3, type4) {
        (MMFFAtomType::C_3, MMFFAtomType::C_3, MMFFAtomType::C_2, MMFFAtomType::O_2) => {
            Some(TorsionParams {
                v1: 0.0,
                v2: 0.15,
                v3: 0.0,
            })
        }
        (_, MMFFAtomType::C_3, MMFFAtomType::C_2, _) => Some(TorsionParams {
            v1: 0.0,
            v2: 0.2,
            v3: 0.15,
        }),

        (MMFFAtomType::C_3, MMFFAtomType::N_3, MMFFAtomType::C_2, MMFFAtomType::O_2) => {
            Some(TorsionParams {
                v1: 0.0,
                v2: 0.15,
                v3: 0.0,
            })
        }
        (_, MMFFAtomType::N_3, MMFFAtomType::C_2, _) => Some(TorsionParams {
            v1: 0.0,
            v2: 0.2,
            v3: 0.15,
        }),

        (_, MMFFAtomType::C_3, MMFFAtomType::N_PL3, _) => Some(TorsionParams {
            v1: 0.0,
            v2: 10.0,
            v3: 0.5,
        }),

        (_, MMFFAtomType::C_3, MMFFAtomType::N_AM, _) => Some(TorsionParams {
            v1: 0.0,
            v2: 10.0,
            v3: 0.5,
        }),

        (_, MMFFAtomType::C_AR, MMFFAtomType::O_R, _) => Some(TorsionParams {
            v1: 0.0,
            v2: 2.0,
            v3: 0.0,
        }),

        (_, MMFFAtomType::C_2, MMFFAtomType::C_2, _) => Some(TorsionParams {
            v1: 0.0,
            v2: 10.0,
            v3: 0.0,
        }),

        (_, MMFFAtomType::C_AR, MMFFAtomType::C_AR, _) => Some(TorsionParams {
            v1: 0.0,
            v2: 0.0,
            v3: 0.0,
        }),

        (_, MMFFAtomType::C_AR, MMFFAtomType::N_AR, _) => Some(TorsionParams {
            v1: 0.0,
            v2: 0.0,
            v3: 0.0,
        }),

        (_, MMFFAtomType::C_3, MMFFAtomType::C_3, _) => Some(TorsionParams {
            v1: 0.0,
            v2: 0.2,
            v3: 0.15,
        }),

        (_, MMFFAtomType::C_3, MMFFAtomType::C_AR, _) => Some(TorsionParams {
            v1: 0.0,
            v2: 0.2,
            v3: 0.15,
        }),

        (_, MMFFAtomType::C_3, MMFFAtomType::N_3, _) => Some(TorsionParams {
            v1: 0.0,
            v2: 0.2,
            v3: 0.0,
        }),

        (_, MMFFAtomType::C_3, MMFFAtomType::O_3, _) => Some(TorsionParams {
            v1: 0.0,
            v2: 0.5,
            v3: 0.0,
        }),

        (_, MMFFAtomType::C_3, MMFFAtomType::O_R, _) => Some(TorsionParams {
            v1: 0.0,
            v2: 0.5,
            v3: 0.0,
        }),

        (_, MMFFAtomType::C_3, MMFFAtomType::S_3, _) => Some(TorsionParams {
            v1: 0.0,
            v2: 0.2,
            v3: 0.0,
        }),

        (_, MMFFAtomType::C_3, MMFFAtomType::P_3, _) => Some(TorsionParams {
            v1: 0.0,
            v2: 0.2,
            v3: 0.0,
        }),

        _ => Some(TorsionParams {
            v1: 0.0,
            v2: 0.0,
            v3: 0.0,
        }),
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

    let m1 = [n1[0] / n1_norm, n1[1] / n1_norm, n1[2] / n1_norm];
    let m2 = [n2[0] / n2_norm, n2[1] / n2_norm, n2[2] / n2_norm];

    let x = m1;

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

    let _z = [
        m1[1] * m2[2] - m1[2] * m2[1],
        m1[2] * m2[0] - m1[0] * m2[2],
        m1[0] * m2[1] - m1[1] * m2[0],
    ];

    let b1_norm = (b1[0].powi(2) + b1[1].powi(2) + b1[2].powi(2)).sqrt();
    let b3_norm = (b3[0].powi(2) + b3[1].powi(2) + b3[2].powi(2)).sqrt();

    if b1_norm == 0.0 || b3_norm == 0.0 {
        return 0.0;
    }

    let b1_x = (b1[0] * x[0] + b1[1] * x[1] + b1[2] * x[2]) / b1_norm;
    let b1_y = (b1[0] * y[0] + b1[1] * y[1] + b1[2] * y[2]) / b1_norm;
    let b3_x = (b3[0] * x[0] + b3[1] * x[1] + b3[2] * x[2]) / b3_norm;
    let b3_y = (b3[0] * y[0] + b3[1] * y[1] + b3[2] * y[2]) / b3_norm;

    (b1_x.atan2(b1_y) - b3_x.atan2(b3_y)).abs()
}

pub fn torsion_gradient(
    coords: &[[f64; 3]],
    atom1: usize,
    atom2: usize,
    atom3: usize,
    atom4: usize,
    params: &TorsionParams,
) -> ([f64; 3], [f64; 3], [f64; 3], [f64; 3]) {
    let eps = 1e-7;
    let e_ref = torsion_energy(coords, atom1, atom2, atom3, atom4, params);
    let mut g1 = [0.0; 3];
    let mut g2 = [0.0; 3];
    let mut g3 = [0.0; 3];
    let mut g4 = [0.0; 3];

    for (atom_idx, grad) in [
        (atom1, &mut g1),
        (atom2, &mut g2),
        (atom3, &mut g3),
        (atom4, &mut g4),
    ] {
        for dim in 0..3 {
            let mut coords_p = coords.to_vec();
            coords_p[atom_idx][dim] += eps;
            let e_plus = torsion_energy(&coords_p, atom1, atom2, atom3, atom4, params);
            grad[dim] = (e_plus - e_ref) / eps;
        }
    }

    (g1, g2, g3, g4)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_torsion_energy() {
        let coords = vec![
            [1.526, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [1.526, 0.934, 0.0],
            [2.526, 0.934, 1.526],
        ];

        let params = TorsionParams {
            v1: 0.0,
            v2: 0.2,
            v3: 0.15,
        };
        let energy = torsion_energy(&coords, 0, 1, 2, 3, &params);

        assert!(energy.is_finite());
    }

    #[test]
    fn test_torsion_gradient_numerical() {
        let coords = vec![
            [1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [1.5, 0.0, 0.0],
            [2.5, 1.0, 0.0],
        ];
        let params = TorsionParams {
            v1: 0.0,
            v2: 0.2,
            v3: 0.15,
        };
        let (g1, g2, g3, g4) = torsion_gradient(&coords, 0, 1, 2, 3, &params);

        let eps = 1e-7;
        for (idx, grad) in [(0usize, g1), (1usize, g2), (2usize, g3), (3usize, g4)] {
            for dim in 0..3 {
                let mut cp = coords.clone();
                cp[idx][dim] += eps;
                let ep = torsion_energy(&cp, 0, 1, 2, 3, &params);
                let num = (ep - torsion_energy(&coords, 0, 1, 2, 3, &params)) / eps;
                assert!(
                    (grad[dim] - num).abs() < 1e-4,
                    "Torsion grad[{}] = {} vs numerical {} for atom {}",
                    dim,
                    grad[dim],
                    num,
                    idx
                );
            }
        }
    }
}
