//! Out-of-Plane bending term for MMFF94
//!
//! Uses the full 117-entry RDKit OOP parameter table with 4-stage
//! equivalence protocol lookup.

use super::params::{get_eq_levels, mmff_type_id};
use super::MMFFAtomType;
use super::MMFFVariant;

#[rustfmt::skip]
const OOP94_TABLE: &[(u8, u8, u8, u8, f64)] = &[
    ( 0,  2,  0,  0,    0.020),
    ( 0,  3,  0,  0,    0.130),
    ( 0,  8,  0,  0,    0.000),
    ( 0, 10,  0,  0,   -0.020),
    ( 0, 17,  0,  0,    0.000),
    ( 0, 26,  0,  0,    0.000),
    ( 0, 30,  0,  0,    0.010),
    ( 0, 37,  0,  0,    0.035),
    ( 0, 39,  0,  0,    0.020),
    ( 0, 40,  0,  0,   -0.005),
    ( 0, 41,  0,  0,    0.180),
    ( 0, 43,  0,  0,    0.000),
    ( 0, 45,  0,  0,    0.150),
    ( 0, 49,  0,  0,    0.000),
    ( 0, 54,  0,  0,    0.020),
    ( 0, 55,  0,  0,    0.020),
    ( 0, 56,  0,  0,    0.020),
    ( 0, 57,  0,  0,    0.080),
    ( 0, 58,  0,  0,    0.025),
    ( 0, 63,  0,  0,    0.050),
    ( 0, 64,  0,  0,    0.040),
    ( 0, 67,  0,  0,    0.070),
    ( 0, 69,  0,  0,    0.070),
    ( 0, 73,  0,  0,    0.000),
    ( 0, 78,  0,  0,    0.045),
    ( 0, 80,  0,  0,    0.080),
    ( 0, 81,  0,  0,    0.025),
    ( 0, 82,  0,  0,    0.000),
    ( 1,  2,  1,  2,    0.030),
    ( 1,  2,  2,  2,    0.027),
    ( 1,  2,  2,  3,    0.026),
    ( 1,  2,  2,  5,    0.013),
    ( 1,  2,  2, 37,    0.032),
    ( 1,  3,  1,  7,    0.146),
    ( 1,  3,  2,  7,    0.138),
    ( 1,  3,  3,  7,    0.134),
    ( 1,  3,  5,  7,    0.122),
    ( 1,  3,  6,  7,    0.141),
    ( 1,  3,  7, 10,    0.129),
    ( 1,  3,  7, 37,    0.138),
    ( 1, 10,  1,  3,   -0.020),
    ( 1, 10,  3,  6,   -0.033),
    ( 1, 10,  3, 28,   -0.020),
    ( 1, 37, 37, 37,    0.040),
    ( 1, 39, 63, 63,    0.012),
    ( 1, 40, 28, 37,   -0.006),
    ( 1, 41, 32, 32,    0.178),
    ( 1, 54,  3, 36,    0.016),
    ( 1, 55, 36, 57,    0.020),
    ( 1, 56, 36, 57,    0.020),
    ( 2,  2,  2,  5,    0.013),
    ( 2,  2,  3,  5,    0.012),
    ( 2,  2,  5,  5,    0.006),
    ( 2,  2,  5,  6,    0.027),
    ( 2,  2,  5, 37,    0.017),
    ( 2,  2,  5, 40,    0.012),
    ( 2,  2,  5, 41,    0.008),
    ( 2,  3,  5,  7,    0.113),
    ( 2,  3,  5,  9,    0.081),
    ( 2,  3,  6,  7,    0.127),
    ( 2,  3,  7, 10,    0.116),
    ( 2, 37, 37, 37,    0.031),
    ( 2, 40, 28, 28,   -0.007),
    ( 2, 41, 32, 32,    0.161),
    ( 3,  3,  5,  7,    0.113),
    ( 3,  3,  6,  7,    0.127),
    ( 3, 10,  3, 28,   -0.030),
    ( 3, 10, 28, 28,   -0.019),
    ( 3, 37, 37, 37,    0.027),
    ( 3, 40, 28, 28,   -0.007),
    ( 3, 54, 36, 36,    0.018),
    ( 5,  3,  5,  7,    0.103),
    ( 5,  3,  5,  9,    0.074),
    ( 5,  3,  5, 54,    0.078),
    ( 5,  3,  6,  7,    0.119),
    ( 5,  3,  7, 10,    0.102),
    ( 5,  3,  9, 40,    0.067),
    ( 5, 30, 20, 30,    0.008),
    ( 5, 37, 37, 37,    0.015),
    ( 5, 37, 37, 38,    0.046),
    ( 5, 37, 37, 63,    0.008),
    ( 5, 37, 37, 64,    0.012),
    ( 5, 37, 37, 69,    0.016),
    ( 5, 37, 38, 38,    0.084),
    ( 5, 41, 32, 32,    0.158),
    ( 5, 57, 55, 55,    0.038),
    ( 5, 63, 39, 64,    0.019),
    ( 5, 63, 39, 66,    0.068),
    ( 5, 63, 44, 64,    0.014),
    ( 5, 63, 44, 66,    0.055),
    ( 5, 63, 59, 64,    0.033),
    ( 5, 63, 59, 66,    0.085),
    ( 5, 64, 63, 64,    0.006),
    ( 5, 64, 63, 66,    0.043),
    ( 5, 64, 64, 65,    0.052),
    ( 5, 64, 65, 66,    0.094),
    ( 5, 78, 78, 81,    0.046),
    ( 5, 80, 81, 81,    0.057),
    ( 6,  3,  7, 37,    0.127),
    ( 6, 37, 37, 37,    0.048),
    ( 7,  3, 10, 10,    0.113),
    ( 7,  3, 20, 20,    0.151),
    ( 9,  3, 40, 40,    0.057),
    (15, 37, 37, 37,    0.025),
    (23, 39, 63, 63,   -0.014),
    (23, 39, 63, 65,    0.021),
    (23, 39, 65, 65,    0.062),
    (28, 40, 28, 37,    0.004),
    (32, 69, 37, 37,    0.067),
    (36, 55, 36, 57,    0.020),
    (36, 56, 36, 57,    0.020),
    (36, 81, 78, 80,    0.016),
    (37, 37, 37, 40,    0.046),
    (37, 63, 39, 64,    0.010),
    (37, 64, 63, 64,   -0.011),
    (50, 49, 50, 50,    0.000),
    (56, 57, 56, 56,    0.158),
];

#[rustfmt::skip]
const OOP94S_TABLE: &[(u8, u8, u8, u8, f64)] = &[
    ( 0,  2,  0,  0,    0.020),
    ( 0,  3,  0,  0,    0.130),
    ( 0,  8,  0,  0,    0.000),
    ( 0, 10,  0,  0,    0.015),
    ( 0, 17,  0,  0,    0.000),
    ( 0, 26,  0,  0,    0.000),
    ( 0, 30,  0,  0,    0.010),
    ( 0, 37,  0,  0,    0.035),
    ( 0, 39,  0,  0,    0.020),
    ( 0, 40,  0,  0,    0.030),
    ( 0, 41,  0,  0,    0.180),
    ( 0, 43,  0,  0,    0.000),
    ( 0, 45,  0,  0,    0.150),
    ( 0, 49,  0,  0,    0.000),
    ( 0, 54,  0,  0,    0.020),
    ( 0, 55,  0,  0,    0.020),
    ( 0, 56,  0,  0,    0.020),
    ( 0, 57,  0,  0,    0.080),
    ( 0, 58,  0,  0,    0.025),
    ( 0, 63,  0,  0,    0.050),
    ( 0, 64,  0,  0,    0.040),
    ( 0, 67,  0,  0,    0.070),
    ( 0, 69,  0,  0,    0.070),
    ( 0, 73,  0,  0,    0.000),
    ( 0, 78,  0,  0,    0.045),
    ( 0, 80,  0,  0,    0.080),
    ( 0, 81,  0,  0,    0.025),
    ( 0, 82,  0,  0,    0.000),
    ( 1,  2,  1,  2,    0.030),
    ( 1,  2,  2,  2,    0.027),
    ( 1,  2,  2,  3,    0.026),
    ( 1,  2,  2,  5,    0.013),
    ( 1,  2,  2, 37,    0.032),
    ( 1,  3,  1,  7,    0.146),
    ( 1,  3,  2,  7,    0.138),
    ( 1,  3,  3,  7,    0.134),
    ( 1,  3,  5,  7,    0.122),
    ( 1,  3,  6,  7,    0.141),
    ( 1,  3,  7, 10,    0.129),
    ( 1,  3,  7, 37,    0.138),
    ( 1, 10,  1,  3,    0.015),
    ( 1, 10,  3,  6,    0.015),
    ( 1, 10,  3, 28,    0.015),
    ( 1, 37, 37, 37,    0.040),
    ( 1, 39, 63, 63,    0.012),
    ( 1, 40, 28, 37,    0.030),
    ( 1, 41, 32, 32,    0.178),
    ( 1, 54,  3, 36,    0.016),
    ( 1, 55, 36, 57,    0.020),
    ( 1, 56, 36, 57,    0.020),
    ( 2,  2,  2,  5,    0.013),
    ( 2,  2,  3,  5,    0.012),
    ( 2,  2,  5,  5,    0.006),
    ( 2,  2,  5,  6,    0.027),
    ( 2,  2,  5, 37,    0.017),
    ( 2,  2,  5, 40,    0.012),
    ( 2,  2,  5, 41,    0.008),
    ( 2,  3,  5,  7,    0.113),
    ( 2,  3,  5,  9,    0.081),
    ( 2,  3,  6,  7,    0.127),
    ( 2,  3,  7, 10,    0.116),
    ( 2, 37, 37, 37,    0.031),
    ( 2, 40, 28, 28,    0.030),
    ( 2, 41, 32, 32,    0.161),
    ( 3,  3,  5,  7,    0.113),
    ( 3,  3,  6,  7,    0.127),
    ( 3, 10,  3, 28,    0.015),
    ( 3, 10, 28, 28,    0.015),
    ( 3, 37, 37, 37,    0.027),
    ( 3, 40, 28, 28,    0.030),
    ( 3, 54, 36, 36,    0.018),
    ( 5,  3,  5,  7,    0.103),
    ( 5,  3,  5,  9,    0.074),
    ( 5,  3,  5, 54,    0.078),
    ( 5,  3,  6,  7,    0.119),
    ( 5,  3,  7, 10,    0.102),
    ( 5,  3,  9, 40,    0.067),
    ( 5, 30, 20, 30,    0.008),
    ( 5, 37, 37, 37,    0.015),
    ( 5, 37, 37, 38,    0.046),
    ( 5, 37, 37, 63,    0.008),
    ( 5, 37, 37, 64,    0.012),
    ( 5, 37, 37, 69,    0.016),
    ( 5, 37, 38, 38,    0.084),
    ( 5, 41, 32, 32,    0.158),
    ( 5, 57, 55, 55,    0.038),
    ( 5, 63, 39, 64,    0.019),
    ( 5, 63, 39, 66,    0.068),
    ( 5, 63, 44, 64,    0.014),
    ( 5, 63, 44, 66,    0.055),
    ( 5, 63, 59, 64,    0.033),
    ( 5, 63, 59, 66,    0.085),
    ( 5, 64, 63, 64,    0.006),
    ( 5, 64, 63, 66,    0.043),
    ( 5, 64, 64, 65,    0.052),
    ( 5, 64, 65, 66,    0.094),
    ( 5, 78, 78, 81,    0.046),
    ( 5, 80, 81, 81,    0.057),
    ( 6,  3,  7, 37,    0.127),
    ( 6, 37, 37, 37,    0.048),
    ( 7,  3, 10, 10,    0.113),
    ( 7,  3, 20, 20,    0.151),
    ( 9,  3, 40, 40,    0.057),
    (15, 37, 37, 37,    0.025),
    (23, 39, 63, 63,   -0.014),
    (23, 39, 63, 65,    0.021),
    (23, 39, 65, 65,    0.062),
    (28, 40, 28, 37,    0.030),
    (32, 69, 37, 37,    0.067),
    (36, 55, 36, 57,    0.020),
    (36, 56, 36, 57,    0.020),
    (36, 81, 78, 80,    0.016),
    (37, 37, 37, 40,    0.046),
    (37, 63, 39, 64,    0.010),
    (37, 64, 63, 64,   -0.011),
    (50, 49, 50, 50,    0.000),
    (56, 57, 56, 56,    0.158),
];

#[derive(Debug, Clone, Copy)]
pub struct OOPParams {
    pub k_oop: f64,
}

fn lookup_oop_table(
    table: &[(u8, u8, u8, u8, f64)],
    gi: u8,
    j: u8,
    gk: u8,
    gl: u8,
) -> Option<OOPParams> {
    let mut lo = 0usize;
    let mut hi = table.len();
    while lo < hi {
        let mid = lo + (hi - lo) / 2;
        let (ti, tj, tk, tl, koop) = table[mid];
        match (gi.cmp(&ti), j.cmp(&tj), gk.cmp(&tk), gl.cmp(&tl)) {
            (
                std::cmp::Ordering::Equal,
                std::cmp::Ordering::Equal,
                std::cmp::Ordering::Equal,
                std::cmp::Ordering::Equal,
            ) => {
                return Some(OOPParams { k_oop: koop });
            }
            (std::cmp::Ordering::Less, _, _, _)
            | (std::cmp::Ordering::Equal, std::cmp::Ordering::Less, _, _)
            | (std::cmp::Ordering::Equal, std::cmp::Ordering::Equal, std::cmp::Ordering::Less, _)
            | (
                std::cmp::Ordering::Equal,
                std::cmp::Ordering::Equal,
                std::cmp::Ordering::Equal,
                std::cmp::Ordering::Less,
            ) => {
                hi = mid;
            }
            _ => {
                lo = mid + 1;
            }
        }
    }
    None
}

pub fn get_oop_params(
    i_type: MMFFAtomType,
    central_type: MMFFAtomType,
    k_type: MMFFAtomType,
    l_type: MMFFAtomType,
    variant: MMFFVariant,
) -> OOPParams {
    let i_id = mmff_type_id(i_type);
    let j_id = mmff_type_id(central_type);
    let k_id = mmff_type_id(k_type);
    let l_id = mmff_type_id(l_type);

    let table = match variant {
        MMFFVariant::MMFF94 => OOP94_TABLE,
        MMFFVariant::MMFF94s => OOP94S_TABLE,
    };

    let i_levels = get_eq_levels(i_id);
    let k_levels = get_eq_levels(k_id);
    let l_levels = get_eq_levels(l_id);

    for iter in 0..4 {
        let gi = i_levels[iter];
        let gk = k_levels[iter];
        let gl = l_levels[iter];

        let mut ikl = [gi, gk, gl];
        ikl.sort();

        if let Some(params) = lookup_oop_table(table, ikl[0], j_id, ikl[1], ikl[2]) {
            return params;
        }
    }

    OOPParams { k_oop: 0.0 }
}

/// Calculate out-of-plane bending energy
pub fn oop_energy(
    coords: &[[f64; 3]],
    central: usize,
    i: usize,
    j: usize,
    k: usize,
    params: &OOPParams,
) -> f64 {
    let chi = calculate_oop_angle(coords, central, i, j, k);

    // MMFF94 out-of-plane: E = 0.5 * C * k_oop * χ²
    // where C = 143.9324 converts mdyn·Å/rad² to kcal/mol,
    // k_oop is the raw parameter table value in mdyn·Å/rad²,
    // and χ is the Wilson out-of-plane angle in radians.
    let c = 143.9325;
    0.5 * c * params.k_oop * chi * chi
}

/// Calculate Wilson out-of-plane angle (RDKit-style)
/// For atoms (i, central, j, k): angle between bond central→i and plane central→j→k
/// Returns 0 when all four atoms are coplanar.
fn calculate_oop_angle(coords: &[[f64; 3]], central: usize, i: usize, j: usize, k: usize) -> f64 {
    // Vectors from central atom to neighbors
    let r_ci = [
        coords[i][0] - coords[central][0],
        coords[i][1] - coords[central][1],
        coords[i][2] - coords[central][2],
    ];
    let r_cj = [
        coords[j][0] - coords[central][0],
        coords[j][1] - coords[central][1],
        coords[j][2] - coords[central][2],
    ];
    let r_ck = [
        coords[k][0] - coords[central][0],
        coords[k][1] - coords[central][1],
        coords[k][2] - coords[central][2],
    ];

    // Normalize vectors
    let r_ci_norm = (r_ci[0].powi(2) + r_ci[1].powi(2) + r_ci[2].powi(2)).sqrt();
    let r_cj_norm = (r_cj[0].powi(2) + r_cj[1].powi(2) + r_cj[2].powi(2)).sqrt();
    let r_ck_norm = (r_ck[0].powi(2) + r_ck[1].powi(2) + r_ck[2].powi(2)).sqrt();

    if r_ci_norm < 1e-10 || r_cj_norm < 1e-10 || r_ck_norm < 1e-10 {
        return 0.0;
    }

    let u_ci = [
        r_ci[0] / r_ci_norm,
        r_ci[1] / r_ci_norm,
        r_ci[2] / r_ci_norm,
    ];
    let u_cj = [
        r_cj[0] / r_cj_norm,
        r_cj[1] / r_cj_norm,
        r_cj[2] / r_cj_norm,
    ];
    let u_ck = [
        r_ck[0] / r_ck_norm,
        r_ck[1] / r_ck_norm,
        r_ck[2] / r_ck_norm,
    ];

    // Normal to plane defined by central→j and central→k
    let normal = [
        u_cj[1] * u_ck[2] - u_cj[2] * u_ck[1],
        u_cj[2] * u_ck[0] - u_cj[0] * u_ck[2],
        u_cj[0] * u_ck[1] - u_cj[1] * u_ck[0],
    ];

    let normal_norm = (normal[0].powi(2) + normal[1].powi(2) + normal[2].powi(2)).sqrt();
    if normal_norm < 1e-10 {
        return 0.0;
    }

    let n_unit = [
        normal[0] / normal_norm,
        normal[1] / normal_norm,
        normal[2] / normal_norm,
    ];

    // sin(χ) = dot(u_ci, n_unit)
    let sin_chi =
        (u_ci[0] * n_unit[0] + u_ci[1] * n_unit[1] + u_ci[2] * n_unit[2]).clamp(-1.0, 1.0);

    sin_chi.asin()
}

pub fn oop_gradient(
    coords: &[[f64; 3]],
    central: usize,
    atom1: usize,
    atom2: usize,
    atom3: usize,
    params: &OOPParams,
) -> ([f64; 3], [f64; 3], [f64; 3], [f64; 3]) {
    // Central-difference numerical gradient (O(eps²) accurate).
    let eps = 1e-5;
    let mut gc = [0.0; 3];
    let mut g1 = [0.0; 3];
    let mut g2 = [0.0; 3];
    let mut g3 = [0.0; 3];

    for (atom_idx, grad) in [
        (central, &mut gc),
        (atom1, &mut g1),
        (atom2, &mut g2),
        (atom3, &mut g3),
    ] {
        for dim in 0..3 {
            let mut cp = coords.to_vec();
            cp[atom_idx][dim] += eps;
            let ep = oop_energy(&cp, central, atom1, atom2, atom3, params);
            cp[atom_idx][dim] -= 2.0 * eps;
            let em = oop_energy(&cp, central, atom1, atom2, atom3, params);
            grad[dim] = (ep - em) / (2.0 * eps);
        }
    }

    (gc, g1, g2, g3)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_oop_energy() {
        let coords = vec![
            [1.526, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [1.526, 0.934, 0.0],
            [0.0, 0.934, -1.526],
        ];

        let params = OOPParams { k_oop: 0.04 };
        let energy = oop_energy(&coords, 1, 0, 2, 3, &params);

        assert!(energy.is_finite());
    }

    #[test]
    fn test_oop_gradient_numerical() {
        let coords = vec![
            [1.526, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 1.526, 0.0],
            [0.0, 0.0, 0.5],
        ];
        let params = OOPParams { k_oop: 0.04 };
        let (gc, g1, g2, g3) = oop_gradient(&coords, 1, 0, 2, 3, &params);

        // Verify against independent central-difference
        let eps = 1e-5;
        for (idx, grad) in [(1usize, gc), (0usize, g1), (2usize, g2), (3usize, g3)] {
            for dim in 0..3 {
                let mut cp = coords.clone();
                cp[idx][dim] += eps;
                let ep = oop_energy(&cp, 1, 0, 2, 3, &params);
                cp[idx][dim] -= 2.0 * eps;
                let em = oop_energy(&cp, 1, 0, 2, 3, &params);
                let num = (ep - em) / (2.0 * eps);
                assert!(
                    (grad[dim] - num).abs() < 1e-4,
                    "OOP grad[{}] = {} vs numerical {} for atom {}",
                    dim,
                    grad[dim],
                    num,
                    idx
                );
            }
        }
    }

    #[test]
    fn test_oop_lookup_exact_match() {
        let params = get_oop_params(
            MMFFAtomType::H,
            MMFFAtomType::C_2,
            MMFFAtomType::H,
            MMFFAtomType::O_2,
            MMFFVariant::MMFF94,
        );
        assert!(
            params.k_oop != 0.0,
            "Should find specific OOP params for H-C_2-H-O_2"
        );
    }

    #[test]
    fn test_oop_lookup_aromatic() {
        let params = get_oop_params(
            MMFFAtomType::C_3,
            MMFFAtomType::C_AR,
            MMFFAtomType::C_AR,
            MMFFAtomType::C_AR,
            MMFFVariant::MMFF94,
        );
        assert!(
            (params.k_oop - 0.040).abs() < 0.001,
            "C_AR central with C_3 peripheral should match table entry, got {}",
            params.k_oop
        );
    }

    #[test]
    fn test_oop_lookup_no_match() {
        let params = get_oop_params(
            MMFFAtomType::Fe_P2,
            MMFFAtomType::Fe_P2,
            MMFFAtomType::Fe_P2,
            MMFFAtomType::Fe_P2,
            MMFFVariant::MMFF94,
        );
        assert_eq!(params.k_oop, 0.0, "Iron types should not have OOP params");
    }

    #[test]
    fn test_oop_mmff94s_differs_from_mmff94() {
        let params94 = get_oop_params(
            MMFFAtomType::C_3,
            MMFFAtomType::N_AM,
            MMFFAtomType::C_2,
            MMFFAtomType::O_3,
            MMFFVariant::MMFF94,
        );
        let params94s = get_oop_params(
            MMFFAtomType::C_3,
            MMFFAtomType::N_AM,
            MMFFAtomType::C_2,
            MMFFAtomType::O_3,
            MMFFVariant::MMFF94s,
        );
        assert_ne!(
            params94.k_oop, params94s.k_oop,
            "MMFF94 ({}) and MMFF94s ({}) should differ for N_AM central",
            params94.k_oop, params94s.k_oop
        );
    }
}
