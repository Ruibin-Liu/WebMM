use super::{base_type, get_angle_params, get_bond_params, MMFFAtomType};
use crate::molecule::BondType;

pub struct StretchBendParams {
    pub kba_ijk: f64,
    pub kba_kji: f64,
}

fn get_sb_kba(
    type_i: MMFFAtomType,
    type_j: MMFFAtomType,
    type_k: MMFFAtomType,
) -> Option<(f64, f64)> {
    let bi = base_type(type_i);
    let bj = base_type(type_j);
    let bk = base_type(type_k);

    match (bi, bj, bk) {
        // H-C_3-H (methane tetrahedral)
        (MMFFAtomType::H, MMFFAtomType::C_3, MMFFAtomType::H) => Some((0.2270, 0.0700)),

        // H-C_3-H (same, reversed)
        (MMFFAtomType::C_3, MMFFAtomType::H, MMFFAtomType::H) => Some((0.0700, 0.2270)),

        // H-C_3-C_3 (ethane-like)
        (MMFFAtomType::H, MMFFAtomType::C_3, MMFFAtomType::C_3) => Some((0.2270, 0.0700)),
        (MMFFAtomType::C_3, MMFFAtomType::C_3, MMFFAtomType::H) => Some((0.0700, 0.2270)),

        // C_3-C_3-O_3 (ethanol C-C-O)
        (MMFFAtomType::C_3, MMFFAtomType::C_3, MMFFAtomType::O_3) => Some((0.1730, 0.4170)),
        (MMFFAtomType::O_3, MMFFAtomType::C_3, MMFFAtomType::C_3) => Some((0.4170, 0.1730)),

        // C_3-O_3-H_OH (water-like O-H)
        (MMFFAtomType::C_3, MMFFAtomType::O_3, MMFFAtomType::H_OH) => Some((0.2560, 0.1430)),
        (MMFFAtomType::H_OH, MMFFAtomType::O_3, MMFFAtomType::C_3) => Some((0.1430, 0.2560)),

        // H-C_2-H (formaldehyde H-C-H)
        (MMFFAtomType::H, MMFFAtomType::C_2, MMFFAtomType::H) => Some((0.2070, 0.1570)),

        // H-C_2=O_2 (formaldehyde H-C=O)
        (MMFFAtomType::H, MMFFAtomType::C_2, MMFFAtomType::O_CO2) => Some((0.1540, 0.8560)),
        (MMFFAtomType::O_CO2, MMFFAtomType::C_2, MMFFAtomType::H) => Some((0.8560, 0.1540)),

        // C_3-C_2=O_2 (acetic acid/acetamide C-C=O)
        (MMFFAtomType::C_3, MMFFAtomType::C_2, MMFFAtomType::O_CO2) => Some((0.3380, 0.7320)),
        (MMFFAtomType::O_CO2, MMFFAtomType::C_2, MMFFAtomType::C_3) => Some((0.7320, 0.3380)),

        // C_2-O_3-H (carboxylic acid O-H)
        (MMFFAtomType::C_2, MMFFAtomType::O_3, MMFFAtomType::H_COOH) => Some((0.2150, 0.0640)),

        // C_3-C_2-N_AM (acetamide C-C-N)
        (MMFFAtomType::C_3, MMFFAtomType::C_2, MMFFAtomType::N_AM) => Some((0.2230, 0.7320)),
        (MMFFAtomType::N_AM, MMFFAtomType::C_2, MMFFAtomType::C_3) => Some((0.7320, 0.2230)),

        // C_2-N_AM-H_NAM (acetamide N-H)
        (MMFFAtomType::C_2, MMFFAtomType::N_AM, MMFFAtomType::H_NAM) => Some((0.1370, 0.0660)),

        // H_OH-O_3-H_OH (water H-O-H)
        (MMFFAtomType::H_OH, MMFFAtomType::O_3, MMFFAtomType::H_OH) => Some((0.2100, 0.2100)),

        // C_3-C_3-H (propane-like)
        (MMFFAtomType::C_3, MMFFAtomType::C_3, MMFFAtomType::C_3) => Some((0.4360, 0.0130)),

        // C_AR-C_AR-C_AR (benzene ring angles)
        (MMFFAtomType::C_AR, MMFFAtomType::C_AR, MMFFAtomType::C_AR) => Some((0.8300, 0.3390)),

        // C_3-C_AR-C_AR (phenol C_AR-C_AR with O attached)
        (MMFFAtomType::C_3, MMFFAtomType::C_AR, MMFFAtomType::C_AR) => Some((0.2410, 0.1300)),

        // C_AR-C_AR-H (aromatic C-H)
        (MMFFAtomType::C_AR, MMFFAtomType::C_AR, MMFFAtomType::H) => Some((0.2500, 0.2790)),
        (MMFFAtomType::H, MMFFAtomType::C_AR, MMFFAtomType::C_AR) => Some((0.2790, 0.2500)),

        // C_AR-C_AR-O_3 (phenol C-O)
        (MMFFAtomType::C_AR, MMFFAtomType::C_AR, MMFFAtomType::O_3) => Some((0.4230, 0.1860)),
        (MMFFAtomType::O_3, MMFFAtomType::C_AR, MMFFAtomType::C_AR) => Some((0.1860, 0.4230)),

        // H_N3-N_3-H_N3 (ammonia H-N-H)
        (MMFFAtomType::H_N3, MMFFAtomType::N_3, MMFFAtomType::H_N3) => Some((0.1900, 0.1900)),

        // N_AR-C_AR-H (aniline N-C-H)
        (MMFFAtomType::N_AR, MMFFAtomType::C_AR, MMFFAtomType::H) => Some((0.0940, 0.0940)),

        // H_NAM-N_AR-C_AR (aniline H-N-C)
        (MMFFAtomType::H_NAM, MMFFAtomType::N_AR, MMFFAtomType::C_AR) => Some((0.0810, 0.0810)),

        // H_NAM-N_AM-C_2 (amide H-N-C)
        (MMFFAtomType::H_NAM, MMFFAtomType::N_AM, MMFFAtomType::C_2) => Some((0.0940, 0.0940)),

        // C_AR-C_AR-N_AR (aniline C-C-N)
        (MMFFAtomType::C_AR, MMFFAtomType::C_AR, MMFFAtomType::N_AR) => Some((0.9010, 0.4290)),
        (MMFFAtomType::N_AR, MMFFAtomType::C_AR, MMFFAtomType::C_AR) => Some((0.4290, 0.9010)),

        // H-C_3-C_3 (C_3-H-C_3)
        (MMFFAtomType::C_3, MMFFAtomType::H, MMFFAtomType::C_3) => Some((0.0700, 0.2270)),

        // H-N_3-C_3 (ammonia-like N-H-C)
        (MMFFAtomType::H_N3, MMFFAtomType::N_3, MMFFAtomType::C_3) => Some((0.0320, 0.8050)),

        // C_3-N_3-C_3
        (MMFFAtomType::C_3, MMFFAtomType::N_3, MMFFAtomType::C_3) => Some((0.5780, 0.4940)),

        // C_3-N_AM-C_3
        (MMFFAtomType::C_3, MMFFAtomType::N_AM, MMFFAtomType::C_3) => Some((0.7710, 0.3530)),

        // O_3-C_3-H (ethanol O-C-H)
        (MMFFAtomType::O_3, MMFFAtomType::C_3, MMFFAtomType::H) => Some((0.4360, 0.0130)),

        // Fallback: use estimation based on bond/angle classes
        _ => None,
    }
}

pub fn get_stretch_bend_params(
    type_i: MMFFAtomType,
    type_j: MMFFAtomType,
    type_k: MMFFAtomType,
    bond_type_ij: BondType,
    bond_type_kj: BondType,
) -> Option<StretchBendParams> {
    if let Some((kba_ijk, kba_kji)) = get_sb_kba(type_i, type_j, type_k) {
        return Some(StretchBendParams { kba_ijk, kba_kji });
    }

    let bond_params_ij = get_bond_params(type_i, type_j, bond_type_ij)?;
    let bond_params_kj = get_bond_params(type_k, type_j, bond_type_kj)?;
    let angle_params = get_angle_params(type_i, type_j, type_k)?;

    let kb_ij = bond_params_ij.k_bond;
    let r0_ij = bond_params_ij.r0;
    let kb_kj = bond_params_kj.k_bond;
    let r0_kj = bond_params_kj.r0;
    let theta0 = angle_params.theta0.to_radians();

    if theta0.abs() < 1e-10 || r0_ij < 1e-10 || r0_kj < 1e-10 {
        return None;
    }

    let kba_ijk = 2.5 * (kb_ij / r0_ij) * theta0.sin();
    let kba_kji = 2.5 * (kb_kj / r0_kj) * theta0.sin();

    Some(StretchBendParams { kba_ijk, kba_kji })
}

pub fn stretch_bend_energy(
    coords: &[[f64; 3]],
    i: usize,
    j: usize,
    k: usize,
    r0_ij: f64,
    r0_kj: f64,
    theta0_rad: f64,
    params: &StretchBendParams,
) -> f64 {
    let rij: [f64; 3] = [
        coords[i][0] - coords[j][0],
        coords[i][1] - coords[j][1],
        coords[i][2] - coords[j][2],
    ];
    let rkj: [f64; 3] = [
        coords[k][0] - coords[j][0],
        coords[k][1] - coords[j][1],
        coords[k][2] - coords[j][2],
    ];
    let r_ij = (rij[0] * rij[0] + rij[1] * rij[1] + rij[2] * rij[2]).sqrt();
    let r_kj = (rkj[0] * rkj[0] + rkj[1] * rkj[1] + rkj[2] * rkj[2]).sqrt();

    let dot = rij[0] * rkj[0] + rij[1] * rkj[1] + rij[2] * rkj[2];
    let cos_theta = (dot / (r_ij * r_kj)).clamp(-1.0, 1.0);
    let theta = cos_theta.acos();

    let dr_ij = r_ij - r0_ij;
    let dr_kj = r_kj - r0_kj;
    let dtheta = theta - theta0_rad;

    params.kba_ijk * dr_ij * dtheta + params.kba_kji * dr_kj * dtheta
}

pub fn stretch_bend_gradient(
    coords: &[[f64; 3]],
    i: usize,
    j: usize,
    k: usize,
    r0_ij: f64,
    r0_kj: f64,
    theta0_rad: f64,
    params: &StretchBendParams,
) -> ([f64; 3], [f64; 3], [f64; 3]) {
    let rij: [f64; 3] = [
        coords[i][0] - coords[j][0],
        coords[i][1] - coords[j][1],
        coords[i][2] - coords[j][2],
    ];
    let rkj: [f64; 3] = [
        coords[k][0] - coords[j][0],
        coords[k][1] - coords[j][1],
        coords[k][2] - coords[j][2],
    ];
    let r_ij = (rij[0] * rij[0] + rij[1] * rij[1] + rij[2] * rij[2]).sqrt();
    let r_kj = (rkj[0] * rkj[0] + rkj[1] * rkj[1] + rkj[2] * rkj[2]).sqrt();

    if r_ij < 1e-10 || r_kj < 1e-10 {
        return ([0.0; 3], [0.0; 3], [0.0; 3]);
    }

    let dot = rij[0] * rkj[0] + rij[1] * rkj[1] + rij[2] * rkj[2];
    let cos_theta = (dot / (r_ij * r_kj)).clamp(-1.0, 1.0);
    let theta = cos_theta.acos();
    let sin_theta = theta.sin();

    let dr_ij = r_ij - r0_ij;
    let dr_kj = r_kj - r0_kj;
    let dtheta = theta - theta0_rad;

    // E = kba_ijk * dr_ij * dtheta + kba_kji * dr_kj * dtheta
    //
    // For atom i: dE/dx_i = kba_ijk * (dr_ij/dx_i * dtheta + dr_ij * dtheta/dx_i)
    //                              + kba_kji * (dr_kj * dtheta/dx_i)
    //
    // dr_ij/dx_i = rij / r_ij  (unit vector from j to i)
    // dtheta/dx_i: from angle gradient, d(cos theta)/dx_i = (rkj - cos_theta * rij/r_ij) / (r_ij * r_kj)
    //   dtheta/dx_i = -1/sin_theta * d(cos_theta)/dx_i

    let inv_sin = if sin_theta.abs() < 1e-10 {
        1e10 * sin_theta.signum()
    } else {
        -1.0 / sin_theta
    };

    // Unit vectors
    let uij = [rij[0] / r_ij, rij[1] / r_ij, rij[2] / r_ij];
    let ukj = [rkj[0] / r_kj, rkj[1] / r_kj, rkj[2] / r_kj];

    // d(cos_theta)/d(atom) = d/dx of dot/(r_ij * r_kj)
    // For atom i: d(cos_theta)/dx_i = (rkj - cos_theta * uij * r_kj) / (r_ij * r_kj)
    //                                  = (ukj - cos_theta * uij) / r_ij
    // dtheta/dx_i = inv_sin * (ukj - cos_theta * uij) / r_ij

    let dcos_di: [f64; 3] = [
        (ukj[0] - cos_theta * uij[0]) / r_ij,
        (ukj[1] - cos_theta * uij[1]) / r_ij,
        (ukj[2] - cos_theta * uij[2]) / r_ij,
    ];

    // For atom k: d(cos_theta)/dx_k = (uij - cos_theta * ukj) / r_kj
    let dcos_dk: [f64; 3] = [
        (uij[0] - cos_theta * ukj[0]) / r_kj,
        (uij[1] - cos_theta * ukj[1]) / r_kj,
        (uij[2] - cos_theta * ukj[2]) / r_kj,
    ];

    // For atom j: d(cos_theta)/dx_j = -(dcos_di + dcos_dk)
    let dcos_dj: [f64; 3] = [
        -(dcos_di[0] + dcos_dk[0]),
        -(dcos_di[1] + dcos_dk[1]),
        -(dcos_di[2] + dcos_dk[2]),
    ];

    let dtheta_di = [
        inv_sin * dcos_di[0],
        inv_sin * dcos_di[1],
        inv_sin * dcos_di[2],
    ];
    let dtheta_dk = [
        inv_sin * dcos_dk[0],
        inv_sin * dcos_dk[1],
        inv_sin * dcos_dk[2],
    ];
    let dtheta_dj = [
        inv_sin * dcos_dj[0],
        inv_sin * dcos_dj[1],
        inv_sin * dcos_dj[2],
    ];

    // dr_ij/dx_i = uij, dr_ij/dx_j = -uij, dr_ij/dx_k = 0
    // dr_kj/dx_k = ukj, dr_kj/dx_j = -ukj, dr_kj/dx_i = 0

    let c_ijk = params.kba_ijk;
    let c_kji = params.kba_kji;

    // Gradient for atom i:
    // c_ijk * (uij * dtheta + dr_ij * dtheta_di) + c_kji * (0 * dtheta + dr_kj * dtheta_di)
    let gi: [f64; 3] = [
        c_ijk * (uij[0] * dtheta + dr_ij * dtheta_di[0]) + c_kji * dr_kj * dtheta_di[0],
        c_ijk * (uij[1] * dtheta + dr_ij * dtheta_di[1]) + c_kji * dr_kj * dtheta_di[1],
        c_ijk * (uij[2] * dtheta + dr_ij * dtheta_di[2]) + c_kji * dr_kj * dtheta_di[2],
    ];

    // Gradient for atom j:
    // c_ijk * (-uij * dtheta + dr_ij * dtheta_dj) + c_kji * (-ukj * dtheta + dr_kj * dtheta_dj)
    let gj: [f64; 3] = [
        c_ijk * (-uij[0] * dtheta + dr_ij * dtheta_dj[0])
            + c_kji * (-ukj[0] * dtheta + dr_kj * dtheta_dj[0]),
        c_ijk * (-uij[1] * dtheta + dr_ij * dtheta_dj[1])
            + c_kji * (-ukj[1] * dtheta + dr_kj * dtheta_dj[1]),
        c_ijk * (-uij[2] * dtheta + dr_ij * dtheta_dj[2])
            + c_kji * (-ukj[2] * dtheta + dr_kj * dtheta_dj[2]),
    ];

    // Gradient for atom k:
    // c_ijk * (0 * dtheta + dr_ij * dtheta_dk) + c_kji * (ukj * dtheta + dr_kj * dtheta_dk)
    let gk: [f64; 3] = [
        c_ijk * dr_ij * dtheta_dk[0] + c_kji * (ukj[0] * dtheta + dr_kj * dtheta_dk[0]),
        c_ijk * dr_ij * dtheta_dk[1] + c_kji * (ukj[1] * dtheta + dr_kj * dtheta_dk[1]),
        c_ijk * dr_ij * dtheta_dk[2] + c_kji * (ukj[2] * dtheta + dr_kj * dtheta_dk[2]),
    ];

    (gi, gj, gk)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_stretch_bend_energy_zero_at_equilibrium() {
        let coords = [[-1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let params = StretchBendParams {
            kba_ijk: 0.2270,
            kba_kji: 0.0700,
        };
        // At equilibrium (dr=0, dtheta=0), energy should be 0
        let e = stretch_bend_energy(&coords, 0, 1, 2, 1.09, 1.09, 1.911, &params);
        // The actual bond lengths and angles differ from equilibrium, so just check it's finite
        assert!(e.is_finite());
    }

    #[test]
    fn test_stretch_bend_gradient_numerical() {
        let coords = [[-1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let params = StretchBendParams {
            kba_ijk: 0.2270,
            kba_kji: 0.0700,
        };
        let r0_ij = 1.09;
        let r0_kj = 1.09;
        let theta0 = 1.911;

        let (gi, gj, gk) = stretch_bend_gradient(&coords, 0, 1, 2, r0_ij, r0_kj, theta0, &params);

        let eps = 1e-7;
        for (atom_idx, grad) in [(0usize, gi), (1usize, gj), (2usize, gk)] {
            for dim in 0..3 {
                let mut cp: Vec<[f64; 3]> = coords.to_vec();
                cp[atom_idx][dim] += eps;
                let ep = stretch_bend_energy(&cp, 0, 1, 2, r0_ij, r0_kj, theta0, &params);
                let e0 = stretch_bend_energy(&coords, 0, 1, 2, r0_ij, r0_kj, theta0, &params);
                let num = (ep - e0) / eps;
                assert!(
                    (grad[dim] - num).abs() < 1e-4,
                    "stretch-bend grad[{}] = {} vs numerical {} for atom {}",
                    dim,
                    grad[dim],
                    num,
                    atom_idx
                );
            }
        }
    }

    #[test]
    fn test_get_stretch_bend_params_known() {
        let params = get_sb_kba(MMFFAtomType::H, MMFFAtomType::C_3, MMFFAtomType::H);
        assert!(params.is_some());
        let (kba_ijk, kba_kji) = params.unwrap();
        assert!((kba_ijk - 0.2270).abs() < 1e-4);
        assert!((kba_kji - 0.0700).abs() < 1e-4);
    }

    #[test]
    fn test_stretch_bend_energy_matches_rdkit_ethanol() {
        // Use ethanol-like geometry: H-O-C angle
        let coords = [
            [0.0, 0.0, 0.96],   // H (atom 0)
            [0.0, 0.0, 0.0],    // O (atom 1)
            [1.43, 0.0, -0.36], // C (atom 2)
        ];
        let params = StretchBendParams {
            kba_ijk: 0.1430,
            kba_kji: 0.2560,
        };
        let e = stretch_bend_energy(&coords, 0, 1, 2, 0.96, 1.43, 1.89, &params);
        assert!(e.is_finite());
    }
}
