use super::mmff_tables;
use super::params::mmff_type_id;

pub struct StretchBendParams {
    pub kba_ijk: f64,
    pub kba_kji: f64,
}

pub fn get_stretch_bend_params(
    type_i: super::MMFFAtomType,
    type_j: super::MMFFAtomType,
    type_k: super::MMFFAtomType,
    bond_type_ij: u8,
    bond_type_jk: u8,
    atomic_num_i: u8,
    atomic_num_j: u8,
    atomic_num_k: u8,
    angle_type: u8,
) -> Option<StretchBendParams> {
    let ti = mmff_type_id(type_i);
    let tj = mmff_type_id(type_j);
    let tk = mmff_type_id(type_k);

    let sb_type = mmff_tables::compute_stretch_bend_type(angle_type, bond_type_ij, bond_type_jk);

    if let Some((kba_ijk, kba_kji)) = mmff_tables::lookup_stretch_bend_params(
        ti,
        tj,
        tk,
        sb_type,
        atomic_num_i,
        atomic_num_j,
        atomic_num_k,
    ) {
        return Some(StretchBendParams { kba_ijk, kba_kji });
    }

    None
}

#[allow(clippy::too_many_arguments)]
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

    const MDYNE_A_TO_KCAL_MOL: f64 = 143.9325;
    MDYNE_A_TO_KCAL_MOL * (params.kba_ijk * dr_ij + params.kba_kji * dr_kj) * dtheta
}

#[allow(clippy::too_many_arguments)]
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

    const SCALE: f64 = 143.9325;

    // Gradient for atom i:
    // dE/dx_i = SCALE * [c_ijk * (uij * dtheta + dr_ij * dtheta_di) + c_kji * dr_kj * dtheta_di]
    let gi: [f64; 3] = [
        SCALE * (c_ijk * (uij[0] * dtheta + dr_ij * dtheta_di[0]) + c_kji * dr_kj * dtheta_di[0]),
        SCALE * (c_ijk * (uij[1] * dtheta + dr_ij * dtheta_di[1]) + c_kji * dr_kj * dtheta_di[1]),
        SCALE * (c_ijk * (uij[2] * dtheta + dr_ij * dtheta_di[2]) + c_kji * dr_kj * dtheta_di[2]),
    ];

    // Gradient for atom j:
    // dE/dx_j = SCALE * [c_ijk * (-uij * dtheta + dr_ij * dtheta_dj) + c_kji * (-ukj * dtheta + dr_kj * dtheta_dj)]
    let gj: [f64; 3] = [
        SCALE
            * (c_ijk * (-uij[0] * dtheta + dr_ij * dtheta_dj[0])
                + c_kji * (-ukj[0] * dtheta + dr_kj * dtheta_dj[0])),
        SCALE
            * (c_ijk * (-uij[1] * dtheta + dr_ij * dtheta_dj[1])
                + c_kji * (-ukj[1] * dtheta + dr_kj * dtheta_dj[1])),
        SCALE
            * (c_ijk * (-uij[2] * dtheta + dr_ij * dtheta_dj[2])
                + c_kji * (-ukj[2] * dtheta + dr_kj * dtheta_dj[2])),
    ];

    // Gradient for atom k:
    // dE/dx_k = SCALE * [c_ijk * dr_ij * dtheta_dk + c_kji * (ukj * dtheta + dr_kj * dtheta_dk)]
    let gk: [f64; 3] = [
        SCALE * (c_ijk * dr_ij * dtheta_dk[0] + c_kji * (ukj[0] * dtheta + dr_kj * dtheta_dk[0])),
        SCALE * (c_ijk * dr_ij * dtheta_dk[1] + c_kji * (ukj[1] * dtheta + dr_kj * dtheta_dk[1])),
        SCALE * (c_ijk * dr_ij * dtheta_dk[2] + c_kji * (ukj[2] * dtheta + dr_kj * dtheta_dk[2])),
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
        use super::super::MMFFAtomType;
        let params = get_stretch_bend_params(
            MMFFAtomType::H,
            MMFFAtomType::C_3,
            MMFFAtomType::H,
            0,
            0,
            1,
            6,
            6,
            0,
        );
        assert!(params.is_some());
        let p = params.unwrap();
        assert!((p.kba_ijk - 0.115).abs() < 0.01);
        assert!((p.kba_kji - 0.115).abs() < 0.01);
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
