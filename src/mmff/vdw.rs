use super::params::mmff_type_id;
use super::MMFFAtomType;

const VDW_POWER: f64 = 0.25;
const VDW_B: f64 = 0.2;
const VDW_BETA: f64 = 12.0;
const VDW_DARAD: f64 = 0.8;
const VDW_DAEPS: f64 = 0.5;

#[derive(Debug, Clone, Copy)]
struct MMFFVdW {
    r_star: f64,
    alpha_i: f64,
    n_i: f64,
    g_i: f64,
    da: u8,
}

#[rustfmt::skip]
const VDW_TABLE: [(u8, f64, f64, f64, f64, u8); 95] = [
    // (type_id, alpha_i, N_i, A_i, G_i, DA)  DA: 0='-', 1='A', 2='D'
    ( 1, 1.050, 2.490, 3.890, 1.282, 0),
    ( 2, 1.350, 2.490, 3.890, 1.282, 0),
    ( 3, 1.100, 2.490, 3.890, 1.282, 0),
    ( 4, 1.300, 2.490, 3.890, 1.282, 0),
    ( 5, 0.250, 0.800, 4.200, 1.209, 0),
    ( 6, 0.700, 3.150, 3.890, 1.282, 1),
    ( 7, 0.650, 3.150, 3.890, 1.282, 1),
    ( 8, 1.150, 2.820, 3.890, 1.282, 1),
    ( 9, 0.900, 2.820, 3.890, 1.282, 1),
    (10, 1.000, 2.820, 3.890, 1.282, 1),
    (11, 0.350, 3.480, 3.890, 1.282, 1),
    (12, 2.300, 5.100, 3.320, 1.345, 1),
    (13, 3.400, 6.000, 3.190, 1.359, 1),
    (14, 5.500, 6.950, 3.080, 1.404, 1),
    (15, 3.000, 4.800, 3.320, 1.345, 1),
    (16, 3.900, 4.800, 3.320, 1.345, 1),
    (17, 2.700, 4.800, 3.320, 1.345, 0),
    (18, 2.100, 4.800, 3.320, 1.345, 0),
    (19, 4.500, 4.200, 3.320, 1.345, 0),
    (20, 1.050, 2.490, 3.890, 1.282, 0),
    (21, 0.150, 0.800, 4.200, 1.209, 2),
    (22, 1.100, 2.490, 3.890, 1.282, 0),
    (23, 0.150, 0.800, 4.200, 1.209, 2),
    (24, 0.150, 0.800, 4.200, 1.209, 2),
    (25, 1.600, 4.500, 3.320, 1.345, 0),
    (26, 3.600, 4.500, 3.320, 1.345, 1),
    (27, 0.150, 0.800, 4.200, 1.209, 2),
    (28, 0.150, 0.800, 4.200, 1.209, 2),
    (29, 0.150, 0.800, 4.200, 1.209, 2),
    (30, 1.350, 2.490, 3.890, 1.282, 0),
    (31, 0.150, 0.800, 4.200, 1.209, 2),
    (32, 0.750, 3.150, 3.890, 1.282, 1),
    (33, 0.150, 0.800, 4.200, 1.209, 2),
    (34, 1.000, 2.820, 3.890, 1.282, 0),
    (35, 1.500, 3.150, 3.890, 1.282, 1),
    (36, 0.150, 0.800, 4.200, 1.209, 2),
    (37, 1.350, 2.490, 3.890, 1.282, 0),
    (38, 0.850, 2.820, 3.890, 1.282, 1),
    (39, 1.100, 2.820, 3.890, 1.282, 0),
    (40, 1.000, 2.820, 3.890, 1.282, 1),
    (41, 1.100, 2.490, 3.890, 1.282, 0),
    (42, 1.000, 2.820, 3.890, 1.282, 1),
    (43, 1.000, 2.820, 3.890, 1.282, 1),
    (44, 3.000, 4.800, 3.320, 1.345, 1),
    (45, 1.150, 2.820, 3.890, 1.282, 0),
    (46, 1.300, 2.820, 3.890, 1.282, 0),
    (47, 1.000, 2.820, 3.890, 1.282, 1),
    (48, 1.200, 2.820, 3.890, 1.282, 1),
    (49, 1.000, 3.150, 3.890, 1.282, 0),
    (50, 0.150, 0.800, 4.200, 1.209, 2),
    (51, 0.400, 3.150, 3.890, 1.282, 0),
    (52, 0.150, 0.800, 4.200, 1.209, 2),
    (53, 1.000, 2.820, 3.890, 1.282, 0),
    (54, 1.300, 2.820, 3.890, 1.282, 0),
    (55, 0.800, 2.820, 3.890, 1.282, 0),
    (56, 0.800, 2.820, 3.890, 1.282, 0),
    (57, 1.000, 2.490, 3.890, 1.282, 0),
    (58, 0.800, 2.820, 3.890, 1.282, 0),
    (59, 0.650, 3.150, 3.890, 1.282, 1),
    (60, 1.800, 2.490, 3.890, 1.282, 1),
    (61, 0.800, 2.820, 3.890, 1.282, 1),
    (62, 1.300, 2.820, 3.890, 1.282, 1),
    (63, 1.350, 2.490, 3.890, 1.282, 0),
    (64, 1.350, 2.490, 3.890, 1.282, 0),
    (65, 1.000, 2.820, 3.890, 1.282, 1),
    (66, 0.750, 2.820, 3.890, 1.282, 1),
    (67, 0.950, 2.820, 3.890, 1.282, 1),
    (68, 0.900, 2.820, 3.890, 1.282, 1),
    (69, 0.950, 2.820, 3.890, 1.282, 1),
    (70, 0.870, 3.150, 3.890, 1.282, 1),
    (71, 0.150, 0.800, 4.200, 1.209, 2),
    (72, 4.000, 4.800, 3.320, 1.345, 1),
    (73, 3.000, 4.800, 3.320, 1.345, 0),
    (74, 3.000, 4.800, 3.320, 1.345, 0),
    (75, 4.000, 4.500, 3.320, 1.345, 1),
    (76, 1.200, 2.820, 3.890, 1.282, 1),
    (77, 1.500, 5.100, 3.320, 1.345, 1),
    (78, 1.350, 2.490, 3.890, 1.282, 0),
    (79, 1.000, 2.820, 3.890, 1.282, 1),
    (80, 1.000, 2.490, 3.890, 1.282, 0),
    (81, 0.800, 2.820, 3.890, 1.282, 0),
    (82, 0.950, 2.820, 3.890, 1.282, 1),
    (87, 0.450, 6.000, 4.000, 1.400, 0),
    (88, 0.550, 6.000, 4.000, 1.400, 0),
    (89, 1.400, 3.480, 3.890, 1.282, 1),
    (90, 4.500, 5.100, 3.320, 1.345, 1),
    (91, 6.000, 6.000, 3.190, 1.359, 1),
    (92, 0.150, 2.000, 4.000, 1.300, 0),
    (93, 0.400, 3.500, 4.000, 1.300, 0),
    (94, 1.000, 5.000, 4.000, 1.300, 0),
    (95, 0.430, 6.000, 4.000, 1.400, 0),
    (96, 0.900, 5.000, 4.000, 1.400, 0),
    (97, 0.350, 6.000, 4.000, 1.400, 0),
    (98, 0.400, 6.000, 4.000, 1.400, 0),
    (99, 0.350, 3.500, 4.000, 1.300, 0),
];

fn lookup_vdw(type_id: u8) -> MMFFVdW {
    VDW_TABLE
        .binary_search_by_key(&type_id, |&(tid, _, _, _, _, _)| tid)
        .map(|idx| {
            let (_, alpha_i, n_i, a_i, g_i, da) = VDW_TABLE[idx];
            MMFFVdW {
                r_star: a_i * alpha_i.powf(VDW_POWER),
                alpha_i,
                n_i,
                g_i,
                da,
            }
        })
        .unwrap_or(MMFFVdW {
            r_star: 3.890,
            alpha_i: 1.0,
            n_i: 5.0,
            g_i: 1.282,
            da: 0,
        })
}

pub struct VDWParams {
    pub r_star: f64,
    pub alpha_i: f64,
    pub n_i: f64,
    pub g_i: f64,
    pub da: u8,
}

pub fn get_vdw_params(atom_type: MMFFAtomType) -> VDWParams {
    let type_id = mmff_type_id(atom_type);
    let v = lookup_vdw(type_id);
    VDWParams {
        r_star: v.r_star,
        alpha_i: v.alpha_i,
        n_i: v.n_i,
        g_i: v.g_i,
        da: v.da,
    }
}

fn calc_r_star_ij(pi: &VDWParams, pj: &VDWParams) -> f64 {
    let gamma = if (pi.r_star + pj.r_star) > 1e-10 {
        (pi.r_star - pj.r_star) / (pi.r_star + pj.r_star)
    } else {
        0.0
    };
    let b_factor = if pi.da == 2 || pj.da == 2 {
        0.0
    } else {
        VDW_B * (1.0 - (-VDW_BETA * gamma * gamma).exp())
    };
    0.5 * (pi.r_star + pj.r_star) * (1.0 + b_factor)
}

fn calc_well_depth(r_star_ij: f64, pi: &VDWParams, pj: &VDWParams) -> f64 {
    let c4 = 181.16;
    let r2 = r_star_ij * r_star_ij;
    let r6 = r2 * r2 * r2;
    c4 * pi.g_i * pj.g_i * pi.alpha_i * pj.alpha_i
        / ((pi.alpha_i / pi.n_i).sqrt() + (pj.alpha_i / pj.n_i).sqrt())
        / r6
}

fn apply_da_scaling(r_star_ij: &mut f64, well_depth: &mut f64, pi: &VDWParams, pj: &VDWParams) {
    if (pi.da == 2 && pj.da == 1) || (pi.da == 1 && pj.da == 2) {
        *r_star_ij *= VDW_DARAD;
        *well_depth *= VDW_DAEPS;
    }
}

fn calc_vdw_energy(dist: f64, r_star_ij: f64, well_depth: f64) -> f64 {
    let vdw1 = 1.07;
    let vdw1m1 = 0.07;
    let vdw2 = 1.12;
    let vdw2m1 = 0.12;

    let dist2 = dist * dist;
    let dist7 = dist2 * dist2 * dist2 * dist;
    let r2 = r_star_ij * r_star_ij;
    let r7 = r2 * r2 * r2 * r_star_ij;

    let a_term = vdw1 * r_star_ij / (dist + vdw1m1 * r_star_ij);
    let a2 = a_term * a_term;
    let a7 = a2 * a2 * a2 * a_term;
    let b_term = vdw2 * r7 / (dist7 + vdw2m1 * r7) - 2.0;

    well_depth * a7 * b_term
}

pub fn vdw_energy_and_gradient(
    coords: &[[f64; 3]],
    i: usize,
    j: usize,
    params_i: &VDWParams,
    params_j: &VDWParams,
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

    let mut r_star_ij = calc_r_star_ij(params_i, params_j);
    let mut well_depth = calc_well_depth(r_star_ij, params_i, params_j);
    apply_da_scaling(&mut r_star_ij, &mut well_depth, params_i, params_j);

    let energy = calc_vdw_energy(r, r_star_ij, well_depth);

    let eps = 1e-7;
    let mut grad_i = [0.0; 3];
    let mut grad_j = [0.0; 3];
    for dim in 0..3 {
        let mut cp: Vec<[f64; 3]> = coords.to_vec();
        cp[i][dim] += eps;
        let e_plus = calc_vdw_energy_for(&cp, i, j, r_star_ij, well_depth);
        let e0 = calc_vdw_energy_for(coords, i, j, r_star_ij, well_depth);
        let num = (e_plus - e0) / eps;
        grad_i[dim] = num;
        grad_j[dim] = -num;
    }

    (energy, grad_i, grad_j)
}

fn calc_vdw_energy_for(
    coords: &[[f64; 3]],
    i: usize,
    j: usize,
    r_star_ij: f64,
    well_depth: f64,
) -> f64 {
    let dx = coords[j][0] - coords[i][0];
    let dy = coords[j][1] - coords[i][1];
    let dz = coords[j][2] - coords[i][2];
    let r = (dx * dx + dy * dy + dz * dz).sqrt();

    if r < 1e-10 {
        return 0.0;
    }

    calc_vdw_energy(r, r_star_ij, well_depth)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_params(type_id: u8) -> VDWParams {
        let v = lookup_vdw(type_id);
        VDWParams {
            r_star: v.r_star,
            alpha_i: v.alpha_i,
            n_i: v.n_i,
            g_i: v.g_i,
            da: v.da,
        }
    }

    #[test]
    fn test_vdw_energy_has_minimum() {
        let pi = make_params(1); // C_3
        let pj = make_params(1);
        let coords_close = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let coords_eq = [[0.0, 0.0, 0.0], [pi.r_star, 0.0, 0.0]];
        let coords_far = [[0.0, 0.0, 0.0], [3.0 * pi.r_star, 0.0, 0.0]];

        let (e_close, _, _) = vdw_energy_and_gradient(&coords_close, 0, 1, &pi, &pj, false);
        let (e_eq, _, _) = vdw_energy_and_gradient(&coords_eq, 0, 1, &pi, &pj, false);
        let (e_far, _, _) = vdw_energy_and_gradient(&coords_far, 0, 1, &pi, &pj, false);

        assert!(e_close > 0.0, "VDW repulsive at short range: {}", e_close);
        assert!(e_eq < 0.0, "VDW attractive well near r_star: {}", e_eq);
        assert!(
            e_far > e_eq,
            "VDW rises toward zero at long range: e_far={} > e_eq={}",
            e_far,
            e_eq
        );
    }

    #[test]
    fn test_vdw_no_14_scaling() {
        let pi = make_params(1);
        let pj = make_params(1);
        let coords = [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]];

        let (e_full, _, _) = vdw_energy_and_gradient(&coords, 0, 1, &pi, &pj, false);
        let (e_14, _, _) = vdw_energy_and_gradient(&coords, 0, 1, &pi, &pj, true);

        assert!(
            (e_14 - e_full).abs() < 1e-10,
            "VDW should NOT be scaled for 1-4: got {}, expected {}",
            e_14,
            e_full
        );
    }

    #[test]
    fn test_slater_kirkwood_well_depth() {
        let pi = make_params(1); // C_3: type_id=1
        let pj = make_params(5); // H: type_id=5
        let r_star_ij = calc_r_star_ij(&pi, &pj);
        let eps = calc_well_depth(r_star_ij, &pi, &pj);
        assert!(eps > 0.0, "Well depth should be positive: {}", eps);
        assert!(
            eps < 1.0,
            "Well depth should be small (kcal/mol range): {}",
            eps
        );
    }

    #[test]
    fn test_vdw_gradient_numerical() {
        let pi = make_params(1);
        let pj = make_params(7);
        let coords = [[0.0, 0.0, 0.0], [1.8, 0.0, 0.0]];

        let (_, gi, gj) = vdw_energy_and_gradient(&coords, 0, 1, &pi, &pj, false);

        let eps = 1e-7;
        for (atom_idx, grad) in [(0usize, gi), (1usize, gj)] {
            for dim in 0..3 {
                let mut cp2: Vec<[f64; 3]> = coords.to_vec();
                cp2[atom_idx][dim] += eps;
                let (ep, _, _) = vdw_energy_and_gradient(&cp2, 0, 1, &pi, &pj, false);
                let (e0, _, _) = vdw_energy_and_gradient(&coords, 0, 1, &pi, &pj, false);
                let num = (ep - e0) / eps;
                assert!(
                    (grad[dim] - num).abs() < 1e-4,
                    "VDW grad[{}] = {} vs numerical {} for atom {}",
                    dim,
                    grad[dim],
                    num,
                    atom_idx
                );
            }
        }
    }

    #[test]
    fn test_da_scaling() {
        let pi = make_params(6); // O_3 (type_id=6, DA=A)
        let pj = make_params(21); // H_ONC (type_id=21, DA=D)
        let mut r_star_ij = calc_r_star_ij(&pi, &pj);
        let mut well_depth = calc_well_depth(r_star_ij, &pi, &pj);
        let r_before = r_star_ij;
        let w_before = well_depth;
        apply_da_scaling(&mut r_star_ij, &mut well_depth, &pi, &pj);
        assert!(
            (r_star_ij - r_before * VDW_DARAD).abs() < 1e-10,
            "DARAD scaling"
        );
        assert!(
            (well_depth - w_before * VDW_DAEPS).abs() < 1e-10,
            "DAEPS scaling"
        );
    }

    #[test]
    fn test_no_da_scaling_non_da_pair() {
        let pi = make_params(1); // C_3 (DA=0)
        let pj = make_params(5); // H (DA=0)
        let mut r_star_ij = calc_r_star_ij(&pi, &pj);
        let mut well_depth = calc_well_depth(r_star_ij, &pi, &pj);
        let r_before = r_star_ij;
        let w_before = well_depth;
        apply_da_scaling(&mut r_star_ij, &mut well_depth, &pi, &pj);
        assert!((r_star_ij - r_before).abs() < 1e-10, "No DARAD scaling");
        assert!((well_depth - w_before).abs() < 1e-10, "No DAEPS scaling");
    }
}
