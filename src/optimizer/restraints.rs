//! Geometric restraints for the Cartesian optimizer (v1.7.0).
//!
//! Adds flat-bottom harmonic terms (distance / angle / dihedral) and frozen
//! atoms on top of any ForceField. Semantics: within ±tol of the target the
//! term is exactly zero; outside, E = ½k·(x − x0 ∓ tol)². Dihedral targets
//! wrap to (−180°, 180°] before the tolerance test so restraints near ±180°
//! behave correctly. Gradients are analytic — the dihedral Jacobian is the
//! same ETKDG kernel metadynamics uses (`etkdg::dihedral_gradient_contrib`),
//! the distance term the same unit-vector chain rule as `DistanceCV`; the
//! angle gradient is the standard cos-form chain rule with a collinearity
//! guard. A restraint with k = 0 contributes nothing and is skipped.

use std::f64::consts::PI;

/// One geometric restraint. Units: lengths Å, angles degrees (converted
/// internally), k kcal/mol/Å² (distance) or kcal/mol/rad² (angle, dihedral),
/// tol in the natural unit of the coordinate (Å / degrees).
#[derive(Debug, Clone, Copy)]
pub enum Restraint {
    /// Restrain |x_i − x_j| to r0.
    Distance {
        i: usize,
        j: usize,
        r0: f64,
        fc: f64,
        tol: f64,
    },
    /// Restrain angle (i, vertex j, k) to a0 degrees.
    Angle {
        i: usize,
        j: usize,
        k: usize,
        a0: f64,
        fc: f64,
        tol: f64,
    },
    /// Restrain dihedral (i, j, k, l) to a0 degrees.
    Dihedral {
        i: usize,
        j: usize,
        k: usize,
        l: usize,
        a0: f64,
        fc: f64,
        tol: f64,
    },
    /// Harmonically pull atom i toward a fixed point in space (POSRES —
    /// the flexible-alignment primitive: map probe atoms onto reference
    /// coordinates).
    Position {
        i: usize,
        target: [f64; 3],
        fc: f64,
        tol: f64,
    },
}

/// A restraint set plus frozen atom indices. Frozen atoms are implemented by
/// zeroing their gradient rows in the optimizer objective — they never move
/// because the search direction stays zero in those components.
#[derive(Debug, Clone, Default)]
pub struct RestraintSet {
    pub restraints: Vec<Restraint>,
    pub frozen: Vec<usize>,
}

impl RestraintSet {
    pub fn is_empty(&self) -> bool {
        self.restraints.is_empty() && self.frozen.is_empty()
    }

    /// Every atom index referenced (restraint atoms + frozen), for bounds
    /// validation against the molecule before optimizing.
    pub fn atom_indices(&self) -> Vec<usize> {
        let mut v: Vec<usize> = self.frozen.clone();
        for r in &self.restraints {
            match *r {
                Restraint::Distance { i, j, .. } => v.extend_from_slice(&[i, j]),
                Restraint::Angle { i, j, k, .. } => v.extend_from_slice(&[i, j, k]),
                Restraint::Dihedral { i, j, k, l, .. } => v.extend_from_slice(&[i, j, k, l]),
                Restraint::Position { i, .. } => v.push(i),
            }
        }
        v
    }

    /// Restraint energy (kcal/mol) and gradient at `coords`. The same
    /// function backs both the objective's `f_and_g` and `energy` paths —
    /// line-search trials cannot drift from accepted-point energies.
    pub fn energy_and_gradient(&self, coords: &[[f64; 3]]) -> (f64, Vec<[f64; 3]>) {
        let n = coords.len();
        let mut e = 0.0f64;
        let mut g = vec![[0.0f64; 3]; n];
        for r in &self.restraints {
            match *r {
                Restraint::Distance { i, j, r0, fc, tol } => {
                    if fc == 0.0 {
                        continue;
                    }
                    let dx = coords[i][0] - coords[j][0];
                    let dy = coords[i][1] - coords[j][1];
                    let dz = coords[i][2] - coords[j][2];
                    let d = (dx * dx + dy * dy + dz * dz).sqrt();
                    if d < 1e-10 {
                        continue; // coincident atoms: distance gradient undefined
                    }
                    let dev = d - r0;
                    if dev.abs() <= tol {
                        continue; // flat bottom: exactly zero inside
                    }
                    let x = dev - dev.signum() * tol;
                    e += 0.5 * fc * x * x;
                    // dE/dx_i = fc·x·(x_i−x_j)/d; opposite on j
                    let f = fc * x / d;
                    let gi = [f * dx, f * dy, f * dz];
                    g[i][0] += gi[0];
                    g[i][1] += gi[1];
                    g[i][2] += gi[2];
                    g[j][0] -= gi[0];
                    g[j][1] -= gi[1];
                    g[j][2] -= gi[2];
                }
                Restraint::Angle {
                    i,
                    j,
                    k,
                    a0,
                    fc,
                    tol,
                } => {
                    if fc == 0.0 {
                        continue;
                    }
                    let u = [
                        coords[i][0] - coords[j][0],
                        coords[i][1] - coords[j][1],
                        coords[i][2] - coords[j][2],
                    ];
                    let v = [
                        coords[k][0] - coords[j][0],
                        coords[k][1] - coords[j][1],
                        coords[k][2] - coords[j][2],
                    ];
                    let un = (u[0] * u[0] + u[1] * u[1] + u[2] * u[2]).sqrt();
                    let vn = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
                    if un < 1e-10 || vn < 1e-10 {
                        continue;
                    }
                    let cos_t = (u[0] * v[0] + u[1] * v[1] + u[2] * v[2]) / (un * vn);
                    let theta = cos_t.clamp(-1.0, 1.0).acos();
                    let dev = theta - a0 * PI / 180.0;
                    if dev.abs() <= tol * PI / 180.0 {
                        continue;
                    }
                    let sin_t = theta.sin();
                    if sin_t.abs() < 1e-8 {
                        continue; // collinear: dtheta/dx singular — skip (documented)
                    }
                    let x = dev - dev.signum() * tol * PI / 180.0;
                    e += 0.5 * fc * x * x;
                    // dE/dx = fc·x·dtheta/dx; dtheta = −dcos/sin
                    let c = -fc * x / sin_t;
                    // dcos/du = v/(|u||v|) − (u·v)u/(|u|³|v|); dcos/dv symmetric
                    let udv = u[0] * v[0] + u[1] * v[1] + u[2] * v[2];
                    let du = [
                        v[0] / (un * vn) - udv * u[0] / (un * un * un * vn),
                        v[1] / (un * vn) - udv * u[1] / (un * un * un * vn),
                        v[2] / (un * vn) - udv * u[2] / (un * un * un * vn),
                    ];
                    let dv = [
                        u[0] / (un * vn) - udv * v[0] / (un * vn * vn * vn),
                        u[1] / (un * vn) - udv * v[1] / (un * vn * vn * vn),
                        u[2] / (un * vn) - udv * v[2] / (un * vn * vn * vn),
                    ];
                    // x_i moves u (+), x_j moves both (−), x_k moves v (+)
                    g[i][0] += c * du[0];
                    g[i][1] += c * du[1];
                    g[i][2] += c * du[2];
                    g[k][0] += c * dv[0];
                    g[k][1] += c * dv[1];
                    g[k][2] += c * dv[2];
                    g[j][0] -= c * (du[0] + dv[0]);
                    g[j][1] -= c * (du[1] + dv[1]);
                    g[j][2] -= c * (du[2] + dv[2]);
                }
                Restraint::Position { i, target, fc, tol } => {
                    if fc == 0.0 {
                        continue;
                    }
                    let d = [
                        coords[i][0] - target[0],
                        coords[i][1] - target[1],
                        coords[i][2] - target[2],
                    ];
                    let r = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
                    let dev = r - tol;
                    if dev <= 0.0 {
                        continue; // inside the flat bottom
                    }
                    e += 0.5 * fc * dev * dev;
                    // dE/dx_i = fc·dev·(x_i−target)/r
                    let f = fc * dev / r;
                    g[i][0] += f * d[0];
                    g[i][1] += f * d[1];
                    g[i][2] += f * d[2];
                }
                Restraint::Dihedral {
                    i,
                    j,
                    k,
                    l,
                    a0,
                    fc,
                    tol,
                } => {
                    if fc == 0.0 {
                        continue;
                    }
                    let phi =
                        crate::etkdg::dihedral_angle4(coords[i], coords[j], coords[k], coords[l]);
                    // wrap the deviation to (−180°, 180°] — targets near ±180°
                    let mut dev = phi - a0 * PI / 180.0;
                    while dev > PI {
                        dev -= 2.0 * PI;
                    }
                    while dev <= -PI {
                        dev += 2.0 * PI;
                    }
                    if dev.abs() <= tol * PI / 180.0 {
                        continue;
                    }
                    let x = dev - dev.signum() * tol * PI / 180.0;
                    e += 0.5 * fc * x * x;
                    // dE/dx = fc·x · dphi/dx via the shared analytic Jacobian
                    let (g0, g1, g2, g3) =
                        crate::etkdg::dihedral_gradient_contrib(coords, i, j, k, l, fc * x);
                    for a in 0..3 {
                        g[i][a] += g0[a];
                        g[j][a] += g1[a];
                        g[k][a] += g2[a];
                        g[l][a] += g3[a];
                    }
                }
            }
        }
        (e, g)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn fd_check(rs: &RestraintSet, coords: &mut [[f64; 3]], label: &str) {
        let (_, g) = rs.energy_and_gradient(coords);
        let h = 1e-6;
        let mut maxerr = 0.0f64;
        for a in 0..coords.len() {
            for c in 0..3 {
                let old = coords[a][c];
                coords[a][c] = old + h;
                let ep = rs.energy_and_gradient(coords).0;
                coords[a][c] = old - h;
                let em = rs.energy_and_gradient(coords).0;
                coords[a][c] = old;
                let fd = (ep - em) / (2.0 * h);
                maxerr = maxerr.max((fd - g[a][c]).abs());
            }
        }
        assert!(maxerr < 1e-5, "{label}: FD err {maxerr}");
    }

    #[test]
    fn distance_gradient_fd_inside_and_outside_tol() {
        // outside the flat bottom (r = 2.9 vs r0 = 3.5, tol = 0.1)
        let rs = RestraintSet {
            restraints: vec![Restraint::Distance {
                i: 0,
                j: 1,
                r0: 3.5,
                fc: 10.0,
                tol: 0.1,
            }],
            frozen: vec![],
        };
        let mut x = [[0.0, 0.0, 0.0], [2.9, 0.3, -0.4], [5.0, 5.0, 5.0]];
        fd_check(&rs, &mut x, "distance-outside");
        // and on the other side (r = 4.2 > r0 + tol)
        let mut x2 = [[0.0, 0.0, 0.0], [4.0, 1.3, 0.2], [5.0, 5.0, 5.0]];
        fd_check(&rs, &mut x2, "distance-far-side");
        // inside the flat bottom: energy AND gradient exactly zero
        let x3 = [[0.0, 0.0, 0.0], [3.5, 0.05, 0.0], [9.0, 9.0, 9.0]];
        let (e, g) = rs.energy_and_gradient(&x3);
        assert_eq!(e, 0.0);
        assert!(g.iter().flatten().all(|&v| v == 0.0));
    }

    #[test]
    fn angle_gradient_fd() {
        let rs = RestraintSet {
            restraints: vec![Restraint::Angle {
                i: 0,
                j: 1,
                k: 2,
                a0: 120.0,
                fc: 10.0,
                tol: 5.0,
            }],
            frozen: vec![],
        };
        let mut x = [[1.0, 0.2, 0.0], [0.0, 0.0, 0.0], [-0.6, 1.1, 0.3]];
        fd_check(&rs, &mut x, "angle-outside");
        // tolerance straddling: 109.5 vs 120 ± 5 → 115 outside by 5.5°, still fine
        let mut x2 = [[1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [-0.33, 0.94, 0.0]];
        fd_check(&rs, &mut x2, "angle-near-tol");
    }

    #[test]
    fn dihedral_gradient_fd_including_wrap() {
        let rs = RestraintSet {
            restraints: vec![Restraint::Dihedral {
                i: 0,
                j: 1,
                k: 2,
                l: 3,
                a0: 170.0,
                fc: 10.0,
                tol: 5.0,
            }],
            frozen: vec![],
        };
        // phi ≈ −168°: dev wraps to +22° (not −338°) — the exact case the
        // wrap exists for
        let mut x = [
            [1.2, 0.1, 0.2],
            [0.0, 0.0, 0.0],
            [-1.0, 0.3, 0.1],
            [-1.6, -0.9, -0.4],
        ];
        fd_check(&rs, &mut x, "dihedral-wrap");
        // plain mid-range case
        let rs2 = RestraintSet {
            restraints: vec![Restraint::Dihedral {
                i: 0,
                j: 1,
                k: 2,
                l: 3,
                a0: 60.0,
                fc: 8.0,
                tol: 0.0,
            }],
            frozen: vec![],
        };
        let mut x2 = [
            [1.1, 0.4, -0.2],
            [0.0, 0.0, 0.0],
            [-0.9, 0.4, 0.2],
            [-1.5, 1.0, -0.5],
        ];
        fd_check(&rs2, &mut x2, "dihedral-plain");
    }

    #[test]
    fn wrap_matches_shortest_arc() {
        let rs = RestraintSet {
            restraints: vec![Restraint::Dihedral {
                i: 0,
                j: 1,
                k: 2,
                l: 3,
                a0: 179.0,
                fc: 10.0,
                tol: 0.0,
            }],
            frozen: vec![],
        };
        // phi ≈ −178° must measure dev ≈ −3° (shortest arc through ±180),
        // i.e. E ≈ ½·10·(0.052)² ≈ 0.014 kcal/mol, NOT a 357° catastrophe
        let x = [
            [1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [-1.0, 0.0, 0.0],
            [-2.0, 0.10, -0.02],
        ];
        let (e, _) = rs.energy_and_gradient(&x);
        assert!(e < 0.2, "wrapped energy should be small, got {e}");
    }

    #[test]
    fn position_gradient_fd_and_flat_bottom() {
        let rs = RestraintSet {
            restraints: vec![Restraint::Position {
                i: 1,
                target: [1.5, -0.4, 0.8],
                fc: 25.0,
                tol: 0.5,
            }],
            frozen: vec![],
        };
        let mut x = [[0.0, 0.0, 0.0], [2.9, 0.3, -0.4], [5.0, 5.0, 5.0]];
        fd_check(&rs, &mut x, "position-outside");
        // inside the flat bottom (|x−target| < tol): exactly zero
        let x2 = [[0.0; 3], [1.6, -0.2, 0.7], [9.0; 3]];
        let (e, g) = rs.energy_and_gradient(&x2);
        assert_eq!(e, 0.0);
        assert!(g.iter().flatten().all(|&v| v == 0.0));
    }

    #[test]
    fn k_zero_skips_and_degenerate_guards_are_finite() {
        let rs = RestraintSet {
            restraints: vec![
                Restraint::Distance {
                    i: 0,
                    j: 1,
                    r0: 1.0,
                    fc: 0.0,
                    tol: 0.0,
                },
                Restraint::Angle {
                    i: 0,
                    j: 1,
                    k: 2,
                    a0: 90.0,
                    fc: 0.0,
                    tol: 0.0,
                },
                Restraint::Dihedral {
                    i: 0,
                    j: 1,
                    k: 2,
                    l: 3,
                    a0: 90.0,
                    fc: 0.0,
                    tol: 0.0,
                },
            ],
            frozen: vec![],
        };
        let x = [[0.0; 3], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0], [3.0, 0.0, 0.0]];
        let (e, g) = rs.energy_and_gradient(&x);
        assert_eq!(e, 0.0);
        assert!(g.iter().flatten().all(|&v| v == 0.0));
        // fully collinear angle with k>0: guarded, finite, no NaN
        let rs2 = RestraintSet {
            restraints: vec![Restraint::Angle {
                i: 0,
                j: 1,
                k: 2,
                a0: 90.0,
                fc: 10.0,
                tol: 0.0,
            }],
            frozen: vec![],
        };
        let (e2, g2) = rs2.energy_and_gradient(&x);
        assert!(e2.is_finite());
        assert!(g2.iter().flatten().all(|v| v.is_finite()));
    }
}
