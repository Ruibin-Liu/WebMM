//! Generalized Born (GB) implicit-solvent model — a composable `ForceField` term.
//!
//! Wraps a vacuum force field (e.g. MMFF) and adds the Onufriev-Bashford-Case
//! (OBC2) GB electrostatic solvation correction.
//!
//! ## Born radii (HCT + OBC2)
//! ψ_ij = exact Coulomb-field integral of 1/r⁴ over atom j's sphere:
//! `ψ_ij = (1/4d)·ln((d-Rj)/(d+Rj)) + Rj/(2(d²-Rj²))`  [d = r_ij, d > Rj]
//! Born radius: `1/R_i = (1/Ri*)·(1 - tanh(αψ_i - βψ_i² + γψ_i³))`  (OBC2).
//!
//! ## GB energy (additive solvation correction on top of vacuum Coulomb)
//! `ΔG_GB = -C·(1/ε_in - 1/ε_out)·[Σ_{i<j} q_i·q_j/f_ij + ½·Σ_i q_i²/R_i]`
//! `f_ij = sqrt(r_ij² + R_i·R_j·exp(-r_ij²/(4·R_i·R_j)))`,  C = 332.053
//!
//! ## Gradient (3 parts): direct pair d(1/f)/dr; born-radius self q_i²/2·d(1/R_i)/dx;
//! born-radius pair q_i·q_j·∂(1/f)/∂R_i·dR_i/dx (f depends on R_i!).
//! The latter two share the structure B_i·chain[i]·dψ_i/dx where
//! B_i = c·(½·q_i² + Σ_j q_i·q_j·∂(1/f_ij)/∂(1/R_i)).

use crate::forces::ForceField;
use crate::molecule::Molecule;

const KE: f64 = 332.053;
const OBC_ALPHA: f64 = 1.0;
const OBC_BETA: f64 = 0.8;
const OBC_GAMMA: f64 = 4.851_713_813_7;

fn bondi_radius(z: u8) -> f64 {
    match z {
        1 => 1.20,
        6 => 1.70,
        7 => 1.55,
        8 => 1.52,
        9 => 1.47,
        14 => 2.10,
        15 => 1.80,
        16 => 1.80,
        17 => 1.75,
        35 => 1.85,
        53 => 1.98,
        _ => 1.5,
    }
}

/// LCPO (Linear Combination of Pairwise Overlaps) SASA coefficients [a1, a2, a3]
/// from Weiser, Tsui, Case 1999. Used for the nonpolar SA term.
fn lcpo_coeffs(z: u8) -> [f64; 3] {
    match z {
        1 => [0.1120, 0.0048, -0.0022],  // H
        6 => [0.1179, 0.0022, -0.0013],  // C (sp3; sp2 similar)
        7 => [0.1560, 0.0111, -0.0048],  // N
        8 => [0.1960, 0.0100, -0.0046],  // O
        16 => [0.2350, 0.0100, -0.0046], // S
        _ => [0.1179, 0.0022, -0.0013],  // default: C
    }
}

#[derive(Clone)]
pub struct GBSAConfig {
    pub solvent_dielectric: f64,
    pub solute_dielectric: f64,
    pub radius_offset: f64,
    pub obc2: bool,
    pub sa_surface_tension: f64,
}

impl Default for GBSAConfig {
    fn default() -> Self {
        Self {
            solvent_dielectric: 78.5,
            solute_dielectric: 1.0,
            radius_offset: 0.195,
            obc2: true,
            sa_surface_tension: 0.0,
        }
    }
}

pub struct GBSA<'a> {
    ff: &'a dyn ForceField,
    charges: Vec<f64>,
    radii: Vec<f64>,     // R_vdw + offset (Born radius boundary)
    sa_radii: Vec<f64>,  // R_vdw + probe (SASA boundary)
    lcpo: Vec<[f64; 3]>, // (a1, a2, a3) LCPO coefficients per atom
    dielectric_factor: f64,
    obc2: bool,
    sa_tension: f64,
    n: usize,
}

impl<'a> GBSA<'a> {
    pub fn new(
        ff: &'a dyn ForceField,
        mol: &Molecule,
        charges: &[f64],
        config: &GBSAConfig,
    ) -> Self {
        let n = mol.atoms.len();
        let radii: Vec<f64> = mol
            .atoms
            .iter()
            .map(|a| bondi_radius(a.atomic_number) + config.radius_offset)
            .collect();
        let sa_radii: Vec<f64> = mol
            .atoms
            .iter()
            .map(|a| bondi_radius(a.atomic_number) + 1.4)
            .collect();
        let lcpo: Vec<[f64; 3]> = mol
            .atoms
            .iter()
            .map(|a| lcpo_coeffs(a.atomic_number))
            .collect();
        Self {
            ff,
            charges: charges.to_vec(),
            radii,
            sa_radii,
            lcpo,
            dielectric_factor: 1.0 / config.solute_dielectric - 1.0 / config.solvent_dielectric,
            obc2: config.obc2,
            sa_tension: config.sa_surface_tension,
            n,
        }
    }

    // ---- Pairwise desolvation: exact HCT integral over V_j outside V_i ----
    // Handles all overlap cases (non-overlap, partial, full containment).
    // Derived analytically: the full-sphere and spherical-cap antiderivatives.

    /// Antiderivative of the full-shell angular integral: ∫ s²/(d²-s²)² ds
    #[inline]
    fn f_full(s: f64, d: f64) -> f64 {
        let diff = d * d - s * s;
        if diff.abs() < 1e-14 {
            return 0.0;
        }
        (1.0 / (4.0 * d)) * (((d - s) / (d + s)).ln() + 2.0 * s * d / diff)
    }
    /// Antiderivative of the cap angular integral: (1/(4d))·∫ s·(1/ri²-1/(s+d)²) ds
    #[inline]
    fn f_cap(s: f64, d: f64, ri: f64) -> f64 {
        (1.0 / (4.0 * d)) * (s * s / (2.0 * ri * ri) - (s + d).ln() - d / (s + d))
    }

    /// ψ_ij: desolvation of atom i (radius ri) by atom j (radius rj) at distance d.
    /// Finite for ALL d (no divergence at d=rj).
    #[inline]
    fn desolv(d: f64, ri: f64, rj: f64) -> f64 {
        let rsum = ri + rj;
        if d >= rsum {
            // Non-overlap: full sphere of j is outside i.
            let diff = d * d - rj * rj;
            if diff <= 1e-12 {
                return 0.0;
            }
            return (0.25 / d) * ((d - rj) / (d + rj)).ln() + rj / (2.0 * diff);
        }
        if d + rj <= ri {
            return 0.0; // j entirely inside i → no displaced dielectric.
        }
        // Partial overlap.
        if d <= ri {
            // Center of j inside i's sphere: cap from s=ri-d to s=rj.
            Self::f_cap(rj, d, ri) - Self::f_cap(ri - d, d, ri)
        } else {
            // Center of j outside i: full-shell [0, d-ri] + cap [d-ri, rj].
            Self::f_full(d - ri, d) + Self::f_cap(rj, d, ri) - Self::f_cap(d - ri, d, ri)
        }
    }

    /// ∂F_full/∂d (holding s constant).
    #[inline]
    fn f_full_dd(s: f64, d: f64) -> f64 {
        let diff = d * d - s * s;
        if diff.abs() < 1e-14 {
            return 0.0;
        }
        let lr = ((d - s) / (d + s)).ln();
        (-1.0 / (4.0 * d * d)) * (lr + 2.0 * s * d / diff) - s.powi(3) / (d * diff * diff)
    }
    /// ∂F_cap/∂d (holding s, ri constant).
    #[inline]
    fn f_cap_dd(s: f64, d: f64, ri: f64) -> f64 {
        let sp = s + d;
        let inner = s * s / (2.0 * ri * ri) - sp.ln() - d / sp;
        (-1.0 / (4.0 * d * d)) * inner + (1.0 / (4.0 * d)) * (-1.0 / sp - s / (sp * sp))
    }

    /// dψ_ij/dd: the desolvation derivative. Analytical for all overlap cases
    /// (Leibniz boundary terms cancel at the transition points).
    #[inline]
    fn desolv_d(d: f64, ri: f64, rj: f64) -> f64 {
        let rsum = ri + rj;
        if d >= rsum {
            // Non-overlap.
            let d2 = d * d;
            let diff = d2 - rj * rj;
            if diff <= 1e-12 {
                return 0.0;
            }
            let lr = ((d - rj) / (d + rj)).ln();
            return -lr / (4.0 * d2) + rj / (2.0 * d * diff) - rj * d / (diff * diff);
        }
        if d + rj <= ri {
            return 0.0;
        }
        // Partial overlap: dψ/dd = Σ ∂F/∂d at limits (boundary terms cancel).
        if d <= ri {
            Self::f_cap_dd(rj, d, ri) - Self::f_cap_dd(ri - d, d, ri)
        } else {
            Self::f_full_dd(d - ri, d) + Self::f_cap_dd(rj, d, ri) - Self::f_cap_dd(d - ri, d, ri)
        }
    }
}

impl<'a> ForceField for GBSA<'a> {
    fn energy_and_gradient(&self, coords: &[[f64; 3]], grad: &mut [[f64; 3]]) -> f64 {
        let e_ff = self.ff.energy_and_gradient(coords, grad);
        let n = self.n;
        let c = -KE * self.dielectric_factor;
        let qi = &self.charges;

        // --- Pass 1: Born radii (desolvation sums) ---
        let mut psi = vec![0.0f64; n];
        for i in 0..n {
            for j in 0..i {
                let dx = coords[i][0] - coords[j][0];
                let dy = coords[i][1] - coords[j][1];
                let dz = coords[i][2] - coords[j][2];
                let d = (dx * dx + dy * dy + dz * dz).sqrt();
                // Coincident atoms (degenerate input): the desolvation integral
                // has no meaningful value; skip rather than propagate NaN.
                if d < 1e-8 {
                    continue;
                }
                // Exact HCT desolvation (handles all overlap cases, no divergence).
                psi[i] += Self::desolv(d, self.radii[i], self.radii[j]);
                psi[j] += Self::desolv(d, self.radii[j], self.radii[i]);
            }
        }
        let mut born = vec![0.0f64; n];
        let mut chain = vec![0.0f64; n]; // d(1/R_i)/dψ_i
        for i in 0..n {
            let inv_ri = 1.0 / self.radii[i];
            let p = psi[i];
            let min_inv = inv_ri * 0.01; // safety floor: R_i ≤ 100·ρ_i
            let (born_inv, d_born_inv_dpsi) = if self.obc2 {
                let arg = OBC_ALPHA * p - OBC_BETA * p * p + OBC_GAMMA * p.powi(3);
                let unclamped = inv_ri * (1.0 - arg.tanh());
                let deriv = -inv_ri
                    * (1.0 - arg.tanh().powi(2))
                    * (OBC_ALPHA - 2.0 * OBC_BETA * p + 3.0 * OBC_GAMMA * p * p);
                if unclamped < min_inv {
                    (min_inv, 0.0) // clamp active: derivative is zero in this region
                } else {
                    (unclamped, deriv)
                }
            } else {
                let unclamped = inv_ri - p;
                if unclamped < min_inv {
                    (min_inv, 0.0)
                } else {
                    (unclamped, -1.0)
                }
            };
            born[i] = 1.0 / born_inv;
            chain[i] = d_born_inv_dpsi;
        }

        // --- Pass 2: GB energy + direct pair gradient + accumulate B_i ---
        let mut e_gb = 0.0;
        let mut b_coeff = vec![0.0f64; n]; // B_i = c·(½q_i² + Σ_j q_iq_j ∂(1/f)/∂(1/R_i))
        for i in 0..n {
            b_coeff[i] = c * 0.5 * qi[i] * qi[i]; // self-energy Born coefficient
            e_gb += 0.5 * qi[i] * qi[i] / born[i];
        }
        for i in 0..n {
            for j in 0..i {
                let dx = coords[i][0] - coords[j][0];
                let dy = coords[i][1] - coords[j][1];
                let dz = coords[i][2] - coords[j][2];
                let r2 = dx * dx + dy * dy + dz * dz;
                let r = r2.sqrt().max(1e-10);
                let a = born[i] * born[j];
                let exp_term = (-r2 / (4.0 * a)).exp();
                let f2 = r2 + a * exp_term;
                let f = f2.sqrt().max(1e-10);
                let qq = qi[i] * qi[j];
                e_gb += qq / f;

                // Direct pair gradient: ∂(1/f)/∂r · dr/dx
                let f3 = f * f2;
                let d_inv_f_du = -(1.0 - exp_term / 4.0) / (2.0 * f3); // d(1/f)/du, u=r²
                let g_dir = c * qq * d_inv_f_du * 2.0; // × du/dx = 2·dx
                grad[i][0] += g_dir * dx;
                grad[i][1] += g_dir * dy;
                grad[i][2] += g_dir * dz;
                grad[j][0] -= g_dir * dx;
                grad[j][1] -= g_dir * dy;
                grad[j][2] -= g_dir * dz;

                // Pair Born-radius coefficient: ∂(1/f)/∂R_i
                // ∂(f²)/∂R_i = R_j·E·(1 + r²/(4a));  ∂(1/f)/∂R_i = -that/(2·f³)
                let ri2 = born[i] * born[i];
                let rj2 = born[j] * born[j];
                let one_plus = 1.0 + r2 / (4.0 * a);
                let d_inv_f_dri = -born[j] * exp_term * one_plus / (2.0 * f3);
                let d_inv_f_drj = -born[i] * exp_term * one_plus / (2.0 * f3);
                // ∂(1/f)/∂(1/R_i) = ∂(1/f)/∂R_i · (-R_i²)
                b_coeff[i] += c * qq * d_inv_f_dri * (-ri2);
                b_coeff[j] += c * qq * d_inv_f_drj * (-rj2);
            }
        }
        e_gb *= c;

        // --- Pass 3: Born-radius gradient (B_i · chain[i] · dψ_i/dx) ---
        for i in 0..n {
            for j in 0..i {
                let dx = coords[i][0] - coords[j][0];
                let dy = coords[i][1] - coords[j][1];
                let dz = coords[i][2] - coords[j][2];
                let r = (dx * dx + dy * dy + dz * dz).sqrt();
                if r < 1e-8 {
                    continue; // coincident atoms: no desolvation gradient
                }
                let dpsi_ij = Self::desolv_d(r, self.radii[i], self.radii[j]); // dψ_i/dr
                let dpsi_ji = Self::desolv_d(r, self.radii[j], self.radii[i]); // dψ_j/dr
                let inv_r = 1.0 / r;
                let g_brn =
                    (b_coeff[i] * chain[i] * dpsi_ij + b_coeff[j] * chain[j] * dpsi_ji) * inv_r;
                grad[i][0] += g_brn * dx;
                grad[i][1] += g_brn * dy;
                grad[i][2] += g_brn * dz;
                grad[j][0] -= g_brn * dx;
                grad[j][1] -= g_brn * dy;
                grad[j][2] -= g_brn * dz;
            }
        }

        // --- SA (surface area) nonpolar term via LCPO ---
        let mut e_sa = 0.0;
        if self.sa_tension > 0.0 {
            let g = self.sa_tension;
            // Pass 1: per-atom LCPO sums (Σ Sij, Σ Sij²) for the energy.
            let mut s_sum = vec![0.0f64; n];
            let mut s_sq_sum = vec![0.0f64; n];
            for i in 0..n {
                for j in 0..i {
                    let dx = coords[i][0] - coords[j][0];
                    let dy = coords[i][1] - coords[j][1];
                    let dz = coords[i][2] - coords[j][2];
                    let d = (dx * dx + dy * dy + dz * dz).sqrt();
                    let ri = self.sa_radii[i];
                    let rj = self.sa_radii[j];
                    if d >= ri + rj {
                        continue;
                    }
                    // Sij = π·Ri·(Ri - xi), xi = (d²+Ri²-Rj²)/(2d)
                    let sij = if d + rj <= ri {
                        0.0
                    } else {
                        let xi = (d * d + ri * ri - rj * rj) / (2.0 * d);
                        std::f64::consts::PI * ri * (ri - xi).max(0.0)
                    };
                    let sji = if d + ri <= rj {
                        0.0
                    } else {
                        let xj = (d * d + rj * rj - ri * ri) / (2.0 * d);
                        std::f64::consts::PI * rj * (rj - xj).max(0.0)
                    };
                    s_sum[i] += sij;
                    s_sq_sum[i] += sij * sij;
                    s_sum[j] += sji;
                    s_sq_sum[j] += sji * sji;
                }
            }
            for i in 0..n {
                let si = 4.0 * std::f64::consts::PI * self.sa_radii[i].powi(2);
                e_sa += self.lcpo[i][0] * si
                    + self.lcpo[i][1] * s_sum[i]
                    + self.lcpo[i][2] * s_sq_sum[i];
            }
            e_sa *= g;

            // Pass 2: SA gradient. dSij/dd = -π·Ri·(d²-Ri²+Rj²)/(2d²)
            for i in 0..n {
                for j in 0..i {
                    let dx = coords[i][0] - coords[j][0];
                    let dy = coords[i][1] - coords[j][1];
                    let dz = coords[i][2] - coords[j][2];
                    let d = (dx * dx + dy * dy + dz * dz).sqrt().max(1e-10);
                    let ri = self.sa_radii[i];
                    let rj = self.sa_radii[j];
                    if d >= ri + rj {
                        continue;
                    }
                    let inv_d = 1.0 / d;
                    // Sij and its derivative
                    let (sij, dsij_dd) = if d + rj <= ri {
                        (0.0, 0.0)
                    } else {
                        // NOTE: S_ij is identically 0 for full containment (d + rj <= ri)
                        // while dS_ij/dd is nonzero just outside that boundary, so the
                        // pairwise gradient has a small kink at d = ri - rj. This is
                        // inherent to the LCPO pairwise overlap model (Weiser-Tsui-Case)
                        // and is preserved for fidelity; atoms only cross it in
                        // unphysical buried configurations.
                        let xi = (d * d + ri * ri - rj * rj) / (2.0 * d);
                        let s = std::f64::consts::PI * ri * (ri - xi).max(0.0);
                        let ds = -std::f64::consts::PI * ri * (d * d - ri * ri + rj * rj)
                            / (2.0 * d * d);
                        (s, ds)
                    };
                    let (sji, dsji_dd) = if d + ri <= rj {
                        (0.0, 0.0)
                    } else {
                        let xj = (d * d + rj * rj - ri * ri) / (2.0 * d);
                        let s = std::f64::consts::PI * rj * (rj - xj).max(0.0);
                        let ds = -std::f64::consts::PI * rj * (d * d - rj * rj + ri * ri)
                            / (2.0 * d * d);
                        (s, ds)
                    };
                    // Gradient coefficient: (a2_k + 2·a3_k·Skj)·dSkj/dd + cross-term
                    let coef_i = g * (self.lcpo[i][1] + 2.0 * self.lcpo[i][2] * sij) * dsij_dd;
                    let coef_j = g * (self.lcpo[j][1] + 2.0 * self.lcpo[j][2] * sji) * dsji_dd;
                    let g_sa = (coef_i + coef_j) * inv_d;
                    grad[i][0] += g_sa * dx;
                    grad[i][1] += g_sa * dy;
                    grad[i][2] += g_sa * dz;
                    grad[j][0] -= g_sa * dx;
                    grad[j][1] -= g_sa * dy;
                    grad[j][2] -= g_sa * dz;
                }
            }
        }

        e_ff + e_gb + e_sa
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mmff::MMFFForceField;
    use crate::molecule::parser::parse_sdf;

    struct ZeroFF;
    impl ForceField for ZeroFF {
        fn energy_and_gradient(&self, _c: &[[f64; 3]], g: &mut [[f64; 3]]) -> f64 {
            for r in g.iter_mut() {
                *r = [0.0; 3];
            }
            0.0
        }
    }

    #[test]
    fn gb_self_energy_matches_born_equation() {
        let q = 1.0_f64;
        let r_born = 1.75 + 0.195;
        let expected = -166.053 * (1.0 - 1.0 / 78.5) * q * q / r_born;
        let zff = ZeroFF;
        let mol = crate::molecule::Molecule {
            atoms: vec![crate::molecule::Atom {
                symbol: "Cl".into(),
                atomic_number: 17,
                mass: 35.0,
                charge: q,
                position: [0.0; 3],
                index: 0,
                stereo_parity: 0,
            }],
            bonds: vec![],
            name: String::new(),
            adjacency: vec![vec![]],
        };
        let gbsa = GBSA::new(&zff, &mol, &[q], &GBSAConfig::default());
        let e = gbsa.energy_and_gradient(&[[0.0; 3]], &mut [[0.0; 3]]);
        assert!(
            (e - expected).abs() / expected.abs() < 0.02,
            "GB self-energy {:.4} vs Born {:.4}",
            e,
            expected
        );
    }

    #[test]
    fn gb_gradient_matches_finite_difference() {
        let sdf = std::fs::read_to_string("scripts/val_set/ethanol.sdf").unwrap();
        let mol = parse_sdf(&sdf).unwrap();
        let ff = MMFFForceField::new(&mol, crate::MMFFVariant::MMFF94s);
        // Test with SA enabled (full GBSA).
        let cfg = GBSAConfig {
            sa_surface_tension: 0.00542,
            ..GBSAConfig::default()
        };
        let gbsa = GBSA::new(&ff, &mol, &ff.charges, &cfg);
        let coords: Vec<[f64; 3]> = mol.atoms.iter().map(|a| a.position).collect();
        let mut grad = vec![[0.0; 3]; mol.atoms.len()];
        gbsa.energy_and_gradient(&coords, &mut grad);
        let h = 5e-5;
        let mut max_err = 0.0f64;
        for i in 0..mol.atoms.len() {
            for k in 0..3 {
                let mut p = coords.clone();
                let mut m = coords.clone();
                p[i][k] += h;
                m[i][k] -= h;
                let mut g0 = vec![[0.0; 3]; mol.atoms.len()];
                let mut g1 = vec![[0.0; 3]; mol.atoms.len()];
                let ep = gbsa.energy_and_gradient(&p, &mut g0);
                let em = gbsa.energy_and_gradient(&m, &mut g1);
                let fd = (ep - em) / (2.0 * h);
                max_err = max_err.max((grad[i][k] - fd).abs());
            }
        }
        eprintln!("GB+SA gradient max FD error: {:.6}", max_err);
        assert!(
            max_err < 1e-4,
            "GB gradient max FD error {:.6} > 1e-4",
            max_err
        );
    }

    #[test]
    fn gb_lowers_energy_for_polar_molecule() {
        let sdf = std::fs::read_to_string("scripts/val_set/ethanol.sdf").unwrap();
        let mol = parse_sdf(&sdf).unwrap();
        let ff = MMFFForceField::new(&mol, crate::MMFFVariant::MMFF94s);
        let coords: Vec<[f64; 3]> = mol.atoms.iter().map(|a| a.position).collect();
        let e_vac = ff.calculate_energy_and_gradient(&coords).0;
        let gbsa = GBSA::new(&ff, &mol, &ff.charges, &GBSAConfig::default());
        let mut g = vec![[0.0; 3]; mol.atoms.len()];
        let e_solv = gbsa.energy_and_gradient(&coords, &mut g);
        assert!(
            e_solv < e_vac,
            "GBSA {:.2} should be < vacuum {:.2}",
            e_solv,
            e_vac
        );
    }

    /// GBSA NVE stability smoke test. The GB gradient is verified correct by
    /// finite-difference (gb_gradient_matches_finite_difference). This test
    /// GBSA NVE energy conservation. The overlap integral eliminates the HCT
    /// divergence and cutoff discontinuity, so NVE should be stable.
    #[test]
    fn gb_nve_stays_finite() {
        use crate::md::{MDConfig, MDRunner};
        let sdf = std::fs::read_to_string("scripts/val_set/ethanol.sdf").unwrap();
        let mol = parse_sdf(&sdf).unwrap();
        let ff = MMFFForceField::new(&mol, crate::MMFFVariant::MMFF94s);
        let gbsa = GBSA::new(&ff, &mol, &ff.charges, &GBSAConfig::default());
        let mut runner = MDRunner::from_molecule(
            std::rc::Rc::new(gbsa),
            &mol,
            MDConfig {
                dt_fs: 0.25,
                temperature_k: 300.0,
                friction_per_ps: 0.0,
                seed: 42,
            },
        );
        let e0 = runner.potential_energy() + runner.kinetic_energy();
        let mut max_drift: f64 = 0.0;
        for _ in 0..2000 {
            runner.step();
            let e = runner.potential_energy() + runner.kinetic_energy();
            assert!(e.is_finite());
            max_drift = max_drift.max((e - e0).abs() / e0.abs().max(1.0));
        }
        eprintln!("GB NVE drift over 2000 steps: {:.2}%", max_drift * 100.0);
        assert!(
            max_drift < 0.01,
            "GB NVE drift {:.2}% > 1%",
            max_drift * 100.0
        );
    }
}
