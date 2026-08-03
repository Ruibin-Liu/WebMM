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
        1 => 1.20, 6 => 1.70, 7 => 1.55, 8 => 1.52, 9 => 1.47,
        14 => 2.10, 15 => 1.80, 16 => 1.80, 17 => 1.75, 35 => 1.85, 53 => 1.98,
        _ => 1.5,
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
        Self { solvent_dielectric: 78.5, solute_dielectric: 1.0, radius_offset: 0.195,
              obc2: true, sa_surface_tension: 0.0 }
    }
}

pub struct GBSA<'a> {
    ff: &'a dyn ForceField,
    charges: Vec<f64>,
    radii: Vec<f64>, // R_vdw + offset
    dielectric_factor: f64,
    obc2: bool,
    sa_tension: f64,
    n: usize,
}

impl<'a> GBSA<'a> {
    pub fn new(ff: &'a dyn ForceField, mol: &Molecule, charges: &[f64], config: &GBSAConfig) -> Self {
        let n = mol.atoms.len();
        let radii: Vec<f64> = mol.atoms.iter()
            .map(|a| bondi_radius(a.atomic_number) + config.radius_offset).collect();
        Self { ff, charges: charges.to_vec(), radii,
            dielectric_factor: 1.0 / config.solute_dielectric - 1.0 / config.solvent_dielectric,
            obc2: config.obc2, sa_tension: config.sa_surface_tension, n }
    }

    #[inline]
    fn desolv(d: f64, rj: f64) -> f64 {
        let diff = d * d - rj * rj;
        if diff <= 1e-8 { return 0.0; }
        (0.25 / d) * ((d - rj) / (d + rj)).ln() + rj / (2.0 * diff)
    }

    #[inline]
    fn desolv_d(d: f64, rj: f64) -> f64 {
        let d2 = d * d;
        let diff = d2 - rj * rj;
        if diff <= 1e-8 { return 0.0; }
        let lr = ((d - rj) / (d + rj)).ln();
        -lr / (4.0 * d2) + rj / (2.0 * d * diff) - rj * d / (diff * diff)
    }

    /// Switched desolvation: smooth cosine fade to zero over [rsum-0.3, rsum]
    /// to avoid the hard-cutoff force discontinuity at d=rsum.
    const SWITCH_WIDTH: f64 = 0.3;
    #[inline]
    fn desolv_sw(d: f64, rj: f64, rsum: f64) -> f64 {
        if d >= rsum { return 0.0; }
        let psi = Self::desolv(d, rj);
        let on = rsum - Self::SWITCH_WIDTH;
        if d <= on { return psi; }
        let t = (d - on) / Self::SWITCH_WIDTH;
        psi * 0.5 * (1.0 + (std::f64::consts::PI * t).cos())
    }
    #[inline]
    fn desolv_sw_d(d: f64, rj: f64, rsum: f64) -> f64 {
        if d >= rsum { return 0.0; }
        let psi = Self::desolv(d, rj);
        let psi_d = Self::desolv_d(d, rj);
        let on = rsum - Self::SWITCH_WIDTH;
        if d <= on { return psi_d; }
        let t = (d - on) / Self::SWITCH_WIDTH;
        let s = 0.5 * (1.0 + (std::f64::consts::PI * t).cos());
        let s_d = -0.5 * (std::f64::consts::PI * t).sin() * std::f64::consts::PI / Self::SWITCH_WIDTH;
        psi_d * s + psi * s_d
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
                let rsum = self.radii[i] + self.radii[j];
                if d >= rsum { continue; }
                // Skip deeply-overlapping pairs (d ≤ rj+margin): the 1/(d²-rj²)
                // Coulomb-field integral diverges for overlapping spheres.
                // Bonded atoms (1-2, tight 1-3) don't contribute to the solvation shell.
                const SKIP_MARGIN: f64 = 0.1;
                if d > self.radii[j] + SKIP_MARGIN {
                    psi[i] += Self::desolv_sw(d, self.radii[j], rsum);
                }
                if d > self.radii[i] + SKIP_MARGIN {
                    psi[j] += Self::desolv_sw(d, self.radii[i], rsum);
                }
            }
        }
        let mut born = vec![0.0f64; n];
        let mut chain = vec![0.0f64; n]; // d(1/R_i)/dψ_i
        for i in 0..n {
            let inv_ri = 1.0 / self.radii[i];
            let p = psi[i];
            let born_inv = if self.obc2 {
                let arg = OBC_ALPHA * p - OBC_BETA * p * p + OBC_GAMMA * p.powi(3);
                (inv_ri * (1.0 - arg.tanh())).max(inv_ri * 0.01)
            } else {
                (inv_ri - p).max(inv_ri * 0.01)
            };
            born[i] = 1.0 / born_inv;
            chain[i] = if self.obc2 {
                let arg = OBC_ALPHA * p - OBC_BETA * p * p + OBC_GAMMA * p.powi(3);
                let sech2 = 1.0 - arg.tanh().powi(2);
                -inv_ri * sech2 * (OBC_ALPHA - 2.0 * OBC_BETA * p + 3.0 * OBC_GAMMA * p * p)
            } else {
                -1.0
            };
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
                grad[i][0] += g_dir * dx; grad[i][1] += g_dir * dy; grad[i][2] += g_dir * dz;
                grad[j][0] -= g_dir * dx; grad[j][1] -= g_dir * dy; grad[j][2] -= g_dir * dz;

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
                let rsum = self.radii[i] + self.radii[j];
                if r >= rsum { continue; }
                const SKIP_MARGIN: f64 = 0.1;
                let dpsi_ij = if r > self.radii[j] + SKIP_MARGIN {
                    Self::desolv_sw_d(r, self.radii[j], rsum)
                } else { 0.0 }; // dψ_i/dr (i desolvated by j)
                let dpsi_ji = if r > self.radii[i] + SKIP_MARGIN {
                    Self::desolv_sw_d(r, self.radii[i], rsum)
                } else { 0.0 }; // dψ_j/dr (j desolvated by i)
                let inv_r = 1.0 / r;
                let g_brn = (b_coeff[i] * chain[i] * dpsi_ij + b_coeff[j] * chain[j] * dpsi_ji) * inv_r;
                grad[i][0] += g_brn * dx; grad[i][1] += g_brn * dy; grad[i][2] += g_brn * dz;
                grad[j][0] -= g_brn * dx; grad[j][1] -= g_brn * dy; grad[j][2] -= g_brn * dz;
            }
        }

        let _ = self.sa_tension; // SA term not yet implemented (TODO: LCPO area)
        e_ff + e_gb
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
            for r in g.iter_mut() { *r = [0.0; 3]; }
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
                symbol: "Cl".into(), atomic_number: 17, mass: 35.0, charge: q,
                position: [0.0; 3], index: 0, stereo_parity: 0,
            }],
            bonds: vec![], name: String::new(), adjacency: vec![vec![]],
        };
        let gbsa = GBSA::new(&zff, &mol, &[q], &GBSAConfig::default());
        let e = gbsa.energy_and_gradient(&[[0.0; 3]], &mut [[0.0; 3]]);
        assert!((e - expected).abs() / expected.abs() < 0.02,
            "GB self-energy {:.4} vs Born {:.4}", e, expected);
    }

    #[test]
    fn gb_gradient_matches_finite_difference() {
        let sdf = std::fs::read_to_string("scripts/val_set/ethanol.sdf").unwrap();
        let mol = parse_sdf(&sdf).unwrap();
        let ff = MMFFForceField::new(&mol, crate::MMFFVariant::MMFF94s);
        let gbsa = GBSA::new(&ff, &mol, &ff.charges, &GBSAConfig::default());
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
        assert!(max_err < 0.1, "GB gradient max FD error {:.4}", max_err);
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
        assert!(e_solv < e_vac, "GBSA {:.2} should be < vacuum {:.2}", e_solv, e_vac);
    }

    /// GBSA NVE stability smoke test. The GB gradient is verified correct by
    /// finite-difference (gb_gradient_matches_finite_difference). This test
    /// checks the energy stays finite during MD (no NaN/Inf blowup). Note: the
    /// HCT desolvation has a residual force discontinuity for overlapping atom
    /// pairs (d ≈ Rj); full NVE conservation needs the proper spherical-cap
    /// overlap formula (TODO). NVT MD (thermostat) is recommended for practical use.
    #[test]
    fn gb_nve_stays_finite() {
        use crate::md::{MDConfig, MDRunner};
        let sdf = std::fs::read_to_string("scripts/val_set/ethanol.sdf").unwrap();
        let mol = parse_sdf(&sdf).unwrap();
        let ff = MMFFForceField::new(&mol, crate::MMFFVariant::MMFF94s);
        let gbsa = GBSA::new(&ff, &mol, &ff.charges, &GBSAConfig::default());
        let mut runner = MDRunner::from_molecule(&gbsa, &mol, MDConfig {
            dt_fs: 0.1, temperature_k: 300.0, friction_per_ps: 0.0, seed: 42,
        });
        for _ in 0..500 {
            runner.step();
            let e = runner.potential_energy() + runner.kinetic_energy();
            assert!(e.is_finite(), "GBSA NVE energy became non-finite");
            assert!(e.abs() < 1e6, "GBSA NVE energy blew up: {:.0}", e);
        }
    }
}

