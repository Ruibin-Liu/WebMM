//! Molecular dynamics engine for gas-phase conformational exploration.
//!
//! Integrators: velocity-Verlet (NVE) and BAOAB Langevin (NVT). Generic over any
//! [`crate::forces::ForceField`]. No nonbonded cutoff (all-pairs, faithful to the
//! underlying potential — important for later metadynamics, which bias-explores the
//! exact surface).
//!
//! # Units
//! Internal unit system: length Å, mass amu, energy kcal/mol, time τ where
//! `1 τ = TIMEUNIT_FS fs = sqrt(amu·Å²·mol/kcal) ≈ 48.885 fs`. In these units the
//! equations of motion need *no* conversion factors: acceleration `a = F/m` comes
//! out in Å/τ², and kinetic energy `0.5 Σ m v²` comes out in kcal/mol, so
//! temperature follows with kB in kcal/(mol·K).

use crate::forces::ForceField;
use std::rc::Rc;

/// fs per internal time unit τ (τ = sqrt(amu·Å²·mol/kcal)).
const TIMEUNIT_FS: f64 = 48.885;
/// Boltzmann constant, kcal/(mol·K).
pub const KB: f64 = 0.0019872041;

/// Atomic mass (amu) by atomic number. Covers the common MMFF elements; falls
/// back to a rough estimate for anything missing.
pub fn atomic_mass(z: u8) -> f64 {
    match z {
        1 => 1.008,    // H
        5 => 10.81,    // B
        6 => 12.011,   // C
        7 => 14.007,   // N
        8 => 15.999,   // O
        9 => 18.998,   // F
        11 => 22.990,  // Na
        12 => 24.305,  // Mg
        13 => 26.982,  // Al
        14 => 28.085,  // Si
        15 => 30.974,  // P
        16 => 32.06,   // S
        17 => 35.45,   // Cl
        19 => 39.098,  // K
        20 => 40.078,  // Ca
        25 => 54.938,  // Mn
        26 => 55.845,  // Fe
        30 => 65.38,   // Zn
        35 => 79.904,  // Br
        53 => 126.904, // I
        _ => 2.0 * z as f64,
    }
}

/// A minimal, deterministic PRNG (splitmix64-seeded xoshiro256**) for Langevin
/// noise. Reproducible given the seed.
struct Rng {
    s: [u64; 4],
}
impl Rng {
    fn new(seed: u64) -> Self {
        // splitmix64 to expand a u64 seed into the 256-bit state.
        let mut z = seed;
        let mut next = || {
            z = z.wrapping_add(0x9E37_79B9_7F4A_7C15);
            let mut x = z;
            x = (x ^ (x >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
            x = (x ^ (x >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
            (x ^ (x >> 31)) | 1 // ensure nonzero state
        };
        Rng {
            s: [next(), next(), next(), next()],
        }
    }
    fn next_u64(&mut self) -> u64 {
        // xoshiro256**
        let result = self.s[1].wrapping_mul(5).rotate_left(7).wrapping_mul(9);
        let t = self.s[1] << 17;
        self.s[2] ^= self.s[0];
        self.s[3] ^= self.s[1];
        self.s[1] ^= self.s[2];
        self.s[0] ^= self.s[3];
        self.s[2] ^= t;
        self.s[3] = self.s[3].rotate_left(45);
        result
    }
    fn next_f64(&mut self) -> f64 {
        // uniform in [0,1)
        (self.next_u64() >> 11) as f64 / ((1u64 << 53) as f64)
    }
    fn next_gaussian(&mut self) -> f64 {
        // Box–Muller → N(0,1)
        let u1 = self.next_f64().max(1e-300);
        let u2 = self.next_f64();
        (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos()
    }
}

/// MD configuration. Sensible defaults: 1 fs step, 300 K, friction 1/ps.
#[derive(Debug, Clone, Copy)]
pub struct MDConfig {
    /// Integration timestep in femtoseconds.
    pub dt_fs: f64,
    /// Target temperature (K). Velocities are initialized at this temperature.
    /// The thermostat runs only if `friction_per_ps > 0`.
    pub temperature_k: f64,
    /// Langevin friction (1/ps). `> 0` → NVT (BAOAB); `<= 0` → NVE.
    pub friction_per_ps: f64,
    /// RNG seed (Langevin noise + velocity initialization). Reproducible.
    pub seed: u64,
}
impl Default for MDConfig {
    fn default() -> Self {
        MDConfig {
            dt_fs: 1.0,
            temperature_k: 300.0,
            friction_per_ps: 1.0,
            seed: 42,
        }
    }
}

/// A molecular dynamics runner. Owns positions, velocities, and the reusable
/// gradient buffer; advances with one force evaluation per step.
pub struct MDRunner<'a> {
    ff: Rc<dyn ForceField + 'a>,
    n: usize,
    masses: Vec<f64>,
    coords: Vec<[f64; 3]>,
    vel: Vec<[f64; 3]>,
    grad: Vec<[f64; 3]>, // last gradient (dE/dx); forces = -grad
    pe: f64,
    dt_tau: f64,
    temperature_k: f64,
    friction_tau: f64,
    use_thermostat: bool,
    dof: usize,
    rng: Rng,
    step: usize,
    time_fs: f64,
}

impl<'a> MDRunner<'a> {
    /// Build a runner from explicit masses + starting coordinates.
    pub fn new(
        ff: Rc<dyn ForceField + 'a>,
        masses: Vec<f64>,
        coords: Vec<[f64; 3]>,
        config: MDConfig,
    ) -> Self {
        let n = coords.len();
        let mut grad = vec![[0.0; 3]; n];
        let pe = ff.energy_and_gradient(&coords, &mut grad);
        let dof = if n > 1 { 3 * n - 3 } else { 3 }; // remove COM translation
        let mut rng = Rng::new(config.seed);
        let mut vel = vec![[0.0; 3]; n];
        if config.temperature_k > 0.0 {
            // Maxwell–Boltzmann initialization
            let s = (KB * config.temperature_k).sqrt();
            for i in 0..n {
                let a = s / masses[i].sqrt();
                vel[i] = [
                    a * rng.next_gaussian(),
                    a * rng.next_gaussian(),
                    a * rng.next_gaussian(),
                ];
            }
            if n > 1 {
                remove_com_velocity(&mut vel, &masses);
            }
            // rescale to exactly the target temperature
            let ke = kinetic_energy(&masses, &vel);
            let t = 2.0 * ke / (dof as f64 * KB);
            if t > 1e-12 {
                let sc = (config.temperature_k / t).sqrt();
                for v in vel.iter_mut() {
                    v[0] *= sc;
                    v[1] *= sc;
                    v[2] *= sc;
                }
            }
        }
        MDRunner {
            ff,
            n,
            masses,
            coords,
            vel,
            grad,
            pe,
            dt_tau: config.dt_fs / TIMEUNIT_FS,
            temperature_k: config.temperature_k,
            friction_tau: config.friction_per_ps * TIMEUNIT_FS / 1000.0,
            use_thermostat: config.friction_per_ps > 0.0 && config.temperature_k > 0.0,
            dof,
            rng,
            step: 0,
            time_fs: 0.0,
        }
    }

    /// Build a runner from a molecule (masses + coordinates taken from atoms).
    pub fn from_molecule(
        ff: Rc<dyn ForceField + 'a>,
        mol: &crate::molecule::Molecule,
        config: MDConfig,
    ) -> Self {
        let masses: Vec<f64> = mol
            .atoms
            .iter()
            .map(|a| atomic_mass(a.atomic_number))
            .collect();
        let coords: Vec<[f64; 3]> = mol.atoms.iter().map(|a| a.position).collect();
        Self::new(ff, masses, coords, config)
    }

    /// Advance one step (NVE velocity-Verlet or NVT BAOAB Langevin).
    /// Returns the current potential energy (kcal/mol). One force evaluation.
    pub fn step(&mut self) -> f64 {
        let n = self.n;
        let dt = self.dt_tau;
        if self.use_thermostat {
            let c1 = (-self.friction_tau * dt).exp();
            let c2 = (1.0 - c1 * c1).sqrt();
            let s = (KB * self.temperature_k).sqrt();
            // B
            for i in 0..n {
                let invm = 1.0 / self.masses[i];
                for c in 0..3 {
                    self.vel[i][c] -= 0.5 * self.grad[i][c] * invm * dt;
                }
            }
            // A
            for i in 0..n {
                for c in 0..3 {
                    self.coords[i][c] += 0.5 * self.vel[i][c] * dt;
                }
            }
            // O
            for i in 0..n {
                let amp = c2 * s * (1.0 / self.masses[i]).sqrt();
                for c in 0..3 {
                    self.vel[i][c] = c1 * self.vel[i][c] + amp * self.rng.next_gaussian();
                }
            }
            // A
            for i in 0..n {
                for c in 0..3 {
                    self.coords[i][c] += 0.5 * self.vel[i][c] * dt;
                }
            }
            self.pe = self.ff.energy_and_gradient(&self.coords, &mut self.grad);
            // B
            for i in 0..n {
                let invm = 1.0 / self.masses[i];
                for c in 0..3 {
                    self.vel[i][c] -= 0.5 * self.grad[i][c] * invm * dt;
                }
            }
        } else {
            // velocity-Verlet
            for i in 0..n {
                let invm = 1.0 / self.masses[i];
                for c in 0..3 {
                    self.vel[i][c] -= 0.5 * self.grad[i][c] * invm * dt;
                    self.coords[i][c] += self.vel[i][c] * dt;
                }
            }
            self.pe = self.ff.energy_and_gradient(&self.coords, &mut self.grad);
            for i in 0..n {
                let invm = 1.0 / self.masses[i];
                for c in 0..3 {
                    self.vel[i][c] -= 0.5 * self.grad[i][c] * invm * dt;
                }
            }
        }
        self.step += 1;
        self.time_fs += dt * TIMEUNIT_FS;
        self.pe
    }

    /// Potential energy from the last force evaluation (kcal/mol).
    pub fn potential_energy(&self) -> f64 {
        self.pe
    }
    /// Instantaneous kinetic energy (kcal/mol).
    pub fn kinetic_energy(&self) -> f64 {
        kinetic_energy(&self.masses, &self.vel)
    }
    /// Instantaneous temperature (K) from `2·KE / (dof·kB)`.
    pub fn temperature(&self) -> f64 {
        2.0 * self.kinetic_energy() / (self.dof as f64 * KB)
    }
    /// Total energy PE + KE (conserved in NVE).
    pub fn total_energy(&self) -> f64 {
        self.pe + self.kinetic_energy()
    }
    pub fn coords(&self) -> &[[f64; 3]] {
        &self.coords
    }
    pub fn velocities(&self) -> &[[f64; 3]] {
        &self.vel
    }
    pub fn step_count(&self) -> usize {
        self.step
    }
    pub fn time_fs(&self) -> f64 {
        self.time_fs
    }
    pub fn dof(&self) -> usize {
        self.dof
    }
}

fn kinetic_energy(masses: &[f64], vel: &[[f64; 3]]) -> f64 {
    let mut ke = 0.0;
    for i in 0..vel.len() {
        let m = masses[i];
        let v = vel[i];
        ke += 0.5 * m * (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    }
    ke
}

/// Subtract the mass-weighted center-of-mass velocity from every atom.
fn remove_com_velocity(vel: &mut [[f64; 3]], masses: &[f64]) {
    let mtot: f64 = masses.iter().sum();
    if mtot < 1e-12 {
        return;
    }
    let mut pcom = [0.0; 3];
    for i in 0..vel.len() {
        let m = masses[i];
        pcom[0] += m * vel[i][0];
        pcom[1] += m * vel[i][1];
        pcom[2] += m * vel[i][2];
    }
    for v in pcom.iter_mut() {
        *v /= mtot;
    }
    for v in vel.iter_mut() {
        v[0] -= pcom[0];
        v[1] -= pcom[1];
        v[2] -= pcom[2];
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::rc::Rc;

    /// Test force field: N independent 3-D harmonic oscillators, E = 0.5·k·|x|².
    /// Exact, so it validates the integrator + units without depending on MMFF.
    struct HarmonicOscillator {
        n: usize,
        k: f64,
    }
    impl ForceField for HarmonicOscillator {
        fn energy_and_gradient(&self, coords: &[[f64; 3]], grad: &mut [[f64; 3]]) -> f64 {
            for g in grad.iter_mut() {
                *g = [0.0; 3];
            }
            let mut e = 0.0;
            for i in 0..self.n {
                for c in 0..3 {
                    e += 0.5 * self.k * coords[i][c] * coords[i][c];
                    grad[i][c] += self.k * coords[i][c]; // dE/dx = k·x
                }
            }
            e
        }
    }

    #[test]
    fn nve_energy_is_conserved() {
        // Velocity-Verlet on a harmonic oscillator conserves energy to O(dt²).
        let ho = HarmonicOscillator { n: 8, k: 1.0 };
        let masses = vec![1.0; 8];
        let coords = (0..8)
            .map(|i| [i as f64, (i as f64) * 0.5, -0.3])
            .collect::<Vec<_>>();
        // temperature_k>0 initializes velocities; friction=0 → NVE dynamics.
        let mut md = MDRunner::new(
            Rc::new(ho),
            masses,
            coords,
            MDConfig {
                dt_fs: 0.5,
                temperature_k: 300.0,
                friction_per_ps: 0.0,
                seed: 1,
            },
        );
        let e0 = md.total_energy();
        for _ in 0..20_000 {
            md.step();
        }
        let e1 = md.total_energy();
        let rel = ((e1 - e0) / e0).abs();
        assert!(
            rel < 1e-3,
            "NVE energy drift too large: {} -> {} (rel {})",
            e0,
            e1,
            rel
        );
    }

    #[test]
    fn langevin_equilibrates_to_target_temperature() {
        let ho = HarmonicOscillator { n: 10, k: 1.0 };
        let masses = vec![1.0; 10];
        let coords = vec![[0.0; 3]; 10];
        let mut md = MDRunner::new(
            Rc::new(ho),
            masses,
            coords,
            MDConfig {
                dt_fs: 1.0,
                temperature_k: 300.0,
                friction_per_ps: 5.0,
                seed: 7,
            },
        );
        let mut t_sum = 0.0;
        let mut n = 0;
        for i in 0..30_000 {
            md.step();
            if i > 8_000 {
                t_sum += md.temperature();
                n += 1;
            }
        }
        let t_avg = t_sum / n as f64;
        assert!(
            (t_avg - 300.0).abs() / 300.0 < 0.10,
            "Langevin avg T {} vs 300",
            t_avg
        );
    }

    #[test]
    fn mmff_md_runs_and_stays_sane() {
        let sdf = "Methane\n     RDKit          3D\n\n  5  4  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0\n    0.6294    0.6294    0.6294 H   0  0  0  0  0  0  0  0  0\n   -0.6294   -0.6294    0.6294 H   0  0  0  0  0  0  0  0  0\n   -0.6294    0.6294   -0.6294 H   0  0  0  0  0  0  0  0  0\n    0.6294   -0.6294   -0.6294 H   0  0  0  0  0  0  0  0  0\n  1  2  1  0\n  1  3  1  0\n  1  4  1  0\n  1  5  1  0\nM  END\n";
        let mol = crate::molecule::parser::parse_sdf(sdf).expect("parse");
        let ff = crate::mmff::MMFFForceField::new(&mol, crate::MMFFVariant::MMFF94s);
        let mut md = MDRunner::from_molecule(
            Rc::new(ff),
            &mol,
            MDConfig {
                dt_fs: 0.5,
                temperature_k: 300.0,
                friction_per_ps: 2.0,
                seed: 3,
            },
        );
        for _ in 0..2000 {
            md.step();
        }
        let t = md.temperature();
        let coords = md.coords();
        assert!(
            t.is_finite() && (50.0..=800.0).contains(&t),
            "MMFF MD temperature {} out of sane range",
            t
        );
        assert!(coords
            .iter()
            .all(|p| p[0].is_finite() && p[1].is_finite() && p[2].is_finite()));
    }

    const METHANE_SDF: &str = "Methane\n     RDKit          3D\n\n  5  4  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0\n    0.6294    0.6294    0.6294 H   0  0  0  0  0  0  0  0  0\n   -0.6294   -0.6294    0.6294 H   0  0  0  0  0  0  0  0  0\n   -0.6294    0.6294   -0.6294 H   0  0  0  0  0  0  0  0  0\n    0.6294   -0.6294   -0.6294 H   0  0  0  0  0  0  0  0  0\n  1  2  1  0\n  1  3  1  0\n  1  4  1  0\n  1  5  1  0\nM  END\n";

    #[test]
    fn wasm_run_md_from_sdf_returns_trajectory() {
        let mut opts = crate::MDOptions::new();
        opts.set_n_steps(200);
        opts.set_snapshot_interval(50);
        opts.set_friction_per_ps(2.0);
        let res = crate::run_md_from_sdf(METHANE_SDF, opts);
        assert!(res.success(), "{}", res.error());
        assert_eq!(res.n_atoms(), 5);
        // frames: initial + steps 50/100/150 + final = 5
        assert_eq!(res.n_frames(), 5);
        assert_eq!(res.coordinates().len(), res.n_frames() * res.n_atoms() * 3);
        assert_eq!(res.energies().len(), res.n_frames());
        assert_eq!(res.temperatures().len(), res.n_frames());
        assert!(res.final_energy().is_finite());
        assert!(
            res.final_temperature() > 0.0,
            "final T = {}",
            res.final_temperature()
        );
    }
}
