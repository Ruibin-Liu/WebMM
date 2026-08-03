//! Well-tempered metadynamics for enhanced conformational sampling.
//!
//! Deposits Gaussian bias hills along a collective variable (CV) during an MD
//! run, filling free-energy minima and forcing the system to explore the full
//! conformational landscape. The accumulated bias approximates the negative
//! free-energy surface (FES) along the CV.
//!
//! Architecture: [`MetaDynamics`] implements [`crate::forces::ForceField`],
//! wrapping the underlying potential (e.g. MMFF) and adding the bias.
//! [`crate::md::MDRunner`] drives it generically — hill deposition happens
//! inside the force evaluation via interior mutability (`RefCell`).

use crate::forces::ForceField;
use std::cell::RefCell;

// ── helpers ──

fn sub(a: &[f64; 3], b: &[f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
fn cross(a: &[f64; 3], b: &[f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn dot(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
fn norm(a: &[f64; 3]) -> f64 {
    dot(a, a).sqrt()
}
fn normalize(a: &[f64; 3]) -> [f64; 3] {
    let n = norm(a);
    if n < 1e-12 {
        return [0.0; 3];
    }
    [a[0] / n, a[1] / n, a[2] / n]
}

/// Dihedral angle (radians) p0–p1–p2–p3, range [−π, π].
pub fn dihedral_angle(p0: &[f64; 3], p1: &[f64; 3], p2: &[f64; 3], p3: &[f64; 3]) -> f64 {
    let b1 = sub(p0, p1); // Praxeolitic convention: b0 = p0 - p1 (negated first bond)
    let b2 = sub(p2, p1);
    let b3 = sub(p3, p2);
    let b2n = normalize(&b2);
    // v = b3 with its b2-component removed (lies in plane ⊥ b2)
    let v = [
        b3[0] - dot(&b3, &b2n) * b2n[0],
        b3[1] - dot(&b3, &b2n) * b2n[1],
        b3[2] - dot(&b3, &b2n) * b2n[2],
    ];
    // w = b1 with its b2-component removed
    let w = [
        b1[0] - dot(&b1, &b2n) * b2n[0],
        b1[1] - dot(&b1, &b2n) * b2n[1],
        b1[2] - dot(&b1, &b2n) * b2n[2],
    ];
    let x = dot(&v, &w);
    let y = dot(&cross(&b2n, &v), &w);
    y.atan2(x)
}

// ── collective variables ──

/// A scalar function of atomic coordinates, with its gradient.
pub trait CollectiveVariable: Send {
    /// Compute the CV value and d(s)/d(x_i) at the given coordinates.
    fn value_and_gradient(&self, coords: &[[f64; 3]]) -> (f64, Vec<[f64; 3]>);
    /// Human-readable name.
    fn name(&self) -> &str;
    /// Typical [min, max] range for FES plotting.
    fn range(&self) -> (f64, f64);
    /// Atom indices that influence this CV (for efficient FD gradient).
    fn atoms(&self) -> Vec<usize>;
    /// Compute the CV value only (no gradient).
    fn value(&self, coords: &[[f64; 3]]) -> f64 {
        self.value_and_gradient(coords).0
    }
}

/// Central finite-difference gradient for a CV that only depends on `atoms`.
fn fd_gradient<F: Fn(&[[f64; 3]]) -> f64>(
    cv_fn: F,
    coords: &[[f64; 3]],
    atoms: &[usize],
) -> Vec<[f64; 3]> {
    let mut grad = vec![[0.0; 3]; coords.len()];
    let h = 1e-6;
    let mut perturbed: Vec<[f64; 3]> = coords.to_vec();
    for &a in atoms {
        for c in 0..3 {
            perturbed[a][c] = coords[a][c] + h;
            let sp = cv_fn(&perturbed);
            perturbed[a][c] = coords[a][c] - h;
            let sm = cv_fn(&perturbed);
            perturbed[a][c] = coords[a][c]; // restore
                                            // handle periodic wraparound for angular CVs
            let mut d = sp - sm;
            if d > std::f64::consts::PI {
                d -= 2.0 * std::f64::consts::PI;
            } else if d < -std::f64::consts::PI {
                d += 2.0 * std::f64::consts::PI;
            }
            grad[a][c] = d / (2.0 * h);
        }
    }
    grad
}

/// A dihedral-angle CV (radians, [−π, π]). Biases rotation around the j–k bond.
pub struct DihedralCV {
    pub i: usize,
    pub j: usize,
    pub k: usize,
    pub l: usize,
    pub label: String,
}

impl DihedralCV {
    pub fn new(i: usize, j: usize, k: usize, l: usize) -> Self {
        DihedralCV {
            i,
            j,
            k,
            l,
            label: format!("Dihedral {}-{}-{}-{}", i, j, k, l),
        }
    }
}

impl CollectiveVariable for DihedralCV {
    fn value_and_gradient(&self, coords: &[[f64; 3]]) -> (f64, Vec<[f64; 3]>) {
        let s = dihedral_angle(
            &coords[self.i],
            &coords[self.j],
            &coords[self.k],
            &coords[self.l],
        );
        let (i, j, k, l) = (self.i, self.j, self.k, self.l);
        let grad = fd_gradient(
            move |c| dihedral_angle(&c[i], &c[j], &c[k], &c[l]),
            coords,
            &[i, j, k, l],
        );
        (s, grad)
    }
    fn name(&self) -> &str {
        &self.label
    }
    fn range(&self) -> (f64, f64) {
        (-std::f64::consts::PI, std::f64::consts::PI)
    }
    fn atoms(&self) -> Vec<usize> {
        vec![self.i, self.j, self.k, self.l]
    }
}

/// A distance CV (Å). Biases the distance between atoms i and j.
pub struct DistanceCV {
    pub i: usize,
    pub j: usize,
    pub label: String,
}

impl DistanceCV {
    pub fn new(i: usize, j: usize) -> Self {
        DistanceCV {
            i,
            j,
            label: format!("Distance {}-{}", i, j),
        }
    }
}

impl CollectiveVariable for DistanceCV {
    fn value_and_gradient(&self, coords: &[[f64; 3]]) -> (f64, Vec<[f64; 3]>) {
        let r = sub(&coords[self.i], &coords[self.j]);
        let d = norm(&r);
        let grad = fd_gradient(
            |c| norm(&sub(&c[self.i], &c[self.j])),
            coords,
            &[self.i, self.j],
        );
        (d, grad)
    }
    fn name(&self) -> &str {
        &self.label
    }
    fn range(&self) -> (f64, f64) {
        (1.0, 15.0)
    }
    fn atoms(&self) -> Vec<usize> {
        vec![self.i, self.j]
    }
}

// ── metadynamics ──

/// A deposited Gaussian hill in CV space.
#[derive(Clone, Copy)]
struct Hill {
    center: f64,
    height: f64,
}

/// Configuration for a metadynamics run.
#[derive(Clone, Copy)]
pub struct MetaDConfig {
    /// Initial hill height (kcal/mol). Typical: 0.1–1.0.
    pub hill_height: f64,
    /// Hill width σ in CV units (radians for dihedral, Å for distance). Typical: 0.1–0.35.
    pub hill_width: f64,
    /// Deposit a hill every N MD steps. Typical: 10–100 (≈ every 10–100 fs at 1 fs/step).
    pub deposit_interval: usize,
    /// Well-tempered bias factor γ = (T+ΔT)/T. 0 = standard (non-well-tempered). Typical: 5–20.
    pub bias_factor: f64,
    /// Temperature (K) — used for well-tempered hill rescaling.
    pub temperature_k: f64,
}

impl Default for MetaDConfig {
    fn default() -> Self {
        MetaDConfig {
            hill_height: 0.3,
            hill_width: 0.2,
            deposit_interval: 50,
            bias_factor: 10.0,
            temperature_k: 300.0,
        }
    }
}

/// Metadynamics bias potential. Wraps an underlying [`ForceField`] and adds
/// Gaussian bias hills along a [`CollectiveVariable`]. Implements `ForceField`
/// so it can be passed directly to [`crate::md::MDRunner`].
pub struct MetaDynamics<'a> {
    ff: &'a dyn ForceField,
    cv: Box<dyn CollectiveVariable>,
    hills: RefCell<Vec<Hill>>,
    config: MetaDConfig,
    call_count: RefCell<usize>,
    last_cv: RefCell<f64>,
}

impl<'a> MetaDynamics<'a> {
    pub fn new(
        ff: &'a dyn ForceField,
        cv: Box<dyn CollectiveVariable>,
        config: MetaDConfig,
    ) -> Self {
        MetaDynamics {
            ff,
            cv,
            hills: RefCell::new(Vec::new()),
            config,
            call_count: RefCell::new(0),
            last_cv: RefCell::new(0.0),
        }
    }

    /// The CV value from the last force evaluation.
    pub fn last_cv(&self) -> f64 {
        *self.last_cv.borrow()
    }

    /// The CV's [min, max] range (for FES grid plotting).
    pub fn cv_range(&self) -> (f64, f64) {
        self.cv.range()
    }

    /// Number of hills deposited so far.
    pub fn hill_count(&self) -> usize {
        self.hills.borrow().len()
    }

    /// Centers of all deposited hills (for visualization).
    pub fn hill_centers(&self) -> Vec<f64> {
        self.hills.borrow().iter().map(|h| h.center).collect()
    }

    /// Reconstruct the free-energy surface (FES) along the CV.
    /// Standard: F(s) = −V_bias(s). Well-tempered: F(s) = −γ/(γ−1) · V_bias(s).
    pub fn free_energy_surface(&self, s_values: &[f64]) -> Vec<f64> {
        let scale = if self.config.bias_factor > 1.0 {
            self.config.bias_factor / (self.config.bias_factor - 1.0)
        } else {
            1.0 // standard metadynamics
        };
        let sigma = self.config.hill_width;
        let inv_2sigma2 = 1.0 / (2.0 * sigma * sigma);
        s_values
            .iter()
            .map(|&s| {
                let v_bias: f64 = self
                    .hills
                    .borrow()
                    .iter()
                    .map(|h| {
                        let ds = s - h.center;
                        h.height * (-(ds * ds) * inv_2sigma2).exp()
                    })
                    .sum();
                -scale * v_bias
            })
            .collect()
    }
}

impl<'a> ForceField for MetaDynamics<'a> {
    fn energy_and_gradient(&self, coords: &[[f64; 3]], grad: &mut [[f64; 3]]) -> f64 {
        // 1. underlying force field (writes MMFF gradient into grad)
        let mut energy = self.ff.energy_and_gradient(coords, grad);

        // 2. CV value + gradient at current coordinates
        let (s, ds_dx) = self.cv.value_and_gradient(coords);
        *self.last_cv.borrow_mut() = s;

        // 3. bias from all deposited hills (energy + gradient via chain rule)
        let sigma = self.config.hill_width;
        let inv_2sigma2 = 1.0 / (2.0 * sigma * sigma);
        for hill in self.hills.borrow().iter() {
            let ds = s - hill.center;
            let gauss = (-(ds * ds) * inv_2sigma2).exp();
            energy += hill.height * gauss;
            // dV/ds = height * gauss * (−ds / σ²)
            let dv_ds = hill.height * gauss * (-ds * inv_2sigma2 * 2.0);
            // chain rule: dV/dx = dV/ds · ds/dx
            for i in 0..coords.len() {
                grad[i][0] += dv_ds * ds_dx[i][0];
                grad[i][1] += dv_ds * ds_dx[i][1];
                grad[i][2] += dv_ds * ds_dx[i][2];
            }
        }

        // 4. deposit a hill if it's time (skip the initialization call at step 0)
        let mut count = self.call_count.borrow_mut();
        *count += 1;
        if *count > 1 && (*count - 1).is_multiple_of(self.config.deposit_interval) {
            let height = if self.config.bias_factor > 1.0 {
                // well-tempered: w_i = w0 · exp(−V_bias(s_i) / (kB · ΔT))
                let v_bias_here: f64 = self
                    .hills
                    .borrow()
                    .iter()
                    .map(|h| {
                        let d = s - h.center;
                        h.height * (-(d * d) * inv_2sigma2).exp()
                    })
                    .sum();
                let delta_t = self.config.temperature_k * (self.config.bias_factor - 1.0);
                let kbt = crate::md::KB * delta_t;
                self.config.hill_height * (-v_bias_here / kbt).exp()
            } else {
                self.config.hill_height
            };
            self.hills.borrow_mut().push(Hill { center: s, height });
        }

        energy
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::forces::ForceField;
    use crate::md::{MDConfig, MDRunner};

    /// Test: dihedral angle for known geometries.
    #[test]
    fn dihedral_angle_known_values() {
        let pi = std::f64::consts::PI;
        // cis (0°): all 4 atoms in a line
        let cis = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [3.0, 0.0, 0.0],
        ];
        let d = dihedral_angle(&cis[0], &cis[1], &cis[2], &cis[3]);
        assert!(d.abs() < 0.01, "cis dihedral should be ~0, got {}", d);
        // trans (180°): p3 on the opposite side
        let trans = [
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, -1.0, 0.0],
        ];
        let d = dihedral_angle(&trans[0], &trans[1], &trans[2], &trans[3]);
        assert!(
            (d.abs() - pi).abs() < 0.01,
            "trans dihedral should be ~±π, got {}",
            d
        );
    }

    /// Test: CV gradient is finite and reasonable (no NaN/Inf, magnitude > 0).
    #[test]
    fn cv_gradient_is_finite() {
        let coords = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [2.0, 1.0, 0.0],
        ];
        let cv = DihedralCV::new(0, 1, 2, 3);
        let (s, grad) = cv.value_and_gradient(&coords);
        assert!(s.is_finite(), "CV value not finite: {}", s);
        for g in grad.iter() {
            assert!(g.iter().all(|&v| v.is_finite()), "CV gradient not finite");
        }
        // atoms 0-3 should have nonzero gradient; others zero
        let total: f64 = grad[0].iter().map(|v| v * v).sum::<f64>().sqrt()
            + grad[1].iter().map(|v| v * v).sum::<f64>().sqrt()
            + grad[2].iter().map(|v| v * v).sum::<f64>().sqrt()
            + grad[3].iter().map(|v| v * v).sum::<f64>().sqrt();
        assert!(total > 0.1, "CV gradient magnitude too small: {}", total);
    }

    /// Harmonic oscillator for testing metadynamics without MMFF.
    struct HO {
        k: f64,
    }
    impl ForceField for HO {
        fn energy_and_gradient(&self, coords: &[[f64; 3]], grad: &mut [[f64; 3]]) -> f64 {
            for g in grad.iter_mut() {
                *g = [0.0; 3];
            }
            let mut e = 0.0;
            for i in 0..coords.len() {
                for c in 0..3 {
                    e += 0.5 * self.k * coords[i][c] * coords[i][c];
                    grad[i][c] += self.k * coords[i][c];
                }
            }
            e
        }
    }

    /// Smoke test: run metadynamics on a harmonic oscillator with a distance CV.
    /// The FES should be non-trivial (the bias fills the harmonic well).
    #[test]
    fn metad_runs_and_produces_fes() {
        let ho = HO { k: 1.0 };
        let cv = Box::new(DistanceCV::new(0, 1));
        let metad = MetaDynamics::new(
            &ho,
            cv,
            MetaDConfig {
                hill_height: 0.2,
                hill_width: 0.3,
                deposit_interval: 10,
                bias_factor: 0.0, // standard
                temperature_k: 300.0,
            },
        );
        let coords = vec![[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]];
        let masses = vec![1.0, 1.0];
        let mut runner = MDRunner::new(
            &metad,
            masses,
            coords,
            MDConfig {
                dt_fs: 0.5,
                temperature_k: 300.0,
                friction_per_ps: 5.0,
                seed: 42,
            },
        );
        // run 2000 steps → ~200 hills deposited
        for _ in 0..2000 {
            runner.step();
        }
        assert!(
            metad.hill_count() > 50,
            "expected >50 hills, got {}",
            metad.hill_count()
        );
        // FES over the CV range
        let grid: Vec<f64> = (0..100).map(|i| 0.5 + i as f64 * 0.05).collect();
        let fes = metad.free_energy_surface(&grid);
        let fes_range = fes.iter().cloned().fold(f64::INFINITY, f64::min)
            ..=fes.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let span = fes_range.end() - fes_range.start();
        assert!(
            span > 0.5,
            "FES should have non-trivial range, got span {}",
            span
        );
    }

    /// Test: metadynamics with MMFF on ethanol — dihedral CV, completes & explores.
    #[test]
    fn metad_on_ethanol_completes() {
        let sdf = "ethanol\n     RDKit          3D\n\n  9  8  0  0  0  0  0  0  0  0999 V2000\n    0.8867    0.1097    0.1398 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.4649   -0.1632   -0.4871 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.4908    0.4341    0.2924 O   0  0  0  0  0  0  0  0  0  0  0  0\n    1.0600    1.1872    0.2278 H   0  0  0  0  0  0  0  0  0  0  0  0\n    1.6897   -0.3301   -0.4581 H   0  0  0  0  0  0  0  0  0  0  0  0\n    0.9343   -0.3042    1.1524 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.5119    0.2599   -1.4947 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.6562   -1.2384   -0.5550 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.4469    0.0448    1.1825 H   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0\n  2  3  1  0\n  1  4  1  0\n  1  5  1  0\n  1  6  1  0\n  2  7  1  0\n  2  8  1  0\n  3  9  1  0\nM  END\n";
        let mol = crate::molecule::parser::parse_sdf(sdf).expect("parse");
        let ff = crate::mmff::MMFFForceField::new(&mol, crate::MMFFVariant::MMFF94s);
        // CV: H(O)–O–C–C dihedral (atoms 8, 2, 1, 0 in 0-indexed)
        let cv = Box::new(DihedralCV::new(8, 2, 1, 0));
        let metad = MetaDynamics::new(
            &ff,
            cv,
            MetaDConfig {
                hill_height: 0.3,
                hill_width: 0.2,
                deposit_interval: 20,
                bias_factor: 8.0,
                temperature_k: 300.0,
            },
        );
        let mut runner = MDRunner::from_molecule(
            &metad,
            &mol,
            MDConfig {
                dt_fs: 1.0,
                temperature_k: 300.0,
                friction_per_ps: 5.0,
                seed: 7,
            },
        );
        for _ in 0..3000 {
            runner.step();
        }
        let n_hills = metad.hill_count();
        assert!(n_hills > 50, "expected >50 hills, got {}", n_hills);
        // FES
        let pi = std::f64::consts::PI;
        let grid: Vec<f64> = (0..72)
            .map(|i| -pi + i as f64 * (2.0 * pi / 71.0))
            .collect();
        let fes = metad.free_energy_surface(&grid);
        let span = fes.iter().cloned().fold(f64::NEG_INFINITY, f64::max)
            - fes.iter().cloned().fold(f64::INFINITY, f64::min);
        assert!(
            span.is_finite() && span > 0.5,
            "FES span should be >0.5, got {}",
            span
        );
    }
}
