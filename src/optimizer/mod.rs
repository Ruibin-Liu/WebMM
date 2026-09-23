//! L-BFGS optimization algorithm

pub mod internal_opt;
pub mod internals;
pub mod jacobi;

use crate::forces::ForceField;
use crate::ConvergenceOptions;

/// Env-gated optimizer trace (OPT_DEBUG=1): per-iteration direction/slope and
/// per-trial Armijo numbers to stderr. Cached — no per-call env lookup cost.
fn opt_debug() -> bool {
    static DBG: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *DBG.get_or_init(|| std::env::var("OPT_DEBUG").is_ok())
}

/// Optimization result
pub struct OptimizationResult {
    pub optimized_coords: Vec<[f64; 3]>,
    pub final_energy: f64,
    /// Force criterion met (max_force AND rms_force below thresholds).
    pub converged: bool,
    /// Energy-resolution floor reached: line searches are failing while the
    /// energy of the last K accepted steps is stationary to <1e-9 kcal/mol —
    /// the surface cannot be resolved further in f64 (GFN-FF large molecules;
    /// real residual forces live in soft torsional modes). Callers may report
    /// this as a converged stop with an explanatory message.
    pub energy_converged: bool,
    pub iterations: usize,
}

/// Optimization variable space for the L-BFGS core.
///
/// The core loop is written once against this trait; the classic Cartesian
/// optimizer wraps a ForceField (z = flattened coordinates, bit-identical to
/// the pre-trait implementation), and the internal-coordinate optimizer maps
/// delocalized internals back and forth to Cartesian for every evaluation.
pub trait Objective {
    fn dim(&self) -> usize;
    /// Energy and dE/dz at `z` (the new accepted point). Also refreshes the
    /// Cartesian force statistics.
    fn f_and_g(&mut self, z: &[f64]) -> (f64, Vec<f64>);
    /// Energy at a trial point; must not advance any internal state.
    fn energy(&mut self, z: &[f64]) -> f64;
    /// (max |g_x|, rms |g_x|) of the last f_and_g — the physical force
    /// criterion (kcal/mol/A), independent of the variable space.
    fn force_stats(&self) -> (f64, f64);
    /// If the last f_and_g rebuilt its coordinate system, consume the flag:
    /// the caller must reset z to zeros (dim() changed) and clear L-BFGS
    /// history.
    fn take_reset(&mut self) -> bool;
    /// Cartesian geometry of the last accepted point.
    fn final_coords(&mut self) -> Vec<[f64; 3]>;
}

/// z = flattened Cartesian coordinates over a ForceField.
struct CartesianObjective<'a> {
    ff: &'a dyn ForceField,
    n_atoms: usize,
    last_z: Vec<f64>,
    last_max: f64,
    last_rms: f64,
}

impl<'a> CartesianObjective<'a> {
    fn new(ff: &'a dyn ForceField, n_atoms: usize) -> Self {
        CartesianObjective {
            ff,
            n_atoms,
            last_z: Vec::new(),
            last_max: 0.0,
            last_rms: 0.0,
        }
    }
}

impl Objective for CartesianObjective<'_> {
    fn dim(&self) -> usize {
        3 * self.n_atoms
    }
    fn f_and_g(&mut self, z: &[f64]) -> (f64, Vec<f64>) {
        let coords_2d = flatten_to_2d(z, self.n_atoms);
        let mut g_2d = vec![[0.0f64; 3]; self.n_atoms];
        let e = self.ff.energy_and_gradient(&coords_2d, &mut g_2d);
        self.last_z = z.to_vec();
        let g: Vec<f64> = g_2d
            .iter()
            .flat_map(|grad| [grad[0], grad[1], grad[2]])
            .collect();
        let nc = g.len();
        self.last_max = g
            .iter()
            .map(|gi| gi.abs())
            .fold(0.0f64, |a, b| a.max(b.abs()));
        self.last_rms = (g.iter().map(|gi| gi * gi).sum::<f64>() / nc as f64).sqrt();
        (e, g)
    }
    fn energy(&mut self, z: &[f64]) -> f64 {
        let coords_2d = flatten_to_2d(z, self.n_atoms);
        self.ff.energy(&coords_2d)
    }
    fn force_stats(&self) -> (f64, f64) {
        (self.last_max, self.last_rms)
    }
    fn take_reset(&mut self) -> bool {
        false
    }
    fn final_coords(&mut self) -> Vec<[f64; 3]> {
        flatten_to_2d(&self.last_z, self.n_atoms)
    }
}

/// L-BFGS optimizer over Cartesian coordinates (the classic entry point).
pub fn optimize(
    ff: &dyn ForceField,
    initial_coords: &[[f64; 3]],
    convergence: &ConvergenceOptions,
) -> OptimizationResult {
    let mut z = Vec::with_capacity(3 * initial_coords.len());
    for coord in initial_coords.iter() {
        z.push(coord[0]);
        z.push(coord[1]);
        z.push(coord[2]);
    }
    let mut obj = CartesianObjective::new(ff, initial_coords.len());
    let r = lbfgs_core(&mut obj, z, convergence);
    // Cartesian z IS the coordinate vector
    r
}

/// The variable-space-agnostic L-BFGS loop (line search: energy-only Armijo
/// trials with quadratic-interpolation backtracking, unit initial step for
/// L-BFGS directions). `initial_z` is the starting point in the objective's
/// variable space.
pub fn lbfgs_core(
    obj: &mut dyn Objective,
    initial_z: Vec<f64>,
    convergence: &ConvergenceOptions,
) -> OptimizationResult {
    let mut z = initial_z;
    let n_coords = z.len();

    // L-BFGS history
    let memory_size = 20;
    let mut s_history = Vec::with_capacity(memory_size);
    let mut y_history = Vec::with_capacity(memory_size);
    let mut rho_history = Vec::with_capacity(memory_size);

    let mut converged = false;
    let mut energy_converged = false;
    let mut final_energy = 0.0;
    let mut final_iter = 0;
    let mut fail_count = 0usize;
    // Energy-resolution floor detection: cumulative line-search failures plus
    // the energy trace of the last 11 accepted points (10 steps).
    let mut total_failures = 0usize;
    let mut energy_trace: Vec<f64> = Vec::new();
    let mut tiny_step_count = 0usize;

    // E+G once per iteration: the gradient evaluated at the updated point is
    // carried into the next iteration instead of being discarded and
    // recomputed at the top of the loop (halves the force-field calls).
    let (mut energy, mut g) = obj.f_and_g(&z);

    for iter in 0..convergence.max_iterations {
        // Check force convergence (physical, Cartesian — objective-supplied)
        let (max_f, rms_f) = obj.force_stats();
        if max_f < convergence.max_force && rms_f < convergence.rms_force {
            converged = true;
            final_energy = energy;
            final_iter = iter;
            break;
        }

        // Compute search direction
        let d = if iter == 0 || s_history.is_empty() {
            g.iter().map(|&gi| -gi).collect()
        } else {
            compute_lbfgs_direction(&g, &s_history, &y_history, &rho_history)
        };

        // Verify descent direction: g^T * d must be negative
        let gt_dot_d: f64 = g.iter().zip(d.iter()).map(|(gi, di)| gi * di).sum();
        let d = if gt_dot_d >= 0.0 {
            // L-BFGS produced a non-descent direction; reset to steepest descent
            s_history.clear();
            y_history.clear();
            rho_history.clear();
            g.iter().map(|&gi| -gi).collect()
        } else {
            d
        };

        // Line search (Armijo backtracking with quadratic interpolation)
        // L-BFGS directions from the two-loop recursion carry the inverse-
        // curvature scaling, so the standard unit initial step applies
        // (Nocedal-Wright). Only raw steepest-descent steps need explicit
        // displacement scaling (~0.5 A-equivalent). Safeguard: when the
        // direction norm is pathological, fall back to a 10 A-equivalent cap
        // so backtracking starts in physical territory.
        let max_component = d
            .iter()
            .map(|di| di.abs())
            .fold(0.0f64, f64::max)
            .max(1e-10);
        let initial_alpha = if s_history.len() < 3 {
            0.5 / max_component // ~0.5 A-equivalent displacement for steepest descent
        } else {
            1.0f64.min(10.0 / max_component) // unit step, displacement cap
        };
        let slope: f64 = g.iter().zip(d.iter()).map(|(gi, di)| gi * di).sum();
        if opt_debug() {
            eprintln!(
                "iter {iter}: hist={} maxd={max_component:.3e} a0={initial_alpha:.3e} slope={slope:.3e} maxf={max_f:.3e}",
                s_history.len()
            );
        }
        let (alpha, _floor_accept) =
            armijo_line_search(obj, &z, &d, energy, slope, initial_alpha, 1e-10);

        // Line-search failure (no Armijo decrease at the floor): reset the
        // L-BFGS memory and retry from steepest descent; after several
        // consecutive failures the surface is numerically flat/broken here —
        // stop rather than wander.
        if alpha == 0.0 {
            fail_count += 1;
            total_failures += 1;
            s_history.clear();
            y_history.clear();
            rho_history.clear();
            if fail_count >= 5 {
                final_energy = energy;
                final_iter = iter;
                // 5 consecutive failures with a stationary recent trace is the
                // same f64 floor, not a broken surface
                if energy_trace.len() >= 6
                    && energy_trace
                        .windows(2)
                        .rev()
                        .take(5)
                        .map(|w| (w[1] - w[0]).abs())
                        .fold(0.0f64, f64::max)
                        < 1e-8
                {
                    energy_converged = true;
                }
                break;
            }
            continue;
        }
        fail_count = 0;

        // Tiny accepted step: reset L-BFGS history (corrupted approximation)
        // and count it — steps below ~1e-8 carry no resolvable energy change;
        // 3 in a row is the f64 floor (GFN-FF large-molecule terminal pattern).
        if alpha * max_component < 1e-8 {
            s_history.clear();
            y_history.clear();
            rho_history.clear();
            tiny_step_count += 1;
        } else {
            tiny_step_count = 0;
        }

        // Update z
        for i in 0..n_coords {
            z[i] += alpha * d[i];
        }

        // E+G at the updated point — reused as the next iteration's input
        // (keeps the invariant final_energy == E(accepted point) even when
        // exiting via max_iterations)
        let g_new: Vec<f64>;
        {
            let (e, gn) = obj.f_and_g(&z);
            energy = e;
            g_new = gn;
        }
        final_energy = energy;

        // Coordinate-system rebuild (internal-coordinate objective): adopt
        // the new basis, restart z at the build point.
        if obj.take_reset() {
            z = vec![0.0f64; obj.dim()];
            s_history.clear();
            y_history.clear();
            rho_history.clear();
            // the returned gradient is already in the new basis at z = 0
            g = g_new;
            continue; // history empty; next iteration re-derives direction
        }

        // Energy-resolution floor: line searches keep failing while accepted
        // steps stopped moving the energy (<1e-9 kcal/mol over the last 10) —
        // the f64 limit of the surface.
        energy_trace.push(energy);
        if energy_trace.len() > 11 {
            energy_trace.remove(0);
        }
        if total_failures >= 3 && energy_trace.len() == 11 {
            let stationary = energy_trace
                .windows(2)
                .map(|w| (w[1] - w[0]).abs())
                .fold(0.0f64, f64::max)
                < 1e-9;
            if stationary {
                energy_converged = true;
                final_iter = iter + 1;
                break;
            }
        }
        if tiny_step_count >= 3 {
            energy_converged = true;
            final_iter = iter + 1;
            break;
        }

        let g_diff: Vec<f64> = g_new.iter().zip(g.iter()).map(|(gn, go)| gn - go).collect();
        let x_diff: Vec<f64> = d.iter().map(|di| alpha * di).collect();

        let y_dot_s = g_diff
            .iter()
            .zip(x_diff.iter())
            .map(|(gi, xi)| gi * xi)
            .sum::<f64>();

        if y_dot_s > 1e-10 {
            if s_history.len() >= memory_size {
                s_history.remove(0);
                y_history.remove(0);
                rho_history.remove(0);
            }
            s_history.push(x_diff);
            y_history.push(g_diff);
            rho_history.push(1.0 / y_dot_s);
        }

        // carry the fresh gradient into the next iteration
        g = g_new;

        final_iter = iter + 1;
    }

    OptimizationResult {
        optimized_coords: obj.final_coords(),
        final_energy,
        converged,
        energy_converged,
        iterations: final_iter,
    }
}

fn flatten_to_2d(x: &[f64], n_atoms: usize) -> Vec<[f64; 3]> {
    let mut coords_2d = Vec::with_capacity(n_atoms);
    for i in 0..n_atoms {
        coords_2d.push([x[i * 3], x[i * 3 + 1], x[i * 3 + 2]]);
    }
    coords_2d
}

/// Compute L-BFGS search direction using two-loop recursion
fn compute_lbfgs_direction(
    g: &[f64],
    s_history: &[Vec<f64>],
    y_history: &[Vec<f64>],
    rho_history: &[f64],
) -> Vec<f64> {
    let m = s_history.len();
    if m == 0 {
        return g.iter().map(|&gi| -gi).collect();
    }

    let n = g.len();
    let mut q = g.to_vec();
    let mut alpha_arr = vec![0.0; m];

    // First loop (backward)
    for i in (0..m).rev() {
        let s_dot_q: f64 = s_history[i]
            .iter()
            .zip(q.iter())
            .map(|(s, qi)| s * qi)
            .sum();
        alpha_arr[i] = rho_history[i] * s_dot_q;
        for j in 0..n {
            q[j] -= alpha_arr[i] * y_history[i][j];
        }
    }

    // H0 scaling
    let last = m - 1;
    let y_y: f64 = y_history[last].iter().map(|y| y * y).sum();
    let s_y_last: f64 = s_history[last]
        .iter()
        .zip(y_history[last].iter())
        .map(|(s, y)| s * y)
        .sum();
    let gamma = if y_y > 1e-20 { s_y_last / y_y } else { 1.0 };

    let mut r = q.iter().map(|qi| gamma * qi).collect::<Vec<f64>>();

    // Second loop (forward)
    for i in 0..m {
        let y_dot_r: f64 = y_history[i]
            .iter()
            .zip(r.iter())
            .map(|(y, ri)| y * ri)
            .sum();
        let beta = rho_history[i] * y_dot_r;
        for j in 0..n {
            r[j] += (alpha_arr[i] - beta) * s_history[i][j];
        }
    }

    r.iter().map(|ri| -ri).collect()
}

/// Armijo line search with quadratic-interpolation backtracking over the
/// objective's variable space.
///
/// Trial points are evaluated with `obj.energy` only (no gradient); `slope`
/// is g(z)·d (< 0). Backtracking shrinks alpha by the minimizer of the
/// one-dimensional quadratic through (0, f0) with slope `slope` and the last
/// trial, safeguarded to [0.1, 0.5]×alpha (plain halving on degenerate
/// interpolation). The Armijo constant is c1 = 1e-4.
fn armijo_line_search(
    obj: &mut dyn Objective,
    z: &[f64],
    d: &[f64],
    f0: f64,
    slope: f64,
    alpha0: f64,
    min_alpha: f64,
) -> (f64, bool) {
    const C1: f64 = 1e-4;
    let mut alpha = alpha0;
    let n = z.len();

    loop {
        // Compute z_new = z + alpha * d
        let mut z_new = z.to_vec();
        for i in 0..n {
            z_new[i] += alpha * d[i];
        }

        // Energy at the trial point (gradient not needed for Armijo)
        let f_new = obj.energy(&z_new);
        let rhs = f0 + C1 * alpha * slope;
        if opt_debug() {
            eprintln!(
                "   trial a={alpha:.3e} f={f_new:+.10} (f0={f0:+.10} rhs={rhs:+.10} ok={})",
                f_new <= rhs
            );
        }

        if f_new <= rhs {
            return (alpha, false);
        }
        if alpha <= min_alpha {
            // Accept the floor step ONLY if it is finite and does not increase
            // the energy (the old unconditional accept once jumped caffeine
            // GFN-FF +248 kcal/mol into an inescapable basin).
            if f_new.is_finite() && f_new <= f0 {
                return (alpha, true); // non-increasing: keep the slither-through
            }
            return (0.0, false);
        }

        // Quadratic interpolation: minimize the parabola through (0, f0) with
        // slope `slope` and (alpha, f_new).
        let denom = 2.0 * (f0 + slope * alpha - f_new);
        let a_interp = if denom > 1e-300 {
            slope * alpha * alpha / denom
        } else {
            0.5 * alpha
        };
        alpha = a_interp.clamp(0.1 * alpha, 0.5 * alpha);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::forces::ForceField;
    use crate::mmff::MMFFForceField;
    use crate::molecule::parser::parse_sdf;
    use crate::MMFFVariant;
    use std::cell::Cell;

    /// Regression for the line-search rework (unit initial step for L-BFGS +
    /// quadratic-interpolation backtracking + energy-only trials): force-field
    /// evaluations per iteration must stay low. Before the rework the fixed
    /// 1.5 A initial displacement overshot ~1000x and burned 10-13 E+G calls
    /// per iteration (25+ on GFN-FF pre-gradient-fix).
    struct Counting<'a> {
        inner: &'a MMFFForceField,
        eg: Cell<usize>,
        e: Cell<usize>,
    }
    impl ForceField for Counting<'_> {
        fn energy_and_gradient(&self, coords: &[[f64; 3]], grad: &mut [[f64; 3]]) -> f64 {
            self.eg.set(self.eg.get() + 1);
            self.inner.energy_and_gradient(coords, grad)
        }
        fn energy(&self, coords: &[[f64; 3]]) -> f64 {
            self.e.set(self.e.get() + 1);
            self.inner.energy(coords)
        }
    }

    #[test]
    fn line_search_call_budget() {
        let sdf = std::fs::read_to_string(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/tests/fixtures/conformers/ethanol.sdf"
        ))
        .unwrap();
        let mol = parse_sdf(&sdf).unwrap();
        let coords: Vec<[f64; 3]> = mol.atoms.iter().map(|a| a.position).collect();
        let ff = MMFFForceField::new(&mol, MMFFVariant::MMFF94s);
        let counting = Counting {
            inner: &ff,
            eg: Cell::new(0),
            e: Cell::new(0),
        };
        let r = optimize(&counting, &coords, &ConvergenceOptions::default());
        assert!(r.converged, "did not converge");
        let calls = counting.eg.get() + counting.e.get();
        let per_iter = calls as f64 / r.iterations.max(1) as f64;
        assert!(
            per_iter <= 5.0,
            "line search too expensive: {calls} calls over {} iters = {per_iter:.2}/iter",
            r.iterations
        );
    }
}
