//! L-BFGS optimization algorithm

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
    pub converged: bool,
    pub iterations: usize,
}

/// L-BFGS optimizer
pub fn optimize(
    ff: &dyn ForceField,
    initial_coords: &[[f64; 3]],
    convergence: &ConvergenceOptions,
) -> OptimizationResult {
    let n_atoms = initial_coords.len();
    let n_coords = n_atoms * 3;

    // Flatten to 1D: [x0, y0, z0, x1, y1, z1, ...]
    let mut x = Vec::with_capacity(n_coords);
    for coord in initial_coords.iter() {
        x.push(coord[0]);
        x.push(coord[1]);
        x.push(coord[2]);
    }

    // L-BFGS history
    let memory_size = 20;
    let mut s_history = Vec::with_capacity(memory_size);
    let mut y_history = Vec::with_capacity(memory_size);
    let mut rho_history = Vec::with_capacity(memory_size);

    let mut converged = false;
    let mut final_energy = 0.0;
    let mut final_iter = 0;
    let mut fail_count = 0usize;

    // E+G once per iteration: the gradient evaluated at the updated point is
    // carried into the next iteration instead of being discarded and
    // recomputed at the top of the loop (halves the force-field calls).
    let mut coords_2d = flatten_to_2d(&x, n_atoms);
    let mut g_2d = vec![[0.0f64; 3]; n_atoms];
    let mut energy = ff.energy_and_gradient(&coords_2d, &mut g_2d);
    let mut g: Vec<f64> = g_2d
        .iter()
        .flat_map(|grad| [grad[0], grad[1], grad[2]])
        .collect();

    for iter in 0..convergence.max_iterations {
        // Check force convergence
        let max_f = g
            .iter()
            .map(|&gi| gi.abs())
            .fold(0.0f64, |a, b| a.max(b.abs()));
        let rms_f = (g.iter().map(|gi| gi * gi).sum::<f64>() / n_coords as f64).sqrt();

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
        // displacement scaling (~0.5 A). Safeguard: when the direction norm
        // is pathological (force-spike curvature pairs can inflate |d| to
        // ~1e12, making the unit step a ~1e12 A displacement), fall back to
        // a 10 A displacement cap so backtracking starts in physical territory
        // (the old fixed 1.5 A cap for ALL L-BFGS steps overshot ~1000x on
        // normal directions and burned 10-26 E+G evaluations per iteration).
        let max_component = d
            .iter()
            .map(|di| di.abs())
            .fold(0.0f64, f64::max)
            .max(1e-10);
        let initial_alpha = if s_history.len() < 3 {
            0.5 / max_component // ~0.5 A max displacement for steepest descent
        } else {
            1.0f64.min(10.0 / max_component) // unit step, 10 A displacement cap
        };
        let slope: f64 = g.iter().zip(d.iter()).map(|(gi, di)| gi * di).sum();
        if opt_debug() {
            eprintln!(
                "iter {iter}: hist={} maxd={max_component:.3e} a0={initial_alpha:.3e} slope={slope:.3e} maxf={:.3e}",
                s_history.len(),
                g.iter().map(|v| v.abs()).fold(0.0, f64::max)
            );
        }
        let alpha = armijo_line_search(ff, &x, &d, n_atoms, energy, slope, initial_alpha, 1e-10);

        // Line-search failure (no Armijo decrease at the floor): reset the
        // L-BFGS memory and retry from steepest descent; after several
        // consecutive failures the surface is numerically flat/broken here —
        // stop rather than wander.
        if alpha == 0.0 {
            fail_count += 1;
            s_history.clear();
            y_history.clear();
            rho_history.clear();
            if fail_count >= 5 {
                final_energy = energy;
                final_iter = iter;
                break;
            }
            continue;
        }
        fail_count = 0;

        // If step is tiny, reset L-BFGS history (corrupted approximation)
        if alpha * max_component < 1e-8 {
            s_history.clear();
            y_history.clear();
            rho_history.clear();
        }

        // Update x
        for i in 0..n_coords {
            x[i] += alpha * d[i];
        }

        // E+G at the updated point — reused as the next iteration's input
        // (keeps the invariant final_energy == E(optimized_coords) even when
        // exiting via max_iterations)
        coords_2d = flatten_to_2d(&x, n_atoms);
        energy = ff.energy_and_gradient(&coords_2d, &mut g_2d);
        final_energy = energy;

        let mut g_new = Vec::with_capacity(n_coords);
        for grad in g_2d.iter() {
            g_new.push(grad[0]);
            g_new.push(grad[1]);
            g_new.push(grad[2]);
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

            s_history.push(x_diff.clone());
            y_history.push(g_diff.clone());
            rho_history.push(1.0 / y_dot_s);
        }

        // carry the fresh gradient into the next iteration (this is what makes
        // the restructure equivalent to the old double-evaluation loop)
        g = g_new;

        final_iter = iter + 1;
    }

    let optimized_coords = flatten_to_2d(&x, n_atoms);

    OptimizationResult {
        optimized_coords,
        final_energy,
        converged,
        iterations: final_iter,
    }
}

/// Flatten 1D coordinates to 2D array
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

/// Armijo line search with quadratic-interpolation backtracking.
///
/// Trial points are evaluated with `ff.energy` only (no gradient); `slope`
/// is g(x)·d (< 0). Backtracking shrinks alpha by the minimizer of the
/// one-dimensional quadratic through (0, f0) with slope `slope` and the
/// last trial, safeguarded to [0.1, 0.5]×alpha (plain halving on
/// degenerate interpolation). The Armijo constant is c1 = 1e-4.
#[allow(clippy::too_many_arguments)]
fn armijo_line_search(
    ff: &dyn ForceField,
    x: &[f64],
    d: &[f64],
    n_atoms: usize,
    f0: f64,
    slope: f64,
    alpha0: f64,
    min_alpha: f64,
) -> f64 {
    const C1: f64 = 1e-4;
    let mut alpha = alpha0;
    let n_coords = n_atoms * 3;

    loop {
        // Compute x_new = x + alpha * d
        let mut x_new = x.to_vec();
        for i in 0..n_coords {
            x_new[i] += alpha * d[i];
        }

        // Energy at the trial point (gradient not needed for Armijo)
        let coords_2d = flatten_to_2d(&x_new, n_atoms);
        let f_new = ff.energy(&coords_2d);

        // Armijo condition: f(x + alpha*d) <= f(x) + c1 * alpha * g(x)^T * d
        let rhs = f0 + C1 * alpha * slope;
        if opt_debug() {
            eprintln!(
                "   trial a={alpha:.3e} f={f_new:+.10} (f0={f0:+.10} rhs={rhs:+.10} ok={})",
                f_new <= rhs
            );
        }

        if f_new <= rhs {
            return alpha;
        }
        if alpha <= min_alpha {
            // The Armijo condition never held. Accept the floor step ONLY if
            // it is finite and strictly decreases the energy — the old
            // behaviour accepted it unconditionally, and with an exploded
            // L-BFGS direction (|d| up to 1e12 after a force spike) the
            // floor step is still macroscopic (~0.1-100 A): caffeine GFN-FF
            // once jumped +248 kcal/mol into a basin it could never leave.
            if f_new.is_finite() && f_new <= f0 {
                return alpha; // non-increasing: keep the old slither-through
            }
            return 0.0;
        }

        // Quadratic interpolation: minimize the parabola through
        // (0, f0) with slope `slope` and (alpha, f_new):
        // alpha* = slope·alpha² / (2·(f0 + slope·alpha − f_new))
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
