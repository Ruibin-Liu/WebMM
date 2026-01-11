//! L-BFGS optimization algorithm

use crate::mmff::MMFFForceField;
use crate::ConvergenceOptions;

/// Optimization result
pub struct OptimizationResult {
    pub optimized_coords: Vec<[f64; 3]>,
    pub final_energy: f64,
    pub converged: bool,
    pub iterations: usize,
}

/// L-BFGS optimizer
pub fn optimize(
    ff: &MMFFForceField,
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

    for iter in 0..convergence.max_iterations {
        // Calculate energy and gradient
        let coords_2d = flatten_to_2d(&x, n_atoms);
        let (energy, g_2d) = ff.calculate_energy_and_gradient(&coords_2d);

        // Flatten gradient to 1D
        let mut g = Vec::with_capacity(n_coords);
        for grad in g_2d.iter() {
            g.push(grad[0]);
            g.push(grad[1]);
            g.push(grad[2]);
        }

        final_energy = energy;

        // Check convergence
        let max_f = g
            .iter()
            .map(|&gi| gi.abs())
            .fold(0.0f64, |a, b| a.max(b.abs()));
        let rms_f = (g.iter().map(|gi| gi * gi).sum::<f64>() / n_coords as f64).sqrt();

        if max_f < convergence.max_force && rms_f < convergence.rms_force {
            converged = true;
            break;
        }

        // Compute search direction
        let d = if iter == 0 || s_history.is_empty() {
            g.iter().map(|&gi| -gi).collect()
        } else {
            compute_lbfgs_direction(&g, &s_history, &y_history, &rho_history)
        };

        // Line search (Armijo backtracking)
        let alpha = armijo_line_search(ff, &x, &d, n_atoms, energy, &g, 0.01, 0.1, 1e-4, 1e-10);

        // Update x
        for i in 0..n_coords {
            x[i] += alpha * d[i];
        }

        // Update history
        let g_new_2d = flatten_to_2d(&x, n_atoms);
        let (_, g_new_2d) = ff.calculate_energy_and_gradient(&g_new_2d);

        let mut g_new = Vec::with_capacity(n_coords);
        for grad in g_new_2d.iter() {
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

        if y_dot_s.abs() > 1e-10 {
            if s_history.len() >= memory_size {
                s_history.remove(0);
                y_history.remove(0);
                rho_history.remove(0);
            }

            s_history.push(x_diff.clone());
            y_history.push(g_diff.clone());
            rho_history.push(1.0 / y_dot_s);
        }
    }

    let optimized_coords = flatten_to_2d(&x, n_atoms);

    OptimizationResult {
        optimized_coords,
        final_energy,
        converged,
        iterations: if converged {
            converged as usize
        } else {
            convergence.max_iterations
        },
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

    let mut q = g.to_vec();

    // Two-loop recursion
    for i in (0..m).rev() {
        let alpha = rho_history[i]
            * y_history[i]
                .iter()
                .zip(s_history[i].iter())
                .map(|(yi, sj)| yi * sj)
                .sum::<f64>();

        let s_g = s_history[i]
            .iter()
            .zip(g.iter())
            .map(|(si, gj)| si * gj)
            .sum::<f64>();

        for j in 0..q.len() {
            q[j] -= alpha * s_g;
        }
    }

    q
}

/// Armijo line search
fn armijo_line_search(
    ff: &MMFFForceField,
    x: &[f64],
    d: &[f64],
    n_atoms: usize,
    f0: f64,
    g0: &[f64],
    alpha0: f64,
    rho: f64,
    c1: f64,
    min_alpha: f64,
) -> f64 {
    let mut alpha = alpha0;
    let n_coords = n_atoms * 3;

    loop {
        // Compute x_new = x + alpha * d
        let mut x_new = x.to_vec();
        for i in 0..n_coords {
            x_new[i] += alpha * d[i];
        }

        // Calculate new energy
        let coords_2d = flatten_to_2d(&x_new, n_atoms);
        let f_new = ff.calculate_energy(&coords_2d);

        // Compute g^T * d
        let g_dot_d = g0.iter().zip(d.iter()).map(|(gi, di)| gi * di).sum::<f64>();

        // Armijo condition: f(x + alpha*d) <= f(x) + c1 * alpha * g(x)^T * d
        let rhs = f0 + c1 * alpha * g_dot_d;

        if f_new <= rhs || alpha <= min_alpha {
            break;
        }

        // Backtrack
        alpha *= rho;
    }

    alpha
}
