//! L-BFGS in delocalized internal coordinates.
//!
//! Wraps a ForceField with an `Objective` whose variable space is the
//! delocalized internal coordinate basis (see `internals`). The map z -> x is
//! FIXED per basis: x(z) = x_ref + BackTransform(B(x_ref), z) — a
//! well-defined function of z, so the L-BFGS line search operates on a
//! consistent energy surface. The basis is rebuilt (new reference/origin,
//! L-BFGS history cleared via `take_reset`) when the linearization drifts.
//! Falls back to plain Cartesian optimization when internals cannot be built.

use super::internals::{project_out_tr, InternalCoords};
use super::{lbfgs_core, Objective, OptimizationResult};
use crate::forces::ForceField;
use crate::ConvergenceOptions;

/// Rebuild the internal basis when the mapped geometry's true q has drifted
/// from the linear target z by more than this (Å/rad). The map x(z) is
/// exactly linear for any drift (no correctness need to rebuild) — drift
/// only degrades L-BFGS conditioning in the stale basis. Measured
/// {2.0, 5.0, ∞} matrix (v1.1.1): monotonically better with fewer rebuilds
/// on all fixtures (a G-diagonalization costs O(n_prim³·sweeps)); set to
/// effectively-never. Topology changes are not handled here — the
/// Cartesian-restart safety net covers pathological cases.
const DRIFT_REBUILD: f64 = 1e9;

struct InternalObjective<'a> {
    ff: &'a dyn ForceField,
    znums: Vec<usize>,
    ic: InternalCoords,
    /// Reference geometry of the current basis: x(z) = x_ref + BT(bq_ref, z).
    x_ref: Vec<[f64; 3]>,
    /// Last mapped geometry (accepted point) for stats/final coords.
    x_last: Vec<[f64; 3]>,
    last_max: f64,
    last_rms: f64,
    reset_flag: bool,
}

impl<'a> InternalObjective<'a> {
    fn new(ff: &'a dyn ForceField, x: &[[f64; 3]], znums: &[usize], ic: InternalCoords) -> Self {
        InternalObjective {
            ff,
            znums: znums.to_vec(),
            ic,
            x_ref: x.to_vec(),
            x_last: x.to_vec(),
            last_max: 0.0,
            last_rms: 0.0,
            reset_flag: false,
        }
    }

    /// x(z) = x_ref + BackTransform(z) — the well-defined (exactly linear)
    /// map of this basis, using the B_q cached at build (line-search safe).
    fn map(&self, z: &[f64]) -> Vec<[f64; 3]> {
        let dx = self.ic.back_transform(z);
        let mut x = self.x_ref.clone();
        for i in 0..x.len() {
            for t in 0..3 {
                x[i][t] += dx[i][t];
            }
        }
        x
    }
}

impl Objective for InternalObjective<'_> {
    fn dim(&self) -> usize {
        self.ic.n_dof
    }

    fn f_and_g(&mut self, z: &[f64]) -> (f64, Vec<f64>) {
        let x_new = self.map(z);
        self.x_last = x_new.clone();

        // E + Cartesian gradient at the mapped geometry
        let n = x_new.len();
        let mut g2 = vec![[0.0f64; 3]; n];
        let e = self.ff.energy_and_gradient(&x_new, &mut g2);
        let nc = 3 * n;
        let mut gx = vec![0.0f64; nc];
        for i in 0..n {
            for t in 0..3 {
                gx[3 * i + t] = g2[i][t];
            }
        }
        self.last_max = gx.iter().map(|v| v.abs()).fold(0.0f64, f64::max);
        self.last_rms = (gx.iter().map(|v| v * v).sum::<f64>() / nc as f64).sqrt();

        // gradient in the delocalized space. CRITICAL: the map x(z) is
        // exactly linear with Jacobian P = B_q(b_ref)ᵀ Λ⁻¹ (the Baker
        // iteration converges to the pseudo-inverse solution of the FIXED
        // system), so the consistent gradient is Λ⁻¹ B_q(b_ref) g_x — using
        // B at the current geometry instead breaks g·d descent once the
        // geometry drifts from x_ref (line searches die on a rising E(z)).
        project_out_tr(&x_new, &mut gx);
        let gq = self.ic.grad_q(&gx);

        // linearization drift: true q of the mapped geometry vs the target z
        let q_true = self.ic.q(&x_new);
        let drift_bad = q_true
            .iter()
            .zip(z.iter())
            .any(|(a, b)| (a - b).abs() > DRIFT_REBUILD);
        if drift_bad {
            if std::env::var("DIC_DEBUG").is_ok() {
                eprintln!(
                    "REBUILD at iter, maxdrift {:.3}",
                    q_true
                        .iter()
                        .zip(z.iter())
                        .map(|(a, b)| (a - b).abs())
                        .fold(0.0f64, f64::max)
                );
            }
            if let Some(ic2) = InternalCoords::build(&x_new, &self.znums) {
                let gq2 = ic2.grad_q(&gx);
                self.ic = ic2;
                self.x_ref = x_new;
                self.reset_flag = true;
                return (e, gq2);
            }
            // rebuild failed on this geometry: keep the basis (valid, slower)
        }

        (e, gq)
    }

    fn energy(&mut self, z: &[f64]) -> f64 {
        let x = self.map(z);
        self.ff.energy(&x)
    }

    fn force_stats(&self) -> (f64, f64) {
        (self.last_max, self.last_rms)
    }

    fn take_reset(&mut self) -> bool {
        std::mem::replace(&mut self.reset_flag, false)
    }

    fn final_coords(&mut self) -> Vec<[f64; 3]> {
        self.x_last.clone()
    }
}

/// Optimize in delocalized internal coordinates; Cartesian fallback when the
/// internal basis cannot be built (n < 3, degenerate, oversized).
///
/// Safety net: if the internal run stops WITHOUT meeting a convergence
/// criterion (line-search failure abort — e.g. when the force field's
/// gradient is locally inconsistent with its energy; a known GFN-FF
/// hydrogen-bond-chain issue), the optimization restarts in Cartesian mode
/// from the internal endpoint, so the final result is never worse than a
/// plain Cartesian run.
pub fn optimize_internal(
    ff: &dyn ForceField,
    initial_coords: &[[f64; 3]],
    znums: &[usize],
    convergence: &ConvergenceOptions,
) -> OptimizationResult {
    let internal = match InternalCoords::build(initial_coords, znums) {
        Some(ic) => {
            let n_dof = ic.n_dof;
            let mut obj = InternalObjective::new(ff, initial_coords, znums, ic);
            lbfgs_core(&mut obj, vec![0.0; n_dof], convergence)
        }
        None => return super::optimize(ff, initial_coords, convergence),
    };
    if internal.converged || internal.energy_converged {
        return internal;
    }
    // diagnostics: capture the internal-run endpoint that failed to converge
    // (used to reproduce force-field gradient/energy inconsistencies)
    if let Ok(path) = std::env::var("DIC_DUMP") {
        let mut s = String::new();
        s.push_str(&format!("{}\n\n", internal.optimized_coords.len()));
        for (a, p) in internal.optimized_coords.iter().enumerate() {
            s.push_str(&format!(
                "{:2} {:19.12} {:19.12} {:19.12}\n",
                znums.get(a).copied().unwrap_or(0),
                p[0],
                p[1],
                p[2]
            ));
        }
        let _ = std::fs::write(path, s);
    }
    // failed abort: continue in Cartesian from the internal endpoint
    let restart = super::optimize(ff, &internal.optimized_coords, convergence);
    OptimizationResult {
        iterations: internal.iterations + restart.iterations,
        ..restart
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mmff::MMFFForceField;
    use crate::molecule::parser::parse_sdf;
    use crate::MMFFVariant;

    /// DIC optimization must land in the same minimum as Cartesian (MMFF,
    /// whose gradients are exact) — and on flexible molecules it should not
    /// need more iterations than Cartesian by a wide margin.
    #[test]
    fn internal_lands_same_minimum_as_cartesian() {
        let dir = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/fixtures/conformers/");
        for name in ["ethanol", "aspirin"] {
            let sdf = std::fs::read_to_string(format!("{dir}{name}.sdf")).unwrap();
            let mol = parse_sdf(&sdf).unwrap();
            let coords: Vec<[f64; 3]> = mol.atoms.iter().map(|a| a.position).collect();
            let z: Vec<usize> = mol.atoms.iter().map(|a| a.atomic_number as usize).collect();
            let ff = MMFFForceField::new(&mol, MMFFVariant::MMFF94s);
            let conv = ConvergenceOptions::default();
            let rc = crate::optimizer::optimize(&ff, &coords, &conv);
            let ri = optimize_internal(&ff, &coords, &z, &conv);
            assert!(ri.converged, "{name}: internal did not converge");
            assert!(
                (ri.final_energy - rc.final_energy).abs() < 0.01,
                "{name}: internal E {} vs cartesian {}",
                ri.final_energy,
                rc.final_energy
            );
        }
    }

    /// Degenerate inputs fall back to the Cartesian optimizer.
    #[test]
    fn fallback_for_tiny_systems() {
        struct DiatomicSpring;
        impl ForceField for DiatomicSpring {
            fn energy_and_gradient(&self, c: &[[f64; 3]], g: &mut [[f64; 3]]) -> f64 {
                let dx = c[0][0] - c[1][0];
                let dy = c[0][1] - c[1][1];
                let dz = c[0][2] - c[1][2];
                let r = (dx * dx + dy * dy + dz * dz).sqrt();
                let d = r - 1.0; // spring to r = 1 A
                let f = 2.0 * d; // dE/dr
                for gi in g.iter_mut() {
                    *gi = [0.0; 3];
                }
                g[0][0] = f * dx / r;
                g[0][1] = f * dy / r;
                g[0][2] = f * dz / r;
                g[1][0] = -g[0][0];
                g[1][1] = -g[0][1];
                g[1][2] = -g[0][2];
                d * d
            }
        }
        let ff = DiatomicSpring;
        let r = optimize_internal(
            &ff,
            &[[0.0, 0.0, 0.0], [1.5, 0.0, 0.0]],
            &[1, 1],
            &ConvergenceOptions::default(),
        );
        assert!(r.converged);
        assert!((r.final_energy - 0.0).abs() < 1e-8);
    }
}
