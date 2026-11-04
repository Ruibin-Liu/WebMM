//! Climbing-image nudged elastic band (CI-NEB) for saddle refinement.
//!
//! Given two minima of the same molecule (any commensurate atom ordering),
//! NEB finds an approximation of the minimum-energy path between them; the
//! climbing image converges to the highest saddle along it. Used by the
//! iMTD-GC landscape to REFINE the observed transition-barrier upper bounds
//! (window PE peaks) into saddle-energy estimates.
//!
//! Implementation: Kabsch-aligned linear initial path, improved (energy-
//! weighted) tangents, perpendicular true force + tangential spring force,
//! climbing image, FIRE integrator. Fully deterministic — no RNG.

use crate::forces::ForceField;
use crate::metad::imtd::kabsch_transform;

/// NEB configuration.
#[derive(Debug, Clone, Copy)]
pub struct NebConfig {
    /// Interior images between the fixed endpoints. Typical: 12–24.
    pub n_images: usize,
    /// Spring constant (kcal/mol/Å²).
    pub spring: f64,
    /// Iteration cap.
    pub max_iterations: usize,
    /// Convergence: max per-image force NORM over interior images
    /// (kcal/mol/Å). The CI-NEB literature standard is ~0.05 eV/Å ≈ 1.2;
    /// 0.3 is comfortably tighter.
    pub force_tol: f64,
    /// Enable the climbing image (saddle refinement) — run as PHASE 2 after
    /// the standard band relaxed (climbing on an unconverged band wastes the
    /// budget on band relaxation at climbing pace).
    pub climbing: bool,
}

impl Default for NebConfig {
    fn default() -> Self {
        NebConfig {
            n_images: 16,
            // ~0.22 eV/Å² — the NEB literature range is 0.1–0.5 eV/Å²
            // (2.3–11.5 kcal/mol/Å²); weaker springs let images slide off
            // barriers and fake a low saddle
            spring: 5.0,
            max_iterations: 20000,
            force_tol: 0.5,
            // OFF by default: the two-phase machinery is in place, but the
            // QuickMin integrator falls into a limit cycle in the climbing
            // phase (observed maxF oscillating 1.66/2.13 with the saddle
            // frozen below the true crossing) — needs a damped-projected
            // LBFGS-style optimizer, deferred past the feature freeze. The
            // reported saddle is the relaxed band PEAK: an upper bound that
            // decreases with image count.
            climbing: false,
        }
    }
}

/// Result of a NEB run.
#[derive(Debug, Clone)]
pub struct NebResult {
    /// Images along the path, endpoints included (index 0 = A, last = B as
    /// aligned onto A).
    pub path: Vec<Vec<[f64; 3]>>,
    /// Energy of every image (kcal/mol).
    pub energies: Vec<f64>,
    /// Highest interior image energy — the saddle estimate. With climbing
    /// disabled (v1 default) this is the relaxed band PEAK: an upper bound
    /// of the true saddle energy that decreases with image count.
    pub saddle_energy: f64,
    /// saddle_energy minus the LOWER endpoint energy (the forward barrier
    /// from the deeper well).
    pub barrier: f64,
    /// Index of the saddle image in `path`.
    pub saddle_image: usize,
    pub converged: bool,
    pub iterations: usize,
    /// Iteration at which the climbing phase activated (0 = never).
    pub climb_started: usize,
    /// Gradient evaluations skipped by the freeze/lazy policy (diagnostics).
    pub evals_skipped: usize,
    /// Number of images frozen at exit.
    pub n_frozen: usize,
}

/// Per-image QuickMin integrator (the VTST-standard NEB optimizer): each
/// image keeps a velocity that is zeroed when it opposes the force and
/// projected onto the force direction otherwise — unconditionally stable for
/// the non-conservative NEB force, deterministic, with a velocity clamp for
/// the strained linear initial path (whose forces can be huge).
struct QuickMin {
    dt: f64,
    /// Max image velocity norm (Å/step) — bounds the first relaxation steps.
    vmax: f64,
}

impl QuickMin {
    fn step(&self, x: &mut [f64], v: &mut [f64], f: &[f64]) {
        let dot =
            |a: &[f64], b: &[f64]| -> f64 { a.iter().zip(b.iter()).map(|(p, q)| p * q).sum() };
        if dot(v, f) < 0.0 {
            for vi in v.iter_mut() {
                *vi = 0.0;
            }
        }
        for (vi, fi) in v.iter_mut().zip(f.iter()) {
            *vi += self.dt * fi;
        }
        // keep only the force-parallel velocity component
        let f2 = dot(f, f);
        if f2 > 1e-20 {
            let a = dot(v, f) / f2;
            for (vi, fi) in v.iter_mut().zip(f.iter()) {
                *vi = a * fi;
            }
        }
        let vn = dot(v, v).sqrt();
        if vn > self.vmax && vn > 1e-20 {
            let sc = self.vmax / vn;
            for vi in v.iter_mut() {
                *vi *= sc;
            }
        }
        for (xi, vi) in x.iter_mut().zip(v.iter()) {
            *xi += self.dt * vi;
        }
    }
}

/// An image whose force norm stays below `force_tol`·FREEZE_FRAC for this
/// many consecutive iterations is FROZEN: its eval and integration stop.
/// Converged endpoint-side images stop moving long before the saddle region
/// settles — freezing them skips the majority of gradient evaluations.
const FREEZE_PATIENCE: usize = 40;
const FREEZE_FRAC: f64 = 0.5;
/// Collective distance between two images: the norm of the concatenated
/// displacement vector (NOT the sum of per-atom distances — the standard NEB
/// spring uses image norms).
fn dist(a: &[[f64; 3]], b: &[[f64; 3]]) -> f64 {
    a.iter()
        .zip(b.iter())
        .map(|(p, q)| {
            let d = [p[0] - q[0], p[1] - q[1], p[2] - q[2]];
            d[0] * d[0] + d[1] * d[1] + d[2] * d[2]
        })
        .sum::<f64>()
        .sqrt()
}

/// Run CI-NEB between two minima of the same molecule. `b` is aligned onto
/// `a` with a proper rotation (Kabsch) before the path is built.
pub fn neb(ff: &dyn ForceField, a: &[[f64; 3]], b: &[[f64; 3]], config: &NebConfig) -> NebResult {
    let n_atoms = a.len();
    debug_assert_eq!(n_atoms, b.len(), "NEB endpoints must be the same molecule");
    // align b onto a (mobile = b, reference = a)
    let (r, t, _) = kabsch_transform(b, a);
    let b_al: Vec<[f64; 3]> = b
        .iter()
        .map(|p| {
            [
                r[0] * p[0] + r[1] * p[1] + r[2] * p[2] + t[0],
                r[3] * p[0] + r[4] * p[1] + r[5] * p[2] + t[1],
                r[6] * p[0] + r[7] * p[1] + r[8] * p[2] + t[2],
            ]
        })
        .collect();

    // endpoints must be true minima: callers hand over search-ensemble
    // coordinates quenched under an iteration cap that can leave residual
    // forces (≈9 kcal/mol/Å observed) — an under-relaxed endpoint pins its
    // neighbor image and the band can never converge. Tight re-quench.
    let quench = crate::ConvergenceOptions {
        max_iterations: 3000,
        ..Default::default()
    };
    let ra = crate::optimizer::optimize(ff, a, &quench);
    let rb = crate::optimizer::optimize(ff, &b_al, &quench);
    let a: Vec<[f64; 3]> = ra.optimized_coords;
    let b_al: Vec<[f64; 3]> = rb.optimized_coords;

    let n_interior = config.n_images.max(2);
    let total = n_interior + 2;
    // initial path: linear interpolation (flattened per image)
    let mut images: Vec<Vec<[f64; 3]>> = Vec::with_capacity(total);
    for i in 0..total {
        let s = i as f64 / (total - 1) as f64;
        images.push(
            a.iter()
                .zip(b_al.iter())
                .map(|(p, q)| {
                    [
                        p[0] + s * (q[0] - p[0]),
                        p[1] + s * (q[1] - p[1]),
                        p[2] + s * (q[2] - p[2]),
                    ]
                })
                .collect(),
        );
    }

    let mut energies: Vec<f64> = Vec::with_capacity(total);
    let mut grads: Vec<Vec<[f64; 3]>> = Vec::with_capacity(total);
    for img in &images {
        let mut g = vec![[0.0f64; 3]; n_atoms];
        let e = ff.energy_and_gradient(img, &mut g);
        energies.push(e);
        grads.push(g);
    }

    let flatten = |ims: &Vec<Vec<[f64; 3]>>| -> Vec<f64> {
        let mut out = Vec::with_capacity(n_interior * n_atoms * 3);
        for img in ims.iter().take(n_interior + 1).skip(1) {
            for p in img {
                out.extend_from_slice(p);
            }
        }
        out
    };
    let mut x = flatten(&images);
    let mut v = vec![0.0f64; x.len()];
    let qm = QuickMin {
        dt: 0.015,
        vmax: 0.06,
    };

    // two-phase: standard band first (climbing off), climbing image after
    // the band relaxes (force tol) or after 60% of the budget
    let mut climb_started = 0usize;
    let mut converged = false;
    let mut iterations = 0usize;
    let mut saddle_image = 0usize;
    let mut prev_max_force = f64::INFINITY;

    // hoisted iteration buffers (the force loop ran 4+ Vec allocs per image
    // per iteration — ~125k allocations on a 2k-iteration band)
    let mut cur: Vec<Vec<[f64; 3]>> = images.clone();
    let mut f: Vec<f64> = vec![0.0; n_interior * n_atoms * 3];
    let mut image_norm2: Vec<f64> = vec![0.0; n_interior + 1];
    let mut tp: Vec<[f64; 3]> = vec![[0.0; 3]; n_atoms];
    let mut tm: Vec<[f64; 3]> = vec![[0.0; 3]; n_atoms];
    let mut tangent: Vec<[f64; 3]> = vec![[0.0; 3]; n_atoms];
    // freeze state per interior image (index k-1)
    let mut frozen: Vec<bool> = vec![false; n_interior];
    let mut calm: Vec<usize> = vec![0; n_interior];
    let mut evals_skipped: usize = 0;

    for it in 0..config.max_iterations {
        iterations = it + 1;
        // unpack interior images into the geometry buffer
        for (k, chunk) in x.chunks(n_atoms * 3).enumerate() {
            for (ai, c) in chunk.iter().enumerate() {
                cur[k + 1][ai / 3][ai % 3] = *c;
            }
        }
        // refresh energies/gradients: frozen images skip entirely; active
        // images always re-evaluate (a stale-gradient coasting experiment
        // made the tail dynamics 6× WORSE — 1937 → 12293 iterations)
        for k in 1..=n_interior {
            if frozen[k - 1] {
                evals_skipped += 1;
                continue;
            }
            energies[k] = ff.energy_and_gradient(&cur[k], &mut grads[k]);
        }

        // climbing image = highest interior energy, activated in phase 2
        let band_relaxed = prev_max_force < config.force_tol * 3.0;
        if config.climbing
            && climb_started == 0
            && (band_relaxed || it >= config.max_iterations * 3 / 5)
        {
            climb_started = it + 1;
        }
        let climb_active = config.climbing && climb_started > 0;
        let climb_idx = if climb_active {
            let ci = (1..=n_interior)
                .max_by(|&i, &j| energies[i].total_cmp(&energies[j]))
                .unwrap_or(1);
            saddle_image = ci;
            Some(ci)
        } else {
            None
        };

        // forces on interior images (flattened, hoisted buffers); convergence
        // tracks the per-image force NORM. Frozen images: spring/tangent
        // terms only need POSITIONS (unchanged), the gradient is cached, so
        // their force is recomputed cheaply from cache — but they are
        // excluded from the convergence max and from integration below.
        let mut max_force = 0.0f64;
        for v in image_norm2.iter_mut() {
            *v = 0.0;
        }
        for k in 1..=n_interior {
            // improved tangent (Henkelman et al. 2000) — written into the
            // hoisted tp/tm/tangent buffers
            for i in 0..n_atoms {
                tp[i] = [
                    cur[k + 1][i][0] - cur[k][i][0],
                    cur[k + 1][i][1] - cur[k][i][1],
                    cur[k + 1][i][2] - cur[k][i][2],
                ];
                tm[i] = [
                    cur[k][i][0] - cur[k - 1][i][0],
                    cur[k][i][1] - cur[k - 1][i][1],
                    cur[k][i][2] - cur[k - 1][i][2],
                ];
            }
            let (e_prev, e_next) = (energies[k - 1], energies[k + 1]);
            if e_next > energies[k] && energies[k] > e_prev {
                tangent.copy_from_slice(&tp);
            } else if e_next < energies[k] && energies[k] < e_prev {
                tangent.copy_from_slice(&tm);
            } else {
                let dmax = (e_next - energies[k])
                    .abs()
                    .max((energies[k] - e_prev).abs());
                let dmin = (e_next - energies[k])
                    .abs()
                    .min((energies[k] - e_prev).abs());
                let (w_plus, w_minus) = if e_next > e_prev {
                    (dmax, dmin)
                } else {
                    (dmin, dmax)
                };
                for i in 0..n_atoms {
                    for c in 0..3 {
                        tangent[i][c] = w_plus * tp[i][c] + w_minus * tm[i][c];
                    }
                }
            }
            // normalize the tangent per-image (collective norm)
            let t_norm: f64 = tangent
                .iter()
                .map(|p| p[0] * p[0] + p[1] * p[1] + p[2] * p[2])
                .sum::<f64>()
                .sqrt();
            if t_norm < 1e-12 {
                continue; // degenerate band
            }
            let inv_t = 1.0 / t_norm;

            let is_climb = climb_idx == Some(k);
            // spring force along the tangent
            let spring_f = if is_climb {
                0.0
            } else {
                config.spring * (dist(&cur[k + 1], &cur[k]) - dist(&cur[k], &cur[k - 1]))
            };
            let coef = if is_climb {
                // climbing image: full gradient with the tangential
                // component INVERTED (drives up to the saddle)
                2.0
            } else {
                1.0
            };
            for i in 0..n_atoms {
                let (tx, ty, tz) = (
                    tangent[i][0] * inv_t,
                    tangent[i][1] * inv_t,
                    tangent[i][2] * inv_t,
                );
                let gdot = grads[k][i][0] * tx + grads[k][i][1] * ty + grads[k][i][2] * tz;
                for c in 0..3 {
                    let t_c = [tx, ty, tz][c];
                    // F = -∇E + coef·(∇E·t̂)t̂ + k·Δℓ·t̂
                    let mut fc = -grads[k][i][c] + coef * gdot * t_c;
                    if !is_climb {
                        fc += spring_f * t_c;
                    }
                    let dst = (k - 1) * n_atoms * 3 + i * 3 + c;
                    f[dst] = fc;
                    image_norm2[k] += fc * fc;
                }
            }
        }
        // freeze bookkeeping + convergence over ACTIVE images only. Iron
        // rule: the highest-energy image (the saddle carrier) NEVER freezes
        // — a prematurely frozen peak locks a distorted band (measured
        // barrier 6.03 vs the physical ~3.3).
        let peak = (1..=n_interior)
            .max_by(|&a, &b| energies[a].total_cmp(&energies[b]))
            .unwrap_or(1);
        let mut all_frozen = true;
        for k in 1..=n_interior {
            let norm = image_norm2[k].sqrt();
            if frozen[k - 1] {
                continue;
            }
            all_frozen = false;
            if norm < config.force_tol * FREEZE_FRAC && k != peak {
                calm[k - 1] += 1;
                if calm[k - 1] >= FREEZE_PATIENCE {
                    frozen[k - 1] = true;
                }
            } else {
                calm[k - 1] = 0;
            }
            max_force = max_force.max(norm);
        }
        if all_frozen {
            converged = true;
            break;
        }
        // phase 2 restarts the relaxation: the climbing image injects force
        // when activated, so only check convergence in a stable phase
        if max_force < config.force_tol && (!climb_active || it > climb_started + 20) {
            converged = true;
            break;
        }
        prev_max_force = max_force;
        // per-image QuickMin steps (the NEB force couples images only
        // through neighbors — independent image steps are the VTST norm)
        let stride = n_atoms * 3;
        for k in 0..n_interior {
            if frozen[k] {
                for vi in &mut v[k * stride..(k + 1) * stride] {
                    *vi = 0.0;
                }
                continue;
            }
            let lo = k * stride;
            let hi = lo + stride;
            qm.step(&mut x[lo..hi], &mut v[lo..hi], &f[lo..hi]);
        }
    }

    // final geometry/energy refresh (x may have moved after the last eval)
    let mut final_path: Vec<Vec<[f64; 3]>> = images.clone();
    for (k, chunk) in x.chunks(n_atoms * 3).enumerate() {
        for (ai, c) in chunk.iter().enumerate() {
            final_path[k + 1][ai / 3][ai % 3] = *c;
        }
    }
    let mut final_energies = energies.clone();
    for k in 1..=n_interior {
        let mut g = vec![[0.0f64; 3]; n_atoms];
        final_energies[k] = ff.energy_and_gradient(&final_path[k], &mut g);
    }
    let (saddle_energy, saddle_idx) = (1..=n_interior)
        .map(|k| (final_energies[k], k))
        .max_by(|a, b| a.0.total_cmp(&b.0))
        .unwrap_or((final_energies[0], 0));
    let e_low = final_energies[0].min(final_energies[total - 1]);
    NebResult {
        path: final_path,
        energies: final_energies,
        saddle_energy,
        barrier: saddle_energy - e_low,
        saddle_image: if saddle_image > 0 {
            saddle_image
        } else {
            saddle_idx
        },
        converged,
        iterations,
        climb_started,
        evals_skipped,
        n_frozen: frozen.iter().filter(|&&b| b).count(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mmff::MMFFForceField;
    use crate::molecule::parser::parse_sdf;

    const BUTANE_HEAVY_SDF: &str = "Butane\n     RDKit          3D\n\n  4  3  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.5500    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    3.0500    0.1000    0.9000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    4.5500    0.1000    0.9000 C   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0\n  2  3  1  0\n  3  4  1  0\nM  END\n";

    fn butane_rotamers() -> (Vec<[f64; 3]>, Vec<[f64; 3]>, crate::molecule::Molecule) {
        // heavy zigzag → add H; gauche built by rotating the last C + its H
        // around the central bond (the same construction the imtd crossing
        // test uses) then both endpoints fully quenched
        let heavy = parse_sdf(BUTANE_HEAVY_SDF).expect("parse");
        let mol = crate::molecule::hydrogens::add_hydrogens(&heavy);
        let anti: Vec<[f64; 3]> = mol.atoms.iter().map(|a| a.position).collect();
        // rotate everything bonded beyond the central C1-C2 bond by 120°
        let (j, k) = (1usize, 2usize);
        let axis = [
            anti[k][0] - anti[j][0],
            anti[k][1] - anti[j][1],
            anti[k][2] - anti[j][2],
        ];
        let norm = (axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2]).sqrt();
        let u = [axis[0] / norm, axis[1] / norm, axis[2] / norm];
        let th = std::f64::consts::PI * 2.0 / 3.0;
        let (s, c) = (th.sin(), th.cos());
        // BFS from k avoiding j
        let mut seen = vec![false; mol.atoms.len()];
        seen[j] = true;
        seen[k] = true;
        let mut f2 = vec![k];
        let mut queue = std::collections::VecDeque::new();
        queue.push_back(k);
        while let Some(z) = queue.pop_front() {
            for &w in &mol.adjacency[z] {
                if !seen[w] {
                    seen[w] = true;
                    f2.push(w);
                    queue.push_back(w);
                }
            }
        }
        let mut gauche = anti.clone();
        for &i in &f2 {
            let d = [
                anti[i][0] - anti[k][0],
                anti[i][1] - anti[k][1],
                anti[i][2] - anti[k][2],
            ];
            let cr = [
                u[1] * d[2] - u[2] * d[1],
                u[2] * d[0] - u[0] * d[2],
                u[0] * d[1] - u[1] * d[0],
            ];
            let dot = u[0] * d[0] + u[1] * d[1] + u[2] * d[2];
            for x in 0..3 {
                gauche[i][x] = anti[k][x] + c * d[x] + s * cr[x] + u[x] * dot * (1.0 - c);
            }
        }
        // quench both endpoints to real minima (H placement is heuristic)
        let conv = crate::ConvergenceOptions {
            max_iterations: 500,
            ..Default::default()
        };
        let ff = MMFFForceField::new(&mol, crate::MMFFVariant::MMFF94s);
        let ra = crate::optimizer::optimize(&ff, &anti, &conv);
        let rb = crate::optimizer::optimize(&ff, &gauche, &conv);
        (ra.optimized_coords, rb.optimized_coords, mol)
    }

    #[test]
    fn neb_butane_anti_gauche_barrier_is_physical() {
        let (a, b, mol) = butane_rotamers();
        let ff = MMFFForceField::new(&mol, crate::MMFFVariant::MMFF94s);
        let cfg = NebConfig::default();
        let res = neb(&ff, &a, &b, &cfg);
        let ea = res.energies[0];
        let eb = *res.energies.last().expect("endpoint");
        assert!(
            res.saddle_energy > ea.max(eb),
            "saddle must bound both wells: saddle {} vs {ea}/{eb}",
            res.saddle_energy
        );
        assert!(
            (2.0..6.0).contains(&res.barrier),
            "butane central-torsion MMFF barrier ~3-4 kcal/mol, got {}",
            res.barrier
        );
        // v1 semantics: converged standard band; the peak is an upper bound
        assert!(
            res.converged,
            "standard band must converge, used {}",
            res.iterations
        );
        assert_eq!(res.climb_started, 0, "climbing off by default");
        // endpoints stay fixed
        for (p, q) in res.path[0].iter().zip(a.iter()) {
            assert_eq!(p, q, "endpoint A fixed");
        }
        // determinism
        let res2 = neb(&ff, &a, &b, &cfg);
        assert_eq!(res.saddle_energy, res2.saddle_energy);
        assert_eq!(res.energies, res2.energies);
    }

    #[test]
    fn neb_relaxes_below_the_linear_interpolation_peak() {
        let (a, b, mol) = butane_rotamers();
        let ff = MMFFForceField::new(&mol, crate::MMFFVariant::MMFF94s);
        // linear-interpolation peak energy
        let n = 17usize;
        let mut lin_peak = f64::NEG_INFINITY;
        for i in 0..n {
            let s = i as f64 / (n - 1) as f64;
            let img: Vec<[f64; 3]> = a
                .iter()
                .zip(b.iter())
                .map(|(p, q)| {
                    [
                        p[0] + s * (q[0] - p[0]),
                        p[1] + s * (q[1] - p[1]),
                        p[2] + s * (q[2] - p[2]),
                    ]
                })
                .collect();
            let mut g = vec![[0.0f64; 3]; mol.atoms.len()];
            lin_peak = lin_peak.max(ff.energy_and_gradient(&img, &mut g));
        }
        let res = neb(&ff, &a, &b, &NebConfig::default());
        assert!(
            res.saddle_energy <= lin_peak + 1e-6,
            "NEB relaxes the band: {} vs linear peak {lin_peak}",
            res.saddle_energy
        );
    }

    #[test]
    fn neb_trivial_path_when_endpoints_coincide() {
        let (a, _b, mol) = butane_rotamers();
        let ff = MMFFForceField::new(&mol, crate::MMFFVariant::MMFF94s);
        let res = neb(&ff, &a, &a, &NebConfig::default());
        let ea = res.energies[0];
        assert!(
            (res.saddle_energy - ea).abs() < 1e-6,
            "identical endpoints → flat band, got {} vs {ea}",
            res.saddle_energy
        );
    }
}
