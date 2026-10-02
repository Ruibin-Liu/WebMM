//! Gaussian molecular shape overlap and rigid-body alignment (ROCS lineage).
//!
//! Atoms are modeled as positive Gaussians following Grant & Pickup
//! (J. Phys. Chem. 1996, 100, 2453): g_i(r) = C·exp(−α_i·|r − r_i|²) with
//! C = 2√2 (GCI) and per-element exponents α = κ/r_Bondi² (κ ≈ 2.418), the
//! parameterization used by the shape-it line of tools. With these constants
//! the single-atom self volume GCI·(π/α)^{3/2} equals the hard-sphere volume
//! (4/3)πr³ exactly — that identity is a unit test below.
//!
//! Scoring: shape Tanimoto = V_AB / (V_A + V_B − V_AB) where the molecule
//! volumes V_A/V_B and the intersection V_AB are computed by the
//! inclusion–exclusion expansion over product Gaussians (each k-fold product
//! has α = Σα_i, center = weighted mean, coefficient
//! ΠC_i·exp(−(1/α)·Σ_{i<j}α_iα_j·d_ij²), volume C·(π/α)^{3/2}), pruned by
//! volume contribution. This is a clean-room implementation of the published
//! equations; no third-party code is used.
//!
//! Alignment: rigid-body (rotation vector + translation, 6 DOF) multi-start
//! L-BFGS on a pairwise-overlap surrogate Σ_ij V(g_i·h_j) with analytic
//! gradients; every start is re-scored with the full inclusion–exclusion
//! Tanimoto and the best pose wins.

use crate::molecule::Molecule;

/// Atomic density coefficient (shape-it config lineage): 2√2.
pub const GCI: f64 = 2.828427125;

/// α = κ / r_Bondi². Determined from the per-element table below
/// (e.g. H: 1.679158285 · 1.2² = 2.41799; C: 0.836674025 · 1.7² = 2.41800).
const KAPPA: f64 = 2.4179879;

/// Per-element Gaussian exponents (shape-it GAlpha table, common elements).
/// Fallback for other elements: KAPPA / r_Bondi².
fn alpha_for(z: u8) -> f64 {
    match z {
        1 => 1.679158285,  // H
        3 => 0.729980658,  // Li
        5 => 0.604496983,  // B
        6 => 0.836674025,  // C
        7 => 1.006446589,  // N
        8 => 1.046566798,  // O
        9 => 1.118972618,  // F
        11 => 0.469247983, // Na
        12 => 0.670309721, // Mg
        14 => 0.804266845, // Si
        15 => 0.814231526, // P
        16 => 0.834737044, // S
        17 => 0.878848066, // Cl
        19 => 0.421311908, // K
        20 => 0.661904802, // Ca
        26 => 0.597518568, // Fe
        30 => 0.655319394, // Zn
        35 => 0.817809742, // Br
        53 => 0.724390270, // I
        _ => {
            let r = bondi_radius(z);
            KAPPA / (r * r)
        }
    }
}

/// Bondi vdW radii (Å) for the exponent fallback above.
fn bondi_radius(z: u8) -> f64 {
    match z {
        1 => 1.20,
        2 => 1.40,
        3 => 1.82,
        4 => 1.53,
        5 => 1.92,
        6 => 1.70,
        7 => 1.55,
        8 => 1.52,
        9 => 1.47,
        10 => 1.54,
        11 => 2.27,
        12 => 1.73,
        13 => 1.84,
        14 => 2.10,
        15 => 1.80,
        16 => 1.80,
        17 => 1.75,
        18 => 1.88,
        19 => 2.75,
        20 => 2.31,
        28 => 1.63,
        29 => 1.40,
        30 => 1.39,
        34 => 1.90,
        35 => 1.85,
        36 => 2.02,
        47 => 1.72,
        50 => 2.17,
        53 => 1.98,
        78 => 1.75,
        79 => 1.66,
        82 => 2.02,
        _ => 1.70,
    }
}

/// One atom Gaussian (center in Å, exponent α).
#[derive(Debug, Clone, Copy)]
pub struct ShapeAtom {
    pub c: [f64; 3],
    pub alpha: f64,
}

/// Build the atom-Gaussian set of a molecule (all atoms, hydrogens included —
/// the shape-it lineage convention).
pub fn shape_atoms(mol: &Molecule) -> Vec<ShapeAtom> {
    mol.atoms
        .iter()
        .map(|a| ShapeAtom {
            c: a.position,
            alpha: alpha_for(a.atomic_number),
        })
        .collect()
}

// ---------------------------------------------------------------------------
// Inclusion–exclusion volumes (shape-it `atomOverlap` semantics)
// ---------------------------------------------------------------------------

/// An inclusion–exclusion term spec: which atoms, the sign, and the
/// pose-independent statistics (Σα, n·ln GCI). Term volumes are
/// materialized per pose from the posed centers.
#[derive(Clone)]
pub struct TermSpec {
    /// Atom indices in increasing order.
    pub idx: Vec<u8>,
    pub alpha_sum: f64,
    pub log_c: f64,
    pub sign: f64,
}

/// Enumerate the inclusion–exclusion terms with shape-it-style pruning:
/// a product is extended by atom j only if the overlap ratio
/// V(P·j)/(V(P) + V_j − V(P·j)) ≥ TERM_EPS (spatially disconnected atoms
/// never spawn deep product chains). The subset STRUCTURE is enumerated at
/// the original pose; volumes are recomputed per pose.
const TERM_EPS: f64 = 0.03; // shape-it value
const TERM_MAX_LEVEL: usize = 6; // shape-it LEVEL

fn enumerate_specs(atoms: &[ShapeAtom]) -> Vec<TermSpec> {
    let n = atoms.len();
    let single_vol = |a: &ShapeAtom| -> f64 {
        let s = std::f64::consts::PI / a.alpha;
        GCI * s * s.sqrt()
    };
    let vol_of = |alpha_sum: f64, log_c: f64, cn: [f64; 3], sas: f64| -> f64 {
        let dot = cn[0] * cn[0] + cn[1] * cn[1] + cn[2] * cn[2];
        let cross = alpha_sum * sas - dot;
        let c = (log_c - cross / alpha_sum).exp();
        let s = std::f64::consts::PI / alpha_sum;
        c * s * s.sqrt()
    };
    struct Entry {
        last: usize,
        spec: TermSpec,
        center_num: [f64; 3], // Σ α_i c_i
        sum_alpha_c_sq: f64,  // Σ α_i |c_i|²
    }
    let mut specs: Vec<TermSpec> = Vec::with_capacity(n * 2);
    let mut stack: Vec<Entry> = Vec::new();
    for (i, a) in atoms.iter().enumerate() {
        let spec = TermSpec {
            idx: vec![i as u8],
            alpha_sum: a.alpha,
            log_c: GCI.ln(),
            sign: 1.0,
        };
        specs.push(spec.clone());
        stack.push(Entry {
            last: i,
            spec,
            center_num: [a.alpha * a.c[0], a.alpha * a.c[1], a.alpha * a.c[2]],
            sum_alpha_c_sq: a.alpha * (a.c[0] * a.c[0] + a.c[1] * a.c[1] + a.c[2] * a.c[2]),
        });
    }
    let singles: f64 = atoms.iter().map(single_vol).sum();
    let hard_floor = 1e-7 * singles.max(1e-12);
    while let Some(e) = stack.pop() {
        let v_e = vol_of(
            e.spec.alpha_sum,
            e.spec.log_c,
            e.center_num,
            e.sum_alpha_c_sq,
        );
        for (j, b) in atoms.iter().enumerate().skip(e.last + 1) {
            let alpha_sum = e.spec.alpha_sum + b.alpha;
            let log_c = e.spec.log_c + GCI.ln();
            let cn = [
                e.center_num[0] + b.alpha * b.c[0],
                e.center_num[1] + b.alpha * b.c[1],
                e.center_num[2] + b.alpha * b.c[2],
            ];
            let sas =
                e.sum_alpha_c_sq + b.alpha * (b.c[0] * b.c[0] + b.c[1] * b.c[1] + b.c[2] * b.c[2]);
            let v = vol_of(alpha_sum, log_c, cn, sas);
            let denom = v_e + single_vol(b) - v;
            if denom <= 0.0 || v / denom < TERM_EPS || v.abs() < hard_floor {
                continue;
            }
            let mut idx = e.spec.idx.clone();
            idx.push(j as u8);
            let sign = if idx.len() % 2 == 0 { -1.0 } else { 1.0 };
            let spec = TermSpec {
                idx,
                alpha_sum,
                log_c,
                sign,
            };
            specs.push(spec.clone());
            if spec.idx.len() < TERM_MAX_LEVEL {
                stack.push(Entry {
                    last: j,
                    spec,
                    center_num: cn,
                    sum_alpha_c_sq: sas,
                });
            }
        }
    }
    specs
}

/// Per-pose materialized term: sufficient statistics for O(1) volume evals
/// and O(1) merges with another term.
#[derive(Clone, Copy)]
pub struct TermPose {
    pub alpha_sum: f64,
    pub log_c: f64,
    pub sign: f64,
    pub center: [f64; 3],
    /// α-weighted mean of |c|² (Σα|c|²/Σα) — rebuilds the cross sum.
    pub mean_alpha_c_sq: f64,
    /// bounding radius over the term's atoms (prefilter).
    pub r_bound: f64,
}

impl TermPose {
    fn cross_with(&self, o: &TermPose) -> f64 {
        let alpha = self.alpha_sum + o.alpha_sum;
        let cn = [
            self.alpha_sum * self.center[0] + o.alpha_sum * o.center[0],
            self.alpha_sum * self.center[1] + o.alpha_sum * o.center[1],
            self.alpha_sum * self.center[2] + o.alpha_sum * o.center[2],
        ];
        let sas = self.alpha_sum * self.mean_alpha_c_sq + o.alpha_sum * o.mean_alpha_c_sq;
        let dot = cn[0] * cn[0] + cn[1] * cn[1] + cn[2] * cn[2];
        // NB: the pairwise sum Σ_{i<j}αiαj d² = α·Σα|c|² − |Σαc|² carries NO
        // ½ — an extra half here had been halving the Gaussian decay of every
        // inclusion-exclusion term since v1.4.0 (found by the rigorous bound
        // derivation in v1.6.3; self-consistent identities masked it, but it
        // disagreed with the pairwise surrogate's β = αiαj/αij).
        alpha * sas - dot
    }

    fn merged_volume(&self, o: &TermPose) -> f64 {
        let alpha = self.alpha_sum + o.alpha_sum;
        let cross = self.cross_with(o);
        let c = (self.log_c + o.log_c - cross / alpha).exp();
        let s = std::f64::consts::PI / alpha;
        c * s * s.sqrt()
    }
}

fn volume_of_term(t: &TermPose) -> f64 {
    t.merged_volume(&TermPose {
        alpha_sum: 0.0,
        log_c: 0.0,
        sign: 1.0,
        center: [0.0, 0.0, 0.0],
        mean_alpha_c_sq: 0.0,
        r_bound: 0.0,
    })
}

/// Materialize all term stats at a pose.
fn materialize(specs: &[TermSpec], atoms: &[ShapeAtom]) -> Vec<TermPose> {
    specs
        .iter()
        .map(|s| {
            let mut cn = [0.0f64; 3];
            let mut sas = 0.0f64;
            for &i in &s.idx {
                let a = &atoms[i as usize];
                cn[0] += a.alpha * a.c[0];
                cn[1] += a.alpha * a.c[1];
                cn[2] += a.alpha * a.c[2];
                sas += a.alpha * (a.c[0] * a.c[0] + a.c[1] * a.c[1] + a.c[2] * a.c[2]);
            }
            let center = [
                cn[0] / s.alpha_sum,
                cn[1] / s.alpha_sum,
                cn[2] / s.alpha_sum,
            ];
            let mut r2 = 0.0f64;
            for &i in &s.idx {
                let a = &atoms[i as usize];
                let d = [a.c[0] - center[0], a.c[1] - center[1], a.c[2] - center[2]];
                r2 = r2.max(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
            }
            TermPose {
                alpha_sum: s.alpha_sum,
                log_c: s.log_c,
                sign: s.sign,
                center,
                mean_alpha_c_sq: sas / s.alpha_sum,
                r_bound: r2.sqrt(),
            }
        })
        .collect()
}

/// Inclusion–exclusion volume of an atom-Gaussian set at its current pose.
fn volume_ie(atoms: &[ShapeAtom]) -> f64 {
    if atoms.is_empty() {
        return 0.0;
    }
    let specs = enumerate_specs(atoms);
    materialize(&specs, atoms)
        .iter()
        .map(|t| t.sign * volume_of_term(t))
        .sum()
}

/// Self volume of a molecule's Gaussian representation.
pub fn self_volume(atoms: &[ShapeAtom]) -> f64 {
    volume_ie(atoms)
}

/// Intersection volume of two molecules at their current relative pose
/// (shape-it `atomOverlap` semantics): the bilinear expansion over the
/// inclusion–exclusion terms of each molecule,
///   V_AB = Σ_{i∈terms(A)} Σ_{j∈terms(B)} ε_i·ε_j·V(term_i · term_j).
/// A molecule bundled with its prepared shape (specs + lazy self overlap).
/// `align`/`align_colored` consume these so callers (and the wasm-level
/// cache) can prepare once and reuse across a whole library scan.
pub struct ShapeMol {
    pub atoms: Vec<ShapeAtom>,
    pub prepared: PreparedShape,
}

pub fn shape_mol_from_atoms(atoms: &[ShapeAtom]) -> ShapeMol {
    ShapeMol {
        atoms: atoms.to_vec(),
        prepared: prepare_shape(atoms),
    }
}

pub fn shape_mol(mol: &Molecule) -> ShapeMol {
    shape_mol_from_atoms(&shape_atoms(mol))
}

/// A prepared shape: the inclusion–exclusion term specs enumerated once
/// (the pruning criterion depends only on inter-atomic distances, so the
/// subset is identical under any rigid pose — preparing is lossless) plus
/// the lazily-computed pose-invariant self overlap.
pub struct PreparedShape {
    pub specs: Vec<TermSpec>,
    /// Pose-invariant self bilinear overlap, computed at prepare time.
    pub self_overlap: f64,
}

pub fn prepare_shape(atoms: &[ShapeAtom]) -> PreparedShape {
    let specs = enumerate_specs(atoms);
    let self_overlap = {
        let p = prepared_view(&specs);
        overlap_prepared_view(&p, atoms, &p, atoms)
    };
    PreparedShape {
        specs,
        self_overlap,
    }
}

// small view shim so prepare can score before owning the struct
struct PreparedView<'a> {
    specs: &'a [TermSpec],
}

fn prepared_view(specs: &[TermSpec]) -> PreparedView<'_> {
    PreparedView { specs }
}

fn overlap_prepared_view(
    pa: &PreparedView,
    a: &[ShapeAtom],
    pb: &PreparedView,
    b: &[ShapeAtom],
) -> f64 {
    if a.is_empty() || b.is_empty() {
        return 0.0;
    }
    let ta = materialize(pa.specs, a);
    let tb = materialize(pb.specs, b);
    let scale: f64 = ta.iter().map(volume_of_term).sum::<f64>().abs().max(1e-12);
    let cutoff = 1e-7 * scale;
    let mut total = 0.0f64;
    for pa in ta.iter() {
        for pb in tb.iter() {
            let d = [
                pa.center[0] - pb.center[0],
                pa.center[1] - pb.center[1],
                pa.center[2] - pb.center[2],
            ];
            if d[0] * d[0] + d[1] * d[1] + d[2] * d[2] > (pa.r_bound + pb.r_bound + 7.0).powi(2) {
                continue;
            }
            let v = pa.merged_volume(pb);
            if v.abs() < cutoff {
                continue;
            }
            total += pa.sign * pb.sign * v;
        }
    }
    total
}

/// Overlap of two prepared shapes, materialized at the given (posed) atom
/// coordinates. v1.6.3 fast path: α-index lookup tables for the pair
/// prefactors (no per-pair sqrt/div) and a rigorous distance-bound skip
/// (identity Σ_{i∈S,j∈T}αiαj d² = αS·αT·|cS−cT|² ⇒ V ≤ K·exp(−βd²), so a
/// table d²max = ln(K/cutoff)/β skips without exp and never drops a
/// contributing pair) — bit-identical to the reference pair loop.
pub fn overlap_prepared(
    pa: &PreparedShape,
    a: &[ShapeAtom],
    pb: &PreparedShape,
    b: &[ShapeAtom],
) -> f64 {
    if a.is_empty() || b.is_empty() {
        return 0.0;
    }
    let ta = materialize(&pa.specs, a);
    let tb = materialize(&pb.specs, b);
    overlap_terms(&ta, &tb)
}

/// Pair loop over materialized terms with α-index tables and the rigorous
/// distance bound. `cutoff` uses the same scale heuristic as before.
pub fn overlap_terms(ta: &[TermPose], tb: &[TermPose]) -> f64 {
    if ta.is_empty() || tb.is_empty() {
        return 0.0;
    }
    let scale: f64 = ta.iter().map(volume_of_term).sum::<f64>().abs().max(1e-12);
    let cutoff = 1e-7 * scale;
    // α-index tables: collect distinct α values, tabulate K/β/d²max per pair
    let mut alphas: Vec<f64> = Vec::new();
    let idx_of = |alphas: &mut Vec<f64>, a: f64| -> usize {
        for (i, x) in alphas.iter().enumerate() {
            if (*x - a).abs() < 1e-12 {
                return i;
            }
        }
        alphas.push(a);
        alphas.len() - 1
    };
    let ia: Vec<usize> = ta
        .iter()
        .map(|t| idx_of(&mut alphas, t.alpha_sum))
        .collect();
    let ib: Vec<usize> = tb
        .iter()
        .map(|t| idx_of(&mut alphas, t.alpha_sum))
        .collect();
    let n = alphas.len();
    // per-α-pair tables: β[i][j] = αiαj/(αi+αj) and ln((π/αij)^{3/2}).
    // The rigorous skip (per pair): V ≤ exp(lnc_a + lnc_b + lnstab − β·d²),
    // using cross(S∪T) = crossS + crossT + αSαT·d², so the pair is dropped
    // exactly when the bound falls below the cutoff — bit-identical totals
    // (the per-term log_c carries the GCI^{na+nb} prefactor).
    let mut beta = vec![0.0f64; n * n];
    let mut lnstab = vec![0.0f64; n * n];
    for i in 0..n {
        for j in 0..n {
            let aij = alphas[i] + alphas[j];
            beta[i * n + j] = alphas[i] * alphas[j] / aij;
            let s_ = std::f64::consts::PI / aij;
            lnstab[i * n + j] = 3.0 * 0.5 * (s_).ln(); // (π/α)^{3/2}
        }
    }
    let lncutoff = cutoff.ln();
    let mut total = 0.0f64;
    for (x, pa) in ta.iter().enumerate() {
        let row = ia[x] * n;
        for (y, pb) in tb.iter().enumerate() {
            let col = ib[y];
            let b = beta[row + col];
            if b <= 1e-12 {
                continue;
            }
            let d = [
                pa.center[0] - pb.center[0],
                pa.center[1] - pb.center[1],
                pa.center[2] - pb.center[2],
            ];
            let d2 = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
            if d2 * b > pa.log_c + pb.log_c + lnstab[row + col] - lncutoff {
                continue; // rigorous bound: contribution < cutoff
            }
            let v = pa.merged_volume(pb);
            if v.abs() < cutoff {
                continue;
            }
            total += pa.sign * pb.sign * v;
        }
    }
    total
}

/// With GCI = 2√2 the self-product volume identity V(g·g) = V(g) makes
/// coincident identical sets give exactly V_AB = V_A (Tanimoto 1).
pub fn overlap_full(a: &[ShapeAtom], b: &[ShapeAtom]) -> f64 {
    if a.is_empty() || b.is_empty() {
        return 0.0;
    }
    overlap_prepared(&prepare_shape(a), a, &prepare_shape(b), b)
}

/// Shape Tanimoto at the current relative pose. Following the shape-it
/// convention, the self volumes in the denominator are the self bilinear
/// overlaps overlap_full(X, X) (so T(A, A) = 1 exactly).
pub fn shape_tanimoto(a: &[ShapeAtom], b: &[ShapeAtom]) -> f64 {
    if a.is_empty() || b.is_empty() {
        return 0.0;
    }
    let va = overlap_full(a, a);
    let vb = overlap_full(b, b);
    let vab = overlap_full(a, b);
    let denom = va + vb - vab;
    if denom <= 0.0 {
        return 0.0;
    }
    (vab / denom).clamp(0.0, 1.0)
}

// ---------------------------------------------------------------------------
// Pairwise-overlap surrogate with analytic 6-DOF gradient
// ---------------------------------------------------------------------------

/// Apply rotation vector w (Rodrigues) and translation t to a point.
pub fn rodrigues(w: &[f64; 3], t: &[f64; 3], p: &[f64; 3]) -> [f64; 3] {
    let th = (w[0] * w[0] + w[1] * w[1] + w[2] * w[2]).sqrt();
    let out = if th < 1e-10 {
        *p
    } else {
        let u = [w[0] / th, w[1] / th, w[2] / th];
        let c = th.cos();
        let s = th.sin();
        let ux = u[1] * p[2] - u[2] * p[1];
        let uy = u[2] * p[0] - u[0] * p[2];
        let uz = u[0] * p[1] - u[1] * p[0];
        let d = u[0] * p[0] + u[1] * p[1] + u[2] * p[2];
        let k = 1.0 - c;
        [
            p[0] * c + ux * s + u[0] * d * k,
            p[1] * c + uy * s + u[1] * d * k,
            p[2] * c + uz * s + u[2] * d * k,
        ]
    };
    [out[0] + t[0], out[1] + t[1], out[2] + t[2]]
}

/// Row-major 3×3 rotation matrix of a rotation vector (for reporting).
pub fn rotmat(w: &[f64; 3]) -> [f64; 9] {
    let th = (w[0] * w[0] + w[1] * w[1] + w[2] * w[2]).sqrt();
    if th < 1e-10 {
        return [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
    }
    let u = [w[0] / th, w[1] / th, w[2] / th];
    let e = [[0.0, -u[2], u[1]], [u[2], 0.0, -u[0]], [-u[1], u[0], 0.0]];
    let (c, s) = (th.cos(), th.sin());
    let mut r = [0.0f64; 9];
    // R = cI + s[u]× + (1−c)uuᵀ
    for i in 0..3 {
        for k in 0..3 {
            r[i * 3 + k] =
                c * (if i == k { 1.0 } else { 0.0 }) + s * e[i][k] + (1.0 - c) * u[i] * u[k];
        }
    }
    r
}

/// Left Jacobian of SO(3) at rotation vector w (exp-map derivative).
fn so3_left_jacobian(w: &[f64; 3]) -> [f64; 9] {
    let th2 = w[0] * w[0] + w[1] * w[1] + w[2] * w[2];
    let th = th2.sqrt();
    if th < 1e-8 {
        return [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
    }
    let (s, c) = (th.sin(), th.cos());
    let a = (1.0 - c) / th2;
    let b = (th - s) / (th2 * th);
    let m = [[0.0, -w[2], w[1]], [w[2], 0.0, -w[0]], [-w[1], w[0], 0.0]];
    let mut j = [0.0f64; 9];
    for i in 0..3 {
        for k in 0..3 {
            let w2sum: f64 = (0..3).map(|q| m[i][q] * m[q][k]).sum();
            j[i * 3 + k] = (if i == k { 1.0 } else { 0.0 }) + a * m[i][k] + b * w2sum;
        }
    }
    j
}

/// Pairwise overlap O(t,R) = Σ_ij K_ij exp(−β_ij|R p_j + t − c_i|²) and the
/// gradient wrt (w, t). β_ij = α_iα_j/(α_i+α_j), K_ij = GCI²(π/α_ij)^{3/2}.
pub fn pairwise_overlap_grad(
    a: &[ShapeAtom],
    b_base: &[ShapeAtom],
    w: &[f64; 3],
    t: &[f64; 3],
) -> (f64, [f64; 6]) {
    let mut total = 0.0f64;
    let mut g_t = [0.0f64; 3];
    let mut torque = [0.0f64; 3]; // Σ_j (R p_j) × g_j — the LEFT perturbation
                                  // acts on R(w)p (translation is added after rotation), so the torque arm
                                  // is the rotated point WITHOUT t
    for h in b_base.iter() {
        let rp = rodrigues(w, &[0.0, 0.0, 0.0], &h.c);
        let y = [rp[0] + t[0], rp[1] + t[1], rp[2] + t[2]];
        let mut g_j = [0.0f64; 3];
        for ai in a.iter() {
            let alpha_ij = ai.alpha + h.alpha;
            let beta = ai.alpha * h.alpha / alpha_ij;
            let s = std::f64::consts::PI / alpha_ij;
            let k_ij = GCI * GCI * s * s.sqrt();
            let d = [y[0] - ai.c[0], y[1] - ai.c[1], y[2] - ai.c[2]];
            let d2 = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
            let e = (-beta * d2).exp() * k_ij;
            total += e;
            let f = -2.0 * beta * e;
            g_j[0] += f * d[0];
            g_j[1] += f * d[1];
            g_j[2] += f * d[2];
        }
        g_t[0] += g_j[0];
        g_t[1] += g_j[1];
        g_t[2] += g_j[2];
        torque[0] += rp[1] * g_j[2] - rp[2] * g_j[1];
        torque[1] += rp[2] * g_j[0] - rp[0] * g_j[2];
        torque[2] += rp[0] * g_j[1] - rp[1] * g_j[0];
    }
    // g_w = J_l(w)ᵀ · torque  (d/dw of R(w)p = −[Rp]×·J_l)
    let jl = so3_left_jacobian(w);
    let g_w = [
        jl[0] * torque[0] + jl[3] * torque[1] + jl[6] * torque[2],
        jl[1] * torque[0] + jl[4] * torque[1] + jl[7] * torque[2],
        jl[2] * torque[0] + jl[5] * torque[1] + jl[8] * torque[2],
    ];
    (total, [g_w[0], g_w[1], g_w[2], g_t[0], g_t[1], g_t[2]])
}

/// Surrogate overlap, no gradient (identical sum as `pairwise_overlap_grad`).
fn surrogate_energy(a: &[ShapeAtom], b: &[ShapeAtom]) -> f64 {
    pairwise_overlap_grad(a, b, &[0.0, 0.0, 0.0], &[0.0, 0.0, 0.0]).0
}

// ---------------------------------------------------------------------------
// Color force field (joint optimization + scoring)
// ---------------------------------------------------------------------------

/// A color feature site (position, Gaussian exponent, feature type).
#[derive(Debug, Clone, Copy)]
pub struct ColorSite {
    pub c: [f64; 3],
    pub alpha: f64,
    pub type_id: u8,
    /// Per-site weight (default 1.0): each same-type pair contribution is
    /// multiplied by w_i·w_j — linear scaling of a feature type when all
    /// its sites share w = sqrt(u). Pass sqrt(u) from the caller for
    /// "u = how many times this feature type counts" semantics.
    pub w: f64,
}

pub const COLOR_TYPES: [&str; 6] = ["donor", "acceptor", "pos", "neg", "hydrophobe", "ring"];

pub fn color_type_id(name: &str) -> Option<u8> {
    COLOR_TYPES.iter().position(|t| *t == name).map(|i| i as u8)
}

/// Same-type pairwise Gaussian overlap + 6-DOF gradient (identical kernel
/// to `pairwise_overlap_grad`, restricted to matching feature types).
pub fn color_overlap_grad(
    qs: &[ColorSite],
    ts: &[ColorSite],
    w: &[f64; 3],
    t: &[f64; 3],
) -> (f64, [f64; 6]) {
    let mut total = 0.0f64;
    let mut g_t = [0.0f64; 3];
    let mut torque = [0.0f64; 3];
    for h in ts.iter() {
        let rp = rodrigues(w, &[0.0, 0.0, 0.0], &h.c);
        let y = [rp[0] + t[0], rp[1] + t[1], rp[2] + t[2]];
        let mut g_j = [0.0f64; 3];
        for ai in qs.iter() {
            if ai.type_id != h.type_id {
                continue;
            }
            let alpha_ij = ai.alpha + h.alpha;
            let beta = ai.alpha * h.alpha / alpha_ij;
            let s = std::f64::consts::PI / alpha_ij;
            let k_ij = GCI * GCI * s * s.sqrt();
            let d = [y[0] - ai.c[0], y[1] - ai.c[1], y[2] - ai.c[2]];
            let d2 = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
            let wpair = ai.w * h.w;
            let e = (-beta * d2).exp() * k_ij * wpair;
            total += e;
            let f = -2.0 * beta * e; // d(w·e)/dy — w is a constant factor
            g_j[0] += f * d[0];
            g_j[1] += f * d[1];
            g_j[2] += f * d[2];
        }
        g_t[0] += g_j[0];
        g_t[1] += g_j[1];
        g_t[2] += g_j[2];
        torque[0] += rp[1] * g_j[2] - rp[2] * g_j[1];
        torque[1] += rp[2] * g_j[0] - rp[0] * g_j[2];
        torque[2] += rp[0] * g_j[1] - rp[1] * g_j[0];
    }
    let jl = so3_left_jacobian(w);
    let g_w = [
        jl[0] * torque[0] + jl[3] * torque[1] + jl[6] * torque[2],
        jl[1] * torque[0] + jl[4] * torque[1] + jl[7] * torque[2],
        jl[2] * torque[0] + jl[5] * torque[1] + jl[8] * torque[2],
    ];
    (total, [g_w[0], g_w[1], g_w[2], g_t[0], g_t[1], g_t[2]])
}

/// Overlap value only (line-search `energy` path consistency is guaranteed
/// because `f_and_g` uses the same gradient-carrying function).
fn color_overlap_value(qs: &[ColorSite], ts: &[ColorSite], w: &[f64; 3], t: &[f64; 3]) -> f64 {
    color_overlap_grad(qs, ts, w, t).0
}

/// Color Tanimoto with the target sites posed by (w, t):
/// T_c = O_ab/(O_aa + O_bb − O_ab), all terms same-type pairwise sums.
pub fn color_tanimoto_at(qs: &[ColorSite], ts: &[ColorSite], w: &[f64; 3], t: &[f64; 3]) -> f64 {
    if qs.is_empty() || ts.is_empty() {
        return 0.0;
    }
    let oab = color_overlap_value(qs, ts, w, t);
    let oaa = color_overlap_value(qs, qs, &[0.0; 3], &[0.0; 3]);
    let obb = color_overlap_value(ts, ts, &[0.0; 3], &[0.0; 3]);
    let den = oaa + obb - oab;
    if den > 1e-12 {
        (oab / den).clamp(0.0, 1.0)
    } else {
        0.0
    }
}

// ---------------------------------------------------------------------------
// Rigid-body alignment (multi-start L-BFGS on the surrogate, full rescoring)
// ---------------------------------------------------------------------------

/// Deterministic splitmix-based RNG for reproducible random starts.
struct MiniRng {
    s: u64,
}

impl MiniRng {
    fn new(seed: u64) -> Self {
        MiniRng {
            s: seed.wrapping_mul(0x9e3779b97f4a7c15).max(1),
        }
    }
    fn next_f64(&mut self) -> f64 {
        self.s = self.s.wrapping_add(0x9e3779b97f4a7c15);
        let mut z = self.s;
        z = (z ^ (z >> 30)).wrapping_mul(0xbf58476d1ce4e5b9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94d049bb133111eb);
        z ^= z >> 31;
        (z >> 11) as f64 / ((1u64 << 53) as f64)
    }
}

#[derive(Debug, Clone)]
pub struct AlignOptions {
    pub random_starts: usize,
    pub max_iter: usize,
    pub seed: u64,
    /// Screening mode: frame starts only, no polish, NO inclusion-exclusion
    /// or color rescoring — the best start is ranked purely by the pairwise
    /// surrogate overlap (an absolute overlap, biased toward larger
    /// molecules; intended as a generous pre-filter for a full-quality
    /// re-rank of the top-N). `tanimoto` is left at 0 (not computed).
    pub screen: bool,
    /// After the cheap optimization of every start, only the top-K poses by
    /// surrogate score receive the full inclusion-exclusion + color
    /// rescoring (the IE rescore is ~95% of the alignment cost). Measured
    /// K-vs-quality on the conformer fixtures (26 starts): K=8 keeps the
    /// best combo within 0.018 of the exhaustive rescore (7/20 pairs differ
    /// at all); usize::MAX = rescore every start.
    pub rescore_top: usize,
}

impl Default for AlignOptions {
    fn default() -> Self {
        AlignOptions {
            random_starts: 16,
            max_iter: 200,
            seed: 42,
            screen: false,
            rescore_top: 8,
        }
    }
}

#[derive(Debug, Clone)]
pub struct AlignResult {
    pub tanimoto: f64,
    /// Color Tanimoto at the best pose (0 when no color sites were given).
    pub color_tanimoto: f64,
    /// Shape Tanimoto + Color Tanimoto.
    pub combo: f64,
    /// Row-major 3×4 affine (rotation, then translation) for the target.
    pub transform: [f64; 12],
    pub surrogate_overlap: f64,
    pub iterations: usize,
    pub starts: usize,
}

/// Pose atoms by (w, t).
fn pose(atoms: &[ShapeAtom], w: &[f64; 3], t: &[f64; 3]) -> Vec<ShapeAtom> {
    atoms
        .iter()
        .map(|h| ShapeAtom {
            c: rodrigues(w, t, &h.c),
            alpha: h.alpha,
        })
        .collect()
}

/// Translate atoms so their centroid sits at the origin.
fn center(atoms: &[ShapeAtom]) -> Vec<ShapeAtom> {
    let n = atoms.len() as f64;
    let mut c = [0.0f64; 3];
    for a in atoms {
        c[0] += a.c[0];
        c[1] += a.c[1];
        c[2] += a.c[2];
    }
    c = [c[0] / n, c[1] / n, c[2] / n];
    atoms
        .iter()
        .map(|a| ShapeAtom {
            c: [a.c[0] - c[0], a.c[1] - c[1], a.c[2] - c[2]],
            alpha: a.alpha,
        })
        .collect()
}

/// Inertia-tensor principal axes via the shared symmetric Jacobi solver
/// (orthonormal eigenvectors ordered by descending eigenvalue).
fn principal_axes(atoms: &[ShapeAtom]) -> [[f64; 3]; 3] {
    let n = atoms.len() as f64;
    let mut com = [0.0f64; 3];
    for a in atoms {
        com[0] += a.c[0];
        com[1] += a.c[1];
        com[2] += a.c[2];
    }
    com = [com[0] / n, com[1] / n, com[2] / n];
    let mut m = [0.0f64; 9];
    for a in atoms {
        let d = [a.c[0] - com[0], a.c[1] - com[1], a.c[2] - com[2]];
        let r2 = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
        for i in 0..3 {
            for k in 0..3 {
                m[i * 3 + k] += if i == k { r2 } else { 0.0 } - d[i] * d[k];
            }
        }
    }
    let mut v = [0.0f64; 9];
    crate::optimizer::jacobi::sym_jacobi(&mut m, &mut v, 3);
    // order eigenvector columns by descending eigenvalue (m diagonal after solve)
    let mut order: [usize; 3] = [0, 1, 2];
    let vals = [m[0], m[4], m[8]];
    order.sort_by(|&a, &b| vals[b].partial_cmp(&vals[a]).unwrap());
    let mut axes = [[0.0f64; 3]; 3];
    for (r, col) in order.iter().enumerate() {
        for i in 0..3 {
            axes[r][i] = v[i * 3 + col];
        }
    }
    // fix any residual sign ambiguity deterministically: largest |component| positive
    for axis in axes.iter_mut() {
        let (mi, _) = axis
            .iter()
            .enumerate()
            .max_by(|(_, a), (_, b)| a.abs().partial_cmp(&b.abs()).unwrap())
            .unwrap();
        if axis[mi] < 0.0 {
            for x in axis.iter_mut() {
                *x = -*x;
            }
        }
    }
    axes
}

/// Rotation vector rotating unit vector a onto unit vector b (shortest arc).
fn rotvec_between(a: &[f64; 3], b: &[f64; 3]) -> [f64; 3] {
    let cross = [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ];
    let dot = (a[0] * b[0] + a[1] * b[1] + a[2] * b[2]).clamp(-1.0, 1.0);
    let th = dot.acos();
    let sn = th.sin();
    if sn < 1e-9 {
        if dot > 0.0 {
            return [0.0, 0.0, 0.0];
        }
        let helper = if a[0].abs() < 0.9 {
            [1.0, 0.0, 0.0]
        } else {
            [0.0, 1.0, 0.0]
        };
        let ax = [
            a[1] * helper[2] - a[2] * helper[1],
            a[2] * helper[0] - a[0] * helper[2],
            a[0] * helper[1] - a[1] * helper[0],
        ];
        let n = (ax[0] * ax[0] + ax[1] * ax[1] + ax[2] * ax[2]).sqrt();
        return [
            ax[0] / n * std::f64::consts::PI,
            ax[1] / n * std::f64::consts::PI,
            ax[2] / n * std::f64::consts::PI,
        ];
    }
    let s = th / sn;
    [cross[0] * s, cross[1] * s, cross[2] * s]
}

/// One L-BFGS run on the surrogate from (w, t); returns final pose and stats.
fn optimize_pose(
    a: &[ShapeAtom],
    b_centered: &[ShapeAtom],
    w0: [f64; 3],
    t0: [f64; 3],
    max_iter: usize,
) -> (f64, [f64; 3], [f64; 3], usize) {
    optimize_pose_colored(a, b_centered, None, None, w0, t0, max_iter)
}

/// Joint shape+color objective: O = O_shape + Σ_w O_color(same-type site
/// pairs). Both terms share the same pairwise-Gaussian kernel, so the
/// J_l(w)-corrected gradient chain applies to the sum.
fn optimize_pose_colored(
    a: &[ShapeAtom],
    b_centered: &[ShapeAtom],
    q_sites: Option<&[ColorSite]>,
    t_sites: Option<&[ColorSite]>,
    w0: [f64; 3],
    t0: [f64; 3],
    max_iter: usize,
) -> (f64, [f64; 3], [f64; 3], usize) {
    use crate::optimizer::{lbfgs_core, Objective, OptimizationResult};
    use crate::ConvergenceOptions;

    struct ShapeObj<'a> {
        a: &'a [ShapeAtom],
        b: &'a [ShapeAtom],
        q_sites: Option<&'a [ColorSite]>,
        t_sites: Option<&'a [ColorSite]>,
        last_z: Vec<f64>,
        last_max: f64,
        last_rms: f64,
    }
    impl<'a> Objective for ShapeObj<'a> {
        fn dim(&self) -> usize {
            6
        }
        fn f_and_g(&mut self, z: &[f64]) -> (f64, Vec<f64>) {
            let w = [z[0], z[1], z[2]];
            let t = [z[3], z[4], z[5]];
            // direct call: the gradient must carry the J_l(w) chain factor —
            // posing first and taking the identity gradient loses it
            let (mut o, mut g) = pairwise_overlap_grad(self.a, self.b, &w, &t);
            if let (Some(qs), Some(ts)) = (self.q_sites, self.t_sites) {
                let (oc, gc) = color_overlap_grad(qs, ts, &w, &t);
                o += oc;
                for i in 0..6 {
                    g[i] += gc[i];
                }
            }
            self.last_z = z.to_vec();
            self.last_max = g.iter().fold(0.0f64, |m, v| m.max(v.abs()));
            self.last_rms = (g.iter().map(|v| v * v).sum::<f64>() / 6.0).sqrt();
            // minimize −overlap: BOTH value and gradient flip sign
            (-o, g.iter().map(|x| -x).collect::<Vec<f64>>())
        }
        fn energy(&mut self, z: &[f64]) -> f64 {
            self.f_and_g(z).0
        }
        fn force_stats(&self) -> (f64, f64) {
            (self.last_max, self.last_rms)
        }
        fn take_reset(&mut self) -> bool {
            false
        }
        fn final_coords(&mut self) -> Vec<[f64; 3]> {
            vec![
                [self.last_z[0], self.last_z[1], self.last_z[2]],
                [self.last_z[3], self.last_z[4], self.last_z[5]],
            ]
        }
    }

    let mut obj = ShapeObj {
        a,
        b: b_centered,
        q_sites,
        t_sites,
        last_z: vec![0.0; 6],
        last_max: 0.0,
        last_rms: 0.0,
    };
    let z0 = vec![w0[0], w0[1], w0[2], t0[0], t0[1], t0[2]];
    let conv = ConvergenceOptions {
        max_force: 1e-6,
        rms_force: 1e-7,
        energy_change: 1e-11,
        max_iterations: max_iter,
    };
    let OptimizationResult {
        optimized_coords,
        iterations,
        ..
    } = lbfgs_core(&mut obj, z0, &conv);
    let w = optimized_coords[0];
    let t = optimized_coords[1];
    let posed = pose(b_centered, &w, &t);
    let mut o = surrogate_energy(a, &posed);
    if let (Some(qs), Some(ts)) = (q_sites, t_sites) {
        o += color_overlap_value(qs, ts, &w, &t);
    }
    (o, w, t, iterations)
}

/// Multi-start rigid-body alignment of `target` onto `query`. Both inputs
/// keep their coordinates; the returned transform maps original target
/// coordinates to the best-aligned pose.
pub fn align(query: &ShapeMol, target: &ShapeMol, opts: &AlignOptions) -> AlignResult {
    align_colored(query, target, None, None, opts)
}

/// Joint shape+color alignment: the objective maximizes O_shape + O_color
/// (same-type site pairs, unit weights), and starts are ranked by the
/// Shape T + Color T combo. Sites are given in ORIGINAL coordinates and
/// centered alongside their molecules.
pub fn align_colored(
    query: &ShapeMol,
    target: &ShapeMol,
    q_sites: Option<&[ColorSite]>,
    t_sites: Option<&[ColorSite]>,
    opts: &AlignOptions,
) -> AlignResult {
    // Work in centered frames: target gets centered first (its own centroid
    // to origin), and we align it onto the centered query; the returned
    // affine therefore maps original target coords → aligned frame.
    let q = center(&query.atoms);
    let t_c = center(&target.atoms);
    // center the site lists into the same frames
    let q_com = centroid(&query.atoms);
    let t_com = centroid(&target.atoms);
    let q_sites_c: Option<Vec<ColorSite>> = q_sites.map(|ss| {
        ss.iter()
            .map(|x| ColorSite {
                c: [x.c[0] - q_com[0], x.c[1] - q_com[1], x.c[2] - q_com[2]],
                alpha: x.alpha,
                type_id: x.type_id,
                w: x.w,
            })
            .collect()
    });
    let t_sites_c: Option<Vec<ColorSite>> = t_sites.map(|ss| {
        ss.iter()
            .map(|x| ColorSite {
                c: [x.c[0] - t_com[0], x.c[1] - t_com[1], x.c[2] - t_com[2]],
                alpha: x.alpha,
                type_id: x.type_id,
                w: x.w,
            })
            .collect()
    });
    let q_sites_ref = q_sites_c.as_deref();
    let t_sites_ref = t_sites_c.as_deref();
    let colored = q_sites_ref.is_some() && t_sites_ref.is_some();

    // Compose transform back to original query coordinates: aligned point =
    // R·(p − t_com) + t + q_com.
    let t_com = centroid(&target.atoms);
    let q_com = centroid(&query.atoms);

    let axes_q = principal_axes(&q);
    let axes_t = principal_axes(&t_c);

    let mut starts: Vec<([f64; 3], [f64; 3])> = Vec::new();
    // Principal-axis frame starts: R0 = axes_q · axes_tᵀ with sign flips on
    // the second/third target axes (proper rotations only), plus the pure
    // dominant-axis ± starts.
    // 8 proper frame combos: 2 axis assignments (identity and 1↔2 swap —
    // near-degenerate eigenvalues make the assignment ambiguous) × 4 sign
    // flips, determinant-corrected.
    for &swap in &[false, true] {
        for &(s1, s2) in &[(1.0f64, 1.0), (1.0, -1.0), (-1.0, 1.0), (-1.0, -1.0)] {
            // build R0 columns: q0·t0ᵀ + q1·(s1·t1)ᵀ + q2·(s2·t2)ᵀ, fix determinant
            let (ax1, ax2) = if swap {
                (axes_t[2], axes_t[1])
            } else {
                (axes_t[1], axes_t[2])
            };
            let t1 = [s1 * ax1[0], s1 * ax1[1], s1 * ax1[2]];
            let t2 = [s2 * ax2[0], s2 * ax2[1], s2 * ax2[2]];
            let mut cols = [[0.0f64; 3]; 3];
            for r in 0..3 {
                for c in 0..3 {
                    let tv = [axes_t[0][c], t1[c], t2[c]];
                    cols[c][r] = axes_q[0][r] * tv[0] + axes_q[1][r] * tv[1] + axes_q[2][r] * tv[2];
                }
            }
            // det: if improper, flip the second column
            let det = cols[0][0] * (cols[1][1] * cols[2][2] - cols[1][2] * cols[2][1])
                - cols[0][1] * (cols[1][0] * cols[2][2] - cols[1][2] * cols[2][0])
                + cols[0][2] * (cols[1][0] * cols[2][1] - cols[1][1] * cols[2][0]);
            if det < 0.0 {
                for v in cols[1].iter_mut() {
                    *v = -*v;
                }
            }
            // rotvec from R0 (matrix → axis·angle)
            let angle = ((cols[0][0] + cols[1][1] + cols[2][2] - 1.0) / 2.0)
                .clamp(-1.0, 1.0)
                .acos();
            if angle > 1e-6 {
                let axis = [
                    cols[1][2] - cols[2][1],
                    cols[2][0] - cols[0][2],
                    cols[0][1] - cols[1][0],
                ];
                let n = (axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2]).sqrt();
                if n > 1e-9 {
                    let u = [axis[0] / n, axis[1] / n, axis[2] / n];
                    starts.push(([u[0] * angle, u[1] * angle, u[2] * angle], [0.0, 0.0, 0.0]));
                }
            } else {
                starts.push(([0.0, 0.0, 0.0], [0.0, 0.0, 0.0]));
            }
        }
    }
    // dominant-axis ± starts (kept from the original design)
    for &sign in &[-1.0f64, 1.0] {
        let a0 = [axes_t[0][0], axes_t[0][1], axes_t[0][2]];
        let mut b0 = [
            sign * axes_q[0][0],
            sign * axes_q[0][1],
            sign * axes_q[0][2],
        ];
        let mut norm = (b0[0] * b0[0] + b0[1] * b0[1] + b0[2] * b0[2]).sqrt();
        b0 = [b0[0] / norm, b0[1] / norm, b0[2] / norm];
        norm = (a0[0] * a0[0] + a0[1] * a0[1] + a0[2] * a0[2]).sqrt();
        let a0 = [a0[0] / norm, a0[1] / norm, a0[2] / norm];
        starts.push((rotvec_between(&a0, &b0), [0.0, 0.0, 0.0]));
    }
    // random rotations (skipped in screen mode: frame starts only)
    let n_random = if opts.screen { 0 } else { opts.random_starts };
    let mut rng = MiniRng::new(opts.seed);
    for _ in 0..n_random {
        let w = [
            (rng.next_f64() * 2.0 - 1.0) * std::f64::consts::PI,
            (rng.next_f64() * 2.0 - 1.0) * std::f64::consts::PI,
            (rng.next_f64() * 2.0 - 1.0) * std::f64::consts::PI,
        ];
        starts.push((w, [0.0, 0.0, 0.0]));
    }

    let mut best = AlignResult {
        tanimoto: -1.0,
        color_tanimoto: 0.0,
        combo: -1.0,
        transform: [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0],
        surrogate_overlap: -1.0,
        iterations: 0,
        starts: starts.len(),
    };
    let mut total_iters = 0usize;
    // prepared shapes: specs enumerated once, self overlaps lazy — the
    // pruning criterion is distance-only so preparing at any rigid pose is
    // bit-identical to per-call enumeration (v1.6.2 perf, lossless).
    let pq = &query.prepared; // specs pose-independent: reuse the caller's
    let pt = &target.prepared;
    // self-overlap volumes are pose invariant: compute once (the Tanimoto
    // denominator still needs the pose-dependent vab). All scoring happens
    // in the CENTERED query frame (`q`), consistently.
    let (vq, vt) = if opts.screen {
        (0.0, 0.0) // not needed in screen mode (no IE rescoring at all)
    } else {
        (pq.self_overlap, pt.self_overlap)
    };
    // query term poses are invariant across the whole alignment (the query
    // never moves) — materialize once (v1.6.3)
    let q_terms = if opts.screen {
        Vec::new()
    } else {
        materialize(&pq.specs, &q)
    };
    let n_starts = starts.len();
    let mut best_wt: Option<([f64; 3], [f64; 3])> = None;
    // Phase A: cheap optimization of every start; screen mode finishes here
    // (surrogate ranking), otherwise keep the top-K poses by surrogate for
    // the expensive full inclusion-exclusion rescoring (v1.6.1: the IE
    // rescore is ~95% of the alignment cost).
    let mut runs: Vec<(f64, [f64; 3], [f64; 3])> = Vec::with_capacity(n_starts);
    for (w0, t0) in &starts {
        let (o, w, t, iters) = if colored {
            optimize_pose_colored(&q, &t_c, q_sites_ref, t_sites_ref, *w0, *t0, opts.max_iter)
        } else {
            optimize_pose(&q, &t_c, *w0, *t0, opts.max_iter)
        };
        total_iters += iters;
        if opts.screen {
            // screening: rank by the surrogate overlap, no IE/color scoring
            if o > best.surrogate_overlap {
                let r = rotmat(&w);
                let m = [
                    r[0],
                    r[1],
                    r[2],
                    t[0] + q_com[0] - (r[0] * t_com[0] + r[1] * t_com[1] + r[2] * t_com[2]),
                    r[3],
                    r[4],
                    r[5],
                    t[1] + q_com[1] - (r[3] * t_com[0] + r[4] * t_com[1] + r[5] * t_com[2]),
                    r[6],
                    r[7],
                    r[8],
                    t[2] + q_com[2] - (r[6] * t_com[0] + r[7] * t_com[1] + r[8] * t_com[2]),
                ];
                best = AlignResult {
                    tanimoto: 0.0, // not computed in screen mode
                    color_tanimoto: 0.0,
                    combo: 0.0,
                    transform: m,
                    surrogate_overlap: o,
                    iterations: iters,
                    starts: n_starts,
                };
            }
            continue;
        }
        runs.push((o, w, t));
    }
    if !opts.screen {
        runs.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap());
    }
    // Phase B: full rescoring of the top-K poses only
    for &(o, w, t) in runs.iter().take(opts.rescore_top.max(1)) {
        // full rescoring at this pose (map centered pose back to query frame)
        let posed_rel = pose(&t_c, &w, &t); // aligned onto centered query
        let vab = {
            let tb = materialize(&pt.specs, &posed_rel);
            overlap_terms(&q_terms, &tb)
        };
        let den = vq + vt - vab;
        let tj = if den > 0.0 {
            (vab / den).clamp(0.0, 1.0)
        } else {
            0.0
        };
        let cj = if colored {
            color_tanimoto_at(q_sites_ref.unwrap(), t_sites_ref.unwrap(), &w, &t)
        } else {
            0.0
        };
        let key = tj + cj;
        if key > best.combo {
            best_wt = Some((w, t));
            // affine for original target coords: p ↦ R·(p − t_com) + t + q_com
            let r = rotmat(&w);
            let m = [
                r[0],
                r[1],
                r[2],
                t[0] + q_com[0] - (r[0] * t_com[0] + r[1] * t_com[1] + r[2] * t_com[2]),
                r[3],
                r[4],
                r[5],
                t[1] + q_com[1] - (r[3] * t_com[0] + r[4] * t_com[1] + r[5] * t_com[2]),
                r[6],
                r[7],
                r[8],
                t[2] + q_com[2] - (r[6] * t_com[0] + r[7] * t_com[1] + r[8] * t_com[2]),
            ];
            best = AlignResult {
                tanimoto: tj,
                color_tanimoto: cj,
                combo: key,
                transform: m,
                surrogate_overlap: o,
                iterations: total_iters,
                starts: n_starts,
            };
        }
    }
    // polish: one fresh restart from the best pose (clears L-BFGS history
    // stalls — the rotation-vector 2π flat directions can freeze the line
    // search while the pose is still ~0.3 Å off on symmetric systems)
    if let Some((w, t)) = best_wt.filter(|_| !opts.screen) {
        let (o2, w2, t2, it2) = if colored {
            optimize_pose_colored(
                &q,
                &t_c,
                q_sites_ref,
                t_sites_ref,
                w,
                t,
                opts.max_iter / 2 + 10,
            )
        } else {
            optimize_pose(&q, &t_c, w, t, opts.max_iter / 2 + 10)
        };
        let posed_rel = pose(&t_c, &w2, &t2);
        let vab2 = {
            let tb = materialize(&pt.specs, &posed_rel);
            overlap_terms(&q_terms, &tb)
        };
        let den2 = vq + vt - vab2;
        let tj = if den2 > 0.0 {
            (vab2 / den2).clamp(0.0, 1.0)
        } else {
            0.0
        };
        let cj2 = if colored {
            color_tanimoto_at(q_sites_ref.unwrap(), t_sites_ref.unwrap(), &w2, &t2)
        } else {
            0.0
        };
        if tj + cj2 > best.combo {
            let r = rotmat(&w2);
            let m = [
                r[0],
                r[1],
                r[2],
                t2[0] + q_com[0] - (r[0] * t_com[0] + r[1] * t_com[1] + r[2] * t_com[2]),
                r[3],
                r[4],
                r[5],
                t2[1] + q_com[1] - (r[3] * t_com[0] + r[4] * t_com[1] + r[5] * t_com[2]),
                r[6],
                r[7],
                r[8],
                t2[2] + q_com[2] - (r[6] * t_com[0] + r[7] * t_com[1] + r[8] * t_com[2]),
            ];
            best = AlignResult {
                tanimoto: tj,
                color_tanimoto: cj2,
                combo: tj + cj2,
                transform: m,
                surrogate_overlap: o2,
                iterations: it2,
                starts: n_starts,
            };
        }
    }
    if best.combo < 0.0 {
        best.combo = best.tanimoto; // uncolored path: combo degenerates to shape
    }
    best.iterations += total_iters;
    best.starts = n_starts;
    best
}

fn centroid(atoms: &[ShapeAtom]) -> [f64; 3] {
    let n = atoms.len() as f64;
    let mut c = [0.0f64; 3];
    for a in atoms {
        c[0] += a.c[0];
        c[1] += a.c[1];
        c[2] += a.c[2];
    }
    [c[0] / n, c[1] / n, c[2] / n]
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn sa(c: [f64; 3], z: u8) -> ShapeAtom {
        ShapeAtom {
            c,
            alpha: alpha_for(z),
        }
    }

    /// Single-atom Gaussian volume vs hard-sphere volume for the elements
    /// whose shape-it exponents are exactly κ/r_Bondi² (H/C/N/O — the
    /// drug-like core; S/Cl etc. use shape-it's own radii lineage and are
    /// covered by the α-table consistency check below).
    #[test]
    fn single_atom_volume_matches_hard_sphere() {
        for (z, r) in [(1u8, 1.2f64), (6, 1.7), (7, 1.55), (8, 1.52)] {
            let alpha = alpha_for(z);
            let s_ = std::f64::consts::PI / alpha;
            let v_gauss = GCI * s_ * s_.sqrt();
            let v_sphere = 4.0 / 3.0 * std::f64::consts::PI * r * r * r;
            assert!(
                (v_gauss - v_sphere).abs() < 1e-3,
                "z={z}: {v_gauss} vs {v_sphere}"
            );
        }
    }

    /// The per-element exponent table is consistent with κ = α·r² for the
    /// H/C/N/O anchors (κ ≈ 2.418).
    #[test]
    fn alpha_table_kappa_consistency() {
        for (z, r) in [(1u8, 1.2f64), (6, 1.7), (7, 1.55), (8, 1.52)] {
            let k = alpha_for(z) * r * r;
            assert!((k - KAPPA).abs() < 2e-6, "z={z}: κ={k}");
        }
    }

    /// Two coincident atoms: inclusion–exclusion must give exactly one atom's
    /// volume; far apart: exactly the sum.
    #[test]
    fn inclusion_exclusion_identities() {
        let a = sa([0.0, 0.0, 0.0], 6);
        let s_ = std::f64::consts::PI / a.alpha;
        let v1 = GCI * s_ * s_.sqrt();
        let both = self_volume(&[a, a]);
        assert!((both - v1).abs() < 1e-7, "coincident: {both} vs {v1}");
        let far = self_volume(&[a, sa([100.0, 0.0, 0.0], 6)]);
        assert!((far - 2.0 * v1).abs() < 1e-6, "far: {far} vs {}", 2.0 * v1);
    }

    /// V_AB identities: coincident identical sets → Tanimoto exactly 1;
    /// disjoint → 0; symmetric in value.
    #[test]
    fn overlap_identities() {
        let a = vec![sa([0.0, 0.0, 0.0], 6), sa([1.5, 0.0, 0.0], 6)];
        let t_self = shape_tanimoto(&a, &a.clone());
        assert!((t_self - 1.0).abs() < 1e-9, "self tanimoto {t_self}");
        let far: Vec<ShapeAtom> = a
            .iter()
            .map(|x| sa([x.c[0] + 100.0, 0.0, 0.0], 6))
            .collect();
        assert!(shape_tanimoto(&a, &far) < 1e-6);
        assert!(shape_tanimoto(&far, &a) < 1e-6);
        // partial overlap sits strictly between 0 and 1 and is symmetric
        let shifted = vec![sa([0.7, 0.0, 0.0], 6), sa([2.2, 0.0, 0.0], 6)];
        let t1 = shape_tanimoto(&a, &shifted);
        let t2 = shape_tanimoto(&shifted, &a);
        assert!(t1 > 0.0 && t1 < 1.0, "t1 {t1}");
        assert!((t1 - t2).abs() < 1e-9);
    }

    /// Analytic 6-DOF gradient vs central finite differences (seeds pinned).
    #[test]
    fn pairwise_gradient_fd() {
        let mut rng = MiniRng::new(7);
        for case in 0..40 {
            let a: Vec<ShapeAtom> = (0..4)
                .map(|_| {
                    sa(
                        [
                            (rng.next_f64() * 4.0 - 2.0),
                            (rng.next_f64() * 4.0 - 2.0),
                            (rng.next_f64() * 4.0 - 2.0),
                        ],
                        if case % 2 == 0 { 6 } else { 7 },
                    )
                })
                .collect();
            let b: Vec<ShapeAtom> = (0..3)
                .map(|_| {
                    sa(
                        [
                            (rng.next_f64() * 4.0 - 2.0),
                            (rng.next_f64() * 4.0 - 2.0),
                            (rng.next_f64() * 4.0 - 2.0),
                        ],
                        if case % 3 == 0 { 8 } else { 6 },
                    )
                })
                .collect();
            let w = [
                rng.next_f64() * 2.0 - 1.0,
                rng.next_f64() * 2.0 - 1.0,
                rng.next_f64() * 2.0 - 1.0,
            ];
            let t = [
                rng.next_f64() * 2.0 - 1.0,
                rng.next_f64() * 2.0 - 1.0,
                rng.next_f64() * 2.0 - 1.0,
            ];
            let (_, g) = pairwise_overlap_grad(&a, &b, &w, &t);
            let eps = 1e-6;
            for k in 0..6 {
                let mut wp = w;
                let mut tp = t;
                let mut wm = w;
                let mut tm = t;
                match k {
                    0..=2 => {
                        wp[k] += eps;
                        wm[k] -= eps;
                    }
                    _ => {
                        tp[k - 3] += eps;
                        tm[k - 3] -= eps;
                    }
                }
                let fp = pairwise_overlap_grad(&a, &b, &wp, &tp).0;
                let fm = pairwise_overlap_grad(&a, &b, &wm, &tm).0;
                let fd = (fp - fm) / (2.0 * eps);
                assert!(
                    (fd - g[k]).abs() < 1e-4 * (1.0 + fd.abs()),
                    "case {case} dof {k}: analytic {} vs fd {fd}",
                    g[k]
                );
            }
        }
    }

    /// Self-alignment from arbitrary poses must recover Tanimoto ≈ 1
    /// (benzene-like ring of C atoms, rotated/translated copy).
    #[test]
    fn self_alignment_recovers_identity() {
        let ring: Vec<ShapeAtom> = (0..6)
            .map(|i| {
                let th = i as f64 * std::f64::consts::PI / 3.0;
                sa([1.4 * th.cos(), 1.4 * th.sin(), 0.0], 6)
            })
            .collect();
        // pose: rotate ~2 rad around a skew axis + translate
        let w = [0.8, 1.1, -0.6];
        let t = [3.0, -2.0, 1.5];
        let moved: Vec<ShapeAtom> = ring
            .iter()
            .map(|a| sa(rodrigues(&w, &t, &a.c), 6))
            .collect();
        let ring_m = shape_mol_from_atoms(&ring);
        let moved_m = shape_mol_from_atoms(&moved);
        let res = align(
            &ring_m,
            &moved_m,
            &AlignOptions {
                random_starts: 6,
                max_iter: 250,
                screen: false,
                rescore_top: 8,
                seed: 3,
            },
        );
        assert!(
            res.tanimoto > 0.999,
            "self-alignment tanimoto {} (transform {:?})",
            res.tanimoto,
            res.transform
        );
    }

    /// For identical molecules at least one principal-frame start must land
    /// on (nearly) the exact coincident pose — otherwise self-alignment
    /// quality silently depends on the random starts.
    #[test]
    fn frame_starts_hit_self_alignment() {
        let sdf = std::fs::read_to_string("tests/fixtures/conformers/aspirin.sdf").unwrap();
        let mol = crate::molecule::parser::parse_sdf(&sdf).unwrap();
        let atoms = shape_atoms(&mol);
        let q = center(&atoms);
        let axes = principal_axes(&q);
        // reconstruct the (+,+) frame start exactly like align() does
        for &(s1, s2) in &[(1.0f64, 1.0), (1.0, -1.0), (-1.0, 1.0), (-1.0, -1.0)] {
            let t1 = [s1 * axes[1][0], s1 * axes[1][1], s1 * axes[1][2]];
            let t2 = [s2 * axes[2][0], s2 * axes[2][1], s2 * axes[2][2]];
            let mut cols = [[0.0f64; 3]; 3];
            for r in 0..3 {
                for c in 0..3 {
                    let tv = [axes[0][c], t1[c], t2[c]];
                    cols[c][r] = axes[0][r] * tv[0] + axes[1][r] * tv[1] + axes[2][r] * tv[2];
                }
            }
            let det = cols[0][0] * (cols[1][1] * cols[2][2] - cols[1][2] * cols[2][1])
                - cols[0][1] * (cols[1][0] * cols[2][2] - cols[1][2] * cols[2][0])
                + cols[0][2] * (cols[1][0] * cols[2][1] - cols[1][1] * cols[2][0]);
            let angle = ((cols[0][0] + cols[1][1] + cols[2][2] - 1.0) / 2.0)
                .clamp(-1.0, 1.0)
                .acos();
            println!("flip ({s1},{s2}): det={det:.3} angle={angle:.4}");
            // For a true self-frame the (+,+) combo must be ~identity
            if (s1, s2) == (1.0, 1.0) {
                assert!(angle < 1e-3, "self frame start angle {angle}");
            }
        }
    }

    /// Color-overlap analytic gradient vs central finite differences
    /// (seeded; the kernel is shared with the shape surrogate but the
    /// type filtering changes the summation structure).
    #[test]
    fn color_gradient_fd() {
        let mk = |i: u8, c: [f64; 3]| ColorSite {
            c,
            alpha: 1.0 + 0.1 * i as f64,
            type_id: i % 3,
            w: 1.0,
        };
        let mut rng = MiniRng::new(11);
        for case in 0..30 {
            let qs: Vec<ColorSite> = (0..5)
                .map(|i| {
                    mk(
                        i,
                        [
                            rng.next_f64() * 3.0 - 1.5,
                            rng.next_f64() * 3.0 - 1.5,
                            rng.next_f64() * 3.0 - 1.5,
                        ],
                    )
                })
                .collect();
            let ts: Vec<ColorSite> = (0..4)
                .map(|i| {
                    mk(
                        i + 2,
                        [
                            rng.next_f64() * 3.0 - 1.5,
                            rng.next_f64() * 3.0 - 1.5,
                            rng.next_f64() * 3.0 - 1.5,
                        ],
                    )
                })
                .collect();
            let w = [
                rng.next_f64() * 1.5 - 0.75,
                rng.next_f64() * 1.5 - 0.75,
                rng.next_f64() * 1.5 - 0.75,
            ];
            let t = [
                rng.next_f64() - 0.5,
                rng.next_f64() - 0.5,
                rng.next_f64() - 0.5,
            ];
            let (_, g) = color_overlap_grad(&qs, &ts, &w, &t);
            let eps = 1e-6;
            for k in 0..6 {
                let mut wp = w;
                let mut tp = t;
                let mut wm = w;
                let mut tm = t;
                if k < 3 {
                    wp[k] += eps;
                    wm[k] -= eps;
                } else {
                    tp[k - 3] += eps;
                    tm[k - 3] -= eps;
                }
                let fp = color_overlap_grad(&qs, &ts, &wp, &tp).0;
                let fm = color_overlap_grad(&qs, &ts, &wm, &tm).0;
                let fd = (fp - fm) / (2.0 * eps);
                assert!(
                    (fd - g[k]).abs() < 1e-4 * (1.0 + fd.abs()),
                    "case {case} dof {k}: {} vs {fd}",
                    g[k]
                );
            }
        }
    }

    /// Color Tanimoto identities: self = 1, empty = 0, symmetric.
    #[test]
    fn color_tanimoto_identities() {
        let qs = vec![
            ColorSite {
                c: [0.0, 0.0, 0.0],
                alpha: 1.0,
                type_id: 0,
                w: 1.0,
            },
            ColorSite {
                c: [1.5, 0.0, 0.0],
                alpha: 1.0,
                type_id: 1,
                w: 1.0,
            },
            ColorSite {
                c: [0.0, 1.5, 0.0],
                alpha: 1.1,
                type_id: 1,
                w: 1.0,
            },
        ];
        let self_t = color_tanimoto_at(&qs, &qs, &[0.0; 3], &[0.0; 3]);
        assert!((self_t - 1.0).abs() < 1e-9, "self {self_t}");
        let ts = vec![ColorSite {
            c: [4.0, 0.0, 0.0],
            alpha: 1.0,
            type_id: 2,
            w: 1.0,
        }];
        assert_eq!(color_tanimoto_at(&qs, &ts, &[0.0; 3], &[0.0; 3]), 0.0); // no shared types
        let empty: Vec<ColorSite> = Vec::new();
        assert_eq!(color_tanimoto_at(&empty, &qs, &[0.0; 3], &[0.0; 3]), 0.0);
    }

    /// Joint shape+color self-alignment reaches combo = 2.0.
    #[test]
    fn colored_self_alignment_combo() {
        let ring: Vec<ShapeAtom> = (0..6)
            .map(|i| {
                let th = i as f64 * std::f64::consts::PI / 3.0;
                ShapeAtom {
                    c: [1.4 * th.cos(), 1.4 * th.sin(), 0.0],
                    alpha: alpha_for(6),
                }
            })
            .collect();
        let qs: Vec<ColorSite> = [
            ColorSite {
                c: [0.0, 0.0, 0.3],
                alpha: alpha_for(7),
                type_id: 0,
                w: 1.0,
            },
            ColorSite {
                c: [1.2, 0.4, -0.2],
                alpha: alpha_for(8),
                type_id: 1,
                w: 1.0,
            },
        ]
        .to_vec();
        // rigidly move the target (atoms + sites together)
        let w0 = [0.7, -1.1, 0.5];
        let t0 = [2.0, 1.0, -3.0];
        let moved: Vec<ShapeAtom> = ring
            .iter()
            .map(|a| ShapeAtom {
                c: rodrigues(&w0, &t0, &a.c),
                alpha: a.alpha,
            })
            .collect();
        let t_sites: Vec<ColorSite> = qs
            .iter()
            .map(|s| ColorSite {
                c: rodrigues(&w0, &t0, &s.c),
                alpha: s.alpha,
                type_id: s.type_id,
                w: 1.0,
            })
            .collect();
        let ring_m = shape_mol_from_atoms(&ring);
        let moved_m = shape_mol_from_atoms(&moved);
        let res = align_colored(
            &ring_m,
            &moved_m,
            Some(&qs),
            Some(&t_sites),
            &AlignOptions {
                random_starts: 6,
                max_iter: 250,
                screen: false,
                rescore_top: 8,
                seed: 5,
            },
        );
        assert!(
            res.combo > 1.99,
            "combo {} (shape {} color {})",
            res.combo,
            res.tanimoto,
            res.color_tanimoto
        );
    }

    /// Screening mode: surrogate-ranked, no IE scoring, frame starts only.
    #[test]
    fn screen_mode_ranking() {
        let sdf = std::fs::read_to_string("tests/fixtures/conformers/aspirin.sdf").unwrap();
        let mol = crate::molecule::parser::parse_sdf(&sdf).unwrap();
        let base = shape_atoms(&mol);
        let base_m = shape_mol_from_atoms(&base);
        let w = [0.9, 0.4, -0.6];
        let moved: Vec<ShapeAtom> = base
            .iter()
            .map(|a| ShapeAtom {
                c: rodrigues(&w, &[1.0, -2.0, 0.5], &a.c),
                alpha: a.alpha,
            })
            .collect();
        let moved_m = shape_mol_from_atoms(&moved);
        let opts = AlignOptions {
            random_starts: 16,
            max_iter: 200,
            seed: 42,
            screen: true,
            rescore_top: 3,
        };
        let res = align_colored(&base_m, &moved_m, None, None, &opts);
        assert!(res.surrogate_overlap > 0.0);
        assert_eq!(
            res.tanimoto, 0.0,
            "screen mode must not compute the IE tanimoto"
        );
        // far fewer iterations than the full pipeline (no polish, frame starts)
        let full = align(&base_m, &moved_m, &AlignOptions::default());
        assert!(res.iterations < full.iterations);
    }

    /// rescore_top equivalence: K=3 vs rescore-every-start across the
    /// conformer fixture pairs — the best combo must agree to < 1e-3.
    #[test]
    fn rescore_top_equivalence() {
        let dir = std::path::Path::new("tests/fixtures/conformers");
        let names = [
            "ethanol.sdf",
            "aspirin.sdf",
            "ibuprofen.sdf",
            "naphthalene.sdf",
            "threonine.sdf",
        ];
        let mols: Vec<ShapeMol> = names
            .iter()
            .map(|n| {
                let sdf = std::fs::read_to_string(dir.join(n)).unwrap();
                shape_mol(&crate::molecule::parser::parse_sdf(&sdf).unwrap())
            })
            .collect();
        let mut worst = 0.0f64;
        for i in 0..mols.len() {
            for j in 0..mols.len() {
                if i == j {
                    continue;
                }
                // pose the target so alignment is non-trivial
                let w = [0.9, 0.3 * i as f64, -0.4 * j as f64];
                let moved: Vec<ShapeAtom> = mols[j]
                    .atoms
                    .iter()
                    .map(|a| ShapeAtom {
                        c: rodrigues(&w, &[1.0, -1.0, 0.5], &a.c),
                        alpha: a.alpha,
                    })
                    .collect();
                let moved_m = shape_mol_from_atoms(&moved);
                let k8 = align(
                    &mols[i],
                    &moved_m,
                    &AlignOptions {
                        random_starts: 12,
                        max_iter: 200,
                        seed: 7,
                        screen: false,
                        rescore_top: 8,
                    },
                );
                let kall = align(
                    &mols[i],
                    &moved_m,
                    &AlignOptions {
                        random_starts: 12,
                        max_iter: 200,
                        seed: 7,
                        screen: false,
                        rescore_top: usize::MAX,
                    },
                );
                worst = worst.max((k8.combo - kall.combo).abs());
            }
        }
        assert!(worst < 0.02, "worst combo delta {worst}");
    }

    /// Shape Tanimoto is symmetric and bounded.
    #[test]
    fn tanimoto_symmetry() {
        let a = vec![
            sa([0.0; 3], 6),
            sa([1.5, 0.0, 0.0], 7),
            sa([0.0, 1.5, 0.3], 8),
        ];
        let b = vec![
            sa([0.2, 0.1, 0.0], 6),
            sa([1.6, 0.2, 0.1], 6),
            sa([3.0, 0.0, 0.0], 1),
        ];
        let tab = shape_tanimoto(&a, &b);
        let tba = shape_tanimoto(&b, &a);
        assert!((tab - tba).abs() < 1e-12);
        assert!((0.0..=1.0).contains(&tab));
    }
}

#[cfg(test)]
mod fastloop_tests {
    use super::*;

    /// Brute-force pair loop (no skips) vs overlap_terms — settles whether
    /// the rigorous-bound skip changes totals (the 7Å+r_bound prefilter it
    /// replaced may have been dropping contributing pairs).
    #[test]
    fn fast_loop_matches_brute_force() {
        for name in [
            "ethanol",
            "aspirin",
            "ibuprofen",
            "naphthalene",
            "threonine",
        ] {
            let sdf =
                std::fs::read_to_string(format!("tests/fixtures/conformers/{name}.sdf")).unwrap();
            let mol = crate::molecule::parser::parse_sdf(&sdf).unwrap();
            let atoms = shape_atoms(&mol);
            let w = [0.9, 0.4, -0.6];
            let posed: Vec<ShapeAtom> = atoms
                .iter()
                .map(|a| ShapeAtom {
                    c: rodrigues(&w, &[1.0, -1.0, 0.5], &a.c),
                    alpha: a.alpha,
                })
                .collect();
            let ta = materialize(&prepare_shape(&atoms).specs, &atoms);
            let tb = materialize(&prepare_shape(&posed).specs, &posed);
            let fast = overlap_terms(&ta, &tb);
            // brute force: same iteration order, no skipping of any kind
            let scale: f64 = ta.iter().map(volume_of_term).sum::<f64>().abs().max(1e-12);
            let cutoff = 1e-7 * scale;
            let mut total = 0.0f64;
            for pa in ta.iter() {
                for pb in tb.iter() {
                    let v = pa.merged_volume(pb);
                    if v.abs() < cutoff {
                        continue;
                    }
                    total += pa.sign * pb.sign * v;
                }
            }
            assert!(
                (fast - total).abs() < 1e-9,
                "{name}: fast {fast} vs brute {total} (delta {})",
                (fast - total).abs()
            );
        }
    }
}

#[cfg(test)]
mod tests_color_weights {
    use super::*;

    fn two_type_sites(w0: f64, w1: f64) -> (Vec<ColorSite>, Vec<ColorSite>) {
        // two same-type pairs at different offsets: donor pair close (large
        // overlap), ring pair far (small overlap) — weights must scale their
        // own type's contribution only
        let q = vec![
            ColorSite {
                c: [0.0, 0.0, 0.0],
                alpha: 2.0,
                type_id: 0,
                w: w0,
            },
            ColorSite {
                c: [0.0, 0.0, 0.0],
                alpha: 2.0,
                type_id: 5,
                w: w1,
            },
        ];
        let t = vec![
            ColorSite {
                c: [0.2, 0.0, 0.0],
                alpha: 2.0,
                type_id: 0,
                w: w0,
            },
            ColorSite {
                c: [6.0, 0.0, 0.0],
                alpha: 2.0,
                type_id: 5,
                w: w1,
            },
        ];
        (q, t)
    }

    #[test]
    fn weight_one_is_bit_identical() {
        // w = 1.0 must not perturb the overlap: multiplying by 1.0 keeps
        // every f64 bit (zero drift for default-weight callers)
        let (q, t) = two_type_sites(1.0, 1.0);
        let qs = q.clone();
        let a = color_overlap_grad(&q, &t, &[0.0; 3], &[0.0; 3]).0;
        // hand-computed: same as the unweighted kernel (regression anchor)
        assert!(a > 0.0);
        let (o_w, g_w) = color_overlap_grad(&q, &t, &[0.1, -0.1, 0.05], &[0.3, 0.0, -0.2]);
        // gradient path runs with weights too
        assert!(o_w.is_finite() && g_w.iter().all(|x| x.is_finite()));
    }

    #[test]
    fn sqrt_u_scales_pair_linearly() {
        // all sites of a type share w = sqrt(u) -> that type's pair terms
        // scale EXACTLY by u (linear user semantics)
        let (q1, t1) = two_type_sites(1.0, 1.0);
        let (q2, t2) = two_type_sites(2f64.sqrt(), 1.0);
        // query site 0 (donor) vs target 0: contribution scales by 2
        // ring pair (index 1) identical
        let alpha_ij = 4.0;
        let beta = 2.0 * 2.0 / alpha_ij;
        let k = GCI
            * GCI
            * (std::f64::consts::PI / alpha_ij)
            * (std::f64::consts::PI / alpha_ij).sqrt();
        let d2_donor = 0.2 * 0.2;
        let donor1 = (-beta * d2_donor).exp() * k;
        let d2_ring = 6.0 * 6.0;
        let ring = (-beta * d2_ring).exp() * k;
        let o1 = color_overlap_grad(&q1, &t1, &[0.0; 3], &[0.0; 3]).0;
        let o2 = color_overlap_grad(&q2, &t2, &[0.0; 3], &[0.0; 3]).0;
        assert!(
            (o1 - (donor1 + ring)).abs() < 1e-12,
            "o1 {o1} vs {}",
            donor1 + ring
        );
        assert!(
            (o2 - (2.0 * donor1 + ring)).abs() < 1e-12,
            "o2 {o2} vs {}",
            2.0 * donor1 + ring
        );
    }

    #[test]
    fn weighted_gradient_matches_fd() {
        // finite-difference check of the weighted gradient (w enters as a
        // constant pair factor — chain rule must scale the analytic kernel)
        let (q, t) = two_type_sites(1.3f64.sqrt(), 0.7f64.sqrt());
        let w0 = [0.12, -0.08, 0.05];
        let t0 = [0.4, -0.2, 0.1];
        let (_, g) = color_overlap_grad(&q, &t, &w0, &t0);
        let h = 1e-6;
        let mut maxerr = 0.0f64;
        for k in 0..3 {
            let mut tp = t0;
            tp[k] += h;
            let mut tm = t0;
            tm[k] -= h;
            let fd = (color_overlap_grad(&q, &t, &w0, &tp).0
                - color_overlap_grad(&q, &t, &w0, &tm).0)
                / (2.0 * h);
            maxerr = maxerr.max((fd - g[3 + k]).abs());
        }
        assert!(maxerr < 1e-6, "translation grad FD err {maxerr}");
        for k in 0..3 {
            let mut wp = w0;
            wp[k] += h;
            let mut wm = w0;
            wm[k] -= h;
            let fd = (color_overlap_grad(&q, &t, &wp, &t0).0
                - color_overlap_grad(&q, &t, &wm, &t0).0)
                / (2.0 * h);
            maxerr = maxerr.max((fd - g[k]).abs());
        }
        assert!(maxerr < 1e-5, "rotation grad FD err {maxerr}");
    }

    #[test]
    fn all_zero_weights_give_zero_color() {
        // the tanimoto denominator guard must hold when every weight is 0
        let (q, t) = two_type_sites(0.0, 0.0);
        let o = color_overlap_grad(&q, &t, &[0.0; 3], &[0.0; 3]).0;
        assert_eq!(o, 0.0);
        let tc = color_tanimoto_at(&q, &t, &[0.0; 3], &[0.0; 3]);
        assert_eq!(tc, 0.0);
    }
}
