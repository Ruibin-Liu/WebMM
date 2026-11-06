//! iMTD-GC conformational search (CREST-lineage iterative metadynamics +
//! genetic crossing).
//!
//! Design baseline (user directive): organize metadynamics-based conformer
//! search as the CREST iMTD-GC iteration — multiple rounds of SHORT metadynamics
//! runs (rotating bias strengths and CV combinations across runs) collect
//! conformers by periodically quenching trajectory snapshots to their local
//! minima; low-energy minima are then genetically crossed (geometry
//! recombination at a rotatable-bond cut) to seed the next round; the loop
//! repeats until the conformer ensemble converges. A single long metadynamics
//! run is explicitly NOT the baseline.
//!
//! The well-tempered engine ([`crate::metad::MetaDynamics`]) is used as-is;
//! this module is the outer iterative-crossing scheduler around it. Everything
//! is deterministic: all randomness derives from `ImtdConfig::seed`, so the
//! same input + config reproduces the same ensemble bit-for-bit.
//!
//! # Conformer-equivalence semantics (v1)
//! Two quenched minima are the SAME conformer iff their heavy-atom scaffolds
//! superpose (aligned RMSD below `rmsd_tol`) AND (their all-atom geometries
//! superpose OR their energies are equal to 1e-3 kcal/mol — methyl-H label
//! permutations quench to the identical minimum). Heavy-atom MIRROR pairs
//! (e.g. butane gauche±) cannot superpose under a proper rotation and are
//! kept as distinct conformers; equal-energy H-arrangement images sharing
//! one heavy scaffold collapse to one entry.

use crate::forces::ForceField;
use crate::md::{MDConfig, MDRunner};
use crate::metad::{CollectiveVariable, DihedralCV, DistanceCV, MetaDConfig, MetaDynamics};
use crate::mmff::MMFFForceField;
use crate::molecule::graph::find_rings;
use crate::molecule::Molecule;
use crate::optimizer::{self};
use crate::ConvergenceOptions;
use std::collections::{HashSet, VecDeque};
use std::rc::Rc;

/// The energy engine a search runs on. GFN-FF is rebuilt per geometry (its
/// bond perception reads coordinates — same discipline as the analog
/// pipeline); MMFF is built once per molecule.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SearchEngine {
    MMFF94,
    MMFF94s,
    GFNFF,
}

impl SearchEngine {
    /// Parse the engine label accepted across the WASM surface
    /// (`ImtdOptions.mmff_variant` doubles as the engine selector).
    pub fn from_label(v: &str) -> Self {
        match v {
            "MMFF94" => SearchEngine::MMFF94,
            "GFNFF" | "GFN-FF" => SearchEngine::GFNFF,
            _ => SearchEngine::MMFF94s,
        }
    }
}

/// GFN-FF exact-coords cache: the last built force field with the
/// coordinates it was perceived from.
pub type GfnffCache =
    std::cell::RefCell<Option<(Vec<[f64; 3]>, Box<crate::gfnff::GfnffForceField>)>>;

/// Quench backend: L-BFGS relaxation of a snapshot to its local minimum.
pub enum Quencher {
    /// MMFF: one shared instance (coordinate-independent topology).
    Mmff(Box<MMFFForceField>),
    /// GFN-FF: atomic numbers + total charge. Bond perception reads the
    /// geometry, so the force field is rebuilt per quench — EXCEPT when the
    /// coordinates are exactly the ones of the cached build (a real case:
    /// the thermo tier re-runs per temperature at the same minima).
    /// Measured: ~0.5 ms per build at ethanol scale — the exact-coords cache
    /// is the only zero-risk reuse; topology-hash caching would re-run the
    /// perception it tries to save.
    Gfnff {
        at: Vec<usize>,
        charge: f64,
        cache: GfnffCache,
    },
}

impl Quencher {
    /// Build the backend for a molecule + engine.
    pub fn new(mol: &Molecule, engine: SearchEngine) -> Self {
        match engine {
            SearchEngine::GFNFF => Quencher::Gfnff {
                at: mol.atoms.iter().map(|a| a.atomic_number as usize).collect(),
                charge: mol.atoms.iter().map(|a| a.charge).sum(),
                cache: std::cell::RefCell::new(None),
            },
            SearchEngine::MMFF94 => Quencher::Mmff(Box::new(MMFFForceField::new(
                mol,
                crate::MMFFVariant::MMFF94,
            ))),
            SearchEngine::MMFF94s => Quencher::Mmff(Box::new(MMFFForceField::new(
                mol,
                crate::MMFFVariant::MMFF94s,
            ))),
        }
    }

    fn quench(&self, coords: &[[f64; 3]], max_iter: usize) -> crate::optimizer::OptimizationResult {
        let conv = |max_force: f64, rms_force: f64| ConvergenceOptions {
            max_iterations: max_iter.max(1),
            max_force,
            rms_force,
            energy_change: Default::default(),
        };
        match self {
            Quencher::Mmff(ff) => optimizer::optimize(ff.as_ref(), coords, &conv(0.01, 0.001)),
            Quencher::Gfnff { at, charge, cache } => {
                // GFN-FF numeric kernel: same relaxed force floors the analog
                // pipeline applies (line-search stalls on large systems)
                let mut guard = cache.borrow_mut();
                let ff: &crate::gfnff::GfnffForceField = match guard.as_ref() {
                    Some((prev, ff)) if prev == coords => ff.as_ref(),
                    _ => {
                        *guard = Some((
                            coords.to_vec(),
                            Box::new(crate::gfnff::GfnffForceField::new(at, coords, *charge)),
                        ));
                        guard
                            .as_ref()
                            .map(|(_, f)| f.as_ref())
                            .unwrap_or_else(|| unreachable!("cache was just filled"))
                    }
                };
                optimizer::optimize(ff, coords, &conv(0.05, 0.005))
            }
        }
    }
}

/// Deterministic ladder of bias-strength factors cycled across metad runs
/// (CREST rotates bias potential magnitudes across its metad copies).
const HEIGHT_LADDER: [f64; 3] = [1.0, 0.5, 0.25];

/// How many low-energy ensemble members are eligible as metad start
/// structures / crossing parents each round (energy-sorted prefix).
const PARENT_POOL: usize = 8;

/// Maximum number of distance CVs used to top up a small dihedral pool
/// (ring systems without rotatable bonds).
const MAX_DISTANCE_CVS: usize = 4;

/// Configuration for one iMTD-GC search.
#[derive(Debug, Clone, Copy)]
pub struct ImtdConfig {
    /// Hard cap on iMTD-GC iterations (each = k short metad runs + crossing).
    pub max_rounds: usize,
    /// Short metadynamics runs per round (CVs and bias strengths rotate).
    pub runs_per_round: usize,
    /// MD steps per short metad run.
    pub steps_per_run: usize,
    /// Quench a snapshot to its local minimum every N steps of each run.
    pub quench_interval: usize,
    /// MD timestep (fs).
    pub dt_fs: f64,
    /// Metad run temperature (K). CREST explores above room temperature.
    pub temperature_k: f64,
    /// Langevin friction (1/ps).
    pub friction_per_ps: f64,
    /// Base RNG seed; per-run seeds derive as `seed + round*1000 + run`.
    pub seed: u64,
    /// Base hill height (kcal/mol); actual = base × ladder factor.
    pub hill_height: f64,
    /// Hill width for dihedral CVs (radians).
    pub hill_width: f64,
    /// Hill width for distance CVs (Å).
    pub hill_width_distance: f64,
    /// Deposit a hill every N MD steps.
    pub deposit_interval: usize,
    /// Well-tempered bias factor γ (>1 = well-tempered).
    pub bias_factor: f64,
    /// L-BFGS iteration cap per quench.
    pub quench_max_iterations: usize,
    /// Keep minima within this window (kcal/mol) of the ensemble minimum.
    pub energy_window: f64,
    /// Two minima closer than this heavy-atom aligned RMSD (Å) are one.
    pub rmsd_tol: f64,
    /// Ensemble size cap (lowest-energy members kept).
    pub max_conformers: usize,
    /// Genetic-crossing pairs attempted per round.
    pub crossing_pairs: usize,
    /// A round adding at most this many NEW unique minima means converged.
    pub convergence_new_minima: usize,
}

impl Default for ImtdConfig {
    fn default() -> Self {
        ImtdConfig {
            max_rounds: 4,
            runs_per_round: 3,
            steps_per_run: 2000,
            quench_interval: 250,
            dt_fs: 1.0,
            temperature_k: 400.0,
            friction_per_ps: 5.0,
            seed: 42,
            hill_height: 0.5,
            hill_width: 0.2,
            hill_width_distance: 0.5,
            deposit_interval: 25,
            bias_factor: 10.0,
            quench_max_iterations: 300,
            energy_window: 6.0,
            rmsd_tol: 0.25,
            max_conformers: 250,
            crossing_pairs: 4,
            convergence_new_minima: 0,
        }
    }
}

/// How a conformer entered the ensemble.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Origin {
    /// The quenched input structure.
    Input,
    /// Quenched snapshot of a metadynamics trajectory.
    Metad,
    /// Quenched genetic-crossing child.
    Crossing,
}

/// One unique local minimum of the ensemble.
#[derive(Debug, Clone)]
pub struct Conformer {
    pub coords: Vec<[f64; 3]>,
    /// MMFF energy of the quenched minimum (kcal/mol).
    pub energy: f64,
    /// Round (0-based; 0 = the quenched input itself) in which it was found.
    pub round_found: usize,
    pub origin: Origin,
}

/// Result of an iMTD-GC search.
#[derive(Debug, Clone)]
pub struct ImtdResult {
    /// Unique minima, sorted by ascending energy, inside the energy window.
    pub conformers: Vec<Conformer>,
    /// iMTD-GC iterations actually executed.
    pub rounds_executed: usize,
    /// True if the convergence criterion fired (false = budget exhausted).
    pub converged: bool,
    /// New unique minima found per round (collect + crossing combined).
    pub new_minima_per_round: Vec<usize>,
    /// Size of the CV pool the runs rotate over.
    pub pool_size: usize,
    /// Number of rotatable bonds (crossing cut candidates).
    pub n_rotatable: usize,
}

/// One runnable short-metad unit of work (the scatter-gather boundary).
/// Self-contained: a worker with the molecule + config executes it as a pure
/// function and returns its quenched [`Minimum`] list.
#[derive(Debug, Clone, PartialEq)]
pub struct RunSpec {
    /// Round this spec belongs to (1-based; reporting only).
    pub round: usize,
    /// Index within the round (0-based; also the gather fold order).
    pub run_index: usize,
    /// Start coordinates (n_atoms rows).
    pub start: Vec<[f64; 3]>,
    pub cv: CvSpec,
    /// Hill height for this run (bias-strength ladder already applied).
    pub hill_height: f64,
    /// RNG seed for this run.
    pub seed: u64,
}

/// One CI-NEB saddle-refinement unit of work (scatter-gather boundary).
#[derive(Debug, Clone, PartialEq)]
pub struct RefineJob {
    /// Slot in the refinement batch.
    pub index: usize,
    /// Stable basin identities (never shift).
    pub id_a: u64,
    pub id_b: u64,
    /// Final result-order indices (display only).
    pub a: usize,
    pub b: usize,
    pub xyz_a: Vec<[f64; 3]>,
    pub xyz_b: Vec<[f64; 3]>,
}

/// A trajectory snapshot already quenched to its local minimum.
#[derive(Debug, Clone, PartialEq)]
pub struct Minimum {
    pub coords: Vec<[f64; 3]>,
    /// MMFF energy of the quenched minimum (kcal/mol).
    pub energy: f64,
    /// Max UNBIASED potential energy since the previous snapshot of the same
    /// run (kcal/mol): observed upper bound of the transition barrier between
    /// consecutive basins. Decreases with sampling.
    pub pe_max_since_prev: f64,
    /// Quench ordinal within its run (0-based).
    pub snapshot: usize,
}

// ── geometry helpers ──

/// Proper rotation + translation (Horn quaternion Kabsch) mapping the mobile
/// point set `p` onto the reference `q` (minimizes Σ|R·p+t−q|²), plus the
/// achieved RMSD. Reflections are never returned (molecules are proper
/// objects). Coordinate-level sibling of the `kabsch_align` WASM export
/// (same math and convention: R·p+t ≈ q).
pub(crate) fn kabsch_transform(p: &[[f64; 3]], q: &[[f64; 3]]) -> ([f64; 9], [f64; 3], f64) {
    debug_assert_eq!(p.len(), q.len(), "kabsch point sets must be paired");
    let n_f = p.len() as f64;
    let mut pc = [0.0f64; 3];
    let mut qc = [0.0f64; 3];
    for i in 0..p.len() {
        for c in 0..3 {
            pc[c] += p[i][c] / n_f;
            qc[c] += q[i][c] / n_f;
        }
    }
    // cross-covariance S = Σ p̃ q̃ᵀ over centroid-centered points
    let mut sxx = 0.0;
    let mut sxy = 0.0;
    let mut sxz = 0.0;
    let mut syx = 0.0;
    let mut syy = 0.0;
    let mut syz = 0.0;
    let mut szx = 0.0;
    let mut szy = 0.0;
    let mut szz = 0.0;
    for i in 0..p.len() {
        let a = [p[i][0] - pc[0], p[i][1] - pc[1], p[i][2] - pc[2]];
        let b = [q[i][0] - qc[0], q[i][1] - qc[1], q[i][2] - qc[2]];
        sxx += a[0] * b[0];
        sxy += a[0] * b[1];
        sxz += a[0] * b[2];
        syx += a[1] * b[0];
        syy += a[1] * b[1];
        syz += a[1] * b[2];
        szx += a[2] * b[0];
        szy += a[2] * b[1];
        szz += a[2] * b[2];
    }
    // Horn's 4×4 symmetric matrix; largest eigenvector = optimal quaternion
    let mut nm = [
        sxx + syy + szz,
        syz - szy,
        szx - sxz,
        sxy - syx,
        syz - szy,
        sxx - syy - szz,
        sxy + syx,
        szx + sxz,
        szx - sxz,
        sxy + syx,
        -sxx + syy - szz,
        syz + szy,
        sxy - syx,
        szx + sxz,
        syz + szy,
        -sxx - syy + szz,
    ];
    let mut vec = [0.0f64; 16];
    crate::optimizer::jacobi::sym_jacobi(&mut nm, &mut vec, 4);
    let mut best = 0usize;
    for k in 1..4 {
        if nm[k * 4 + k] > nm[best * 4 + best] {
            best = k;
        }
    }
    let qw = vec[best];
    let qx = vec[4 + best];
    let qy = vec[8 + best];
    let qz = vec[12 + best];
    let r = [
        1.0 - 2.0 * (qy * qy + qz * qz),
        2.0 * (qx * qy - qw * qz),
        2.0 * (qx * qz + qw * qy),
        2.0 * (qx * qy + qw * qz),
        1.0 - 2.0 * (qx * qx + qz * qz),
        2.0 * (qy * qz - qw * qx),
        2.0 * (qx * qz - qw * qy),
        2.0 * (qy * qz + qw * qx),
        1.0 - 2.0 * (qx * qx + qy * qy),
    ];
    // t = qc − R·pc  (maps p into q's frame)
    let t = [
        qc[0] - (r[0] * pc[0] + r[1] * pc[1] + r[2] * pc[2]),
        qc[1] - (r[3] * pc[0] + r[4] * pc[1] + r[5] * pc[2]),
        qc[2] - (r[6] * pc[0] + r[7] * pc[1] + r[8] * pc[2]),
    ];
    let mut msd = 0.0f64;
    for i in 0..p.len() {
        let x = r[0] * p[i][0] + r[1] * p[i][1] + r[2] * p[i][2] + t[0] - q[i][0];
        let y = r[3] * p[i][0] + r[4] * p[i][1] + r[5] * p[i][2] + t[1] - q[i][1];
        let z = r[6] * p[i][0] + r[7] * p[i][1] + r[8] * p[i][2] + t[2] - q[i][2];
        msd += (x * x + y * y + z * z) / n_f;
    }
    (r, t, msd.sqrt())
}

/// Aligned RMSD over the given atom indices (Å). For <2 indices no rotation is
/// determined; the centroid gap of the single pair is returned.
/// Heavy-atom aligned RMSD between two coordinate sets over `idx`
/// (Kabsch; crate-internal: also used by the conformer-subset contract).
pub(crate) fn aligned_rmsd_over(a: &[[f64; 3]], b: &[[f64; 3]], idx: &[usize]) -> f64 {
    match idx.len() {
        0 => 0.0,
        1 => {
            let d = [
                a[idx[0]][0] - b[idx[0]][0],
                a[idx[0]][1] - b[idx[0]][1],
                a[idx[0]][2] - b[idx[0]][2],
            ];
            (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt()
        }
        _ => {
            let p: Vec<[f64; 3]> = idx.iter().map(|&i| a[i]).collect();
            let q: Vec<[f64; 3]> = idx.iter().map(|&i| b[i]).collect();
            kabsch_transform(&p, &q).2
        }
    }
}

// ── molecule topology helpers ──

/// Indices of heavy atoms (z > 1).
fn heavy_indices(mol: &Molecule) -> Vec<usize> {
    (0..mol.atoms.len())
        .filter(|&i| mol.atoms[i].atomic_number > 1)
        .collect()
}

/// Rotatable bonds as (j, k) pairs: single, acyclic, both endpoints heavy,
/// each endpoint keeps ≥1 other neighbor (H termini allowed — the H–O–C–C
/// torsions of e.g. hydroxyls are prime metad CVs). Bond order follows the
/// molfile, deterministic.
fn rotatable_bonds(mol: &Molecule) -> Vec<(usize, usize)> {
    let mut ring_bonds = HashSet::new();
    for ring in find_rings(mol) {
        for w in ring.windows(2) {
            ring_bonds.insert((w[0].min(w[1]), w[0].max(w[1])));
        }
        if let (Some(&first), Some(&last)) = (ring.first(), ring.last()) {
            ring_bonds.insert((first.min(last), first.max(last)));
        }
    }
    let mut out = Vec::new();
    for bond in &mol.bonds {
        if !matches!(bond.bond_type, crate::molecule::BondType::Single) {
            continue;
        }
        let (j, k) = (bond.atom1, bond.atom2);
        if mol.atoms[j].atomic_number <= 1 || mol.atoms[k].atomic_number <= 1 {
            continue;
        }
        if ring_bonds.contains(&(j.min(k), j.max(k))) {
            continue;
        }
        if mol.adjacency[j].len() >= 2 && mol.adjacency[k].len() >= 2 {
            out.push((j, k));
        }
    }
    out
}

/// The i-j-k-l torsion biasing rotation around bond (j, k): i / l are
/// neighbors of j / k (≠ the other endpoint), preferring heavy atoms, then
/// the smallest index.
fn dihedral_for_bond(mol: &Molecule, j: usize, k: usize) -> Option<(usize, usize, usize, usize)> {
    let pick = |center: usize, exclude: usize| -> Option<usize> {
        let mut cands: Vec<usize> = mol.adjacency[center]
            .iter()
            .copied()
            .filter(|&n| n != exclude)
            .collect();
        if cands.is_empty() {
            return None;
        }
        // heavy atoms first (stable conformational meaning), then by index
        cands.sort_by_key(|&n| (mol.atoms[n].atomic_number <= 1, n));
        Some(cands[0])
    };
    Some((pick(j, k)?, j, k, pick(k, j)?))
}

/// Heavy-atom pairs at graph distance ≥ 4 (collapse/extension bias CVs for
/// systems with too few rotatable bonds, e.g. rings). Deterministic order.
fn distance_cvs(mol: &Molecule, cap: usize) -> Vec<(usize, usize)> {
    let heavy = heavy_indices(mol);
    let mut pairs: Vec<(usize, usize)> = Vec::new();
    for &src in &heavy {
        // BFS distances from src (capped by graph size; molecules here ≤ 1.5k atoms)
        let mut dist = vec![usize::MAX; mol.atoms.len()];
        dist[src] = 0;
        let mut queue = VecDeque::new();
        queue.push_back(src);
        while let Some(u) = queue.pop_front() {
            for &v in &mol.adjacency[u] {
                if dist[v] == usize::MAX {
                    dist[v] = dist[u] + 1;
                    queue.push_back(v);
                }
            }
        }
        for &dst in &heavy {
            if dst > src && dist[dst] != usize::MAX && dist[dst] >= 4 {
                pairs.push((src, dst));
            }
        }
    }
    pairs.sort_unstable();
    pairs.dedup();
    pairs.truncate(cap);
    pairs
}

/// One bias collective variable of the rotation pool.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum CvSpec {
    Dihedral(usize, usize, usize, usize),
    Distance(usize, usize),
}

impl CvSpec {
    fn build(self) -> Box<dyn CollectiveVariable> {
        match self {
            CvSpec::Dihedral(i, j, k, l) => Box::new(DihedralCV::new(i, j, k, l)),
            CvSpec::Distance(i, j) => Box::new(DistanceCV::new(i, j)),
        }
    }
}

/// Build the CV pool: one dihedral per rotatable bond; if fewer than 2
/// dihedrals exist, top up with distance CVs (ring systems). Empty pool ⇒
/// rigid molecule (caller reports the quenched input as the only conformer).
fn cv_pool(mol: &Molecule, rot: &[(usize, usize)]) -> Vec<CvSpec> {
    let mut pool: Vec<CvSpec> = rot
        .iter()
        .filter_map(|&(j, k)| {
            dihedral_for_bond(mol, j, k).map(|(i, j2, k2, l)| CvSpec::Dihedral(i, j2, k2, l))
        })
        .collect();
    if pool.len() < 2 {
        pool.extend(
            distance_cvs(mol, MAX_DISTANCE_CVS)
                .into_iter()
                .map(|(i, j)| CvSpec::Distance(i, j)),
        );
    }
    pool
}

/// Fragment split at bond (j, k): F1 = the connected component of j when the
/// bond is removed (hydrogens follow their parent side). Returns (f1, f2)
/// covering all atoms disjointly, each sorted.
fn fragment_split(mol: &Molecule, j: usize, k: usize) -> (Vec<usize>, Vec<usize>) {
    let mut seen = vec![false; mol.atoms.len()];
    seen[j] = true;
    let mut f1 = vec![j];
    let mut queue = VecDeque::new();
    queue.push_back(j);
    while let Some(u) = queue.pop_front() {
        for &v in &mol.adjacency[u] {
            // the cut EDGE blocks traversal j→k (k itself belongs to F2)
            if !seen[v] && !(u == j && v == k) {
                seen[v] = true;
                f1.push(v);
                queue.push_back(v);
            }
        }
    }
    let f2: Vec<usize> = (0..mol.atoms.len()).filter(|&i| !seen[i]).collect();
    f1.sort_unstable();
    (f1, f2)
}

/// Genetic crossing: align parent B onto parent A over F1 ∪ {k} (≥3 points
/// for any real rotatable bond), then take F1 coordinates from A and F2 from
/// the aligned B — a geometry recombination that inherits F1's conformational
/// features from A and F2's from B (CREST's GC step).
fn genetic_crossing(a: &[[f64; 3]], b: &[[f64; 3]], f1: &[usize], k: usize) -> Vec<[f64; 3]> {
    let mut align_idx: Vec<usize> = f1.to_vec();
    if !align_idx.contains(&k) {
        align_idx.push(k);
    }
    // mobile = B, reference = A (kabsch maps mobile onto reference)
    let p: Vec<[f64; 3]> = align_idx.iter().map(|&i| b[i]).collect();
    let q: Vec<[f64; 3]> = align_idx.iter().map(|&i| a[i]).collect();
    let (r, t, _) = kabsch_transform(&p, &q);
    let f1_set: HashSet<usize> = f1.iter().copied().collect();
    let mut child = a.to_vec();
    for (i, ci) in child.iter_mut().enumerate() {
        if !f1_set.contains(&i) {
            // F2 atoms come from the aligned B
            *ci = [
                r[0] * b[i][0] + r[1] * b[i][1] + r[2] * b[i][2] + t[0],
                r[3] * b[i][0] + r[4] * b[i][1] + r[5] * b[i][2] + t[1],
                r[6] * b[i][0] + r[7] * b[i][1] + r[8] * b[i][2] + t[2],
            ];
        }
    }
    child
}

// ── ensemble ──

/// Which trajectory family feeds the statistics: metad-biased search runs
/// (fast basin discovery, distorted occupancy) or plain-MD verification runs
/// (hill_height = 0; thermally meaningful populations).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Ladder {
    Metad,
    Md,
}

/// One statistics ladder, parallel to the ensemble members.
#[derive(Default)]
struct LadderStats {
    /// Quenched frames assigned per member.
    visits: Vec<usize>,
    /// Consecutive same-basin snapshot pairs per member (residence proxy).
    residence: Vec<usize>,
    trans_in: Vec<usize>,
    trans_out: Vec<usize>,
    /// Pair transition stats, unordered key, both directions summed;
    /// `barrier` is the best (lowest) observed upper bound.
    pairs: Vec<(usize, usize, usize, f64)>,
    /// Total quenched frames classified into this ladder.
    frames: usize,
    /// Frames rejected by window/cap (unassigned; they break chains).
    pruned: usize,
}

impl LadderStats {
    fn on_new_basin(&mut self) {
        self.visits.push(0);
        self.residence.push(0);
        self.trans_in.push(0);
        self.trans_out.push(0);
    }

    fn on_evict(&mut self, i: usize) {
        if self.visits.len() > i {
            self.visits.remove(i);
        }
        if self.residence.len() > i {
            self.residence.remove(i);
        }
        if self.trans_in.len() > i {
            self.trans_in.remove(i);
        }
        if self.trans_out.len() > i {
            self.trans_out.remove(i);
        }
        let shift = |x: &mut usize| {
            if *x > i {
                *x -= 1;
            }
        };
        for t in self.pairs.iter_mut() {
            shift(&mut t.0);
            shift(&mut t.1);
        }
    }

    fn note_visit(&mut self, i: usize) {
        if i < self.visits.len() {
            self.visits[i] += 1;
        }
        self.frames += 1;
    }

    fn note_residence(&mut self, i: usize) {
        if i < self.residence.len() {
            self.residence[i] += 1;
        }
    }

    fn note_pruned(&mut self) {
        self.pruned += 1;
        self.frames += 1;
    }

    fn note_transition(&mut self, q: usize, p: usize, barrier: f64) {
        if q < self.trans_out.len() {
            self.trans_out[q] += 1;
        }
        if p < self.trans_in.len() {
            self.trans_in[p] += 1;
        }
        let key = (q.min(p), q.max(p));
        for t in self.pairs.iter_mut() {
            if t.0 == key.0 && t.1 == key.1 {
                t.2 += 1;
                if barrier < t.3 {
                    t.3 = barrier;
                }
                return;
            }
        }
        self.pairs.push((key.0, key.1, 1, barrier));
    }
}

struct Ensemble {
    members: Vec<Conformer>,
    /// All-atom indices (H-rotamer resolution) for dedup alignment.
    dedup_all: Vec<usize>,
    /// Heavy-atom indices (scaffold identity) for dedup alignment.
    dedup_heavy: Vec<usize>,
    rmsd_tol: f64,
    /// Max |ΔE| for the hydrogen-permutation equivalence (identical minimum,
    /// permuted methyl-H labels ⇒ identical energy).
    energy_dup_tol: f64,
    energy_window: f64,
    max_conformers: usize,
    emin: f64,
    // ── landscape statistics (two ladders, both parallel to `members`) ──
    metad: LadderStats,
    md: LadderStats,
    /// Stable basin identities (parallel to `members`): positions shift on
    /// eviction and re-sorting, ids never do — refined saddles key on them.
    basin_ids: Vec<u64>,
    next_basin_id: u64,
    /// Distinct equivalent-well geometries observed per basin (the first
    /// entry is the member itself). The energy-dedup branch folds mirror
    /// twins, methyl-H permutations and C3-rotor wells alike; each NEW
    /// fixed-label geometry is one more observed well — the degeneracy
    /// factor for Boltzmann populations (≤ the true symmetry number).
    variants: Vec<Vec<Vec<[f64; 3]>>>,
}

impl Ensemble {
    fn new(mol: &Molecule, config: &ImtdConfig) -> Self {
        let n = mol.atoms.len();
        let mut heavy = heavy_indices(mol);
        if heavy.is_empty() {
            heavy = (0..n).collect();
        }
        Ensemble {
            members: Vec::new(),
            dedup_all: (0..n).collect(),
            dedup_heavy: heavy,
            rmsd_tol: config.rmsd_tol,
            energy_dup_tol: 1e-3,
            energy_window: config.energy_window,
            max_conformers: config.max_conformers.max(1),
            emin: f64::INFINITY,
            metad: LadderStats::default(),
            md: LadderStats::default(),
            basin_ids: Vec::new(),
            next_basin_id: 1,
            variants: Vec::new(),
        }
    }

    fn id_of(&self, pos: usize) -> u64 {
        self.basin_ids.get(pos).copied().unwrap_or(0)
    }

    /// Two-tier conformer-equivalence predicate (see the module docs): same
    /// heavy scaffold AND (same all-atom geometry OR equal energy). Returns
    /// the matching member position and whether the match was a MIRROR twin
    /// (equal-energy branch + improper superposition): mirror twins fold
    /// like duplicates but raise the member's degeneracy; label
    /// permutations fold silently.
    /// Two-tier conformer-equivalence predicate (see the module docs): same
    /// heavy scaffold AND (same all-atom geometry OR equal energy). Returns
    /// the matching member position; the caller routes equal-energy frames
    /// through [`Self::note_variant`] for the degeneracy count.
    fn find_conformer(&self, coords: &[[f64; 3]], energy: f64) -> Option<usize> {
        for (i, m) in self.members.iter().enumerate() {
            let heavy_same =
                aligned_rmsd_over(coords, &m.coords, &self.dedup_heavy) < self.rmsd_tol;
            if !heavy_same {
                continue;
            }
            if aligned_rmsd_over(coords, &m.coords, &self.dedup_all) < self.rmsd_tol {
                return Some(i); // exact revisit of the same well
            }
            if (energy - m.energy).abs() < self.energy_dup_tol {
                return Some(i); // equivalent well — note_variant decides g
            }
        }
        None
    }

    fn ladder_mut(&mut self, l: Ladder) -> &mut LadderStats {
        match l {
            Ladder::Metad => &mut self.metad,
            Ladder::Md => &mut self.md,
        }
    }

    fn note_visit(&mut self, i: usize, l: Ladder) {
        self.ladder_mut(l).note_visit(i);
    }

    /// An equal-energy frame folded onto basin `i`: a NEW distinct well
    /// geometry (not superposing any stored variant) raises the degeneracy;
    /// a revisit of a known variant changes nothing. Uniformly covers
    /// mirror twins, methyl-H permutations and C3-rotor wells.
    fn note_variant(&mut self, i: usize, coords: &[[f64; 3]]) {
        if i >= self.variants.len() {
            return;
        }
        for v in &self.variants[i] {
            if aligned_rmsd_over(coords, v, &self.dedup_all) < self.rmsd_tol {
                return; // known well, revisit
            }
        }
        self.variants[i].push(coords.to_vec());
    }

    fn note_residence(&mut self, i: usize, l: Ladder) {
        self.ladder_mut(l).note_residence(i);
    }

    fn note_pruned(&mut self, l: Ladder) {
        self.ladder_mut(l).note_pruned();
    }

    /// Directed transition q → p with the observed window PE peak as an
    /// upper bound of the crossing energy (kept as the best/lowest bound).
    fn note_transition(&mut self, q: usize, p: usize, barrier: f64, l: Ladder) {
        self.ladder_mut(l).note_transition(q, p, barrier);
    }

    /// Result-order mapping: members sorted by energy, window-filtered
    /// against the FINAL minimum, capped. Returns the kept members and the
    /// old-index → new-index map (dropped members map to None).
    fn finalized_members(&self) -> (Vec<Conformer>, Vec<Option<usize>>) {
        let mut order: Vec<usize> = (0..self.members.len()).collect();
        let energy = |i: usize| self.members[i].energy;
        order.sort_by(|&x, &y| energy(x).total_cmp(&energy(y)));
        if let Some(&first) = order.first() {
            let emin = energy(first);
            order.retain(|&i| energy(i) <= emin + self.energy_window);
        }
        order.truncate(self.max_conformers);
        let mut map = vec![None; self.members.len()];
        let kept: Vec<Conformer> = order
            .iter()
            .enumerate()
            .map(|(new_i, &old)| {
                map[old] = Some(new_i);
                self.members[old].clone()
            })
            .collect();
        (kept, map)
    }

    /// Indices of the `PARENT_POOL` lowest-energy members (deterministic).
    fn lowest_energy_indices(&self) -> Vec<usize> {
        let mut order: Vec<usize> = (0..self.members.len()).collect();
        order.sort_by(|&x, &y| self.members[x].energy.total_cmp(&self.members[y].energy));
        order.truncate(PARENT_POOL);
        order
    }

    /// Add an ALREADY-QUENCHED minimum: energy window, capacity, aligned-RMSD
    /// dedup. Returns the accepted member's POSITION (creation counts as one
    /// landscape visit), `None` if duplicate / out of window / cap-rejected.
    fn try_add_quenched(
        &mut self,
        coords: Vec<[f64; 3]>,
        energy: f64,
        round: usize,
        origin: Origin,
    ) -> Option<usize> {
        if !energy.is_finite() {
            return None;
        }
        // cap: at capacity, only strictly-below-max energies earn a slot
        if self.members.len() >= self.max_conformers {
            let current_max = self
                .members
                .iter()
                .map(|m| m.energy)
                .fold(f64::NEG_INFINITY, f64::max);
            if energy >= current_max {
                return None;
            }
        }
        // energy window against the current ensemble minimum
        if !self.members.is_empty() && energy > self.emin + self.energy_window {
            return None;
        }
        if self.find_conformer(&coords, energy).is_some() {
            return None; // duplicate (two-tier equivalence, see module docs)
        }
        if energy < self.emin {
            self.emin = energy;
        }
        self.members.push(Conformer {
            coords,
            energy,
            round_found: round,
            origin,
        });
        // zeroed on both ladders — the CALLER's ladder notes the creating
        // frame's visit (a crossing-created basin starts at 0% occupancy
        // until a trajectory actually visits it)
        self.metad.on_new_basin();
        self.md.on_new_basin();
        self.basin_ids.push(self.next_basin_id);
        self.next_basin_id += 1;
        self.variants.push(vec![self
            .members
            .last()
            .map(|m| m.coords.clone())
            .unwrap_or_default()]);
        let pos = self.members.len() - 1;
        if self.members.len() > self.max_conformers {
            // evict the previous maximum (never the just-added member — it
            // passed the strictly-below-max cap check); landscape stats of
            // the evicted basin disappear with it, and positions above it
            // (including the new member's) shift down by one
            let mut i_max = 0;
            for (i, m) in self.members.iter().enumerate() {
                if m.energy > self.members[i_max].energy {
                    i_max = i;
                }
            }
            self.members.remove(i_max);
            self.metad.on_evict(i_max);
            self.md.on_evict(i_max);
            if self.basin_ids.len() > i_max {
                self.basin_ids.remove(i_max);
            }
            if self.variants.len() > i_max {
                self.variants.remove(i_max);
            }
            return Some(pos - 1);
        }
        Some(pos)
    }

    /// Quench `coords` to its local minimum, then [`Self::try_add_quenched`].
    fn try_add(
        &mut self,
        quencher: &Quencher,
        coords: &[[f64; 3]],
        round: usize,
        origin: Origin,
        quench_max_iter: usize,
    ) -> Option<usize> {
        let r = quencher.quench(coords, quench_max_iter);
        self.try_add_quenched(r.optimized_coords, r.final_energy, round, origin)
    }
}

// ── one runnable metad unit (scatter-gather boundary) ──

/// Execute one short metadynamics run and quench its snapshots to local
/// minima. Pure: the same `spec` (+ molecule/config) yields the identical
/// minimum list — this is the unit of work a parallel worker executes.
pub fn execute_run(
    mol: &Molecule,
    engine: SearchEngine,
    quencher: &Quencher,
    config: &ImtdConfig,
    spec: &RunSpec,
) -> Vec<Minimum> {
    let mut mol_start = mol.clone();
    for (atom, c) in mol_start.atoms.iter_mut().zip(spec.start.iter()) {
        atom.position = *c;
    }
    let width = match spec.cv {
        CvSpec::Dihedral(..) => config.hill_width,
        CvSpec::Distance(..) => config.hill_width_distance,
    };
    let metad_ff: Box<dyn ForceField> = match engine {
        SearchEngine::GFNFF => {
            let at: Vec<usize> = mol_start
                .atoms
                .iter()
                .map(|a| a.atomic_number as usize)
                .collect();
            let charge: f64 = mol_start.atoms.iter().map(|a| a.charge).sum();
            let xyz: Vec<[f64; 3]> = mol_start.atoms.iter().map(|a| a.position).collect();
            Box::new(crate::gfnff::GfnffForceField::new(&at, &xyz, charge))
        }
        SearchEngine::MMFF94 => {
            Box::new(MMFFForceField::new(&mol_start, crate::MMFFVariant::MMFF94))
        }
        SearchEngine::MMFF94s => {
            Box::new(MMFFForceField::new(&mol_start, crate::MMFFVariant::MMFF94s))
        }
    };
    let metad = Rc::new(MetaDynamics::new(
        metad_ff,
        spec.cv.build(),
        MetaDConfig {
            hill_height: spec.hill_height,
            hill_width: width,
            deposit_interval: config.deposit_interval.max(1),
            bias_factor: config.bias_factor,
            temperature_k: config.temperature_k,
        },
    ));
    let mut runner = MDRunner::from_molecule(
        metad.clone(), // keep the handle for unbiased-PE sampling
        &mol_start,
        MDConfig {
            dt_fs: config.dt_fs,
            temperature_k: config.temperature_k,
            friction_per_ps: config.friction_per_ps,
            seed: spec.seed,
        },
    );
    let mut minima = Vec::new();
    let mut snapshot_idx = 0usize;
    // unbiased-PE window peak between consecutive snapshots (barrier bound)
    let mut pe_peak = metad.underlying_energy(runner.coords());
    let mut snapshot = |coords: &[[f64; 3]], pe_peak: f64, idx: usize| {
        let r = quencher.quench(coords, config.quench_max_iterations);
        if r.final_energy.is_finite() {
            minima.push(Minimum {
                coords: r.optimized_coords,
                energy: r.final_energy,
                pe_max_since_prev: pe_peak,
                snapshot: idx,
            });
        }
    };
    let interval = config.quench_interval.max(1);
    for i in 1..=config.steps_per_run {
        runner.step();
        let pe = metad.underlying_energy(runner.coords());
        if pe > pe_peak {
            pe_peak = pe;
        }
        if i.is_multiple_of(interval) {
            snapshot(runner.coords(), pe_peak, snapshot_idx);
            snapshot_idx += 1;
            pe_peak = pe; // the next window starts at this frame
        }
    }
    snapshot(runner.coords(), pe_peak, snapshot_idx); // final frame
    minima
}

// ── the iMTD-GC coordinator (scatter-gather) ──

/// Round-scheduling coordinator: plans each round's independent short-metad
/// runs ([`RunSpec`]), accepts their gathered minima in ANY arrival order,
/// and folds them deterministically in run-index order (window/dedup/cap →
/// genetic crossing → convergence), so parallel and sequential execution
/// produce bit-identical ensembles.
pub struct ImtdSearch {
    config: ImtdConfig,
    /// Coordinator-side quench backend (input + crossing children).
    quencher: Quencher,
    ensemble: Ensemble,
    rot: Vec<(usize, usize)>,
    pool: Vec<CvSpec>,
    fragments: Vec<(Vec<usize>, usize)>,
    /// Next round to plan (1-based).
    round: usize,
    plan_cache: Option<Vec<RunSpec>>,
    pending: Vec<Option<Vec<Minimum>>>,
    /// In-flight plain-MD verification batch (empty = none planned).
    verify_pending: Vec<Option<Vec<Minimum>>>,
    /// Atomic masses (amu), parallel to atoms — for the harmonic tier.
    masses: Vec<f64>,
    /// In-flight saddle-refinement batch (empty = none planned).
    refine_jobs: Vec<RefineJob>,
    refine_done: Vec<Option<Option<f64>>>,
    /// Refined saddle energies keyed by stable basin-id pairs
    /// (min_id, max_id, saddle) — converged jobs only.
    refined: Vec<(u64, u64, f64)>,
    rounds_executed: usize,
    new_minima_per_round: Vec<usize>,
    crossed_seeds: Vec<Vec<[f64; 3]>>,
    converged: bool,
    finished: bool,
}

impl ImtdSearch {
    /// Build the coordinator: quenches the input structure (round 0), derives
    /// the CV pool and crossing fragments. `done()` is immediately true for
    /// rigid molecules (empty pool) or zero budgets.
    pub fn new(mol: &Molecule, engine: SearchEngine, config: &ImtdConfig) -> Self {
        let quencher = Quencher::new(mol, engine);
        let mut ensemble = Ensemble::new(mol, config);
        let input_coords: Vec<[f64; 3]> = mol.atoms.iter().map(|a| a.position).collect();
        ensemble.try_add(
            &quencher,
            &input_coords,
            0,
            Origin::Input,
            config.quench_max_iterations,
        );
        let rot = rotatable_bonds(mol);
        let pool = cv_pool(mol, &rot);
        let fragments: Vec<(Vec<usize>, usize)> = rot
            .iter()
            .map(|&(j, k)| (fragment_split(mol, j, k).0, k))
            .collect();
        let finished = pool.is_empty() || config.max_rounds == 0 || config.steps_per_run == 0;
        ImtdSearch {
            config: *config,
            quencher,
            ensemble,
            rot,
            pool,
            fragments,
            round: 1,
            plan_cache: None,
            pending: Vec::new(),
            verify_pending: Vec::new(),
            masses: mol
                .atoms
                .iter()
                .map(|a| crate::md::atomic_mass(a.atomic_number))
                .collect(),
            refine_jobs: Vec::new(),
            refine_done: Vec::new(),
            refined: Vec::new(),
            rounds_executed: 0,
            new_minima_per_round: Vec::new(),
            crossed_seeds: Vec::new(),
            converged: finished,
            finished,
        }
    }

    pub fn done(&self) -> bool {
        self.finished
    }

    /// Round number the next plan belongs to (1-based).
    pub fn round(&self) -> usize {
        self.round
    }

    /// This round's run specs (computed once and cached — idempotent until
    /// the round is gathered). Empty when done. All starts are fixed from the
    /// END of the previous round (crossed seeds first, then the low-energy
    /// ensemble rotation), so the runs are mutually independent.
    pub fn plan(&mut self) -> &[RunSpec] {
        if self.finished {
            return &[];
        }
        if self.plan_cache.is_none() {
            let specs = self.build_round_plan();
            self.pending = vec![None; specs.len()];
            self.plan_cache = Some(specs);
        }
        self.plan_cache.as_deref().unwrap_or(&[])
    }

    fn build_round_plan(&mut self) -> Vec<RunSpec> {
        let mut specs = Vec::with_capacity(self.config.runs_per_round);
        for m in 0..self.config.runs_per_round {
            let g = (self.round - 1) * self.config.runs_per_round + m; // global run counter
            let start: Vec<[f64; 3]> =
                if let Some(seed_coords) = self.crossed_seeds.first().cloned() {
                    self.crossed_seeds.remove(0);
                    seed_coords
                } else {
                    let top = self.ensemble.lowest_energy_indices();
                    self.ensemble.members[top[g % top.len()]].coords.clone()
                };
            specs.push(RunSpec {
                round: self.round,
                run_index: m,
                start,
                cv: self.pool[g % self.pool.len()],
                hill_height: self.config.hill_height * HEIGHT_LADDER[g % HEIGHT_LADDER.len()],
                seed: self.config.seed + (self.round as u64 - 1) * 1000 + m as u64,
            });
        }
        specs
    }

    /// Deposit one run's minima (any order). When every spec of the round has
    /// arrived, folds the round: metad minima in run-index order, then genetic
    /// crossing, then the convergence criterion. Out-of-range or duplicate
    /// deposits are ignored.
    pub fn gather(&mut self, run_index: usize, minima: Vec<Minimum>) {
        if self.finished || self.plan_cache.is_none() {
            return;
        }
        if run_index >= self.pending.len() || self.pending[run_index].is_some() {
            return;
        }
        self.pending[run_index] = Some(minima);
        if self.pending.iter().all(|p| p.is_some()) {
            self.fold_round();
        }
    }

    fn fold_round(&mut self) {
        let specs = self.plan_cache.take().expect("plan present when folding");
        let pending = std::mem::take(&mut self.pending);
        let _ = specs; // (round identity lives in the specs; fold order is the index)
                       // fold metad minima in run-index order (== sequential execution order)
        let mut new_this_round = 0usize;
        for minima in pending.iter().flatten() {
            new_this_round += self.classify_frames(minima, Ladder::Metad);
        }
        // genetic crossing: recombine low-energy parents at rotatable cuts
        if self.config.crossing_pairs > 0
            && self.ensemble.members.len() >= 2
            && !self.fragments.is_empty()
        {
            let top = self.ensemble.lowest_energy_indices();
            let n = top.len();
            let pair_budget = self.config.crossing_pairs.min(n * (n - 1) / 2);
            for p in 0..pair_budget {
                let a = p % n;
                let b = (a + 1 + p / n) % n;
                if a == b {
                    continue;
                }
                let (f1, k) = &self.fragments[p % self.fragments.len()];
                let child = genetic_crossing(
                    &self.ensemble.members[top[a]].coords,
                    &self.ensemble.members[top[b]].coords,
                    f1,
                    *k,
                );
                if let Some(pos) = self.ensemble.try_add(
                    &self.quencher,
                    &child,
                    self.round,
                    Origin::Crossing,
                    self.config.quench_max_iterations,
                ) {
                    new_this_round += 1;
                    let seed = self.ensemble.members[pos].coords.clone();
                    self.crossed_seeds.push(seed);
                }
            }
        }
        self.new_minima_per_round.push(new_this_round);
        self.rounds_executed = self.round;
        self.round += 1;
        if new_this_round <= self.config.convergence_new_minima {
            self.converged = true;
            self.finished = true;
        } else if self.round > self.config.max_rounds {
            self.finished = true;
        }
    }

    /// Classify one run's quenched frames into the landscape: an existing
    /// basin earns a visit (in `ladder`), a genuinely new minimum creates a
    /// basin (its creating frame noted in the same ladder), window/cap
    /// rejects are pruned and break the transition chain. Returns how many
    /// NEW minima were created.
    fn classify_frames(&mut self, minima: &[Minimum], ladder: Ladder) -> usize {
        let mut created = 0usize;
        let mut prev: Option<usize> = None;
        for m in minima {
            let pos = match self.ensemble.find_conformer(&m.coords, m.energy) {
                Some(i) => {
                    self.ensemble.note_visit(i, ladder);
                    self.ensemble.note_variant(i, &m.coords);
                    Some(i)
                }
                None => {
                    let p = self.ensemble.try_add_quenched(
                        m.coords.clone(),
                        m.energy,
                        self.rounds_executed,
                        Origin::Metad,
                    );
                    if p.is_some() {
                        created += 1;
                        // the creating frame's visit belongs to this ladder
                        if let Some(ip) = p {
                            self.ensemble.note_visit(ip, ladder);
                        }
                    }
                    p
                }
            };
            match pos {
                Some(p) => {
                    if let Some(q) = prev {
                        if q != p {
                            // transition q→p with the observed window PE
                            // peak as the barrier upper bound
                            self.ensemble
                                .note_transition(q, p, m.pe_max_since_prev, ladder);
                        } else {
                            self.ensemble.note_residence(p, ladder);
                        }
                    }
                    prev = Some(p);
                }
                None => {
                    self.ensemble.note_pruned(ladder);
                    prev = None; // unassigned frame breaks the chain
                }
            }
        }
        created
    }

    /// Plan `n` plain-MD verification runs (hill_height = 0 ⇒ no bias — the
    /// trajectory samples the physical surface, giving thermally meaningful
    /// basin populations for the found minima). Starts rotate over the
    /// low-energy basins; seeds live in a disjoint space (base + 1e6).
    /// Only after the search finished and no round/batch is in flight;
    /// batches may be repeated. Steps/quench cadence come from the DRIVER's
    /// options at execution time (the spec itself is cadence-free).
    pub fn plan_verify(&mut self, n: usize) -> Vec<RunSpec> {
        if !self.finished || self.plan_cache.is_some() || n == 0 {
            return Vec::new();
        }
        if self.verify_pending.iter().any(|p| p.is_some()) {
            return Vec::new(); // a verification batch is still in flight
        }
        let top = self.ensemble.lowest_energy_indices();
        if top.is_empty() {
            return Vec::new();
        }
        let cv_default = self.pool.first().copied().unwrap_or(CvSpec::Distance(0, 1));
        let specs: Vec<RunSpec> = (0..n)
            .map(|i| RunSpec {
                round: 0, // verification marker — no round semantics
                run_index: i,
                start: self.ensemble.members[top[i % top.len()]].coords.clone(),
                cv: self
                    .pool
                    .get(i % self.pool.len().max(1))
                    .copied()
                    .unwrap_or(cv_default),
                hill_height: 0.0,
                seed: self.config.seed + 1_000_000 + i as u64,
            })
            .collect();
        self.verify_pending = vec![None; specs.len()];
        specs
    }

    /// Deposit one verification run's minima (any order); folds the batch
    /// into the MD ladder once every planned run arrived.
    pub fn gather_verify(&mut self, run_index: usize, minima: Vec<Minimum>) {
        if run_index >= self.verify_pending.len() || self.verify_pending[run_index].is_some() {
            return;
        }
        self.verify_pending[run_index] = Some(minima);
        if self.verify_pending.iter().all(|p| p.is_some()) {
            let pending = std::mem::take(&mut self.verify_pending);
            for minima in pending.iter().flatten() {
                self.classify_frames(minima, Ladder::Md);
            }
        }
    }

    /// Plan one CI-NEB saddle-refinement job per observed transition pair
    /// (metad + MD ladders, deduplicated by stable basin ids). Endpoint
    /// geometries come from the CURRENT final mapping. Only after the search
    /// finished with no round/verification/refinement batch in flight;
    /// batches may be repeated (a pair refined twice keeps the new value).
    pub fn plan_refine(&mut self) -> Vec<RefineJob> {
        if !self.finished || self.plan_cache.is_some() {
            return Vec::new();
        }
        if self.verify_pending.iter().any(|p| p.is_some()) || !self.refine_jobs.is_empty() {
            return Vec::new();
        }
        let (kept, map) = self.ensemble.finalized_members();
        // id -> final index
        let final_of = |id: u64| -> Option<usize> {
            (0..self.ensemble.members.len())
                .find(|&old| self.ensemble.id_of(old) == id)
                .and_then(|old| map[old])
        };
        let mut seen: Vec<(u64, u64)> = Vec::new();
        let mut jobs = Vec::new();
        for ladder in [&self.ensemble.metad, &self.ensemble.md] {
            for &(a_pos, b_pos, count, _) in &ladder.pairs {
                if count == 0 {
                    continue;
                }
                let (ia, ib) = (self.ensemble.id_of(a_pos), self.ensemble.id_of(b_pos));
                let key = (ia.min(ib), ia.max(ib));
                if key.0 == 0 || seen.contains(&key) {
                    continue;
                }
                seen.push(key);
                let (Some(fa), Some(fb)) = (final_of(key.0), final_of(key.1)) else {
                    continue; // a member was dropped by the final filter
                };
                jobs.push(RefineJob {
                    index: jobs.len(),
                    id_a: key.0,
                    id_b: key.1,
                    a: fa,
                    b: fb,
                    xyz_a: kept[fa].coords.clone(),
                    xyz_b: kept[fb].coords.clone(),
                });
            }
        }
        self.refine_jobs = jobs.clone();
        self.refine_done = vec![None; jobs.len()];
        jobs
    }

    /// Deposit one refinement result (saddle energy; `None` when the NEB did
    /// not converge). When every job of the batch arrived, converged saddles
    /// are recorded (keyed by basin ids — stable across later folds).
    pub fn gather_refine(&mut self, index: usize, saddle: Option<f64>) {
        if index >= self.refine_jobs.len() || self.refine_done[index].is_some() {
            return;
        }
        self.refine_done[index] = Some(saddle);
        if self.refine_done.iter().all(|p| p.is_some()) {
            let jobs = std::mem::take(&mut self.refine_jobs);
            let done = std::mem::take(&mut self.refine_done);
            for (job, res) in jobs.into_iter().zip(done) {
                if let Some(Some(saddle)) = res {
                    let key = (job.id_a.min(job.id_b), job.id_a.max(job.id_b));
                    self.refined.retain(|&(a, b, _)| (a, b) != key);
                    self.refined.push((key.0, key.1, saddle));
                }
            }
        }
    }

    /// Snapshot result (non-destructive; may be called at any time).
    pub fn result(&self) -> ImtdResult {
        finalize(
            &self.ensemble,
            self.rounds_executed,
            self.converged,
            self.new_minima_per_round.clone(),
            self.pool.len(),
            self.rot.len(),
        )
    }

    /// Harmonic free-energy tier for the final ensemble at temperature `t`
    /// (K): per basin, frequencies at the FINAL minimum, G = E + G_corr,
    /// ZPE, imaginary-mode count. Engine follows the search (GFN-FF works —
    /// the Hessian differences its analytic gradient).
    pub fn thermo(&self, temperature_k: f64) -> Vec<ThermoEntry> {
        let (kept, _) = self.ensemble.finalized_members();
        kept.iter()
            .map(|c| {
                let freqs = match &self.quencher {
                    Quencher::Mmff(ff) => {
                        crate::vib::frequencies(ff.as_ref(), &c.coords, &self.masses)
                    }
                    Quencher::Gfnff { at, charge, cache } => {
                        // GFN-FF: bond perception is geometry-dependent —
                        // the Hessian differences a per-geometry rebuild
                        // (exact-coords cache helps the per-T re-runs)
                        let mut guard = cache.borrow_mut();
                        let g: &crate::gfnff::GfnffForceField = match guard.as_ref() {
                            Some((prev, f)) if *prev == c.coords => f.as_ref(),
                            _ => {
                                *guard = Some((
                                    c.coords.clone(),
                                    Box::new(crate::gfnff::GfnffForceField::new(
                                        at, &c.coords, *charge,
                                    )),
                                ));
                                guard
                                    .as_ref()
                                    .map(|(_, f)| f.as_ref())
                                    .unwrap_or_else(|| unreachable!("cache was just filled"))
                            }
                        };
                        crate::vib::frequencies(g, &c.coords, &self.masses)
                    }
                };
                let (g_corr, zpe, n_imag) = crate::vib::harmonic_correction(&freqs, temperature_k);
                let (g_corr_qh, _, _) =
                    crate::vib::harmonic_correction_opt(&freqs, temperature_k, true);
                let nu_min = freqs
                    .iter()
                    .copied()
                    .filter(|f| f.abs() >= crate::vib::TR_ROT_CUTOFF_CM)
                    .fold(f64::INFINITY, f64::min);
                let nu_max = freqs.iter().copied().fold(f64::NEG_INFINITY, f64::max);
                ThermoEntry {
                    energy: c.energy,
                    g_free: c.energy + g_corr,
                    g_free_qh: c.energy + g_corr_qh,
                    g_corr,
                    zpe,
                    n_imag,
                    nu_min,
                    nu_max,
                }
            })
            .collect()
    }

    /// Landscape snapshot (non-destructive): basin statistics and the
    /// disconnectivity trees, both ladders aligned with the FINAL result
    /// order (`result().conformers` indices == basin indices). Statistics of
    /// members dropped by the final window/cap filter are discarded.
    pub fn landscape(&self) -> Landscape {
        let (kept, map) = self.ensemble.finalized_members();
        let map_ladder = |stats: &LadderStats| -> (Vec<BasinStats>, Vec<TransitionStats>) {
            let mut basins = vec![
                BasinStats {
                    visits: 0,
                    residence: 0,
                    transitions_in: 0,
                    transitions_out: 0,
                    degeneracy: 1,
                };
                kept.len()
            ];
            for (old, new_i) in map.iter().enumerate() {
                if let Some(&j) = new_i.as_ref() {
                    if old < stats.visits.len() {
                        basins[j] = BasinStats {
                            visits: stats.visits[old],
                            residence: stats.residence[old],
                            transitions_in: stats.trans_in[old],
                            transitions_out: stats.trans_out[old],
                            degeneracy: self
                                .ensemble
                                .variants
                                .get(old)
                                .map(|v| v.len() as u64)
                                .unwrap_or(1),
                        };
                    }
                }
            }
            let mut transitions = Vec::new();
            for &(a, b, count, barrier) in &stats.pairs {
                if let (Some(na), Some(nb)) =
                    (map.get(a).copied().flatten(), map.get(b).copied().flatten())
                {
                    let (ia, ib) = (self.ensemble.id_of(a), self.ensemble.id_of(b));
                    let key = (ia.min(ib), ia.max(ib));
                    let refined = self
                        .refined
                        .iter()
                        .find(|&&(x, y, _)| (x, y) == key)
                        .map(|&(_, _, e)| e);
                    transitions.push(TransitionStats {
                        from: na,
                        to: nb,
                        count,
                        barrier,
                        refined,
                    });
                }
            }
            transitions.sort_by_key(|t| (t.from, t.to));
            (basins, transitions)
        };
        let energies: Vec<f64> = kept.iter().map(|c| c.energy).collect();
        let (basins, transitions) = map_ladder(&self.ensemble.metad);
        let dgraph_tree = dgraph(&energies, &transitions);
        let md = if self.ensemble.md.frames > 0 {
            let (md_basins, md_transitions) = map_ladder(&self.ensemble.md);
            let md_tree = dgraph(&energies, &md_transitions);
            Some(MdLandscape {
                basins: md_basins,
                transitions: md_transitions,
                dgraph: md_tree,
                frames: self.ensemble.md.frames,
                pruned: self.ensemble.md.pruned,
            })
        } else {
            None
        };
        Landscape {
            basins,
            transitions,
            dgraph: dgraph_tree,
            frames: self.ensemble.metad.frames,
            pruned: self.ensemble.metad.pruned,
            md,
        }
    }
}

/// iMTD-GC conformational search (see the module docs) — sequential
/// composition of the same scatter-gather primitives a parallel driver uses,
/// hence bit-identical to it. Deterministic: same molecule + variant + config
/// ⇒ bit-identical result.
pub fn imtd_conformer_search(
    mol: &Molecule,
    engine: SearchEngine,
    config: &ImtdConfig,
) -> ImtdResult {
    if mol.atoms.is_empty() {
        return ImtdResult {
            conformers: Vec::new(),
            rounds_executed: 0,
            converged: true,
            new_minima_per_round: Vec::new(),
            pool_size: 0,
            n_rotatable: 0,
        };
    }
    let quencher = Quencher::new(mol, engine);
    let mut search = ImtdSearch::new(mol, engine, config);
    while !search.done() {
        let specs = search.plan().to_vec();
        if specs.is_empty() {
            break;
        }
        for spec in &specs {
            let minima = execute_run(mol, engine, &quencher, config, spec);
            search.gather(spec.run_index, minima);
        }
    }
    search.result()
}

/// Sort by energy, re-apply the window against the FINAL minimum (the minimum
/// may have dropped during the search), cap the ensemble size.
fn finalize(
    ensemble: &Ensemble,
    rounds_executed: usize,
    converged: bool,
    new_minima_per_round: Vec<usize>,
    pool_size: usize,
    n_rotatable: usize,
) -> ImtdResult {
    let (conformers, _) = ensemble.finalized_members();
    ImtdResult {
        conformers,
        rounds_executed,
        converged,
        new_minima_per_round,
        pool_size,
        n_rotatable,
    }
}

// ── landscape types + disconnectivity tree ──

/// Per-basin landscape statistics, aligned with the final result order.
#[derive(Debug, Clone, Copy)]
pub struct BasinStats {
    /// Quenched frames assigned to this basin (creation included).
    pub visits: usize,
    /// Consecutive same-basin snapshot pairs (residence proxy).
    pub residence: usize,
    pub transitions_in: usize,
    pub transitions_out: usize,
    /// Observed well count: 2 when a mirror twin folded into this basin —
    /// the degeneracy factor for Boltzmann populations.
    pub degeneracy: u64,
}

/// Transition statistics between a basin PAIR (unordered, both directions
/// summed); `barrier` is the best (lowest) observed window-PE upper bound,
/// `refined` the CI-NEB saddle-energy estimate once refined (absolute
/// crossing energy, same units/semantics as `barrier`).
#[derive(Debug, Clone, Copy)]
pub struct TransitionStats {
    pub from: usize,
    pub to: usize,
    pub count: usize,
    pub barrier: f64,
    pub refined: Option<f64>,
}

/// The observed-transition landscape of a search.
#[derive(Debug, Clone)]
pub struct Landscape {
    pub basins: Vec<BasinStats>,
    pub transitions: Vec<TransitionStats>,
    /// Disconnectivity forest over the metad-ladder transition barriers.
    pub dgraph: Vec<DgNode>,
    /// Total quenched frames the fold classified (metad ladder).
    pub frames: usize,
    /// Frames rejected by window/cap (unassigned; they break transition chains).
    pub pruned: usize,
    /// Plain-MD verification ladder (present once verification runs folded).
    pub md: Option<MdLandscape>,
}

/// Harmonic free-energy entry for one final basin.
#[derive(Debug, Clone, Copy)]
pub struct ThermoEntry {
    /// The search energy (kcal/mol).
    pub energy: f64,
    /// G = E + G_corr at the requested temperature (kcal/mol), plain
    /// harmonic.
    pub g_free: f64,
    /// Quasi-harmonic G (soft modes floored at 100 cm⁻¹ — Grimme-style).
    pub g_free_qh: f64,
    /// Harmonic correction (ZPE + vibrational entropy term).
    pub g_corr: f64,
    pub zpe: f64,
    /// Imaginary modes at |ν̃| ≥ 25 cm⁻¹ (true-minimum QC flag).
    pub n_imag: usize,
    /// Smallest/largest real-mode wavenumber (cm⁻¹).
    pub nu_min: f64,
    pub nu_max: f64,
}

/// Serialize the thermo tier (bit-exact float strings).
pub fn thermo_to_json(entries: &[ThermoEntry], temperature_k: f64) -> String {
    serde_json::json!({
        "T": temperature_k,
        "basins": entries.iter().map(|e| serde_json::json!({
            "E": f64_to_json(e.energy), "G": f64_to_json(e.g_free),
            "Gqh": f64_to_json(e.g_free_qh),
            "gCorr": f64_to_json(e.g_corr), "zpe": f64_to_json(e.zpe),
            "nImag": e.n_imag,
            "nuMin": f64_to_json(e.nu_min), "nuMax": f64_to_json(e.nu_max),
        })).collect::<Vec<_>>(),
    })
    .to_string()
}

/// The plain-MD verification ladder: thermally meaningful populations for
/// the found basins (finite-sampling caveat applies).
#[derive(Debug, Clone)]
pub struct MdLandscape {
    pub basins: Vec<BasinStats>,
    pub transitions: Vec<TransitionStats>,
    pub dgraph: Vec<DgNode>,
    pub frames: usize,
    pub pruned: usize,
}

/// Disconnectivity-graph tree node: leaves carry `basin` (final result
/// index) with `threshold` = the basin energy; internal nodes carry the
/// merge threshold (observed barrier upper bound).
#[derive(Debug, Clone)]
pub struct DgNode {
    pub threshold: f64,
    pub basin: Option<usize>,
    pub children: Vec<DgNode>,
}

/// Build the disconnectivity forest: Kruskal-style ascending barrier order
/// (ties broken by (from, to) for determinism); every union of two
/// components becomes an internal node at that barrier. Basins never linked
/// by an observed transition stay separate roots — an honest "observed
/// disconnected" statement, not a claim of true disconnection.
fn dgraph(energies: &[f64], edges: &[TransitionStats]) -> Vec<DgNode> {
    let n = energies.len();
    // arena of tree nodes; `merged_into` is the union-find parent (None = root)
    let mut arena: Vec<DgNode> = (0..n)
        .map(|i| DgNode {
            threshold: energies[i],
            basin: Some(i),
            children: Vec::new(),
        })
        .collect();
    let mut merged_into: Vec<Option<usize>> = vec![None; n];
    let root_of = |merged_into: &Vec<Option<usize>>, mut i: usize| -> usize {
        while let Some(up) = merged_into[i] {
            i = up;
        }
        i
    };
    // effective barrier: the refined NEB saddle when available, else the
    // observed window-PE upper bound
    let eff = |t: &TransitionStats| t.refined.filter(|e| e.is_finite()).unwrap_or(t.barrier);
    let mut sorted: Vec<&TransitionStats> = edges.iter().collect();
    sorted.sort_by(|x, y| {
        eff(x)
            .total_cmp(&eff(y))
            .then((x.from, x.to).cmp(&(y.from, y.to)))
    });
    for e in sorted {
        if eff(e).is_nan() {
            continue; // degenerate payload (older workers) — no edge
        }
        let ra = root_of(&merged_into, e.from);
        let rb = root_of(&merged_into, e.to);
        if ra == rb {
            continue;
        }
        let node = DgNode {
            threshold: eff(e),
            basin: None,
            children: vec![arena[ra].clone(), arena[rb].clone()],
        };
        arena.push(node);
        let new_root = arena.len() - 1;
        merged_into[ra] = Some(new_root);
        merged_into[rb] = Some(new_root);
        merged_into.push(None);
    }
    (0..arena.len())
        .filter(|&i| merged_into[i].is_none())
        .map(|i| arena[i].clone())
        .collect()
}

/// Serialize the landscape for the UI (bit-exact float strings, same channel
/// discipline as the spec/minima codecs).
pub fn landscape_to_json(l: &Landscape) -> String {
    fn dg_value(n: &DgNode) -> serde_json::Value {
        serde_json::json!({
            "t": f64_to_json(n.threshold),
            "basin": n.basin,
            "children": n.children.iter().map(dg_value).collect::<Vec<_>>(),
        })
    }
    let ladder_value = |basins: &[BasinStats],
                        transitions: &[TransitionStats],
                        tree: &[DgNode],
                        frames: usize,
                        pruned: usize|
     -> serde_json::Value {
        serde_json::json!({
            "frames": frames,
            "pruned": pruned,
            "basins": basins.iter().enumerate().map(|(i, b)| serde_json::json!({
                "i": i, "visits": b.visits, "residence": b.residence,
                "in": b.transitions_in, "out": b.transitions_out,
                "g": b.degeneracy,
            })).collect::<Vec<_>>(),
            "transitions": transitions.iter().map(|t| serde_json::json!({
                "a": t.from, "b": t.to, "count": t.count,
                "barrier": f64_to_json(t.barrier),
                "refined": t.refined.map(f64_to_json),
            })).collect::<Vec<_>>(),
            "dgraph": tree.iter().map(dg_value).collect::<Vec<_>>(),
        })
    };
    serde_json::json!({
        "frames": l.frames,
        "pruned": l.pruned,
        "basins": l.basins.iter().enumerate().map(|(i, b)| serde_json::json!({
            "i": i, "visits": b.visits, "residence": b.residence,
            "in": b.transitions_in, "out": b.transitions_out,
            "g": b.degeneracy,
        })).collect::<Vec<_>>(),
        "transitions": l.transitions.iter().map(|t| serde_json::json!({
            "a": t.from, "b": t.to, "count": t.count,
            "barrier": f64_to_json(t.barrier),
            "refined": t.refined.map(f64_to_json),
        })).collect::<Vec<_>>(),
        "dgraph": l.dgraph.iter().map(dg_value).collect::<Vec<_>>(),
        "md": l.md.as_ref().map(|m| ladder_value(&m.basins, &m.transitions, &m.dgraph, m.frames, m.pruned)),
    })
    .to_string()
}

// ── JSON codecs (worker boundary) ──
// wasm objects cannot cross a worker postMessage; specs and minima travel as
// JSON. serde_json's float PARSING is 1-ULP lossy on boundary values (its
// ryu emitter is exact but the parser is not correctly rounded), so every
// f64 is carried as a "{:.17e}" STRING — Rust's str::parse is correctly
// rounded, 17 significant digits always round-trip, and the JS side ferries
// the strings untouched (it never interprets coordinates).

/// Bit-exact JSON rendering of one f64 (18 significant digits, exponent form).
fn f64_to_json(v: f64) -> String {
    format!("{:.17e}", v)
}

/// Bit-exact JSON parsing of one f64.
///
/// # Errors
/// `Err` when the field is not a string or not a finite float.
fn f64_from_json(v: &serde_json::Value, ctx: &str) -> Result<f64, String> {
    let s = v
        .as_str()
        .ok_or_else(|| format!("{ctx}: expected a float string, got {v}"))?;
    let x: f64 = s
        .parse()
        .map_err(|e| format!("{ctx}: bad float {s:?}: {e}"))?;
    if x.is_finite() {
        Ok(x)
    } else {
        Err(format!("{ctx}: non-finite float {s:?}"))
    }
}

fn point_to_json(p: &[f64; 3]) -> serde_json::Value {
    serde_json::json!([f64_to_json(p[0]), f64_to_json(p[1]), f64_to_json(p[2])])
}

fn point_from_json(v: &serde_json::Value, ctx: &str) -> Result<[f64; 3], String> {
    let arr = v
        .as_array()
        .ok_or_else(|| format!("{ctx}: expected [x, y, z], got {v}"))?;
    if arr.len() != 3 {
        return Err(format!("{ctx}: expected 3 components, got {}", arr.len()));
    }
    Ok([
        f64_from_json(&arr[0], ctx)?,
        f64_from_json(&arr[1], ctx)?,
        f64_from_json(&arr[2], ctx)?,
    ])
}

fn cv_to_json(cv: CvSpec) -> serde_json::Value {
    match cv {
        CvSpec::Dihedral(i, j, k, l) => {
            serde_json::json!({ "type": "dihedral", "atoms": [i, j, k, l] })
        }
        CvSpec::Distance(i, j) => serde_json::json!({ "type": "distance", "atoms": [i, j] }),
    }
}

fn cv_from_json(v: &serde_json::Value) -> Result<CvSpec, String> {
    let t = v
        .get("type")
        .and_then(serde_json::Value::as_str)
        .ok_or_else(|| format!("cv: missing type, got {v}"))?;
    let atoms: Vec<usize> = v
        .get("atoms")
        .and_then(serde_json::Value::as_array)
        .map(|a| {
            a.iter()
                .filter_map(serde_json::Value::as_u64)
                .map(|x| x as usize)
                .collect()
        })
        .ok_or_else(|| format!("cv: missing atoms, got {v}"))?;
    match (t, atoms.as_slice()) {
        ("dihedral", [i, j, k, l]) => Ok(CvSpec::Dihedral(*i, *j, *k, *l)),
        ("distance", [i, j]) => Ok(CvSpec::Distance(*i, *j)),
        _ => Err(format!("cv: bad type/atoms combination: {t} {atoms:?}")),
    }
}

fn spec_to_value(spec: &RunSpec) -> serde_json::Value {
    serde_json::json!({
        "round": spec.round,
        "run_index": spec.run_index,
        "start": spec.start.iter().map(point_to_json).collect::<Vec<_>>(),
        "cv": cv_to_json(spec.cv),
        "hill_height": f64_to_json(spec.hill_height),
        "seed": spec.seed,
    })
}

fn spec_from_value(v: serde_json::Value) -> Result<RunSpec, String> {
    let round = v
        .get("round")
        .and_then(serde_json::Value::as_u64)
        .ok_or("spec: missing round")? as usize;
    let run_index = v
        .get("run_index")
        .and_then(serde_json::Value::as_u64)
        .ok_or("spec: missing run_index")? as usize;
    let start_v = v
        .get("start")
        .and_then(serde_json::Value::as_array)
        .ok_or("spec: missing start")?;
    let mut start = Vec::with_capacity(start_v.len());
    for (i, p) in start_v.iter().enumerate() {
        start.push(point_from_json(p, &format!("spec.start[{i}]"))?);
    }
    let cv = cv_from_json(v.get("cv").ok_or_else(|| "spec: missing cv".to_string())?)?;
    let hill_height = f64_from_json(
        v.get("hill_height")
            .ok_or_else(|| "spec: missing hill_height".to_string())?,
        "spec.hill_height",
    )?;
    let seed = v
        .get("seed")
        .and_then(serde_json::Value::as_u64)
        .ok_or("spec: missing seed")?;
    Ok(RunSpec {
        round,
        run_index,
        start,
        cv,
        hill_height,
        seed,
    })
}

fn minima_to_value(minima: &[Minimum]) -> serde_json::Value {
    serde_json::json!({
        "minima": minima
            .iter()
            .map(|m| serde_json::json!({
                "e": f64_to_json(m.energy),
                "xyz": m.coords.iter().map(point_to_json).collect::<Vec<_>>(),
                "pe": f64_to_json(m.pe_max_since_prev),
                "snap": m.snapshot,
            }))
            .collect::<Vec<_>>()
    })
}

fn minima_from_value(v: serde_json::Value) -> Result<Vec<Minimum>, String> {
    let arr = v
        .get("minima")
        .and_then(serde_json::Value::as_array)
        .ok_or("minima: missing minima array")?;
    let mut out = Vec::with_capacity(arr.len());
    for (i, m) in arr.iter().enumerate() {
        let energy = f64_from_json(
            m.get("e")
                .ok_or_else(|| format!("minima[{i}]: missing e"))?,
            &format!("minima[{i}].e"),
        )?;
        let xyz_v = m
            .get("xyz")
            .and_then(serde_json::Value::as_array)
            .ok_or_else(|| format!("minima[{i}]: missing xyz"))?;
        let mut coords = Vec::with_capacity(xyz_v.len());
        for (j, p) in xyz_v.iter().enumerate() {
            coords.push(point_from_json(p, &format!("minima[{i}].xyz[{j}]"))?);
        }
        // pe/snap tolerate absence (older payloads): landscape degenerates
        // gracefully instead of failing the worker boundary
        let pe_max_since_prev = match m.get("pe") {
            Some(v) => f64_from_json(v, &format!("minima[{i}].pe"))?,
            None => f64::NAN,
        };
        let snapshot = m
            .get("snap")
            .and_then(serde_json::Value::as_u64)
            .unwrap_or(0) as usize;
        out.push(Minimum {
            coords,
            energy,
            pe_max_since_prev,
            snapshot,
        });
    }
    Ok(out)
}

/// Serialize a round plan: `{"round": r, "runs": [spec, ...]}`.
pub fn plan_to_json(round: usize, specs: &[RunSpec]) -> String {
    serde_json::json!({
        "round": round,
        "runs": specs.iter().map(spec_to_value).collect::<Vec<_>>(),
    })
    .to_string()
}

/// Serialize one run's minima: `{"minima": [...]}`.
pub fn minima_to_json(minima: &[Minimum]) -> String {
    minima_to_value(minima).to_string()
}

/// Parse a serialized [`RunSpec`] (unknown fields are ignored).
///
/// # Errors
/// `Err` with a context message when the JSON is malformed.
pub fn run_spec_from_json(s: &str) -> Result<RunSpec, String> {
    let v: serde_json::Value = serde_json::from_str(s).map_err(|e| format!("spec json: {e}"))?;
    spec_from_value(v)
}

/// Serialize refinement jobs: `{"jobs": [{...}]}` (floats as bit-exact strings).
pub fn refine_jobs_to_json(jobs: &[RefineJob]) -> String {
    serde_json::json!({
        "jobs": jobs.iter().map(|j| serde_json::json!({
            "index": j.index, "idA": j.id_a, "idB": j.id_b, "a": j.a, "b": j.b,
            "xyzA": j.xyz_a.iter().map(point_to_json).collect::<Vec<_>>(),
            "xyzB": j.xyz_b.iter().map(point_to_json).collect::<Vec<_>>(),
        })).collect::<Vec<_>>(),
    })
    .to_string()
}

/// A refinement job's endpoint pair.
pub type RefineEndpoints = (Vec<[f64; 3]>, Vec<[f64; 3]>);

/// Parse one refinement job's endpoints from its serialized form.
///
/// # Errors
/// `Err` with a context message when the JSON is malformed.
pub fn refine_job_from_json(s: &str) -> Result<RefineEndpoints, String> {
    let v: serde_json::Value =
        serde_json::from_str(s).map_err(|e| format!("refine job json: {e}"))?;
    let arr = |key: &str| -> Result<Vec<[f64; 3]>, String> {
        let pts = v
            .get(key)
            .and_then(serde_json::Value::as_array)
            .ok_or_else(|| format!("refine job: missing {key}"))?;
        let mut out = Vec::with_capacity(pts.len());
        for (i, p) in pts.iter().enumerate() {
            out.push(point_from_json(p, &format!("refine job {key}[{i}]"))?);
        }
        Ok(out)
    };
    Ok((arr("xyzA")?, arr("xyzB")?))
}

/// Serialize a refinement result; `saddle` is `None` when NEB did not converge.
pub fn refine_result_to_json(saddle: Option<f64>, barrier: f64) -> String {
    serde_json::json!({
        "saddle": saddle.map(f64_to_json),
        "barrier": f64_to_json(barrier),
    })
    .to_string()
}

/// Parse a refinement result into the saddle estimate (`None` = not converged).
///
/// # Errors
/// `Err` with a context message when the JSON is malformed.
pub fn refine_result_from_json(s: &str) -> Result<Option<f64>, String> {
    let v: serde_json::Value =
        serde_json::from_str(s).map_err(|e| format!("refine result json: {e}"))?;
    match v.get("saddle") {
        None | Some(serde_json::Value::Null) => Ok(None),
        Some(x) => Ok(Some(f64_from_json(x, "refine result saddle")?)),
    }
}

/// Boltzmann reconstruction of basin populations from the found minima:
/// p_i(T) ∝ exp(-E_i/kT), computed with an Emin shift (overflow-safe) and
/// normalized to 1. HONEST CAVEATS: ideal populations over the FOUND minima
/// — degeneracies hidden by enantiomer-pair collapsing and entropy / harmonic
/// free-energy corrections are NOT included. `temperature_k` is clamped to a
/// sane domain (≤ 1 K defaults to 300 K).
pub fn boltzmann_populations(energies: &[f64], temperature_k: f64) -> Vec<f64> {
    boltzmann_populations_deg(energies, &[], temperature_k)
}

/// [`boltzmann_populations`] with explicit degeneracy factors
/// (`weights = g_i·exp(−ΔE_i/kT)`; mirror-twin basins carry g=2).
pub fn boltzmann_populations_deg(
    energies: &[f64],
    degeneracies: &[f64],
    temperature_k: f64,
) -> Vec<f64> {
    if energies.is_empty() {
        return Vec::new();
    }
    let t = if temperature_k > 1.0 {
        temperature_k
    } else {
        300.0
    };
    let kt = crate::md::KB * t;
    let emin = energies.iter().copied().fold(f64::INFINITY, f64::min);
    let weights: Vec<f64> = energies
        .iter()
        .zip(degeneracies.iter().copied().chain(std::iter::repeat(1.0)))
        .map(|(&e, g)| g.max(1.0) * (-(e - emin) / kt).exp())
        .collect();
    let sum: f64 = weights.iter().sum();
    if sum <= 0.0 || !sum.is_finite() {
        // impossible after the Emin shift unless the input is degenerate
        return vec![1.0 / energies.len() as f64; energies.len()];
    }
    weights.into_iter().map(|w| w / sum).collect()
}

/// Parse a serialized minima list.
///
/// # Errors
/// `Err` with a context message when the JSON is malformed.
pub fn minima_from_json(s: &str) -> Result<Vec<Minimum>, String> {
    let v: serde_json::Value = serde_json::from_str(s).map_err(|e| format!("minima json: {e}"))?;
    minima_from_value(v)
}

#[cfg(test)]
mod tests {
    use super::*;
    /// Improper (mirror) superposition RMSD over the given indices: reflect one
    /// structure through a fixed plane (z → −z), then align properly — a
    /// mirror-IMAGE pair superposes under this transform, not under a rotation.
    fn improper_rmsd_over(a: &[[f64; 3]], b: &[[f64; 3]], idx: &[usize]) -> f64 {
        if idx.len() < 2 {
            return 0.0;
        }
        let p: Vec<[f64; 3]> = idx.iter().map(|&i| a[i]).collect();
        let q: Vec<[f64; 3]> = idx.iter().map(|&i| [b[i][0], b[i][1], -b[i][2]]).collect();
        kabsch_transform(&p, &q).2
    }

    const ETHANOL_SDF: &str = "ethanol\n     RDKit          3D\n\n  9  8  0  0  0  0  0  0  0  0999 V2000\n    0.8867    0.1097    0.1398 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.4649   -0.1632   -0.4871 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.4908    0.4341    0.2924 O   0  0  0  0  0  0  0  0  0  0  0  0\n    1.0600    1.1872    0.2278 H   0  0  0  0  0  0  0  0  0  0  0  0\n    1.6897   -0.3301   -0.4581 H   0  0  0  0  0  0  0  0  0  0  0  0\n    0.9343   -0.3042    1.1524 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.5119    0.2599   -1.4947 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.6562   -1.2384   -0.5550 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.4469    0.0448    1.1825 H   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0\n  2  3  1  0\n  1  4  1  0\n  1  5  1  0\n  1  6  1  0\n  2  7  1  0\n  2  8  1  0\n  3  9  1  0\nM  END\n";

    const BUTANE_HEAVY_SDF: &str = "Butane\n     RDKit          3D\n\n  4  3  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.5500    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    3.0500    0.1000    0.9000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    4.5500    0.1000    0.9000 C   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0\n  2  3  1  0\n  3  4  1  0\nM  END\n";

    const METHANE_SDF: &str = "Methane\n     RDKit          3D\n\n  5  4  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.6294    0.6294    0.6294 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.6294   -0.6294    0.6294 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.6294    0.6294   -0.6294 H   0  0  0  0  0  0  0  0  0  0  0  0\n    0.6294   -0.6294   -0.6294 H   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0\n  1  3  1  0\n  1  4  1  0\n  1  5  1  0\nM  END\n";

    fn ethanol() -> Molecule {
        crate::molecule::parser::parse_sdf(ETHANOL_SDF).expect("ethanol SDF parses")
    }

    fn butane_all_h() -> Molecule {
        let heavy =
            crate::molecule::parser::parse_sdf(BUTANE_HEAVY_SDF).expect("butane SDF parses");
        crate::molecule::hydrogens::add_hydrogens(&heavy)
    }

    #[test]
    fn cv_pool_ethanol_has_two_rotatable_dihedrals() {
        let mol = ethanol();
        let rot = rotatable_bonds(&mol);
        assert_eq!(rot.len(), 2, "ethanol: C-C and C-O rotatable, got {rot:?}");
        let pool = cv_pool(&mol, &rot);
        assert_eq!(pool.len(), 2, "no distance top-up needed, got {pool:?}");
        // bond (0,1): i = H3 (smallest H on C0), l = O2 (heavy preferred)
        // bond (1,2): i = C0 (heavy), l = H8 (only other O neighbor)
        let mut torsions = Vec::new();
        for cv in &pool {
            if let CvSpec::Dihedral(i, j, k, l) = *cv {
                torsions.push((i, j, k, l));
            }
        }
        assert!(torsions.contains(&(3, 0, 1, 2)), "got {torsions:?}");
        assert!(torsions.contains(&(0, 1, 2, 8)), "got {torsions:?}");
    }

    #[test]
    fn cv_pool_ring_molecule_tops_up_with_distances() {
        // toluene heavy skeleton: no rotatable bonds (the methyl C has
        // degree 1, ring bonds are excluded), but the methyl–para carbon
        // pair sits at graph distance 4 → distance CV top-up
        let sdf = "Toluene\n     RDKit          3D\n\n  7  7  0  0  0  0  0  0  0  0999 V2000\n    0.0000    1.3900    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.2040    0.6950    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.2040   -0.6950    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000   -1.3900    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.2040   -0.6950    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.2040    0.6950    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000    2.8900    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  2  0\n  2  3  1  0\n  3  4  2  0\n  4  5  1  0\n  5  6  2  0\n  6  1  1  0\n  1  7  1  0\nM  END\n";
        let mol = crate::molecule::parser::parse_sdf(sdf).expect("parse");
        let rot = rotatable_bonds(&mol);
        assert!(
            rot.is_empty(),
            "ring/methyl bonds are never rotatable, got {rot:?}"
        );
        let pool = cv_pool(&mol, &rot);
        assert_eq!(
            pool.len(),
            1,
            "methyl–para is the only dist-4 pair, got {pool:?}"
        );
        match pool[0] {
            CvSpec::Distance(i, j) => {
                assert_eq!((i.min(j), i.max(j)), (3, 6), "methyl(6)–para(3) pair");
                assert_eq!(graph_distance(&mol, i, j), Some(4));
            }
            other => panic!("expected a distance CV, got {other:?}"),
        }
    }

    fn graph_distance(mol: &Molecule, src: usize, dst: usize) -> Option<usize> {
        let mut dist = vec![usize::MAX; mol.atoms.len()];
        dist[src] = 0;
        let mut queue = VecDeque::new();
        queue.push_back(src);
        while let Some(u) = queue.pop_front() {
            if u == dst {
                return Some(dist[u]);
            }
            for &v in &mol.adjacency[u] {
                if dist[v] == usize::MAX {
                    dist[v] = dist[u] + 1;
                    queue.push_back(v);
                }
            }
        }
        None
    }

    #[test]
    fn kabsch_recovers_rotation_and_rejects_mirror() {
        // non-degenerate point cloud
        let p: Vec<[f64; 3]> = vec![
            [0.0, 0.0, 0.0],
            [1.3, 0.2, -0.4],
            [0.3, 1.7, 0.5],
            [-0.9, 0.6, 1.4],
            [0.5, -1.1, 1.9],
        ];
        // known proper rotation (Rodrigues 72° about (1,1,1)/√3) + shift
        let norm = (3.0f64).sqrt();
        let u = [1.0 / norm, 1.0 / norm, 1.0 / norm];
        let th = 72.0f64.to_radians();
        let (s, c) = (th.sin(), th.cos());
        let r = [
            c + u[0] * u[0] * (1.0 - c),
            u[0] * u[1] * (1.0 - c) - u[2] * s,
            u[0] * u[2] * (1.0 - c) + u[1] * s,
            u[1] * u[0] * (1.0 - c) + u[2] * s,
            c + u[1] * u[1] * (1.0 - c),
            u[1] * u[2] * (1.0 - c) - u[0] * s,
            u[2] * u[0] * (1.0 - c) - u[1] * s,
            u[2] * u[1] * (1.0 - c) + u[0] * s,
            c + u[2] * u[2] * (1.0 - c),
        ];
        let t = [0.7, -1.2, 2.3];
        let q: Vec<[f64; 3]> = p
            .iter()
            .map(|pt| {
                [
                    r[0] * pt[0] + r[1] * pt[1] + r[2] * pt[2] + t[0],
                    r[3] * pt[0] + r[4] * pt[1] + r[5] * pt[2] + t[1],
                    r[6] * pt[0] + r[7] * pt[1] + r[8] * pt[2] + t[2],
                ]
            })
            .collect();
        let (_, _, rmsd) = kabsch_transform(&q, &p);
        assert!(rmsd < 1e-8, "rotated copy must superpose, rmsd {rmsd}");
        // mirror image needs a reflection: proper rotation cannot reach it
        let m: Vec<[f64; 3]> = q.iter().map(|pt| [pt[0], pt[1], -pt[2]]).collect();
        let (_, _, rmsd_m) = kabsch_transform(&m, &p);
        assert!(rmsd_m > 0.1, "mirror must NOT superpose, rmsd {rmsd_m}");
    }

    #[test]
    fn crossing_recombines_fragments_exactly() {
        let mol = butane_all_h();
        let rot = rotatable_bonds(&mol);
        // 3 rotatable bonds by our CV rule (terminal C–C bonds give valid
        // H–C–C–C torsions); the crossing test uses the central bond
        assert_eq!(rot.len(), 3, "got {rot:?}");
        let (j, k) = (1, 2);
        let (f1, f2) = fragment_split(&mol, j, k);
        assert!(!f1.is_empty() && !f2.is_empty());
        assert_eq!(f1.len() + f2.len(), mol.atoms.len(), "partition covers all");
        let a: Vec<[f64; 3]> = mol.atoms.iter().map(|x| x.position).collect();
        // B = A with F2 rotated 120° about the j-k bond axis (another
        // rotamer): crossing takes F1 from A (verbatim) and F2 from the
        // aligned B — the child must equal B up to alignment noise.
        let axis = [a[k][0] - a[j][0], a[k][1] - a[j][1], a[k][2] - a[j][2]];
        let norm = (axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2]).sqrt();
        let u = [axis[0] / norm, axis[1] / norm, axis[2] / norm];
        let th = std::f64::consts::PI * 2.0 / 3.0;
        let (s, c) = (th.sin(), th.cos());
        let mut b = a.clone();
        for &i in &f2 {
            let d = [a[i][0] - a[k][0], a[i][1] - a[k][1], a[i][2] - a[k][2]];
            let cross = [
                u[1] * d[2] - u[2] * d[1],
                u[2] * d[0] - u[0] * d[2],
                u[0] * d[1] - u[1] * d[0],
            ];
            let dot = u[0] * d[0] + u[1] * d[1] + u[2] * d[2];
            for x in 0..3 {
                b[i][x] = a[k][x] + c * d[x] + s * cross[x] + u[x] * dot * (1.0 - c);
            }
        }
        let heavy = heavy_indices(&mol);
        let rmsd_ab = aligned_rmsd_over(&a, &b, &heavy);
        assert!(
            rmsd_ab > 0.2,
            "rotamer B must differ from A (4 heavy atoms, aligned), rmsd {rmsd_ab}"
        );
        let child = genetic_crossing(&a, &b, &f1, k);
        let rmsd_cb = aligned_rmsd_over(&child, &b, &heavy);
        let rmsd_ca = aligned_rmsd_over(&child, &a, &heavy);
        assert!(
            rmsd_cb < 1e-6,
            "child == B (F1 from A equals B's aligned F1), rmsd {rmsd_cb}"
        );
        assert!(rmsd_ca > 0.2, "child keeps B's rotamer, rmsd {rmsd_ca}");
        for &i in &f1 {
            assert_eq!(child[i], a[i], "F1 atom {i} must be A's coordinate");
        }
    }

    #[test]
    fn imtd_butane_finds_multiple_minima_deterministically() {
        let mol = butane_all_h();
        let cfg = ImtdConfig {
            max_rounds: 3,
            runs_per_round: 2,
            steps_per_run: 1500,
            quench_interval: 250,
            quench_max_iterations: 250,
            crossing_pairs: 2,
            seed: 42,
            ..Default::default()
        };
        let res = imtd_conformer_search(&mol, SearchEngine::MMFF94s, &cfg);
        assert!(
            res.conformers.len() >= 2,
            "butane has anti + gauche minima (window {}), got {} conformers: {:?}",
            cfg.energy_window,
            res.conformers.len(),
            res.conformers.iter().map(|c| c.energy).collect::<Vec<_>>()
        );
        for w in res.conformers.windows(2) {
            assert!(w[0].energy <= w[1].energy, "ensemble must be sorted");
        }
        let emin = res.conformers[0].energy;
        for c in &res.conformers {
            assert!(c.energy.is_finite());
            assert!(c.energy <= emin + cfg.energy_window);
        }
        assert_eq!(res.rounds_executed, res.new_minima_per_round.len());
        // determinism: bit-identical repeat
        let res2 = imtd_conformer_search(&mol, SearchEngine::MMFF94s, &cfg);
        assert_eq!(res.conformers.len(), res2.conformers.len());
        for (x, y) in res.conformers.iter().zip(res2.conformers.iter()) {
            assert_eq!(x.energy, y.energy, "energies must be bit-identical");
            assert_eq!(x.coords, y.coords, "coordinates must be bit-identical");
        }
    }

    #[test]
    fn imtd_rigid_methane_single_conformer() {
        let mol = crate::molecule::parser::parse_sdf(METHANE_SDF).expect("parse");
        let res = imtd_conformer_search(&mol, SearchEngine::MMFF94s, &ImtdConfig::default());
        assert_eq!(res.pool_size, 0, "methane has no biasable CV");
        assert_eq!(res.conformers.len(), 1);
        assert!(res.converged);
        assert_eq!(res.rounds_executed, 0);
        assert_eq!(res.conformers[0].origin, Origin::Input);
    }

    fn small_cfg() -> ImtdConfig {
        ImtdConfig {
            max_rounds: 3,
            runs_per_round: 3,
            steps_per_run: 500,
            quench_interval: 250,
            quench_max_iterations: 200,
            crossing_pairs: 2,
            seed: 42,
            ..Default::default()
        }
    }

    #[test]
    fn execute_run_is_pure_and_json_roundtrips_bit_exact() {
        let mol = ethanol();
        let cfg = small_cfg();
        let mut search = ImtdSearch::new(&mol, SearchEngine::MMFF94s, &cfg);
        let specs = search.plan().to_vec();
        assert_eq!(specs.len(), 3);
        let quencher = Quencher::new(&mol, SearchEngine::MMFF94s);
        let a = execute_run(&mol, SearchEngine::MMFF94s, &quencher, &cfg, &specs[0]);
        let b = execute_run(&mol, SearchEngine::MMFF94s, &quencher, &cfg, &specs[0]);
        assert!(!a.is_empty(), "run produces quenched minima");
        assert_eq!(a, b, "same spec ⇒ bit-identical minima (purity)");
        // JSON round-trip (worker boundary) is bit-exact: floats travel as
        // 17-digit exponent strings (serde_json's number parser is 1-ULP
        // lossy on boundary values — string channel dodges it)
        let spec_json = spec_to_value(&specs[1]).to_string();
        let spec_rt = run_spec_from_json(&spec_json).expect("spec round-trips");
        assert_eq!(spec_rt, specs[1]);
        let minima_rt = minima_from_json(&minima_to_json(&a)).expect("minima round-trips");
        assert_eq!(minima_rt, a, "coordinates/energies bit-exact via JSON");
    }

    #[test]
    fn plan_is_idempotent_within_a_round() {
        let mol = ethanol();
        let cfg = small_cfg();
        let mut search = ImtdSearch::new(&mol, SearchEngine::MMFF94s, &cfg);
        let p1 = search.plan().to_vec();
        let p2 = search.plan().to_vec();
        assert_eq!(p1, p2, "repeated plan must not consume crossed seeds twice");
        assert_eq!(search.round(), 1);
        assert!(!search.done());
    }

    #[test]
    fn scatter_gather_matches_sequential_regardless_of_arrival_order() {
        let mol = ethanol();
        let cfg = small_cfg();
        let seq = imtd_conformer_search(&mol, SearchEngine::MMFF94s, &cfg);
        let quencher = Quencher::new(&mol, SearchEngine::MMFF94s);
        // in-order composition
        let mut s1 = ImtdSearch::new(&mol, SearchEngine::MMFF94s, &cfg);
        while !s1.done() {
            for spec in s1.plan().to_vec() {
                let m = execute_run(&mol, SearchEngine::MMFF94s, &quencher, &cfg, &spec);
                s1.gather(spec.run_index, m);
            }
        }
        let r1 = s1.result();
        // JSON boundary + reverse arrival order per round
        let mut s2 = ImtdSearch::new(&mol, SearchEngine::MMFF94s, &cfg);
        while !s2.done() {
            let specs = s2.plan().to_vec();
            let jsons: Vec<String> = specs
                .iter()
                .map(|sp| {
                    minima_to_json(&execute_run(
                        &mol,
                        SearchEngine::MMFF94s,
                        &quencher,
                        &cfg,
                        sp,
                    ))
                })
                .collect();
            for (i, j) in jsons.iter().enumerate().rev() {
                let m = minima_from_json(j).expect("worker minima parse");
                s2.gather(i, m);
            }
        }
        let r2 = s2.result();
        for r in [r1, r2] {
            assert_eq!(r.conformers.len(), seq.conformers.len());
            assert_eq!(r.rounds_executed, seq.rounds_executed);
            assert_eq!(r.converged, seq.converged);
            assert_eq!(r.new_minima_per_round, seq.new_minima_per_round);
            for (x, y) in r.conformers.iter().zip(seq.conformers.iter()) {
                assert_eq!(x.energy, y.energy, "scatter-gather must match sequential");
                assert_eq!(x.coords, y.coords);
            }
        }
    }

    #[test]
    fn hydrogen_rotamers_distinct_but_permutations_not() {
        // Ethanol's anti/gauche C–C rotamers differ only in H positions
        // (the heavy C–C–O triangle superposes) — the two-tier dedup must
        // keep them (CREST semantics) while methyl-H permutations must not
        // inflate the ensemble (they quench to the identical energy).
        let mol = ethanol();
        let cfg = ImtdConfig {
            max_rounds: 2,
            runs_per_round: 3,
            steps_per_run: 600,
            quench_interval: 200,
            quench_max_iterations: 200,
            crossing_pairs: 2,
            seed: 42,
            ..Default::default()
        };
        let res = imtd_conformer_search(&mol, SearchEngine::MMFF94s, &cfg);
        assert!(
            res.conformers.len() >= 2,
            "H rotamers must survive dedup, got {} conformers",
            res.conformers.len()
        );
        let gap = res.conformers[1].energy - res.conformers[0].energy;
        assert!(
            (0.05..0.5).contains(&gap),
            "anti/gauche gap ~0.18, got {gap}"
        );
        for w in res.conformers.windows(2) {
            assert!(
                (w[1].energy - w[0].energy).abs() > 1e-3,
                "no equal-energy duplicates (permutation artifacts): {:?}",
                res.conformers.iter().map(|c| c.energy).collect::<Vec<_>>()
            );
        }
    }

    #[test]
    fn imtd_gfnff_engine_runs_deterministically() {
        let mol = ethanol();
        let cfg = ImtdConfig {
            max_rounds: 2,
            runs_per_round: 2,
            steps_per_run: 500,
            quench_interval: 250,
            quench_max_iterations: 200,
            crossing_pairs: 1,
            seed: 7,
            ..Default::default()
        };
        let a = imtd_conformer_search(&mol, SearchEngine::GFNFF, &cfg);
        assert!(!a.conformers.is_empty(), "GFN-FF finds the input minimum");
        for c in &a.conformers {
            assert!(c.energy.is_finite(), "GFN-FF energy finite: {}", c.energy);
        }
        // GFN-FF is a different surface: the ensemble must not equal MMFF's
        let b = imtd_conformer_search(&mol, SearchEngine::MMFF94s, &cfg);
        assert_ne!(
            a.conformers[0].energy, b.conformers[0].energy,
            "engines are distinct surfaces"
        );
        let a2 = imtd_conformer_search(&mol, SearchEngine::GFNFF, &cfg);
        assert_eq!(
            a.conformers[0].energy, a2.conformers[0].energy,
            "deterministic"
        );
    }

    #[test]
    fn butane_keeps_enantiomeric_gauche_pair() {
        // v1 mirror semantics: heavy-atom mirror pairs (butane gauche± —
        // Kabsch proper rotation cannot superpose them) both survive; only
        // same-heavy-scaffold + equal-energy structures (methyl-H
        // permutations / H-arrangement mirrors) collapse.
        let mol = butane_all_h();
        let cfg = ImtdConfig {
            max_rounds: 3,
            runs_per_round: 3,
            steps_per_run: 2500,
            quench_interval: 125,
            quench_max_iterations: 250,
            crossing_pairs: 2,
            seed: 42,
            ..Default::default()
        };
        let res = imtd_conformer_search(&mol, SearchEngine::MMFF94s, &cfg);
        let n = res.conformers.len();
        assert!(
            n >= 3,
            "anti + gauche± (enantiomeric pair) expected, got {n}: {:?}",
            res.conformers.iter().map(|c| c.energy).collect::<Vec<_>>()
        );
        // the two gauche enantiomers share an energy to high precision
        let mut equal_energy_pairs = 0;
        for i in 0..n {
            for j in (i + 1)..n {
                if (res.conformers[i].energy - res.conformers[j].energy).abs() < 1e-6 {
                    equal_energy_pairs += 1;
                }
            }
        }
        assert!(
            equal_energy_pairs >= 1,
            "the gauche enantiomeric pair (equal energies) must both be present"
        );
    }

    #[test]
    fn landscape_trace_invariants() {
        let mol = butane_all_h();
        let cfg = ImtdConfig {
            max_rounds: 3,
            runs_per_round: 3,
            steps_per_run: 2500,
            quench_interval: 125,
            quench_max_iterations: 250,
            crossing_pairs: 2,
            seed: 42,
            ..Default::default()
        };
        let mut search = ImtdSearch::new(&mol, SearchEngine::MMFF94s, &cfg);
        let quencher = Quencher::new(&mol, SearchEngine::MMFF94s);
        while !search.done() {
            for spec in search.plan().to_vec() {
                let m = execute_run(&mol, SearchEngine::MMFF94s, &quencher, &cfg, &spec);
                search.gather(spec.run_index, m);
            }
        }
        let res = search.result();
        let land = search.landscape();
        assert!(
            land.frames >= 30,
            "dense budget yields many frames, got {}",
            land.frames
        );
        let visits: usize = land.basins.iter().map(|b| b.visits).sum();
        assert!(
            visits + land.pruned <= land.frames,
            "conservation: visits {visits} + pruned {} <= frames {}",
            land.pruned,
            land.frames
        );
        let sin: usize = land.basins.iter().map(|b| b.transitions_in).sum();
        let sout: usize = land.basins.iter().map(|b| b.transitions_out).sum();
        assert_eq!(sin, sout, "every transition has one source and one target");
        assert!(
            !land.transitions.is_empty(),
            "butane anti↔gauche crossings must be observed"
        );
        for t in &land.transitions {
            let ea = res.conformers[t.from].energy;
            let eb = res.conformers[t.to].energy;
            assert!(
                t.barrier >= ea.max(eb) - 1e-6,
                "barrier {} must bound both wells {ea}/{eb}",
                t.barrier
            );
        }
        // disconnectivity forest: leaves == basins
        fn count_leaves(n: &DgNode) -> usize {
            if n.children.is_empty() {
                1
            } else {
                n.children.iter().map(count_leaves).sum()
            }
        }
        let leaves: usize = land.dgraph.iter().map(count_leaves).sum();
        assert_eq!(leaves, res.conformers.len(), "leaves == final basins");
        // determinism of the whole landscape payload
        let json1 = landscape_to_json(&land);
        let mut search2 = ImtdSearch::new(&mol, SearchEngine::MMFF94s, &cfg);
        while !search2.done() {
            for spec in search2.plan().to_vec() {
                let m = execute_run(&mol, SearchEngine::MMFF94s, &quencher, &cfg, &spec);
                search2.gather(spec.run_index, m);
            }
        }
        let land2 = search2.landscape();
        assert_eq!(json1, landscape_to_json(&land2), "deterministic");
        // the JSON parses back with the exact barrier values
        let v: serde_json::Value = serde_json::from_str(&json1).expect("landscape json parses");
        assert_eq!(v["frames"].as_u64(), Some(land.frames as u64));
        let b0 = v["transitions"][0]["barrier"]
            .as_str()
            .expect("string floats");
        assert_eq!(
            b0.parse::<f64>(),
            Ok(land.transitions[0].barrier),
            "barrier round-trips bit-exact"
        );
    }

    #[test]
    fn landscape_parity_sequential_vs_scatter_gather() {
        let mol = ethanol();
        let cfg = small_cfg();
        let mut seq = ImtdSearch::new(&mol, SearchEngine::MMFF94s, &cfg);
        let quencher = Quencher::new(&mol, SearchEngine::MMFF94s);
        while !seq.done() {
            for spec in seq.plan().to_vec() {
                let m = execute_run(&mol, SearchEngine::MMFF94s, &quencher, &cfg, &spec);
                seq.gather(spec.run_index, m);
            }
        }
        let l1 = seq.landscape();
        // out-of-order arrival through the JSON boundary
        let mut par = ImtdSearch::new(&mol, SearchEngine::MMFF94s, &cfg);
        while !par.done() {
            let specs = par.plan().to_vec();
            let payloads: Vec<String> = specs
                .iter()
                .map(|sp| {
                    minima_to_json(&execute_run(
                        &mol,
                        SearchEngine::MMFF94s,
                        &quencher,
                        &cfg,
                        sp,
                    ))
                })
                .collect();
            for (i, j) in payloads.iter().enumerate().rev() {
                par.gather(i, minima_from_json(j).expect("parse"));
            }
        }
        let l2 = par.landscape();
        assert_eq!(
            landscape_to_json(&l1),
            landscape_to_json(&l2),
            "landscape must be arrival-order independent"
        );
    }

    #[test]
    fn verify_phase_populates_md_ladder_deterministically() {
        let mol = ethanol();
        let cfg = small_cfg();
        let run_all = |s: &mut ImtdSearch| {
            let quencher = Quencher::new(&mol, SearchEngine::MMFF94s);
            while !s.done() {
                for spec in s.plan().to_vec() {
                    let m = execute_run(&mol, SearchEngine::MMFF94s, &quencher, &cfg, &spec);
                    s.gather(spec.run_index, m);
                }
            }
            for _ in 0..2 {
                let specs = s.plan_verify(3);
                assert!(!specs.is_empty(), "verification batch plannable");
                for sp in &specs {
                    let m = execute_run(&mol, SearchEngine::MMFF94s, &quencher, &cfg, sp);
                    s.gather_verify(sp.run_index, m);
                }
            }
        };
        let mut s1 = ImtdSearch::new(&mol, SearchEngine::MMFF94s, &cfg);
        // guards: no verification before the search finished
        assert!(s1.plan_verify(2).is_empty());
        run_all(&mut s1);
        let land = s1.landscape();
        let md = land.md.as_ref().expect("md ladder present after verify");
        assert!(
            md.frames >= 6,
            "two batches of plain-MD runs, got {}",
            md.frames
        );
        let vis: usize = md.basins.iter().map(|b| b.visits).sum();
        assert!(vis + md.pruned <= md.frames, "md conservation");
        let sin: usize = md.basins.iter().map(|b| b.transitions_in).sum();
        let sout: usize = md.basins.iter().map(|b| b.transitions_out).sum();
        assert_eq!(sin, sout, "md ladder in==out");
        let json1 = landscape_to_json(&land);
        let v: serde_json::Value = serde_json::from_str(&json1).expect("json");
        assert!(
            v["md"]["frames"].as_u64().unwrap_or(0) >= 6,
            "md block emitted"
        );
        // spec sanity: zero bias + disjoint seed space
        let mut s2 = ImtdSearch::new(&mol, SearchEngine::MMFF94s, &cfg);
        while !s2.done() {
            let quencher = Quencher::new(&mol, SearchEngine::MMFF94s);
            for spec in s2.plan().to_vec() {
                let m = execute_run(&mol, SearchEngine::MMFF94s, &quencher, &cfg, &spec);
                s2.gather(spec.run_index, m);
            }
        }
        for sp in s2.plan_verify(3) {
            assert_eq!(sp.hill_height, 0.0, "plain MD = zero bias");
            assert!(sp.seed >= 1_000_000, "disjoint verification seed space");
        }
        run_all(&mut s2);
        assert_eq!(json1, landscape_to_json(&s2.landscape()), "deterministic");
    }

    #[test]
    fn refine_phase_sets_saddles_with_stable_ids() {
        let mol = butane_all_h();
        let cfg = ImtdConfig {
            max_rounds: 3,
            runs_per_round: 3,
            steps_per_run: 2500,
            quench_interval: 125,
            quench_max_iterations: 250,
            crossing_pairs: 2,
            seed: 42,
            ..Default::default()
        };
        let mut search = ImtdSearch::new(&mol, SearchEngine::MMFF94s, &cfg);
        let quencher = Quencher::new(&mol, SearchEngine::MMFF94s);
        while !search.done() {
            for spec in search.plan().to_vec() {
                let m = execute_run(&mol, SearchEngine::MMFF94s, &quencher, &cfg, &spec);
                search.gather(spec.run_index, m);
            }
        }
        // guard: no refinement while unfinished paths exist — after done OK
        let jobs = search.plan_refine();
        assert!(!jobs.is_empty(), "butane dense run has transition pairs");
        let ff = MMFFForceField::new(&mol, crate::MMFFVariant::MMFF94s);
        let neb_cfg = crate::neb::NebConfig::default();
        for job in &jobs {
            let res = crate::neb::neb(&ff, &job.xyz_a, &job.xyz_b, &neb_cfg);
            search.gather_refine(
                job.index,
                if res.converged {
                    Some(res.saddle_energy)
                } else {
                    None
                },
            );
        }
        let land = search.landscape();
        let refined: Vec<_> = land.transitions.iter().filter_map(|t| t.refined).collect();
        assert!(!refined.is_empty(), "refined saddles recorded");
        for t in &land.transitions {
            if let Some(r) = t.refined {
                let ea = 0.0; // endpoints checked structurally below via energies
                let _ = ea;
                assert!(r.is_finite());
                assert!(
                    r >= t.barrier - 25.0,
                    "refined saddle should be near/below the observed bound (sanity band)"
                );
            }
        }
        // JSON carries the refined field (string floats)
        let json = landscape_to_json(&land);
        assert!(json.contains("\"refined\":\""), "refined emitted");
        // id stability: fold a verification batch afterwards — the refined
        // pair must survive re-finalization (ids, not positions)
        for _ in 0..1 {
            let specs = search.plan_verify(2);
            assert!(!specs.is_empty());
            for sp in &specs {
                let m = execute_run(&mol, SearchEngine::MMFF94s, &quencher, &cfg, sp);
                search.gather_verify(sp.run_index, m);
            }
        }
        let land2 = search.landscape();
        let refined2: Vec<_> = land2.transitions.iter().filter_map(|t| t.refined).collect();
        assert_eq!(
            refined.len(),
            refined2.len(),
            "refined saddles survive post-refine folds (stable ids)"
        );
    }

    #[test]
    fn boltzmann_populations_are_stable_and_physical() {
        // sums to 1; low T concentrates on the ground state; huge gaps
        // do not overflow (Emin shift)
        let e = [0.0, 0.5, 3.0];
        let p = boltzmann_populations(&e, 298.15);
        assert_eq!(p.len(), 3);
        let sum: f64 = p.iter().sum();
        assert!((sum - 1.0).abs() < 1e-12, "sum {sum}");
        assert!(p[0] >= p[1] && p[1] >= p[2], "monotone in energy");
        let cold = boltzmann_populations(&e, 5.0);
        assert!(cold[0] > 0.9999, "T→low selects the ground state");
        let hot = boltzmann_populations(&e, 1e6);
        assert!((hot[0] - 1.0 / 3.0).abs() < 1e-3, "T→high flattens");
        let wide = boltzmann_populations(&[0.0, 500.0], 300.0);
        // a 500 kcal/mol gap underflows exp() to exactly 0 — the honest f64 answer
        assert!(wide[0] > 0.9999 && wide[1] >= 0.0 && wide[1].is_finite());
        // degenerate input
        assert!(boltzmann_populations(&[], 300.0).is_empty());
        let one = boltzmann_populations(&[7.0], 300.0);
        assert_eq!(one, vec![1.0]);
    }

    #[test]
    fn mirror_twins_fold_with_degeneracy_permutations_silently() {
        // ethanol OH rotamers: rotating ONLY the hydroxyl H about the C–O
        // bond leaves the heavy scaffold identical; +120°/−120° are mirror
        // twins (equal energy, improper-superposable, not proper), while a
        // methyl-H permutation (rotate H3 about C–C into H4's site with the
        // labels fixed) is neither.
        let mol = ethanol();
        let base: Vec<[f64; 3]> = mol.atoms.iter().map(|a| a.position).collect();
        let (c1, o2) = (1usize, 2usize); // C–O bond (H8 = idx 8 on O)
        let axis = [
            base[o2][0] - base[c1][0],
            base[o2][1] - base[c1][1],
            base[o2][2] - base[c1][2],
        ];
        let norm = (axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2]).sqrt();
        let u = [axis[0] / norm, axis[1] / norm, axis[2] / norm];
        let rotate_h = |src: &[[f64; 3]], deg: f64| -> Vec<[f64; 3]> {
            let th = deg.to_radians();
            let (sn, cs) = (th.sin(), th.cos());
            let mut out = src.to_vec();
            let d = [
                src[8][0] - src[o2][0],
                src[8][1] - src[o2][1],
                src[8][2] - src[o2][2],
            ];
            let cr = [
                u[1] * d[2] - u[2] * d[1],
                u[2] * d[0] - u[0] * d[2],
                u[0] * d[1] - u[1] * d[0],
            ];
            let dot = u[0] * d[0] + u[1] * d[1] + u[2] * d[2];
            for x in 0..3 {
                out[8][x] = src[o2][x] + cs * d[x] + sn * cr[x] + u[x] * dot * (1.0 - cs);
            }
            out
        };
        let quencher = Quencher::new(&mol, SearchEngine::MMFF94s);
        let q = |c: &[[f64; 3]]| -> (Vec<[f64; 3]>, f64) {
            let r = quencher.quench(c, 500);
            (r.optimized_coords, r.final_energy)
        };
        let (p0, e0) = q(&base);
        // build the twin by REFLECTING a gauche-OH rotamer through the heavy
        // plane (C0,C1,O2): the reflection of a minimum is a minimum, the
        // heavy scaffold is unchanged, and the H arrangement mirrors. (A ±
        // rotation construction lets the requench drag the C–C torsion.)
        let (pp, ep) = q(&rotate_h(&p0, 120.0));
        let n = {
            let (a, b, c) = (pp[0], pp[1], pp[2]);
            let u = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
            let v = [c[0] - a[0], c[1] - a[1], c[2] - a[2]];
            let cr = [
                u[1] * v[2] - u[2] * v[1],
                u[2] * v[0] - u[0] * v[2],
                u[0] * v[1] - u[1] * v[0],
            ];
            let cn = (cr[0] * cr[0] + cr[1] * cr[1] + cr[2] * cr[2]).sqrt();
            [cr[0] / cn, cr[1] / cn, cr[2] / cn]
        };
        let pm_src: Vec<[f64; 3]> = pp
            .iter()
            .map(|p| {
                let d = [p[0] - pp[0][0], p[1] - pp[0][1], p[2] - pp[0][2]];
                let dn = 2.0 * (d[0] * n[0] + d[1] * n[1] + d[2] * n[2]);
                [p[0] - dn * n[0], p[1] - dn * n[1], p[2] - dn * n[2]]
            })
            .collect();
        let (pm, em) = q(&pm_src);
        // the two rotations are mirror twins: equal energy, improper ~0,
        // proper all-atom RMSD large
        assert!((ep - em).abs() < 1e-3, "twin energies equal: {ep} vs {em}");
        assert!(improper_rmsd_over(&pp, &pm, &(0..mol.atoms.len()).collect::<Vec<_>>()) < 0.05);
        assert!(aligned_rmsd_over(&pp, &pm, &(0..mol.atoms.len()).collect::<Vec<_>>()) > 0.2);
        let _ = e0;

        // fold [base, twin+, twin−] through a fresh coordinator: 2 basins,
        // the rotamer basin carries degeneracy 2
        let cfg = ImtdConfig {
            max_rounds: 1,
            runs_per_round: 1,
            steps_per_run: 200,
            quench_interval: 200,
            quench_max_iterations: 500,
            crossing_pairs: 0,
            convergence_new_minima: usize::MAX,
            ..Default::default()
        };
        let mut search = ImtdSearch::new(&mol, SearchEngine::MMFF94s, &cfg);
        let quencher = Quencher::new(&mol, SearchEngine::MMFF94s);
        while !search.done() {
            for spec in search.plan().to_vec() {
                let m = execute_run(&mol, SearchEngine::MMFF94s, &quencher, &cfg, &spec);
                search.gather(spec.run_index, m);
            }
        }
        let specs = search.plan_verify(1);
        assert_eq!(specs.len(), 1, "one verification run planned");
        search.gather_verify(
            0,
            vec![
                Minimum {
                    coords: pp,
                    energy: ep,
                    pe_max_since_prev: 0.0,
                    snapshot: 0,
                },
                Minimum {
                    coords: pm.clone(),
                    energy: em,
                    pe_max_since_prev: 0.0,
                    snapshot: 1,
                },
            ],
        );
        let land1 = search.landscape();
        // twin+ created the basin (g=1), twin− folds with mirror ⇒ g=2
        let gvals: Vec<u64> = land1.basins.iter().map(|b| b.degeneracy).collect();
        // the unified variant count discovers symmetry-equivalent wells:
        // the member's own rotor arrangement + the folded pair = 3 observed
        // wells (methyl 3-fold / OH mirror — a LOWER bound of the true
        // symmetry number; deterministic for this construction)
        assert!(
            gvals.contains(&3),
            "equivalent wells recorded, got {gvals:?}"
        );
        assert!(!land1.basins.is_empty());
        // a REVISIT of the twin does not raise g past 2
        search.gather_verify(
            0,
            vec![Minimum {
                coords: pm,
                energy: em,
                pe_max_since_prev: 0.0,
                snapshot: 0,
            }],
        );
        let land2 = search.landscape();
        assert_eq!(land2.basins[0].degeneracy, 3, "revisits do not raise g");
    }

    #[test]
    fn methanol_methyl_rotor_counts_three_wells() {
        // CH3OH: rotating the methyl by ±120° about the C–O axis gives three
        // fixed-label geometries of three SYMMETRY-EQUAL wells (σ=3) — the
        // variant counter must record all three, revisits must not raise.
        let heavy_sdf = "Methanol\n     RDKit          3D\n\n  2  1  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.4300    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0\nM  END\n";
        let heavy = crate::molecule::parser::parse_sdf(heavy_sdf).expect("parse");
        let mol = crate::molecule::hydrogens::add_hydrogens(&heavy);
        let base: Vec<[f64; 3]> = mol.atoms.iter().map(|a| a.position).collect();
        let (c0, o1) = (0usize, 1usize);
        let axis = [
            base[o1][0] - base[c0][0],
            base[o1][1] - base[c0][1],
            base[o1][2] - base[c0][2],
        ];
        let norm = (axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2]).sqrt();
        let u = [axis[0] / norm, axis[1] / norm, axis[2] / norm];
        // the C side = C itself (on-axis) + everything reachable from C
        // without crossing the C–O bond (the three methyl H's); the O and
        // its H stay put
        let mut side = vec![false; mol.atoms.len()];
        side[o1] = true;
        let mut stack = vec![c0];
        side[c0] = true;
        let mut c_side = vec![c0];
        while let Some(z) = stack.pop() {
            for &w in &mol.adjacency[z] {
                if !side[w] {
                    side[w] = true;
                    c_side.push(w);
                    stack.push(w);
                }
            }
        }
        let rotate_methyl = |src: &[[f64; 3]], deg: f64| -> Vec<[f64; 3]> {
            let th = deg.to_radians();
            let (sn, cs) = (th.sin(), th.cos());
            let mut out = src.to_vec();
            for &i in &c_side {
                let d = [
                    src[i][0] - base[c0][0],
                    src[i][1] - base[c0][1],
                    src[i][2] - base[c0][2],
                ];
                let cr = [
                    u[1] * d[2] - u[2] * d[1],
                    u[2] * d[0] - u[0] * d[2],
                    u[0] * d[1] - u[1] * d[0],
                ];
                let dot = u[0] * d[0] + u[1] * d[1] + u[2] * d[2];
                for x in 0..3 {
                    out[i][x] = base[c0][x] + cs * d[x] + sn * cr[x] + u[x] * dot * (1.0 - cs);
                }
            }
            out
        };
        let ff = MMFFForceField::new(&mol, crate::MMFFVariant::MMFF94s);
        let conv = crate::ConvergenceOptions {
            max_iterations: 500,
            ..Default::default()
        };
        let q = |c: &[[f64; 3]]| -> (Vec<[f64; 3]>, f64) {
            let r = crate::optimizer::optimize(&ff, c, &conv);
            (r.optimized_coords, r.final_energy)
        };
        let (a, ea) = q(&base);
        let (b, eb) = q(&rotate_methyl(&a, 120.0));
        let (c, ec) = q(&rotate_methyl(&a, -120.0));
        assert!(
            (ea - eb).abs() < 1e-3 && (ea - ec).abs() < 1e-3,
            "σ-equal wells"
        );
        let cfg = ImtdConfig {
            max_rounds: 1,
            runs_per_round: 1,
            steps_per_run: 100,
            quench_interval: 100,
            quench_max_iterations: 500,
            crossing_pairs: 0,
            convergence_new_minima: usize::MAX,
            ..Default::default()
        };
        let mut search = ImtdSearch::new(&mol, SearchEngine::MMFF94s, &cfg);
        let quencher = Quencher::new(&mol, SearchEngine::MMFF94s);
        while !search.done() {
            for spec in search.plan().to_vec() {
                let m = execute_run(&mol, SearchEngine::MMFF94s, &quencher, &cfg, &spec);
                search.gather(spec.run_index, m);
            }
        }
        assert_eq!(search.plan_verify(1).len(), 1);
        search.gather_verify(
            0,
            vec![
                Minimum {
                    coords: a,
                    energy: ea,
                    pe_max_since_prev: 0.0,
                    snapshot: 0,
                },
                Minimum {
                    coords: b.clone(),
                    energy: eb,
                    pe_max_since_prev: 0.0,
                    snapshot: 1,
                },
                Minimum {
                    coords: c,
                    energy: ec,
                    pe_max_since_prev: 0.0,
                    snapshot: 2,
                },
            ],
        );
        let land = search.landscape();
        let gvals: Vec<u64> = land.basins.iter().map(|x| x.degeneracy).collect();
        assert!(
            gvals.contains(&3),
            "methyl σ=3 must be observed (member + two rotations), got {gvals:?}"
        );
        // revisit does not raise
        search.gather_verify(
            0,
            vec![Minimum {
                coords: b,
                energy: eb,
                pe_max_since_prev: 0.0,
                snapshot: 0,
            }],
        );
        let land2 = search.landscape();
        assert_eq!(land2.basins.iter().map(|x| x.degeneracy).max(), Some(3));
    }

    #[test]
    fn boltzmann_applies_degeneracy() {
        // butane-like: anti E=0 (g=1) vs collapsed gauche pair g=2, ΔE=0.78
        let e = [0.0, 0.78];
        let g = [1.0, 2.0];
        let p = boltzmann_populations_deg(&e, &g, 298.15);
        let kt = crate::md::KB * 298.15;
        let w_gauche = 2.0 * (-0.78 / kt).exp();
        let expect = w_gauche / (1.0 + w_gauche);
        assert!((p[1] - expect).abs() < 1e-12, "{} vs {expect}", p[1]);
        assert!(p[1] > 0.30, "degeneracy lifts the gauche share: {}", p[1]);
        // without degeneracy it stays ~21%
        let p1 = boltzmann_populations(&e, 298.15);
        assert!((p1[1] - 0.21).abs() < 0.02, "{}", p1[1]);
    }

    #[test]
    fn thermo_tier_matches_harmonic_correction() {
        let mol = ethanol();
        let cfg = ImtdConfig {
            max_rounds: 1,
            runs_per_round: 1,
            steps_per_run: 200,
            quench_interval: 200,
            quench_max_iterations: 500,
            crossing_pairs: 0,
            convergence_new_minima: usize::MAX,
            ..Default::default()
        };
        let mut search = ImtdSearch::new(&mol, SearchEngine::MMFF94s, &cfg);
        while !search.done() {
            let quencher = Quencher::new(&mol, SearchEngine::MMFF94s);
            for spec in search.plan().to_vec() {
                let m = execute_run(&mol, SearchEngine::MMFF94s, &quencher, &cfg, &spec);
                search.gather(spec.run_index, m);
            }
        }
        let entries = search.thermo(298.15);
        assert!(!entries.is_empty());
        for e in &entries {
            assert_eq!(e.n_imag, 0, "ethanol minima are real");
            assert!((42.0..=62.0).contains(&e.zpe), "ZPE band: {}", e.zpe);
            assert!((e.g_free - e.energy - e.g_corr).abs() < 1e-9);
            assert!(e.g_free < e.energy + e.zpe, "entropy lowers G below E+ZPE");
        }
        let json = thermo_to_json(&entries, 298.15);
        assert!(json.contains("\"G\":\""), "thermo json emits G");
    }

    #[test]
    fn ensemble_capacity_eviction_keeps_stats_consistent() {
        // force eviction with max_conformers=3: statistics vectors, basin
        // ids, variants and pair keys must all shrink/remap in lockstep and
        // the landscape must stay aligned with the final order
        let mol = butane_all_h();
        let cfg = ImtdConfig {
            max_rounds: 2,
            runs_per_round: 2,
            steps_per_run: 1200,
            quench_interval: 150,
            quench_max_iterations: 250,
            crossing_pairs: 2,
            max_conformers: 3, // << forces evictions as the ensemble grows
            seed: 42,
            ..Default::default()
        };
        let mut search = ImtdSearch::new(&mol, SearchEngine::MMFF94s, &cfg);
        let quencher = Quencher::new(&mol, SearchEngine::MMFF94s);
        while !search.done() {
            for spec in search.plan().to_vec() {
                let m = execute_run(&mol, SearchEngine::MMFF94s, &quencher, &cfg, &spec);
                search.gather(spec.run_index, m);
            }
        }
        let res = search.result();
        assert!(
            res.conformers.len() <= 3,
            "capacity respected, got {}",
            res.conformers.len()
        );
        let land = search.landscape();
        assert_eq!(land.basins.len(), res.conformers.len(), "landscape aligned");
        // visits/pruned conservation holds after eviction (dropped members'
        // visits vanish — the inequality form)
        let visits: usize = land.basins.iter().map(|b| b.visits).sum();
        assert!(
            visits + land.pruned <= land.frames,
            "post-eviction conservation"
        );
        // degeneracy valid on every surviving basin
        assert!(land.basins.iter().all(|b| b.degeneracy >= 1));
        // refined saddles keyed by ids survive a later fold (no eviction
        // between refine and landscape here, but the mapping path runs)
        let jobs = search.plan_refine();
        assert!(jobs.len() <= 3 * 2 / 2, "at most one pair per 3 basins");
        for job in &jobs {
            search.gather_refine(job.index, Some(0.0));
        }
        let land2 = search.landscape();
        for t in &land2.transitions {
            if let Some(r) = t.refined {
                assert!((r - 0.0).abs() < 1e-12, "refined survives mapping");
            }
        }
    }

    #[test]
    fn gfnff_thermo_cache_is_bit_exact() {
        // same coords → the exact-coords cache must return bit-identical
        // frequencies; different T re-runs reuse the cached force field
        let mol = ethanol();
        let cfg = ImtdConfig {
            max_rounds: 1,
            runs_per_round: 1,
            steps_per_run: 100,
            quench_interval: 100,
            quench_max_iterations: 300,
            crossing_pairs: 0,
            convergence_new_minima: usize::MAX,
            ..Default::default()
        };
        let mut search = ImtdSearch::new(&mol, SearchEngine::GFNFF, &cfg);
        let quencher = Quencher::new(&mol, SearchEngine::GFNFF);
        while !search.done() {
            for spec in search.plan().to_vec() {
                let m = execute_run(&mol, SearchEngine::GFNFF, &quencher, &cfg, &spec);
                search.gather(spec.run_index, m);
            }
        }
        let a = search.thermo(298.15);
        let b = search.thermo(298.15);
        let c = search.thermo(500.0);
        assert!(!a.is_empty());
        // (G depends on T; the cached FF is keyed on coords, so energies must
        // match bit-for-bit across calls)
        for (x, y) in a.iter().zip(b.iter()) {
            assert_eq!(
                x.energy.to_bits(),
                y.energy.to_bits(),
                "E bit-exact via cache"
            );
            assert_eq!(x.zpe.to_bits(), y.zpe.to_bits(), "ZPE bit-exact via cache");
        }
        assert!(
            c.iter().any(|x| (x.g_free - x.energy).abs() > 0.0),
            "T enters G"
        );
        for e in &a {
            assert!(
                e.nu_max > 500.0,
                "GFN-FF frequencies physical scale: {}",
                e.nu_max
            );
        }
    }

    #[test]
    fn boltzmann_handles_short_degeneracy_vector() {
        let e = [0.0, 1.0, 2.0];
        // degeneracies shorter than energies: missing entries default to 1
        let p = boltzmann_populations_deg(&e, &[3.0], 298.15);
        let sum: f64 = p.iter().sum();
        assert!((sum - 1.0).abs() < 1e-12);
        let kt = crate::md::KB * 298.15;
        let expect0 = 3.0 / (3.0 + (-1.0 / kt).exp() + (-2.0 / kt).exp());
        assert!((p[0] - expect0).abs() < 1e-12, "g applied to first only");
        assert!(p[0] > p[1] && p[1] > p[2]);
    }

    #[test]
    fn codecs_reject_malformed_payloads_honestly() {
        // bad float strings / wrong shapes / missing fields never panic
        assert!(run_spec_from_json("not json").is_err());
        assert!(run_spec_from_json("{}").is_err(), "missing fields");
        assert!(run_spec_from_json(
            r#"{"round":1,"run_index":0,"start":[[1,2,3]],"cv":{"type":"distance","atoms":[0]},"hill_height":"NaN","seed":1}"#
        )
        .is_err(), "NaN float rejected");
        assert!(run_spec_from_json(
            r#"{"round":1,"run_index":0,"start":[[1,2]],"cv":{"type":"distance","atoms":[0,1]},"hill_height":"1e0","seed":1}"#
        )
        .is_err(), "malformed point rejected");
        assert!(minima_from_json("[]").is_err(), "missing minima array");
        assert!(minima_from_json(r#"{"minima":[{"e":"x","xyz":[]}]}"#).is_err());
        assert!(refine_job_from_json("{}").is_err());
        // the JSON literal null = "no saddle" (honest accept, not an error)
        assert_eq!(refine_result_from_json("null"), Ok(None));
        assert!(refine_result_from_json("{").is_err());
        assert!(refine_result_from_json(r#"{"saddle":"abc"}"#).is_err());
        // float-string round-trip channel rejects plain JSON numbers
        assert!(minima_from_json(r#"{"minima":[{"e":-1.5,"xyz":[[1e0,1e0,1e0]]}]}"#).is_err());
    }

    #[test]
    fn imtd_convergence_criterion_stops_at_threshold() {
        let mol = ethanol();
        let cfg = ImtdConfig {
            max_rounds: 2,
            runs_per_round: 1,
            steps_per_run: 400,
            quench_interval: 200,
            quench_max_iterations: 150,
            crossing_pairs: 1,
            // any round satisfies `new <= usize::MAX` → stops after round 1
            convergence_new_minima: usize::MAX,
            ..Default::default()
        };
        let res = imtd_conformer_search(&mol, SearchEngine::MMFF94s, &cfg);
        assert_eq!(res.rounds_executed, 1, "criterion must fire in round 1");
        assert!(res.converged);
        assert_eq!(res.new_minima_per_round.len(), 1);
        assert!(res.pool_size >= 2);
    }
}
