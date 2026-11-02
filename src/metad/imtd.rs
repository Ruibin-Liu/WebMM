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

use crate::md::{MDConfig, MDRunner};
use crate::metad::{CollectiveVariable, DihedralCV, DistanceCV, MetaDConfig, MetaDynamics};
use crate::mmff::MMFFForceField;
use crate::molecule::graph::find_rings;
use crate::molecule::Molecule;
use crate::optimizer::{self};
use crate::ConvergenceOptions;
use std::collections::{HashSet, VecDeque};
use std::rc::Rc;

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

// ── geometry helpers ──

/// Proper rotation + translation (Horn quaternion Kabsch) mapping the mobile
/// point set `p` onto the reference `q` (minimizes Σ|R·p+t−q|²), plus the
/// achieved RMSD. Reflections are never returned (molecules are proper
/// objects). Coordinate-level sibling of the `kabsch_align` WASM export
/// (same math and convention: R·p+t ≈ q).
fn kabsch_transform(p: &[[f64; 3]], q: &[[f64; 3]]) -> ([f64; 9], [f64; 3], f64) {
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
fn aligned_rmsd_over(a: &[[f64; 3]], b: &[[f64; 3]], idx: &[usize]) -> f64 {
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
#[derive(Debug, Clone, Copy)]
enum CvSpec {
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

struct Ensemble {
    members: Vec<Conformer>,
    /// Atom indices used for dedup alignment (heavy atoms; all if none).
    dedup_idx: Vec<usize>,
    rmsd_tol: f64,
    energy_window: f64,
    max_conformers: usize,
    emin: f64,
}

impl Ensemble {
    fn new(mol: &Molecule, config: &ImtdConfig) -> Self {
        let mut dedup_idx = heavy_indices(mol);
        if dedup_idx.is_empty() {
            dedup_idx = (0..mol.atoms.len()).collect();
        }
        Ensemble {
            members: Vec::new(),
            dedup_idx,
            rmsd_tol: config.rmsd_tol,
            energy_window: config.energy_window,
            max_conformers: config.max_conformers.max(1),
            emin: f64::INFINITY,
        }
    }

    /// Indices of the `PARENT_POOL` lowest-energy members (deterministic).
    fn lowest_energy_indices(&self) -> Vec<usize> {
        let mut order: Vec<usize> = (0..self.members.len()).collect();
        order.sort_by(|&x, &y| self.members[x].energy.total_cmp(&self.members[y].energy));
        order.truncate(PARENT_POOL);
        order
    }

    /// Quench `coords` to its local minimum and try to add it. Returns the
    /// accepted conformer's coordinates (`None` if duplicate, out of energy
    /// window, cap-rejected, or the quench produced a non-finite energy).
    fn try_add(
        &mut self,
        ff: &MMFFForceField,
        coords: &[[f64; 3]],
        round: usize,
        origin: Origin,
        quench_max_iter: usize,
    ) -> Option<Vec<[f64; 3]>> {
        let conv = ConvergenceOptions {
            max_iterations: quench_max_iter.max(1),
            ..Default::default()
        };
        let r = optimizer::optimize(ff, coords, &conv);
        let energy = r.final_energy;
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
        // duplicate minimum? (aligned heavy-atom RMSD below tolerance)
        let opt = r.optimized_coords;
        for m in &self.members {
            if aligned_rmsd_over(&opt, &m.coords, &self.dedup_idx) < self.rmsd_tol {
                return None;
            }
        }
        if energy < self.emin {
            self.emin = energy;
        }
        self.members.push(Conformer {
            coords: opt,
            energy,
            round_found: round,
            origin,
        });
        if self.members.len() > self.max_conformers {
            // evict the previous maximum (never the just-added member: it
            // passed the strictly-below-max check above)
            let mut i_max = 0;
            for (i, m) in self.members.iter().enumerate() {
                if m.energy > self.members[i_max].energy {
                    i_max = i;
                }
            }
            self.members.remove(i_max);
        }
        self.members.last().map(|m| m.coords.clone())
    }
}

// ── short metad run ──

/// Bundles the immutable per-search context for one short metad run.
struct MetadRun<'a> {
    mol: &'a Molecule,
    variant: crate::MMFFVariant,
    config: &'a ImtdConfig,
}

impl MetadRun<'_> {
    /// Run one short metadynamics trajectory from `start`, calling
    /// `on_snapshot` with the coordinates every `quench_interval` steps and
    /// at the final frame.
    fn run(
        &self,
        start: &[[f64; 3]],
        cv: CvSpec,
        hill_height: f64,
        seed: u64,
        mut on_snapshot: impl FnMut(&[[f64; 3]]),
    ) {
        let config = self.config;
        let mut mol_start = self.mol.clone();
        for (atom, c) in mol_start.atoms.iter_mut().zip(start.iter()) {
            atom.position = *c;
        }
        let width = match cv {
            CvSpec::Dihedral(..) => config.hill_width,
            CvSpec::Distance(..) => config.hill_width_distance,
        };
        let ff = MMFFForceField::new(&mol_start, self.variant);
        let metad = Rc::new(MetaDynamics::new(
            Box::new(ff),
            cv.build(),
            MetaDConfig {
                hill_height,
                hill_width: width,
                deposit_interval: config.deposit_interval.max(1),
                bias_factor: config.bias_factor,
                temperature_k: config.temperature_k,
            },
        ));
        let mut runner = MDRunner::from_molecule(
            metad,
            &mol_start,
            MDConfig {
                dt_fs: config.dt_fs,
                temperature_k: config.temperature_k,
                friction_per_ps: config.friction_per_ps,
                seed,
            },
        );
        let interval = config.quench_interval.max(1);
        for i in 1..=config.steps_per_run {
            runner.step();
            if i.is_multiple_of(interval) {
                on_snapshot(runner.coords());
            }
        }
        on_snapshot(runner.coords()); // final frame
    }
}

// ── the iMTD-GC scheduler ──

/// iMTD-GC conformational search (see the module docs). Deterministic: same
/// molecule + variant + config ⇒ bit-identical result.
pub fn imtd_conformer_search(
    mol: &Molecule,
    variant: crate::MMFFVariant,
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
    let ff = MMFFForceField::new(mol, variant);
    let mut ensemble = Ensemble::new(mol, config);

    // round 0 seed: the quenched input structure
    let input_coords: Vec<[f64; 3]> = mol.atoms.iter().map(|a| a.position).collect();
    ensemble.try_add(
        &ff,
        &input_coords,
        0,
        Origin::Input,
        config.quench_max_iterations,
    );

    let rot = rotatable_bonds(mol);
    let pool = cv_pool(mol, &rot);
    let n_rotatable = rot.len();
    if pool.is_empty() || config.max_rounds == 0 || config.steps_per_run == 0 {
        return finish(ensemble, 0, true, Vec::new(), pool.len(), n_rotatable);
    }

    let fragments: Vec<(Vec<usize>, usize)> = rot
        .iter()
        .map(|&(j, k)| (fragment_split(mol, j, k).0, k))
        .collect();

    let mut new_minima_per_round = Vec::new();
    let mut crossed_seeds: Vec<Vec<[f64; 3]>> = Vec::new();
    let mut rounds_executed = 0;
    let mut converged = false;

    for r in 0..config.max_rounds {
        rounds_executed = r + 1;
        let mut new_this_round = 0usize;

        // — collect: short metad runs, rotating CV and bias strength —
        for m in 0..config.runs_per_round {
            let g = r * config.runs_per_round + m; // global run counter
            let cv = pool[g % pool.len()];
            let hill = config.hill_height * HEIGHT_LADDER[g % HEIGHT_LADDER.len()];
            // start structure: fresh crossed seeds first, then rotate
            // through the low-energy ensemble
            let start: Vec<[f64; 3]> = if let Some(seed_coords) = crossed_seeds.first().cloned() {
                crossed_seeds.remove(0);
                seed_coords
            } else {
                let top = ensemble.lowest_energy_indices();
                ensemble.members[top[g % top.len()]].coords.clone()
            };
            let seed = config.seed + (r as u64) * 1000 + m as u64;
            MetadRun {
                mol,
                variant,
                config,
            }
            .run(&start, cv, hill, seed, |coords| {
                if ensemble
                    .try_add(
                        &ff,
                        coords,
                        r + 1,
                        Origin::Metad,
                        config.quench_max_iterations,
                    )
                    .is_some()
                {
                    new_this_round += 1;
                }
            });
        }

        // — genetic crossing: recombine low-energy parents at rotatable cuts —
        if config.crossing_pairs > 0 && ensemble.members.len() >= 2 && !fragments.is_empty() {
            let top = ensemble.lowest_energy_indices();
            let n = top.len();
            let pair_budget = config.crossing_pairs.min(n * (n - 1) / 2);
            for p in 0..pair_budget {
                let a = p % n;
                let b = (a + 1 + p / n) % n;
                if a == b {
                    continue;
                }
                let (f1, k) = &fragments[p % fragments.len()];
                let child = genetic_crossing(
                    &ensemble.members[top[a]].coords,
                    &ensemble.members[top[b]].coords,
                    f1,
                    *k,
                );
                if let Some(seed_coords) = ensemble.try_add(
                    &ff,
                    &child,
                    r + 1,
                    Origin::Crossing,
                    config.quench_max_iterations,
                ) {
                    new_this_round += 1;
                    crossed_seeds.push(seed_coords);
                }
            }
        }

        new_minima_per_round.push(new_this_round);
        if new_this_round <= config.convergence_new_minima {
            converged = true;
            break;
        }
    }

    finish(
        ensemble,
        rounds_executed,
        converged,
        new_minima_per_round,
        pool.len(),
        n_rotatable,
    )
}

/// Sort by energy, re-apply the window against the FINAL minimum (the minimum
/// may have dropped during the search), cap the ensemble size.
fn finish(
    mut ensemble: Ensemble,
    rounds_executed: usize,
    converged: bool,
    new_minima_per_round: Vec<usize>,
    pool_size: usize,
    n_rotatable: usize,
) -> ImtdResult {
    let window = ensemble.energy_window;
    let cap = ensemble.max_conformers;
    ensemble
        .members
        .sort_by(|x, y| x.energy.total_cmp(&y.energy));
    if let Some(first) = ensemble.members.first() {
        let emin = first.energy;
        ensemble.members.retain(|m| m.energy <= emin + window);
    }
    ensemble.members.truncate(cap);
    ImtdResult {
        conformers: ensemble.members,
        rounds_executed,
        converged,
        new_minima_per_round,
        pool_size,
        n_rotatable,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

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
        let res = imtd_conformer_search(&mol, crate::MMFFVariant::MMFF94s, &cfg);
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
        let res2 = imtd_conformer_search(&mol, crate::MMFFVariant::MMFF94s, &cfg);
        assert_eq!(res.conformers.len(), res2.conformers.len());
        for (x, y) in res.conformers.iter().zip(res2.conformers.iter()) {
            assert_eq!(x.energy, y.energy, "energies must be bit-identical");
            assert_eq!(x.coords, y.coords, "coordinates must be bit-identical");
        }
    }

    #[test]
    fn imtd_rigid_methane_single_conformer() {
        let mol = crate::molecule::parser::parse_sdf(METHANE_SDF).expect("parse");
        let res = imtd_conformer_search(&mol, crate::MMFFVariant::MMFF94s, &ImtdConfig::default());
        assert_eq!(res.pool_size, 0, "methane has no biasable CV");
        assert_eq!(res.conformers.len(), 1);
        assert!(res.converged);
        assert_eq!(res.rounds_executed, 0);
        assert_eq!(res.conformers[0].origin, Origin::Input);
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
        let res = imtd_conformer_search(&mol, crate::MMFFVariant::MMFF94s, &cfg);
        assert_eq!(res.rounds_executed, 1, "criterion must fire in round 1");
        assert!(res.converged);
        assert_eq!(res.new_minima_per_round.len(), 1);
        assert!(res.pool_size >= 2);
    }
}
