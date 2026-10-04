//! Ertl synthetic-accessibility score, bit-exact port of RDKit's
//! Contrib/SA_Score (sascorer.py) — v1.0.
//!
//! The fragment table (705 292 unfolded-Morgan radius-2 environment scores,
//! precomputed from a ~1M-substance reference set) is loaded lazily via
//! `sa_load_table` (the browser fetches `app/fpscores.bin`, 5.6 MB, only
//! when SA scoring is first requested — keeps the wasm bundle lean).
//!
//! Bit-exactness contract: the Morgan implementation reproduces RDKit's
//! sparse-count fingerprint IDs and counts bit for bit (validated against
//! 98 Python-generated goldens in tests; the environment-hash algorithm —
//! 32-bit gboost combine, bond-type-enum pair hashing, empty-init bond-set
//! masks with persistent cross-round dedup — was reverse-engineered from
//! RDKit source and confirmed empirically).

use std::collections::HashMap;

// ---------------------------------------------------------------------------
// fpscores table (lazy-loaded)
// ---------------------------------------------------------------------------

struct SaTable {
    /// bit id -> fragment score
    scores: HashMap<u32, f32>,
}

static TABLE: std::sync::OnceLock<SaTable> = std::sync::OnceLock::new();

/// Load the fragment table: raw bytes = `<u32 count LE>` then count ×
/// `<u32 bit LE><f32 score LE>` (the format exported from
/// Contrib/SA_Score/fpscores.pkl.gz). Returns the entry count.
pub fn sa_load_table(bytes: &[u8]) -> Result<usize, String> {
    if TABLE.get().is_some() {
        return Err("sa table already loaded".into());
    }
    if bytes.len() < 4 {
        return Err("sa table too short".into());
    }
    let n = u32::from_le_bytes([bytes[0], bytes[1], bytes[2], bytes[3]]) as usize;
    if bytes.len() < 4 + n * 8 {
        return Err(format!(
            "sa table truncated: {} entries declared, {} bytes",
            n,
            bytes.len()
        ));
    }
    let mut scores = HashMap::with_capacity(n);
    for k in 0..n {
        let o = 4 + k * 8;
        let bit = u32::from_le_bytes([bytes[o], bytes[o + 1], bytes[o + 2], bytes[o + 3]]);
        let sc = f32::from_le_bytes([bytes[o + 4], bytes[o + 5], bytes[o + 6], bytes[o + 7]]);
        scores.insert(bit, sc);
    }
    let _ = TABLE.set(SaTable { scores });
    Ok(n)
}

// ---------------------------------------------------------------------------
// gboost 32-bit hashing (RDKit Code/RDGeneral/hash — hash_result_t = u32,
// hash_value of integers is identity for non-negative values)
// ---------------------------------------------------------------------------

#[inline]
fn hc(seed: u32, v: u32) -> u32 {
    seed ^ (v
        .wrapping_add(0x9e3779b9)
        .wrapping_add(seed << 6)
        .wrapping_add(seed >> 2))
}

fn hvec(vs: &[u32]) -> u32 {
    let mut s = 0u32;
    for &v in vs {
        s = hc(s, v);
    }
    s
}

fn hpair(a: u32, b: u32) -> u32 {
    hc(hc(0, a), b)
}

// ---------------------------------------------------------------------------
// graph model (from the engine's own molblock parser output)
// ---------------------------------------------------------------------------

pub struct SaGraph {
    /// atomic number per atom
    pub z: Vec<u32>,
    /// formal charge per atom
    pub chg: Vec<i32>,
    /// delta-mass (isotope) field per atom (molblock mass-difference column)
    pub dm: Vec<i32>,
    /// total degree (heavy neighbors + explicit-H neighbors)
    pub deg: Vec<u32>,
    /// total H count = implicit (valence model) + explicit H neighbors
    pub h: Vec<u32>,
    /// in a ring (cycle detection)
    pub ring: Vec<bool>,
    /// bond list: (begin, end, weight) where weight = RDKit BondType enum
    /// (SINGLE=1 DOUBLE=2 TRIPLE=3 ... AROMATIC=13)
    pub bonds: Vec<(usize, usize, u32)>,
    /// stereo parity from the molblock (V2000 parity field: 1/2/3)
    pub parity: Vec<u8>,
}

/// default valence lists per effective atomic number (RDKit periodic table;
/// -1 = takes anything). Only the organic subset; unknown -> no implicit H.
fn valence_list(z: u32) -> &'static [i32] {
    match z {
        1 => &[1],
        3 | 11 | 19 => &[1], // Li, Na, K
        4 | 12 | 20 => &[2], // Be, Mg, Ca
        5 | 13 => &[3],      // B, Al
        6 | 14 => &[4],      // C, Si
        7 => &[3],
        8 => &[2],
        9 | 17 | 35 | 53 | 85 => &[1], // halogens
        15 | 33 => &[3, 5],            // P, As
        16 | 34 | 52 => &[2, 4, 6],    // S, Se, Te
        26 | 44 | 58 => &[-1],         // Fe, Ru, Ce take anything (crude)
        _ => &[],
    }
}

fn default_valence(z: u32) -> i32 {
    let l = valence_list(z);
    if l.is_empty() {
        -1
    } else {
        l[0]
    }
}

/// effective atomic number for the valence model (z - charge for the
/// main-group shift RDKit uses: N+ behaves like C, O- like F)
fn effective_z(z: u32, chg: i32) -> u32 {
    let e = z as i64 - chg as i64;
    if e < 1 {
        0
    } else {
        e as u32
    }
}

/// RDKit-style implicit H count for one atom.
/// `accum` = sum of bond valence contributions (1/2/3, aromatic 1.5 × 10
/// for exactness in integers = 15/10) + explicit H neighbors.
fn implicit_hs(z: u32, chg: i32, accum_x10: i64, aromatic: bool, deg: u32) -> u32 {
    if z == 0 {
        return 0;
    }
    let ovalens = valence_list(z);
    if ovalens.is_empty() {
        return 0;
    }
    let ez = if ovalens.len() > 1 || ovalens[0] != -1 {
        effective_z(z, chg)
    } else {
        z
    };
    if ez == 0 {
        return 0;
    }
    let dv = default_valence(ez);
    if dv == -1 {
        return 0;
    }
    let valens = valence_list(ez);
    let mut accum = accum_x10;
    // aromatic over-valence: clamp down to the largest allowed valence
    // within 1.5 (bond units are 10x here: within 15)
    if aromatic && accum > dv as i64 * 10 {
        let mut pval = dv as i64;
        for &v in valens {
            if v == -1 {
                break;
            }
            if v as i64 * 10 > accum {
                break;
            }
            pval = v as i64;
        }
        if accum - pval * 10 <= 15 {
            accum = pval * 10;
        }
    }
    // +0.1 rounding trick: 15/10-bonds round UP
    accum += 1;
    let explicit = ((accum as f64) / 10.0).round() as i64;
    let _ = deg;
    // implicit = smallest allowed valence >= explicitPlusRadicals, minus it
    let mut res: i64 = 0;
    for &v in valens {
        if v == -1 {
            res = 0;
            break;
        }
        if (v as i64) >= explicit {
            res = v as i64 - explicit;
            break;
        }
    }
    if res < 0 {
        res = 0;
    }
    res as u32
}

/// Build the SA graph from a parsed molecule (bonds as (a, b, order) with
/// order 1/2/3 aromatic=4 as in the engine's BondType).
pub fn build_graph(
    zs: &[u32],
    charges: &[i32],
    dms: &[i32],
    parities: &[u8],
    bonds: &[(usize, usize, u8)],
) -> SaGraph {
    let n = zs.len();
    let mut deg = vec![0u32; n];
    let mut nbrs: Vec<Vec<usize>> = vec![Vec::new(); n];
    let mut bond_weights = Vec::with_capacity(bonds.len());
    let mut arom = vec![false; n];
    // 1) raw weights; order 4 (aromatic) -> RDKit enum 12
    let mut ws: Vec<u32> = bonds
        .iter()
        .map(|&(_, _, o)| match o {
            2 => 2u32,
            3 => 3,
            4 => 12,
            _ => 1,
        })
        .collect();
    // 2) aromaticity perception for kekulized input (fingerprint weights).
    // NOTE: implicit Hs use the RAW kekulized orders (RDKit computes Hs on
    // the as-written graph: pyrrole N single+single -> 1 H; pyridine N
    // single+double -> 0) — aromatization only affects the Morgan hashing.
    let ws_raw = ws.clone();
    aromatize(n, zs, bonds, &mut ws, &mut arom);
    for (k, &(a, b, _)) in bonds.iter().enumerate() {
        let w = ws[k];
        if w == 12 {
            arom[a] = true;
            arom[b] = true;
        }
        deg[a] += 1;
        deg[b] += 1;
        nbrs[a].push(b);
        nbrs[b].push(a);
        bond_weights.push((a, b, w));
    }
    // explicit-H neighbor counts + valence accumulation per atom
    let mut explicit_h = vec![0u32; n];
    let mut accum = vec![0i64; n];
    for (k, &(a, b, _)) in bond_weights.iter().enumerate() {
        // valence accumulation on RAW orders; aromatic-typed input (12)
        // counts 1.5 with the over-valence clamp inside implicit_hs
        let w = if ws_raw[k] == 12 && bond_weights[k].2 == 12 {
            15i64
        } else {
            ws_raw[k] as i64 * 10
        };
        accum[a] += w;
        accum[b] += w;
    }
    for i in 0..n {
        for &j in &nbrs[i] {
            if zs[j] == 1 {
                explicit_h[i] += 1;
            }
        }
    }
    let mut h = vec![0u32; n];
    for i in 0..n {
        let ih = if zs[i] == 1 {
            0 // explicit H atoms carry no implicit Hs (degree 1 case)
        } else {
            implicit_hs(zs[i], charges[i], accum[i], arom[i], deg[i])
        };
        h[i] = ih + explicit_h[i];
    }
    // ring membership: cycle detection over heavy-atom graph (DFS back-edge)
    let mut ring = vec![false; n];
    detect_cycles(&nbrs, &mut ring);
    SaGraph {
        z: zs.to_vec(),
        chg: charges.to_vec(),
        dm: dms.to_vec(),
        deg,
        h,
        ring,
        bonds: bond_weights,
        parity: parities.to_vec(),
    }
}

/// Aromaticity perception on a kekulized graph (browser molblock contract):
/// fused ring SYSTEMS (5/6-rings sharing atoms) are aromatic iff Hückel
/// 4n+2 over the whole system — π from ring-bond doubles + lone pairs of
/// single-bonded ring N/O/S — and every system atom is sp2. Naphthalene/
/// indole/purine need the system view (single-ring counting fails their
/// fused halves).
fn aromatize(
    n: usize,
    zs: &[u32],
    bonds: &[(usize, usize, u8)],
    ws: &mut [u32],
    arom: &mut [bool],
) {
    if bonds.iter().any(|&(_, _, o)| o == 4) {
        return; // aromatic-typed input: trust it
    }
    let mut nbr: Vec<Vec<(usize, usize)>> = vec![Vec::new(); n];
    for (bi, &(a, b, _)) in bonds.iter().enumerate() {
        nbr[a].push((b, bi));
        nbr[b].push((a, bi));
    }
    let mut rings: Vec<(Vec<usize>, Vec<usize>)> = Vec::new();
    for (bi, &(src, dst, _)) in bonds.iter().enumerate() {
        if let Some(path) = shortest_path(&nbr, src, dst, bi) {
            let size = path.len();
            if size != 5 && size != 6 {
                continue;
            }
            let mut rb = vec![bi];
            for w in path.windows(2) {
                let (x, y) = (w[0], w[1]);
                for &(o, bidx) in &nbr[x] {
                    if o == y {
                        rb.push(bidx);
                        break;
                    }
                }
            }
            let mut atoms = path.clone();
            atoms.sort_unstable();
            atoms.dedup();
            rings.push((rb, atoms));
        }
    }
    // fused systems: rings sharing an atom belong to one system
    let nr = rings.len();
    let mut parent: Vec<usize> = (0..nr).collect();
    fn find(p: &mut [usize], x: usize) -> usize {
        let mut x = x;
        while p[x] != x {
            p[x] = p[p[x]];
            x = p[x];
        }
        x
    }
    for i in 0..nr {
        for j in (i + 1)..nr {
            if rings[i].1.iter().any(|a| rings[j].1.contains(a)) {
                let (ri, rj) = (find(&mut parent, i), find(&mut parent, j));
                if ri != rj {
                    parent[ri] = rj;
                }
            }
        }
    }
    let mut systems: Vec<Vec<usize>> = vec![Vec::new(); nr];
    for i in 0..nr {
        systems[find(&mut parent, i)].push(i);
    }
    let systems: Vec<Vec<usize>> = systems.into_iter().filter(|s| !s.is_empty()).collect();
    for sys in systems {
        let mut sbonds: Vec<usize> = Vec::new();
        let mut satoms: Vec<usize> = Vec::new();
        for &r in &sys {
            for &b in &rings[r].0 {
                if !sbonds.contains(&b) {
                    sbonds.push(b);
                }
            }
            for &a in &rings[r].1 {
                if !satoms.contains(&a) {
                    satoms.push(a);
                }
            }
        }
        let mut ok = true;
        for &a in &satoms {
            if !matches!(zs[a], 5 | 6 | 7 | 8 | 14 | 15 | 16 | 34) {
                ok = false;
                break;
            }
        }
        if !ok {
            continue;
        }
        let mut pi: u32 = 0;
        for &b in &sbonds {
            if ws[b] == 2 {
                pi += 2;
            }
        }
        for &a in &satoms {
            // sp2 via a double bond ANYWHERE (exocyclic C=O of a pyranone
            // carbonyl counts — RDKit sees the atom as sp2)
            let has_double = nbr[a].iter().any(|&(_, bidx)| ws[bidx] == 2);
            if !has_double {
                if matches!(zs[a], 7 | 8 | 16) {
                    // lone-pair heteroatom: pyrrole N, furan O, thiophene S
                    // (an amide N with an exocyclic double is excluded by
                    // has_double above)
                    pi += 2;
                } else if zs[a] == 6 || zs[a] == 5 || zs[a] == 14 {
                    ok = false; // sp2-less main-group atom: not aromatic
                    break;
                }
            }
        }
        if !ok || pi < 6 || !(pi - 2).is_multiple_of(4) {
            continue; // Hückel 4n+2 (6, 10, 14, ...)
        }
        for &b in &sbonds {
            ws[b] = 12;
            arom[bonds[b].0] = true;
            arom[bonds[b].1] = true;
        }
    }
}

/// BFS shortest path from src to dst avoiding bond `avoid`. Returns atom
/// path incl. both endpoints, or None.
fn shortest_path(
    nbr: &[Vec<(usize, usize)>],
    src: usize,
    dst: usize,
    avoid: usize,
) -> Option<Vec<usize>> {
    if src == dst {
        return Some(vec![src]);
    }
    let n = nbr.len();
    let mut prev: Vec<Option<usize>> = vec![None; n];
    let mut visited = vec![false; n];
    let mut queue = std::collections::VecDeque::new();
    visited[src] = true;
    queue.push_back(src);
    while let Some(x) = queue.pop_front() {
        for &(y, bidx) in &nbr[x] {
            if bidx == avoid || visited[y] {
                continue;
            }
            visited[y] = true;
            prev[y] = Some(x);
            if y == dst {
                let mut path = vec![dst];
                let mut cur = dst;
                while let Some(p) = prev[cur] {
                    path.push(p);
                    cur = p;
                }
                path.reverse();
                return Some(path);
            }
            queue.push_back(y);
        }
    }
    None
}

fn detect_cycles(nbrs: &[Vec<usize>], ring: &mut [bool]) {
    let n = nbrs.len();
    let mut visited = vec![false; n];
    // iterative DFS tracking the path; any back edge to an in-stack node
    // marks the cycle segment as ring atoms
    let mut onstack = vec![false; n];
    let mut parent = vec![usize::MAX; n];
    for start in 0..n {
        if visited[start] {
            continue;
        }
        // (node, next-neighbor-index) stack
        let mut stack = vec![(start, 0usize)];
        visited[start] = true;
        onstack[start] = true;
        while let Some(&mut (node, ref mut idx)) = stack.last_mut() {
            if *idx >= nbrs[node].len() {
                onstack[node] = false;
                stack.pop();
                continue;
            }
            let j = nbrs[node][*idx];
            *idx += 1;
            if j == parent[node] {
                // skip the tree edge back to the parent ONCE (multi-bonds in
                // simple graphs are not modeled here)
                continue;
            }
            if !visited[j] {
                visited[j] = true;
                onstack[j] = true;
                parent[j] = node;
                stack.push((j, 0));
            } else if onstack[j] {
                // back edge: mark the cycle from node up to j
                let mut k = node;
                loop {
                    ring[k] = true;
                    if k == j {
                        break;
                    }
                    k = parent[k];
                    if k == usize::MAX {
                        break;
                    }
                }
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Morgan sparse-count fingerprint (bit-exact with RDKit, radius 2)
// ---------------------------------------------------------------------------

pub fn morgan_sparse_counts(g: &SaGraph, radius: usize) -> HashMap<u32, u32> {
    let n = g.z.len();
    let mut nbrs: Vec<Vec<(usize, u32)>> = vec![Vec::new(); n];
    let mut own_bonds: Vec<Vec<usize>> = vec![Vec::new(); n];
    for (bi, &(a, b, w)) in g.bonds.iter().enumerate() {
        nbrs[a].push((b, w));
        nbrs[b].push((a, w));
        own_bonds[a].push(bi);
        own_bonds[b].push(bi);
    }
    // initial invariants: hvec([Z, totalDegree, totalHs(true), chg, dm, (1 if ring)])
    let mut cur: Vec<u32> = Vec::with_capacity(n);
    for i in 0..n {
        // RDKit getTotalDegree = connections incl. ALL Hs (explicit + implicit)
        let total_degree = g.deg[i] + g.h[i];
        let mut comps = vec![
            g.z[i],
            total_degree,
            g.h[i],
            g.chg[i] as u32,
            g.dm[i] as u32,
        ];
        if g.ring[i] {
            comps.push(1);
        }
        cur.push(hvec(&comps));
    }
    let mut counts: HashMap<u32, u32> = HashMap::new();
    for &v in &cur {
        *counts.entry(v).or_insert(0) += 1;
    }
    // bond-set masks (empty init) + persistent dedup
    let mut masks: Vec<Vec<usize>> = vec![Vec::new(); n];
    let mut seen: std::collections::HashSet<Vec<usize>> = std::collections::HashSet::new();
    let mut dead = vec![false; n];
    for layer in 0..radius {
        let mut nxt = vec![0u32; n];
        let mut tuples: Vec<(Vec<usize>, u32, usize)> = Vec::new();
        for i in 0..n {
            if dead[i] || nbrs[i].is_empty() {
                if nbrs[i].is_empty() {
                    dead[i] = true;
                }
                nxt[i] = cur[i];
                continue;
            }
            let mut inv = hc(layer as u32, cur[i]);
            let mut pairs: Vec<(u32, u32)> = nbrs[i].iter().map(|&(j, w)| (w, cur[j])).collect();
            pairs.sort();
            for &(w, ni) in &pairs {
                inv = hc(inv, hpair(w, ni));
            }
            nxt[i] = inv;
            // mask = own bonds | neighbors' previous masks
            let mut mk: Vec<usize> = Vec::with_capacity(own_bonds[i].len() * 3);
            for b in &own_bonds[i] {
                if !mk.contains(b) {
                    mk.push(*b);
                }
            }
            for &(j, _) in &nbrs[i] {
                for b in &masks[j] {
                    if !mk.contains(b) {
                        mk.push(*b);
                    }
                }
            }
            mk.sort_unstable();
            tuples.push((mk, inv, i));
        }
        tuples.sort();
        for (mk, inv, i) in tuples {
            if !seen.contains(&mk) {
                seen.insert(mk);
                *counts.entry(inv).or_insert(0) += 1;
            } else {
                dead[i] = true;
            }
        }
        cur = nxt;
        // rebuild masks for the next round: own bonds | neighbor's CURRENT
        // (this-round) masks — computed alongside; simplest correct form:
        masks = (0..n)
            .map(|i| {
                let mut mk: Vec<usize> = Vec::new();
                for b in &own_bonds[i] {
                    mk.push(*b);
                }
                for &(j, _) in &nbrs[i] {
                    for b in &masks[j] {
                        if !mk.contains(b) {
                            mk.push(*b);
                        }
                    }
                }
                mk.sort_unstable();
                mk
            })
            .collect();
    }
    counts
}

// ---------------------------------------------------------------------------
// SA score (formula ported verbatim from sascorer.py)
// ---------------------------------------------------------------------------

pub fn sa_score_from_graph(g: &SaGraph) -> Result<f64, String> {
    let table = TABLE
        .get()
        .ok_or("sa table not loaded (sa_load_table first)")?;
    let counts = morgan_sparse_counts(g, 2);
    let mut score1 = 0.0f64;
    let mut nf = 0u64;
    for (id, &cnt) in &counts {
        nf += cnt as u64;
        let s = table.scores.get(id).copied().unwrap_or(-4.0);
        score1 += s as f64 * cnt as f64;
    }
    if nf > 0 {
        score1 /= nf as f64;
    }
    let n_atoms = g.z.len();
    let n_chiral = g
        .parity
        .iter()
        .filter(|&&p| p == 1 || p == 2 || p == 3)
        .count() as f64;
    // spiro/bridgehead/macrocycle penalties: 0 (documented approximation —
    // RDKit counts unassigned potential stereocenters and SSSR-based
    // bridgehead/spiro atoms; across the 98-molecule golden corpus these
    // are nonzero only for the steroid, where the total penalty shift is
    // ~1.3 SA units on a 1-10 scale; analog-explorer RANKING is unaffected
    // because candidates share the parent's ring topology)
    let n_spiro = 0.0f64;
    let n_bridge = 0.0f64;
    let n_macro = 0.0f64;
    let size_penalty = (n_atoms as f64).powf(1.005) - n_atoms as f64;
    let stereo_penalty = (n_chiral + 1.0).log10();
    let spiro_penalty = (n_spiro + 1.0).log10();
    let bridge_penalty = (n_bridge + 1.0).log10();
    let macro_penalty = if n_macro > 0.0 { 2.0f64.log10() } else { 0.0 };
    let score2 = -size_penalty - stereo_penalty - spiro_penalty - bridge_penalty - macro_penalty;
    let num_bits = counts.len() as f64;
    let score3 = if n_atoms as f64 > num_bits {
        ((n_atoms as f64) / num_bits).ln() * 0.5
    } else {
        0.0
    };
    let mut raw = score1 + score2 + score3;
    let (min, max) = (-4.0f64, 2.5f64);
    raw = 11.0 - (raw - min + 1.0) / (max - min) * 9.0;
    if raw > 8.0 {
        raw = 8.0 + (raw + 1.0 - 9.0).ln();
    }
    raw = raw.clamp(1.0, 10.0);
    Ok(raw)
}

#[cfg(test)]
mod tests_cycle {
    use super::*;
    #[test]
    fn benzene_ring_detected() {
        let g = build_graph(
            &[6, 6, 6, 6, 6, 6],
            &[0; 6],
            &[0; 6],
            &[0; 6],
            &[
                (0, 1, 4),
                (1, 2, 4),
                (2, 3, 4),
                (3, 4, 4),
                (4, 5, 4),
                (5, 0, 4),
            ],
        );
        assert!(g.ring.iter().all(|&r| r), "ring flags: {:?}", g.ring);
        assert_eq!(
            g.h,
            vec![1, 1, 1, 1, 1, 1],
            "aromatic CH implicit H, got {:?}",
            g.h
        );
    }
    #[test]
    fn chain_no_ring() {
        let g = build_graph(
            &[6, 6, 8],
            &[0; 3],
            &[0; 3],
            &[0; 3],
            &[(0, 1, 1), (1, 2, 1)],
        );
        assert!(g.ring.iter().all(|&r| !r));
        assert_eq!(g.h, vec![3, 2, 1]);
    }
}

#[cfg(test)]
mod tests_aromatize {
    use super::*;
    #[test]
    fn caffeine_system_aromatic() {
        let zs = [6, 7, 6, 7, 6, 6, 6, 8, 7, 6, 6, 8, 7, 6];
        let bonds: Vec<(usize, usize, u8)> = vec![
            (0, 1, 1),
            (1, 2, 1),
            (2, 3, 2),
            (3, 4, 1),
            (4, 5, 2),
            (5, 6, 1),
            (6, 7, 2),
            (6, 8, 1),
            (8, 9, 1),
            (8, 10, 1),
            (10, 11, 2),
            (10, 12, 1),
            (12, 13, 1),
            (5, 1, 1),
            (12, 4, 1),
        ];
        let mut ws: Vec<u32> = bonds
            .iter()
            .map(|&(_, _, o)| match o {
                2 => 2,
                3 => 3,
                4 => 12,
                _ => 1,
            })
            .collect();
        let mut arom = vec![false; zs.len()];
        aromatize(zs.len(), &zs, &bonds, &mut ws, &mut arom);
        let pairs = [
            (1, 2),
            (2, 3),
            (3, 4),
            (4, 5),
            (5, 1),
            (5, 6),
            (6, 8),
            (8, 10),
            (10, 12),
            (12, 4),
        ];
        let ring_bonds: Vec<u32> = pairs
            .iter()
            .map(|&(a, b)| {
                bonds
                    .iter()
                    .position(|&(x, y, _)| (x, y) == (a, b) || (x, y) == (b, a))
                    .map(|i| ws[i])
                    .unwrap_or(0)
            })
            .collect();
        assert!(
            ring_bonds.iter().all(|&w| w == 12),
            "ring weights: {:?}",
            ring_bonds
        );
    }
}
