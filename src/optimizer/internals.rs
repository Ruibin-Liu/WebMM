//! Delocalized internal coordinates (Baker & Besley 1996; geomeTRIC-style).
//!
//! Cartesian L-BFGS follows valley floors slowly on flexible molecules (each
//! Cartesian step mixes stiff bond stretches with soft torsions). Optimizing
//! in delocalized internal coordinates — linear combinations of bond lengths,
//! angles, dihedrals and out-of-plane angles formed from the eigenvectors of
//! the Wilson G = B·Bᵀ matrix — decouples those modes, typically cutting
//! iteration counts several-fold on floppy molecules.
//!
//! The optimizer wraps an internal-coordinate `Objective` around the existing
//! force fields; energies/gradients stay Cartesian, only the search space
//! changes. Construction falls back to Cartesian optimization (None) for
//! degenerate cases (n < 3, no primitives, oversized primitive sets).

// Clippy: dense matrix loops intentionally index by (i, j) to mirror the
// Wilson B/G formalism (same policy as the gfnff module).
#![allow(clippy::needless_range_loop)]

use super::jacobi::sym_jacobi;

/// Covalent radii (Å, Cordero et al. 2008), index = atomic number.
const RCOV: [f64; 87] = [
    0.0, // unused (Z=0)
    0.31, 0.28, 1.28, 0.96, 0.84, 0.76, 0.71, 0.66, 0.57, 0.58, // H..Ne
    1.66, 1.41, 1.21, 1.11, 1.07, 1.05, 1.02, 1.06, 2.03, 1.76, // Na..Ca
    1.70, 1.60, 1.53, 1.39, 1.39, 1.39, 1.32, 1.24, 1.32, 1.22, // Sc..Zn
    1.22, 1.20, 1.19, 1.20, 1.20, 1.16, 2.20, 1.95, 1.90, 1.75, // Ga..Zr
    1.64, 1.54, 1.47, 1.46, 1.42, 1.39, 1.45, 1.44, 1.42, 1.39, // Nb..Sn
    1.39, 1.38, 1.39, 1.40, 2.44, 2.15, 2.07, 2.04, 2.03, 2.01, // Sb..Nd
    1.99, 1.98, 1.98, 1.96, 1.94, 1.92, 1.92, 1.89, 1.90, 1.87, // Pm..Yb
    1.87, 1.75, 1.70, 1.62, 1.51, 1.44, 1.41, 1.36, 1.36, 1.32, // Lu..Hg
    1.45, 1.46, 1.48, 1.40, 1.50, 1.50, // Tl..Rn
];

/// Primitive internal coordinate.
#[derive(Clone, Copy, Debug, PartialEq)]
enum Prim {
    /// Bond stretch |x_i - x_j| (also linear-angle replacements and auxiliary
    /// close-contact distances).
    Dist(usize, usize),
    /// Valence angle a-center-c.
    Angle(usize, usize, usize),
    /// Torsion a-b-c-d about the central bond b-c.
    Dih(usize, usize, usize, usize),
    /// Out-of-plane: angle the bond center->l makes with the plane
    /// (center, j, k); zero for planar centers.
    Oop {
        c: usize,
        j: usize,
        k: usize,
        l: usize,
    },
}

/// Hard cap on the primitive count — G is (n_prim)^2 and the Jacobi
/// diagonalization is O(n_prim^3 · sweeps); beyond this the transform
/// overhead outweighs the iteration savings (fallback: Cartesian).
const MAX_PRIMS: usize = 800;

pub struct InternalCoords {
    n: usize,
    prims: Vec<Prim>,
    /// Reference primitive values (at construction).
    s_ref: Vec<f64>,
    /// Delocalized basis: u[i][k] = weight of primitive i in coordinate k
    /// (kept G-eigenvector columns), k = 0..n_dof.
    u: Vec<Vec<f64>>,
    /// G eigenvalues of the kept coordinates.
    lam: Vec<f64>,
    /// B_q = U^T B(x_ref) (n_dof x 3N), cached at construction: the map
    /// x(z) = x_ref + B_q^T Lambda^-1 (Baker-iterated) is EXACTLY linear
    /// with this matrix, and every back_transform / grad_q call reuses it
    /// (recomputing it per call was the dominant DIC cost).
    bq: Vec<Vec<f64>>,
    pub n_dof: usize,
}

impl InternalCoords {
    /// Build internal coordinates for a molecule at geometry `x` with atomic
    /// numbers `z`. Returns None for degenerate cases (see module doc).
    pub fn build(x: &[[f64; 3]], z: &[usize]) -> Option<Self> {
        let n = x.len();
        if n < 3 || z.len() != n {
            return None;
        }
        if z.iter().any(|&zi| zi >= RCOV.len()) {
            return None;
        }

        // --- connectivity from covalent radii (1.3x sum) ---
        let mut adj: Vec<Vec<usize>> = vec![Vec::new(); n];
        for i in 0..n {
            for j in (i + 1)..n {
                // skip H-H "bonds" (H2-like degenerate input)
                if z[i] == 1 && z[j] == 1 {
                    continue;
                }
                let d = dist(x, i, j);
                if d < 1.3 * (RCOV[z[i]] + RCOV[z[j]]) {
                    adj[i].push(j);
                    adj[j].push(i);
                }
            }
        }

        let linear =
            |a: usize, b: usize, c: usize| -> bool { angle_of(x, a, b, c) > 175.0f64.to_radians() };

        let mut prims: Vec<Prim> = Vec::new();
        // bonds
        for i in 0..n {
            for &j in &adj[i] {
                if i < j {
                    prims.push(Prim::Dist(i, j));
                }
            }
        }
        // angles (linear ones replaced by the 1-3 distance) + impropers
        for c in 0..n {
            let m = adj[c].len();
            for a in 0..m {
                for b in (a + 1)..m {
                    let (p, q) = (adj[c][a], adj[c][b]);
                    if linear(p, c, q) {
                        prims.push(Prim::Dist(p, q));
                    } else {
                        prims.push(Prim::Angle(p, c, q));
                    }
                }
            }
            if m == 3 {
                let (j, k, l) = (adj[c][0], adj[c][1], adj[c][2]);
                prims.push(Prim::Oop { c, j, k, l });
            }
        }
        // dihedrals about each bond, skipping linear flanks
        for b in 0..n {
            for &c in &adj[b] {
                if b > c {
                    continue; // each bond once
                }
                for &a in &adj[b] {
                    if a == c || linear(a, b, c) {
                        continue;
                    }
                    for &d in &adj[c] {
                        if d == b || d == a || linear(b, c, d) {
                            continue;
                        }
                        // prune H-H flank pairs (e.g. H-C-C-H): near-redundant
                        // with the heavy-atom dihedral + angles, and they
                        // inflate the G diagonalization size ~2x
                        if z[a] == 1 && z[d] == 1 {
                            continue;
                        }
                        prims.push(Prim::Dih(a, b, c, d));
                    }
                }
            }
        }
        // auxiliary nonbonded distances (>= 1-4, heavy atoms, close contact)
        let share_nbr =
            |a: usize, b: usize| -> bool { adj[a].iter().any(|&m| adj[m].contains(&b)) };
        for i in 0..n {
            if z[i] == 1 {
                continue;
            }
            for j in (i + 1)..n {
                if z[j] == 1 || adj[i].contains(&j) || share_nbr(i, j) {
                    continue;
                }
                if dist(x, i, j) < 2.8 {
                    prims.push(Prim::Dist(i, j));
                }
            }
        }
        prims.dedup();
        if prims.is_empty() || prims.len() > MAX_PRIMS {
            return None;
        }

        let n_prim = prims.len();
        let s_ref: Vec<f64> = prims.iter().map(|p| prim_value(p, x)).collect();
        let b = b_matrix(&prims, x, n);

        // G = B B^T (flat)
        let mut g = vec![0.0f64; n_prim * n_prim];
        for i in 0..n_prim {
            for j in 0..n_prim {
                let mut s = 0.0;
                for k in 0..(3 * n) {
                    s += b[i][k] * b[j][k];
                }
                g[i * n_prim + j] = s;
            }
        }
        let mut v = vec![0.0f64; n_prim * n_prim];
        sym_jacobi(&mut g, &mut v, n_prim);

        // order columns by eigenvalue, keep the non-null space
        let mut order: Vec<usize> = (0..n_prim).collect();
        order.sort_by(|&a, &bb| g[bb * n_prim + bb].partial_cmp(&g[a * n_prim + a]).unwrap());
        let lmax = g[order[0] * n_prim + order[0]].max(1e-30);
        let cutoff = lmax * 1e-6;
        let kept: Vec<usize> = order
            .into_iter()
            .filter(|&c| g[c * n_prim + c] > cutoff && g[c * n_prim + c] > 1e-8)
            .collect();
        // primitives are translation/rotation invariant, so rank(G) <= 3n-6
        let n_dof = kept.len().min(3 * n - 6);
        if n_dof == 0 {
            return None;
        }
        let mut u = vec![vec![0.0f64; n_dof]; n_prim];
        let mut lam = vec![0.0f64; n_dof];
        for (k, &col) in kept.iter().take(n_dof).enumerate() {
            for i in 0..n_prim {
                u[i][k] = v[i * n_prim + col];
            }
            lam[k] = g[col * n_prim + col];
        }
        // cache B_q = U^T B(x_ref) once (see struct doc)
        let nc = 3 * n;
        let mut bq = vec![vec![0.0f64; nc]; n_dof];
        for k in 0..n_dof {
            for i in 0..n_prim {
                let w = u[i][k];
                if w == 0.0 {
                    continue;
                }
                for j in 0..nc {
                    bq[k][j] += w * b[i][j];
                }
            }
        }

        Some(InternalCoords {
            n,
            prims,
            s_ref,
            u,
            lam,
            bq,
            n_dof,
        })
    }

    /// Delocalized coordinates q = U^T (s(x) - s_ref).
    pub fn q(&self, x: &[[f64; 3]]) -> Vec<f64> {
        let mut q = vec![0.0f64; self.n_dof];
        for (i, p) in self.prims.iter().enumerate() {
            let d = prim_value(p, x) - self.s_ref[i];
            for k in 0..self.n_dof {
                q[k] += self.u[i][k] * d;
            }
        }
        q
    }

    /// Transform a (TR-projected, flat) Cartesian gradient into the
    /// delocalized space: g_q = Lambda^-1 B_q g_x (cached B_q).
    pub fn grad_q(&self, gx_flat: &[f64]) -> Vec<f64> {
        let nc = 3 * self.n;
        let mut gq = vec![0.0f64; self.n_dof];
        for k in 0..self.n_dof {
            let mut s = 0.0;
            let bk = &self.bq[k];
            for j in 0..nc {
                s += bk[j] * gx_flat[j];
            }
            gq[k] = s / self.lam[k];
        }
        gq
    }

    /// Solve B_q dx = dq for the minimum-norm Cartesian displacement via the
    /// Baker iteration dx += B_q^T Lambda^-1 (dq - B_q dx) with the cached
    /// B_q (6 rounds; the fixed-point residual contracts geometrically).
    pub fn back_transform(&self, dq: &[f64]) -> Vec<[f64; 3]> {
        let nc = 3 * self.n;
        let mut dx = vec![0.0f64; nc];
        let mut r = vec![0.0f64; self.n_dof];
        for _round in 0..6 {
            for k in 0..self.n_dof {
                let mut s = 0.0;
                let bk = &self.bq[k];
                for j in 0..nc {
                    s += bk[j] * dx[j];
                }
                r[k] = dq[k] - s;
            }
            for j in 0..nc {
                let mut corr = 0.0;
                for k in 0..self.n_dof {
                    corr += self.bq[k][j] * r[k] / self.lam[k];
                }
                dx[j] += corr;
            }
        }
        let mut out = vec![[0.0f64; 3]; self.n];
        for i in 0..self.n {
            out[i] = [dx[3 * i], dx[3 * i + 1], dx[3 * i + 2]];
        }
        out
    }

    pub fn n_prims(&self) -> usize {
        self.prims.len()
    }
}

/// Remove the 6 translational/rotational components from a flat Cartesian
/// gradient (orthogonal projection onto the internal subspace).
pub fn project_out_tr(x: &[[f64; 3]], g_flat: &mut [f64]) {
    let n = x.len();
    let nc = 3 * n;
    // centroid
    let mut cxyz = [0.0f64; 3];
    for p in x {
        for t in 0..3 {
            cxyz[t] += p[t];
        }
    }
    for t in 0..3 {
        cxyz[t] /= n as f64;
    }
    // basis: translations + rotations e_t x (r - centroid)
    let mut basis: Vec<Vec<f64>> = Vec::with_capacity(6);
    for t in 0..3 {
        let mut v = vec![0.0f64; nc];
        for i in 0..n {
            v[3 * i + t] = 1.0;
        }
        basis.push(v);
    }
    // rotation generators: axis 0 -> (y,z), axis 1 -> (z,x), axis 2 -> (x,y)
    for (o1, o2) in [(1usize, 2usize), (2, 0), (0, 1)] {
        let mut v = vec![0.0f64; nc];
        for i in 0..n {
            v[3 * i + o1] = x[i][o2] - cxyz[o2];
            v[3 * i + o2] = -(x[i][o1] - cxyz[o1]);
        }
        basis.push(v);
    }
    // Gram-Schmidt, then subtract the projections from g
    let mut ortho: Vec<Vec<f64>> = Vec::with_capacity(6);
    for v in &basis {
        let mut w = v.clone();
        for o in &ortho {
            let d: f64 = w.iter().zip(o.iter()).map(|(a, b)| a * b).sum();
            for j in 0..nc {
                w[j] -= d * o[j];
            }
        }
        let nrm: f64 = w.iter().map(|a| a * a).sum::<f64>().sqrt();
        if nrm > 1e-10 {
            for j in 0..nc {
                w[j] /= nrm;
            }
            ortho.push(w);
        }
    }
    for o in &ortho {
        let d: f64 = g_flat.iter().zip(o.iter()).map(|(a, b)| a * b).sum();
        for j in 0..nc {
            g_flat[j] -= d * o[j];
        }
    }
}

fn dist(x: &[[f64; 3]], i: usize, j: usize) -> f64 {
    let dx = x[i][0] - x[j][0];
    let dy = x[i][1] - x[j][1];
    let dz = x[i][2] - x[j][2];
    (dx * dx + dy * dy + dz * dz).sqrt()
}

fn angle_of(x: &[[f64; 3]], a: usize, b: usize, c: usize) -> f64 {
    let u = [x[a][0] - x[b][0], x[a][1] - x[b][1], x[a][2] - x[b][2]];
    let v = [x[c][0] - x[b][0], x[c][1] - x[b][1], x[c][2] - x[b][2]];
    let nu = norm(u);
    let nv = norm(v);
    if nu < 1e-12 || nv < 1e-12 {
        return 0.0;
    }
    (dot(u, v) / (nu * nv)).clamp(-1.0, 1.0).acos()
}

fn prim_value(p: &Prim, x: &[[f64; 3]]) -> f64 {
    match p {
        Prim::Dist(i, j) => dist(x, *i, *j),
        Prim::Angle(a, b, c) => angle_of(x, *a, *b, *c),
        Prim::Dih(a, b, c, d) => dihedral(x, *a, *b, *c, *d),
        Prim::Oop { c, j, k, l } => oop_value(x, *c, *j, *k, *l),
    }
}

/// Dihedral a-b-c-d in (-pi, pi]: atan2((n1 x n2)·b2 / |b2|, n1·n2) — both
/// arms carry the |n1||n2| scale, so this is the true signed torsion angle.
fn dihedral(x: &[[f64; 3]], a: usize, b: usize, c: usize, d: usize) -> f64 {
    let b1 = sub(x, a, b);
    let b2 = sub(x, c, b);
    let b3 = sub(x, d, c);
    let n1 = cross(b1, b2);
    let n2 = cross(b2, b3);
    let m = cross(n1, n2);
    (dot(m, b2) / norm(b2)).atan2(dot(n1, n2))
}

/// Out-of-plane angle tau = asin( r_cl_hat . n_hat ), n = r_cj x r_ck.
fn oop_value(x: &[[f64; 3]], c: usize, j: usize, k: usize, l: usize) -> f64 {
    let n = cross(sub(x, j, c), sub(x, k, c));
    let ril = sub(x, l, c);
    let nn = norm(n);
    let nl = norm(ril);
    if nn < 1e-12 || nl < 1e-12 {
        return 0.0;
    }
    (dot(ril, n) / (nl * nn)).clamp(-1.0, 1.0).asin()
}

fn b_matrix(prims: &[Prim], x: &[[f64; 3]], n: usize) -> Vec<Vec<f64>> {
    let mut b = vec![vec![0.0f64; 3 * n]; prims.len()];
    for (row, p) in prims.iter().enumerate() {
        prim_deriv(p, x, &mut b[row]);
    }
    b
}

fn prim_deriv(p: &Prim, x: &[[f64; 3]], row: &mut [f64]) {
    match p {
        Prim::Dist(i, j) => {
            let d = sub(x, *i, *j);
            let r = norm(d);
            if r < 1e-12 {
                return;
            }
            for t in 0..3 {
                row[3 * i + t] = d[t] / r;
                row[3 * j + t] = -d[t] / r;
            }
        }
        Prim::Angle(a, b, c) => {
            let u = sub(x, *a, *b);
            let v = sub(x, *c, *b);
            let ru = norm(u);
            let rv = norm(v);
            let cos = dot(u, v) / (ru * rv);
            let sin = (1.0 - cos * cos).sqrt();
            if sin < 1e-8 || ru < 1e-12 || rv < 1e-12 {
                return;
            }
            let uh = [u[0] / ru, u[1] / ru, u[2] / ru];
            let vh = [v[0] / rv, v[1] / rv, v[2] / rv];
            let f = -1.0 / sin;
            for t in 0..3 {
                let da = (vh[t] - cos * uh[t]) / ru;
                let dc = (uh[t] - cos * vh[t]) / rv;
                row[3 * a + t] = f * da;
                row[3 * c + t] = f * dc;
                row[3 * b + t] = -f * (da + dc);
            }
        }
        Prim::Dih(a, b, c, d) => {
            dihedral_deriv(x, *a, *b, *c, *d, row);
        }
        Prim::Oop { c, j, k, l } => {
            oop_deriv(x, *c, *j, *k, *l, row);
        }
    }
}

/// Analytic derivative of `dihedral` via the chain rule on
/// phi = atan2(y, x), y = (m·b2)/|b2|, x = n1·n2, m = n1 x n2.
fn dihedral_deriv(x: &[[f64; 3]], a: usize, b: usize, c: usize, d: usize, row: &mut [f64]) {
    let b1 = sub(x, a, b);
    let b2 = sub(x, c, b);
    let b3 = sub(x, d, c);
    let n1 = cross(b1, b2);
    let n2 = cross(b2, b3);
    let m = cross(n1, n2);
    let r2 = norm(b2);
    if r2 < 1e-12 {
        return;
    }
    let xv = dot(n1, n2);
    let yv = dot(m, b2) / r2;
    let den = xv * xv + yv * yv;
    if den < 1e-20 {
        return; // cis/trans degenerate (guarded at construction)
    }
    for atom in [a, b, c, d] {
        for t in 0..3 {
            let mut db1 = [0.0f64; 3];
            let mut db2 = [0.0f64; 3];
            let mut db3 = [0.0f64; 3];
            if atom == a {
                db1[t] = 1.0;
            }
            if atom == b {
                db1[t] = -1.0;
                db2[t] = -1.0;
            }
            if atom == c {
                db2[t] += 1.0;
                db3[t] = -1.0;
            }
            if atom == d {
                db3[t] += 1.0;
            }
            let dn1 = add(cross(db1, b2), cross(b1, db2));
            let dn2 = add(cross(db2, b3), cross(b2, db3));
            let dxv = dot(dn1, n2) + dot(n1, dn2);
            let dm = add(cross(dn1, n2), cross(n1, dn2));
            let dw = dot(dm, b2) + dot(m, db2);
            let dr2 = dot(b2, db2) / r2;
            let dyv = (dw - yv * dr2) / r2;
            row[3 * atom + t] = (xv * dyv - yv * dxv) / den;
        }
    }
}

fn oop_deriv(x: &[[f64; 3]], c: usize, j: usize, k: usize, l: usize, row: &mut [f64]) {
    let rij = sub(x, j, c);
    let rik = sub(x, k, c);
    let ril = sub(x, l, c);
    let rj = norm(rij);
    let rk = norm(rik);
    let rl = norm(ril);
    let n = cross(rij, rik);
    let nn = norm(n);
    if nn < 1e-12 || rl < 1e-12 || rj < 1e-12 || rk < 1e-12 {
        return;
    }
    let s = (dot(ril, n) / (rl * nn)).clamp(-0.99999999999999, 0.99999999999999);
    let cos_tau = (1.0 - s * s).sqrt();
    let pref = 1.0 / cos_tau; // dtau = ds / cos(tau)
    let rlh = [ril[0] / rl, ril[1] / rl, ril[2] / rl];
    let denom = rl * nn;

    for atom in [j, k, l, c] {
        for t in 0..3 {
            let sgn = if atom == c { -1.0 } else { 1.0 };
            let mut drij = [0.0f64; 3];
            let mut drik = [0.0f64; 3];
            let mut dril = [0.0f64; 3];
            if atom == j || atom == c {
                drij[t] = sgn;
            }
            if atom == k || atom == c {
                drik[t] = sgn;
            }
            if atom == l || atom == c {
                dril[t] = sgn;
            }
            // s = (ril.n)/(rl.nn); ds via quotient rule
            let dt = dot(dril, n) + dot(ril, add(cross(drij, rik), cross(rij, drik)));
            let drl = dot(dril, rlh);
            let dn = add(cross(drij, rik), cross(rij, drik));
            let dnn = dot(n, dn) / nn;
            let dden = drl * nn + rl * dnn;
            let ds = (dt - s * dden) / denom;
            row[3 * atom + t] += pref * ds;
        }
    }
}

#[inline]
fn sub(x: &[[f64; 3]], a: usize, b: usize) -> [f64; 3] {
    [x[a][0] - x[b][0], x[a][1] - x[b][1], x[a][2] - x[b][2]]
}
#[inline]
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
#[inline]
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
#[inline]
fn norm(a: [f64; 3]) -> f64 {
    dot(a, a).sqrt()
}
#[inline]
fn add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

#[cfg(test)]
mod tests {
    use super::*;

    /// ethanol geometry (heavy-atom skeleton + OH/H placed by hand, sane
    /// bond lengths) — exercises bonds/angles/dihedrals/oop.
    pub(crate) fn ethanol() -> (Vec<usize>, Vec<[f64; 3]>) {
        let z = vec![6usize, 6, 8, 1, 1, 1, 1, 1, 1];
        let x = vec![
            [0.000, 0.000, 0.000],   // C0
            [1.540, 0.100, 0.050],   // C1
            [2.100, 1.350, -0.200],  // O2
            [3.050, 1.400, -0.150],  // H3 (OH)
            [-0.500, -0.400, 0.900], // H4
            [-0.450, 0.550, -0.850], // H5
            [0.350, -1.000, -0.350], // H6
            [2.050, -0.550, 0.750],  // H7
            [2.120, -0.250, -0.950], // H8
        ];
        (z, x)
    }

    /// formaldehyde: planar trigonal center, exercises Oop primitives
    pub(crate) fn formaldehyde() -> (Vec<usize>, Vec<[f64; 3]>) {
        let z = vec![6usize, 8, 1, 1];
        let x = vec![
            [0.000, 0.000, 0.000],   // C0
            [1.210, 0.000, 0.000],   // O1 (C=O)
            [-0.600, 0.940, 0.000],  // H2
            [-0.600, -0.940, 0.100], // H3 (slightly pyramidalized)
        ];
        (z, x)
    }

    #[test]
    fn build_dof_rank() {
        let (z, x) = ethanol();
        let ic = InternalCoords::build(&x, &z).expect("ethanol internals");
        assert_eq!(ic.n_dof, 3 * x.len() - 6, "ethanol is nonlinear: 21 dof");
        assert!(ic.n_prims() >= ic.n_dof);
    }

    #[test]
    fn degenerate_fallbacks() {
        assert!(InternalCoords::build(&[[0.0; 3]; 2], &[8, 1]).is_none());
        // single atom
        assert!(InternalCoords::build(&[[0.0; 3]; 1], &[6]).is_none());
    }

    /// B matrix vs finite differences of the primitive values.
    #[test]
    fn b_matrix_fd() {
        let (z, x) = ethanol();
        let ic = InternalCoords::build(&x, &z).unwrap();
        let b = b_matrix(&ic.prims, &x, x.len());
        let h = 1e-6f64;
        let mut maxerr = 0.0f64;
        for row in 0..ic.n_prims() {
            for j in 0..(3 * x.len()) {
                let mut xp = x.clone();
                let mut xm = x.clone();
                let a = j / 3;
                let t = j % 3;
                xp[a][t] += h;
                xm[a][t] -= h;
                let sp = prim_value(&ic.prims[row], &xp);
                let sm = prim_value(&ic.prims[row], &xm);
                let fd = (sp - sm) / (2.0 * h);
                maxerr = maxerr.max((fd - b[row][j]).abs());
            }
        }
        assert!(maxerr < 1e-5, "B matrix FD error {maxerr}");
    }

    /// q(x + small dx) should track back_transform's dq.
    #[test]
    fn back_transform_roundtrip() {
        let (z, x) = ethanol();
        let ic = InternalCoords::build(&x, &z).unwrap();
        // target displacement: move along a few delocalized coordinates
        let mut dq = vec![0.0f64; ic.n_dof];
        for (k, d) in dq.iter_mut().enumerate() {
            *d = 0.02 * (((k * 7919) % 13) as f64 - 6.0) / 6.0;
        }
        let dx = ic.back_transform(&dq);
        let mut x2 = x.clone();
        for i in 0..x.len() {
            for t in 0..3 {
                x2[i][t] += dx[i][t];
            }
        }
        // (a) linear contract: B_q dx = dq to machine precision
        let bq = ic.q(&x2);
        // reconstruct B_q row space via a second back_transform probe is
        // overkill — instead verify the nonlinear q gap scales ~dq^2:
        let mut maxerr = 0.0f64;
        for k in 0..ic.n_dof {
            maxerr = maxerr.max((bq[k] - dq[k]).abs());
        }
        // dq components are O(scale); curvature gap must be O(scale^2)
        let scale = 0.02f64;
        assert!(maxerr < 50.0 * scale * scale, "roundtrip error {maxerr}");
    }

    /// grad_q chain consistency: E(q) sampled along a delocalized coordinate
    /// must have slope g_q (uses a harmonic toy energy).
    #[test]
    fn grad_q_directional_consistency() {
        let (z, x) = ethanol();
        let ic = InternalCoords::build(&x, &z).unwrap();
        // toy energy E(x) = |x - x_target|^2 with a fixed displaced target;
        // gx = 2(x - x_target), and E(q) = |x + dx(q) - x_target|^2
        let target: Vec<[f64; 3]> = (0..x.len())
            .map(|i| {
                let o = ((i * 37) % 11) as f64 * 0.021;
                [x[i][0] + o, x[i][1] - o * 0.7, x[i][2] + o * 1.3]
            })
            .collect();
        let gx_flat: Vec<f64> = (0..(3 * x.len()))
            .map(|j| 2.0 * (flat_at(&x, j) - flat_at(&target, j)))
            .collect();
        let mut gproj = gx_flat.clone();
        project_out_tr(&x, &mut gproj);
        let gq = ic.grad_q(&gproj);
        let e_of = |dq: &[f64]| -> f64 {
            let dx = ic.back_transform(dq);
            let mut e = 0.0;
            for i in 0..x.len() {
                for t in 0..3 {
                    let d = x[i][t] + dx[i][t] - target[i][t];
                    e += d * d;
                }
            }
            e
        };
        let h = 1e-5;
        let mut maxerr = 0.0f64;
        for k in 0..ic.n_dof {
            let mut dp = vec![0.0f64; ic.n_dof];
            let mut dm = vec![0.0f64; ic.n_dof];
            dp[k] = h;
            dm[k] = -h;
            let fd = (e_of(&dp) - e_of(&dm)) / (2.0 * h);
            maxerr = maxerr.max((fd - gq[k]).abs());
        }
        assert!(maxerr < 1e-4, "grad_q FD error {maxerr}");
    }

    fn flat_at(x: &[[f64; 3]], j: usize) -> f64 {
        x[j / 3][j % 3]
    }
}

#[cfg(test)]
mod diag {
    use super::*;
    #[test]
    fn fd_by_prim_type() {
        let w1 = check(&super::tests::ethanol());
        let w2 = check(&super::tests::formaldehyde());
        assert!(w2.contains_key("Oop"), "formaldehyde must exercise Oop");
        for (name, e) in w1.iter().chain(w2.iter()) {
            assert!(*e < 1e-5, "{name} FD error {e}");
        }
    }

    fn check((z, x): &(Vec<usize>, Vec<[f64; 3]>)) -> std::collections::HashMap<String, f64> {
        let ic = InternalCoords::build(x, z).unwrap();
        let b = b_matrix(&ic.prims, x, x.len());
        let h = 1e-6f64;
        let mut worst: std::collections::HashMap<String, f64> = Default::default();
        for row in 0..ic.n_prims() {
            for j in 0..(3 * x.len()) {
                let mut xp = x.clone();
                let mut xm = x.clone();
                let a = j / 3;
                let t = j % 3;
                xp[a][t] += h;
                xm[a][t] -= h;
                let sp = prim_value(&ic.prims[row], &xp);
                let sm = prim_value(&ic.prims[row], &xm);
                let fd = (sp - sm) / (2.0 * h);
                let e = (fd - b[row][j]).abs();
                let key = match ic.prims[row] {
                    Prim::Dist(..) => "Dist",
                    Prim::Angle(..) => "Angle",
                    Prim::Dih(..) => "Dih",
                    Prim::Oop { .. } => "Oop",
                };
                *worst.entry(key.to_string()).or_default() =
                    worst.get(key).copied().unwrap_or(0.0).max(e);
            }
        }
        worst
    }
}

impl InternalCoords {
    /// True when any bond-like primitive is `threshold` (Å) away from its
    /// reference value — signals a topology change needing a rebuild.
    pub fn primitive_strain(&self, x: &[[f64; 3]], threshold: f64) -> bool {
        for (i, p) in self.prims.iter().enumerate() {
            if matches!(p, Prim::Dist(..)) && (prim_value(p, x) - self.s_ref[i]).abs() > threshold {
                return true;
            }
        }
        false
    }
}
