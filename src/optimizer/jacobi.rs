//! Symmetric-matrix eigendecomposition via cyclic Jacobi rotations.
//!
//! Used by the internal-coordinate optimizer (G = B·Bᵀ eigendecomposition for
//! the delocalized-coordinate basis and pseudo-inverse). Hand-rolled to keep
//! the dependency set unchanged; sizes here are the primitive count (a few
//! hundred at most), where O(n^3 · sweeps) Jacobi is entirely adequate.
//! Matrices are flat row-major `&mut [f64]` (length n^2) for locality.

/// Diagonalize a symmetric matrix in place.
///
/// On entry `a` holds the symmetric matrix (only the upper triangle is
/// read). On exit the diagonal holds the eigenvalues (unordered) and `v` the
/// accumulated eigenvector matrix: column j of `v` is the eigenvector for
/// a[j*n+j], i.e. A = V·diag(a)·Vᵀ.
///
/// Convergence: off-diagonal Frobenius norm below `1e-9 × (diag norm + 1)`,
/// or 8 sweeps. The delocalized-coordinate consumer needs the basis to
/// ~1e-6, so 1e-9 buys margin at a fraction of the sweep cost (the G
/// rebuild dominates DIC overhead; full 1e-13 took 2-3x longer).
pub fn sym_jacobi(a: &mut [f64], v: &mut [f64], n: usize) {
    for i in 0..n {
        for j in 0..n {
            v[i * n + j] = if i == j { 1.0 } else { 0.0 };
        }
    }
    if n < 2 {
        return;
    }
    let mut diag_scale: f64 = 0.0;
    for i in 0..n {
        diag_scale += a[i * n + i].abs();
    }
    let tol = 1e-9 * (diag_scale / n as f64 + 1.0);

    for _sweep in 0..8 {
        let mut off = 0.0f64;
        for p in 0..n {
            for q in (p + 1)..n {
                let apq = a[p * n + q];
                off += apq * apq;
            }
        }
        if off.sqrt() < tol {
            break;
        }
        for p in 0..(n - 1) {
            for q in (p + 1)..n {
                let apq = a[p * n + q];
                if apq.abs() < 1e-300 {
                    continue;
                }
                let app = a[p * n + p];
                let aqq = a[q * n + q];
                // rotation angle: tan(2φ) = 2 apq / (aqq - app); pick the
                // smaller-|t| root for stability
                let theta = (aqq - app) / (2.0 * apq);
                let t = if theta >= 0.0 {
                    1.0 / (theta + (1.0 + theta * theta).sqrt())
                } else {
                    1.0 / (theta - (1.0 + theta * theta).sqrt())
                };
                let c = 1.0 / (1.0 + t * t).sqrt();
                let s = t * c;
                // update rows/cols p and q of a (symmetric)
                for k in 0..n {
                    if k != p && k != q {
                        let akp = a[k * n + p];
                        let akq = a[k * n + q];
                        let akp_new = c * akp - s * akq;
                        let akq_new = s * akp + c * akq;
                        a[k * n + p] = akp_new;
                        a[p * n + k] = akp_new;
                        a[k * n + q] = akq_new;
                        a[q * n + k] = akq_new;
                    }
                }
                let app_new = c * c * app - 2.0 * s * c * apq + s * s * aqq;
                let aqq_new = s * s * app + 2.0 * s * c * apq + c * c * aqq;
                a[p * n + p] = app_new;
                a[q * n + q] = aqq_new;
                a[p * n + q] = 0.0;
                a[q * n + p] = 0.0;
                // accumulate eigenvectors: columns p, q of v
                for k in 0..n {
                    let vkp = v[k * n + p];
                    let vkq = v[k * n + q];
                    v[k * n + p] = c * vkp - s * vkq;
                    v[k * n + q] = s * vkp + c * vkq;
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn mat(rows: &[&[f64]]) -> Vec<f64> {
        let n = rows.len();
        let mut m = vec![0.0f64; n * n];
        for (i, r) in rows.iter().enumerate() {
            for (j, v) in r.iter().enumerate() {
                m[i * n + j] = *v;
            }
        }
        m
    }

    #[test]
    fn diagonal_matrix_unchanged() {
        let n = 3;
        let mut a = mat(&[&[2.0, 0.0, 0.0], &[0.0, 5.0, 0.0], &[0.0, 0.0, 9.0]]);
        let mut v = vec![0.0f64; n * n];
        sym_jacobi(&mut a, &mut v, n);
        let mut vals: Vec<f64> = (0..n).map(|i| a[i * n + i]).collect();
        vals.sort_by(|p, q| p.partial_cmp(q).unwrap());
        assert!((vals[0] - 2.0).abs() < 1e-12);
        assert!((vals[1] - 5.0).abs() < 1e-12);
        assert!((vals[2] - 9.0).abs() < 1e-12);
    }

    #[test]
    fn known_2x2() {
        // [[2, 1], [1, 2]] -> eigenvalues 1 and 3, eigenvectors (1,-1)/√2, (1,1)/√2
        let n = 2;
        let mut a = mat(&[&[2.0, 1.0], &[1.0, 2.0]]);
        let mut v = vec![0.0f64; n * n];
        sym_jacobi(&mut a, &mut v, n);
        let mut vals = [a[0], a[3]];
        vals.sort_by(|p, q| p.partial_cmp(q).unwrap());
        assert!((vals[0] - 1.0).abs() < 1e-9);
        assert!((vals[1] - 3.0).abs() < 1e-9);
        // residual A·v = λ·v
        for j in 0..2 {
            for i in 0..2 {
                let avi = 2.0 * v[i * n + j] + if i == 0 { v[n + j] } else { v[j] };
                assert!((avi - a[j * n + j] * v[i * n + j]).abs() < 1e-7);
            }
        }
    }

    #[test]
    fn residual_and_orthogonality() {
        // deterministic pseudo-random symmetric 12x12
        let n = 12;
        let seed = |i: usize, j: usize| -> f64 {
            let x = (i * 73 + j * 151 + 13) as f64;
            ((x * 0.1013).fract() - 0.5) * 2.0
        };
        let mut a0 = vec![0.0f64; n * n];
        for i in 0..n {
            for j in 0..n {
                a0[i * n + j] = if i <= j { seed(i, j) } else { a0[j * n + i] };
            }
        }
        for i in 0..n {
            a0[i * n + i] += 3.0;
        }
        let mut a = a0.clone();
        let mut v = vec![0.0f64; n * n];
        sym_jacobi(&mut a, &mut v, n);
        // A0 = V diag(A) V^T -> residual per entry
        for i in 0..n {
            for j in 0..n {
                let mut s = 0.0;
                for k in 0..n {
                    s += v[i * n + k] * a[k * n + k] * v[j * n + k];
                }
                assert!(
                    (s - a0[i * n + j]).abs() < 1e-7,
                    "({i},{j}): {s} vs {}",
                    a0[i * n + j]
                );
            }
        }
        // orthogonality V^T V = I
        for i in 0..n {
            for j in 0..n {
                let mut s = 0.0;
                for k in 0..n {
                    s += v[k * n + i] * v[k * n + j];
                }
                let expect = if i == j { 1.0 } else { 0.0 };
                assert!((s - expect).abs() < 1e-8, "V^T V ({i},{j}) = {s}");
            }
        }
    }
}
