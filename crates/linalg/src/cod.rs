//! Complete orthogonal decomposition `A·P = Q·[T 0; 0 0]·Z` (Golub & Van
//! Loan, *Matrix Computations*, 4th ed., §5.4.2 and §5.5.4; the `R = [T 0]·Z`
//! step is LAPACK's `dtzrzf`). It gives the minimum-norm least-squares
//! solution `A⁺ b` and an orthonormal basis of the null space of a
//! rank-deficient matrix.

use crate::householder::make_reflector;
use crate::matrix::Matrix;
use crate::qr::QrPivoted;
use crate::vector::{axpy, dot};

/// The complete orthogonal decomposition of an `m × n` matrix of numerical
/// rank `r`: `A·P = Q·[T 0; 0 0]·Z` with `T` upper triangular `r × r` and
/// nonsingular, `Q` and `Z` orthogonal, `P` a permutation. The trailing block
/// of the pivoted QR factor (below the rank threshold) is treated as zero.
#[derive(Clone, Debug)]
pub struct CompleteOrthogonal {
    qr: QrPivoted,
    rank: usize,
    /// `r × n`: `T` in its leading `r × r` upper triangle, the tails of the
    /// `Z` reflectors in columns `r..n`.
    t: Matrix,
    /// Coefficients of the `Z` reflectors (`H_k` acts on coordinates
    /// `{k} ∪ {r, …, n−1}`).
    tau_z: Vec<f64>,
}

impl CompleteOrthogonal {
    /// Decomposes `a`, deciding the rank with [`QrPivoted::rank`]`(rel_tol)`.
    pub fn new(a: &Matrix, rel_tol: f64) -> Self {
        let qr = QrPivoted::new(a);
        let rank = qr.rank(rel_tol);
        Self::from_qr(qr, rank)
    }

    /// Completes an existing pivoted QR factorization, treating rows
    /// `rank..` of `R` as zero. The rank is reduced, if necessary, to the
    /// number of leading nonzero diagonal entries of `R`.
    pub fn from_qr(qr: QrPivoted, rank: usize) -> Self {
        let n = qr.cols();
        let mut r = rank.min(qr.steps());
        let diag = qr.r_diagonal();
        if let Some(z) = diag[..r].iter().position(|d| *d == 0.0) {
            r = z;
        }
        let packed = qr.packed();
        let mut t = Matrix::from_fn(r, n, |i, j| if i <= j { packed[(i, j)] } else { 0.0 });
        let mut tau_z = vec![0.0; r];
        if r < n {
            let tail_len = n - r;
            let mut x = vec![0.0; 1 + tail_len];
            let mut w = vec![0.0; r];
            for k in (0..r).rev() {
                x[0] = t[(k, k)];
                for j in 0..tail_len {
                    x[1 + j] = t[(k, r + j)];
                }
                let tau = make_reflector(&mut x);
                tau_z[k] = tau;
                t[(k, k)] = x[0];
                for j in 0..tail_len {
                    t[(k, r + j)] = x[1 + j];
                }
                if tau == 0.0 || k == 0 {
                    continue;
                }
                // Apply H_k from the right to rows 0..k:
                // w = T[0..k, k] + Σ_j v_j T[0..k, r+j].
                w[..k].copy_from_slice(&t.col(k)[..k]);
                for j in 0..tail_len {
                    let vj = x[1 + j];
                    if vj != 0.0 {
                        axpy(vj, &t.col(r + j)[..k], &mut w[..k]);
                    }
                }
                axpy(-tau, &w[..k], &mut t.col_mut(k)[..k]);
                for j in 0..tail_len {
                    let vj = x[1 + j];
                    if vj != 0.0 {
                        axpy(-tau * vj, &w[..k], &mut t.col_mut(r + j)[..k]);
                    }
                }
            }
        }
        CompleteOrthogonal {
            qr,
            rank: r,
            t,
            tau_z,
        }
    }

    /// The numerical rank `r`.
    pub fn rank(&self) -> usize {
        self.rank
    }

    /// The underlying pivoted QR factorization.
    pub fn qr(&self) -> &QrPivoted {
        &self.qr
    }

    /// Applies `Zᵀ = H_{r−1} ⋯ H_0` to `w` (length `n`, in permuted order).
    fn apply_zt(&self, w: &mut [f64]) {
        let r = self.rank;
        let n = self.qr.cols();
        if r == n {
            return;
        }
        for k in 0..r {
            let tau = self.tau_z[k];
            if tau == 0.0 {
                continue;
            }
            // v = e_k + Σ_j tail_j e_{r+j}; tail_j = T[k, r+j].
            let mut s = w[k];
            for j in 0..n - r {
                s += self.t[(k, r + j)] * w[r + j];
            }
            let f = tau * s;
            w[k] -= f;
            for j in 0..n - r {
                w[r + j] -= f * self.t[(k, r + j)];
            }
        }
    }

    /// The minimum-norm least-squares solution `x = A⁺ b`: among the
    /// minimizers of `‖A x − b‖₂`, the one of least `‖x‖₂`.
    ///
    /// # Panics
    /// If `b.len()` differs from the row count.
    pub fn solve(&self, b: &[f64]) -> Vec<f64> {
        let n = self.qr.cols();
        let r = self.rank;
        let mut c = b.to_vec();
        self.qr.apply_qt(&mut c);
        // Back substitution T y = c[0..r], column oriented.
        let mut y = vec![0.0; n];
        for j in (0..r).rev() {
            let yj = c[j] / self.t[(j, j)];
            y[j] = yj;
            if yj != 0.0 {
                axpy(-yj, &self.t.col(j)[..j], &mut c[..j]);
            }
        }
        self.apply_zt(&mut y);
        let mut x = vec![0.0; n];
        for (j, &p) in self.qr.permutation().iter().enumerate() {
            x[p] = y[j];
        }
        x
    }

    /// An orthonormal basis of the null space, as the columns of an
    /// `n × (n − r)` matrix.
    pub fn null_space(&self) -> Matrix {
        let n = self.qr.cols();
        let r = self.rank;
        let perm = self.qr.permutation();
        let mut out = Matrix::zeros(n, n - r);
        let mut w = vec![0.0; n];
        for j in 0..n - r {
            w.iter_mut().for_each(|v| *v = 0.0);
            w[r + j] = 1.0;
            self.apply_zt(&mut w);
            let col = out.col_mut(j);
            for (i, &p) in perm.iter().enumerate() {
                col[p] = w[i];
            }
        }
        out
    }
}

/// The minimum-norm least-squares solution of `A x ≈ b` with the rank decided
/// by `rel_tol` (see [`QrPivoted::rank`]).
pub fn lstsq_min_norm(a: &Matrix, b: &[f64], rel_tol: f64) -> Vec<f64> {
    CompleteOrthogonal::new(a, rel_tol).solve(b)
}

/// Result of [`weighted_min_norm`].
#[derive(Clone, Debug, PartialEq)]
pub struct WeightedSolution {
    /// The solution.
    pub x: Vec<f64>,
    /// Numerical rank of `A·W^{−1/2}`.
    pub rank: usize,
}

/// Among the least-squares solutions of `A x ≈ b`, the one minimizing the
/// weighted norm `‖W^{1/2} x‖₂` with `W = diag(weights)` (all weights
/// positive): `x = W^{−1/2} (A W^{−1/2})⁺ b`. The rank of `A W^{−1/2}` is
/// decided by `rel_tol`. Heavier unknowns move less.
///
/// # Panics
/// If the dimensions disagree or a weight is not positive and finite.
pub fn weighted_min_norm(a: &Matrix, b: &[f64], weights: &[f64], rel_tol: f64) -> WeightedSolution {
    assert_eq!(weights.len(), a.cols(), "weighted_min_norm: weights length");
    assert!(
        weights.iter().all(|w| *w > 0.0 && w.is_finite()),
        "weighted_min_norm: weights must be positive and finite"
    );
    let s: Vec<f64> = weights.iter().map(|w| 1.0 / w.sqrt()).collect();
    let mut aw = a.clone();
    aw.scale_cols(&s);
    let cod = CompleteOrthogonal::new(&aw, rel_tol);
    let y = cod.solve(b);
    WeightedSolution {
        x: y.iter().zip(&s).map(|(yi, si)| yi * si).collect(),
        rank: cod.rank(),
    }
}

/// An orthonormal basis of the null space of `a` (`n × (n − rank)`), with the
/// rank decided by `rel_tol`.
pub fn null_space(a: &Matrix, rel_tol: f64) -> Matrix {
    CompleteOrthogonal::new(a, rel_tol).null_space()
}

/// The least-squares residual `b − A x` for the solution `x`.
pub fn residual(a: &Matrix, x: &[f64], b: &[f64]) -> Vec<f64> {
    let ax = a.mul_vec(x);
    b.iter().zip(&ax).map(|(bi, ai)| bi - ai).collect()
}

/// `‖x‖₂` restricted to the components of `x` in the column space of `n`,
/// i.e. `‖nᵀ x‖₂` for a matrix `n` with orthonormal columns.
pub fn projection_norm(n: &Matrix, x: &[f64]) -> f64 {
    let p: Vec<f64> = (0..n.cols()).map(|j| dot(n.col(j), x)).collect();
    crate::vector::norm2(&p)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn full_rank_square_solve() {
        let a = Matrix::from_row_major(2, 2, &[2.0, 1.0, 1.0, 3.0]);
        let x = lstsq_min_norm(&a, &[3.0, 5.0], 1e-12);
        assert!((x[0] - 0.8).abs() < 1e-14 && (x[1] - 1.4).abs() < 1e-14);
    }

    #[test]
    fn underdetermined_min_norm() {
        // x + y = 2 → minimum norm (1, 1).
        let a = Matrix::from_row_major(1, 2, &[1.0, 1.0]);
        let x = lstsq_min_norm(&a, &[2.0], 1e-12);
        assert!((x[0] - 1.0).abs() < 1e-14 && (x[1] - 1.0).abs() < 1e-14);
        let n = null_space(&a, 1e-12);
        assert_eq!((n.rows(), n.cols()), (2, 1));
        assert!((n[(0, 0)] + n[(1, 0)]).abs() < 1e-14);
        assert!((n[(0, 0)].abs() - 0.5f64.sqrt()).abs() < 1e-14);
    }

    #[test]
    fn weighted_moves_light_unknown() {
        // x + y = 2 with y 1e6 times heavier: x takes nearly everything.
        let a = Matrix::from_row_major(1, 2, &[1.0, 1.0]);
        let s = weighted_min_norm(&a, &[2.0], &[1.0, 1e6], 1e-12);
        assert!((s.x[0] + s.x[1] - 2.0).abs() < 1e-14);
        assert!(s.x[1].abs() < 1e-5 && s.rank == 1);
    }

    #[test]
    fn rank_deficient_least_squares() {
        // Two identical rows with inconsistent right-hand sides.
        let a = Matrix::from_row_major(2, 2, &[1.0, 0.0, 1.0, 0.0]);
        let x = lstsq_min_norm(&a, &[1.0, 3.0], 1e-12);
        assert!((x[0] - 2.0).abs() < 1e-14 && x[1].abs() < 1e-14);
        let r = residual(&a, &x, &[1.0, 3.0]);
        assert!((r[0] + 1.0).abs() < 1e-14 && (r[1] - 1.0).abs() < 1e-14);
        let n = null_space(&a, 1e-12);
        assert!((projection_norm(&n, &[0.0, 1.0]) - 1.0).abs() < 1e-14);
    }

    #[test]
    fn zero_matrix() {
        let a = Matrix::zeros(2, 3);
        let x = lstsq_min_norm(&a, &[1.0, 2.0], 1e-12);
        assert_eq!(x, vec![0.0; 3]);
        let n = null_space(&a, 1e-12);
        assert_eq!(n.cols(), 3);
    }
}
