//! One-sided Jacobi singular value decomposition (Hestenes 1958; Golub & Van
//! Loan, *Matrix Computations*, 4th ed., §8.6.4; Demmel & Veselić,
//! "Jacobi's method is more accurate than QR", *SIAM J. Matrix Anal. Appl.*
//! 13(4), 1992). Tall matrices are first reduced by pivoted QR, which also
//! speeds up convergence (Drmač & Veselić, "New fast and accurate Jacobi SVD
//! algorithm", *SIAM J. Matrix Anal. Appl.* 29(4), 2008).

use crate::matrix::Matrix;
use crate::qr::QrPivoted;
use crate::vector::{dot, norm2};

/// Sweep limit; cyclic one-sided Jacobi converges quadratically and needs
/// far fewer sweeps in practice.
const MAX_SWEEPS: usize = 80;

/// `A = U·diag(s)·Vᵀ` for an `m × n` matrix `A`.
#[derive(Clone, Debug)]
pub struct Svd {
    /// `m × n`: the left singular vector for each singular value (a zero
    /// column where the singular value is zero).
    pub u: Matrix,
    /// The `n` singular values, sorted in decreasing order (ties keep the
    /// order of the input columns). When `m < n`, at least `n − m` are zero
    /// up to rounding.
    pub s: Vec<f64>,
    /// `n × n` orthogonal: the right singular vectors, one per column.
    pub v: Matrix,
}

impl Svd {
    /// Numerical rank: the number of singular values `> rel_tol · s[0]`.
    pub fn rank(&self, rel_tol: f64) -> usize {
        match self.s.first() {
            Some(&s0) if s0 > 0.0 && s0.is_finite() => {
                self.s.iter().take_while(|&&x| x > rel_tol * s0).count()
            }
            _ => 0,
        }
    }

    /// An orthonormal basis of the null space (`n × (n − rank)`): the right
    /// singular vectors of the singular values at or below the threshold.
    pub fn null_space(&self, rel_tol: f64) -> Matrix {
        let r = self.rank(rel_tol);
        let cols: Vec<usize> = (r..self.v.cols()).collect();
        self.v.select_cols(&cols)
    }
}

/// Orthogonalizes the columns of `b` by plane rotations, accumulating the
/// rotations into `v`: on return `b_in · v = b_out` with mutually orthogonal
/// columns.
fn jacobi_sweeps(b: &mut Matrix, v: &mut Matrix) {
    let n = b.cols();
    let m = b.rows();
    let tol = f64::EPSILON * (m.max(1) as f64).sqrt();
    // Columns shorter than this are zero to working precision: rotating
    // them changes nothing above rounding level, and without the cut-off a
    // column that must vanish (more columns than rows, or rank deficiency)
    // would only shrink geometrically, sweep after sweep.
    let negligible = f64::EPSILON * b.frobenius_norm() * 1e-2;
    for _ in 0..MAX_SWEEPS {
        let mut rotated = false;
        for i in 0..n {
            for j in i + 1..n {
                let (bi, bj) = b.two_cols_mut(i, j);
                let ni = norm2(bi);
                let nj = norm2(bj);
                if ni <= negligible || nj <= negligible {
                    continue;
                }
                let gamma = dot(bi, bj);
                if gamma.abs() <= tol * ni * nj {
                    continue;
                }
                rotated = true;
                // Rotation (bi, bj) ← (c bi − s bj, s bi + c bj) making them
                // orthogonal; t is the smaller root of t² + 2ζt − 1 = 0.
                let zeta = ((nj - ni) * (nj + ni)) / (2.0 * gamma);
                let sign = if zeta >= 0.0 { 1.0 } else { -1.0 };
                let t = sign / (zeta.abs() + norm2(&[1.0, zeta]));
                let c = 1.0 / (1.0 + t * t).sqrt();
                let s = c * t;
                rotate(bi, bj, c, s);
                let (vi, vj) = v.two_cols_mut(i, j);
                rotate(vi, vj, c, s);
            }
        }
        if !rotated {
            break;
        }
    }
}

fn rotate(x: &mut [f64], y: &mut [f64], c: f64, s: f64) {
    for (a, b) in x.iter_mut().zip(y.iter_mut()) {
        let (p, q) = (*a, *b);
        *a = c * p - s * q;
        *b = s * p + c * q;
    }
}

/// The singular value decomposition of `a` by one-sided Jacobi.
pub fn svd(a: &Matrix) -> Svd {
    let (m, n) = (a.rows(), a.cols());
    let (b, v, u_left): (Matrix, Matrix, Option<Matrix>) = if m > n {
        // A·P = Q·R; Jacobi on R: R·V_R = B ⇒ A·(P·V_R) = Q_thin·B.
        let qr = QrPivoted::new(a);
        let mut b = qr.r();
        let mut vr = Matrix::identity(n);
        jacobi_sweeps(&mut b, &mut vr);
        let perm = qr.permutation();
        let mut v = Matrix::zeros(n, n);
        for j in 0..n {
            for (i, &p) in perm.iter().enumerate() {
                v[(p, j)] = vr[(i, j)];
            }
        }
        (b, v, Some(qr.thin_q()))
    } else {
        let mut b = a.clone();
        let mut v = Matrix::identity(n);
        jacobi_sweeps(&mut b, &mut v);
        (b, v, None)
    };
    let norms: Vec<f64> = (0..n).map(|j| norm2(b.col(j))).collect();
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&x, &y| {
        norms[y]
            .partial_cmp(&norms[x])
            .unwrap_or(std::cmp::Ordering::Equal)
            .then(x.cmp(&y))
    });
    let s: Vec<f64> = order.iter().map(|&j| norms[j]).collect();
    let v_sorted = v.select_cols(&order);
    // Normalized columns of B, in sorted order.
    let mut ub = Matrix::zeros(b.rows(), n);
    for (k, &j) in order.iter().enumerate() {
        let sj = norms[j];
        if sj > 0.0 {
            for (dst, src) in ub.col_mut(k).iter_mut().zip(b.col(j)) {
                *dst = src / sj;
            }
        }
    }
    let u = match u_left {
        Some(q) => q.matmul(&ub),
        None => ub,
    };
    debug_assert_eq!(u.rows(), m);
    Svd { u, s, v: v_sorted }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn diagonal_matrix() {
        let a = Matrix::from_row_major(3, 3, &[1.0, 0.0, 0.0, 0.0, -3.0, 0.0, 0.0, 0.0, 2.0]);
        let d = svd(&a);
        assert_eq!(d.s, vec![3.0, 2.0, 1.0]);
        assert_eq!(d.rank(1e-10), 3);
        assert_eq!(d.null_space(1e-10).cols(), 0);
    }

    #[test]
    fn reconstructs_small_shapes() {
        for (m, n) in [(3, 2), (2, 3), (4, 4), (1, 3), (3, 1)] {
            let a = Matrix::from_fn(m, n, |i, j| ((i * 5 + j * 3) % 7) as f64 - 3.0);
            let d = svd(&a);
            let mut usv = d.u.clone();
            usv.scale_cols(&d.s);
            let rec = usv.matmul(&d.v.transpose());
            for i in 0..m {
                for j in 0..n {
                    assert!((rec[(i, j)] - a[(i, j)]).abs() < 1e-13, "{m}x{n} {rec:?} {a:?} {d:?}");
                }
            }
            for w in d.s.windows(2) {
                assert!(w[0] >= w[1]);
            }
        }
    }

    #[test]
    fn zero_and_empty() {
        let d = svd(&Matrix::zeros(2, 2));
        assert_eq!(d.rank(1e-10), 0);
        assert_eq!(d.null_space(1e-10).cols(), 2);
        let e = svd(&Matrix::zeros(0, 2));
        assert_eq!(e.s, vec![0.0, 0.0]);
    }
}
