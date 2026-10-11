//! Householder QR with column pivoting (Businger & Golub 1965; Golub & Van
//! Loan, *Matrix Computations*, 4th ed., Algorithm 5.4.1), with the column
//! norm downdating and recomputation safeguard of LAPACK `dlaqp2`
//! (Drmač & Bujanović, "On the failure of rank-revealing QR factorization
//! software", *ACM TOMS* 35(2), 2008).

use crate::householder::{apply_reflector, make_reflector};
use crate::matrix::Matrix;
use crate::vector::norm2;

/// `A·P = Q·R` for an `m × n` matrix `A`: `Q` is orthogonal (`m × m`), `R`
/// is upper trapezoidal (`min(m, n) × n`) with `|R₀₀| ≥ |R₁₁| ≥ …`
/// (up to the norm downdating), and `P` is a permutation chosen greedily so
/// that each step eliminates the remaining column of largest norm (ties go
/// to the lowest index). `Q` is stored as Householder reflectors.
#[derive(Clone, Debug)]
pub struct QrPivoted {
    /// `R` in and above the diagonal, reflector tails below it.
    qr: Matrix,
    /// Reflector coefficients, one per step.
    tau: Vec<f64>,
    /// Column `j` of `A·P` is column `perm[j]` of `A`.
    perm: Vec<usize>,
}

impl QrPivoted {
    /// Factorizes `a`.
    pub fn new(a: &Matrix) -> Self {
        let (m, n) = (a.rows(), a.cols());
        let k = m.min(n);
        let mut qr = a.clone();
        let mut perm: Vec<usize> = (0..n).collect();
        let mut tau = vec![0.0; k];
        let mut vn1: Vec<f64> = (0..n).map(|j| norm2(qr.col(j))).collect();
        let mut vn2 = vn1.clone();
        let tol3z = f64::EPSILON.sqrt();
        for i in 0..k {
            // Pivot: first column of largest (downdated) norm.
            let mut p = i;
            for j in i + 1..n {
                if vn1[j] > vn1[p] {
                    p = j;
                }
            }
            if p != i {
                qr.swap_cols(i, p);
                perm.swap(i, p);
                vn1.swap(i, p);
                vn2.swap(i, p);
            }
            tau[i] = make_reflector(&mut qr.col_mut(i)[i..]);
            for j in i + 1..n {
                let (ci, cj) = qr.two_cols_mut(i, j);
                apply_reflector(tau[i], &ci[i + 1..], &mut cj[i..]);
            }
            // Downdate the remaining column norms (LAPACK dlaqp2).
            for j in i + 1..n {
                if vn1[j] != 0.0 {
                    let r = qr[(i, j)].abs() / vn1[j];
                    let temp = (1.0 - r * r).max(0.0);
                    let ratio = vn1[j] / vn2[j];
                    let temp2 = temp * ratio * ratio;
                    if temp2 <= tol3z {
                        let fresh = if i + 1 < m {
                            norm2(&qr.col(j)[i + 1..])
                        } else {
                            0.0
                        };
                        vn1[j] = fresh;
                        vn2[j] = fresh;
                    } else {
                        vn1[j] *= temp.sqrt();
                    }
                }
            }
        }
        QrPivoted { qr, tau, perm }
    }

    /// Rows of the factorized matrix.
    pub fn rows(&self) -> usize {
        self.qr.rows()
    }

    /// Columns of the factorized matrix.
    pub fn cols(&self) -> usize {
        self.qr.cols()
    }

    /// Number of Householder steps, `min(m, n)`.
    pub fn steps(&self) -> usize {
        self.tau.len()
    }

    /// The column permutation: column `j` of `A·P` is column
    /// `permutation()[j]` of `A`.
    pub fn permutation(&self) -> &[usize] {
        &self.perm
    }

    /// Entry `(i, j)` of `R` (zero below the diagonal).
    pub fn r_entry(&self, i: usize, j: usize) -> f64 {
        if i > j || i >= self.steps() {
            0.0
        } else {
            self.qr[(i, j)]
        }
    }

    /// The diagonal of `R`.
    pub fn r_diagonal(&self) -> Vec<f64> {
        (0..self.steps()).map(|i| self.qr[(i, i)]).collect()
    }

    /// `R` as a `min(m, n) × n` upper-trapezoidal matrix.
    pub fn r(&self) -> Matrix {
        Matrix::from_fn(self.steps(), self.cols(), |i, j| self.r_entry(i, j))
    }

    /// The full orthogonal factor `Q` (`m × m`).
    pub fn q(&self) -> Matrix {
        let m = self.rows();
        let mut q = Matrix::identity(m);
        for j in 0..m {
            self.apply_q(q.col_mut(j));
        }
        q
    }

    /// The first `min(m, n)` columns of `Q`.
    pub fn thin_q(&self) -> Matrix {
        let m = self.rows();
        let k = self.steps();
        let mut q = Matrix::zeros(m, k);
        for j in 0..k {
            let c = q.col_mut(j);
            c[j] = 1.0;
            self.apply_q(c);
        }
        q
    }

    /// `b ← Qᵀ b` (`b.len() == rows()`).
    pub fn apply_qt(&self, b: &mut [f64]) {
        assert_eq!(b.len(), self.rows(), "apply_qt: dimension mismatch");
        for i in 0..self.steps() {
            apply_reflector(self.tau[i], &self.qr.col(i)[i + 1..], &mut b[i..]);
        }
    }

    /// `b ← Q b` (`b.len() == rows()`).
    pub fn apply_q(&self, b: &mut [f64]) {
        assert_eq!(b.len(), self.rows(), "apply_q: dimension mismatch");
        for i in (0..self.steps()).rev() {
            apply_reflector(self.tau[i], &self.qr.col(i)[i + 1..], &mut b[i..]);
        }
    }

    /// Numerical rank: the number of leading diagonal entries of `R` with
    /// `|Rᵢᵢ| > rel_tol · |R₀₀|` (counting stops at the first that fails).
    /// A zero matrix has rank 0.
    pub fn rank(&self, rel_tol: f64) -> usize {
        let k = self.steps();
        if k == 0 {
            return 0;
        }
        let r00 = self.qr[(0, 0)].abs();
        if r00 == 0.0 || !r00.is_finite() {
            return 0;
        }
        let thr = rel_tol * r00;
        (0..k).take_while(|&i| self.qr[(i, i)].abs() > thr).count()
    }

    /// The compact factor storage (`R` on and above the diagonal, reflector
    /// tails below), for the complete orthogonal decomposition.
    pub(crate) fn packed(&self) -> &Matrix {
        &self.qr
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn check_factorization(a: &Matrix) {
        let qr = QrPivoted::new(a);
        let q = qr.q();
        let r = qr.r();
        // Pad R to m × n for the product.
        let rp = Matrix::from_fn(a.rows(), a.cols(), |i, j| {
            if i < r.rows() {
                r[(i, j)]
            } else {
                0.0
            }
        });
        let ap = a.select_cols(qr.permutation());
        let qrp = q.matmul(&rp);
        for i in 0..a.rows() {
            for j in 0..a.cols() {
                assert!((ap[(i, j)] - qrp[(i, j)]).abs() < 1e-13, "{i} {j}");
            }
        }
        let qtq = q.tr_matmul(&q);
        for i in 0..q.rows() {
            for j in 0..q.cols() {
                let e = if i == j { 1.0 } else { 0.0 };
                assert!((qtq[(i, j)] - e).abs() < 1e-14);
            }
        }
    }

    #[test]
    fn small_shapes() {
        check_factorization(&Matrix::from_row_major(
            3,
            2,
            &[1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
        ));
        check_factorization(&Matrix::from_row_major(
            2,
            3,
            &[1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
        ));
        check_factorization(&Matrix::zeros(3, 3));
        check_factorization(&Matrix::zeros(0, 3));
        check_factorization(&Matrix::zeros(3, 0));
    }

    #[test]
    fn rank_of_rank_one() {
        let a = Matrix::from_fn(4, 3, |i, j| ((i + 1) * (j + 2)) as f64);
        let qr = QrPivoted::new(&a);
        assert_eq!(qr.rank(1e-10), 1);
        assert_eq!(QrPivoted::new(&Matrix::zeros(2, 2)).rank(1e-10), 0);
        assert_eq!(QrPivoted::new(&Matrix::identity(3)).rank(1e-10), 3);
    }

    #[test]
    fn apply_q_inverts_apply_qt() {
        let a = Matrix::from_fn(5, 3, |i, j| ((i * 7 + j * 3) % 5) as f64 - 2.0);
        let qr = QrPivoted::new(&a);
        let b0 = vec![1.0, -2.0, 0.5, 3.0, 4.0];
        let mut b = b0.clone();
        qr.apply_qt(&mut b);
        qr.apply_q(&mut b);
        for (x, y) in b.iter().zip(&b0) {
            assert!((x - y).abs() < 1e-14);
        }
        let tq = qr.thin_q();
        let fq = qr.q();
        for i in 0..5 {
            for j in 0..3 {
                assert_eq!(tq[(i, j)], fq[(i, j)]);
            }
        }
    }
}
