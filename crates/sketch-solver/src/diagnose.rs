//! Diagnosis at the solution (design note 0005 §5 step 4): numerical rank,
//! degrees of freedom, free unknowns, redundant and conflicting equations.

use openshape_linalg::{
    lstsq_min_norm, norm2, select_columns_in_order, svd, CompleteOrthogonal, Matrix, QrPivoted,
};

use crate::newton::{CompSystem, RANK_TOL};

/// Unknowns whose row of an orthonormal null-space basis is longer than
/// this are free (design note §5).
pub(crate) const FREE_TOL: f64 = 1e-8;
/// Least-squares residuals above this (scaled) mark conflicting equations.
pub(crate) const CONFLICT_TOL: f64 = 1e-9;

/// Diagnosis of one component.
#[derive(Clone, Debug)]
pub(crate) struct CompDiag {
    /// Numerical rank of the Jacobian.
    pub rank: usize,
    /// Degrees of freedom of the visible unknowns.
    pub dof: usize,
    /// Per local unknown: free.
    pub free: Vec<bool>,
    /// Local indices of redundant equations.
    pub redundant: Vec<usize>,
    /// An orthonormal null-space basis, one row per local unknown (the
    /// identity when nothing could be said: everything free).
    pub null: Matrix,
}

/// The numerical rank and an orthonormal null-space basis of `j`: rank from
/// pivoted QR (`|R_kk| ≤ 1e-10·|R_00|` is zero), confirmed by rigorous
/// singular-value bounds (`σ_{r+1} ≤ ‖R₂₂‖_F`, `σ_r ≥ 1/‖R₁₁⁻¹‖_F`; Golub &
/// Van Loan §5.4.1) that prove the SVD criterion `σ_k > 1e-10·σ_1` gives
/// the same rank. Only when the bounds cannot decide is the Jacobi SVD
/// computed, and then its rank and null space are used. (The design note
/// asked for the SVD on every component of up to 200 unknowns; the bounds
/// give the same guarantee at a fraction of the cost: the SVD takes about
/// 0.35 s at 200 unknowns.)
pub(crate) fn rank_and_null_space(j: &Matrix) -> (usize, Matrix) {
    let qr = QrPivoted::new(j);
    let r_qr = qr.rank(RANK_TOL);
    let use_svd = !bounds_confirm(&qr, r_qr, j);
    if use_svd {
        let d = svd(j);
        let r = d.rank(RANK_TOL);
        (r, d.null_space(RANK_TOL))
    } else {
        let cod = CompleteOrthogonal::from_qr(qr, r_qr);
        let r = cod.rank();
        (r, cod.null_space())
    }
}

/// True if the QR rank `r` provably agrees with the SVD criterion
/// `σ_k > 1e-10·σ_1`.
fn bounds_confirm(qr: &QrPivoted, r: usize, j: &Matrix) -> bool {
    let k = qr.steps();
    let n = qr.cols();
    let r00 = qr.r_entry(0, 0).abs();
    // σ_{r+1} ≤ ‖R₂₂‖_F and σ_1 ≥ |R₀₀|.
    let mut r22 = Vec::new();
    for c in r..n {
        for i in r..k.min(c + 1) {
            r22.push(qr.r_entry(i, c));
        }
    }
    if norm2(&r22) > RANK_TOL * r00 {
        return false;
    }
    if r == 0 {
        return true;
    }
    // σ_r ≥ σ_min(R₁₁) ≥ 1/‖R₁₁⁻¹‖_F and σ_1 ≤ ‖J‖_F.
    let mut inv_sq = 0.0f64;
    let mut x = vec![0.0; r];
    for col in 0..r {
        // Solve R₁₁ x = e_col by back substitution.
        x.iter_mut().for_each(|v| *v = 0.0);
        x[col] = 1.0;
        for i in (0..=col).rev() {
            let mut s = x[i];
            for (jj, xj) in x.iter().enumerate().take(col + 1).skip(i + 1) {
                s -= qr.r_entry(i, jj) * xj;
            }
            x[i] = s / qr.r_entry(i, i);
        }
        let nx = norm2(&x);
        inv_sq += nx * nx;
    }
    let lower = 1.0 / inv_sq.sqrt();
    lower.is_finite() && lower > RANK_TOL * j.frobenius_norm()
}

/// Diagnoses a component at `u` (a point where it is admissible).
/// `hidden[k]` marks hidden curve parameters.
pub(crate) fn diagnose(sys: &mut CompSystem<'_>, u: &[f64], hidden: &[bool]) -> CompDiag {
    let n = sys.n();
    let Some((_, j, _)) = sys.jacobian(u) else {
        // Not admissible: nothing can be said; report everything free.
        return CompDiag {
            rank: 0,
            dof: hidden.iter().filter(|h| !**h).count(),
            free: vec![true; n],
            redundant: Vec::new(),
            null: Matrix::identity(n),
        };
    };
    let m = j.rows();
    let (rank, nsp) = rank_and_null_space(&j);
    let free: Vec<bool> = (0..n).map(|k| norm2(&nsp.row(k)) > FREE_TOL).collect();
    let dof = if hidden.iter().any(|h| *h) {
        // Dimension of the null space seen through the visible unknowns.
        let vis: Vec<usize> = (0..n).filter(|&k| !hidden[k]).collect();
        if nsp.cols() == 0 || vis.is_empty() {
            0
        } else {
            let nv = nsp.select_rows(&vis);
            svd(&nv).s.iter().filter(|s| **s > FREE_TOL).count()
        }
    } else {
        n - rank
    };
    let redundant = if m > rank {
        select_columns_in_order(&j.transpose(), RANK_TOL, Some(rank)).rejected
    } else {
        Vec::new()
    };
    CompDiag {
        rank,
        dof,
        free,
        redundant,
        null: nsp,
    }
}

/// Local indices of the equations that carry a least-squares residual above
/// [`CONFLICT_TOL`] at `u`, plus angle equations on the antipodal branch.
pub(crate) fn conflicts(sys: &mut CompSystem<'_>, u: &[f64]) -> Vec<usize> {
    let mut out = sys.wrong_branch(u);
    if let Some((r, j, _)) = sys.jacobian(u) {
        let neg: Vec<f64> = r.iter().map(|x| -x).collect();
        let d = lstsq_min_norm(&j, &neg, RANK_TOL);
        let jd = j.mul_vec(&d);
        for (i, (ri, ji)) in r.iter().zip(&jd).enumerate() {
            let rs = ri + ji;
            if rs.abs() > CONFLICT_TOL || !rs.is_finite() {
                out.push(i);
            }
        }
    }
    out.sort_unstable();
    out.dedup();
    out
}
