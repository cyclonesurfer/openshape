//! Property tests for the dense kernels: reconstruction identities,
//! orthogonality, least-squares and minimum-norm optimality, rank decisions
//! and bit-for-bit determinism, on random, tall, wide, scaled and
//! rank-deficient matrices.

use openshape_linalg::{
    lstsq_min_norm, null_space, select_columns_in_order, svd, weighted_min_norm, Matrix,
    QrPivoted,
};
use proptest::prelude::*;

const RANK_TOL: f64 = 1e-10;

/// Matrix shapes: small square-ish, tall and wide, including empty.
fn dims() -> impl Strategy<Value = (usize, usize)> {
    prop_oneof![
        (0usize..8, 0usize..8),
        (10usize..26, 1usize..8),
        (1usize..8, 10usize..26),
    ]
}

/// A random matrix with entries in [-1, 1] times 10^e.
fn dense() -> impl Strategy<Value = Matrix> {
    (dims(), -8i32..=8).prop_flat_map(|((m, n), e)| {
        prop::collection::vec(-1.0f64..1.0, m * n).prop_map(move |v| {
            let s = 10f64.powi(e);
            Matrix::from_col_major(m, n, v.into_iter().map(|x| x * s).collect())
        })
    })
}

/// A matrix of known rank r = rank(L·R) with random L (m×r) and R (r×n),
/// scaled by 10^e, optionally with a duplicated and a zeroed column.
fn low_rank() -> impl Strategy<Value = (Matrix, usize)> {
    (dims(), -6i32..=6, any::<bool>())
        .prop_flat_map(|((m, n), e, structured)| {
            let kmax = m.min(n);
            (Just((m, n, e, structured)), 0..=kmax)
        })
        .prop_flat_map(|((m, n, e, structured), r)| {
            (
                Just((m, n, e, structured, r)),
                prop::collection::vec(-1.0f64..1.0, m * r),
                prop::collection::vec(-1.0f64..1.0, r * n),
            )
        })
        .prop_map(|((m, n, e, structured, r), l, rr)| {
            let l = Matrix::from_col_major(m, r, l);
            let rr = Matrix::from_col_major(r, n, rr);
            let mut a = l.matmul(&rr);
            let s = 10f64.powi(e);
            for j in 0..n {
                for v in a.col_mut(j) {
                    *v *= s;
                }
            }
            if structured && n >= 3 {
                let c0 = a.col(0).to_vec();
                a.col_mut(n - 1).copy_from_slice(&c0);
                a.col_mut(1).iter_mut().for_each(|v| *v = 0.0);
            }
            // Rank of the structured variant can drop below r; recompute it
            // from an exact property: the product rank is r unless the
            // factors are degenerate, which uniform random factors are not
            // (with probability one). Duplicating a column keeps the column
            // space; zeroing one can only lower it.
            (a, r)
        })
}

fn frob_diff(a: &Matrix, b: &Matrix) -> f64 {
    assert_eq!((a.rows(), a.cols()), (b.rows(), b.cols()));
    let d = Matrix::from_fn(a.rows(), a.cols(), |i, j| a[(i, j)] - b[(i, j)]);
    d.frobenius_norm()
}

fn orthonormality_error(q: &Matrix) -> f64 {
    let g = q.tr_matmul(q);
    let mut e = 0.0f64;
    for i in 0..g.rows() {
        for j in 0..g.cols() {
            let target = if i == j { 1.0 } else { 0.0 };
            e = e.max((g[(i, j)] - target).abs());
        }
    }
    e
}

fn bits(m: &Matrix) -> Vec<u64> {
    m.as_col_major().iter().map(|v| v.to_bits()).collect()
}

/// The rank of the structured low-rank matrices is at most r; the exact
/// value is the SVD rank of the unscaled factors, so compare QR and SVD
/// against each other and against the upper bound r.
fn check_rank_agreement(a: &Matrix, r: usize) -> Result<(), TestCaseError> {
    let qr_rank = QrPivoted::new(a).rank(RANK_TOL);
    let svd_rank = svd(a).rank(RANK_TOL);
    prop_assert_eq!(qr_rank, svd_rank, "QR and SVD disagree on the rank");
    prop_assert!(qr_rank <= r);
    Ok(())
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(256))]

    #[test]
    fn qr_reconstructs_and_is_orthogonal(a in dense()) {
        let qr = QrPivoted::new(&a);
        let (m, n) = (a.rows(), a.cols());
        let q = qr.q();
        prop_assert!(orthonormality_error(&q) <= 1e-14 * (m.max(1) as f64));
        let r = qr.r();
        let rp = Matrix::from_fn(m, n, |i, j| if i < r.rows() { r[(i, j)] } else { 0.0 });
        for i in 0..r.rows() {
            for j in 0..i.min(n) {
                prop_assert_eq!(r[(i, j)], 0.0);
            }
        }
        let ap = a.select_cols(qr.permutation());
        let err = frob_diff(&ap, &q.matmul(&rp));
        prop_assert!(err <= 1e-12 * a.frobenius_norm(), "‖AP − QR‖ = {err:e}");
        // The permutation is a permutation.
        let mut p = qr.permutation().to_vec();
        p.sort_unstable();
        prop_assert_eq!(p, (0..n).collect::<Vec<_>>());
        // Pivoting makes |R_kk| non-increasing (up to the norm downdating).
        let d = qr.r_diagonal();
        for w in d.windows(2) {
            prop_assert!(w[1].abs() <= w[0].abs() * (1.0 + 1e-8) + 1e-300);
        }
    }

    #[test]
    fn rank_decisions_agree((a, r) in low_rank()) {
        check_rank_agreement(&a, r)?;
    }

    #[test]
    fn unstructured_low_rank_is_exact((a, r) in low_rank()) {
        // Without the structured modification (detected by a zero column 1
        // in a matrix with n ≥ 3 that is not all zero), the rank is exactly r.
        let n = a.cols();
        let structured = n >= 3 && a.col(1).iter().all(|v| *v == 0.0) && r > 0;
        if !structured {
            prop_assert_eq!(QrPivoted::new(&a).rank(RANK_TOL), r);
            prop_assert_eq!(svd(&a).rank(RANK_TOL), r);
        }
    }

    #[test]
    fn svd_reconstructs((a, _r) in low_rank()) {
        let d = svd(&a);
        let (m, n) = (a.rows(), a.cols());
        prop_assert_eq!(d.s.len(), n);
        prop_assert!(orthonormality_error(&d.v) <= 1e-13 * (n.max(1) as f64));
        for w in d.s.windows(2) {
            prop_assert!(w[0] >= w[1]);
        }
        prop_assert!(d.s.iter().all(|s| *s >= 0.0 && s.is_finite()));
        let mut us = d.u.clone();
        us.scale_cols(&d.s);
        let err = frob_diff(&a, &us.matmul(&d.v.transpose()));
        prop_assert!(err <= 1e-12 * a.frobenius_norm().max(f64::MIN_POSITIVE), "‖A − UΣVᵀ‖ = {err:e}");
        // Left singular vectors of the nonzero singular values are orthonormal.
        let k = d.rank(RANK_TOL);
        let cols: Vec<usize> = (0..k).collect();
        prop_assert!(orthonormality_error(&d.u.select_cols(&cols)) <= 1e-12 * (m.max(1) as f64));
    }

    #[test]
    fn svd_of_dense(a in dense()) {
        let d = svd(&a);
        let mut us = d.u.clone();
        us.scale_cols(&d.s);
        let err = frob_diff(&a, &us.matmul(&d.v.transpose()));
        prop_assert!(err <= 1e-12 * a.frobenius_norm().max(f64::MIN_POSITIVE));
        prop_assert!(orthonormality_error(&d.v) <= 1e-13 * (a.cols().max(1) as f64));
        // Largest singular value equals the 2-norm: it bounds every column norm.
        for j in 0..a.cols() {
            let c = openshape_linalg::norm2(a.col(j));
            prop_assert!(c <= d.s[0] * (1.0 + 1e-12));
        }
    }

    #[test]
    fn min_norm_least_squares((a, _r) in low_rank(), seed in prop::collection::vec(-1.0f64..1.0, 26)) {
        let (m, n) = (a.rows(), a.cols());
        let b: Vec<f64> = seed[..m].to_vec();
        let x = lstsq_min_norm(&a, &b, RANK_TOL);
        prop_assert_eq!(x.len(), n);
        prop_assert!(x.iter().all(|v| v.is_finite()));
        let an = a.frobenius_norm();
        // Least squares: Aᵀ(Ax − b) = 0.
        let ax = a.mul_vec(&x);
        let res: Vec<f64> = ax.iter().zip(&b).map(|(p, q)| p - q).collect();
        let g = a.tr_mul_vec(&res);
        let xn = openshape_linalg::norm2(&x);
        let bn = openshape_linalg::norm2(&b);
        let gn = openshape_linalg::norm2(&g);
        prop_assert!(gn <= 1e-10 * (an * an * xn + an * bn) + 1e-300, "‖Aᵀr‖ = {gn:e}");
        // Minimum norm: x is orthogonal to the null space.
        let nsp = null_space(&a, RANK_TOL);
        let proj = openshape_linalg::projection_norm(&nsp, &x);
        prop_assert!(proj <= 1e-10 * xn + 1e-300, "null-space component {proj:e}");
        // Agrees with the SVD pseudo-inverse.
        let d = svd(&a);
        let k = d.rank(RANK_TOL);
        let mut xs = vec![0.0; n];
        for i in 0..k {
            let ui = d.u.col(i);
            let c = openshape_linalg::dot(ui, &b) / d.s[i];
            openshape_linalg::axpy(c, d.v.col(i), &mut xs);
        }
        let diff: Vec<f64> = x.iter().zip(&xs).map(|(p, q)| p - q).collect();
        prop_assert!(openshape_linalg::norm2(&diff) <= 1e-8 * xn.max(1e-300) + 1e-300,
            "x differs from the SVD pseudo-inverse solution");
    }

    #[test]
    fn null_space_is_orthonormal_and_annihilated((a, _r) in low_rank()) {
        let nsp = null_space(&a, RANK_TOL);
        let rank = QrPivoted::new(&a).rank(RANK_TOL);
        prop_assert_eq!(nsp.cols(), a.cols() - rank);
        prop_assert!(orthonormality_error(&nsp) <= 1e-13 * (a.cols().max(1) as f64));
        if nsp.cols() > 0 && a.rows() > 0 {
            let an = nsp.cols() as f64 * a.frobenius_norm();
            let prod = a.matmul(&nsp);
            prop_assert!(prod.frobenius_norm() <= 1e-11 * an + 1e-300);
        }
    }

    #[test]
    fn weighted_min_norm_optimality(
        (a, _r) in low_rank(),
        seed in prop::collection::vec(-1.0f64..1.0, 26),
        wexp in prop::collection::vec(-3i32..=6, 26),
    ) {
        let (m, n) = (a.rows(), a.cols());
        let b: Vec<f64> = seed[..m].to_vec();
        let w: Vec<f64> = wexp[..n].iter().map(|e| 10f64.powi(*e)).collect();
        let sol = weighted_min_norm(&a, &b, &w, RANK_TOL);
        let x = &sol.x;
        prop_assert!(x.iter().all(|v| v.is_finite()));
        // Least squares in A.
        let ax = a.mul_vec(x);
        let res: Vec<f64> = ax.iter().zip(&b).map(|(p, q)| p - q).collect();
        let g = a.tr_mul_vec(&res);
        let an = a.frobenius_norm();
        let wmax = w.iter().cloned().fold(0.0, f64::max);
        let wmin = w.iter().cloned().fold(f64::INFINITY, f64::min);
        let cond_w = (wmax / wmin).sqrt();
        let xn = openshape_linalg::norm2(x);
        let bn = openshape_linalg::norm2(&b);
        prop_assert!(openshape_linalg::norm2(&g) <= 1e-9 * cond_w * (an * an * xn + an * bn) + 1e-300);
        // Optimality of ‖W^{1/2}x‖ over x + null(A): W x ⟂ null(A W^{-1/2})
        // mapped back, i.e. W x ⟂ null(A) when null(A) is computed for the
        // weighted matrix's rank decision. Use the weighted null space.
        let s: Vec<f64> = w.iter().map(|v| 1.0 / v.sqrt()).collect();
        let mut aw = a.clone();
        aw.scale_cols(&s);
        let nw = null_space(&aw, RANK_TOL);
        // y = W^{1/2} x must be orthogonal to null(A W^{-1/2}).
        let y: Vec<f64> = x.iter().zip(&w).map(|(xi, wi)| xi * wi.sqrt()).collect();
        let yn = openshape_linalg::norm2(&y);
        prop_assert!(openshape_linalg::projection_norm(&nw, &y) <= 1e-10 * yn + 1e-300);
    }

    #[test]
    fn ordered_selection_spans_in_order((a, _r) in low_rank()) {
        let rank = svd(&a).rank(RANK_TOL);
        let sel = select_columns_in_order(&a, RANK_TOL, None);
        prop_assert_eq!(sel.selected.len(), rank);
        let forced = select_columns_in_order(&a, RANK_TOL, Some(rank));
        prop_assert_eq!(&forced, &sel);
        // The selected columns are independent.
        let s = a.select_cols(&sel.selected);
        prop_assert_eq!(svd(&s).rank(RANK_TOL), sel.selected.len());
        // Each rejected column lies in the span of the selected ones that
        // precede it (the greedy, lexicographically first basis).
        let an = a.frobenius_norm();
        for &j in &sel.rejected {
            let before: Vec<usize> = sel.selected.iter().cloned().filter(|&k| k < j).collect();
            let sub = a.select_cols(&before);
            let coef = lstsq_min_norm(&sub, a.col(j), RANK_TOL);
            let fit = sub.mul_vec(&coef);
            let d: Vec<f64> = fit.iter().zip(a.col(j)).map(|(p, q)| p - q).collect();
            prop_assert!(openshape_linalg::norm2(&d) <= 1e-9 * an + 1e-300);
        }
    }

    #[test]
    fn results_are_bit_identical((a, _r) in low_rank(), seed in prop::collection::vec(-1.0f64..1.0, 26)) {
        let b: Vec<f64> = seed[..a.rows()].to_vec();
        let q1 = QrPivoted::new(&a);
        let q2 = QrPivoted::new(&a);
        prop_assert_eq!(bits(&q1.r()), bits(&q2.r()));
        prop_assert_eq!(bits(&q1.q()), bits(&q2.q()));
        prop_assert_eq!(q1.permutation(), q2.permutation());
        let x1 = lstsq_min_norm(&a, &b, RANK_TOL);
        let x2 = lstsq_min_norm(&a, &b, RANK_TOL);
        prop_assert_eq!(
            x1.iter().map(|v| v.to_bits()).collect::<Vec<_>>(),
            x2.iter().map(|v| v.to_bits()).collect::<Vec<_>>()
        );
        let d1 = svd(&a);
        let d2 = svd(&a);
        prop_assert_eq!(bits(&d1.u), bits(&d2.u));
        prop_assert_eq!(bits(&d1.v), bits(&d2.v));
        prop_assert_eq!(
            d1.s.iter().map(|v| v.to_bits()).collect::<Vec<_>>(),
            d2.s.iter().map(|v| v.to_bits()).collect::<Vec<_>>()
        );
        prop_assert_eq!(bits(&null_space(&a, RANK_TOL)), bits(&null_space(&a, RANK_TOL)));
    }
}
