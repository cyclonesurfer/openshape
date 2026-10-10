//! Basis function properties (acceptance test 1).

use openshape_nurbs::KnotVector;
use proptest::prelude::*;

/// Random clamped knot vector of degree `p` with interior multiplicities in
/// `1..=p`, plus a parameter inside the domain.
fn arb_knots_and_param() -> impl Strategy<Value = (usize, KnotVector, f64)> {
    (1usize..=6)
        .prop_flat_map(|p| {
            (
                Just(p),
                prop::collection::vec((0.05..1.0f64, 1usize..=p), 0..8),
                -3.0..3.0f64,
                0.1..20.0f64,
                0.0..=1.0f64,
            )
        })
        .prop_map(|(p, interior, start, width, rel)| {
            let total: f64 = interior.iter().map(|(d, _)| d).sum::<f64>() + 0.5;
            let mut knots = vec![start; p + 1];
            let mut acc = 0.0;
            for (d, mult) in &interior {
                acc += d;
                let k = start + width * acc / total;
                knots.extend(std::iter::repeat_n(k, *mult));
            }
            knots.extend(std::iter::repeat_n(start + width, p + 1));
            let u = start + width * rel;
            (p, KnotVector::new(knots).unwrap(), u)
        })
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(1024))]

    #[test]
    fn partition_of_unity_and_non_negativity((p, kv, u) in arb_knots_and_param()) {
        let n = kv.len() - p - 2;
        let span = kv.find_span(p, n, u);
        let k = kv.as_slice();
        prop_assert!(span >= p && span <= n);
        // A2.1 contract: u_span <= u < u_{span+1}, except at the domain end.
        prop_assert!(k[span] <= u);
        prop_assert!(u < k[span + 1] || (u == k[n + 1] && k[span] < k[span + 1]));
        let basis = kv.basis_funs(span, u, p);
        prop_assert_eq!(basis.len(), p + 1);
        let sum: f64 = basis.iter().sum();
        prop_assert!((sum - 1.0).abs() < 1e-13, "sum = {}", sum);
        for b in &basis {
            prop_assert!(*b >= 0.0, "negative basis value {}", b);
        }
    }

    #[test]
    fn derivative_basis_consistency((p, kv, u) in arb_knots_and_param()) {
        let n = kv.len() - p - 2;
        let span = kv.find_span(p, n, u);
        let ders = kv.ders_basis_funs(span, u, p, p + 2);
        let basis = kv.basis_funs(span, u, p);
        prop_assert_eq!(ders.len(), p + 3);
        for j in 0..=p {
            prop_assert!((ders[0][j] - basis[j]).abs() < 1e-15);
        }
        // Every derivative of a partition of unity sums to zero.
        for (k, row) in ders.iter().enumerate().skip(1) {
            let scale: f64 = row.iter().map(|d| d.abs()).sum::<f64>().max(1.0);
            let sum: f64 = row.iter().sum();
            prop_assert!(sum.abs() <= 1e-11 * scale, "order {} sums to {}", k, sum);
        }
        // Orders above the degree vanish.
        for row in &ders[p + 1..] {
            prop_assert!(row.iter().all(|&d| d == 0.0));
        }
    }

    #[test]
    fn basis_derivatives_match_finite_differences((p, kv, u) in arb_knots_and_param()) {
        let n = kv.len() - p - 2;
        let span = kv.find_span(p, n, u);
        let k = kv.as_slice();
        let (lo, hi) = (k[span], k[span + 1]);
        // Evaluate inside the span so the basis polynomials are smooth.
        let t = lo + (hi - lo) * (0.25 + 0.5 * ((u - lo) / (hi - lo)).clamp(0.0, 1.0));
        let h = 1e-6 * (hi - lo);
        let ders = kv.ders_basis_funs(span, t, p, 1);
        let plus = kv.basis_funs(span, t + h, p);
        let minus = kv.basis_funs(span, t - h, p);
        let scale = ders[1].iter().map(|d| d.abs()).fold(1.0, f64::max);
        for j in 0..=p {
            let fd = (plus[j] - minus[j]) / (2.0 * h);
            prop_assert!((fd - ders[1][j]).abs() <= 1e-5 * scale, "j={} fd={} exact={}", j, fd, ders[1][j]);
        }
    }
}

#[test]
fn unclamped_knot_vector_basis() {
    // Uniform unclamped cubic: domain [3, 5].
    let kv = KnotVector::new((0..=8).map(f64::from).collect()).unwrap();
    assert_eq!(kv.domain(3), (3.0, 5.0));
    for i in 0..=40 {
        let u = 3.0 + 2.0 * i as f64 / 40.0;
        let span = kv.find_span(3, 4, u);
        let sum: f64 = kv.basis_funs(span, u, 3).iter().sum();
        assert!((sum - 1.0).abs() < 1e-14);
    }
    // Uniform cubic B-spline basis at a knot: 1/6, 2/3, 1/6, 0.
    let b = kv.basis_funs(4, 4.0, 3);
    assert!((b[0] - 1.0 / 6.0).abs() < 1e-15);
    assert!((b[1] - 2.0 / 3.0).abs() < 1e-15);
    assert!((b[2] - 1.0 / 6.0).abs() < 1e-15);
    assert!(b[3].abs() < 1e-15);
}
