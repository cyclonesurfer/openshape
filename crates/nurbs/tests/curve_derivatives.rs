//! Rational curve derivatives against finite differences (acceptance test 2).

mod common;

use common::{arb_curve, close, span_width};
use glam::DVec3;
use openshape_nurbs::NurbsCurve;
use proptest::prelude::*;

/// Picks a parameter strictly inside a knot span (relative position 0.1..0.9),
/// where the curve is smooth.
fn param_in_span(curve: &NurbsCurve, which: f64, rel: f64) -> f64 {
    let k = curve.knots().as_slice();
    let p = curve.degree();
    let mut breaks: Vec<f64> = k[p..k.len() - p].to_vec();
    breaks.dedup();
    let span = ((which * (breaks.len() - 1) as f64) as usize).min(breaks.len() - 2);
    breaks[span] + (breaks[span + 1] - breaks[span]) * (0.1 + 0.8 * rel)
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(512))]

    #[test]
    fn derivatives_match_central_differences(
        curve in arb_curve(),
        which in 0.0..1.0f64,
        rel in 0.0..1.0f64,
    ) {
        let t = param_in_span(&curve, which, rel);
        let h = 1e-6 * span_width(curve.knots().as_slice(), t);
        let d = curve.derivatives(t, 3);
        prop_assert_eq!(d.len(), 4);
        prop_assert!((d[0] - curve.point(t)).length() <= 1e-12 * d[0].length().max(1.0));

        let plus = curve.derivatives(t + h, 2);
        let minus = curve.derivatives(t - h, 2);
        let fd1 = (curve.point(t + h) - curve.point(t - h)) / (2.0 * h);
        let fd2 = (plus[1] - minus[1]) / (2.0 * h);
        let fd3 = (plus[2] - minus[2]) / (2.0 * h);
        prop_assert!(close(d[1], fd1, 1e-5), "C' {:?} vs {:?}", d[1], fd1);
        prop_assert!(close(d[2], fd2, 1e-5), "C'' {:?} vs {:?}", d[2], fd2);
        prop_assert!(close(d[3], fd3, 1e-5), "C''' {:?} vs {:?}", d[3], fd3);

        // Second derivative also from point values (looser: cancellation).
        let h2 = 1e-4 * span_width(curve.knots().as_slice(), t);
        let fd2_points = (curve.point(t + h2) - 2.0 * curve.point(t) + curve.point(t - h2)) / (h2 * h2);
        prop_assert!(close(d[2], fd2_points, 1e-3), "C'' {:?} vs {:?}", d[2], fd2_points);
    }

    #[test]
    fn tangent_is_unit_and_parallel(curve in arb_curve(), which in 0.0..1.0f64, rel in 0.0..1.0f64) {
        let t = param_in_span(&curve, which, rel);
        let d1 = curve.derivatives(t, 1)[1];
        let tan = curve.tangent(t);
        if d1.length() > 1e-6 {
            prop_assert!((tan.length() - 1.0).abs() < 1e-12);
            prop_assert!(tan.dot(d1.normalize()) > 1.0 - 1e-12);
        }
    }
}

#[test]
fn derivatives_of_a_polynomial_curve_are_exact() {
    // C(t) = (t, t^2, t^3) on [0, 1] as a cubic Bezier.
    let c = NurbsCurve::new(
        3,
        vec![0., 0., 0., 0., 1., 1., 1., 1.],
        vec![
            DVec3::new(0., 0., 0.),
            DVec3::new(1. / 3., 0., 0.),
            DVec3::new(2. / 3., 1. / 3., 0.),
            DVec3::new(1., 1., 1.),
        ],
        None,
    )
    .unwrap();
    for i in 0..=10 {
        let t = i as f64 / 10.0;
        let d = c.derivatives(t, 4);
        assert!((d[0] - DVec3::new(t, t * t, t * t * t)).length() < 1e-15);
        assert!((d[1] - DVec3::new(1., 2. * t, 3. * t * t)).length() < 1e-14);
        assert!((d[2] - DVec3::new(0., 2., 6. * t)).length() < 1e-13);
        assert!((d[3] - DVec3::new(0., 0., 6.)).length() < 1e-12);
        assert_eq!(d[4], DVec3::ZERO);
    }
}

#[test]
fn rational_circle_derivatives_are_exact() {
    // Quarter circle: |C| = 1, C . C' = 0, |C'|^2 + C . C'' = 0.
    let w = std::f64::consts::FRAC_1_SQRT_2;
    let c = NurbsCurve::new(
        2,
        vec![0., 0., 0., 1., 1., 1.],
        vec![DVec3::X, DVec3::new(1., 1., 0.), DVec3::Y],
        Some(vec![1., w, 1.]),
    )
    .unwrap();
    for i in 0..=10 {
        let t = i as f64 / 10.0;
        let d = c.derivatives(t, 2);
        assert!(d[0].dot(d[1]).abs() < 1e-14);
        assert!((d[1].length_squared() + d[0].dot(d[2])).abs() < 1e-13);
    }
}
