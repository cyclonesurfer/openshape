//! Knot insertion, splitting, reversal and transforms (acceptance test 4).

mod common;

use common::arb_curve;
use glam::{DMat4, DVec3};
use opendatum_nurbs::{NurbsCurve, NurbsError};
use proptest::prelude::*;

fn samples(curve: &NurbsCurve, count: usize) -> Vec<f64> {
    let (a, b) = curve.domain();
    (0..=count)
        .map(|i| a + (b - a) * i as f64 / count as f64)
        .collect()
}

fn scale(curve: &NurbsCurve) -> f64 {
    let (lo, hi) = curve.bounding_box();
    lo.abs().max(hi.abs()).max_element().max(1.0)
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(512))]

    #[test]
    fn knot_insertion_preserves_the_curve(
        curve in arb_curve(),
        rel in 0.0..1.0f64,
        times in 1usize..6,
    ) {
        let (a, b) = curve.domain();
        let t = a + (b - a) * rel;
        let refined = curve.insert_knot(t, times);
        let s0 = curve.knots().multiplicity(t);
        let expected = (s0 + times).min(curve.degree()).max(s0);
        prop_assert_eq!(refined.knots().multiplicity(t), expected);
        prop_assert_eq!(
            refined.homogeneous_points().len(),
            curve.homogeneous_points().len() + (expected - s0)
        );
        prop_assert_eq!(refined.domain(), curve.domain());
        let tol = 1e-12 * scale(&curve);
        for u in samples(&curve, 200) {
            prop_assert!((curve.point(u) - refined.point(u)).length() <= tol);
        }
    }

    #[test]
    fn split_preserves_the_curve_and_pieces_join(curve in arb_curve(), rel in 0.01..0.99f64) {
        let (a, b) = curve.domain();
        let t = a + (b - a) * rel;
        let (left, right) = curve.split(t);
        prop_assert_eq!(left.domain(), (a, t));
        prop_assert_eq!(right.domain(), (t, b));
        let tol = 1e-12 * scale(&curve);
        // The pieces meet at the split point, which is on the curve.
        let at = curve.point(t);
        prop_assert!((left.point(t) - at).length() <= tol);
        prop_assert!((right.point(t) - at).length() <= tol);
        // Clamped ends interpolate the end control points.
        let lp = left.control_points();
        let rp = right.control_points();
        prop_assert!((lp[lp.len() - 1] - at).length() <= tol);
        prop_assert!((rp[0] - at).length() <= tol);
        prop_assert_eq!(lp[lp.len() - 1], rp[0]);
        for u in samples(&curve, 200) {
            let piece = if u <= t { &left } else { &right };
            prop_assert!((curve.point(u) - piece.point(u)).length() <= tol);
        }
        // Lengths add up.
        let total = curve.length();
        prop_assert!((left.length() + right.length() - total).abs() <= 1e-9 * total.max(1.0));
    }

    #[test]
    fn reversal_preserves_the_curve(curve in arb_curve()) {
        let rev = curve.reversed();
        let (a, b) = curve.domain();
        prop_assert_eq!(rev.domain(), (a, b));
        let tol = 1e-12 * scale(&curve);
        // The mirrored parameter a + b - u, and the mirrored knots, are
        // rounded by up to two units in the last place of a + b; a steep
        // curve moves by its speed times that.
        let du = 2.0 * f64::EPSILON * (a + b).abs();
        for u in samples(&curve, 100) {
            let v = a + b - u;
            let slack = rev.derivatives(v, 1)[1].length() * du;
            prop_assert!((curve.point(u) - rev.point(v)).length() <= tol + slack);
        }
        let back = rev.reversed();
        for u in samples(&curve, 50) {
            let slack = curve.derivatives(u, 1)[1].length() * du;
            prop_assert!((curve.point(u) - back.point(u)).length() <= tol + slack);
        }
    }

    #[test]
    fn affine_transform_commutes_with_evaluation(
        curve in arb_curve(),
        angle in -3.0..3.0f64,
        s in 0.1..5.0f64,
        tx in -10.0..10.0f64,
    ) {
        let m = DMat4::from_scale_rotation_translation(
            DVec3::new(s, 1.0 / s, 2.0),
            glam::DQuat::from_axis_angle(DVec3::new(1.0, -2.0, 0.5).normalize(), angle),
            DVec3::new(tx, 1.0, -2.0),
        );
        let moved = curve.transformed(&m);
        prop_assert_eq!(moved.weights(), curve.weights());
        for u in samples(&curve, 50) {
            let expected = m.transform_point3(curve.point(u));
            prop_assert!((moved.point(u) - expected).length() <= 1e-11 * expected.length().max(1.0) * 10.0);
        }
    }
}

#[test]
fn split_at_existing_knot_and_repeated_splits() {
    let curve = NurbsCurve::new(
        3,
        vec![0., 0., 0., 0., 1., 2., 2., 3., 4., 4., 4., 4.],
        (0..8)
            .map(|i| DVec3::new(i as f64, ((i * 7) % 5) as f64, ((i * 3) % 4) as f64))
            .collect(),
        Some(vec![1., 2., 1., 0.5, 1., 3., 1., 1.]),
    )
    .unwrap();
    for &t in &[1.0, 2.0, 3.0, 0.5, 3.5] {
        let (l, r) = curve.split(t);
        for i in 0..=80 {
            let u = 4.0 * i as f64 / 80.0;
            let piece = if u <= t { &l } else { &r };
            assert!((curve.point(u) - piece.point(u)).length() < 1e-12);
        }
        // Split a piece again.
        let (ll, lr) = l.split(t / 2.0);
        assert!((ll.point(t / 2.0) - lr.point(t / 2.0)).length() < 1e-12);
        assert!((ll.point(t / 4.0) - curve.point(t / 4.0)).length() < 1e-12);
    }
}

#[test]
fn a_split_a_rounding_step_from_a_knot_splits_at_the_knot() {
    // A knot one rounding step below 0.5 (as chord-length parameters
    // make: 110/220 summed in pieces), and a split at the midpoint 0.5.
    let knot: f64 = 0.5 - f64::EPSILON / 4.0;
    let t: f64 = 0.5;
    assert!(t != knot && (t - knot).abs() < 4.0 * f64::EPSILON);
    let curve = NurbsCurve::new(
        3,
        vec![0., 0., 0., 0., 0.25, knot, 0.75, 1., 1., 1., 1.],
        (0..7)
            .map(|i| DVec3::new(i as f64, ((i * 7) % 5) as f64, ((i * 3) % 4) as f64))
            .collect(),
        None,
    )
    .unwrap();
    for u in [t, knot - f64::EPSILON / 2.0, knot + 2.0 * f64::EPSILON] {
        let (l, r) = curve.split(u);
        let (el, er) = curve.split(knot);
        assert_eq!(l, el, "{u}");
        assert_eq!(r, er, "{u}");
        // No sliver spans: every nonzero span is a real one.
        for piece in [&l, &r] {
            for w in piece.knots().as_slice().windows(2) {
                assert!(w[1] == w[0] || w[1] - w[0] >= 0.25 - 1e-15, "{u}: {w:?}");
            }
        }
        // End derivatives agree across the cut (the curve is C² there).
        let (dl, dr) = (l.derivatives(knot, 1)[1], r.derivatives(knot, 1)[1]);
        assert!((dl - dr).length() <= 1e-12 * dl.length(), "{dl} {dr}");
    }
    // Away from knots nothing moves.
    let (l, _) = curve.split(0.3);
    assert_eq!(l.domain(), (0.0, 0.3));
}

#[test]
fn split_outside_domain_is_an_error() {
    let curve = NurbsCurve::line(DVec3::ZERO, DVec3::X).unwrap();
    for t in [0.0, 1.0, -1.0, 2.0, f64::NAN] {
        assert!(matches!(
            curve.try_split(t),
            Err(NurbsError::InvalidArgument(_))
        ));
    }
    let (l, r) = curve.split(0.25);
    assert_eq!(l.point(0.25), DVec3::new(0.25, 0.0, 0.0));
    assert_eq!(r.point(1.0), DVec3::X);
}

#[test]
fn unclamped_curve_operations() {
    // Uniform unclamped cubic B-spline (domain [3, 6]).
    let pts: Vec<DVec3> = (0..6)
        .map(|i| DVec3::new(i as f64, (i as f64 * 1.3).sin(), (i as f64 * 0.7).cos()))
        .collect();
    let curve = NurbsCurve::new(3, (0..10).map(f64::from).collect(), pts, None).unwrap();
    assert_eq!(curve.domain(), (3.0, 6.0));
    let refined = curve.insert_knot(3.0, 2);
    let (l, r) = curve.split(4.5);
    let rev = curve.reversed();
    for i in 0..=60 {
        let u = 3.0 + 3.0 * i as f64 / 60.0;
        let p = curve.point(u);
        assert!((refined.point(u) - p).length() < 1e-12);
        let piece = if u <= 4.5 { &l } else { &r };
        assert!((piece.point(u) - p).length() < 1e-12);
        assert!((rev.point(9.0 - u) - p).length() < 1e-12);
    }
}

#[test]
fn insert_knot_ignores_parameters_outside_the_domain() {
    let curve = NurbsCurve::line(DVec3::ZERO, DVec3::X).unwrap();
    assert_eq!(curve.insert_knot(-0.5, 1), curve);
    assert_eq!(curve.insert_knot(1.5, 1), curve);
    assert_eq!(curve.insert_knot(f64::NAN, 1), curve);
    assert_eq!(curve.insert_knot(0.5, 0), curve);
    // Degree 1: at most one copy of an interior knot.
    let refined = curve.insert_knot(0.5, 3);
    assert_eq!(refined.knots().as_slice(), &[0.0, 0.0, 0.5, 1.0, 1.0]);
}
