//! Global curve interpolation (acceptance test 5).

use glam::DVec3;
use opendatum_nurbs::{NurbsCurve, NurbsError};
use proptest::prelude::*;

/// Chord-length parameters, as used by `interpolate`.
fn chord_params(points: &[DVec3]) -> Vec<f64> {
    let chords: Vec<f64> = points.windows(2).map(|w| (w[1] - w[0]).length()).collect();
    let total: f64 = chords.iter().sum();
    let mut out = vec![0.0];
    let mut acc = 0.0;
    for c in &chords[..chords.len() - 1] {
        acc += c;
        out.push(acc / total);
    }
    out.push(1.0);
    out
}

fn arb_points() -> impl Strategy<Value = (usize, Vec<DVec3>)> {
    (1usize..=5).prop_flat_map(|p| {
        (
            Just(p),
            prop::collection::vec(
                (-10.0..10.0f64, -10.0..10.0f64, -10.0..10.0f64)
                    .prop_map(|(x, y, z)| DVec3::new(x, y, z)),
                (p + 1)..=(p + 25),
            ),
        )
    })
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(512))]

    #[test]
    fn interpolant_passes_through_all_points((p, points) in arb_points()) {
        // Skip degenerate draws with (nearly) repeated consecutive points.
        prop_assume!(points.windows(2).all(|w| (w[1] - w[0]).length() > 1e-3));
        let curve = NurbsCurve::interpolate(&points, p).unwrap();
        prop_assert_eq!(curve.degree(), p);
        prop_assert_eq!(curve.domain(), (0.0, 1.0));
        prop_assert_eq!(curve.control_points().len(), points.len());
        prop_assert!(!curve.is_rational());
        let params = chord_params(&points);
        for (q, &t) in points.iter().zip(&params) {
            let c = curve.point(t);
            prop_assert!((c - *q).length() <= 1e-10 * q.length().max(1.0), "at t={}: {:?} vs {:?}", t, c, q);
        }
        prop_assert!((curve.point(0.0) - points[0]).length() <= 1e-12 * 10.0);
        prop_assert!((curve.point(1.0) - points[points.len() - 1]).length() <= 1e-12 * 10.0);
    }
}

#[test]
fn interpolating_many_points_on_a_helix() {
    let points: Vec<DVec3> = (0..=1000)
        .map(|i| {
            let s = i as f64 * 0.05;
            DVec3::new(s.cos(), s.sin(), 0.1 * s)
        })
        .collect();
    for p in 1..=5 {
        let curve = NurbsCurve::interpolate(&points, p).unwrap();
        let params = chord_params(&points);
        for (q, &t) in points.iter().zip(&params) {
            assert!((curve.point(t) - *q).length() <= 1e-10);
        }
        if p == 3 {
            // A cubic through dense helix samples approximates the helix closely.
            let mid = curve.point(0.5 * (params[500] + params[501]));
            let s: f64 = 500.5 * 0.05;
            assert!((mid - DVec3::new(s.cos(), s.sin(), 0.1 * s)).length() < 1e-6);
        }
    }
}

#[test]
fn degree_one_interpolation_is_the_polyline() {
    let points = [
        DVec3::ZERO,
        DVec3::new(1., 0., 0.),
        DVec3::new(1., 2., 0.),
        DVec3::new(1., 2., 2.),
    ];
    let curve = NurbsCurve::interpolate(&points, 1).unwrap();
    for (a, b) in curve.control_points().iter().zip(&points) {
        assert!((*a - *b).length() < 1e-14);
    }
    assert!((curve.length() - 5.0).abs() < 1e-12);
}

#[test]
fn interpolation_errors() {
    let pts = [DVec3::ZERO, DVec3::X, DVec3::Y];
    assert_eq!(
        NurbsCurve::interpolate(&pts, 0),
        Err(NurbsError::ZeroDegree)
    );
    assert_eq!(
        NurbsCurve::interpolate(&pts, 3),
        Err(NurbsError::TooFewPoints {
            required: 4,
            found: 3
        })
    );
    assert!(NurbsCurve::interpolate(&[], 1).is_err());
    assert!(matches!(
        NurbsCurve::interpolate(&[DVec3::X, DVec3::X, DVec3::X], 2),
        Err(NurbsError::InvalidArgument(_))
    ));
    assert!(matches!(
        NurbsCurve::interpolate(&[DVec3::ZERO, DVec3::X, DVec3::X, DVec3::Y], 2),
        Err(NurbsError::InvalidArgument(_))
    ));
    assert!(matches!(
        NurbsCurve::interpolate(&[DVec3::ZERO, DVec3::splat(f64::NAN), DVec3::Y], 2),
        Err(NurbsError::NonFiniteControlPoint { index: 1 })
    ));
}
