//! Curve point inversion (acceptance test 6).

mod common;

use common::arb_curve;
use glam::DVec3;
use opendatum_nurbs::NurbsCurve;
use proptest::prelude::*;

fn dense_min_distance(curve: &NurbsCurve, p: DVec3, count: usize) -> f64 {
    let (a, b) = curve.domain();
    (0..count)
        .map(|i| {
            curve
                .point(a + (b - a) * i as f64 / (count - 1) as f64)
                .distance(p)
        })
        .fold(f64::INFINITY, f64::min)
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(1024))]

    #[test]
    fn closest_point_beats_dense_sampling(
        curve in arb_curve(),
        x in -20.0..20.0f64,
        y in -20.0..20.0f64,
        z in -20.0..20.0f64,
    ) {
        let p = DVec3::new(x, y, z);
        let (t, q) = curve.closest_point(p);
        let (a, b) = curve.domain();
        prop_assert!(t >= a && t <= b);
        prop_assert!((curve.point(t) - q).length() <= 1e-12 * q.length().max(1.0));
        let dense = dense_min_distance(&curve, p, 2000);
        let got = q.distance(p);
        prop_assert!(got <= dense + 1e-9, "closest {} > dense {}", got, dense);
    }

    /// Within a range of the parameter, the closest point keeps to it and
    /// is as near as dense samples of that part.
    #[test]
    fn closest_point_in_a_range_beats_dense_sampling_there(
        curve in arb_curve(),
        r0 in 0.0..=1.0f64,
        r1 in 0.0..=1.0f64,
        x in -20.0..20.0f64,
        y in -20.0..20.0f64,
        z in -20.0..20.0f64,
    ) {
        let p = DVec3::new(x, y, z);
        let (a, b) = curve.domain();
        let (t0, t1) = (a + (b - a) * r0, a + (b - a) * r1);
        let (t, q) = curve.closest_point_in(p, t0, t1);
        prop_assert!(t >= t0.min(t1) && t <= t0.max(t1), "{t} outside [{t0}, {t1}]");
        prop_assert!((curve.point(t) - q).length() <= 1e-12 * q.length().max(1.0));
        let dense = (0..2000)
            .map(|i| curve.point(t0 + (t1 - t0) * i as f64 / 1999.0).distance(p))
            .fold(f64::INFINITY, f64::min);
        let got = q.distance(p);
        prop_assert!(got <= dense + 1e-9, "closest {} > dense {}", got, dense);
    }

    #[test]
    fn points_on_the_curve_map_back(curve in arb_curve(), rel in 0.0..=1.0f64) {
        let (a, b) = curve.domain();
        let t0 = a + (b - a) * rel;
        let p = curve.point(t0);
        let (t, q) = curve.closest_point(p);
        let scale = p.length().max(1.0);
        prop_assert!(q.distance(p) <= 1e-10 * scale, "distance {}", q.distance(p));
        // Unless the curve passes through p again elsewhere (self-intersection),
        // the parameter is recovered accurately too.
        let same_branch = (t - t0).abs() < 1e-3 * (b - a);
        if same_branch && curve.derivatives(t0, 1)[1].length() > 1e-3 {
            prop_assert!((t - t0).abs() <= 1e-9 * (b - a), "t {} vs t0 {}", t, t0);
        }
    }

    #[test]
    fn off_curve_points_along_the_normal_map_back(
        curve in arb_curve(),
        rel in 0.05..0.95f64,
        offset in 1e-6..1e-2f64,
    ) {
        // A point slightly off the curve in a direction normal to the tangent
        // projects back to (very nearly) the foot point, when the curve is
        // locally nice.
        let (a, b) = curve.domain();
        let t0 = a + (b - a) * rel;
        let d = curve.derivatives(t0, 2);
        prop_assume!(d[1].length() > 1e-2);
        let tangent = d[1].normalize();
        let normal = tangent.any_orthonormal_vector();
        let p = d[0] + normal * offset;
        let (_, q) = curve.closest_point(p);
        prop_assert!(q.distance(p) <= offset * (1.0 + 1e-9) + 1e-12);
    }
}

#[test]
fn closest_point_on_circle_and_line() {
    let circle = NurbsCurve::circular_arc(
        DVec3::ZERO,
        DVec3::X,
        DVec3::Y,
        2.0,
        0.0,
        std::f64::consts::TAU,
    )
    .unwrap();
    for i in 0..36 {
        let ang = i as f64 * 0.17 + 0.01;
        let dir = DVec3::new(ang.cos(), ang.sin(), 0.0);
        for &r in &[0.5, 1.9, 2.0, 3.0, 50.0] {
            let p = dir * r + DVec3::new(0.0, 0.0, 0.3);
            let (t, q) = circle.closest_point(p);
            assert!(
                (q - dir * 2.0).length() < 1e-12,
                "angle {ang}, r {r}: {q:?}"
            );
            assert!((circle.point(t) - q).length() < 1e-15);
        }
    }
    let line = NurbsCurve::line(DVec3::ZERO, DVec3::new(2.0, 0.0, 0.0)).unwrap();
    let (t, q) = line.closest_point(DVec3::new(0.5, 1.0, 1.0));
    assert!((t - 0.25).abs() < 1e-15);
    assert!((q - DVec3::new(0.5, 0.0, 0.0)).length() < 1e-15);
    // Beyond the ends: clamped to the end points.
    let (t, q) = line.closest_point(DVec3::new(-3.0, 1.0, 0.0));
    assert_eq!(t, 0.0);
    assert_eq!(q, DVec3::ZERO);
    let (t, q) = line.closest_point(DVec3::new(5.0, -1.0, 0.0));
    assert_eq!(t, 1.0);
    assert_eq!(q, DVec3::new(2.0, 0.0, 0.0));
}

#[test]
fn closest_point_parameter_accuracy() {
    // Off-curve point with a known foot: the parameter converges tightly.
    let curve = NurbsCurve::interpolate(
        &[
            DVec3::ZERO,
            DVec3::new(1.0, 1.0, 0.0),
            DVec3::new(2.0, 0.5, 0.3),
            DVec3::new(3.0, 2.0, -0.2),
            DVec3::new(4.0, 0.0, 0.0),
        ],
        3,
    )
    .unwrap();
    for i in 1..20 {
        let t0 = i as f64 / 20.0;
        let d = curve.derivatives(t0, 1);
        let n = d[1].normalize().any_orthonormal_vector();
        let p = d[0] + n * 0.01;
        let (t, _) = curve.closest_point(p);
        assert!((t - t0).abs() < 1e-12, "t0 {t0}: t {t}");
    }
}
