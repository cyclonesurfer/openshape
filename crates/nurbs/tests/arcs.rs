//! Exact circular arcs (acceptance test 3).

use std::f64::consts::{FRAC_PI_2, PI, TAU};

use glam::{DMat4, DVec3};
use openshape_nurbs::{NurbsCurve, NurbsError};
use proptest::prelude::*;

const SWEEPS: [f64; 12] = [
    1e-3,
    0.3,
    FRAC_PI_2,
    FRAC_PI_2 + 1e-9,
    2.0,
    PI,
    PI + 0.5,
    1.5 * PI,
    1.5 * PI + 0.01,
    5.0,
    TAU - 1e-6,
    TAU,
];

fn frames() -> Vec<(DVec3, DVec3, DVec3, f64)> {
    vec![
        (DVec3::ZERO, DVec3::X, DVec3::Y, 1.0),
        (DVec3::new(3.0, -2.0, 7.5), DVec3::Z, DVec3::X, 0.25),
        (
            DVec3::new(-1e3, 2e3, 5e2),
            DVec3::new(1.0, 1.0, 0.0),
            DVec3::new(-1.0, 1.0, 1.0),
            1234.5,
        ),
        (
            DVec3::new(0.1, 0.2, 0.3),
            DVec3::new(0.3, -0.5, 0.8),
            DVec3::new(0.2, 0.9, 0.4),
            1e-3,
        ),
    ]
}

/// Orthonormal frame matching what `circular_arc` does with its axes.
fn orthonormal(x: DVec3, y: DVec3) -> (DVec3, DVec3) {
    let x = x.normalize();
    let y = (y - x * y.dot(x)).normalize();
    (x, y)
}

#[test]
fn points_lie_on_the_circle_and_ends_are_exact() {
    for (center, xa, ya, r) in frames() {
        let (x, y) = orthonormal(xa, ya);
        let normal = x.cross(y);
        for &sweep in &SWEEPS {
            for &start in &[0.0, 0.7, -2.0, 3.0 * PI] {
                let arc =
                    NurbsCurve::circular_arc(center, xa, ya, r, start, start + sweep).unwrap();
                assert_eq!(arc.degree(), 2);
                assert_eq!(arc.domain(), (0.0, 1.0));
                for i in 0..=400 {
                    let t = i as f64 / 400.0;
                    let p = arc.point(t);
                    let dist = (p - center).length();
                    assert!(
                        (dist - r).abs() <= 1e-12 * r,
                        "sweep {sweep}: |p - c| = {dist}, r = {r}"
                    );
                    assert!((p - center).dot(normal).abs() <= 1e-12 * r);
                }
                let expected_start = center + (x * start.cos() + y * start.sin()) * r;
                let end = start + sweep;
                let expected_end = center + (x * end.cos() + y * end.sin()) * r;
                assert!(
                    (arc.point(0.0) - expected_start).length() <= 1e-12 * r.max(center.length())
                );
                assert!((arc.point(1.0) - expected_end).length() <= 1e-12 * r.max(center.length()));
                if sweep == TAU {
                    assert_eq!(arc.point(0.0), arc.point(1.0));
                }
            }
        }
    }
}

#[test]
fn arc_length_is_radius_times_sweep() {
    for (center, xa, ya, r) in frames() {
        for &sweep in &SWEEPS {
            let arc = NurbsCurve::circular_arc(center, xa, ya, r, 0.4, 0.4 + sweep).unwrap();
            let len = arc.length();
            assert!(
                (len - r * sweep).abs() <= 1e-9 * r * sweep.max(1.0),
                "sweep {sweep}: length {len}, expected {}",
                r * sweep
            );
            // Partial lengths add up.
            let a = arc.arc_length(0.0, 0.37);
            let b = arc.arc_length(0.37, 1.0);
            assert!((a + b - len).abs() <= 1e-12 * len.max(1.0));
            assert!((arc.arc_length(0.37, 0.0) - a).abs() <= 1e-15 * len.max(1.0));
        }
    }
}

#[test]
fn arc_tangents_are_perpendicular_to_radius() {
    let arc = NurbsCurve::circular_arc(DVec3::new(1., 2., 3.), DVec3::X, DVec3::Y, 2.0, 0.0, 5.0)
        .unwrap();
    for i in 0..=50 {
        let t = i as f64 / 50.0;
        let tan = arc.tangent(t);
        let radial = (arc.point(t) - DVec3::new(1., 2., 3.)).normalize();
        assert!(tan.dot(radial).abs() < 1e-12);
        // Counter-clockwise about +Z.
        assert!(radial.cross(tan).z > 0.999_999);
    }
}

#[test]
fn end_before_start_wraps_like_a7_1() {
    let a =
        NurbsCurve::circular_arc(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 0.0, -FRAC_PI_2).unwrap();
    // 0 -> 270 degrees counter-clockwise.
    assert!((a.length() - 1.5 * PI).abs() < 1e-12);
    assert!((a.point(1.0) - DVec3::new(0.0, -1.0, 0.0)).length() < 1e-15);
}

#[test]
fn transformed_arc_stays_circular() {
    let arc = NurbsCurve::circular_arc(DVec3::ZERO, DVec3::X, DVec3::Y, 2.0, 0.0, 4.0).unwrap();
    let m = DMat4::from_translation(DVec3::new(5., -1., 2.))
        * DMat4::from_axis_angle(DVec3::new(1., 2., 3.).normalize(), 0.8);
    let moved = arc.transformed(&m);
    let center = DVec3::new(5., -1., 2.);
    for i in 0..=100 {
        let t = i as f64 / 100.0;
        assert!(((moved.point(t) - center).length() - 2.0).abs() < 1e-12);
        assert!((moved.point(t) - m.transform_point3(arc.point(t))).length() < 1e-12);
    }
}

#[test]
fn invalid_arcs_are_rejected() {
    let c = DVec3::ZERO;
    let bad = [
        NurbsCurve::circular_arc(c, DVec3::X, DVec3::Y, 0.0, 0.0, 1.0),
        NurbsCurve::circular_arc(c, DVec3::X, DVec3::Y, -1.0, 0.0, 1.0),
        NurbsCurve::circular_arc(c, DVec3::X, DVec3::Y, f64::NAN, 0.0, 1.0),
        NurbsCurve::circular_arc(c, DVec3::ZERO, DVec3::Y, 1.0, 0.0, 1.0),
        NurbsCurve::circular_arc(c, DVec3::X, DVec3::X * 2.0, 1.0, 0.0, 1.0),
        NurbsCurve::circular_arc(c, DVec3::X, DVec3::ZERO, 1.0, 0.0, 1.0),
        NurbsCurve::circular_arc(c, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0),
        NurbsCurve::circular_arc(c, DVec3::X, DVec3::Y, 1.0, 0.0, 7.0),
        NurbsCurve::circular_arc(c, DVec3::X, DVec3::Y, 1.0, 0.0, f64::INFINITY),
        NurbsCurve::circular_arc(DVec3::splat(f64::NAN), DVec3::X, DVec3::Y, 1.0, 0.0, 1.0),
    ];
    for (i, result) in bad.iter().enumerate() {
        assert!(
            matches!(result, Err(NurbsError::InvalidArgument(_))),
            "case {i}: {result:?}"
        );
    }
}

proptest! {
    #[test]
    fn random_arcs_are_exact(
        start in -10.0..10.0f64,
        sweep in 0.01..TAU,
        r in 1e-3..1e3f64,
        cx in -100.0..100.0f64,
        ax in -1.0..1.0f64,
        ay in -1.0..1.0f64,
        t in 0.0..=1.0f64,
    ) {
        let center = DVec3::new(cx, -cx * 0.5, 3.0);
        let xa = DVec3::new(1.0, ax, ay);
        let ya = DVec3::new(-ax, 1.0, ay * 0.5);
        let arc = NurbsCurve::circular_arc(center, xa, ya, r, start, start + sweep).unwrap();
        let p = arc.point(t);
        prop_assert!((((p - center).length()) - r).abs() <= 1e-12 * r.max(cx.abs()));
        prop_assert!((arc.length() - r * sweep).abs() <= 1e-9 * r * sweep.max(1.0));
    }
}
