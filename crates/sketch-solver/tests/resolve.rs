//! `resolve_curves` and `measure` (design note 0005 §2, §4): explicit curve
//! geometry for the kernel, and dimension values measured from geometry
//! with the solver's conventions.

mod common;

use common::{de_boor, B};
use openshape_sketch_solver::{
    measure, resolve_curves, Constraint, DVec2, DistanceKind, Geometry, ResolvedCurve, Sketch,
    SolveError,
};
use std::f64::consts::{FRAC_PI_2, PI, TAU};

fn curves(s: &Sketch<u32>) -> Vec<(u32, bool, ResolvedCurve)> {
    resolve_curves(s).expect("resolvable")
}

fn arc_of(c: &ResolvedCurve) -> (DVec2, f64, f64, f64) {
    match c {
        ResolvedCurve::Arc {
            center,
            radius,
            start_angle,
            sweep_angle,
        } => (*center, *radius, *start_angle, *sweep_angle),
        other => panic!("not an arc: {other:?}"),
    }
}

fn spline_of(c: &ResolvedCurve) -> (usize, Vec<f64>, Vec<DVec2>) {
    match c {
        ResolvedCurve::Spline {
            degree,
            knots,
            points,
        } => (*degree, knots.clone(), points.clone()),
        other => panic!("not a spline: {other:?}"),
    }
}

#[test]
fn lines_and_circles_in_entity_order() {
    let mut b = B::new();
    let (c, _) = b.circle(1.0, 2.0, 0.5);
    let lone = b.point(9.0, 9.0);
    let (l, ..) = b.line(0.0, 0.0, 3.0, 4.0);
    b.sketch.entities.last_mut().unwrap().construction = true;
    let out = curves(&b.build());
    assert_eq!(out.len(), 2, "points are skipped");
    assert_eq!(
        out[0],
        (
            c,
            false,
            ResolvedCurve::Circle {
                center: DVec2::new(1.0, 2.0),
                radius: 0.5
            }
        )
    );
    assert_eq!(
        out[1],
        (
            l,
            true,
            ResolvedCurve::Line {
                start: DVec2::ZERO,
                end: DVec2::new(3.0, 4.0)
            }
        )
    );
    assert!(out.iter().all(|(k, _, _)| *k != lone));
}

#[test]
fn arcs_run_counter_clockwise_from_start_to_end() {
    let c = DVec2::new(1.0, -1.0);
    for (a0, a1, sweep) in [
        (0.1, 1.9, 1.8),
        (1.9, 0.1, TAU - 1.8),
        (3.0, -3.0, TAU - 6.0),
        (-0.5, 0.5, 1.0),
        (-3.0, 3.0, 6.0),
        (0.0, PI, PI),
    ] {
        let mut b = B::new();
        b.arc(c, 2.0, a0, a1);
        let out = curves(&b.build());
        let (center, radius, start, sw) = arc_of(&out[0].2);
        assert_eq!(center, c);
        assert!((radius - 2.0).abs() < 1e-15);
        assert!((start - a0).abs() < 1e-14, "start {start} for {a0}");
        assert!((sw - sweep).abs() < 1e-14, "sweep {sw} for {a0} → {a1}");
        assert!(sw > 0.0 && sw < TAU);
    }
    // A start point straight left of the centre is at π, not −π.
    let mut b = B::new();
    let cc = b.point(0.0, 0.0);
    let s = b.point(-1.0, 0.0);
    let e = b.point(0.0, -1.0);
    b.entity(Geometry::Arc {
        center: cc,
        start: s,
        end: e,
    });
    let (_, _, start, sweep) = arc_of(&curves(&b.build())[0].2);
    assert_eq!(start, PI);
    assert!((sweep - FRAC_PI_2).abs() < 1e-15);
    // The radius is the start point's; the end only gives the angle.
    let mut b = B::new();
    let cc = b.point(0.0, 0.0);
    let s = b.point(2.0, 0.0);
    let e = b.point(0.0, 5.0);
    b.entity(Geometry::Arc {
        center: cc,
        start: s,
        end: e,
    });
    let (_, radius, start, sweep) = arc_of(&curves(&b.build())[0].2);
    assert_eq!((radius, start), (2.0, 0.0));
    assert!((sweep - FRAC_PI_2).abs() < 1e-15);
}

#[test]
fn degenerate_curves_are_reported() {
    let degenerate = |s: &Sketch<u32>, what: &str| match resolve_curves(s) {
        Err(SolveError::Degenerate(msg)) => assert!(msg.contains(what), "{msg}"),
        other => panic!("expected Degenerate({what}), got {other:?}"),
    };
    let mut b = B::new();
    b.line(1.0, 1.0, 1.0, 1.0);
    degenerate(&b.build(), "zero length");

    let arc = |c: (f64, f64), s: (f64, f64), e: (f64, f64)| {
        let mut b = B::new();
        let cc = b.point(c.0, c.1);
        let sp = b.point(s.0, s.1);
        let ep = b.point(e.0, e.1);
        b.entity(Geometry::Arc {
            center: cc,
            start: sp,
            end: ep,
        });
        b.build()
    };
    degenerate(&arc((0.0, 0.0), (0.0, 0.0), (1.0, 0.0)), "zero radius");
    degenerate(&arc((0.0, 0.0), (1.0, 0.0), (0.0, 0.0)), "zero radius");
    degenerate(
        &arc((0.0, 0.0), (1.0, 0.0), (2.0, 0.0)),
        "ends where it starts",
    );
}

#[test]
fn control_splines_get_clamped_uniform_knots() {
    for (n, degree, knots) in [
        (4, 3, vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0]),
        (5, 3, vec![0.0, 0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0, 1.0]),
        (6, 2, vec![0.0, 0.0, 0.0, 0.25, 0.5, 0.75, 1.0, 1.0, 1.0]),
        (3, 1, vec![0.0, 0.0, 0.5, 1.0, 1.0]),
        (6, 5, vec![0.0; 6].into_iter().chain([1.0; 6]).collect()),
    ] {
        let mut b = B::new();
        let pts: Vec<DVec2> = (0..n)
            .map(|i| DVec2::new(i as f64, ((i * 7) % 3) as f64 - 1.0))
            .collect();
        let keys: Vec<u32> = pts.iter().map(|p| b.point(p.x, p.y)).collect();
        b.entity(Geometry::ControlSpline {
            control: keys,
            degree,
        });
        let (d, k, points) = spline_of(&curves(&b.build())[0].2);
        assert_eq!(d, degree as usize);
        assert_eq!(k, knots, "{n} points, degree {degree}");
        assert_eq!(points, pts);
    }
}

#[test]
fn interpolated_splines_pass_through_their_points() {
    // (through points, parameters, expected degree, expected knots)
    type Shape = (Vec<(f64, f64)>, Vec<f64>, usize, Vec<f64>);
    let shapes: Vec<Shape> = vec![
        // Cubic, knots by averaging: (0.2 + 0.5 + 0.7) / 3.
        (
            vec![(0.0, 0.0), (1.0, 1.0), (2.0, 0.5), (3.0, 1.5), (4.0, 0.0)],
            vec![0.0, 0.2, 0.5, 0.7, 1.0],
            3,
            vec![0.0, 0.0, 0.0, 0.0, 1.4 / 3.0, 1.0, 1.0, 1.0, 1.0],
        ),
        (
            vec![(0.0, 0.0), (1.0, 2.0), (3.0, -1.0), (4.0, 0.5)],
            vec![0.0, 0.3, 0.8, 1.0],
            3,
            vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0],
        ),
        // Three points: quadratic.
        (
            vec![(0.0, 0.0), (1.0, 1.0), (2.0, 0.0)],
            vec![0.0, 0.25, 1.0],
            2,
            vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
        ),
        // Two points: a straight segment.
        (
            vec![(0.5, 0.5), (2.0, -1.0)],
            vec![0.0, 1.0],
            1,
            vec![0.0, 0.0, 1.0, 1.0],
        ),
    ];
    for (pts, params, degree, knots) in shapes {
        let mut b = B::new();
        let keys: Vec<u32> = pts.iter().map(|&(x, y)| b.point(x, y)).collect();
        b.entity(Geometry::InterpolatedSpline {
            through: keys,
            params: params.clone(),
        });
        let (d, k, control) = spline_of(&curves(&b.build())[0].2);
        assert_eq!(d, degree);
        assert_eq!(k.len(), knots.len());
        for (x, y) in k.iter().zip(&knots) {
            assert!((x - y).abs() < 1e-15, "{k:?} vs {knots:?}");
        }
        assert_eq!(control.len(), pts.len());
        for (&(x, y), &t) in pts.iter().zip(&params) {
            let at = de_boor(d, &k, &control, t);
            assert!(
                (at - DVec2::new(x, y)).length() < 1e-12,
                "degree {d}: {at} at t = {t}, want ({x}, {y})"
            );
        }
        // The end control points are the end points, exactly.
        assert_eq!(control[0], DVec2::new(pts[0].0, pts[0].1));
        let last = pts[pts.len() - 1];
        assert_eq!(control[control.len() - 1], DVec2::new(last.0, last.1));
    }
}

#[test]
fn bad_references_are_invalid() {
    let mut b = B::new();
    let s = b.point(0.0, 0.0);
    b.entity(Geometry::Line { start: s, end: 77 });
    assert!(matches!(
        resolve_curves(&b.build()),
        Err(SolveError::Invalid(_))
    ));
}

#[test]
fn measured_values_follow_the_solver_conventions() {
    let mut b = B::new();
    let p = b.point(1.0, 2.0);
    let q = b.point(4.0, 6.0);
    let (h, ..) = b.line(0.0, 0.0, 1.0, 0.0);
    let (top, ..) = b.line(0.0, 3.0, 2.0, 3.0);
    let (diag, ..) = b.line(0.0, 0.0, 3.0, 4.0);
    let (up, ..) = b.line(5.0, 5.0, 5.0, 7.0);
    let (back, ..) = b.line(1.0, 1.0, -2.0, 1.0);
    let (circle, _) = b.circle(0.0, 0.0, 1.25);
    let (arc, ..) = b.arc(DVec2::new(1.0, 1.0), 0.75, 0.0, 2.0);
    let far = b.point(5.0, 2.0);
    let s = b.build();
    let m = |c: Constraint<u32>| measure(&s, &c).expect("measurable");
    use Constraint as C;
    let dist = |a, b, kind| C::Distance {
        a,
        b,
        kind,
        value: 1.0,
    };
    assert_eq!(m(dist(p, q, DistanceKind::Aligned)), 5.0);
    assert_eq!(m(dist(q, p, DistanceKind::Aligned)), 5.0);
    assert_eq!(m(dist(p, q, DistanceKind::Horizontal)), 3.0);
    assert_eq!(m(dist(q, p, DistanceKind::Vertical)), 4.0);
    // Point–line: to the infinite line, in either order.
    assert_eq!(m(dist(far, h, DistanceKind::Aligned)), 2.0);
    assert_eq!(m(dist(h, far, DistanceKind::Aligned)), 2.0);
    // Line–line: from the first line's start to the second line.
    assert_eq!(m(dist(h, top, DistanceKind::Aligned)), 3.0);
    assert_eq!(m(dist(top, h, DistanceKind::Aligned)), 3.0);
    assert_eq!(
        m(C::Length {
            line: diag,
            value: 1.0
        }),
        5.0
    );
    // Signed angle from a to b, in (−π, π].
    let angle = |a, b| m(C::Angle { a, b, value: 0.0 });
    assert_eq!(angle(h, up), FRAC_PI_2);
    assert_eq!(angle(up, h), -FRAC_PI_2);
    assert_eq!(angle(h, back), PI);
    assert_eq!(angle(back, h), PI);
    assert_eq!(angle(h, top), 0.0);
    assert!((angle(h, diag) - (4.0f64).atan2(3.0)).abs() < 1e-15);
    assert_eq!(
        m(C::Radius {
            curve: circle,
            value: 1.0
        }),
        1.25
    );
    assert_eq!(
        m(C::Diameter {
            curve: circle,
            value: 1.0
        }),
        2.5
    );
    assert!(
        (m(C::Radius {
            curve: arc,
            value: 1.0
        }) - 0.75)
            .abs()
            < 1e-15
    );
    // The stored value plays no part, not even when it is invalid.
    assert_eq!(
        m(C::Length {
            line: diag,
            value: f64::NAN
        }),
        5.0
    );
    assert_eq!(
        m(C::Length {
            line: diag,
            value: -3.0
        }),
        5.0
    );
}

#[test]
fn measure_errors() {
    let mut b = B::new();
    let (l, ..) = b.line(0.0, 0.0, 1.0, 0.0);
    let (z, ..) = b.line(2.0, 2.0, 2.0, 2.0);
    let p = b.point(0.5, 0.5);
    let s = b.build();
    use Constraint as C;
    let invalid = |c: C<u32>| {
        assert!(
            matches!(measure(&s, &c), Err(SolveError::Invalid(_))),
            "{c:?}"
        )
    };
    invalid(C::Horizontal { line: l });
    invalid(C::Coincident { a: p, b: p });
    invalid(C::Length {
        line: 99,
        value: 1.0,
    });
    invalid(C::Radius {
        curve: l,
        value: 1.0,
    });
    invalid(C::Distance {
        a: p,
        b: l,
        kind: DistanceKind::Vertical,
        value: 1.0,
    });
    let degenerate = |c: C<u32>| {
        assert!(
            matches!(measure(&s, &c), Err(SolveError::Degenerate(_))),
            "{c:?}"
        )
    };
    degenerate(C::Angle {
        a: l,
        b: z,
        value: 0.0,
    });
    degenerate(C::Distance {
        a: p,
        b: z,
        kind: DistanceKind::Aligned,
        value: 1.0,
    });
    // A zero-length line's length is simply 0.
    assert_eq!(
        measure(
            &s,
            &C::Length {
                line: z,
                value: 1.0
            }
        ),
        Ok(0.0)
    );
}
