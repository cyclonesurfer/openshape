//! Serialization (design note 0005 §2): sketches, solutions and resolved
//! curves round-trip through serde (JSON here) exactly, floating-point
//! values to the bit.

mod common;

use common::{bits, B};
use openshape_sketch_solver::{
    resolve_curves, solve, Constraint, DVec2, DistanceKind, EntityState, Geometry, ResolvedCurve,
    Sketch, Solution, SolveError, SolveOptions, SolveStatus, TangentKind,
};
use serde::{de::DeserializeOwned, Serialize};

fn round_trip<T: Serialize + DeserializeOwned>(x: &T) -> T {
    let text = serde_json::to_string(x).expect("serializes");
    serde_json::from_str(&text).expect("deserializes")
}

/// A sketch using every entity kind, every constraint kind and every
/// tangency and distance variant, with awkward floating-point values.
fn everything() -> Sketch<u32> {
    let mut b = B::new();
    let (l0, a0, a1) = b.line(0.1, 0.2, 3.0, 1e-17);
    let (l1, b0, b1) = b.line(3.0, -0.0, 3.000_000_000_000_000_4, 2.0);
    let (l2, ..) = b.line(-1.0, 5.0, 1.0 / 3.0, 5.0);
    let (circle, cc) = b.circle(-2.0, -2.0, std::f64::consts::PI);
    let (arc, ac, as_, ae) = b.arc(DVec2::new(5.0, 5.0), 1.5, 0.25, 2.5);
    let ctrl: Vec<u32> = [(6.0, 0.0), (7.0, 1.0), (8.0, -1.0), (9.0, 0.5)]
        .iter()
        .map(|&(x, y)| b.point(x, y))
        .collect();
    let cs = b.entity(Geometry::ControlSpline {
        control: ctrl.clone(),
        degree: 3,
    });
    let thr: Vec<u32> = [(0.0, -5.0), (1.0, -4.0), (2.0, -5.5)]
        .iter()
        .map(|&(x, y)| b.point(x, y))
        .collect();
    b.entity(Geometry::InterpolatedSpline {
        through: thr.clone(),
        params: vec![0.0, 0.1 + 0.2, 1.0],
    });
    let p = b.point(1e-300, -1e300 / 1e294);
    b.sketch.entities.last_mut().unwrap().construction = true;
    use Constraint as C;
    for c in [
        C::Coincident { a: a1, b: b0 },
        C::PointOnCurve {
            point: p,
            curve: cs,
        },
        C::Horizontal { line: l0 },
        C::Vertical { line: l1 },
        C::HorizontalPoints { a: a0, b: p },
        C::VerticalPoints { a: b1, b: ae },
        C::Parallel { a: l0, b: l2 },
        C::Perpendicular { a: l0, b: l1 },
        C::Equal { a: l1, b: l2 },
        C::Midpoint { point: p, line: l2 },
        C::Concentric { a: circle, b: arc },
        C::Symmetric {
            a: a0,
            b: b1,
            axis: l2,
        },
        C::Tangent {
            a: l1,
            b: arc,
            kind: TangentKind::AtPoint { point: b1 },
        },
        C::Tangent {
            a: l2,
            b: circle,
            kind: TangentKind::LineCircle { center_left: true },
        },
        C::Tangent {
            a: circle,
            b: arc,
            kind: TangentKind::CircleCircle { internal: false },
        },
        C::Normal {
            line: l0,
            curve: circle,
        },
        C::Fix {
            point: cc,
            at: DVec2::new(-2.0, f64::MIN_POSITIVE),
        },
        C::Distance {
            a: a0,
            b: as_,
            kind: DistanceKind::Aligned,
            value: 0.1 + 0.2,
        },
        C::Distance {
            a: ac,
            b: thr[0],
            kind: DistanceKind::Horizontal,
            value: 1e-9,
        },
        C::Distance {
            a: ctrl[0],
            b: thr[2],
            kind: DistanceKind::Vertical,
            value: 123_456.789,
        },
        C::Length {
            line: l2,
            value: 2.0f64.sqrt(),
        },
        C::Angle {
            a: l0,
            b: l1,
            value: -std::f64::consts::FRAC_PI_3,
        },
        C::Radius {
            curve: arc,
            value: 1.5,
        },
        C::Diameter {
            curve: circle,
            value: std::f64::consts::TAU,
        },
    ] {
        b.c(c);
    }
    b.build()
}

#[test]
fn sketches_round_trip_exactly() {
    let s = everything();
    let back = round_trip(&s);
    assert_eq!(back, s);
    assert_eq!(bits(&back), bits(&s));
    // Constraint values to the bit too.
    for (x, y) in s.constraints.iter().zip(&back.constraints) {
        assert_eq!(
            x.constraint.value().map(f64::to_bits),
            y.constraint.value().map(f64::to_bits)
        );
    }
    // Other key types.
    let mut t: Sketch<u64> = Sketch::default();
    t.entities.push(openshape_sketch_solver::Entity {
        id: u64::MAX,
        geometry: Geometry::Point { at: DVec2::ONE },
        construction: false,
    });
    assert_eq!(round_trip(&t), t);
}

#[test]
fn solutions_round_trip_exactly() {
    let s = everything();
    let solved = solve(&s, &SolveOptions::default());
    let back: Solution<u32> = round_trip(&solved);
    assert_eq!(back, solved);
    assert_eq!(bits(&back.sketch), bits(&solved.sketch));
    assert_eq!(back.max_residual.to_bits(), solved.max_residual.to_bits());

    // A conflicting solve, with its sets.
    let mut b = B::new();
    let q = b.point(0.0, 0.0);
    b.c(Constraint::Fix {
        point: q,
        at: DVec2::ZERO,
    });
    b.c(Constraint::Fix {
        point: q,
        at: DVec2::X,
    });
    let failed = solve(&b.build(), &SolveOptions::default());
    assert_eq!(failed.status, SolveStatus::Failed(SolveError::Conflicting));
    assert_eq!(round_trip(&failed), failed);
}

#[test]
fn statuses_and_errors_round_trip() {
    for st in [
        SolveStatus::Solved,
        SolveStatus::Redundant,
        SolveStatus::Failed(SolveError::Invalid("entity 3: «bad» \"input\"".into())),
        SolveStatus::Failed(SolveError::Degenerate("line 4 has zero length".into())),
        SolveStatus::Failed(SolveError::Conflicting),
        SolveStatus::Failed(SolveError::NotConverged),
    ] {
        assert_eq!(round_trip(&st), st);
    }
    for e in [EntityState::Free, EntityState::Constrained] {
        assert_eq!(round_trip(&e), e);
    }
}

#[test]
fn resolved_curves_round_trip_exactly() {
    let s = everything();
    let curves = resolve_curves(&s).expect("resolvable");
    assert!(curves
        .iter()
        .any(|(_, _, c)| matches!(c, ResolvedCurve::Spline { .. })));
    assert!(curves
        .iter()
        .any(|(_, _, c)| matches!(c, ResolvedCurve::Arc { .. })));
    let back: Vec<(u32, bool, ResolvedCurve)> = round_trip(&curves);
    assert_eq!(back, curves);
}
