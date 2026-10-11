//! Degrees of freedom and free/constrained entity states (design note 0005
//! §8: known DOF counts).

mod common;

use common::{assert_satisfied, B};
use openshape_sketch_solver::{
    solve, Constraint, DVec2, DistanceKind, EntityState, Geometry, SolveOptions, SolveStatus,
};

fn opts() -> SolveOptions<u32> {
    SolveOptions::default()
}

/// A rectangle from four lines with their own end points, joined by
/// coincident constraints, with horizontal/vertical sides.
/// Returns (builder, lines [bottom, right, top, left], corner points).
pub fn rectangle(w: f64, h: f64) -> (B, [u32; 4], [u32; 8]) {
    let mut b = B::new();
    let (l0, p0, p1) = b.line(0.0, 0.0, w, 0.0);
    let (l1, p2, p3) = b.line(w, 0.0, w, h);
    let (l2, p4, p5) = b.line(w, h, 0.0, h);
    let (l3, p6, p7) = b.line(0.0, h, 0.0, 0.0);
    b.coincident(p1, p2);
    b.coincident(p3, p4);
    b.coincident(p5, p6);
    b.coincident(p7, p0);
    b.c(Constraint::Horizontal { line: l0 });
    b.c(Constraint::Vertical { line: l1 });
    b.c(Constraint::Horizontal { line: l2 });
    b.c(Constraint::Vertical { line: l3 });
    (b, [l0, l1, l2, l3], [p0, p1, p2, p3, p4, p5, p6, p7])
}

#[test]
fn lone_point_has_two() {
    let mut b = B::new();
    let p = b.point(1.0, 2.0);
    let s = solve(&b.build(), &opts());
    assert_eq!(s.status, SolveStatus::Solved);
    assert_eq!(s.dof, 2);
    assert_eq!(s.entities[&p], EntityState::Free);
    assert_eq!(s.sketch, b.build(), "nothing to do: geometry unchanged");
    assert_eq!(s.iterations, 0);
}

#[test]
fn line_has_four() {
    let mut b = B::new();
    let (l, p, q) = b.line(0.0, 0.0, 1.0, 0.5);
    let s = solve(&b.build(), &opts());
    assert_eq!(s.dof, 4);
    for k in [l, p, q] {
        assert_eq!(s.entities[&k], EntityState::Free);
    }
}

#[test]
fn circle_has_three_arc_has_five() {
    let mut b = B::new();
    let (c, cc) = b.circle(0.0, 0.0, 1.0);
    let s = solve(&b.build(), &opts());
    assert_eq!(s.dof, 3);
    assert_eq!(s.entities[&c], EntityState::Free);
    assert_eq!(s.entities[&cc], EntityState::Free);

    let mut b = B::new();
    let (a, ..) = b.arc(DVec2::new(1.0, 1.0), 2.0, 0.1, 1.9);
    let s = solve(&b.build(), &opts());
    assert_eq!(s.status, SolveStatus::Solved);
    assert_eq!(s.dof, 5);
    assert_eq!(s.entities[&a], EntityState::Free);
}

#[test]
fn triangle_with_coincident_corners_has_six() {
    let mut b = B::new();
    let (_, a0, a1) = b.line(0.0, 0.0, 3.0, 0.0);
    let (_, b0, b1) = b.line(3.0, 0.0, 1.0, 2.0);
    let (_, c0, c1) = b.line(1.0, 2.0, 0.0, 0.0);
    b.coincident(a1, b0);
    b.coincident(b1, c0);
    b.coincident(c1, a0);
    let s = solve(&b.build(), &opts());
    assert_eq!(s.status, SolveStatus::Solved);
    assert_eq!(s.dof, 6);
}

#[test]
fn rectangle_has_four_then_zero() {
    let (b, lines, pts) = rectangle(3.0, 2.0);
    let s = solve(&b.build(), &opts());
    assert_eq!(s.status, SolveStatus::Solved);
    assert_eq!(s.dof, 4);
    assert_eq!(s.entities[&lines[0]], EntityState::Free);

    let mut b2 = b.clone();
    b2.c(Constraint::Length {
        line: lines[0],
        value: 4.0,
    });
    b2.c(Constraint::Distance {
        a: pts[2],
        b: pts[3],
        kind: DistanceKind::Aligned,
        value: 1.5,
    });
    b2.fix(pts[0]);
    let s = solve(&b2.build(), &opts());
    assert_eq!(s.status, SolveStatus::Solved, "{:?}", s.status);
    assert_eq!(s.dof, 0);
    for (k, st) in &s.entities {
        assert_eq!(*st, EntityState::Constrained, "entity {k}");
    }
    assert_satisfied(&s.sketch);
    let width = (s.sketch.point(pts[1]).unwrap() - s.sketch.point(pts[0]).unwrap()).length();
    assert!((width - 4.0).abs() < 1e-12);
}

#[test]
fn partially_constrained_rectangle_reports_free_entities() {
    // Fix one corner and the width: the top edge can still move up and down.
    let (mut b, lines, pts) = rectangle(3.0, 2.0);
    b.fix(pts[0]);
    b.c(Constraint::Length {
        line: lines[0],
        value: 3.0,
    });
    let s = solve(&b.build(), &opts());
    assert_eq!(s.status, SolveStatus::Solved);
    assert_eq!(s.dof, 1);
    assert_eq!(s.entities[&lines[0]], EntityState::Constrained);
    assert_eq!(s.entities[&pts[0]], EntityState::Constrained);
    assert_eq!(s.entities[&pts[1]], EntityState::Constrained);
    assert_eq!(s.entities[&lines[2]], EntityState::Free);
    assert_eq!(s.entities[&lines[1]], EntityState::Free);
    assert_eq!(s.entities[&pts[4]], EntityState::Free);
}

#[test]
fn point_on_circle_and_concentric() {
    let mut b = B::new();
    let (c, cc) = b.circle(0.0, 0.0, 1.0);
    let p = b.point(1.0, 0.0);
    b.c(Constraint::PointOnCurve { point: p, curve: c });
    let s = solve(&b.build(), &opts());
    assert_eq!(s.dof, 4);

    let (c2, _) = b.circle(0.1, 0.0, 2.0);
    b.c(Constraint::Concentric { a: c, b: c2 });
    b.fix(cc);
    b.c(Constraint::Radius {
        curve: c,
        value: 1.0,
    });
    let s = solve(&b.build(), &opts());
    assert_eq!(s.status, SolveStatus::Solved);
    // c2's radius and the point's angle remain.
    assert_eq!(s.dof, 2);
    assert_eq!(s.entities[&c], EntityState::Constrained);
    assert_eq!(s.entities[&c2], EntityState::Free);
    assert_eq!(s.entities[&p], EntityState::Free);
    assert_satisfied(&s.sketch);
}

#[test]
fn point_on_spline_slides() {
    let mut b = B::new();
    let ctrl: Vec<u32> = [(0.0, 0.0), (1.0, 1.0), (2.0, -1.0), (3.0, 0.0)]
        .iter()
        .map(|&(x, y)| b.point(x, y))
        .collect();
    for &k in &ctrl {
        b.fix(k);
    }
    let sp = b.entity(Geometry::ControlSpline {
        control: ctrl.clone(),
        degree: 3,
    });
    let p = b.point(1.5, 0.3);
    b.c(Constraint::PointOnCurve {
        point: p,
        curve: sp,
    });
    let s = solve(&b.build(), &opts());
    assert_eq!(s.status, SolveStatus::Solved, "{:?}", s.status);
    assert_eq!(s.dof, 1, "the point slides along the fixed spline");
    assert_eq!(s.entities[&p], EntityState::Free);
    assert_eq!(s.entities[&sp], EntityState::Constrained);
    assert_satisfied(&s.sketch);
}

#[test]
fn interpolated_spline_has_two_per_point() {
    let mut b = B::new();
    let pts: Vec<u32> = [(0.0, 0.0), (1.0, 1.0), (2.0, 0.5), (3.0, 1.5), (4.0, 0.0)]
        .iter()
        .map(|&(x, y)| b.point(x, y))
        .collect();
    let sp = b.entity(Geometry::InterpolatedSpline {
        through: pts.clone(),
        params: vec![0.0, 0.2, 0.5, 0.7, 1.0],
    });
    let s = solve(&b.build(), &opts());
    assert_eq!(s.dof, 10);
    assert_eq!(s.entities[&sp], EntityState::Free);
}
