//! Dragging (design note 0005 §7, §8): a free point ends exactly on the
//! cursor, a constrained one at the nearest feasible position, a fixed one
//! stays, and bad drag requests are rejected.

mod common;

use common::{assert_satisfied, bits, circ, pos, B};
use opendatum_sketch_solver::{
    solve, Constraint, DVec2, EntityState, Geometry, SolveError, SolveOptions, SolveStatus,
};
use std::f64::consts::TAU;

fn drag(points: &[(u32, DVec2)]) -> SolveOptions<u32> {
    SolveOptions {
        drag: points.to_vec(),
        ..SolveOptions::default()
    }
}

fn close(a: DVec2, b: DVec2, tol: f64) -> bool {
    (a - b).length() <= tol
}

#[test]
fn a_free_point_lands_exactly_on_the_target() {
    let mut b = B::new();
    let p = b.point(1.0, 2.0);
    let (_, s, e) = b.line(0.0, 0.0, 1.0, 1.0);
    let sk = b.build();
    for target in [
        DVec2::new(-3.25, 7.5),
        DVec2::new(1e-7, -1e-7),
        DVec2::new(123.456, 0.001),
    ] {
        for k in [p, s, e] {
            let out = solve(&sk, &drag(&[(k, target)]));
            assert_eq!(out.status, SolveStatus::Solved);
            assert_eq!(pos(&out.sketch, k).to_array().map(f64::to_bits), {
                target.to_array().map(f64::to_bits)
            });
            // Nothing else moves.
            for other in [p, s, e].into_iter().filter(|&o| o != k) {
                assert_eq!(pos(&out.sketch, other), pos(&sk, other));
            }
        }
    }
}

#[test]
fn coincident_points_move_together() {
    // Two lines joined at a corner: dragging the corner point of one line
    // drags the other line's end with it.
    let mut b = B::new();
    let (_, a0, a1) = b.line(0.0, 0.0, 1.0, 0.0);
    let (_, b0, b1) = b.line(1.0, 0.0, 1.0, 1.0);
    b.coincident(a1, b0);
    let sk = b.build();
    let target = DVec2::new(2.0, -0.5);
    let out = solve(&sk, &drag(&[(a1, target)]));
    assert_eq!(out.status, SolveStatus::Solved);
    assert!(close(pos(&out.sketch, a1), target, 1e-14));
    assert!(close(pos(&out.sketch, b0), target, 1e-14));
    assert_eq!(pos(&out.sketch, a0), pos(&sk, a0));
    assert_eq!(pos(&out.sketch, b1), pos(&sk, b1));
}

#[test]
fn a_point_on_a_fixed_circle_goes_to_the_nearest_point() {
    let mut b = B::new();
    let (c, cc) = b.circle(0.5, -0.25, 1.5);
    b.fix(cc);
    b.c(Constraint::Radius {
        curve: c,
        value: 1.5,
    });
    let p = b.point(2.0, -0.25);
    b.c(Constraint::PointOnCurve { point: p, curve: c });
    let sk = b.build();
    let centre = DVec2::new(0.5, -0.25);
    for target in [
        DVec2::new(3.0, 4.0),
        DVec2::new(-2.0, 0.1),
        DVec2::new(0.6, -0.2),
        DVec2::new(0.5, -3.0),
        DVec2::new(1.0, 1.0),
    ] {
        let out = solve(&sk, &drag(&[(p, target)]));
        assert_eq!(out.status, SolveStatus::Solved, "{target}");
        let want = centre + (target - centre).normalize() * 1.5;
        let got = pos(&out.sketch, p);
        assert!(close(got, want, 1e-12), "{target}: {got} vs {want}");
        assert_eq!(circ(&out.sketch, c), (centre, 1.5));
        assert_satisfied(&out.sketch);
    }
}

#[test]
fn a_point_dragged_round_a_circle_follows_the_cursor() {
    // As in the UI: one solve per pointer move, each from the last result.
    let mut b = B::new();
    let (c, cc) = b.circle(0.0, 0.0, 2.0);
    b.fix(cc);
    b.c(Constraint::Diameter {
        curve: c,
        value: 4.0,
    });
    let p = b.point(2.0, 0.0);
    b.c(Constraint::PointOnCurve { point: p, curve: c });
    let mut sk = b.build();
    for i in 1..=72 {
        let a = i as f64 * TAU / 72.0;
        let target = DVec2::new(3.0 * a.cos(), 3.0 * a.sin());
        let out = solve(&sk, &drag(&[(p, target)]));
        assert_eq!(out.status, SolveStatus::Solved);
        let got = pos(&out.sketch, p);
        assert!(close(got, DVec2::new(2.0 * a.cos(), 2.0 * a.sin()), 1e-12));
        sk = out.sketch;
    }
}

#[test]
fn a_line_end_with_a_fixed_length_points_at_the_cursor() {
    let mut b = B::new();
    let (l, s, e) = b.line(1.0, 1.0, 2.0, 1.0);
    b.fix(s);
    b.c(Constraint::Length {
        line: l,
        value: 1.0,
    });
    let sk = b.build();
    let start = DVec2::new(1.0, 1.0);
    for target in [
        DVec2::new(4.0, 5.0),
        DVec2::new(1.0, 0.0),
        DVec2::new(-1.0, 1.2),
    ] {
        let out = solve(&sk, &drag(&[(e, target)]));
        assert_eq!(out.status, SolveStatus::Solved, "{target}");
        let want = start + (target - start).normalize();
        assert!(close(pos(&out.sketch, e), want, 1e-12), "{target}");
        assert_eq!(pos(&out.sketch, s), start);
    }
}

#[test]
fn a_fixed_point_does_not_move() {
    let mut b = B::new();
    let (_, s, e) = b.line(0.0, 0.0, 1.0, 1.0);
    b.fix(s);
    let sk = b.build();
    let out = solve(&sk, &drag(&[(s, DVec2::new(5.0, 5.0))]));
    assert_eq!(out.status, SolveStatus::Solved);
    assert_eq!(pos(&out.sketch, s), DVec2::ZERO);
    assert_eq!(pos(&out.sketch, e), DVec2::ONE);
    assert_eq!(out.entities[&s], EntityState::Constrained);
}

#[test]
fn dragging_a_rectangle_corner_resizes_it() {
    let mut b = B::new();
    let (l0, p0, p1) = b.line(0.0, 0.0, 3.0, 0.0);
    let (l1, p2, p3) = b.line(3.0, 0.0, 3.0, 2.0);
    let (l2, p4, p5) = b.line(3.0, 2.0, 0.0, 2.0);
    let (l3, p6, p7) = b.line(0.0, 2.0, 0.0, 0.0);
    b.coincident(p1, p2);
    b.coincident(p3, p4);
    b.coincident(p5, p6);
    b.coincident(p7, p0);
    b.c(Constraint::Horizontal { line: l0 });
    b.c(Constraint::Vertical { line: l1 });
    b.c(Constraint::Horizontal { line: l2 });
    b.c(Constraint::Vertical { line: l3 });
    b.fix(p0);
    let sk = b.build();
    let target = DVec2::new(5.0, 4.0);
    let out = solve(&sk, &drag(&[(p4, target)]));
    assert_eq!(out.status, SolveStatus::Solved);
    assert_satisfied(&out.sketch);
    assert!(close(pos(&out.sketch, p4), target, 1e-14));
    assert!(close(pos(&out.sketch, p1), DVec2::new(5.0, 0.0), 1e-14));
    assert!(close(pos(&out.sketch, p6), DVec2::new(0.0, 4.0), 1e-14));
    assert_eq!(pos(&out.sketch, p0), DVec2::ZERO);

    // With the width dimensioned, the corner goes as far as it can: the
    // height follows the cursor, the width stays.
    let mut b2 = b.clone();
    b2.c(Constraint::Length {
        line: l0,
        value: 3.0,
    });
    let out = solve(&b2.build(), &drag(&[(p4, target)]));
    assert_eq!(out.status, SolveStatus::Solved);
    assert!(close(pos(&out.sketch, p4), DVec2::new(3.0, 4.0), 1e-12));
    assert_satisfied(&out.sketch);
}

#[test]
fn a_fully_constrained_sketch_does_not_follow() {
    let mut b = B::new();
    let (l, s, e) = b.line(0.0, 0.0, 2.0, 0.0);
    b.fix(s);
    b.c(Constraint::Horizontal { line: l });
    b.c(Constraint::Length {
        line: l,
        value: 2.0,
    });
    let sk = b.build();
    let out = solve(&sk, &drag(&[(e, DVec2::new(2.5, 0.7))]));
    assert_eq!(out.status, SolveStatus::Solved);
    assert!(close(pos(&out.sketch, e), DVec2::new(2.0, 0.0), 1e-12));
    assert_eq!(out.dof, 0);
}

#[test]
fn a_point_on_a_fixed_spline_goes_to_the_nearest_point() {
    let mut b = B::new();
    let ctrl: Vec<u32> = [(0.0, 0.0), (1.0, 2.0), (2.0, -1.0), (3.0, 1.0), (4.0, 0.0)]
        .iter()
        .map(|&(x, y)| b.point(x, y))
        .collect();
    for &k in &ctrl {
        b.fix(k);
    }
    let sp = b.entity(Geometry::ControlSpline {
        control: ctrl,
        degree: 3,
        periodic: false,
    });
    let p = b.point(2.0, 0.2);
    b.c(Constraint::PointOnCurve {
        point: p,
        curve: sp,
    });
    let sk = b.build();
    for target in [
        DVec2::new(1.0, 1.5),
        DVec2::new(2.2, -1.0),
        DVec2::new(3.5, 0.9),
        DVec2::new(-0.5, -0.3),
    ] {
        let out = solve(&sk, &drag(&[(p, target)]));
        assert_eq!(out.status, SolveStatus::Solved, "{target}");
        assert_satisfied(&out.sketch);
        let got = pos(&out.sketch, p);
        let nearest = common::dist_to_spline(&out.sketch, sp, target);
        assert!(
            ((got - target).length() - nearest).abs() <= 1e-9,
            "{target}: {} vs {nearest}",
            (got - target).length()
        );
    }
}

#[test]
fn several_points_dragged_at_once() {
    let mut b = B::new();
    let (l, s, e) = b.line(0.0, 0.0, 1.0, 0.0);
    b.c(Constraint::Horizontal { line: l });
    let sk = b.build();
    // Both ends dragged to targets at different heights: the horizontal
    // constraint holds and the ends split the difference.
    let out = solve(
        &sk,
        &drag(&[(s, DVec2::new(0.0, 1.0)), (e, DVec2::new(2.0, 2.0))]),
    );
    assert_eq!(out.status, SolveStatus::Solved);
    assert_satisfied(&out.sketch);
    assert!(close(pos(&out.sketch, s), DVec2::new(0.0, 1.5), 1e-14));
    assert!(close(pos(&out.sketch, e), DVec2::new(2.0, 1.5), 1e-14));
}

#[test]
fn bad_drags_are_invalid_and_change_nothing() {
    let mut b = B::new();
    let (l, s, e) = b.line(0.0, 0.0, 1.0, 1.0);
    let sk = b.build();
    let ok = DVec2::new(1.0, 2.0);
    let cases: Vec<(Vec<(u32, DVec2)>, &str)> = vec![
        (vec![(l, ok)], "not a point"),
        (vec![(999, ok)], "does not exist"),
        (vec![(s, ok), (s, ok)], "twice"),
        (vec![(s, DVec2::new(f64::NAN, 0.0))], "finite"),
        (vec![(e, DVec2::new(0.0, f64::INFINITY))], "finite"),
        (vec![(e, DVec2::new(2e6, 0.0))], "finite"),
    ];
    for (d, what) in cases {
        let out = solve(&sk, &drag(&d));
        match &out.status {
            SolveStatus::Failed(SolveError::Invalid(msg)) => {
                assert!(msg.contains(what), "{d:?}: {msg}")
            }
            other => panic!("{d:?}: {other:?}"),
        }
        assert_eq!(out.sketch, sk);
        assert_eq!(bits(&out.sketch), bits(&sk));
        assert_eq!(out.dof, 4);
        assert!(out.entities.values().all(|&st| st == EntityState::Free));
    }
}
