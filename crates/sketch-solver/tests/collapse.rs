//! Collapsed curves (design note 0005 §9, final acceptance check): a solve
//! never succeeds with a zero-length line, a zero-radius circle or arc, or
//! an arc that closes on itself. When one constraint causes the collapse it
//! is reported as conflicting (the latest of those holding the curve
//! collapsed whose removal avoids it), the geometry is left at the input,
//! and a curve already degenerate in the input is reported `Degenerate`.
//! Only constraints that hold the collapse are tried, so blaming stays
//! cheap in large sketches.

mod common;

use common::{bits, circ, line, pos, B};
use openshape_sketch_solver::{
    solve, Constraint, DVec2, EntityState, Geometry, Sketch, Solution, SolveError, SolveOptions,
    SolveStatus, TangentKind,
};
use std::f64::consts::FRAC_PI_2;

fn solve_default(s: &Sketch<u32>) -> Solution<u32> {
    solve(s, &SolveOptions::default())
}

/// Expects `Conflicting` blaming exactly `culprits`, with the input back.
fn expect_blamed(s: &Sketch<u32>, culprits: &[u32]) -> Solution<u32> {
    let out = solve_default(s);
    assert_eq!(
        out.status,
        SolveStatus::Failed(SolveError::Conflicting),
        "{out:?}"
    );
    assert_eq!(out.conflicting, culprits);
    assert!(out.redundant.iter().all(|k| !culprits.contains(k)));
    assert_eq!(bits(&out.sketch), bits(s), "geometry left at the input");
    assert!(out.max_residual.is_finite());
    // Without the culprits the rest solves, with nothing collapsed.
    let mut rest = s.clone();
    rest.constraints.retain(|c| !culprits.contains(&c.id));
    let again = solve_default(&rest);
    assert!(
        matches!(again.status, SolveStatus::Solved | SolveStatus::Redundant),
        "{:?}",
        again.status
    );
    out
}

/// Expects `Degenerate` naming `what`, reported like an early failure.
fn expect_degenerate(s: &Sketch<u32>, what: &str) {
    let out = solve_default(s);
    match &out.status {
        SolveStatus::Failed(SolveError::Degenerate(msg)) => {
            assert!(msg.contains(what), "expected {what:?} in {msg:?}")
        }
        other => panic!("expected Degenerate ({what}), got {other:?}"),
    }
    assert_eq!(bits(&out.sketch), bits(s));
    assert!(out.entities.values().all(|&e| e == EntityState::Free));
    assert!(out.redundant.is_empty() && out.conflicting.is_empty());
    assert_eq!((out.max_residual, out.iterations), (0.0, 0));
}

#[test]
fn horizontal_then_vertical_blames_the_vertical() {
    let mut b = B::new();
    let (l, ..) = b.line(0.0, 0.02, 0.04, 0.02);
    b.c(Constraint::Horizontal { line: l });
    let v = b.c(Constraint::Vertical { line: l });
    let out = expect_blamed(&b.build(), &[v]);
    // The diagnosis is that of the state reached: one point's worth free.
    assert_eq!(out.dof, 2);

    // On a slanted line either constraint alone is fine, and the later one
    // is blamed.
    let mut b = B::new();
    let (l, ..) = b.line(0.0, 0.0, 0.04, 0.03);
    b.c(Constraint::Vertical { line: l });
    let h = b.c(Constraint::Horizontal { line: l });
    expect_blamed(&b.build(), &[h]);

    // On a horizontal line, vertical alone collapses it (the nearest
    // vertical line has no length), so it is blamed even when earlier.
    let mut b = B::new();
    let (l, ..) = b.line(0.0, 0.02, 0.04, 0.02);
    let v = b.c(Constraint::Vertical { line: l });
    b.c(Constraint::Horizontal { line: l });
    expect_blamed(&b.build(), &[v]);
}

#[test]
fn coincident_end_points_of_one_line_are_blamed() {
    let mut b = B::new();
    let (_, p, q) = b.line(0.0, 0.02, 0.04, 0.02);
    let c = b.coincident(p, q);
    expect_blamed(&b.build(), &[c]);
}

#[test]
fn horizontal_and_vertical_points_blame_the_later() {
    let mut b = B::new();
    let (_, p, q) = b.line(0.0, 0.02, 0.04, 0.02);
    b.c(Constraint::HorizontalPoints { a: p, b: q });
    let v = b.c(Constraint::VerticalPoints { a: p, b: q });
    expect_blamed(&b.build(), &[v]);
}

#[test]
fn the_latest_constraint_whose_removal_helps_is_blamed() {
    // A fixed rectangle, then a line collapsed by a vertical constraint,
    // then a later fix that does not cause the collapse.
    let mut b = B::new();
    let (bottom, p0, p1) = b.line(0.0, 0.0, 0.04, 0.0);
    let (right, p2, p3) = b.line(0.04, 0.0, 0.04, 0.03);
    let (top, p4, p5) = b.line(0.04, 0.03, 0.0, 0.03);
    let (left, p6, p7) = b.line(0.0, 0.03, 0.0, 0.0);
    for (x, y) in [(p1, p2), (p3, p4), (p5, p6), (p7, p0)] {
        b.coincident(x, y);
    }
    b.c(Constraint::Horizontal { line: bottom });
    b.c(Constraint::Horizontal { line: top });
    b.c(Constraint::Vertical { line: right });
    b.c(Constraint::Vertical { line: left });
    b.fix(p0);
    // A loose line, horizontal, then vertical, then fixed at one end.
    let (l, s, _) = b.line(0.1, 0.02, 0.14, 0.02);
    b.c(Constraint::Horizontal { line: l });
    let v = b.c(Constraint::Vertical { line: l });
    b.fix(s);
    let out = expect_blamed(&b.build(), &[v]);
    assert_eq!(out.entities[&bottom], EntityState::Free);

    // Collapsing a rectangle side: the coincident of its ends is blamed.
    let mut b2 = b.clone();
    b2.sketch.constraints.retain(|c| c.id != v);
    let c = b2.coincident(p2, p3);
    expect_blamed(&b2.build(), &[c]);
}

#[test]
fn each_collapsed_part_gets_its_culprit() {
    let mut b = B::new();
    let (l, ..) = b.line(0.0, 0.02, 0.04, 0.02);
    let (_, p, q) = b.line(0.1, 0.0, 0.1, 0.03);
    b.c(Constraint::Horizontal { line: l });
    let c = b.coincident(p, q);
    let v = b.c(Constraint::Vertical { line: l });
    expect_blamed(&b.build(), &[c, v]);
}

#[test]
fn an_arc_closed_on_itself_is_blamed() {
    let mut b = B::new();
    let (_, _, s, e) = b.arc(DVec2::new(0.01, 0.01), 0.02, 0.0, FRAC_PI_2);
    let c = b.coincident(s, e);
    expect_blamed(&b.build(), &[c]);
}

#[test]
fn an_arc_with_its_start_on_its_centre_is_degenerate() {
    // The elimination of the coincidence already leaves no radius.
    let mut b = B::new();
    let (a, cc, s, _) = b.arc(DVec2::new(0.01, 0.01), 0.02, 0.0, FRAC_PI_2);
    b.coincident(s, cc);
    expect_degenerate(&b.build(), &format!("arc {a:?} has zero radius"));
}

#[test]
fn a_circle_pulled_to_zero_radius_is_blamed() {
    // Tangent to a line with its centre on the line: the radius must be 0.
    // It cannot reach 0, so the solve fails heading there; the collapse
    // is blamed on the later constraint.
    let mut b = B::new();
    let (l, ..) = b.line(-0.05, 0.0, 0.05, 0.0);
    let (k, cc) = b.circle(0.0, 0.02, 0.01);
    b.c(Constraint::Tangent {
        a: l,
        b: k,
        kind: TangentKind::LineCircle { center_left: true },
    });
    let on = b.c(Constraint::PointOnCurve {
        point: cc,
        curve: l,
    });
    let out = expect_blamed(&b.build(), &[on]);
    assert_eq!(circ(&out.sketch, k).1, 0.01);

    // Concentric and touching from outside: both radii must be 0.
    let mut b = B::new();
    let (k1, c1) = b.circle(0.0, 0.0, 0.01);
    let (k2, c2) = b.circle(0.03, 0.0, 0.01);
    b.c(Constraint::Tangent {
        a: k1,
        b: k2,
        kind: TangentKind::CircleCircle { internal: false },
    });
    let c = b.coincident(c1, c2);
    expect_blamed(&b.build(), &[c]);
}

#[test]
fn a_degenerate_input_line_never_solves() {
    // No constraints at all.
    let mut b = B::new();
    let (l, ..) = b.line(0.01, 0.02, 0.01, 0.02);
    b.line(0.0, 0.0, 0.04, 0.0);
    expect_degenerate(&b.build(), &format!("line {l:?} has zero length"));

    // End points coincident-constrained and already together.
    let mut b = B::new();
    let (l, p, q) = b.line(0.01, 0.02, 0.01, 0.02);
    b.coincident(p, q);
    expect_degenerate(&b.build(), &format!("line {l:?} has zero length"));

    // End points coincident-constrained and together within the model
    // tolerance: the input satisfies the constraint but is no solution.
    let mut b = B::new();
    let (_, p, q) = b.line(0.01, 0.02, 0.01 + 5e-8, 0.02);
    let c = b.coincident(p, q);
    expect_blamed(&b.build(), &[c]);

    // Pulled apart by later constraints: a valid result.
    let mut b = B::new();
    let (l, p, q) = b.line(0.01, 0.02, 0.01, 0.02);
    b.c(Constraint::Fix {
        point: p,
        at: DVec2::new(0.0, 0.0),
    });
    b.c(Constraint::Fix {
        point: q,
        at: DVec2::new(0.03, 0.0),
    });
    let out = solve_default(&b.build());
    assert_eq!(out.status, SolveStatus::Solved);
    let (a, z) = line(&out.sketch, l);
    assert_eq!((a, z), (DVec2::ZERO, DVec2::new(0.03, 0.0)));
}

#[test]
fn a_drag_onto_the_other_end_fails() {
    let mut b = B::new();
    let (l, p, q) = b.line(0.0, 0.0, 0.04, 0.0);
    b.c(Constraint::Horizontal { line: l });
    let s = b.build();
    let options = SolveOptions {
        drag: vec![(q, pos(&s, p))],
        ..SolveOptions::default()
    };
    let out = solve(&s, &options);
    assert!(
        matches!(out.status, SolveStatus::Failed(SolveError::Degenerate(_))),
        "{:?}",
        out.status
    );
    assert_eq!(bits(&out.sketch), bits(&s));
    // Short of it, the drag is followed.
    let options = SolveOptions {
        drag: vec![(q, DVec2::new(1e-5, 0.0))],
        ..SolveOptions::default()
    };
    let out = solve(&s, &options);
    assert_eq!(out.status, SolveStatus::Solved);
    assert!((pos(&out.sketch, q).x - 1e-5).abs() < 1e-12);
}

#[test]
fn very_short_lines_still_solve() {
    for len in [1e-5, 1e-6, 1e-7] {
        // A 1 m sketch with a dimensioned short line hanging off a corner.
        let mut b = B::new();
        let (base, p0, p1) = b.line(0.0, 0.0, 1.0, 0.0);
        let (short, s, e) = b.line(1.0, 0.0, 1.0 + 2.0 * len, 0.1 * len);
        b.coincident(p1, s);
        b.c(Constraint::Horizontal { line: base });
        b.c(Constraint::Horizontal { line: short });
        b.c(Constraint::Length {
            line: short,
            value: len,
        });
        b.fix(p0);
        let out = solve_default(&b.build());
        assert_eq!(out.status, SolveStatus::Solved, "{len}");
        let (a, z) = line(&out.sketch, short);
        assert!(((z - a).length() - len).abs() <= 1e-12, "{len}");
        assert!(pos(&out.sketch, e).x > pos(&out.sketch, s).x);
        // And again from the result.
        let again = solve_default(&out.sketch);
        assert_eq!(again.status, SolveStatus::Solved, "{len}");
    }
    // A tiny circle and arc.
    let mut b = B::new();
    let (k, _) = b.circle(0.5, 0.5, 1e-6);
    let (a, ..) = b.arc(DVec2::new(0.0, 1.0), 1e-6, 0.0, 1.0);
    b.line(0.0, 0.0, 1.0, 0.0);
    b.c(Constraint::Radius {
        curve: k,
        value: 2e-6,
    });
    b.c(Constraint::Equal { a, b: k });
    let out = solve_default(&b.build());
    assert_eq!(out.status, SolveStatus::Solved);
    assert!((circ(&out.sketch, a).1 - 2e-6).abs() < 1e-12);
    assert!(matches!(
        out.sketch.entity(k).unwrap().geometry,
        Geometry::Circle { radius, .. } if radius == 2e-6
    ));
}

/// A chain of `n` lines from point `start` (at `at`), each with a length
/// (2 % longer than drawn, so that every solve of it iterates) and an angle
/// to the previous one.
fn chain(b: &mut B, n: usize, start: u32, mut at: DVec2) {
    let mut heading = 0.0f64;
    let mut prev: Option<(u32, u32)> = None;
    for i in 0..n {
        let len = 0.01 + 0.005 * ((i * 7 % 5) as f64);
        let turn = if i == 0 {
            0.0
        } else {
            0.3 + 0.3 * (i % 4) as f64
        };
        let turn = if i % 2 == 0 { turn } else { -turn };
        heading += turn;
        let end = at + DVec2::new(heading.cos(), heading.sin()) * len;
        let (l, s, e) = b.line(at.x, at.y, end.x, end.y);
        b.coincident(prev.map_or(start, |(_, pe)| pe), s);
        if let Some((pl, _)) = prev {
            b.c(Constraint::Angle {
                a: pl,
                b: l,
                value: turn,
            });
        }
        b.c(Constraint::Length {
            line: l,
            value: 1.02 * len,
        });
        prev = Some((l, e));
        at = end;
    }
}

#[test]
fn blaming_tries_only_constraints_holding_the_collapse() {
    // An early horizontal and vertical on a line, then a 40-line chain
    // hanging off it: only the line's own constraints hold the collapse,
    // so one trial solve finds the vertical (trying the chain's 120
    // constraints first would take hundreds of iterations).
    let mut b = B::new();
    let (l, s, e) = b.line(0.0, 0.0, 0.01, 0.0);
    b.fix(s);
    b.c(Constraint::Horizontal { line: l });
    let v = b.c(Constraint::Vertical { line: l });
    chain(&mut b, 40, e, DVec2::new(0.01, 0.0));
    let sketch = b.build();
    let out = expect_blamed(&sketch, &[v]);
    let mut without = sketch.clone();
    without.constraints.retain(|c| c.id != v);
    let alone = solve_default(&without);
    assert!(
        out.iterations <= 2 * alone.iterations + 10,
        "{} iterations",
        out.iterations
    );

    // A zero-length line hanging off a chain from a fixed point: nothing
    // holds it collapsed, so nothing is tried.
    let mut b = B::new();
    let (l, s, _) = b.line(0.0, 0.0, 0.0, 0.0);
    b.fix(s);
    chain(&mut b, 40, s, DVec2::ZERO);
    expect_degenerate(&b.build(), &format!("line {l:?} has zero length"));
}
