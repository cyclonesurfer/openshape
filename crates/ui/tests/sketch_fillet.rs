//! Sketch fillet (design note 0007 §8), driven through the editor: exact
//! tangent points and centres, degrees of freedom unchanged, the corner's
//! constraints carried by a virtual sharp, refusals, and a fuzz of random
//! fillets that must never break a sketch.

mod common;

use common::*;
use opendatum_app::sketch_mode::{SketchEditor, Snap, SnapTarget, Tool};
use opendatum_app::sketch_ops::fillet::fillet;
use opendatum_app::sketch_ops::{shapes, Batch};
use opendatum_app::Model;
use opendatum_core_model::sketch::{
    Constraint, ConstraintEntry, DistanceKind, Geometry, TangentKind,
};
use opendatum_core_model::{DVec2, FeatureId, PlaneRef, SketchEdit, SketchId};
use opendatum_sketch_solver::{ResolvedCurve, SolveStatus};
use proptest::prelude::*;
use std::f64::consts::PI;

fn pos(m: &Model, f: FeatureId, id: SketchId) -> DVec2 {
    m.sketch_result(f)
        .unwrap()
        .solution
        .sketch
        .point(id)
        .expect("a point")
}

/// The sketch's arcs as (id, centre, radius).
fn arcs(m: &Model, f: FeatureId) -> Vec<(SketchId, DVec2, f64)> {
    curves(m, f, true)
        .into_iter()
        .filter_map(|(id, c)| match c {
            ResolvedCurve::Arc { center, radius, .. } => Some((id, center, radius)),
            _ => None,
        })
        .collect()
}

fn redundant(m: &Model, f: FeatureId) -> usize {
    m.sketch_result(f).unwrap().solution.redundant.len()
}

/// Fillets the corner at a plane point (millimetres) with radius `r` mm.
fn round(m: &mut Model, ed: &mut SketchEditor, x: f64, y: f64, r: f64) -> Option<String> {
    let map = view();
    ed.set_tool(Tool::SketchFillet);
    ed.fillet_radius = r * MM;
    click(m, ed, &map, x, y);
    ed.set_tool(Tool::Select);
    ed.message.clone()
}

fn near(p: DVec2, x: f64, y: f64) -> bool {
    (p - DVec2::new(x, y) * MM).length() < 1e-6
}

#[test]
fn a_rectangle_corner_gets_a_tangent_arc() {
    let (mut m, mut ed, map) = setup();
    let f = ed.feature;
    draw(
        &mut m,
        &mut ed,
        &map,
        Tool::Rectangle,
        &[(-20.0, -10.0), (20.0, 10.0)],
    );
    let before = dof(&m, f);
    // The hover preview is the fillet arc.
    ed.set_tool(Tool::SketchFillet);
    ed.fillet_radius = 5.0 * MM;
    hover(&m, &mut ed, &map, 20.0, 10.0);
    match ed
        .curve_edit
        .clone()
        .expect("a corner")
        .plan
        .expect("a plan")
        .preview
    {
        opendatum_app::sketch_ops::trim::Preview::Piece(ResolvedCurve::Arc {
            center,
            radius,
            sweep_angle,
            ..
        }) => {
            assert!(near(center, 15.0, 5.0), "{center}");
            assert!((radius - 5.0 * MM).abs() < 1e-12);
            assert!((sweep_angle - PI / 2.0).abs() < 1e-9);
        }
        p => panic!("{p:?}"),
    }
    let message = round(&mut m, &mut ed, 20.0, 10.0, 5.0);
    assert_eq!(message.as_deref(), Some("Filleted"));
    let a = arcs(&m, f);
    assert_eq!(a.len(), 1);
    assert!(near(a[0].1, 15.0, 5.0));
    assert!(has_radius(&m, f, a[0].0, 0.005));
    assert_eq!(
        constraints(&m, f, |c| matches!(
            c,
            Constraint::Tangent {
                kind: TangentKind::AtPoint { .. },
                ..
            }
        )),
        2
    );
    assert_eq!(status(&m, f), SolveStatus::Solved);
    assert_eq!(redundant(&m, f), 0);
    assert_eq!(dof(&m, f), before);
    // The rounded rectangle: 800 mm² less the corner outside the arc.
    let areas = region_areas(&m, f);
    assert_eq!(areas.len(), 1);
    let expected = (800.0 - 25.0 + PI * 25.0 / 4.0) * MM * MM;
    assert!((areas[0] - expected).abs() < 1e-12, "{areas:?}");
    // Every corner.
    for (x, y) in [(-20.0, 10.0), (-20.0, -10.0), (20.0, -10.0)] {
        round(&mut m, &mut ed, x, y, 5.0);
    }
    assert_eq!(arcs(&m, f).len(), 4);
    assert_eq!(dof(&m, f), before);
    let expected = (800.0 - 4.0 * (25.0 - PI * 25.0 / 4.0)) * MM * MM;
    assert!((region_areas(&m, f)[0] - expected).abs() < 1e-12);
    let x = extrude(&mut m, &map, &[(0.0, 0.0)], "10 mm");
    assert!((volume(&m, x) - expected * 0.01).abs() < 1e-14);
}

fn has_radius(m: &Model, f: FeatureId, arc: SketchId, value: f64) -> bool {
    constraints(
        m,
        f,
        |c| matches!(c, Constraint::Radius { curve, value: v } if *curve == arc && (*v - value).abs() < 1e-15),
    ) == 1
}

#[test]
fn a_centre_rectangle_keeps_its_diagonals_on_the_virtual_sharps() {
    let (mut m, mut ed, map) = setup();
    let f = ed.feature;
    draw(
        &mut m,
        &mut ed,
        &map,
        Tool::CenterRectangle,
        &[(0.0, 0.0), (20.0, 10.0)],
    );
    let before = dof(&m, f);
    for (x, y) in [(20.0, 10.0), (-20.0, 10.0), (-20.0, -10.0), (20.0, -10.0)] {
        assert_eq!(
            round(&mut m, &mut ed, x, y, 4.0).as_deref(),
            Some("Filleted")
        );
    }
    assert_eq!(arcs(&m, f).len(), 4);
    assert_eq!(status(&m, f), SolveStatus::Solved);
    assert_eq!(redundant(&m, f), 0);
    assert_eq!(dof(&m, f), before);
    // The diagonals still run corner to (virtual) corner.
    let diagonals: Vec<(DVec2, DVec2)> = m
        .sketch_result(f)
        .unwrap()
        .curves
        .iter()
        .filter_map(|(_, construction, c)| match c {
            ResolvedCurve::Line { start, end } if *construction => Some((*start, *end)),
            _ => None,
        })
        .collect();
    assert_eq!(diagonals.len(), 2);
    for (s, e) in diagonals {
        assert!((s.abs() - DVec2::new(0.02, 0.01)).length() < 1e-9, "{s}");
        assert!((e.abs() - DVec2::new(0.02, 0.01)).length() < 1e-9, "{e}");
    }
    let expected = (800.0 - 4.0 * (16.0 - PI * 16.0 / 4.0)) * MM * MM;
    assert!((region_areas(&m, f)[0] - expected).abs() < 1e-12);
}

#[test]
fn a_dimensioned_corner_keeps_its_size_through_a_virtual_sharp() {
    let (mut m, mut ed, map) = setup();
    let f = ed.feature;
    draw(
        &mut m,
        &mut ed,
        &map,
        Tool::Rectangle,
        &[(-20.0, -10.0), (20.0, 10.0)],
    );
    let top = curves(&m, f, false)
        .into_iter()
        .find(|(_, c)| matches!(c, ResolvedCurve::Line { start, end } if (start.y - 0.01).abs() < 1e-9 && (end.y - 0.01).abs() < 1e-9))
        .unwrap()
        .0;
    m.edit_sketch(
        f,
        vec![SketchEdit::AddConstraint(ConstraintEntry {
            id: SketchId::new_in_batch(0),
            constraint: Constraint::Length {
                line: top,
                value: 0.04,
            },
        })],
    )
    .unwrap();
    let before = dof(&m, f);
    round(&mut m, &mut ed, 20.0, 10.0, 5.0);
    assert_eq!(
        constraints(&m, f, |c| matches!(c, Constraint::Length { .. })),
        0
    );
    assert_eq!(
        constraints(
            &m,
            f,
            |c| matches!(c, Constraint::Distance { kind: DistanceKind::Aligned, value, .. } if *value == 0.04)
        ),
        1
    );
    assert_eq!(status(&m, f), SolveStatus::Solved);
    assert_eq!(redundant(&m, f), 0);
    assert_eq!(dof(&m, f), before);
    // A larger radius leaves the 40 mm side where it was.
    let radius = m
        .sketch_def(f)
        .unwrap()
        .sketch
        .constraints
        .iter()
        .find(|c| matches!(c.constraint, Constraint::Radius { .. }))
        .unwrap()
        .id;
    m.edit_sketch(
        f,
        vec![SketchEdit::SetValue {
            constraint: radius,
            value: 0.008,
        }],
    )
    .unwrap();
    assert_eq!(status(&m, f), SolveStatus::Solved);
    assert!((arcs(&m, f)[0].2 - 0.008).abs() < 1e-12);
    // The sharp is still where the two sides would meet, 40 mm from the
    // top side's far end.
    let (sharp, far) = m
        .sketch_def(f)
        .unwrap()
        .sketch
        .constraints
        .iter()
        .find_map(|c| match c.constraint {
            Constraint::Distance { a, b, .. } => Some((a, b)),
            _ => None,
        })
        .unwrap();
    assert!(((pos(&m, f, sharp) - pos(&m, f, far)).length() - 0.04).abs() < 1e-12);
    let lines: Vec<(DVec2, DVec2)> = curves(&m, f, false)
        .into_iter()
        .filter_map(|(_, c)| match c {
            ResolvedCurve::Line { start, end } => Some((start, end)),
            _ => None,
        })
        .collect();
    let on = |p: DVec2, (a, b): (DVec2, DVec2)| (b - a).normalize().perp_dot(p - a).abs() < 1e-12;
    let s = pos(&m, f, sharp);
    assert_eq!(
        lines.iter().filter(|l| on(s, **l)).count(),
        2,
        "on the two filleted sides"
    );
}

#[test]
fn a_line_and_an_arc_are_filleted_on_the_inside() {
    let (mut m, mut ed, map) = setup();
    let f = ed.feature;
    draw(
        &mut m,
        &mut ed,
        &map,
        Tool::Line,
        &[(-30.0, 0.0), (0.0, 0.0)],
    );
    ed.set_tool(Tool::CenterArc);
    click(&mut m, &mut ed, &map, 20.0, 0.0);
    click(&mut m, &mut ed, &map, 0.0, 0.0);
    hover(&m, &mut ed, &map, 5.86, -14.14);
    click(&mut m, &mut ed, &map, 20.0, -20.0);
    ed.set_tool(Tool::Select);
    assert_eq!(arcs(&m, f).len(), 1);
    let before = dof(&m, f);
    round(&mut m, &mut ed, 0.0, 0.0, 5.0);
    let a = arcs(&m, f);
    assert_eq!(a.len(), 2);
    let fillet = a
        .iter()
        .find(|x| (x.2 - 0.005).abs() < 1e-9)
        .expect("the fillet");
    let cx = 20.0 - (625.0f64 - 25.0).sqrt();
    assert!(near(fillet.1, cx, -5.0), "{}", fillet.1);
    assert_eq!(status(&m, f), SolveStatus::Solved);
    assert_eq!(redundant(&m, f), 0);
    assert_eq!(dof(&m, f), before);
}

#[test]
fn two_arcs_at_a_corner_are_filleted_outside_both() {
    let (mut m, mut ed, _) = setup();
    let f = ed.feature;
    // Arcs about (−20, 0) and (0, 20) meeting at right angles at the origin.
    add_cornered(
        &mut m,
        f,
        &Cornered::Joined(
            DVec2::ZERO,
            Piece::Arc(DVec2::new(-20.0, 0.0), -1.0),
            Piece::Arc(DVec2::new(0.0, 20.0), 1.0),
        ),
    );
    let before = dof(&m, f);
    assert_eq!(
        round(&mut m, &mut ed, 0.0, 0.0, 5.0).as_deref(),
        Some("Filleted")
    );
    let t = (-20.0 + 850.0f64.sqrt()) / 2.0;
    let fillet = arcs(&m, f)
        .into_iter()
        .find(|a| (a.2 - 0.005).abs() < 1e-9)
        .expect("the fillet");
    assert!(near(fillet.1, t, -t), "{}", fillet.1);
    assert_eq!(status(&m, f), SolveStatus::Solved);
    assert_eq!(redundant(&m, f), 0);
    assert_eq!(dof(&m, f), before);
}

#[test]
fn fillets_that_cannot_be_made_are_refused() {
    let (mut m, mut ed, map) = setup();
    let f = ed.feature;
    // A straight joint.
    draw(
        &mut m,
        &mut ed,
        &map,
        Tool::Line,
        &[(-40.0, 30.0), (0.0, 30.0), (40.0, 30.0)],
    );
    let message = round(&mut m, &mut ed, 0.0, 30.0, 5.0).unwrap();
    assert!(message.contains("smoothly"), "{message}");
    // Too large a radius for a 40 × 20 rectangle.
    draw(
        &mut m,
        &mut ed,
        &map,
        Tool::Rectangle,
        &[(-20.0, -10.0), (20.0, 10.0)],
    );
    let message = round(&mut m, &mut ed, 20.0, 10.0, 50.0).unwrap();
    assert!(message.contains("too large"), "{message}");
    // Three lines at a point.
    draw(
        &mut m,
        &mut ed,
        &map,
        Tool::Line,
        &[(20.0, 10.0), (40.0, 0.0)],
    );
    let message = round(&mut m, &mut ed, 20.0, 10.0, 2.0).unwrap();
    assert!(message.contains("exactly two"), "{message}");
    // Not a point.
    let message = round(&mut m, &mut ed, 0.0, -40.0, 2.0).unwrap();
    assert!(message.contains("corner"), "{message}");
    assert_eq!(arcs(&m, f).len(), 0);
}

/// A random closed or open shape with corners (millimetres).
#[derive(Clone, Debug)]
enum Cornered {
    Rectangle(DVec2, DVec2),
    Polygon(DVec2, DVec2, u32),
    Chain(Vec<DVec2>),
    /// Two pieces meeting at a corner.
    Joined(DVec2, Piece, Piece),
}

/// A line or arc leaving a corner.
#[derive(Clone, Debug)]
enum Piece {
    /// To its far end.
    Line(DVec2),
    /// About a centre, turning by an angle (counter-clockwise positive).
    Arc(DVec2, f64),
}

fn cornered() -> impl Strategy<Value = Cornered> {
    let p = || (-50.0..50.0f64, -50.0..50.0f64).prop_map(|(x, y)| DVec2::new(x, y));
    prop_oneof![
        (p(), p()).prop_map(|(c, k)| Cornered::Rectangle(c, k)),
        (p(), p(), 3u32..9).prop_map(|(c, v, n)| Cornered::Polygon(c, v, n)),
        proptest::collection::vec(p(), 3..6).prop_map(Cornered::Chain),
        (p(), piece(), piece()).prop_map(|(v, a, b)| Cornered::Joined(v, a, b)),
    ]
}

fn piece() -> impl Strategy<Value = Piece> {
    let p = || (-50.0..50.0f64, -50.0..50.0f64).prop_map(|(x, y)| DVec2::new(x, y));
    prop_oneof![
        p().prop_map(Piece::Line),
        (p(), prop_oneof![-5.5..-0.3f64, 0.3..5.5f64]).prop_map(|(c, t)| Piece::Arc(c, t)),
    ]
}

fn add_cornered(m: &mut Model, f: FeatureId, s: &Cornered) {
    let free = |p: DVec2| Snap {
        at: p * MM,
        target: SnapTarget::Free,
    };
    let mut b = Batch::default();
    let made = match s {
        Cornered::Rectangle(c, k) => shapes::centre_rectangle(&mut b, &free(*c), &free(*k), false),
        Cornered::Polygon(c, v, n) => shapes::polygon(&mut b, &free(*c), &free(*v), *n, false),
        Cornered::Chain(points) => {
            if points.windows(2).any(|w| (w[0] - w[1]).length() < 2.0) {
                return;
            }
            let mut last: Option<SketchId> = None;
            for w in points.windows(2) {
                let l = b.line(w[0] * MM, w[1] * MM, false);
                if let Some(e) = last {
                    b.coincident(e, l.start);
                }
                last = Some(l.end);
            }
            Ok(())
        }
        Cornered::Joined(v, a, e) => {
            let mut corner = Vec::new();
            for piece in [a, e] {
                match piece {
                    Piece::Line(far) => {
                        if (*far - *v).length() < 2.0 {
                            return;
                        }
                        corner.push(b.line(*v * MM, *far * MM, false).start);
                    }
                    Piece::Arc(c, turn) => {
                        if (*c - *v).length() < 2.0 {
                            return;
                        }
                        let other = *c + DVec2::from_angle(*turn).rotate(*v - *c);
                        let arc = if *turn > 0.0 {
                            b.arc(*c * MM, *v * MM, other * MM, false)
                        } else {
                            b.arc(*c * MM, other * MM, *v * MM, false)
                        };
                        corner.push(if *turn > 0.0 { arc.start } else { arc.end });
                    }
                }
            }
            b.coincident(corner[0], corner[1]);
            Ok(())
        }
    };
    if made.is_ok() {
        m.edit_sketch(f, b.edits).unwrap();
    }
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(256))]

    /// Random fillets on random cornered shapes: each one either is refused
    /// with a reason or commits through the checked path with the sketch
    /// solved, no redundancy, degrees of freedom unchanged (but for the
    /// equal-length and midpoint constraints a fillet removes), and
    /// nothing but the two filleted ends moving.
    #[test]
    fn random_fillets_keep_the_sketch_solved(
        shapes in proptest::collection::vec(cornered(), 1..4),
        steps in proptest::collection::vec((0usize..64, 0.2..15.0f64), 1..6),
    ) {
        let mut m = Model::new();
        let f = m.add_sketch(PlaneRef::Top);
        for s in &shapes {
            add_cornered(&mut m, f, s);
        }
        prop_assume!(status(&m, f) == SolveStatus::Solved);
        prop_assume!(redundant(&m, f) == 0);
        for (which, r) in steps {
            let res = m.sketch_result(f).unwrap();
            let curves: Vec<(SketchId, ResolvedCurve)> =
                res.curves.iter().map(|(id, _, c)| (*id, c.clone())).collect();
            let sketch = res.solution.sketch.clone();
            let points: Vec<SketchId> = sketch
                .entities
                .iter()
                .filter(|e| matches!(e.geometry, Geometry::Point { .. }))
                .map(|e| e.id)
                .collect();
            if points.is_empty() {
                break;
            }
            let corner = points[which % points.len()];
            let before = dof(&m, f);
            let measures = |m: &Model| {
                (
                    constraints(m, f, |c| matches!(c, Constraint::Equal { .. })),
                    constraints(m, f, |c| matches!(c, Constraint::Midpoint { .. })),
                )
            };
            let (equal, mid) = measures(&m);
            let Ok(plan) = fillet(&sketch, &curves, corner, r * MM) else { continue };
            let moved: Vec<(SketchId, DVec2)> = plan
                .edits
                .iter()
                .flat_map(|e| match e {
                    SketchEdit::MovePoints { points } => points.clone(),
                    _ => Vec::new(),
                })
                .collect();
            if let Err(e) = m.edit_sketch_checked(f, plan.edits.clone()) {
                return Err(TestCaseError::fail(format!("fillet r={r} at {corner:?}: {e}")));
            }
            prop_assert_eq!(status(&m, f), SolveStatus::Solved);
            prop_assert_eq!(redundant(&m, f), 0);
            // Only removed equal lengths (one equation each) and midpoints
            // (two) free anything.
            let (equal_after, mid_after) = measures(&m);
            let (de, dm) = (equal - equal_after, mid - mid_after);
            prop_assert_eq!(plan.removed, de + dm);
            prop_assert_eq!(dof(&m, f), before + de + 2 * dm);
            let after = m.sketch_result(f).unwrap();
            for e in &sketch.entities {
                let Geometry::Point { at } = e.geometry else { continue };
                let Some(now) = after.solution.sketch.point(e.id) else { continue };
                let to = moved.iter().find(|(id, _)| *id == e.id).map_or(at, |(_, p)| *p);
                prop_assert!((now - to).length() < 1e-9, "{:?} moved from {} to {}", e.id, to, now);
            }
        }
    }
}
