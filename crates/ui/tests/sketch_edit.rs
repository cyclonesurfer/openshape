//! Trim, extend and split (design note 0007 §7), driven through the editor
//! as the viewport drives it: where ends land, which constraints carry
//! over, degrees of freedom, no new redundancy, exact regions, and a fuzz of
//! random trims that must never break a sketch.

mod common;

use common::*;
use openshape_app::sketch_mode::{SketchEditor, Snap, SnapTarget, Tool};
use openshape_app::sketch_ops::trim::{self, Preview};
use openshape_app::sketch_ops::{shapes, Batch};
use openshape_app::Model;
use openshape_core_model::sketch::{Constraint, ConstraintEntry, Geometry, TangentKind};
use openshape_core_model::{DVec2, FeatureId, PlaneRef, SketchEdit, SketchId};
use openshape_sketch_solver::{ResolvedCurve, SolveStatus};
use proptest::prelude::*;
use std::f64::consts::{PI, TAU};

/// The sketch's solved position of point `id`.
fn pos(m: &Model, f: FeatureId, id: SketchId) -> DVec2 {
    m.sketch_result(f)
        .unwrap()
        .solution
        .sketch
        .point(id)
        .expect("a point")
}

/// Every line as (id, start point, end point).
fn line_ids(m: &Model, f: FeatureId) -> Vec<(SketchId, SketchId, SketchId)> {
    m.sketch_def(f)
        .unwrap()
        .sketch
        .entities
        .iter()
        .filter_map(|e| match e.geometry {
            Geometry::Line { start, end } => Some((e.id, start, end)),
            _ => None,
        })
        .collect()
}

/// Every arc as (id, centre, start, end).
fn arc_ids(m: &Model, f: FeatureId) -> Vec<(SketchId, SketchId, SketchId, SketchId)> {
    m.sketch_def(f)
        .unwrap()
        .sketch
        .entities
        .iter()
        .filter_map(|e| match e.geometry {
            Geometry::Arc { center, start, end } => Some((e.id, center, start, end)),
            _ => None,
        })
        .collect()
}

/// The line whose end points are near `a` and `b` (millimetres), either way
/// round, as (id, start, end).
fn line_near(
    m: &Model,
    f: FeatureId,
    a: (f64, f64),
    b: (f64, f64),
) -> (SketchId, SketchId, SketchId) {
    let (a, b) = (DVec2::new(a.0, a.1) * MM, DVec2::new(b.0, b.1) * MM);
    let near = |p: DVec2, q: DVec2| (p - q).length() < 1e-6;
    *line_ids(m, f)
        .iter()
        .find(|(_, s, e)| {
            let (s, e) = (pos(m, f, *s), pos(m, f, *e));
            (near(s, a) && near(e, b)) || (near(s, b) && near(e, a))
        })
        .unwrap_or_else(|| panic!("no line from {a} to {b}"))
}

fn has(m: &Model, f: FeatureId, c: &Constraint<SketchId>) -> bool {
    constraints(m, f, |x| x == c) == 1
}

fn redundant(m: &Model, f: FeatureId) -> usize {
    m.sketch_result(f).unwrap().solution.redundant.len()
}

/// Uses `tool` once at a plane point (millimetres), then returns to Select.
fn act(m: &mut Model, ed: &mut SketchEditor, tool: Tool, x: f64, y: f64) -> Option<String> {
    let map = view();
    ed.set_tool(tool);
    click(m, ed, &map, x, y);
    ed.set_tool(Tool::Select);
    ed.message.clone()
}

/// Adds one constraint to the sketch.
fn constrain(m: &mut Model, f: FeatureId, c: Constraint<SketchId>) {
    m.edit_sketch(
        f,
        vec![SketchEdit::AddConstraint(ConstraintEntry {
            id: SketchId::new_in_batch(0),
            constraint: c,
        })],
    )
    .unwrap();
}

/// A counter-clockwise centre-point arc about `c` from `s`, turning past
/// `via` to `e` (millimetres).
fn centre_arc(
    m: &mut Model,
    ed: &mut SketchEditor,
    c: (f64, f64),
    s: (f64, f64),
    via: (f64, f64),
    e: (f64, f64),
) {
    let map = view();
    ed.set_tool(Tool::CenterArc);
    click(m, ed, &map, c.0, c.1);
    click(m, ed, &map, s.0, s.1);
    hover(m, ed, &map, via.0, via.1);
    click(m, ed, &map, e.0, e.1);
    ed.set_tool(Tool::Select);
}

/// A horizontal line from (−40, 0) to (40, 0) crossed by vertical lines at
/// the given x (from y = −30 to 30).
fn crossed_line(xs: &[f64]) -> (Model, SketchEditor) {
    let (mut m, mut ed, map) = setup();
    draw(
        &mut m,
        &mut ed,
        &map,
        Tool::Line,
        &[(-40.0, 0.0), (40.0, 0.0)],
    );
    for x in xs {
        draw(
            &mut m,
            &mut ed,
            &map,
            Tool::Line,
            &[(*x, -30.0), (*x, 30.0)],
        );
    }
    (m, ed)
}

#[test]
fn trimming_an_end_piece_moves_the_end_onto_the_cutter() {
    let (mut m, mut ed) = crossed_line(&[0.0]);
    let f = ed.feature;
    let (a, _, a_end) = line_near(&m, f, (-40.0, 0.0), (40.0, 0.0));
    let (b, _, _) = line_near(&m, f, (0.0, -30.0), (0.0, 30.0));
    assert_eq!(dof(&m, f), 6);
    // The hover preview is the piece that goes.
    let map = view();
    ed.set_tool(Tool::Trim);
    hover(&m, &mut ed, &map, 30.0, 0.0);
    let edit = ed.curve_edit.clone().expect("a curve under the pointer");
    assert_eq!(edit.target, a);
    match edit.plan.expect("a plan").preview {
        Preview::Piece(ResolvedCurve::Line { start, end }) => {
            assert!((start - DVec2::ZERO).length() < 1e-9, "{start}");
            assert!((end - DVec2::new(40.0, 0.0) * MM).length() < 1e-6, "{end}");
        }
        p => panic!("{p:?}"),
    }
    let message = act(&mut m, &mut ed, Tool::Trim, 30.0, 0.0);
    assert_eq!(message.as_deref(), Some("Trimmed"));
    assert_eq!(count(&m, f, |g| matches!(g, Geometry::Line { .. })), 2);
    assert!(pos(&m, f, a_end).length() < 1e-9);
    assert!(has(
        &m,
        f,
        &Constraint::PointOnCurve {
            point: a_end,
            curve: b
        }
    ));
    assert!(has(&m, f, &Constraint::Horizontal { line: a }));
    assert_eq!(status(&m, f), SolveStatus::Solved);
    // The end can no longer slide along the line.
    assert_eq!(dof(&m, f), 5);
}

#[test]
fn trimming_a_middle_piece_leaves_two_collinear_lines() {
    let (mut m, mut ed) = crossed_line(&[-10.0, 10.0]);
    let f = ed.feature;
    let map = view();
    // A line hanging off the far end, and a length on the long line.
    draw(
        &mut m,
        &mut ed,
        &map,
        Tool::Line,
        &[(40.0, 0.0), (40.0, 20.0)],
    );
    let (a, a_start, a_end) = line_near(&m, f, (-40.0, 0.0), (40.0, 0.0));
    let (_, hang, _) = line_near(&m, f, (40.0, 0.0), (40.0, 20.0));
    let (b, _, _) = line_near(&m, f, (-10.0, -30.0), (-10.0, 30.0));
    let (c, _, _) = line_near(&m, f, (10.0, -30.0), (10.0, 30.0));
    constrain(
        &mut m,
        f,
        Constraint::Length {
            line: a,
            value: 0.08,
        },
    );
    let before = dof(&m, f);
    let message = act(&mut m, &mut ed, Tool::Trim, 0.0, 0.0);
    assert_eq!(
        message.as_deref(),
        Some("Trimmed; removed 1 constraint that no longer applied")
    );
    assert_eq!(
        constraints(&m, f, |c| matches!(c, Constraint::Length { .. })),
        0
    );
    let (d, d_start, d_end) = line_near(&m, f, (10.0, 0.0), (40.0, 0.0));
    assert_ne!(d, a);
    assert!((pos(&m, f, a_start) - DVec2::new(-40.0, 0.0) * MM).length() < 1e-6);
    assert!((pos(&m, f, a_end) - DVec2::new(-10.0, 0.0) * MM).length() < 1e-6);
    // The new line lies on the old one, starts on its cutter, and its end
    // took over the old end's coincidence with the hanging line.
    for p in [d_start, d_end] {
        assert!(has(&m, f, &Constraint::PointOnCurve { point: p, curve: a }));
    }
    assert!(has(
        &m,
        f,
        &Constraint::PointOnCurve {
            point: a_end,
            curve: b
        }
    ));
    assert!(has(
        &m,
        f,
        &Constraint::PointOnCurve {
            point: d_start,
            curve: c
        }
    ));
    assert_eq!(
        constraints(&m, f, |x| matches!(x, Constraint::Coincident { a, b }
            if (*a == d_end && *b == hang) || (*a == hang && *b == d_end))),
        1
    );
    assert_eq!(status(&m, f), SolveStatus::Solved);
    // Losing the length frees one; the new line adds 4 and takes 3 back
    // (on the line twice, on its cutter once); the old end loses 1.
    assert_eq!(dof(&m, f), before + 1 + 4 - 3 - 1);
    // One undo step brings the whole line back.
    assert!(m.undo());
    assert_eq!(count(&m, f, |g| matches!(g, Geometry::Line { .. })), 4);
    assert_eq!(
        constraints(&m, f, |c| matches!(c, Constraint::Length { .. })),
        1
    );
}

#[test]
fn trimming_an_uncut_curve_deletes_it_and_its_points() {
    let (mut m, mut ed, map) = setup();
    let f = ed.feature;
    draw(
        &mut m,
        &mut ed,
        &map,
        Tool::Line,
        &[(-40.0, 0.0), (40.0, 0.0)],
    );
    draw(
        &mut m,
        &mut ed,
        &map,
        Tool::Circle,
        &[(0.0, 50.0), (10.0, 50.0)],
    );
    act(&mut m, &mut ed, Tool::Trim, 20.0, 0.0);
    assert_eq!(
        count(&m, f, |_| true),
        2,
        "the circle and its centre remain"
    );
    act(&mut m, &mut ed, Tool::Trim, 10.0, 50.0);
    assert_eq!(count(&m, f, |_| true), 0);
}

#[test]
fn a_circle_cut_twice_becomes_an_arc_with_its_constraints() {
    let (mut m, mut ed, map) = setup();
    let f = ed.feature;
    draw(
        &mut m,
        &mut ed,
        &map,
        Tool::Circle,
        &[(0.0, 0.0), (20.0, 0.0)],
    );
    draw(
        &mut m,
        &mut ed,
        &map,
        Tool::Line,
        &[(-40.0, 0.0), (40.0, 0.0)],
    );
    let circle = m
        .sketch_def(f)
        .unwrap()
        .sketch
        .entities
        .iter()
        .find(|e| matches!(e.geometry, Geometry::Circle { .. }))
        .unwrap()
        .id;
    constrain(
        &mut m,
        f,
        Constraint::Radius {
            curve: circle,
            value: 0.02,
        },
    );
    let before = dof(&m, f);
    // The upper half goes; the lower half stays as an arc.
    let message = act(&mut m, &mut ed, Tool::Trim, 0.0, 20.0);
    assert_eq!(message.as_deref(), Some("Trimmed"));
    assert_eq!(count(&m, f, |g| matches!(g, Geometry::Circle { .. })), 0);
    let arcs = arc_ids(&m, f);
    assert_eq!(arcs.len(), 1);
    let (arc, _, s, e) = arcs[0];
    assert!((pos(&m, f, s) - DVec2::new(-20.0, 0.0) * MM).length() < 1e-6);
    assert!((pos(&m, f, e) - DVec2::new(20.0, 0.0) * MM).length() < 1e-6);
    assert!(has(
        &m,
        f,
        &Constraint::Radius {
            curve: arc,
            value: 0.02
        }
    ));
    let (line, _, _) = line_ids(&m, f)[0];
    assert!(has(
        &m,
        f,
        &Constraint::PointOnCurve {
            point: s,
            curve: line
        }
    ));
    assert!(has(
        &m,
        f,
        &Constraint::PointOnCurve {
            point: e,
            curve: line
        }
    ));
    assert_eq!(status(&m, f), SolveStatus::Solved);
    assert_eq!(dof(&m, f), before);
    // The line and the lower half bound a half disc.
    let areas = region_areas(&m, f);
    assert_eq!(areas.len(), 1);
    assert!(
        (areas[0] - PI * 400.0 * MM * MM / 2.0).abs() < 1e-12,
        "{areas:?}"
    );
}

#[test]
fn a_circle_crossed_once_cannot_be_trimmed() {
    let (mut m, mut ed, map) = setup();
    let f = ed.feature;
    draw(
        &mut m,
        &mut ed,
        &map,
        Tool::Circle,
        &[(0.0, 0.0), (20.0, 0.0)],
    );
    draw(
        &mut m,
        &mut ed,
        &map,
        Tool::Line,
        &[(5.0, 5.0), (40.0, 5.0)],
    );
    let before = m.sketch_def(f).unwrap().clone();
    let message = act(&mut m, &mut ed, Tool::Trim, 0.0, -20.0).unwrap();
    assert!(message.contains("two crossings"), "{message}");
    assert_eq!(*m.sketch_def(f).unwrap(), before);
}

#[test]
fn trimming_the_middle_of_an_arc_leaves_two_arcs_on_one_centre() {
    let (mut m, mut ed, map) = setup();
    let f = ed.feature;
    // A half circle above the x axis, cut by two vertical lines.
    centre_arc(
        &mut m,
        &mut ed,
        (0.0, 0.0),
        (30.0, 0.0),
        (0.0, 30.0),
        (-30.0, 0.0),
    );
    assert_eq!(arc_ids(&m, f).len(), 1);
    draw(
        &mut m,
        &mut ed,
        &map,
        Tool::Line,
        &[(-10.0, -10.0), (-10.0, 40.0)],
    );
    draw(
        &mut m,
        &mut ed,
        &map,
        Tool::Line,
        &[(10.0, -10.0), (10.0, 40.0)],
    );
    let before = dof(&m, f);
    act(&mut m, &mut ed, Tool::Trim, 0.0, 30.0);
    let arcs = arc_ids(&m, f);
    assert_eq!(arcs.len(), 2);
    assert_eq!(arcs[0].1, arcs[1].1, "one centre point");
    assert!(has(
        &m,
        f,
        &Constraint::Equal {
            a: arcs[0].0,
            b: arcs[1].0
        }
    ));
    let y = (900.0f64 - 100.0).sqrt();
    let ends: Vec<DVec2> = arcs
        .iter()
        .flat_map(|a| [pos(&m, f, a.2), pos(&m, f, a.3)])
        .collect();
    for p in [(30.0, 0.0), (10.0, y), (-10.0, y), (-30.0, 0.0)] {
        let p = DVec2::new(p.0, p.1) * MM;
        assert!(
            ends.iter().any(|q| (*q - p).length() < 1e-6),
            "{p} in {ends:?}"
        );
    }
    assert_eq!(status(&m, f), SolveStatus::Solved);
    // The new arc's start and end add 4 and take back the circle (1), the
    // equal radius (1) and its cutter (1); the old end loses 1 to its cutter.
    assert_eq!(dof(&m, f), before + 4 - 3 - 1);
}

#[test]
fn trimming_overhangs_leaves_two_regions_to_extrude() {
    let (mut m, mut ed, map) = setup();
    let f = ed.feature;
    draw(
        &mut m,
        &mut ed,
        &map,
        Tool::Rectangle,
        &[(-20.0, -10.0), (20.0, 10.0)],
    );
    draw(
        &mut m,
        &mut ed,
        &map,
        Tool::Line,
        &[(0.0, -30.0), (0.0, 30.0)],
    );
    act(&mut m, &mut ed, Tool::Trim, 0.0, 25.0);
    act(&mut m, &mut ed, Tool::Trim, 0.0, -25.0);
    let (_, s, e) = line_near(&m, f, (0.0, -10.0), (0.0, 10.0));
    let top = line_near(&m, f, (-20.0, 10.0), (20.0, 10.0)).0;
    let bottom = line_near(&m, f, (-20.0, -10.0), (20.0, -10.0)).0;
    let on = |p, c| Constraint::PointOnCurve { point: p, curve: c };
    assert!(
        (has(&m, f, &on(e, top)) && has(&m, f, &on(s, bottom)))
            || (has(&m, f, &on(s, top)) && has(&m, f, &on(e, bottom)))
    );
    assert_eq!(status(&m, f), SolveStatus::Solved);
    let areas = region_areas(&m, f);
    assert_eq!(areas.len(), 2, "{areas:?}");
    for a in &areas {
        assert!((a - 400.0 * MM * MM).abs() < 1e-12, "{areas:?}");
    }
    let x = extrude(&mut m, &map, &[(-10.0, 0.0), (10.0, 0.0)], "10 mm");
    assert!((volume(&m, x) - 8000.0 * MM * MM * MM).abs() < 1e-15);
}

#[test]
fn a_tangency_at_the_old_end_moves_to_the_new_line() {
    let (mut m, mut ed, map) = setup();
    let f = ed.feature;
    draw(
        &mut m,
        &mut ed,
        &map,
        Tool::Line,
        &[(-40.0, 0.0), (0.0, 0.0)],
    );
    draw(
        &mut m,
        &mut ed,
        &map,
        Tool::TangentArc,
        &[(0.0, 0.0), (20.0, 20.0)],
    );
    draw(
        &mut m,
        &mut ed,
        &map,
        Tool::Line,
        &[(-20.0, -10.0), (-20.0, 10.0)],
    );
    draw(
        &mut m,
        &mut ed,
        &map,
        Tool::Line,
        &[(-10.0, -10.0), (-10.0, 10.0)],
    );
    let (a, _, _) = line_near(&m, f, (-40.0, 0.0), (0.0, 0.0));
    let arc = arc_ids(&m, f)[0].0;
    let tangent_to = |m: &Model, line: SketchId| {
        constraints(m, f, |c| {
            matches!(c, Constraint::Tangent { a, b, kind: TangentKind::AtPoint { .. } }
            if (*a == line && *b == arc) || (*a == arc && *b == line))
        })
    };
    assert_eq!(tangent_to(&m, a), 1);
    let before = dof(&m, f);
    act(&mut m, &mut ed, Tool::Trim, -15.0, 0.0);
    let (d, _, _) = line_near(&m, f, (-10.0, 0.0), (0.0, 0.0));
    assert_eq!(tangent_to(&m, a), 0);
    assert_eq!(tangent_to(&m, d), 1);
    assert_eq!(status(&m, f), SolveStatus::Solved);
    assert_eq!(redundant(&m, f), 0);
    assert_eq!(dof(&m, f), before + 4 - 3 - 1);
    // Splitting the new line hands the tangency on again.
    act(&mut m, &mut ed, Tool::Split, -5.0, 0.0);
    let (e, _, _) = line_near(&m, f, (-5.0, 0.0), (0.0, 0.0));
    assert_eq!(tangent_to(&m, d), 0);
    assert_eq!(tangent_to(&m, e), 1);
    assert_eq!(status(&m, f), SolveStatus::Solved);
}

#[test]
fn extending_a_line_reaches_the_next_curve() {
    let (mut m, mut ed, map) = setup();
    let f = ed.feature;
    draw(
        &mut m,
        &mut ed,
        &map,
        Tool::Line,
        &[(-40.0, 0.0), (-10.0, 0.0)],
    );
    draw(
        &mut m,
        &mut ed,
        &map,
        Tool::Line,
        &[(20.0, -30.0), (20.0, 30.0)],
    );
    draw(
        &mut m,
        &mut ed,
        &map,
        Tool::Line,
        &[(50.0, -30.0), (50.0, 30.0)],
    );
    let (a, a_start, a_end) = line_near(&m, f, (-40.0, 0.0), (-10.0, 0.0));
    let (b, _, _) = line_near(&m, f, (20.0, -30.0), (20.0, 30.0));
    constrain(
        &mut m,
        f,
        Constraint::Length {
            line: a,
            value: 0.03,
        },
    );
    // The preview is the added piece.
    ed.set_tool(Tool::Extend);
    hover(&m, &mut ed, &map, -12.0, 0.0);
    match ed.curve_edit.clone().unwrap().plan.unwrap().preview {
        Preview::Piece(ResolvedCurve::Line { start, end }) => {
            assert!((start - DVec2::new(-10.0, 0.0) * MM).length() < 1e-6);
            assert!((end - DVec2::new(20.0, 0.0) * MM).length() < 1e-6);
        }
        p => panic!("{p:?}"),
    }
    let message = act(&mut m, &mut ed, Tool::Extend, -12.0, 0.0);
    assert_eq!(
        message.as_deref(),
        Some("Extended; removed 1 constraint that no longer applied")
    );
    assert!((pos(&m, f, a_end) - DVec2::new(20.0, 0.0) * MM).length() < 1e-6);
    assert!((pos(&m, f, a_start) - DVec2::new(-40.0, 0.0) * MM).length() < 1e-6);
    assert!(has(
        &m,
        f,
        &Constraint::PointOnCurve {
            point: a_end,
            curve: b
        }
    ));
    assert_eq!(status(&m, f), SolveStatus::Solved);
    // Again: on to the next line.
    act(&mut m, &mut ed, Tool::Extend, 15.0, 0.0);
    assert!((pos(&m, f, a_end) - DVec2::new(50.0, 0.0) * MM).length() < 1e-6);
    // And then nothing is left beyond it.
    let message = act(&mut m, &mut ed, Tool::Extend, 45.0, 0.0).unwrap();
    assert!(message.contains("Nothing to extend to"), "{message}");
    // The start end extends backwards: nothing there either.
    let message = act(&mut m, &mut ed, Tool::Extend, -38.0, 0.0).unwrap();
    assert!(message.contains("Nothing to extend to"), "{message}");
}

#[test]
fn extending_an_arc_follows_its_circle() {
    let (mut m, mut ed, map) = setup();
    let f = ed.feature;
    // A quarter arc from +x to +y, and a vertical line at x = −10.
    centre_arc(
        &mut m,
        &mut ed,
        (0.0, 0.0),
        (20.0, 0.0),
        (14.0, 14.0),
        (0.0, 20.0),
    );
    draw(
        &mut m,
        &mut ed,
        &map,
        Tool::Line,
        &[(-10.0, -40.0), (-10.0, 40.0)],
    );
    let (_, _, s, e) = arc_ids(&m, f)[0];
    act(&mut m, &mut ed, Tool::Extend, 2.0, 19.9);
    let y = (400.0f64 - 100.0).sqrt();
    assert!((pos(&m, f, e) - DVec2::new(-10.0, y) * MM).length() < 1e-6);
    assert!((pos(&m, f, s) - DVec2::new(20.0, 0.0) * MM).length() < 1e-6);
    // The start runs clockwise round to the line's other crossing.
    act(&mut m, &mut ed, Tool::Extend, 19.97, 1.0);
    assert!((pos(&m, f, s) - DVec2::new(-10.0, -y) * MM).length() < 1e-6);
    assert_eq!(status(&m, f), SolveStatus::Solved);
}

#[test]
fn extending_keeps_curves_that_met_at_the_end_joined() {
    let (mut m, mut ed, map) = setup();
    let f = ed.feature;
    draw(
        &mut m,
        &mut ed,
        &map,
        Tool::Line,
        &[(-40.0, 0.0), (0.0, 0.0)],
    );
    draw(
        &mut m,
        &mut ed,
        &map,
        Tool::Line,
        &[(0.0, 0.0), (0.0, 30.0)],
    );
    draw(
        &mut m,
        &mut ed,
        &map,
        Tool::Line,
        &[(0.0, 0.0), (30.0, 30.0)],
    );
    draw(
        &mut m,
        &mut ed,
        &map,
        Tool::Line,
        &[(20.0, -30.0), (20.0, 30.0)],
    );
    let (_, _, a_end) = line_near(&m, f, (-40.0, 0.0), (0.0, 0.0));
    let (_, b, _) = line_near(&m, f, (0.0, 0.0), (0.0, 30.0));
    let (_, c, _) = line_near(&m, f, (0.0, 0.0), (30.0, 30.0));
    act(&mut m, &mut ed, Tool::Extend, -3.0, 0.0);
    assert!((pos(&m, f, a_end) - DVec2::new(20.0, 0.0) * MM).length() < 1e-6);
    assert!((pos(&m, f, b) - pos(&m, f, c)).length() < 1e-9);
    assert_eq!(
        constraints(&m, f, |x| matches!(x, Constraint::Coincident { .. })),
        1
    );
    assert_eq!(status(&m, f), SolveStatus::Solved);
    // The two stay joined when one moves.
    constrain(
        &mut m,
        f,
        Constraint::Fix {
            point: b,
            at: DVec2::new(0.005, 0.0),
        },
    );
    assert!((pos(&m, f, c) - DVec2::new(0.005, 0.0)).length() < 1e-9);
}

#[test]
fn splitting_a_line_frees_the_split_point_along_it() {
    let (mut m, mut ed, map) = setup();
    let f = ed.feature;
    draw(
        &mut m,
        &mut ed,
        &map,
        Tool::Line,
        &[(-40.0, 0.0), (40.0, 0.0)],
    );
    let (a, _, a_end) = line_near(&m, f, (-40.0, 0.0), (40.0, 0.0));
    ed.set_tool(Tool::Split);
    hover(&m, &mut ed, &map, 10.0, 0.5);
    match ed.curve_edit.clone().unwrap().plan.unwrap().preview {
        Preview::Point(q) => assert!((q - DVec2::new(10.0, 0.0) * MM).length() < 1e-6),
        p => panic!("{p:?}"),
    }
    let before = dof(&m, f);
    act(&mut m, &mut ed, Tool::Split, 10.0, 0.5);
    let (d, d_start, d_end) = line_near(&m, f, (10.0, 0.0), (40.0, 0.0));
    assert!((pos(&m, f, a_end) - DVec2::new(10.0, 0.0) * MM).length() < 1e-6);
    assert!(has(
        &m,
        f,
        &Constraint::Coincident {
            a: a_end,
            b: d_start
        }
    ));
    assert!(has(
        &m,
        f,
        &Constraint::PointOnCurve {
            point: d_end,
            curve: a
        }
    ));
    assert_ne!(d, a);
    assert_eq!(status(&m, f), SolveStatus::Solved);
    assert_eq!(dof(&m, f), before + 1);
    // Too close to an end, and circles, are refused.
    let message = act(&mut m, &mut ed, Tool::Split, 10.0, 0.0).unwrap();
    assert!(message.contains("Too close"), "{message}");
    draw(
        &mut m,
        &mut ed,
        &map,
        Tool::Circle,
        &[(0.0, 50.0), (10.0, 50.0)],
    );
    let message = act(&mut m, &mut ed, Tool::Split, 10.0, 50.0).unwrap();
    assert!(message.contains("circle"), "{message}");
}

#[test]
fn splitting_an_arc_shares_its_centre() {
    let (mut m, mut ed, _) = setup();
    let f = ed.feature;
    centre_arc(
        &mut m,
        &mut ed,
        (0.0, 0.0),
        (20.0, 0.0),
        (0.0, 20.0),
        (-20.0, 0.0),
    );
    let before = dof(&m, f);
    act(&mut m, &mut ed, Tool::Split, 0.0, 20.0);
    let arcs = arc_ids(&m, f);
    assert_eq!(arcs.len(), 2);
    assert_eq!(arcs[0].1, arcs[1].1);
    assert_eq!(
        constraints(&m, f, |c| matches!(c, Constraint::Equal { .. })),
        0
    );
    assert_eq!(status(&m, f), SolveStatus::Solved);
    assert_eq!(dof(&m, f), before + 1);
}

#[test]
fn checked_edits_refuse_to_break_a_sketch() {
    let (mut m, mut ed, map) = setup();
    let f = ed.feature;
    draw(
        &mut m,
        &mut ed,
        &map,
        Tool::Line,
        &[(-40.0, 0.0), (40.0, 0.0)],
    );
    let (a, s, _) = line_near(&m, f, (-40.0, 0.0), (40.0, 0.0));
    let ops = m.session().log().len();
    let add = |c| {
        vec![SketchEdit::AddConstraint(ConstraintEntry {
            id: SketchId::new_in_batch(0),
            constraint: c,
        })]
    };
    // Already horizontal: a second one is redundant.
    let e1 = m
        .edit_sketch_checked(f, add(Constraint::Horizontal { line: a }))
        .unwrap_err();
    assert!(e1.contains("redundant"), "{e1}");
    // A vertical horizontal line cannot be.
    let e2 = m
        .edit_sketch_checked(f, add(Constraint::Vertical { line: a }))
        .unwrap_err();
    assert!(e2.contains("no longer solve"), "{e2}");
    assert_eq!(m.session().log().len(), ops);
    m.edit_sketch_checked(
        f,
        add(Constraint::Fix {
            point: s,
            at: DVec2::ZERO,
        }),
    )
    .unwrap();
    assert_eq!(m.session().log().len(), ops + 1);
}

/// A random curve (millimetres).
#[derive(Clone, Debug)]
enum Shape {
    Line(DVec2, DVec2),
    Circle(DVec2, f64),
    /// Centre, radius, start angle, sweep.
    Arc(DVec2, f64, f64, f64),
}

fn shape() -> impl Strategy<Value = Shape> {
    let p = || (-50.0..50.0f64, -50.0..50.0f64).prop_map(|(x, y)| DVec2::new(x, y));
    prop_oneof![
        3 => (p(), p())
            .prop_filter("a visible line", |(a, b)| (*a - *b).length() > 2.0)
            .prop_map(|(a, b)| Shape::Line(a, b)),
        1 => (p(), 3.0..30.0f64).prop_map(|(c, r)| Shape::Circle(c, r)),
        1 => (p(), 3.0..30.0f64, 0.0..TAU, 0.3..6.0f64)
            .prop_map(|(c, r, a, s)| Shape::Arc(c, r, a, s)),
    ]
}

/// Adds the shapes as one batch.
fn add_shapes(m: &mut Model, f: FeatureId, shapes: &[Shape]) {
    let mut b = Batch::default();
    for s in shapes {
        match s {
            Shape::Line(a, e) => {
                b.line(*a * MM, *e * MM, false);
            }
            Shape::Circle(c, r) => {
                let center = b.point(*c * MM);
                b.entity(
                    Geometry::Circle {
                        center,
                        radius: r * MM,
                    },
                    false,
                );
            }
            Shape::Arc(c, r, a, sweep) => {
                let at = |t: f64| (*c + DVec2::from_angle(t) * *r) * MM;
                b.arc(*c * MM, at(*a), at(a + sweep), false);
            }
        }
    }
    m.edit_sketch(f, b.edits).unwrap();
}

/// The point a fraction `along` of the way along a curve (in its
/// parameter, for ellipses its eccentric anomaly, for conics t).
fn along(c: &ResolvedCurve, along: f64) -> DVec2 {
    match c {
        ResolvedCurve::Line { start, end } => start.lerp(*end, along),
        ResolvedCurve::Arc {
            center,
            radius,
            start_angle,
            sweep_angle,
        } => *center + DVec2::from_angle(start_angle + along * sweep_angle) * *radius,
        ResolvedCurve::Circle { center, radius } => {
            *center + DVec2::from_angle(along * TAU) * *radius
        }
        ResolvedCurve::Spline { points, .. } => points[0],
        ResolvedCurve::Ellipse {
            center,
            major,
            minor,
        } => {
            let (s, c) = (along * TAU).sin_cos();
            *center + *major * c + *minor * s
        }
        ResolvedCurve::EllipticalArc {
            center,
            major,
            minor,
            start_param,
            sweep_param,
        } => {
            let (s, c) = (start_param + along * sweep_param).sin_cos();
            *center + *major * c + *minor * s
        }
        ResolvedCurve::Conic {
            start,
            control,
            end,
            weight,
        } => openshape_sketch_solver::conic::point([*start, *control, *end], *weight, along),
    }
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(256))]

    /// Random trims, extends and splits of random lines, arcs and circles
    /// never break the sketch: every plan commits through the checked path
    /// (no failure, no new redundancy), the sketch stays solved with its
    /// regions, points the plan did not move stay put, and moved points
    /// land where planned.
    #[test]
    fn random_edits_keep_the_sketch_solved(
        shapes in proptest::collection::vec(shape(), 3..8),
        steps in proptest::collection::vec((0usize..3, 0usize..64, 0.02..0.98f64), 1..7),
    ) {
        let mut m = Model::new();
        let f = m.add_sketch(PlaneRef::Top);
        add_shapes(&mut m, f, &shapes);
        prop_assume!(status(&m, f) == SolveStatus::Solved);
        edit_randomly(&mut m, f, &steps)?;
    }
}

/// A random shape built by the sketch tools' own builders (millimetres).
#[derive(Clone, Debug)]
enum Built {
    Rectangle(DVec2, DVec2),
    Slot(DVec2, DVec2, DVec2),
    Polygon(DVec2, DVec2, u32),
    /// Lines joined end to end by coincidences.
    Chain(Vec<DVec2>),
    Circle(DVec2, f64),
}

fn built() -> impl Strategy<Value = Built> {
    let p = || (-50.0..50.0f64, -50.0..50.0f64).prop_map(|(x, y)| DVec2::new(x, y));
    prop_oneof![
        (p(), p()).prop_map(|(c, k)| Built::Rectangle(c, k)),
        (p(), p(), p()).prop_map(|(a, b, w)| Built::Slot(a, b, w)),
        (p(), p(), 3u32..9).prop_map(|(c, v, n)| Built::Polygon(c, v, n)),
        proptest::collection::vec(p(), 2..5).prop_map(Built::Chain),
        (p(), 3.0..30.0f64).prop_map(|(c, r)| Built::Circle(c, r)),
    ]
}

/// Adds a shape as one batch; degenerate shapes are skipped.
fn add_built(m: &mut Model, f: FeatureId, s: &Built) {
    let free = |p: DVec2| Snap {
        at: p * MM,
        target: SnapTarget::Free,
    };
    let mut b = Batch::default();
    let made = match s {
        Built::Rectangle(c, k) => shapes::centre_rectangle(&mut b, &free(*c), &free(*k), false),
        Built::Slot(a, e, w) => shapes::slot(&mut b, &free(*a), &free(*e), *w * MM, false),
        Built::Polygon(c, v, n) => shapes::polygon(&mut b, &free(*c), &free(*v), *n, false),
        Built::Chain(points) => {
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
        Built::Circle(c, r) => {
            let center = b.point(*c * MM);
            b.entity(
                Geometry::Circle {
                    center,
                    radius: r * MM,
                },
                false,
            );
            Ok(())
        }
    };
    if made.is_ok() {
        m.edit_sketch(f, b.edits).unwrap();
    }
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(256))]

    /// As `random_edits_keep_the_sketch_solved`, on shapes that carry
    /// constraints (rectangles, slots with tangent joints, polygons on a
    /// circle, joined chains), so the constraint hand-over is exercised.
    #[test]
    fn random_edits_of_constrained_shapes_keep_the_sketch_solved(
        shapes in proptest::collection::vec(built(), 2..5),
        steps in proptest::collection::vec((0usize..3, 0usize..64, 0.02..0.98f64), 1..7),
    ) {
        let mut m = Model::new();
        let f = m.add_sketch(PlaneRef::Top);
        for s in &shapes {
            add_built(&mut m, f, s);
        }
        prop_assume!(status(&m, f) == SolveStatus::Solved);
        prop_assume!(redundant(&m, f) == 0);
        edit_randomly(&mut m, f, &steps)?;
    }
}

/// Applies random trims, extends and splits (`(op, curve, fraction along
/// it)`) and checks each one (see `random_edits_keep_the_sketch_solved`).
fn edit_randomly(
    m: &mut Model,
    f: FeatureId,
    steps: &[(usize, usize, f64)],
) -> Result<(), TestCaseError> {
    for &(op, which, t) in steps {
        let r = m.sketch_result(f).unwrap();
        let curves: Vec<(SketchId, ResolvedCurve)> =
            r.curves.iter().map(|(id, _, c)| (*id, c.clone())).collect();
        if curves.is_empty() {
            break;
        }
        let (target, curve) = curves[which % curves.len()].clone();
        let near = along(&curve, t);
        let sketch = r.solution.sketch.clone();
        let planned = match op {
            0 => {
                trim::cuts(m.kernel(), &curves, target).and_then(|k| trim::trim(&sketch, &k, near))
            }
            1 => trim::extend(m.kernel(), &sketch, &curves, target, near),
            _ => trim::split(&sketch, &curves, target, near),
        };
        let plan = match planned {
            Ok(p) => p,
            // Refusals are fine (nothing to extend to, a circle cut
            // once, too close to an end); they change nothing.
            Err(_) => continue,
        };
        let moved: Vec<(SketchId, DVec2)> = plan
            .edits
            .iter()
            .flat_map(|e| match e {
                SketchEdit::MovePoints { points } => points.clone(),
                _ => Vec::new(),
            })
            .collect();
        if let Err(e) = m.edit_sketch_checked(f, plan.edits.clone()) {
            return Err(TestCaseError::fail(format!("op {op} on {curve:?}: {e}")));
        }
        prop_assert_eq!(status(m, f), SolveStatus::Solved);
        prop_assert_eq!(redundant(m, f), 0);
        let after = m.sketch_result(f).unwrap();
        // Regions are found, unless curves touch so closely that the
        // kernel refuses them by design (random geometry can make a sliver
        // thinner than the tolerance, design note 0006).
        if let Err(e) = &after.regions {
            prop_assert!(e.contains("within the model tolerance"), "{}", e);
        }
        for (_, _, c) in &after.curves {
            let size = match c {
                ResolvedCurve::Line { start, end } => (*end - *start).length(),
                ResolvedCurve::Arc {
                    radius,
                    sweep_angle,
                    ..
                } => radius * sweep_angle,
                ResolvedCurve::Circle { radius, .. } => *radius,
                ResolvedCurve::Spline { .. } => 1.0,
                ResolvedCurve::Ellipse { major, minor, .. } => major.length().min(minor.length()),
                ResolvedCurve::EllipticalArc {
                    major,
                    minor,
                    sweep_param,
                    ..
                } => major.length().min(minor.length()) * sweep_param.abs(),
                ResolvedCurve::Conic { start, end, .. } => (*end - *start).length(),
            };
            prop_assert!(size > trim::CUT_GAP, "{:?}", c);
        }
        for e in &sketch.entities {
            let Geometry::Point { at } = e.geometry else {
                continue;
            };
            let Some(now) = after.solution.sketch.point(e.id) else {
                continue;
            };
            let to = moved
                .iter()
                .find(|(id, _)| *id == e.id)
                .map_or(at, |(_, p)| *p);
            prop_assert!(
                (now - to).length() < 1e-9,
                "{:?} moved from {} to {}",
                e.id,
                to,
                now
            );
        }
    }
    Ok(())
}
