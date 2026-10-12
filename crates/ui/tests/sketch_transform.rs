//! Move, rotate and scale (design note 0007 §13, ADR 0014), driven through
//! the editor: exact typed and clicked transforms, what happens to
//! constraints among the selection, to the frame and to the rest of the
//! sketch, copies, refusals, and a fuzz of random transforms.

mod common;

use std::collections::BTreeSet;
use std::f64::consts::PI;

use common::*;
use opendatum_app::sketch_mode::{ScreenMap, SketchEditor, Snap, SnapTarget, Tool};
use opendatum_app::sketch_ops::transform::{moved_set, Similarity};
use opendatum_app::sketch_ops::trim::Preview;
use opendatum_app::sketch_ops::{mirror, Batch};
use opendatum_app::Model;
use opendatum_core_model::sketch::{Constraint, DistanceKind, Geometry};
use opendatum_core_model::{DVec2, FeatureId, SketchData, SketchEdit, SketchId};
use opendatum_sketch_solver::SolveStatus;
use proptest::prelude::*;

fn solved(m: &Model, f: FeatureId) -> SketchData {
    m.sketch_result(f).unwrap().solution.sketch.clone()
}

fn redundant(m: &Model, f: FeatureId) -> usize {
    m.sketch_result(f).unwrap().solution.redundant.len()
}

/// The ids of the sketch's entities matching `pred`.
fn ids<F: Fn(&Geometry<SketchId>) -> bool>(m: &Model, f: FeatureId, pred: F) -> Vec<SketchId> {
    m.sketch_def(f)
        .unwrap()
        .sketch
        .entities
        .iter()
        .filter(|e| pred(&e.geometry))
        .map(|e| e.id)
        .collect()
}

fn lines(m: &Model, f: FeatureId) -> Vec<SketchId> {
    ids(m, f, |g| matches!(g, Geometry::Line { .. }))
}

/// A `w` × `h` mm rectangle with its lower left corner at `(x, y)`, drawn
/// with the rectangle tool (a corner at the origin is fixed there); its
/// four lines.
fn rectangle(
    m: &mut Model,
    ed: &mut SketchEditor,
    map: &ScreenMap,
    (x, y): (f64, f64),
    (w, h): (f64, f64),
) -> Vec<SketchId> {
    let before: BTreeSet<SketchId> = lines(m, ed.feature).into_iter().collect();
    draw(m, ed, map, Tool::Rectangle, &[(x, y), (x + w, y + h)]);
    lines(m, ed.feature)
        .into_iter()
        .filter(|l| !before.contains(l))
        .collect()
}

/// The line among `lines` whose midpoint is nearest `(x, y)` mm.
fn line_near(m: &Model, f: FeatureId, lines: &[SketchId], x: f64, y: f64) -> SketchId {
    let s = solved(m, f);
    let mid = |l: SketchId| {
        let Geometry::Line { start, end } = s.entity(l).unwrap().geometry else {
            unreachable!()
        };
        (s.point(start).unwrap() + s.point(end).unwrap()) * 0.5
    };
    let q = DVec2::new(x, y) * MM;
    *lines
        .iter()
        .min_by(|a, b| (mid(**a) - q).length().total_cmp(&(mid(**b) - q).length()))
        .unwrap()
}

/// Adds constraints through the model (one undo step).
fn add(m: &mut Model, f: FeatureId, cs: Vec<Constraint<SketchId>>) {
    let mut b = Batch::default();
    for c in cs {
        b.constraint(c);
    }
    m.edit_sketch(f, b.edits).unwrap();
    assert!(
        !matches!(status(m, f), SolveStatus::Failed(_)),
        "{:?}",
        status(m, f)
    );
}

/// Every point of `before` moved by `t` in `after` (within `tol` metres),
/// for the points in `moved`; every other point within `still` metres of
/// where it was, if given.
fn moved_by(
    before: &SketchData,
    after: &SketchData,
    moved: &BTreeSet<SketchId>,
    t: &Similarity,
    tol: f64,
    still: Option<f64>,
) {
    for e in &before.entities {
        let Geometry::Point { at } = e.geometry else {
            continue;
        };
        let now = after.point(e.id).unwrap();
        if moved.contains(&e.id) {
            let want = t.point(at);
            assert!(
                (now - want).length() <= tol,
                "{:?} at {now:?}, not {want:?}",
                e.id
            );
        } else if let Some(still) = still {
            assert!(
                (now - at).length() <= still,
                "{:?} moved by {:e}",
                e.id,
                (now - at).length()
            );
        }
    }
}

fn count_of<F: Fn(&Constraint<SketchId>) -> bool>(m: &Model, f: FeatureId, pred: F) -> usize {
    constraints(m, f, pred)
}

fn horizontal(c: &Constraint<SketchId>) -> bool {
    matches!(c, Constraint::Horizontal { .. })
}

fn vertical(c: &Constraint<SketchId>) -> bool {
    matches!(c, Constraint::Vertical { .. })
}

/// The dimension values of the sketch, in constraint order.
fn values(m: &Model, f: FeatureId) -> Vec<f64> {
    m.sketch_def(f)
        .unwrap()
        .sketch
        .constraints
        .iter()
        .filter_map(|c| c.constraint.value())
        .collect()
}

/// Chooses `tool` with `selection` selected.
fn choose(ed: &mut SketchEditor, tool: Tool, selection: &[SketchId]) {
    ed.set_tool(Tool::Select);
    ed.selection = selection.to_vec();
    ed.set_tool(tool);
}

/// A click the transform tools take as a centre or base point.
fn at_point(x: f64, y: f64) -> Snap {
    Snap {
        at: DVec2::new(x, y) * MM,
        target: SnapTarget::Free,
    }
}

#[test]
fn a_typed_move_is_exact_and_keeps_every_constraint() {
    let (mut m, mut ed, map) = setup();
    let f = ed.feature;
    let square = rectangle(&mut m, &mut ed, &map, (0.0, 0.0), (10.0, 10.0));
    assert_eq!(dof(&m, f), 2);
    let n = m.sketch_def(f).unwrap().sketch.constraints.len();
    let before = solved(&m, f);
    choose(&mut ed, Tool::Move, &square);
    assert_eq!(
        ed.hint(),
        "Click a point to move the 4 selected items from, or type the offsets and Apply"
    );
    ed.transform_offset = DVec2::new(25.0, 10.0) * MM;
    ed.apply_typed_transform(&mut m);
    assert_eq!(ed.message.as_deref(), Some("Moved 4 items"));
    assert_eq!(ed.tool, Tool::Select);
    assert_eq!(ed.selection, square);
    let after = solved(&m, f);
    let t = Similarity::translation(DVec2::new(25.0, 10.0) * MM);
    let moved = moved_set(&before, &square);
    moved_by(&before, &after, &moved, &t, 0.0, Some(0.0));
    // The same constraints; the corner fixed at the origin is fixed where
    // it went.
    assert_eq!(m.sketch_def(f).unwrap().sketch.constraints.len(), n);
    let fixes: Vec<DVec2> = m
        .sketch_def(f)
        .unwrap()
        .sketch
        .constraints
        .iter()
        .filter_map(|c| match c.constraint {
            Constraint::Fix { at, .. } => Some(at),
            _ => None,
        })
        .collect();
    assert_eq!(fixes, vec![DVec2::new(25.0, 10.0) * MM]);
    assert_eq!(dof(&m, f), 2);
    assert_eq!(status(&m, f), SolveStatus::Solved);
    // One undo step puts it back.
    assert!(m.undo());
    moved_by(
        &before,
        &solved(&m, f),
        &BTreeSet::new(),
        &t,
        0.0,
        Some(0.0),
    );
}

#[test]
fn clicked_moves_snap_to_points_and_preview_the_result() {
    let (mut m, mut ed, map) = setup();
    let f = ed.feature;
    let square = rectangle(&mut m, &mut ed, &map, (10.0, 5.0), (10.0, 10.0));
    draw(&mut m, &mut ed, &map, Tool::Point, &[(50.0, 20.0)]);
    let before = solved(&m, f);
    choose(&mut ed, Tool::Move, &square);
    // From the upper right corner to the point.
    click(&mut m, &mut ed, &map, 20.0, 15.0);
    assert_eq!(ed.clicks.len(), 1);
    assert_eq!(ed.hint(), "Click where that point goes; Esc starts again");
    hover(&m, &mut ed, &map, 50.0, 20.0);
    let plan = ed.curve_edit.clone().expect("a plan").plan.expect("ok");
    let Preview::Shapes { curves: shown, .. } = plan.preview else {
        panic!("{:?}", plan.preview)
    };
    assert_eq!(shown.len(), 4);
    assert_eq!(ed.rubber_band().len(), 1);
    // Esc starts again; the next two clicks move it.
    ed.escape();
    assert!(ed.clicks.is_empty() && ed.tool == Tool::Move);
    click(&mut m, &mut ed, &map, 20.0, 15.0);
    click(&mut m, &mut ed, &map, 50.0, 20.0);
    assert_eq!(ed.message.as_deref(), Some("Moved 4 items"));
    let corner = |s: &SketchData, x: f64, y: f64| {
        s.entities
            .iter()
            .filter_map(|e| match e.geometry {
                Geometry::Point { at } => Some(at),
                _ => None,
            })
            .min_by(|a, b| {
                let q = DVec2::new(x, y) * MM;
                (*a - q).length().total_cmp(&(*b - q).length())
            })
            .unwrap()
    };
    let target = corner(&before, 50.0, 20.0);
    let v = target - corner(&before, 20.0, 15.0);
    let after = solved(&m, f);
    moved_by(
        &before,
        &after,
        &moved_set(&before, &square),
        &Similarity::translation(v),
        0.0,
        Some(0.0),
    );
    assert!((corner(&after, 50.0, 20.0) - target).length() < 1e-15);
    assert_eq!(dof(&m, f), 4 + 2);
}

/// A circle at (60, 30) mm of radius 5 mm, dimensioned, and a fixed point
/// at (20, 30) mm it is dimensioned from horizontally; returns the circle,
/// its centre and the two dimensions' ids.
fn dimensioned_circle(
    m: &mut Model,
    ed: &mut SketchEditor,
    map: &ScreenMap,
) -> (SketchId, SketchId, SketchId, SketchId) {
    let f = ed.feature;
    draw(m, ed, map, Tool::Point, &[(20.0, 30.0)]);
    draw(m, ed, map, Tool::Circle, &[(60.0, 30.0), (65.0, 30.0)]);
    let s = solved(m, f);
    let circle = ids(m, f, |g| matches!(g, Geometry::Circle { .. }))[0];
    let Geometry::Circle { center, .. } = s.entity(circle).unwrap().geometry else {
        unreachable!()
    };
    let p = ids(m, f, |g| matches!(g, Geometry::Point { .. }))
        .into_iter()
        .find(|p| *p != center)
        .unwrap();
    let at = s.point(p).unwrap();
    let r = match s.entity(circle).unwrap().geometry {
        Geometry::Circle { radius, .. } => radius,
        _ => unreachable!(),
    };
    let dx = s.point(center).unwrap().x - at.x;
    add(
        m,
        f,
        vec![
            Constraint::Fix { point: p, at },
            Constraint::Distance {
                a: p,
                b: center,
                kind: DistanceKind::Horizontal,
                value: dx,
            },
            Constraint::Radius {
                curve: circle,
                value: r,
            },
        ],
    );
    let dims: Vec<SketchId> = m
        .sketch_def(f)
        .unwrap()
        .sketch
        .constraints
        .iter()
        .filter(|c| c.constraint.value().is_some())
        .map(|c| c.id)
        .collect();
    (circle, center, dims[0], dims[1])
}

fn value_of(m: &Model, f: FeatureId, id: SketchId) -> f64 {
    m.sketch_def(f)
        .unwrap()
        .sketch
        .constraint(id)
        .unwrap()
        .constraint
        .value()
        .unwrap()
}

#[test]
fn dimensions_to_the_rest_are_measured_anew_and_scaled_ones_exactly() {
    let (mut m, mut ed, map) = setup();
    let f = ed.feature;
    let (circle, center, distance, radius) = dimensioned_circle(&mut m, &mut ed, &map);
    let d0 = value_of(&m, f, distance);
    let r0 = value_of(&m, f, radius);
    assert_eq!(dof(&m, f), 1);
    choose(&mut ed, Tool::Move, &[circle]);
    ed.transform_offset = DVec2::new(10.0, 5.0) * MM;
    ed.apply_typed_transform(&mut m);
    assert_eq!(ed.message.as_deref(), Some("Moved 1 item"));
    assert_eq!(value_of(&m, f, distance), d0 + 0.01);
    assert_eq!(value_of(&m, f, radius), r0);
    assert_eq!(dof(&m, f), 1);
    // Scaling about the centre (clicked: it snaps to the point) doubles the
    // radius and its dimension exactly; the distance stays.
    let c = solved(&m, f).point(center).unwrap();
    choose(&mut ed, Tool::Scale, &[circle]);
    assert_eq!(ed.transform_factor, 2.0);
    click(&mut m, &mut ed, &map, c.x / MM, c.y / MM);
    assert_eq!(ed.clicks[0].target, SnapTarget::Point(center));
    ed.apply_typed_transform(&mut m);
    assert_eq!(ed.message.as_deref(), Some("Scaled 1 item"));
    assert_eq!(value_of(&m, f, radius), 2.0 * r0);
    assert_eq!(value_of(&m, f, distance), d0 + 0.01);
    let Geometry::Circle { radius: r, .. } = solved(&m, f).entity(circle).unwrap().geometry else {
        unreachable!()
    };
    assert_eq!(r, 2.0 * r0);
    assert_eq!(solved(&m, f).point(center).unwrap(), c);
}

#[test]
fn kept_constraints_stretch_the_rest_and_detached_ones_free_it() {
    let (mut m, mut ed, map) = setup();
    let f = ed.feature;
    let rect = rectangle(&mut m, &mut ed, &map, (10.0, 5.0), (30.0, 20.0));
    let top = line_near(&m, f, &rect, 25.0, 25.0);
    let coincident = |m: &Model| count_of(m, f, |c| matches!(c, Constraint::Coincident { .. }));
    assert_eq!(coincident(&m), 4);
    assert!(close(region_areas(&m, f)[0], 600.0 * MM * MM));
    choose(&mut ed, Tool::Move, &[top]);
    assert!(ed.transform_keep);
    ed.transform_offset = DVec2::new(0.0, 10.0) * MM;
    ed.apply_typed_transform(&mut m);
    assert_eq!(ed.message.as_deref(), Some("Moved 1 item"));
    // The sides stretched: a 30 × 30 rectangle.
    assert_eq!(coincident(&m), 4);
    assert!(close(region_areas(&m, f)[0], 900.0 * MM * MM));
    // Detached instead, the top line leaves the sides where they were.
    assert!(m.undo());
    let before = solved(&m, f);
    choose(&mut ed, Tool::Move, &[top]);
    ed.transform_keep = false;
    ed.apply_typed_transform(&mut m);
    assert_eq!(
        ed.message.as_deref(),
        Some("Moved 1 item; removed 2 constraints to the rest of the sketch")
    );
    assert_eq!(coincident(&m), 2);
    moved_by(
        &before,
        &solved(&m, f),
        &moved_set(&before, &[top]),
        &Similarity::translation(DVec2::new(0.0, 10.0) * MM),
        0.0,
        Some(0.0),
    );
    // The top line's position and length, and the sides' lengths.
    assert_eq!(dof(&m, f), 4 + 4);
}

#[test]
fn a_transform_the_rest_cannot_follow_is_refused_until_detached() {
    let (mut m, mut ed, map) = setup();
    let f = ed.feature;
    let rect = rectangle(&mut m, &mut ed, &map, (0.0, 0.0), (30.0, 20.0));
    let bottom = line_near(&m, f, &rect, 15.0, 0.0);
    let left = line_near(&m, f, &rect, 0.0, 10.0);
    let top = line_near(&m, f, &rect, 15.0, 20.0);
    add(
        &mut m,
        f,
        vec![
            Constraint::Length {
                line: bottom,
                value: 30.0 * MM,
            },
            Constraint::Length {
                line: left,
                value: 20.0 * MM,
            },
        ],
    );
    assert_eq!(dof(&m, f), 0);
    let before = m.sketch_def(f).unwrap().clone();
    choose(&mut ed, Tool::Move, &[top]);
    click(&mut m, &mut ed, &map, 30.0, 20.0);
    hover(&m, &mut ed, &map, 30.0, 30.0);
    let why = ed.curve_edit.clone().expect("a plan").plan.unwrap_err();
    assert_eq!(
        why,
        "Constraints to the rest of the sketch hold the selection: turn off Keep constraints \
         to detach it"
    );
    click(&mut m, &mut ed, &map, 30.0, 30.0);
    assert_eq!(
        ed.message.as_deref(),
        Some(format!("Not done: {why}").as_str())
    );
    assert_eq!(m.sketch_def(f).unwrap(), &before);
    // Detached, the top line moves and is free to slide.
    ed.transform_keep = false;
    ed.transform_offset = DVec2::new(0.0, 10.0) * MM;
    ed.apply_typed_transform(&mut m);
    assert_eq!(
        ed.message.as_deref(),
        Some("Moved 1 item; removed 2 constraints to the rest of the sketch")
    );
    // The top line's position and length, and the right side's length.
    assert_eq!(dof(&m, f), 4);
}

#[test]
fn quarter_turns_swap_horizontal_and_vertical_exactly() {
    let (mut m, mut ed, map) = setup();
    let f = ed.feature;
    let square = rectangle(&mut m, &mut ed, &map, (0.0, 0.0), (10.0, 10.0));
    let (h, v) = (count_of(&m, f, horizontal), count_of(&m, f, vertical));
    assert_eq!((h, v), (2, 2));
    let before = solved(&m, f);
    let n = before.constraints.len();
    choose(&mut ed, Tool::Rotate, &square);
    assert_eq!(
        ed.hint(),
        "Click the centre to turn the 4 selected items about"
    );
    // Typing needs the centre first.
    ed.apply_typed_transform(&mut m);
    assert_eq!(
        ed.message.as_deref(),
        Some("Click the centre to turn about first")
    );
    ed.clicks.push(at_point(0.0, 0.0));
    assert_eq!(ed.transform_angle, PI / 2.0);
    ed.apply_typed_transform(&mut m);
    assert_eq!(ed.message.as_deref(), Some("Rotated 4 items"));
    let after = solved(&m, f);
    let t = Similarity::rotation(DVec2::ZERO, PI / 2.0);
    moved_by(
        &before,
        &after,
        &moved_set(&before, &square),
        &t,
        0.0,
        Some(0.0),
    );
    assert_eq!(
        (count_of(&m, f, horizontal), count_of(&m, f, vertical)),
        (2, 2)
    );
    assert_eq!(m.sketch_def(f).unwrap().sketch.constraints.len(), n);
    assert_eq!(dof(&m, f), 2);
    // Half a turn keeps them.
    let before = after;
    choose(&mut ed, Tool::Rotate, &square);
    ed.clicks.push(at_point(0.0, 0.0));
    ed.transform_angle = PI;
    ed.apply_typed_transform(&mut m);
    moved_by(
        &before,
        &solved(&m, f),
        &moved_set(&before, &square),
        &Similarity::rotation(DVec2::ZERO, PI),
        0.0,
        Some(0.0),
    );
    assert_eq!(dof(&m, f), 2);
}

#[test]
fn other_turns_keep_right_angles_and_free_only_the_rotation() {
    let (mut m, mut ed, map) = setup();
    let f = ed.feature;
    let square = rectangle(&mut m, &mut ed, &map, (0.0, 0.0), (30.0, 20.0));
    let bottom = line_near(&m, f, &square, 15.0, 0.0);
    let left = line_near(&m, f, &square, 0.0, 10.0);
    add(
        &mut m,
        f,
        vec![
            Constraint::Length {
                line: bottom,
                value: 30.0 * MM,
            },
            Constraint::Length {
                line: left,
                value: 20.0 * MM,
            },
        ],
    );
    assert_eq!(dof(&m, f), 0);
    let before = solved(&m, f);
    choose(&mut ed, Tool::Rotate, &square);
    ed.clicks.push(at_point(0.0, 0.0));
    ed.transform_angle = 30f64.to_radians();
    ed.apply_typed_transform(&mut m);
    assert_eq!(
        ed.message.as_deref(),
        Some(
            "Rotated 4 items; the angle is now free, with 3 horizontal or vertical \
             constraints made parallel or perpendicular"
        )
    );
    let t = Similarity::rotation(DVec2::ZERO, 30f64.to_radians());
    moved_by(
        &before,
        &solved(&m, f),
        &moved_set(&before, &square),
        &t,
        1e-15,
        Some(0.0),
    );
    assert_eq!(count_of(&m, f, |c| horizontal(c) || vertical(c)), 0);
    assert_eq!(
        count_of(&m, f, |c| matches!(c, Constraint::Parallel { .. })),
        1
    );
    assert_eq!(
        count_of(&m, f, |c| matches!(c, Constraint::Perpendicular { .. })),
        2
    );
    assert_eq!(dof(&m, f), 1);
    assert_eq!(status(&m, f), SolveStatus::Solved);
    assert!(close(region_areas(&m, f)[0], 600.0 * MM * MM));
    assert_eq!(values(&m, f), vec![30.0 * MM, 20.0 * MM]);
}

#[test]
fn rotations_drop_frame_constraints_they_cannot_keep() {
    let (mut m, mut ed, map) = setup();
    let f = ed.feature;
    draw(&mut m, &mut ed, &map, Tool::Point, &[(10.0, 10.0)]);
    draw(&mut m, &mut ed, &map, Tool::Point, &[(40.0, 10.0)]);
    draw(&mut m, &mut ed, &map, Tool::Point, &[(40.0, 30.0)]);
    let pts = ids(&m, f, |g| matches!(g, Geometry::Point { .. }));
    let s = solved(&m, f);
    let (a, b, c) = (pts[0], pts[1], pts[2]);
    let p = |id: SketchId| s.point(id).unwrap();
    add(
        &mut m,
        f,
        vec![
            Constraint::HorizontalPoints { a, b },
            // Level the other way: becomes an aligned distance.
            Constraint::Distance {
                a,
                b,
                kind: DistanceKind::Horizontal,
                value: (p(b) - p(a)).x,
            },
            // Not level: removed.
            Constraint::Distance {
                a,
                b: c,
                kind: DistanceKind::Vertical,
                value: (p(c) - p(a)).y,
            },
        ],
    );
    choose(&mut ed, Tool::Rotate, &pts);
    ed.clicks.push(at_point(10.0, 10.0));
    ed.transform_angle = 0.4;
    ed.apply_typed_transform(&mut m);
    assert_eq!(
        ed.message.as_deref(),
        Some("Rotated 3 items; removed 2 constraints it could not keep")
    );
    let kept: Vec<Constraint<SketchId>> = m
        .sketch_def(f)
        .unwrap()
        .sketch
        .constraints
        .iter()
        .map(|c| c.constraint.clone())
        .collect();
    assert_eq!(
        kept,
        vec![Constraint::Distance {
            a,
            b,
            kind: DistanceKind::Aligned,
            value: (p(b) - p(a)).x,
        }]
    );
}

#[test]
fn scaling_scales_dimensions_and_areas_exactly() {
    let (mut m, mut ed, map) = setup();
    let f = ed.feature;
    let square = rectangle(&mut m, &mut ed, &map, (0.0, 0.0), (10.0, 10.0));
    draw(
        &mut m,
        &mut ed,
        &map,
        Tool::Circle,
        &[(5.0, 5.0), (7.0, 5.0)],
    );
    let circle = ids(&m, f, |g| matches!(g, Geometry::Circle { .. }))[0];
    let bottom = line_near(&m, f, &square, 5.0, 0.0);
    let left = line_near(&m, f, &square, 0.0, 5.0);
    let r = match solved(&m, f).entity(circle).unwrap().geometry {
        Geometry::Circle { radius, .. } => radius,
        _ => unreachable!(),
    };
    add(
        &mut m,
        f,
        vec![
            Constraint::Length {
                line: bottom,
                value: 10.0 * MM,
            },
            Constraint::Length {
                line: left,
                value: 10.0 * MM,
            },
            Constraint::Radius {
                curve: circle,
                value: r,
            },
        ],
    );
    let area = region_areas(&m, f)[0];
    let mut all = square.clone();
    all.push(circle);
    choose(&mut ed, Tool::Scale, &all);
    click(&mut m, &mut ed, &map, 0.0, 0.0);
    assert_eq!(
        ed.hint(),
        "Click a reference point, or type the factor and Apply"
    );
    ed.transform_factor = 3.0;
    ed.apply_typed_transform(&mut m);
    assert_eq!(ed.message.as_deref(), Some("Scaled 5 items"));
    assert_eq!(values(&m, f), vec![30.0 * MM, 30.0 * MM, 3.0 * r]);
    let a = region_areas(&m, f)[0];
    assert!((a - 9.0 * area).abs() <= 1e-12 * a, "{a} vs {}", 9.0 * area);
    assert_eq!(dof(&m, f), 2);
}

#[test]
fn clicked_rotations_and_scales_use_three_points() {
    let (mut m, mut ed, map) = setup();
    let f = ed.feature;
    let square = rectangle(&mut m, &mut ed, &map, (10.0, 10.0), (20.0, 20.0));
    draw(&mut m, &mut ed, &map, Tool::Point, &[(10.0, 50.0)]);
    let before = solved(&m, f);
    choose(&mut ed, Tool::Rotate, &square);
    // About the lower left corner, from the lower right corner to the
    // point straight above: a quarter turn, exact from the snapped points.
    click(&mut m, &mut ed, &map, 10.0, 10.0);
    click(&mut m, &mut ed, &map, 30.0, 10.0);
    assert_eq!(
        ed.hint(),
        "Click where the reference point turns to; Esc starts again"
    );
    hover(&m, &mut ed, &map, 10.0, 50.0);
    assert!(ed.curve_edit.clone().unwrap().plan.is_ok());
    assert_eq!(ed.rubber_band().len(), 1);
    click(&mut m, &mut ed, &map, 10.0, 50.0);
    assert_eq!(ed.message.as_deref(), Some("Rotated 4 items"));
    let corner = DVec2::new(10.0, 10.0) * MM;
    moved_by(
        &before,
        &solved(&m, f),
        &moved_set(&before, &square),
        &Similarity::rotation(corner, PI / 2.0),
        1e-15,
        Some(0.0),
    );
    // Scale about the same corner, from (10, 30) to the point: × 2.
    let before = solved(&m, f);
    choose(&mut ed, Tool::Scale, &square);
    click(&mut m, &mut ed, &map, 10.0, 10.0);
    click(&mut m, &mut ed, &map, 10.0, 30.0);
    click(&mut m, &mut ed, &map, 10.0, 50.0);
    assert_eq!(ed.message.as_deref(), Some("Scaled 4 items"));
    moved_by(
        &before,
        &solved(&m, f),
        &moved_set(&before, &square),
        &Similarity::scaling(corner, 2.0),
        1e-15,
        Some(0.0),
    );
    assert!(close(region_areas(&m, f)[0], 1600.0 * MM * MM));
}

#[test]
fn copies_get_the_constraints_among_the_selection_only() {
    let (mut m, mut ed, map) = setup();
    let f = ed.feature;
    let square = rectangle(&mut m, &mut ed, &map, (0.0, 0.0), (10.0, 10.0));
    let bottom = line_near(&m, f, &square, 5.0, 0.0);
    let left = line_near(&m, f, &square, 0.0, 5.0);
    add(
        &mut m,
        f,
        vec![
            Constraint::Length {
                line: bottom,
                value: 10.0 * MM,
            },
            Constraint::Length {
                line: left,
                value: 10.0 * MM,
            },
        ],
    );
    let n = m.sketch_def(f).unwrap().sketch.constraints.len();
    let before = solved(&m, f);
    choose(&mut ed, Tool::Move, &square);
    ed.transform_copy = true;
    ed.transform_offset = DVec2::new(30.0, 0.0) * MM;
    ed.apply_typed_transform(&mut m);
    assert_eq!(ed.message.as_deref(), Some("Copied 4 items"));
    // The originals stay; the copies are selected, fixed where they went
    // and fully defined like them.
    moved_by(
        &before,
        &solved(&m, f),
        &BTreeSet::new(),
        &Similarity::translation(DVec2::ZERO),
        0.0,
        Some(0.0),
    );
    assert_eq!(lines(&m, f).len(), 8);
    assert_eq!(ed.selection.len(), 4);
    assert!(ed.selection.iter().all(|id| !square.contains(id)));
    assert_eq!(m.sketch_def(f).unwrap().sketch.constraints.len(), 2 * n);
    assert_eq!(dof(&m, f), 0);
    let areas = region_areas(&m, f);
    assert_eq!(areas.len(), 2);
    assert!(areas.iter().all(|a| close(*a, 100.0 * MM * MM)));
    // A turned copy keeps its right angles and is free to turn.
    choose(&mut ed, Tool::Rotate, &square);
    ed.transform_copy = true;
    ed.clicks.push(at_point(0.0, 0.0));
    ed.transform_angle = 0.5;
    ed.apply_typed_transform(&mut m);
    assert_eq!(
        ed.message.as_deref(),
        Some(
            "Copied 4 items; the angle is now free, with 3 horizontal or vertical \
             constraints made parallel or perpendicular"
        )
    );
    assert_eq!(dof(&m, f), 1);
    // A copy onto the original is refused.
    choose(&mut ed, Tool::Move, &square);
    ed.transform_copy = true;
    ed.transform_offset = DVec2::ZERO;
    ed.apply_typed_transform(&mut m);
    assert_eq!(
        ed.message.as_deref(),
        Some("Not done: The copy would lie on the original")
    );
}

#[test]
fn placed_labels_turn_with_their_dimensions() {
    let (mut m, mut ed, map) = setup();
    let f = ed.feature;
    draw(
        &mut m,
        &mut ed,
        &map,
        Tool::Line,
        &[(10.0, 10.0), (40.0, 10.0)],
    );
    let line = lines(&m, f)[0];
    add(
        &mut m,
        f,
        vec![Constraint::Length {
            line,
            value: 30.0 * MM,
        }],
    );
    let dim = m
        .sketch_def(f)
        .unwrap()
        .sketch
        .constraints
        .last()
        .unwrap()
        .id;
    m.edit_sketch(
        f,
        vec![SketchEdit::PlaceLabel {
            dimension: dim,
            offset: Some(DVec2::new(0.0, 8.0) * MM),
        }],
    )
    .unwrap();
    choose(&mut ed, Tool::Rotate, &[line]);
    ed.clicks.push(at_point(10.0, 10.0));
    ed.apply_typed_transform(&mut m);
    let label = m.sketch_def(f).unwrap().label(dim).unwrap();
    assert_eq!(label, DVec2::new(-8.0, 0.0) * MM);
    choose(&mut ed, Tool::Scale, &[line]);
    ed.clicks.push(at_point(10.0, 10.0));
    ed.transform_factor = 0.5;
    ed.apply_typed_transform(&mut m);
    assert_eq!(
        m.sketch_def(f).unwrap().label(dim).unwrap(),
        DVec2::new(-4.0, 0.0) * MM
    );
    assert_eq!(values(&m, f), vec![15.0 * MM]);
}

#[test]
fn curves_of_every_kind_turn_and_scale() {
    let (mut m, mut ed, map) = setup();
    let f = ed.feature;
    draw(
        &mut m,
        &mut ed,
        &map,
        Tool::Ellipse,
        &[(30.0, 30.0), (45.0, 30.0), (30.0, 38.0)],
    );
    draw(
        &mut m,
        &mut ed,
        &map,
        Tool::Conic,
        &[(10.0, 5.0), (70.0, 5.0), (30.0, 45.0), (45.0, 25.0)],
    );
    draw(
        &mut m,
        &mut ed,
        &map,
        Tool::ControlSpline,
        &[(10.0, 60.0), (20.0, 70.0), (35.0, 55.0), (50.0, 65.0)],
    );
    ed.set_tool(Tool::Select);
    let curves = ids(&m, f, |g| g.is_curve());
    let rho = |m: &Model| {
        solved(m, f)
            .entities
            .iter()
            .find_map(|e| match e.geometry {
                Geometry::Conic { rho, .. } => Some(rho),
                _ => None,
            })
            .unwrap()
    };
    let rho0 = rho(&m);
    let n = m.sketch_def(f).unwrap().sketch.constraints.len();
    let dof0 = dof(&m, f);
    let before = solved(&m, f);
    choose(&mut ed, Tool::Rotate, &curves);
    ed.clicks.push(at_point(20.0, 20.0));
    ed.transform_angle = 1.0;
    ed.apply_typed_transform(&mut m);
    assert!(ed.message.as_deref().unwrap().starts_with("Rotated"));
    let t = Similarity::rotation(DVec2::new(20.0, 20.0) * MM, 1.0);
    let after = solved(&m, f);
    moved_by(
        &before,
        &after,
        &moved_set(&before, &curves),
        &t,
        1e-15,
        Some(0.0),
    );
    assert_eq!(rho(&m), rho0);
    assert_eq!(dof(&m, f), dof0);
    assert_eq!(redundant(&m, f), 0);
    let before = after;
    choose(&mut ed, Tool::Scale, &curves);
    ed.clicks.push(at_point(0.0, 50.0));
    ed.transform_factor = 1.5;
    ed.apply_typed_transform(&mut m);
    assert!(ed.message.as_deref().unwrap().starts_with("Scaled"));
    let t = Similarity::scaling(DVec2::new(0.0, 50.0) * MM, 1.5);
    moved_by(
        &before,
        &solved(&m, f),
        &moved_set(&before, &curves),
        &t,
        1e-15,
        Some(0.0),
    );
    assert_eq!(rho(&m), rho0);
    assert_eq!(m.sketch_def(f).unwrap().sketch.constraints.len(), n);
}

/// The scene of the fuzz: rectangles and circles at random places, some
/// connected by lines between their corners, the first rectangle perhaps
/// mirrored, the first circle perhaps patterned about a centre, and
/// distances between random points.
#[derive(Clone, Debug)]
struct Scene {
    rects: Vec<((f64, f64), (f64, f64))>,
    circles: Vec<((f64, f64), f64)>,
    links: Vec<(usize, usize)>,
    mirror: bool,
    pattern: bool,
    dims: Vec<(usize, usize)>,
}

fn scene() -> impl Strategy<Value = Scene> {
    let rect = (
        (-90.0..40.0f64, -80.0..30.0f64),
        (8.0..40.0f64, 8.0..40.0f64),
    );
    let circle = ((-90.0..90.0f64, -80.0..80.0f64), 3.0..15.0f64);
    (
        prop::collection::vec(rect, 1..3),
        prop::collection::vec(circle, 0..3),
        prop::collection::vec((0usize..8, 0usize..8), 0..2),
        any::<bool>(),
        any::<bool>(),
        prop::collection::vec((0usize..64, 0usize..64), 0..3),
    )
        .prop_map(|(rects, circles, links, mirror, pattern, dims)| Scene {
            rects,
            circles,
            links,
            mirror,
            pattern,
            dims,
        })
}

#[derive(Clone, Copy, Debug)]
enum Kind {
    Move(f64, f64),
    Rotate(f64, f64, f64),
    Scale(f64, f64, f64),
}

fn kind() -> impl Strategy<Value = Kind> {
    prop_oneof![
        (-50.0..50.0f64, -50.0..50.0f64).prop_map(|(x, y)| Kind::Move(x, y)),
        (-50.0..50.0f64, -50.0..50.0f64, -3.0..3.0f64).prop_map(|(x, y, a)| Kind::Rotate(x, y, a)),
        (-50.0..50.0f64, -50.0..50.0f64, 0.3..3.0f64).prop_map(|(x, y, s)| Kind::Scale(x, y, s)),
        (-50.0..50.0f64, -50.0..50.0f64, 0i32..4).prop_map(|(x, y, k)| Kind::Rotate(
            x,
            y,
            f64::from(k) * PI / 2.0
        )),
    ]
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(200))]

    /// Random transforms of random selections either are refused with the
    /// sketch untouched, or land the moved geometry exactly where the
    /// similarity puts it, leave the sketch solved without new redundancy,
    /// and (detached) leave everything else where it was.
    #[test]
    fn random_transforms_land_exactly_or_are_refused(
        scene in scene(),
        pick in prop::collection::vec(any::<bool>(), 24),
        kind in kind(),
        keep in any::<bool>(),
        copy in any::<bool>(),
    ) {
        let (mut m, mut ed, map) = setup();
        let f = ed.feature;
        let mut first = Vec::new();
        for &((x, y), (w, h)) in &scene.rects {
            let r = rectangle(&mut m, &mut ed, &map, (x, y), (w, h));
            if first.is_empty() {
                first = r;
            }
        }
        for &((x, y), r) in &scene.circles {
            draw(&mut m, &mut ed, &map, Tool::Circle, &[(x, y), (x + r, y)]);
        }
        let corners: Vec<DVec2> = scene
            .rects
            .iter()
            .flat_map(|&((x, y), (w, h))| [(x, y), (x + w, y), (x + w, y + h), (x, y + h)])
            .map(|(x, y)| DVec2::new(x, y))
            .collect();
        for &(a, b) in &scene.links {
            let (a, b) = (corners[a % corners.len()], corners[b % corners.len()]);
            if (a - b).length() > 1.0 {
                draw(&mut m, &mut ed, &map, Tool::Line, &[(a.x, a.y), (b.x, b.y)]);
            }
        }
        if scene.mirror {
            draw(&mut m, &mut ed, &map, Tool::Line, &[(95.0, -90.0), (95.0, 90.0)]);
            choose(&mut ed, Tool::Mirror, &first);
            click(&mut m, &mut ed, &map, 95.0, 0.0);
        }
        if let (true, Some(circle)) = (
            scene.pattern,
            ids(&m, f, |g| matches!(g, Geometry::Circle { .. })).first().copied(),
        ) {
            choose(&mut ed, Tool::CircularPattern, &[circle]);
            click(&mut m, &mut ed, &map, 110.0, -70.0);
        }
        let points = ids(&m, f, |g| matches!(g, Geometry::Point { .. }));
        for &(a, b) in &scene.dims {
            let (a, b) = (points[a % points.len()], points[b % points.len()]);
            let mut c = Constraint::Distance { a, b, kind: DistanceKind::Aligned, value: 1.0 };
            if let Ok(v) = opendatum_sketch_solver::measure(&solved(&m, f), &c) {
                if a != b && v > 1e-3 {
                    c.set_value(v);
                    let mut batch = Batch::default();
                    batch.constraint(c);
                    let _ = m.edit_sketch_checked(f, batch.edits);
                }
            }
        }
        ed.set_tool(Tool::Select);
        prop_assume!(!matches!(status(&m, f), SolveStatus::Failed(_)));
        let curves = ids(&m, f, |g| g.is_curve());
        let selection: Vec<SketchId> = curves
            .iter()
            .zip(&pick)
            .filter(|(_, p)| **p)
            .map(|(c, _)| *c)
            .collect();
        prop_assume!(!selection.is_empty());
        let before_def = m.sketch_def(f).unwrap().clone();
        let before = solved(&m, f);
        let redundant0 = redundant(&m, f);
        let (tool, t) = match kind {
            Kind::Move(x, y) => (Tool::Move, Similarity::translation(DVec2::new(x, y) * MM)),
            Kind::Rotate(x, y, a) => (Tool::Rotate, Similarity::rotation(DVec2::new(x, y) * MM, a)),
            Kind::Scale(x, y, s) => (Tool::Scale, Similarity::scaling(DVec2::new(x, y) * MM, s)),
        };
        choose(&mut ed, tool, &selection);
        ed.transform_keep = keep;
        ed.transform_copy = copy;
        ed.transform_offset = t.shift;
        ed.transform_angle = t.angle;
        ed.transform_factor = t.scale;
        if tool != Tool::Move {
            ed.clicks.push(Snap { at: t.centre, target: SnapTarget::Free });
        }
        ed.apply_typed_transform(&mut m);
        let message = ed.message.clone().unwrap_or_default();
        let sources = mirror::sources(&before, &selection);
        let moved = moved_set(&before, &sources);
        if message.starts_with("Not done") {
            prop_assert_eq!(m.sketch_def(f).unwrap(), &before_def, "{}", message);
            // Only the identity, relations the rest cannot follow, or a
            // curve joined to the selection collapsing are reasons to
            // refuse.
            let tied = before.constraints.iter().any(|c| {
                let refs = c.constraint.referenced();
                c.constraint.value().is_none()
                    && refs.iter().any(|k| moved.contains(k))
                    && refs.iter().any(|k| !moved.contains(k))
            }) || before.entities.iter().any(|e| {
                !moved.contains(&e.id)
                    && e.geometry.referenced_points().iter().any(|p| moved.contains(p))
            });
            let collapse = message.contains("collapse a curve joined to the selection");
            prop_assert!(
                t.is_identity() || (!copy && tied && (keep || collapse)),
                "{}",
                message
            );
            return Ok(());
        }
        prop_assert!(!matches!(status(&m, f), SolveStatus::Failed(_)), "{}", message);
        prop_assert!(redundant(&m, f) <= redundant0);
        let after = solved(&m, f);
        // Rounding: coordinates of up to ~0.2 m through a turn and a factor.
        let tol = 1e-15 * 64.0;
        if copy {
            moved_by(&before, &after, &BTreeSet::new(), &t, 0.0, Some(tol));
            prop_assert_eq!(ed.selection.len(), sources.len());
            for (src, copy) in sources.iter().zip(&ed.selection) {
                let (a, b) = (before.entity(*src).unwrap(), after.entity(*copy).unwrap());
                let (pa, pb) = (a.geometry.referenced_points(), b.geometry.referenced_points());
                prop_assert_eq!(pa.len(), pb.len());
                for (p, q) in pa.into_iter().zip(pb) {
                    let want = t.point(before.point(p).unwrap());
                    prop_assert!((after.point(q).unwrap() - want).length() <= tol);
                }
            }
        } else {
            moved_by(&before, &after, &moved, &t, tol, (!keep).then_some(tol));
            if !keep {
                prop_assert!(!message.contains("made parallel") || tool == Tool::Rotate);
            }
        }
        if !message.contains("removed") && !message.contains("free") && !copy {
            prop_assert_eq!(dof(&m, f), before_dof(&before_def));
        }
    }
}

/// The degrees of freedom of a sketch definition, solved on its own.
fn before_dof(def: &opendatum_core_model::SketchDef) -> usize {
    opendatum_sketch_solver::solve(&def.sketch, &Default::default()).dof
}
