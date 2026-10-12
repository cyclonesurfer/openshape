//! Sketch mirror (design note 0007 §9), driven through the editor: images
//! of a selection or of curves clicked one at a time, joined where the
//! sources meet the mirror line, following their sources, with no degrees
//! of freedom or redundancy added; refusals; and a fuzz of random mirrors.

mod common;

use common::*;
use opendatum_app::sketch_mode::{Snap, SnapTarget, Tool};
use opendatum_app::sketch_ops::mirror::{self, Axis};
use opendatum_app::sketch_ops::trim::Preview;
use opendatum_app::sketch_ops::{shapes, Batch};
use opendatum_app::Model;
use opendatum_core_model::sketch::{Constraint, Geometry};
use opendatum_core_model::{DVec2, FeatureId, PlaneRef, SketchId};
use opendatum_sketch_solver::{ResolvedCurve, SolveStatus};
use proptest::prelude::*;

fn redundant(m: &Model, f: FeatureId) -> usize {
    m.sketch_result(f).unwrap().solution.redundant.len()
}

fn pos(m: &Model, f: FeatureId, id: SketchId) -> DVec2 {
    m.sketch_result(f)
        .unwrap()
        .solution
        .sketch
        .point(id)
        .expect("a point")
}

/// The mirror constraints of a sketch: (axis, pairs).
fn mirrors(m: &Model, f: FeatureId) -> Vec<(SketchId, Vec<(SketchId, SketchId)>)> {
    m.sketch_def(f)
        .unwrap()
        .sketch
        .constraints
        .iter()
        .filter_map(|c| match &c.constraint {
            Constraint::Mirror { axis, pairs } => Some((*axis, pairs.clone())),
            _ => None,
        })
        .collect()
}

/// A fixed vertical construction line along x = 0 from y = −50 to 50 mm,
/// and half a 40 × 20 mm rectangle to its right: lines A→B, B→C, C→D with
/// A (0, −10) and D (0, 10) held on the line by point-on-line constraints.
/// Returns the axis and the three lines.
fn half_rectangle(m: &mut Model, f: FeatureId) -> (SketchId, [SketchId; 3]) {
    let mut b = Batch::default();
    let (bottom, top) = (DVec2::new(0.0, -50.0) * MM, DVec2::new(0.0, 50.0) * MM);
    let axis = b.line(bottom, top, true);
    b.constraint(Constraint::Fix {
        point: axis.start,
        at: bottom,
    });
    b.constraint(Constraint::Fix {
        point: axis.end,
        at: top,
    });
    let pts = [(0.0, -10.0), (20.0, -10.0), (20.0, 10.0), (0.0, 10.0)]
        .map(|(x, y)| DVec2::new(x, y) * MM);
    let ab = b.line(pts[0], pts[1], false);
    let bc = b.line(pts[1], pts[2], false);
    let cd = b.line(pts[2], pts[3], false);
    b.coincident(ab.end, bc.start);
    b.coincident(bc.end, cd.start);
    b.constraint(Constraint::Horizontal { line: ab.curve });
    b.constraint(Constraint::Vertical { line: bc.curve });
    b.constraint(Constraint::Horizontal { line: cd.curve });
    b.constraint(Constraint::PointOnCurve {
        point: ab.start,
        curve: axis.curve,
    });
    b.constraint(Constraint::PointOnCurve {
        point: cd.end,
        curve: axis.curve,
    });
    m.edit_sketch(f, b.edits).unwrap();
    // The batch ids resolve to the edit's ids: find the curves by position.
    let lines = curves(m, f, true);
    let find = |a: DVec2, z: DVec2| {
        lines
            .iter()
            .find(|(_, c)| matches!(c, ResolvedCurve::Line { start, end } if (*start - a * MM).length() < 1e-9 && (*end - z * MM).length() < 1e-9))
            .expect("the line")
            .0
    };
    (
        find(DVec2::new(0.0, -50.0), DVec2::new(0.0, 50.0)),
        [
            find(DVec2::new(0.0, -10.0), DVec2::new(20.0, -10.0)),
            find(DVec2::new(20.0, -10.0), DVec2::new(20.0, 10.0)),
            find(DVec2::new(20.0, 10.0), DVec2::new(0.0, 10.0)),
        ],
    )
}

fn points(m: &Model, f: FeatureId) -> usize {
    count(m, f, |g| matches!(g, Geometry::Point { .. }))
}

fn on_axis(m: &Model, f: FeatureId) -> usize {
    constraints(m, f, |c| matches!(c, Constraint::PointOnCurve { .. }))
}

#[test]
fn a_selection_is_mirrored_about_the_line_clicked() {
    let (mut m, mut ed, map) = setup();
    let f = ed.feature;
    let (axis, lines) = half_rectangle(&mut m, f);
    let (before, n_points) = (dof(&m, f), points(&m, f));
    assert_eq!(before, 3);
    ed.selection = lines.to_vec();
    ed.set_tool(Tool::Mirror);
    assert_eq!(
        ed.hint(),
        "Click the line to mirror the 3 selected items about"
    );
    // Hovering the line previews the images.
    hover(&m, &mut ed, &map, 0.0, 30.0);
    let plan = ed
        .curve_edit
        .clone()
        .expect("the line")
        .plan
        .expect("a plan");
    let Preview::Shapes {
        curves: shown,
        points: lone,
    } = plan.preview
    else {
        panic!("{:?}", plan.preview)
    };
    assert_eq!(shown.len(), 3);
    assert!(lone.is_empty());
    assert!(matches!(shown[1], ResolvedCurve::Line { start, end }
        if (start - DVec2::new(-20.0, -10.0) * MM).length() < 1e-9
            && (end - DVec2::new(-20.0, 10.0) * MM).length() < 1e-9));
    click(&mut m, &mut ed, &map, 0.0, 30.0);
    assert_eq!(ed.message.as_deref(), Some("Mirrored 3 items"));
    assert!(ed.selection.is_empty());
    // The images share A and D (on the line); the joints at B and C are
    // pairs of coincident points, so their images are four new points.
    assert_eq!(points(&m, f), n_points + 4);
    let def = &m.sketch_def(f).unwrap().sketch;
    assert_eq!(shared_points(def, axis).len(), 2);
    let mirrored = mirrors(&m, f);
    assert_eq!(mirrored.len(), 1);
    assert_eq!(mirrored[0].1.len(), 3);
    // The mirror holds A and D on the line now.
    assert_eq!(on_axis(&m, f), 0);
    assert_eq!(status(&m, f), SolveStatus::Solved);
    assert_eq!(redundant(&m, f), 0);
    assert_eq!(dof(&m, f), before);
    // One closed 40 × 20 mm rectangle.
    let areas = region_areas(&m, f);
    assert_eq!(areas.len(), 1);
    assert!((areas[0] - 800.0 * MM * MM).abs() < 1e-12, "{areas:?}");
    let x = extrude(&mut m, &map, &[(-5.0, 0.0)], "10 mm");
    assert!((volume(&m, x) - 8e-6).abs() < 1e-15);
    // The tool stays on the line for more.
    assert!(ed.mirror_axis.is_some());
}

#[test]
fn the_images_follow_their_sources() {
    let (mut m, mut ed, map) = setup();
    let f = ed.feature;
    let (_, lines) = half_rectangle(&mut m, f);
    ed.selection = lines.to_vec();
    ed.set_tool(Tool::Mirror);
    click(&mut m, &mut ed, &map, 0.0, 30.0);
    // Dimension the source's top: the image's top follows.
    let mut b = Batch::default();
    b.constraint(Constraint::Length {
        line: lines[2],
        value: 0.03,
    });
    b.constraint(Constraint::Length {
        line: lines[1],
        value: 0.02,
    });
    m.edit_sketch(f, b.edits).unwrap();
    assert_eq!(status(&m, f), SolveStatus::Solved);
    let areas = region_areas(&m, f);
    assert!(
        (areas[0] - 60.0 * 20.0 * MM * MM).abs() < 1e-12,
        "{areas:?}"
    );
    // Every image point is the reflection of its source.
    let (axis, pairs) = mirrors(&m, f).remove(0);
    let sketch = m.sketch_result(f).unwrap().solution.sketch.clone();
    let line = Axis::of(&sketch, axis).unwrap();
    for (src, img) in pairs {
        let (Geometry::Line { start: s0, end: s1 }, Geometry::Line { start: i0, end: i1 }) = (
            sketch.entity(src).unwrap().geometry.clone(),
            sketch.entity(img).unwrap().geometry.clone(),
        ) else {
            panic!("lines")
        };
        for (s, i) in [(s0, i0), (s1, i1)] {
            let want = line.reflect(pos(&m, f, s));
            assert!((pos(&m, f, i) - want).length() < 1e-12);
        }
    }
}

#[test]
fn curves_clicked_one_at_a_time_join_up() {
    let (mut m, mut ed, map) = setup();
    let f = ed.feature;
    let (axis, _) = half_rectangle(&mut m, f);
    let (before, n_points) = (dof(&m, f), points(&m, f));
    ed.set_tool(Tool::Mirror);
    assert_eq!(
        ed.hint(),
        "Click the mirror line, then the curves to mirror"
    );
    click(&mut m, &mut ed, &map, 0.0, 30.0);
    assert_eq!(ed.message.as_deref(), Some("Mirror line chosen"));
    // The line itself is not mirrored.
    hover(&m, &mut ed, &map, 0.0, 30.0);
    assert_eq!(
        ed.curve_edit.clone().unwrap().plan.unwrap_err(),
        "This is the mirror line"
    );
    for (x, y) in [(10.0, -10.0), (20.0, 0.0), (10.0, 10.0)] {
        click(&mut m, &mut ed, &map, x, y);
        assert_eq!(ed.message.as_deref(), Some("Mirrored 1 item"));
    }
    // One constraint about the line, the images sharing A and D and
    // joining each other through the images of the joints.
    let mirrored = mirrors(&m, f);
    assert_eq!(mirrored.len(), 1);
    assert_eq!(mirrored[0].1.len(), 3);
    assert_eq!(points(&m, f), n_points + 4);
    let def = &m.sketch_def(f).unwrap().sketch;
    assert_eq!(shared_points(def, axis).len(), 2);
    assert_eq!(status(&m, f), SolveStatus::Solved);
    assert_eq!(redundant(&m, f), 0);
    assert_eq!(dof(&m, f), before);
    assert_eq!(region_areas(&m, f).len(), 1);
    // Clicking a source or an image again does nothing new.
    click(&mut m, &mut ed, &map, -20.0, 0.0);
    assert_eq!(
        ed.message.as_deref(),
        Some("Already mirrored about this line, or on it (their own mirror image)")
    );
    // Escape chooses another line, then leaves the tool.
    ed.escape();
    assert!(ed.mirror_axis.is_none());
    assert_eq!(ed.tool, Tool::Mirror);
    ed.escape();
    assert_eq!(ed.tool, Tool::Select);
}

#[test]
fn deleting_an_image_keeps_the_others_mirrored() {
    let (mut m, mut ed, map) = setup();
    let f = ed.feature;
    let (_, lines) = half_rectangle(&mut m, f);
    ed.selection = lines.to_vec();
    ed.set_tool(Tool::Mirror);
    click(&mut m, &mut ed, &map, 0.0, 30.0);
    let image = mirrors(&m, f)[0].1[1].1;
    ed.set_tool(Tool::Select);
    ed.selection = vec![image];
    ed.delete_selection(&mut m);
    let mirrored = mirrors(&m, f);
    assert_eq!(mirrored.len(), 1);
    assert_eq!(mirrored[0].1.len(), 2);
    assert_eq!(status(&m, f), SolveStatus::Solved);
    assert_eq!(redundant(&m, f), 0);
    // Deleting a source frees its image the same way.
    ed.selection = vec![lines[0]];
    ed.delete_selection(&mut m);
    assert_eq!(mirrors(&m, f)[0].1.len(), 1);
    assert_eq!(status(&m, f), SolveStatus::Solved);
}

#[test]
fn arcs_circles_and_points_mirror_about_a_slanted_line() {
    let (mut m, mut ed, map) = setup();
    let f = ed.feature;
    let mut b = Batch::default();
    // The line y = x / 2 through the origin.
    b.line(
        DVec2::new(-40.0, -20.0) * MM,
        DVec2::new(40.0, 20.0) * MM,
        true,
    );
    b.arc(
        DVec2::new(10.0, 30.0) * MM,
        DVec2::new(20.0, 30.0) * MM,
        DVec2::new(10.0, 40.0) * MM,
        false,
    );
    let center = b.point(DVec2::new(-20.0, 30.0) * MM);
    b.entity(
        Geometry::Circle {
            center,
            radius: 5.0 * MM,
        },
        false,
    );
    b.point(DVec2::new(30.0, 40.0) * MM);
    m.edit_sketch(f, b.edits).unwrap();
    let before = dof(&m, f);
    let sketch = m.sketch_def(f).unwrap().sketch.clone();
    // Everything but the line and its end points (whose images would be
    // themselves): the arc and circle with their points, and a lone point.
    let of_line = sketch
        .entities
        .iter()
        .find(|e| matches!(e.geometry, Geometry::Line { .. }))
        .unwrap()
        .geometry
        .referenced_points();
    let all: Vec<SketchId> = sketch
        .entities
        .iter()
        .filter(|e| !matches!(e.geometry, Geometry::Line { .. }) && !of_line.contains(&e.id))
        .map(|e| e.id)
        .collect();
    ed.selection = all;
    ed.set_tool(Tool::Mirror);
    click(&mut m, &mut ed, &map, 30.0, 15.0);
    // The arc, the circle and the lone point (the arc's and circle's own
    // points come with them).
    assert_eq!(ed.message.as_deref(), Some("Mirrored 3 items"));
    assert_eq!(status(&m, f), SolveStatus::Solved);
    assert_eq!(redundant(&m, f), 0);
    assert_eq!(dof(&m, f), before);
    let line = Axis {
        at: DVec2::ZERO,
        dir: DVec2::new(2.0, 1.0).normalize(),
    };
    let all = curves(&m, f, true);
    let arcs: Vec<&ResolvedCurve> = all
        .iter()
        .map(|(_, c)| c)
        .filter(|c| matches!(c, ResolvedCurve::Arc { .. }))
        .collect();
    assert_eq!(arcs.len(), 2);
    let expected = line.reflect_curve(arcs[0]);
    let (
        ResolvedCurve::Arc {
            center: c1,
            radius: r1,
            start_angle: a1,
            sweep_angle: s1,
        },
        ResolvedCurve::Arc {
            center: c2,
            radius: r2,
            start_angle: a2,
            sweep_angle: s2,
        },
    ) = (expected, arcs[1].clone())
    else {
        unreachable!()
    };
    assert!((c1 - c2).length() < 1e-9);
    assert!((r1 - r2).abs() < 1e-9);
    assert!((s1 - s2).abs() < 1e-9);
    assert!((DVec2::from_angle(a1) - DVec2::from_angle(a2)).length() < 1e-9);
    let circles: Vec<DVec2> = all
        .iter()
        .filter_map(|(_, c)| match c {
            ResolvedCurve::Circle { center, .. } => Some(*center),
            _ => None,
        })
        .collect();
    assert_eq!(circles.len(), 2);
    assert!((circles[1] - line.reflect(DVec2::new(-20.0, 30.0) * MM)).length() < 1e-9);
    let lone = m
        .sketch_result(f)
        .unwrap()
        .solution
        .sketch
        .entities
        .iter()
        .filter_map(|e| match e.geometry {
            Geometry::Point { at } => Some(at),
            _ => None,
        })
        .any(|at| (at - line.reflect(DVec2::new(30.0, 40.0) * MM)).length() < 1e-9);
    assert!(lone);
}

#[test]
fn mirrors_that_cannot_be_made_are_refused_or_skipped() {
    let (mut m, mut ed, map) = setup();
    let f = ed.feature;
    let (axis, lines) = half_rectangle(&mut m, f);
    // A circle centred on the line is its own image.
    let mut b = Batch::default();
    let center = b.point(DVec2::new(0.0, 30.0) * MM);
    b.entity(
        Geometry::Circle {
            center,
            radius: 5.0 * MM,
        },
        false,
    );
    m.edit_sketch(f, b.edits).unwrap();
    // A circle is not a mirror line.
    ed.selection = lines.to_vec();
    ed.set_tool(Tool::Mirror);
    hover(&m, &mut ed, &map, 5.0, 30.0);
    assert_eq!(
        ed.curve_edit.clone().unwrap().plan.unwrap_err(),
        "The mirror line must be a line"
    );
    // The line alone, or what lies on it, mirrors to nothing.
    let circle = curves(&m, f, true)
        .into_iter()
        .find(|(_, c)| matches!(c, ResolvedCurve::Circle { .. }))
        .unwrap()
        .0;
    let sketch = m.sketch_def(f).unwrap().sketch.clone();
    let res = m.sketch_result(f).unwrap();
    let shown: Vec<(SketchId, ResolvedCurve)> = res
        .curves
        .iter()
        .map(|(id, _, c)| (*id, c.clone()))
        .collect();
    for sources in [vec![axis], vec![circle], vec![]] {
        let why = mirror::mirror(&sketch, &shown, &sources, axis, true).unwrap_err();
        assert!(
            why.starts_with("Already mirrored") || why.starts_with("Nothing to mirror"),
            "{why}"
        );
    }
    // With the circle in the selection, the rest is mirrored and the
    // circle is counted as skipped.
    ed.selection = [lines.to_vec(), vec![circle]].concat();
    ed.set_tool(Tool::Mirror);
    click(&mut m, &mut ed, &map, 0.0, -30.0);
    assert_eq!(
        ed.message.as_deref(),
        Some("Mirrored 3 items; skipped 1 already mirrored about this line or lying on it")
    );
}

#[test]
fn a_point_held_at_the_lines_midpoint_gets_its_own_image() {
    // Sharing D would repeat its midpoint constraint (which already holds
    // it on the line); its image gets a point of its own on the same spot.
    let (mut m, ed, _) = setup();
    let f = ed.feature;
    let mut b = Batch::default();
    let axis = b.line(
        DVec2::new(0.0, -40.0) * MM,
        DVec2::new(0.0, 40.0) * MM,
        true,
    );
    let l = b.line(DVec2::new(20.0, 20.0) * MM, DVec2::ZERO, false);
    b.constraint(Constraint::Midpoint {
        point: l.end,
        line: axis.curve,
    });
    m.edit_sketch(f, b.edits).unwrap();
    let before = (dof(&m, f), points(&m, f));
    let sketch = m.sketch_def(f).unwrap().sketch.clone();
    let res = m.sketch_result(f).unwrap();
    let shown: Vec<(SketchId, ResolvedCurve)> = res
        .curves
        .iter()
        .map(|(id, _, c)| (*id, c.clone()))
        .collect();
    let axis_id = sketch.entities.iter().find(|e| e.construction).unwrap().id;
    let line_id = sketch
        .entities
        .iter()
        .find(|e| e.geometry.is_curve() && !e.construction)
        .unwrap()
        .id;
    let made = mirror::mirror(&sketch, &shown, &[line_id], axis_id, true).unwrap();
    m.edit_sketch_checked(f, made.plan.edits).unwrap();
    assert_eq!(status(&m, f), SolveStatus::Solved);
    assert_eq!(redundant(&m, f), 0);
    assert_eq!(dof(&m, f), before.0);
    assert_eq!(points(&m, f), before.1 + 2);
}

// ---------------------------------------------------------------------------
// Fuzz
// ---------------------------------------------------------------------------

#[derive(Clone, Debug)]
enum Shape {
    Rectangle(DVec2, DVec2),
    Polygon(DVec2, DVec2, u32),
    Slot(DVec2, DVec2, DVec2),
    Circle(DVec2, f64),
    Arc(DVec2, f64, f64, f64),
    Point(DVec2),
}

fn shape() -> impl Strategy<Value = Shape> {
    let p = || (-50.0..50.0f64, -50.0..50.0f64).prop_map(|(x, y)| DVec2::new(x, y));
    prop_oneof![
        (p(), p()).prop_map(|(c, k)| Shape::Rectangle(c, k)),
        (p(), p(), 3u32..8).prop_map(|(c, v, n)| Shape::Polygon(c, v, n)),
        (p(), p(), p()).prop_map(|(a, b, s)| Shape::Slot(a, b, s)),
        (p(), 1.0..30.0f64).prop_map(|(c, r)| Shape::Circle(c, r)),
        (p(), 1.0..30.0f64, -3.0..3.0f64, 0.3..5.5f64)
            .prop_map(|(c, r, a, s)| Shape::Arc(c, r, a, s)),
        p().prop_map(Shape::Point),
    ]
}

fn add_shape(m: &mut Model, f: FeatureId, s: &Shape) {
    let free = |p: DVec2| Snap {
        at: p * MM,
        target: SnapTarget::Free,
    };
    let mut b = Batch::default();
    let made = match s {
        Shape::Rectangle(c, k) => shapes::centre_rectangle(&mut b, &free(*c), &free(*k), false),
        Shape::Polygon(c, v, n) => shapes::polygon(&mut b, &free(*c), &free(*v), *n, false),
        Shape::Slot(a, e, side) => shapes::slot(&mut b, &free(*a), &free(*e), *side * MM, false),
        Shape::Circle(c, r) => {
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
        Shape::Arc(c, r, a, sweep) => {
            let s0 = *c + DVec2::from_angle(*a) * *r;
            let s1 = *c + DVec2::from_angle(a + sweep) * *r;
            b.arc(*c * MM, s0 * MM, s1 * MM, false);
            Ok(())
        }
        Shape::Point(p) => {
            b.point(*p * MM);
            Ok(())
        }
    };
    if made.is_ok() {
        m.edit_sketch(f, b.edits).unwrap();
    }
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(128))]

    /// Random shapes mirrored in random batches about a random line (often
    /// through one of their points): every mirror commits through the
    /// checked path, the sketch stays solved without redundancy, only the
    /// shared points lose a degree of freedom (held on the line), sources
    /// do not move, and every image point is its source's reflection.
    #[test]
    fn random_mirrors_keep_the_sketch_solved(
        shapes_in in proptest::collection::vec(shape(), 1..4),
        through in proptest::option::of(0usize..64),
        angle in -3.2..3.2f64,
        fixed in any::<bool>(),
        batches in proptest::collection::vec(proptest::collection::vec(any::<bool>(), 12), 1..4),
    ) {
        let mut m = Model::new();
        let f = m.add_sketch(PlaneRef::Top);
        for s in &shapes_in {
            add_shape(&mut m, f, s);
        }
        prop_assume!(status(&m, f) == SolveStatus::Solved);
        prop_assume!(redundant(&m, f) == 0);
        // The line: through a point of the shapes, or near the middle.
        let sketch = m.sketch_result(f).unwrap().solution.sketch.clone();
        let pts: Vec<DVec2> = sketch
            .entities
            .iter()
            .filter_map(|e| match e.geometry {
                Geometry::Point { at } => Some(at),
                _ => None,
            })
            .collect();
        let at = match through {
            Some(k) if !pts.is_empty() => pts[k % pts.len()],
            _ => DVec2::new(3.0, -2.0) * MM,
        };
        let dir = DVec2::from_angle(angle);
        let mut b = Batch::default();
        let axis = b.line(at - dir * 0.2, at + dir * 0.2, true);
        if fixed {
            b.constraint(Constraint::Fix { point: axis.start, at: at - dir * 0.2 });
            b.constraint(Constraint::Fix { point: axis.end, at: at + dir * 0.2 });
        }
        m.edit_sketch(f, b.edits).unwrap();
        prop_assume!(status(&m, f) == SolveStatus::Solved);
        let axis = m
            .sketch_def(f)
            .unwrap()
            .sketch
            .entities
            .iter()
            .rev()
            .find(|e| e.construction)
            .unwrap()
            .id;
        for pick in batches {
            let res = m.sketch_result(f).unwrap();
            let shown: Vec<(SketchId, ResolvedCurve)> =
                res.curves.iter().map(|(id, _, c)| (*id, c.clone())).collect();
            let sketch = res.solution.sketch.clone();
            let candidates: Vec<SketchId> = sketch
                .entities
                .iter()
                .filter(|e| e.id != axis)
                .map(|e| e.id)
                .collect();
            let selection: Vec<SketchId> = candidates
                .iter()
                .zip(pick.iter().cycle())
                .filter(|(_, &p)| p)
                .map(|(id, _)| *id)
                .collect();
            let sources = mirror::sources(&sketch, &selection);
            let before = dof(&m, f);
            let mut done = None;
            let mut why = String::new();
            for share in [true, false] {
                let Ok(made) = mirror::mirror(&sketch, &shown, &sources, axis, share) else {
                    break;
                };
                match m.edit_sketch_checked(f, made.plan.edits.clone()) {
                    Ok(_) => {
                        done = Some(made);
                        break;
                    }
                    Err(e) => why = e,
                }
            }
            if done.is_none() && !why.is_empty() {
                return Err(TestCaseError::fail(format!("mirror refused: {why}")));
            }
            let Some(_) = done else { continue };
            prop_assert_eq!(status(&m, f), SolveStatus::Solved);
            prop_assert_eq!(redundant(&m, f), 0);
            let after = m.sketch_result(f).unwrap().solution.sketch.clone();
            // One equation per coincidence class of shared points.
            let newly_held = held_classes(&after, axis) - held_classes(&sketch, axis);
            prop_assert_eq!(dof(&m, f) + newly_held, before);
            // Sources stay put (shared points may settle onto the line).
            for e in &sketch.entities {
                let Geometry::Point { at } = e.geometry else { continue };
                let now = after.point(e.id).unwrap();
                prop_assert!((now - at).length() < 1e-6, "{:?} moved by {}", e.id, (now - at).length());
            }
            // Images are reflections.
            let line = Axis::of(&after, axis).unwrap();
            for (_, pairs) in after.constraints.iter().filter_map(|c| match &c.constraint {
                Constraint::Mirror { axis: a, pairs } if *a == axis => Some((a, pairs)),
                _ => None,
            }) {
                for &(s, i) in pairs {
                    let gs = &after.entity(s).unwrap().geometry;
                    let gi = &after.entity(i).unwrap().geometry;
                    let src = match gs {
                        Geometry::Arc { center, start, end } => vec![*center, *end, *start],
                        Geometry::Point { .. } => vec![s],
                        g => g.referenced_points(),
                    };
                    let img = match gi {
                        Geometry::Point { .. } => vec![i],
                        g => g.referenced_points(),
                    };
                    prop_assert_eq!(src.len(), img.len());
                    for (p, q) in src.iter().zip(&img) {
                        let want = line.reflect(after.point(*p).unwrap());
                        prop_assert!((after.point(*q).unwrap() - want).length() < 1e-9);
                    }
                }
            }
            // Regions still resolve (or the kernel names a tolerance issue).
            if let Err(e) = &m.sketch_result(f).unwrap().regions {
                prop_assert!(e.contains("within the model tolerance"), "{}", e);
            }
        }
    }
}

/// The points the mirrors about `axis` share between a source and its
/// image (held on the line).
fn shared_points(sketch: &opendatum_core_model::SketchData, axis: SketchId) -> Vec<SketchId> {
    let mut out = Vec::new();
    for c in &sketch.constraints {
        let Constraint::Mirror { axis: a, pairs } = &c.constraint else {
            continue;
        };
        if *a != axis {
            continue;
        }
        for &(s, i) in pairs {
            let gs = sketch.entity(s).unwrap().geometry.referenced_points();
            let gi = sketch.entity(i).unwrap().geometry.referenced_points();
            out.extend(gs.into_iter().filter(|p| gi.contains(p)));
        }
    }
    out.sort();
    out.dedup();
    out
}

/// The coincidence classes of the points the mirrors about `axis` share:
/// the mirror holds one point of each on the line.
fn held_classes(sketch: &opendatum_core_model::SketchData, axis: SketchId) -> usize {
    let mut root: std::collections::BTreeMap<SketchId, SketchId> = Default::default();
    fn find(root: &mut std::collections::BTreeMap<SketchId, SketchId>, x: SketchId) -> SketchId {
        let mut r = x;
        while let Some(&p) = root.get(&r) {
            if p == r {
                break;
            }
            r = p;
        }
        r
    }
    for c in &sketch.constraints {
        if let Constraint::Coincident { a, b } = c.constraint {
            let (ra, rb) = (find(&mut root, a), find(&mut root, b));
            if ra != rb {
                root.insert(ra, rb);
            }
        }
    }
    let mut classes: Vec<SketchId> = shared_points(sketch, axis)
        .into_iter()
        .map(|p| find(&mut root, p))
        .collect();
    classes.sort();
    classes.dedup();
    classes.len()
}
