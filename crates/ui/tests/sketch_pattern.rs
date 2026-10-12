//! Sketch patterns (design note 0007 §11), driven through the editor:
//! linear copies along a step the user types, circular copies about a hub
//! they share, copies that follow their sources and keep their places when
//! one is deleted, refusals; and a fuzz of random patterns.

mod common;

use std::f64::consts::{PI, TAU};

use common::*;
use openshape_app::sketch_mode::{annotations, Snap, SnapTarget, Tool};
use openshape_app::sketch_ops::pattern::{self, Anchor, Motion};
use openshape_app::sketch_ops::trim::Preview;
use openshape_app::sketch_ops::{mirror, shapes, Batch};
use openshape_app::Model;
use openshape_core_model::sketch::{Constraint, Geometry};
use openshape_core_model::{DVec2, FeatureId, PlaneRef, SketchData, SketchEdit, SketchId};
use openshape_sketch_solver::{ResolvedCurve, SolveStatus};
use proptest::prelude::*;

fn redundant(m: &Model, f: FeatureId) -> usize {
    m.sketch_result(f).unwrap().solution.redundant.len()
}

fn solved(m: &Model, f: FeatureId) -> SketchData {
    m.sketch_result(f).unwrap().solution.sketch.clone()
}

/// A pattern constraint: its id, its motion's anchor entity (line or
/// centre), the angle of a circular one, and its copies.
type PatternEntry = (
    SketchId,
    SketchId,
    Option<f64>,
    Vec<(SketchId, SketchId, u32)>,
);

fn patterns(m: &Model, f: FeatureId) -> Vec<PatternEntry> {
    m.sketch_def(f)
        .unwrap()
        .sketch
        .constraints
        .iter()
        .filter_map(|c| match &c.constraint {
            Constraint::LinearPattern { direction, copies } => {
                Some((c.id, *direction, None, copies.clone()))
            }
            Constraint::CircularPattern {
                centre,
                angle,
                copies,
            } => Some((c.id, *centre, Some(*angle), copies.clone())),
            _ => None,
        })
        .collect()
}

/// The motion a pattern constraint holds its copies to, in a solved sketch.
fn motion(s: &SketchData, anchor: SketchId, angle: Option<f64>) -> Motion {
    match angle {
        Some(angle) => Motion::Rotate {
            centre: s.point(anchor).unwrap(),
            angle,
        },
        None => {
            let Geometry::Line { start, end } = s.entity(anchor).unwrap().geometry else {
                panic!("the direction is a line")
            };
            Motion::Translate(s.point(end).unwrap() - s.point(start).unwrap())
        }
    }
}

/// The largest distance of a copy point from where its pattern puts it.
fn pattern_error(
    s: &SketchData,
    anchor: SketchId,
    angle: Option<f64>,
    copies: &[(SketchId, SketchId, u32)],
) -> f64 {
    let mv = motion(s, anchor, angle);
    let points = |id: SketchId| match &s.entity(id).unwrap().geometry {
        Geometry::Point { .. } => vec![id],
        g => g.referenced_points(),
    };
    let mut worst: f64 = 0.0;
    for &(src, img, k) in copies {
        let (p, q) = (points(src), points(img));
        assert_eq!(p.len(), q.len());
        for (p, q) in p.into_iter().zip(q) {
            let want = mv.point(s.point(p).unwrap(), k);
            worst = worst.max((s.point(q).unwrap() - want).length());
        }
    }
    worst
}

/// A 10 × 10 mm square with its lower left corner at the origin (fixed
/// there), drawn with the rectangle tool; returns its four lines.
fn square(
    m: &mut Model,
    ed: &mut openshape_app::sketch_mode::SketchEditor,
    map: &openshape_app::sketch_mode::ScreenMap,
) -> Vec<SketchId> {
    draw(m, ed, map, Tool::Rectangle, &[(0.0, 0.0), (10.0, 10.0)]);
    m.sketch_def(ed.feature)
        .unwrap()
        .sketch
        .entities
        .iter()
        .filter(|e| matches!(e.geometry, Geometry::Line { .. }))
        .map(|e| e.id)
        .collect()
}

#[test]
fn a_linear_pattern_copies_the_selection_and_its_step_is_typed() {
    let (mut m, mut ed, map) = setup();
    let f = ed.feature;
    let lines = square(&mut m, &mut ed, &map);
    assert_eq!(lines.len(), 4);
    assert_eq!(dof(&m, f), 2);
    ed.selection = lines.clone();
    ed.set_tool(Tool::LinearPattern);
    assert_eq!(
        ed.hint(),
        "Click where the step between copies of the 4 selected items starts"
    );
    assert_eq!(ed.pattern_count, 3);
    // The step starts at the corner; hovering its end previews the copies.
    click(&mut m, &mut ed, &map, 0.0, 0.0);
    assert_eq!(ed.clicks.len(), 1);
    assert_eq!(
        ed.hint(),
        "Click where the step to the first copy ends; Esc starts it again"
    );
    hover(&m, &mut ed, &map, 20.0, 0.0);
    let plan = ed.curve_edit.clone().expect("a plan").plan.expect("ok");
    let Preview::Shapes { curves: shown, .. } = plan.preview else {
        panic!("{:?}", plan.preview)
    };
    assert_eq!(shown.len(), 8);
    assert!(shown
        .iter()
        .any(|c| matches!(c, ResolvedCurve::Line { start, end }
        if (*start - DVec2::new(40.0, 0.0) * MM).length() < 1e-6
            && (*end - DVec2::new(50.0, 0.0) * MM).length() < 1e-6)));
    assert_eq!(ed.rubber_band().len(), 1);
    click(&mut m, &mut ed, &map, 20.0, 0.0);
    assert_eq!(
        ed.message.as_deref(),
        Some("Patterned 4 items into 3 instances")
    );
    // Back to Select, with the step's length open for typing.
    assert_eq!(ed.tool, Tool::Select);
    assert!(ed.selection.is_empty());
    let prompt = ed.prompt.clone().expect("the step's length");
    assert!(matches!(prompt.constraint, Constraint::Length { .. }));
    assert_eq!(prompt.text, "20");
    // One constraint, eight copies; the direction line's free end is the
    // only new freedom.
    let pats = patterns(&m, f);
    assert_eq!(pats.len(), 1);
    assert_eq!(pats[0].3.len(), 8);
    assert_eq!(construction(&m, f), 1);
    assert_eq!(status(&m, f), SolveStatus::Solved);
    assert_eq!(redundant(&m, f), 0);
    assert_eq!(dof(&m, f), 4);
    let areas = region_areas(&m, f);
    assert_eq!(areas.len(), 3);
    assert!(
        areas.iter().all(|a| close(*a, 100.0 * MM * MM)),
        "{areas:?}"
    );
    // Typing 15 mm moves the copies only.
    ed.prompt.as_mut().unwrap().text = "15".into();
    ed.commit_prompt(&mut m);
    assert!(ed.prompt.is_none());
    assert_eq!(status(&m, f), SolveStatus::Solved);
    assert_eq!(redundant(&m, f), 0);
    assert_eq!(dof(&m, f), 3);
    let s = solved(&m, f);
    let (_, line, _, copies) = pats[0].clone();
    assert!(pattern_error(&s, line, None, &copies) < 1e-12);
    let Geometry::Line { start, end } = s.entity(line).unwrap().geometry else {
        panic!()
    };
    assert!(((s.point(end).unwrap() - s.point(start).unwrap()).length() - 0.015).abs() < 1e-12);
    for l in &lines {
        let Geometry::Line { start, .. } = s.entity(*l).unwrap().geometry else {
            panic!()
        };
        let p = s.point(start).unwrap();
        assert!(p.x.abs() < 1e-9 || (p.x - 0.01).abs() < 1e-9, "{p:?}");
    }
    let areas = region_areas(&m, f);
    assert_eq!(areas.len(), 3);
    let x = extrude(&mut m, &map, &[(35.0, 5.0)], "10 mm");
    assert!((volume(&m, x) - 1e-6).abs() < 1e-14);
}

#[test]
fn a_circular_pattern_turns_spokes_about_their_hub() {
    let (mut m, mut ed, map) = setup();
    let f = ed.feature;
    // A spoke from a hub point held at the origin.
    draw(
        &mut m,
        &mut ed,
        &map,
        Tool::Line,
        &[(0.0, 0.0), (20.0, 0.0)],
    );
    let spoke = m
        .sketch_def(f)
        .unwrap()
        .sketch
        .entities
        .iter()
        .find(|e| matches!(e.geometry, Geometry::Line { .. }))
        .unwrap()
        .id;
    let Geometry::Line { start: hub, .. } = m
        .sketch_def(f)
        .unwrap()
        .sketch
        .entity(spoke)
        .unwrap()
        .geometry
    else {
        panic!()
    };
    // The line was drawn level, so it is horizontal: its length is free.
    assert_eq!(dof(&m, f), 1);
    ed.selection = vec![spoke];
    ed.pattern_count = 4;
    ed.set_tool(Tool::CircularPattern);
    assert_eq!(
        ed.hint(),
        "Click the centre to turn the selected item about"
    );
    hover(&m, &mut ed, &map, 0.0, 0.0);
    assert_eq!(ed.cursor.unwrap().target, SnapTarget::Point(hub));
    let plan = ed.curve_edit.clone().unwrap().plan.unwrap();
    let Preview::Shapes { curves: shown, .. } = plan.preview else {
        panic!()
    };
    assert_eq!(shown.len(), 3);
    click(&mut m, &mut ed, &map, 0.0, 0.0);
    assert_eq!(
        ed.message.as_deref(),
        Some("Patterned 1 item into 4 instances")
    );
    // The copies share the hub: no new points there, no new freedom.
    let pats = patterns(&m, f);
    assert_eq!(pats.len(), 1);
    let (id, centre, angle, copies) = pats[0].clone();
    assert_eq!(centre, hub);
    assert!((angle.unwrap() - TAU / 4.0).abs() < 1e-15);
    assert_eq!(copies.iter().map(|c| c.2).collect::<Vec<_>>(), [1, 2, 3]);
    let def = m.sketch_def(f).unwrap().sketch.clone();
    for &(_, img, _) in &copies {
        let Geometry::Line { start, .. } = def.entity(img).unwrap().geometry else {
            panic!()
        };
        assert_eq!(start, hub);
    }
    assert_eq!(status(&m, f), SolveStatus::Solved);
    assert_eq!(redundant(&m, f), 0);
    assert_eq!(dof(&m, f), 1);
    // The angle opens for typing; 60° moves the copies only.
    let prompt = ed.prompt.clone().expect("the angle");
    assert!(prompt.is_angle());
    assert_eq!(prompt.existing, Some(id));
    assert_eq!(prompt.text, "90");
    ed.prompt.as_mut().unwrap().text = "60".into();
    ed.commit_prompt(&mut m);
    assert!(ed.prompt.is_none(), "{:?}", ed.prompt);
    let s = solved(&m, f);
    assert!(pattern_error(&s, hub, Some(PI / 3.0), &copies) < 1e-12);
    let Geometry::Line { end, .. } = s.entity(spoke).unwrap().geometry else {
        panic!()
    };
    assert!((s.point(end).unwrap() - DVec2::new(0.02, 0.0)).length() < 1e-9);
    let ends: Vec<DVec2> = copies
        .iter()
        .map(|&(_, img, _)| {
            let Geometry::Line { end, .. } = s.entity(img).unwrap().geometry else {
                panic!()
            };
            s.point(end).unwrap()
        })
        .collect();
    let want = [60.0f64, 120.0, 180.0].map(|a| DVec2::from_angle(a.to_radians()) * 0.02);
    for (got, want) in ends.iter().zip(want) {
        assert!((*got - want).length() < 1e-9, "{got:?} {want:?}");
    }
    // The angle is labelled.
    let d = ed.display(&m).unwrap();
    let labels = annotations(&d, &map);
    assert!(labels.iter().any(|a| a.constraint == id && a.text == "60°"));
}

#[test]
fn a_circular_pattern_of_a_circle_about_the_origin() {
    let (mut m, mut ed, map) = setup();
    let f = ed.feature;
    draw(
        &mut m,
        &mut ed,
        &map,
        Tool::Circle,
        &[(30.0, 0.0), (35.0, 0.0)],
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
    assert_eq!(dof(&m, f), 3);
    ed.selection = vec![circle];
    ed.pattern_count = 6;
    ed.set_tool(Tool::CircularPattern);
    click(&mut m, &mut ed, &map, 0.0, 0.0);
    assert_eq!(
        ed.message.as_deref(),
        Some("Patterned 1 item into 6 instances")
    );
    // A new centre point, fixed at the origin.
    assert_eq!(status(&m, f), SolveStatus::Solved);
    assert_eq!(redundant(&m, f), 0);
    assert_eq!(dof(&m, f), 3);
    let areas = region_areas(&m, f);
    assert_eq!(areas.len(), 6);
    let r = 5.0 * MM;
    assert!(
        areas
            .iter()
            .all(|a| (a - PI * r * r).abs() < 1e-3 * PI * r * r),
        "{areas:?}"
    );
    // The copies keep the source's radius when it changes.
    let mut b = Batch::default();
    b.constraint(Constraint::Radius {
        curve: circle,
        value: 0.004,
    });
    m.edit_sketch(f, b.edits).unwrap();
    assert_eq!(status(&m, f), SolveStatus::Solved);
    let s = solved(&m, f);
    let (_, centre, angle, copies) = patterns(&m, f).remove(0);
    assert!(pattern_error(&s, centre, angle, &copies) < 1e-12);
    for &(_, img, _) in &copies {
        let Geometry::Circle { radius, .. } = s.entity(img).unwrap().geometry else {
            panic!()
        };
        assert!((radius - 0.004).abs() < 1e-12);
    }
}

#[test]
fn deleting_a_copy_keeps_the_rest_and_deleting_the_line_frees_them() {
    let (mut m, mut ed, map) = setup();
    let f = ed.feature;
    let lines = square(&mut m, &mut ed, &map);
    ed.selection = lines;
    ed.set_tool(Tool::LinearPattern);
    click(&mut m, &mut ed, &map, 0.0, 0.0);
    click(&mut m, &mut ed, &map, 20.0, 0.0);
    ed.escape();
    let (_, line, _, copies) = patterns(&m, f).remove(0);
    let before = dof(&m, f);
    ed.selection = vec![copies[0].1];
    ed.delete_selection(&mut m);
    let (_, _, _, kept) = patterns(&m, f).remove(0);
    assert_eq!(kept.len(), copies.len() - 1);
    assert_eq!(status(&m, f), SolveStatus::Solved);
    assert_eq!(redundant(&m, f), 0);
    assert_eq!(dof(&m, f), before);
    ed.selection = vec![line];
    ed.delete_selection(&mut m);
    assert!(patterns(&m, f).is_empty());
    assert_eq!(status(&m, f), SolveStatus::Solved);
    // Seven free lines, their points no longer held (two each, four
    // coordinates, less the direction line's free end).
    assert_eq!(dof(&m, f), before - 2 + 7 * 4);
}

#[test]
fn patterns_that_cannot_be_made_are_refused_or_skipped() {
    let (mut m, mut ed, map) = setup();
    let f = ed.feature;
    // Nothing selected.
    ed.set_tool(Tool::LinearPattern);
    assert_eq!(
        ed.hint(),
        "Select the curves and points to pattern, then choose the tool"
    );
    click(&mut m, &mut ed, &map, 0.0, 0.0);
    assert_eq!(
        ed.message.as_deref(),
        Some("Select the curves and points to pattern first")
    );
    assert!(ed.clicks.is_empty());
    // A point at the centre is its own copy.
    draw(&mut m, &mut ed, &map, Tool::Point, &[(5.0, 5.0)]);
    let point = m.sketch_def(f).unwrap().sketch.entities[0].id;
    ed.selection = vec![point];
    ed.set_tool(Tool::CircularPattern);
    hover(&m, &mut ed, &map, 5.0, 5.0);
    assert_eq!(
        ed.curve_edit.clone().unwrap().plan.unwrap_err(),
        "The selection lies on the centre: turning it leaves it where it is"
    );
    // Counts outside 2 to 100.
    let s = m.sketch_def(f).unwrap().sketch.clone();
    let free = |x: f64, y: f64| Snap {
        at: DVec2::new(x, y) * MM,
        target: SnapTarget::Free,
    };
    for n in [0, 1, 101] {
        let err = pattern::pattern(
            &s,
            &[],
            &[point],
            Anchor::Circular {
                centre: free(0.0, 0.0),
            },
            n,
        )
        .unwrap_err();
        assert!(err.starts_with("A pattern has 2 to 100 instances"), "{err}");
    }
    let err = pattern::pattern(
        &s,
        &[],
        &[point],
        Anchor::Linear {
            from: free(1.0, 1.0),
            to: free(1.0, 1.0),
        },
        3,
    )
    .unwrap_err();
    assert_eq!(err, "The step between copies has no length");
    // An angle that turns a copy onto its source, or a turn or more.
    draw(
        &mut m,
        &mut ed,
        &map,
        Tool::Line,
        &[(20.0, 0.0), (30.0, 0.0)],
    );
    let line = m
        .sketch_def(f)
        .unwrap()
        .sketch
        .entities
        .iter()
        .find(|e| matches!(e.geometry, Geometry::Line { .. }))
        .unwrap()
        .id;
    ed.selection = vec![line];
    ed.pattern_count = 3;
    ed.set_tool(Tool::CircularPattern);
    click(&mut m, &mut ed, &map, -20.0, 0.0);
    for (typed, why) in [
        ("180", "copy 2 would be turned onto its source (360°)"),
        ("360", "a pattern angle is more than 0° and less than 360°"),
        ("0", "a pattern angle is more than 0° and less than 360°"),
    ] {
        let p = ed.prompt.as_mut().expect("the angle prompt");
        p.text = typed.into();
        p.error = None;
        ed.commit_prompt(&mut m);
        let p = ed.prompt.as_ref().expect("still open");
        assert_eq!(p.error.as_deref(), Some(why));
    }
    assert_eq!(status(&m, f), SolveStatus::Solved);
    assert_eq!(redundant(&m, f), 0);
}

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

    /// Random shapes patterned in random batches, linearly along a random
    /// step or circularly about a random centre (often one of their
    /// points): every pattern commits through the checked path, the sketch
    /// stays solved without redundancy, the only new freedom is the
    /// anchor's (a linear pattern's line, a new free centre), sources do
    /// not move, and every copy point is where its pattern puts it, also
    /// after the step or angle changes.
    #[test]
    fn random_patterns_keep_the_sketch_solved(
        shapes_in in proptest::collection::vec(shape(), 1..4),
        batches in proptest::collection::vec(
            (
                proptest::collection::vec(any::<bool>(), 12),
                any::<bool>(),
                proptest::option::of(0usize..64),
                (-60.0..60.0f64, -60.0..60.0f64),
                (-60.0..60.0f64, -60.0..60.0f64),
                2u32..5,
                0.2..1.8f64,
            ),
            1..4,
        ),
    ) {
        let mut m = Model::new();
        let f = m.add_sketch(PlaneRef::Top);
        for s in &shapes_in {
            add_shape(&mut m, f, s);
        }
        prop_assume!(status(&m, f) == SolveStatus::Solved);
        prop_assume!(redundant(&m, f) == 0);
        // Sources come from the shapes (copies of copies are the solver
        // tests' business; here they would only grow the sketch).
        let originals: Vec<SketchId> = m
            .sketch_def(f)
            .unwrap()
            .sketch
            .entities
            .iter()
            .map(|e| e.id)
            .collect();
        for (pick, linear, snap_to, a, b, count, scale) in batches {
            let res = m.sketch_result(f).unwrap();
            let shown: Vec<(SketchId, ResolvedCurve)> =
                res.curves.iter().map(|(id, _, c)| (*id, c.clone())).collect();
            let sketch = res.solution.sketch.clone();
            let selection: Vec<SketchId> = originals
                .iter()
                .zip(pick.iter().cycle())
                .filter(|(_, &p)| p)
                .map(|(id, _)| *id)
                .collect();
            let sources = mirror::sources(&sketch, &selection);
            let pts: Vec<(SketchId, DVec2)> = sketch
                .entities
                .iter()
                .filter_map(|e| match e.geometry {
                    Geometry::Point { at } => Some((e.id, at)),
                    _ => None,
                })
                .collect();
            let snapped = |k: Option<usize>, at: (f64, f64)| match k {
                Some(k) if !pts.is_empty() => {
                    let (id, p) = pts[k % pts.len()];
                    Snap { at: p, target: SnapTarget::Point(id) }
                }
                _ => Snap { at: DVec2::new(at.0, at.1) * MM, target: SnapTarget::Free },
            };
            let anchor = if linear {
                Anchor::Linear { from: snapped(snap_to, a), to: snapped(None, b) }
            } else {
                Anchor::Circular { centre: snapped(snap_to, a) }
            };
            let Ok(made) = pattern::pattern(&sketch, &shown, &sources, anchor, count) else {
                continue;
            };
            let before = dof(&m, f);
            if let Err(e) = m.edit_sketch_checked(f, made.plan.edits.clone()) {
                return Err(TestCaseError::fail(format!("pattern refused: {e}")));
            }
            prop_assert_eq!(status(&m, f), SolveStatus::Solved);
            prop_assert_eq!(redundant(&m, f), 0);
            // The anchor's freedom: a line's free ends, a new free centre.
            let anchor_dof = match anchor {
                Anchor::Linear { from, to } => [from, to]
                    .iter()
                    .filter(|s| s.target == SnapTarget::Free)
                    .count() * 2,
                Anchor::Circular { centre } => {
                    if centre.target == SnapTarget::Free { 2 } else { 0 }
                }
            };
            prop_assert_eq!(dof(&m, f), before + anchor_dof);
            let after = solved(&m, f);
            for e in &sketch.entities {
                let Geometry::Point { at } = e.geometry else { continue };
                let now = after.point(e.id).unwrap();
                prop_assert!((now - at).length() < 1e-9, "{:?} moved by {}", e.id, (now - at).length());
            }
            let (id, anchor_id, angle, copies) = patterns(&m, f).pop().unwrap();
            prop_assert!(pattern_error(&after, anchor_id, angle, &copies) < 1e-9);
            // A new angle or step moves the copies, not the sources.
            let edit = match angle {
                Some(angle) => SketchEdit::SetValue { constraint: id, value: angle * scale },
                None => {
                    let mut b = Batch::default();
                    let Geometry::Line { start, end } = after.entity(anchor_id).unwrap().geometry else {
                        unreachable!()
                    };
                    let len = (after.point(end).unwrap() - after.point(start).unwrap()).length();
                    b.constraint(Constraint::Length { line: anchor_id, value: len * scale });
                    b.edits.remove(0)
                }
            };
            if angle.is_some_and(|a| {
                (1..count).any(|k| {
                    let t = a * scale * f64::from(k) / TAU;
                    (t - t.round()).abs() < 1e-6
                })
            }) {
                continue;
            }
            m.edit_sketch(f, vec![edit]).unwrap();
            prop_assert_eq!(status(&m, f), SolveStatus::Solved);
            let moved = solved(&m, f);
            prop_assert!(pattern_error(&moved, anchor_id, angle.map(|a| a * scale), &copies) < 1e-9);
            if angle.is_some() {
                for e in &sketch.entities {
                    let Geometry::Point { at } = e.geometry else { continue };
                    let now = moved.point(e.id).unwrap();
                    prop_assert!((now - at).length() < 1e-9, "{:?} moved by {}", e.id, (now - at).length());
                }
            }
        }
    }
}
