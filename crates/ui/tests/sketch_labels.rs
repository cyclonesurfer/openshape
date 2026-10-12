//! Placing dimension labels (design note 0007 §12), driven through the
//! editor: a label dragged to a spot stays there in one undo step, its
//! dimension line follows it, it moves with what it measures, and every
//! dimension kind accepts a placed label.

mod common;

use common::*;
use egui::Pos2;
use opendatum_app::sketch_mode::{annotations, Annotation, Pick, ScreenMap, SketchEditor, Tool};
use opendatum_app::Model;
use opendatum_core_model::sketch::{Constraint, DistanceKind};
use opendatum_core_model::{DVec2, SketchId};

/// The label of dimension `id` as drawn.
fn label(m: &Model, ed: &SketchEditor, map: &ScreenMap, id: SketchId) -> Annotation {
    let d = ed.display(m).unwrap();
    annotations(&d, map)
        .into_iter()
        .find(|a| a.constraint == id && a.dimension)
        .expect("the dimension has a label")
}

/// Drags dimension `id`'s label by its centre to `to` (millimetres).
fn drag_label(m: &mut Model, ed: &mut SketchEditor, map: &ScreenMap, id: SketchId, to: (f64, f64)) {
    let from = label(m, ed, map, id).rect.center();
    let to = at(map, to.0, to.1);
    ed.drag_start(m, map, from);
    for i in 1..=4 {
        ed.drag_to(m, map, from + (to - from) * (i as f32 / 4.0));
    }
    ed.drag_end(m, map, to);
}

/// Opens the dimension prompt for the selection, optionally switches a
/// distance's kind, and commits the value shown; returns the new dimension.
fn dimension(
    m: &mut Model,
    ed: &mut SketchEditor,
    map: &ScreenMap,
    kind: Option<DistanceKind>,
) -> SketchId {
    let before: Vec<SketchId> = constraint_ids(m, ed);
    ed.dimension(m, at(map, 90.0, 90.0));
    if let Some(kind) = kind {
        ed.set_prompt_kind(m, kind);
    }
    ed.commit_prompt(m);
    assert!(ed.prompt.is_none(), "{:?}", ed.prompt);
    let after = constraint_ids(m, ed);
    let new: Vec<SketchId> = after.into_iter().filter(|c| !before.contains(c)).collect();
    let def = m.sketch_def(ed.feature).unwrap();
    *new.iter()
        .find(|c| {
            def.sketch
                .constraint(**c)
                .unwrap()
                .constraint
                .value()
                .is_some()
        })
        .expect("a new dimension")
}

fn constraint_ids(m: &Model, ed: &SketchEditor) -> Vec<SketchId> {
    m.sketch_def(ed.feature)
        .unwrap()
        .sketch
        .constraints
        .iter()
        .map(|c| c.id)
        .collect()
}

fn near(a: Pos2, b: Pos2) -> bool {
    a.distance(b) < 0.01
}

/// A horizontal line from (-20, 0) to (20, 0) with its length dimension.
fn dimensioned_line() -> (Model, SketchEditor, ScreenMap, SketchId) {
    let (mut m, mut ed, map) = setup();
    draw(
        &mut m,
        &mut ed,
        &map,
        Tool::Line,
        &[(-20.0, 0.0), (20.0, 0.0)],
    );
    ed.escape();
    ed.set_tool(Tool::Select);
    click(&mut m, &mut ed, &map, 5.0, 0.0);
    let id = dimension(&mut m, &mut ed, &map, None);
    ed.selection.clear();
    (m, ed, map, id)
}

#[test]
fn a_dragged_length_label_stays_where_it_is_put_in_one_undo_step() {
    let (mut m, mut ed, map, id) = dimensioned_line();
    let f = ed.feature;
    let geometry = m.sketch_def(f).unwrap().sketch.clone();
    assert_eq!(m.sketch_def(f).unwrap().label(id), None);
    let ops = m.session().log().len();
    drag_label(&mut m, &mut ed, &map, id, (0.0, 30.0));
    assert_eq!(m.session().log().len(), ops + 1, "one op");
    let def = m.sketch_def(f).unwrap();
    assert_eq!(def.sketch, geometry, "moving a label moves no geometry");
    let offset = def.label(id).expect("placed");
    assert!(
        (offset - DVec2::new(0.0, 30.0) * MM).length() < 1e-6,
        "{offset}"
    );
    // Drawn there, with the dimension line through it, parallel to the
    // line and as long, and extension lines up from the line's ends.
    let a = label(&m, &ed, &map, id);
    assert!(near(a.rect.center(), at(&map, 0.0, 30.0)), "{:?}", a.rect);
    let y = at(&map, 0.0, 30.0).y;
    assert!(a.lines.iter().any(|[p, q]| {
        (p.y - y).abs() < 0.01
            && (q.y - y).abs() < 0.01
            && near(*p, at(&map, -20.0, 30.0))
            && near(*q, at(&map, 20.0, 30.0))
    }));
    for end in [-20.0, 20.0] {
        assert!(a
            .lines
            .iter()
            .any(|[p, q]| near(*p, at(&map, end, 0.0)) && q.y < y - 1.0));
    }
    // It is picked and edited where it now is.
    ed.hover(&m, &map, Some(at(&map, 0.0, 30.0)));
    assert_eq!(ed.hover, Some(Pick::Constraint(id)));
    ed.double_click(&mut m, &map, at(&map, 0.0, 30.0));
    assert_eq!(ed.prompt.as_ref().map(|p| p.text.as_str()), Some("40"));
    ed.prompt = None;
    assert!(m.undo());
    assert_eq!(m.sketch_def(f).unwrap().label(id), None);
    assert!(m.redo());
    assert_eq!(m.sketch_def(f).unwrap().label(id), Some(offset));
}

#[test]
fn past_an_end_the_dimension_line_reaches_the_label() {
    let (mut m, mut ed, map, id) = dimensioned_line();
    drag_label(&mut m, &mut ed, &map, id, (45.0, -15.0));
    let a = label(&m, &ed, &map, id);
    assert!(near(a.rect.center(), at(&map, 45.0, -15.0)));
    assert!(a
        .lines
        .iter()
        .any(|[p, q]| near(*p, at(&map, -20.0, -15.0)) && near(*q, at(&map, 45.0, -15.0))));
    // The extension lines run down from the line, past the dimension line.
    let y = at(&map, 0.0, -15.0).y;
    assert!(a
        .lines
        .iter()
        .any(|[p, q]| near(*p, at(&map, 20.0, 0.0)) && q.y > y + 1.0));
}

#[test]
fn a_placed_label_moves_with_what_it_measures() {
    let (mut m, mut ed, map, id) = dimensioned_line();
    drag_label(&mut m, &mut ed, &map, id, (0.0, 30.0));
    // Drag the line by its middle: it moves, the label with it.
    ed.drag_start(&mut m, &map, at(&map, 5.0, 0.0));
    ed.drag_to(&m, &map, at(&map, 15.0, 5.0));
    ed.drag_end(&mut m, &map, at(&map, 15.0, 5.0));
    let a = label(&m, &ed, &map, id);
    assert!(
        a.rect.center().distance(at(&map, 10.0, 35.0)) < 0.05,
        "{:?}",
        a.rect.center()
    );
    // While a label is dragged the drawing follows the pointer; releasing
    // commits where it ends.
    let from = a.rect.center();
    ed.drag_start(&mut m, &map, from);
    ed.drag_to(&m, &map, at(&map, -30.0, 50.0));
    let mid = label(&m, &ed, &map, id);
    assert!(near(mid.rect.center(), at(&map, -30.0, 50.0)));
    ed.drag_to(&m, &map, at(&map, 0.0, -40.0));
    ed.drag_end(&mut m, &map, at(&map, 0.0, -40.0));
    let a = label(&m, &ed, &map, id);
    assert!(near(a.rect.center(), at(&map, 0.0, -40.0)));
}

#[test]
fn a_drag_that_ends_where_it_started_places_nothing() {
    let (mut m, mut ed, map, id) = dimensioned_line();
    let from = label(&m, &ed, &map, id).rect.center();
    let ops = m.session().log().len();
    ed.drag_start(&mut m, &map, from);
    ed.drag_to(&m, &map, from + egui::vec2(30.0, 0.0));
    ed.drag_to(&m, &map, from);
    ed.drag_end(&mut m, &map, from);
    assert_eq!(m.session().log().len(), ops);
    assert_eq!(m.sketch_def(ed.feature).unwrap().label(id), None);
}

#[test]
fn deleting_a_dimension_drops_its_placed_label() {
    let (mut m, mut ed, map, id) = dimensioned_line();
    drag_label(&mut m, &mut ed, &map, id, (0.0, 30.0));
    ed.selection = vec![id];
    ed.delete_selection(&mut m);
    let def = m.sketch_def(ed.feature).unwrap();
    assert!(def.sketch.constraint(id).is_none());
    assert!(def.labels.is_empty());
}

#[test]
fn every_dimension_kind_takes_a_placed_label() {
    let (mut m, mut ed, map) = setup();
    let f = ed.feature;
    let mut dims = Vec::new();

    // A circle's diameter.
    draw(
        &mut m,
        &mut ed,
        &map,
        Tool::Circle,
        &[(40.0, 40.0), (55.0, 40.0)],
    );
    click(&mut m, &mut ed, &map, 40.0, 55.0);
    dims.push((dimension(&mut m, &mut ed, &map, None), (70.0, 70.0)));
    // An arc's radius (centre (-50, 0)).
    draw(
        &mut m,
        &mut ed,
        &map,
        Tool::Arc,
        &[(-40.0, 0.0), (-60.0, 0.0), (-50.0, 10.0)],
    );
    ed.selection.clear();
    click(&mut m, &mut ed, &map, -50.0, 10.0);
    dims.push((dimension(&mut m, &mut ed, &map, None), (-50.0, -30.0)));
    // The angle between two lines from (0, -60), placed across the vertex.
    ed.set_tool(Tool::Line);
    click(&mut m, &mut ed, &map, 0.0, -60.0);
    click(&mut m, &mut ed, &map, 40.0, -60.0);
    ed.escape();
    click(&mut m, &mut ed, &map, 0.0, -60.0);
    click(&mut m, &mut ed, &map, 30.0, -30.0);
    ed.escape();
    ed.set_tool(Tool::Select);
    ed.selection.clear();
    click(&mut m, &mut ed, &map, 20.0, -60.0);
    ed.click(&mut m, &map, at(&map, 15.0, -45.0), true);
    dims.push((dimension(&mut m, &mut ed, &map, None), (-10.0, -65.0)));
    // Horizontal and vertical distances between two points.
    draw(
        &mut m,
        &mut ed,
        &map,
        Tool::Point,
        &[(-80.0, 60.0), (-50.0, 80.0)],
    );
    for (kind, to) in [
        (DistanceKind::Horizontal, (-60.0, 95.0)),
        (DistanceKind::Vertical, (-100.0, 75.0)),
    ] {
        ed.selection.clear();
        click(&mut m, &mut ed, &map, -80.0, 60.0);
        ed.click(&mut m, &map, at(&map, -50.0, 80.0), true);
        dims.push((dimension(&mut m, &mut ed, &map, Some(kind)), to));
    }
    // A point's distance from a line.
    draw(&mut m, &mut ed, &map, Tool::Point, &[(100.0, -20.0)]);
    ed.set_tool(Tool::Line);
    click(&mut m, &mut ed, &map, 80.0, -60.0);
    click(&mut m, &mut ed, &map, 120.0, -60.0);
    ed.escape();
    ed.set_tool(Tool::Select);
    ed.selection.clear();
    click(&mut m, &mut ed, &map, 100.0, -20.0);
    ed.click(&mut m, &map, at(&map, 110.0, -60.0), true);
    dims.push((dimension(&mut m, &mut ed, &map, None), (115.0, -40.0)));
    ed.selection.clear();

    let def = m.sketch_def(f).unwrap();
    let kinds: Vec<&str> = dims
        .iter()
        .map(
            |(id, _)| match def.sketch.constraint(*id).unwrap().constraint {
                Constraint::Diameter { .. } => "diameter",
                Constraint::Radius { .. } => "radius",
                Constraint::Angle { .. } => "angle",
                Constraint::Distance {
                    kind: DistanceKind::Horizontal,
                    ..
                } => "horizontal",
                Constraint::Distance {
                    kind: DistanceKind::Vertical,
                    ..
                } => "vertical",
                Constraint::Distance { .. } => "point to line",
                _ => "other",
            },
        )
        .collect();
    assert_eq!(
        kinds,
        [
            "diameter",
            "radius",
            "angle",
            "horizontal",
            "vertical",
            "point to line"
        ]
    );
    let geometry = m.sketch_result(f).unwrap().solution.sketch.clone();
    for (id, to) in &dims {
        let ops = m.session().log().len();
        drag_label(&mut m, &mut ed, &map, *id, *to);
        assert_eq!(m.session().log().len(), ops + 1, "{id:?}");
        let a = label(&m, &ed, &map, *id);
        assert!(
            near(a.rect.center(), at(&map, to.0, to.1)),
            "{} at {:?}",
            a.text,
            a.rect.center()
        );
        assert!(a
            .lines
            .iter()
            .flatten()
            .all(|p| p.x.is_finite() && p.y.is_finite()));
        ed.hover(&m, &map, Some(at(&map, to.0, to.1)));
        assert_eq!(ed.hover, Some(Pick::Constraint(*id)), "{}", a.text);
    }
    assert_eq!(
        m.sketch_result(f).unwrap().solution.sketch,
        geometry,
        "no label moved any geometry"
    );
    assert_eq!(m.sketch_def(f).unwrap().labels.len(), dims.len());
    // The radius points from the arc's centre towards its label.
    let radius = label(&m, &ed, &map, dims[1].0);
    let c = at(&map, -50.0, 0.0);
    assert!(radius.lines.iter().any(|[p, q]| near(*p, c)
        && ((*q - *p).normalized() - (at(&map, -50.0, -30.0) - c).normalized()).length() < 1e-4));
    // The angle's arc lies across the vertex, on the label's side.
    let angle = label(&m, &ed, &map, dims[2].0);
    let vertex = at(&map, 0.0, -60.0);
    let r = at(&map, -10.0, -65.0).distance(vertex);
    let on_arc: Vec<Pos2> = angle
        .lines
        .iter()
        .flatten()
        .copied()
        .filter(|p| (p.distance(vertex) - r).abs() < 0.05)
        .collect();
    assert!(!on_arc.is_empty());
    assert!(on_arc.iter().all(|p| p.x < vertex.x), "{on_arc:?}");
}

#[test]
fn labels_return_to_automatic_placement_selected_or_all() {
    let (mut m, mut ed, map, length) = dimensioned_line();
    // A second dimension: the distance from the origin's side to a point.
    draw(&mut m, &mut ed, &map, Tool::Point, &[(0.0, 40.0)]);
    click(&mut m, &mut ed, &map, 0.0, 40.0);
    ed.click(&mut m, &map, at(&map, 5.0, 0.0), true);
    let distance = dimension(&mut m, &mut ed, &map, None);
    ed.selection.clear();
    let auto = label(&m, &ed, &map, length).rect.center();
    assert!(ed.labels_to_reset(&m).is_empty(), "nothing placed yet");
    drag_label(&mut m, &mut ed, &map, length, (0.0, -30.0));
    drag_label(&mut m, &mut ed, &map, distance, (30.0, 20.0));
    let f = ed.feature;
    // With one selected, only it.
    ed.selection = vec![distance];
    assert_eq!(ed.labels_to_reset(&m), vec![distance]);
    ed.reset_labels(&mut m);
    let def = m.sketch_def(f).unwrap();
    assert_eq!(def.label(distance), None);
    assert!(def.label(length).is_some());
    // With no dimension selected (a curve does not count), every one; in
    // one undo step, and drawn where it was before it was placed.
    ed.selection.clear();
    click(&mut m, &mut ed, &map, 5.0, 0.0);
    assert_eq!(ed.labels_to_reset(&m), vec![length]);
    drag_label(&mut m, &mut ed, &map, distance, (-30.0, 20.0));
    ed.selection.clear();
    let ops = m.session().log().len();
    ed.reset_labels(&mut m);
    assert_eq!(m.session().log().len(), ops + 1);
    assert!(m.sketch_def(f).unwrap().labels.is_empty());
    assert!(near(label(&m, &ed, &map, length).rect.center(), auto));
    assert!(m.undo());
    assert_eq!(m.sketch_def(f).unwrap().labels.len(), 2);
}
