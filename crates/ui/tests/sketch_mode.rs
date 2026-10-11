//! Sketch mode end to end, without a window: the editor is driven with
//! screen-space clicks and drags through a top-view camera, exactly as the
//! app's viewport drives it, and the results are checked in the document,
//! the solver's diagnostics, the regions and the bodies built from them.

use egui::{Pos2, Rect};
use openshape_app::profile_pick::{ProfilePick, SweepKind};
use openshape_app::sketch_mode::{normal_view, ConstraintTool, ScreenMap, SketchEditor, Tool};
use openshape_app::Model;
use openshape_core_model::sketch::{Constraint, Geometry};
use openshape_core_model::{DVec2, DVec3, FeatureId, PlaneRef};
use openshape_features::FeatureStatus;
use openshape_kernel::{Frame, Kernel};
use openshape_render::{Camera, Projection};
use openshape_sketch_solver::{SolveError, SolveStatus};

const MM: f64 = 1e-3;
/// Relative accuracy of geometry placed by clicking: screen positions are
/// `f32`, so clicked coordinates carry errors around 1e-8 m.
const CLICKED: f64 = 1e-6;

fn close(actual: f64, expected: f64) -> bool {
    (actual - expected).abs() <= CLICKED * expected.abs()
}

/// A top view showing about ±133 mm × ±100 mm in an 800 × 600 viewport
/// (one screen point is a third of a millimetre).
fn view() -> ScreenMap {
    let (yaw, pitch) = normal_view(&Frame::top());
    let camera = Camera {
        target: DVec3::ZERO,
        distance: 0.1 / (openshape_render::DEFAULT_FOV_Y_RAD as f64 * 0.5).tan(),
        yaw,
        pitch,
        projection: Projection::Orthographic,
        scene_bounds: None,
    };
    ScreenMap::new(
        &camera,
        Frame::top(),
        Rect::from_min_size(Pos2::new(10.0, 20.0), egui::vec2(800.0, 600.0)),
    )
}

/// Screen position of a plane point given in millimetres.
fn at(map: &ScreenMap, x: f64, y: f64) -> Pos2 {
    map.to_screen(DVec2::new(x, y) * MM).expect("on screen")
}

fn setup() -> (Model, SketchEditor, ScreenMap) {
    let mut model = Model::new();
    let f = model.add_sketch(PlaneRef::Top);
    (model, SketchEditor::new(f), view())
}

fn click(m: &mut Model, ed: &mut SketchEditor, map: &ScreenMap, x: f64, y: f64) {
    let p = at(map, x, y);
    ed.hover(m, map, Some(p));
    ed.click(m, map, p, false);
}

fn rectangle(m: &mut Model, ed: &mut SketchEditor, map: &ScreenMap, a: (f64, f64), b: (f64, f64)) {
    ed.set_tool(Tool::Rectangle);
    click(m, ed, map, a.0, a.1);
    click(m, ed, map, b.0, b.1);
    ed.set_tool(Tool::Select);
}

fn count<F: Fn(&Geometry<openshape_core_model::SketchId>) -> bool>(
    m: &Model,
    f: FeatureId,
    pred: F,
) -> usize {
    m.sketch_def(f)
        .unwrap()
        .sketch
        .entities
        .iter()
        .filter(|e| pred(&e.geometry))
        .count()
}

fn constraints<F: Fn(&Constraint<openshape_core_model::SketchId>) -> bool>(
    m: &Model,
    f: FeatureId,
    pred: F,
) -> usize {
    m.sketch_def(f)
        .unwrap()
        .sketch
        .constraints
        .iter()
        .filter(|c| pred(&c.constraint))
        .count()
}

fn status(m: &Model, f: FeatureId) -> SolveStatus {
    m.sketch_result(f).unwrap().solution.status.clone()
}

fn volume(m: &Model, feature: FeatureId) -> f64 {
    let b = m
        .bodies()
        .bodies
        .iter()
        .find(|b| b.feature == feature)
        .expect("the feature made a body");
    m.kernel().mass_properties(&b.body, 1e-9).unwrap().volume
}

fn extrude(m: &mut Model, map: &ScreenMap, inside: &[(f64, f64)], distance: &str) -> FeatureId {
    let mut pick = ProfilePick::new(SweepKind::Extrude);
    for (x, y) in inside {
        pick.click(m, map, at(map, *x, *y));
    }
    pick.distance = distance.into();
    let feature = pick.build().expect("extrude builds");
    m.add_feature("Extrude 1", feature)
}

#[test]
fn a_rectangle_is_closed_constrained_and_extrudes() {
    let (mut m, mut ed, map) = setup();
    let f = ed.feature;
    rectangle(&mut m, &mut ed, &map, (-20.0, -10.0), (20.0, 10.0));
    assert_eq!(count(&m, f, |g| matches!(g, Geometry::Point { .. })), 8);
    assert_eq!(count(&m, f, |g| matches!(g, Geometry::Line { .. })), 4);
    assert_eq!(
        constraints(&m, f, |c| matches!(c, Constraint::Coincident { .. })),
        4
    );
    assert_eq!(
        constraints(&m, f, |c| matches!(c, Constraint::Horizontal { .. })),
        2
    );
    assert_eq!(
        constraints(&m, f, |c| matches!(c, Constraint::Vertical { .. })),
        2
    );
    let r = m.sketch_result(f).unwrap();
    assert_eq!(r.solution.status, SolveStatus::Solved);
    // Two corners' positions are free.
    assert_eq!(r.solution.dof, 4);
    let regions = r.regions.as_ref().expect("regions");
    assert_eq!(regions.arrangement.regions.len(), 1);
    let area = regions.arrangement.regions[0].area;
    assert!(close(area, 800.0 * MM * MM), "{area}");

    let e = extrude(&mut m, &map, &[(0.0, 0.0)], "10");
    assert_eq!(m.status(e), Some(&FeatureStatus::Ok));
    let v = volume(&m, e);
    assert!(close(v, 8000.0 * MM.powi(3)), "{v}");
}

#[test]
fn a_dimension_drives_the_geometry_and_the_extrude_follows() {
    let (mut m, mut ed, map) = setup();
    let f = ed.feature;
    rectangle(&mut m, &mut ed, &map, (-20.0, -10.0), (20.0, 10.0));
    let e = extrude(&mut m, &map, &[(0.0, 0.0)], "10");
    let dof = m.sketch_result(f).unwrap().solution.dof;
    // Select the bottom edge and give it a length of 60 mm.
    click(&mut m, &mut ed, &map, 0.0, -10.0);
    assert_eq!(ed.selection.len(), 1);
    ed.dimension(&m, at(&map, 0.0, -20.0));
    let prompt = ed.prompt.as_ref().expect("dimension prompt");
    assert_eq!(prompt.text, "40");
    ed.prompt.as_mut().unwrap().text = "60".into();
    let ops = m.session().log().len();
    ed.commit_prompt(&mut m);
    assert!(ed.prompt.is_none());
    // One op: the dimension and the moved geometry together.
    assert_eq!(m.session().log().len(), ops + 1);
    let r = m.sketch_result(f).unwrap();
    assert_eq!(r.solution.status, SolveStatus::Solved);
    assert_eq!(r.solution.dof, dof - 1);
    let area = r.regions.as_ref().unwrap().arrangement.regions[0].area;
    assert!(close(area, 1200.0 * MM * MM), "{area}");
    // The extrude picked the region before the edit and still finds it.
    let v = volume(&m, e);
    assert!(close(v, 12000.0 * MM.powi(3)), "{v}");
    // Units: "2 cm" is accepted; a bad value is refused with a message.
    let dim = m
        .sketch_def(f)
        .unwrap()
        .sketch
        .constraints
        .iter()
        .find(|c| matches!(c.constraint, Constraint::Length { .. }))
        .unwrap()
        .id;
    ed.edit_dimension(&m, dim, at(&map, 0.0, 0.0));
    ed.prompt.as_mut().unwrap().text = "-3".into();
    ed.commit_prompt(&mut m);
    assert!(ed.prompt.as_ref().unwrap().error.is_some());
    ed.prompt.as_mut().unwrap().text = "5 cm".into();
    ed.commit_prompt(&mut m);
    let area = m
        .sketch_result(f)
        .unwrap()
        .regions
        .as_ref()
        .unwrap()
        .arrangement
        .regions[0]
        .area;
    assert!(close(area, 1000.0 * MM * MM), "{area}");
}

#[test]
fn lines_chain_infer_horizontal_and_snap_to_the_origin() {
    let (mut m, mut ed, map) = setup();
    let f = ed.feature;
    ed.set_tool(Tool::Line);
    click(&mut m, &mut ed, &map, 0.0, 0.0);
    // 1.9° off horizontal: inferred horizontal.
    click(&mut m, &mut ed, &map, 30.0, 1.0);
    click(&mut m, &mut ed, &map, 30.0, 25.0);
    ed.escape();
    assert_eq!(count(&m, f, |g| matches!(g, Geometry::Line { .. })), 2);
    assert_eq!(
        constraints(&m, f, |c| matches!(c, Constraint::Horizontal { .. })),
        1
    );
    // The second line was 0.0° from vertical.
    assert_eq!(
        constraints(&m, f, |c| matches!(c, Constraint::Vertical { .. })),
        1
    );
    // The chain joined the lines; the first click snapped to the origin.
    assert_eq!(
        constraints(&m, f, |c| matches!(c, Constraint::Coincident { .. })),
        1
    );
    assert_eq!(
        constraints(
            &m,
            f,
            |c| matches!(c, Constraint::Fix { at, .. } if *at == DVec2::ZERO)
        ),
        1
    );
    let r = m.sketch_result(f).unwrap();
    assert_eq!(r.solution.status, SolveStatus::Solved);
    let corner = r
        .solution
        .sketch
        .entities
        .iter()
        .filter_map(|e| match e.geometry {
            Geometry::Point { at } => Some(at),
            _ => None,
        })
        .find(|p| (p.x - 30.0 * MM).abs() < 1e-3 && p.y.abs() < 1e-3)
        .expect("corner");
    assert!(
        corner.y.abs() < 1e-12,
        "the first line is exactly horizontal"
    );
}

#[test]
fn dragging_moves_free_geometry_in_one_undo_step() {
    let (mut m, mut ed, map) = setup();
    let f = ed.feature;
    ed.set_tool(Tool::Line);
    click(&mut m, &mut ed, &map, 10.0, 10.0);
    click(&mut m, &mut ed, &map, 40.0, 30.0);
    ed.set_tool(Tool::Select);
    let before = m.sketch_def(f).unwrap().clone();
    let ops = m.session().log().len();
    ed.drag_start(&mut m, &map, at(&map, 40.0, 30.0));
    ed.drag_to(&m, &map, at(&map, 50.0, 35.0));
    ed.drag_to(&m, &map, at(&map, 60.0, 40.0));
    ed.drag_end(&mut m, &map, at(&map, 60.0, 40.0));
    assert_eq!(m.session().log().len(), ops + 1);
    let moved = m.sketch_result(f).unwrap().solution.sketch.clone();
    let has = |x: f64, y: f64| {
        moved.entities.iter().any(|e| {
            matches!(e.geometry, Geometry::Point { at }
                if (at - DVec2::new(x, y) * MM).length() < 0.5 * MM)
        })
    };
    assert!(has(60.0, 40.0) && has(10.0, 10.0));
    assert!(m.undo());
    assert_eq!(m.sketch_def(f).unwrap(), &before);
}

#[test]
fn deleting_a_line_removes_its_orphaned_points_and_constraints() {
    let (mut m, mut ed, map) = setup();
    let f = ed.feature;
    rectangle(&mut m, &mut ed, &map, (-20.0, -10.0), (20.0, 10.0));
    click(&mut m, &mut ed, &map, 0.0, 10.0);
    ed.delete_selection(&mut m);
    assert_eq!(count(&m, f, |g| matches!(g, Geometry::Line { .. })), 3);
    assert_eq!(count(&m, f, |g| matches!(g, Geometry::Point { .. })), 6);
    assert_eq!(
        constraints(&m, f, |c| matches!(c, Constraint::Coincident { .. })),
        2
    );
    assert_eq!(
        constraints(&m, f, |c| matches!(c, Constraint::Horizontal { .. })),
        1
    );
    // An open outline has no regions.
    let r = m.sketch_result(f).unwrap();
    assert!(r.regions.as_ref().unwrap().arrangement.regions.is_empty());
}

#[test]
fn a_conflicting_constraint_is_reported_and_undo_recovers() {
    let (mut m, mut ed, map) = setup();
    let f = ed.feature;
    ed.set_tool(Tool::Line);
    click(&mut m, &mut ed, &map, 0.0, 20.0);
    click(&mut m, &mut ed, &map, 40.0, 20.0);
    ed.escape();
    ed.set_tool(Tool::Select);
    click(&mut m, &mut ed, &map, 20.0, 20.0);
    assert!(ed.can_constrain(&m, ConstraintTool::Vertical));
    ed.constrain(&mut m, ConstraintTool::Vertical);
    assert_eq!(status(&m, f), SolveStatus::Failed(SolveError::Conflicting));
    assert!(matches!(m.status(f), Some(FeatureStatus::Failed(_))));
    assert!(ed.status(&m).contains("Conflicting"));
    assert!(ed.message.as_deref().unwrap_or("").contains("conflict"));
    let d = ed.display(&m).unwrap();
    assert!(!d.conflicting.is_empty());
    assert!(m.undo());
    assert_eq!(status(&m, f), SolveStatus::Solved);
}

#[test]
fn constraint_tools_need_the_right_selection() {
    let (mut m, mut ed, map) = setup();
    rectangle(&mut m, &mut ed, &map, (-20.0, -10.0), (20.0, 10.0));
    // Nothing selected: nothing applies.
    for t in ConstraintTool::ALL {
        assert!(!ed.can_constrain(&m, t), "{t:?}");
    }
    // Two lines: the line-pair tools apply (whether the result is sensible
    // is the solver's call); the circle tools do not.
    click(&mut m, &mut ed, &map, 0.0, -10.0);
    ed.click(&mut m, &map, at(&map, 0.0, 10.0), true);
    assert_eq!(ed.selection.len(), 2);
    assert!(ed.can_constrain(&m, ConstraintTool::Parallel));
    assert!(ed.can_constrain(&m, ConstraintTool::Equal));
    assert!(ed.can_constrain(&m, ConstraintTool::Perpendicular));
    assert!(!ed.can_constrain(&m, ConstraintTool::Concentric));
    assert!(!ed.can_constrain(&m, ConstraintTool::Tangent));
    // Opposite sides of a rectangle are already equal: redundant.
    ed.constrain(&mut m, ConstraintTool::Equal);
    assert_eq!(status(&m, ed.feature), SolveStatus::Redundant);
    assert!(matches!(
        m.status(ed.feature),
        Some(FeatureStatus::Warning(_))
    ));
}

#[test]
fn circles_arcs_and_splines() {
    let (mut m, mut ed, map) = setup();
    let f = ed.feature;
    ed.set_tool(Tool::Circle);
    click(&mut m, &mut ed, &map, 50.0, 50.0);
    click(&mut m, &mut ed, &map, 60.0, 50.0);
    let r = m.sketch_result(f).unwrap();
    let regions = r.regions.as_ref().unwrap();
    assert_eq!(regions.arrangement.regions.len(), 1);
    let area = regions.arrangement.regions[0].area;
    assert!(
        close(area, std::f64::consts::PI * 100.0 * MM * MM),
        "{area}"
    );

    // A 3-point arc through the top runs counter-clockwise from the right.
    ed.set_tool(Tool::Arc);
    click(&mut m, &mut ed, &map, -40.0, 0.0);
    click(&mut m, &mut ed, &map, -60.0, 0.0);
    click(&mut m, &mut ed, &map, -50.0, 10.0);
    let sketch = &m.sketch_result(f).unwrap().solution.sketch;
    let arc = sketch
        .entities
        .iter()
        .find_map(|e| match e.geometry {
            Geometry::Arc { center, start, end } => Some((center, start, end)),
            _ => None,
        })
        .expect("arc");
    let p = |id| sketch.point(id).unwrap() / MM;
    assert!((p(arc.0) - DVec2::new(-50.0, 0.0)).length() < 1e-3);
    assert!((p(arc.1) - DVec2::new(-40.0, 0.0)).length() < 1e-3);
    assert!((p(arc.2) - DVec2::new(-60.0, 0.0)).length() < 1e-3);

    ed.set_tool(Tool::Spline);
    click(&mut m, &mut ed, &map, 0.0, -40.0);
    click(&mut m, &mut ed, &map, 10.0, -30.0);
    click(&mut m, &mut ed, &map, 20.0, -40.0);
    ed.finish(&mut m);
    let spline = m
        .sketch_def(f)
        .unwrap()
        .sketch
        .entities
        .iter()
        .find_map(|e| match &e.geometry {
            Geometry::InterpolatedSpline { through, params } => {
                Some((through.len(), params.clone()))
            }
            _ => None,
        })
        .expect("spline");
    assert_eq!(spline.0, 3);
    assert_eq!(spline.1.first(), Some(&0.0));
    assert_eq!(spline.1.last(), Some(&1.0));
    assert!((spline.1[1] - 0.5).abs() < 1e-3);
    assert_eq!(status(&m, f), SolveStatus::Solved);
}

#[test]
fn revolve_about_a_construction_line() {
    let (mut m, mut ed, map) = setup();
    let f = ed.feature;
    rectangle(&mut m, &mut ed, &map, (10.0, 0.0), (20.0, 30.0));
    // A construction line on the sketch's y axis.
    ed.construction = true;
    ed.set_tool(Tool::Line);
    click(&mut m, &mut ed, &map, 0.0, -40.0);
    click(&mut m, &mut ed, &map, 0.0, 40.0);
    ed.escape();
    ed.construction = false;
    let r = m.sketch_result(f).unwrap();
    // Construction geometry bounds no region.
    assert_eq!(r.regions.as_ref().unwrap().arrangement.regions.len(), 1);
    let mut pick = ProfilePick::new(SweepKind::Revolve);
    pick.click(&m, &map, at(&map, 15.0, 15.0));
    pick.picking_axis = true;
    pick.click(&m, &map, at(&map, 0.0, 20.0));
    assert!(pick.axis.is_some() && !pick.picking_axis);
    let feature = pick.build().expect("revolve builds");
    let rev = m.add_feature("Revolve 1", feature);
    assert_eq!(m.status(rev), Some(&FeatureStatus::Ok));
    // Pappus: 2π × centroid radius × area.
    let expected = std::f64::consts::TAU * 15.0 * 300.0 * MM.powi(3);
    let v = volume(&m, rev);
    assert!(close(v, expected), "{v} vs {expected}");
}

#[test]
fn undoing_the_sketch_removes_it() {
    let (mut m, mut ed, map) = setup();
    rectangle(&mut m, &mut ed, &map, (-20.0, -10.0), (20.0, 10.0));
    assert!(m.undo());
    assert!(m.undo());
    assert!(m.sketch_result(ed.feature).is_none());
}

/// Opens the dimension prompt for the current selection, types `value` and
/// commits it; panics if the prompt reports an error.
fn dimension(m: &mut Model, ed: &mut SketchEditor, map: &ScreenMap, shown: &str, value: &str) {
    ed.dimension(m, at(map, 90.0, 90.0));
    let prompt = ed.prompt.as_mut().expect("dimension prompt");
    assert_eq!(prompt.text, shown);
    prompt.text = value.into();
    ed.commit_prompt(m);
    assert!(ed.prompt.is_none(), "{:?}", ed.prompt);
}

fn solved(m: &Model, f: FeatureId) -> openshape_core_model::SketchData {
    let r = m.sketch_result(f).unwrap();
    assert_eq!(r.solution.status, SolveStatus::Solved);
    r.solution.sketch.clone()
}

#[test]
fn every_dimension_kind_is_placed_from_a_selection_and_drives_the_geometry() {
    let (mut m, mut ed, map) = setup();
    let f = ed.feature;

    // A circle: diameter.
    ed.set_tool(Tool::Circle);
    click(&mut m, &mut ed, &map, 40.0, 40.0);
    click(&mut m, &mut ed, &map, 55.0, 40.0);
    ed.set_tool(Tool::Select);
    click(&mut m, &mut ed, &map, 40.0, 55.0);
    dimension(&mut m, &mut ed, &map, "30", "50");
    let s = solved(&m, f);
    let radius = s
        .entities
        .iter()
        .find_map(|e| match e.geometry {
            Geometry::Circle { radius, .. } => Some(radius),
            _ => None,
        })
        .unwrap();
    assert!(close(radius, 25.0 * MM), "{radius}");
    assert_eq!(
        constraints(&m, f, |c| matches!(c, Constraint::Diameter { .. })),
        1
    );

    // An arc (centre (-50, 0), radius 10): radius.
    ed.set_tool(Tool::Arc);
    click(&mut m, &mut ed, &map, -40.0, 0.0);
    click(&mut m, &mut ed, &map, -60.0, 0.0);
    click(&mut m, &mut ed, &map, -50.0, 10.0);
    ed.set_tool(Tool::Select);
    click(&mut m, &mut ed, &map, -50.0, 10.0);
    dimension(&mut m, &mut ed, &map, "10", "12 mm");
    let s = solved(&m, f);
    let (c, a) = s
        .entities
        .iter()
        .find_map(|e| match e.geometry {
            Geometry::Arc { center, start, .. } => Some((center, start)),
            _ => None,
        })
        .unwrap();
    let r = (s.point(a).unwrap() - s.point(c).unwrap()).length();
    assert!(close(r, 12.0 * MM), "{r}");

    // Two lines from the origin at 45° to each other: angle.
    ed.set_tool(Tool::Line);
    click(&mut m, &mut ed, &map, 0.0, -60.0);
    click(&mut m, &mut ed, &map, 40.0, -60.0);
    ed.escape();
    click(&mut m, &mut ed, &map, 0.0, -60.0);
    click(&mut m, &mut ed, &map, 30.0, -30.0);
    ed.escape();
    ed.set_tool(Tool::Select);
    click(&mut m, &mut ed, &map, 20.0, -60.0);
    ed.click(&mut m, &map, at(&map, 15.0, -45.0), true);
    assert_eq!(ed.selection.len(), 2);
    dimension(&mut m, &mut ed, &map, "45", "60");
    let s = solved(&m, f);
    let lines: Vec<DVec2> = s
        .entities
        .iter()
        .filter_map(|e| match e.geometry {
            Geometry::Line { start, end } => Some(s.point(end)? - s.point(start)?),
            _ => None,
        })
        .collect();
    assert_eq!(lines.len(), 2);
    let angle = lines[0].angle_to(lines[1]).abs().to_degrees();
    assert!((angle - 60.0).abs() < 1e-6, "{angle}");

    // Two points: aligned by default, switched to a horizontal distance.
    ed.set_tool(Tool::Point);
    click(&mut m, &mut ed, &map, -80.0, 60.0);
    click(&mut m, &mut ed, &map, -50.0, 80.0);
    ed.set_tool(Tool::Select);
    click(&mut m, &mut ed, &map, -80.0, 60.0);
    ed.click(&mut m, &map, at(&map, -50.0, 80.0), true);
    ed.dimension(&m, at(&map, 90.0, 90.0));
    ed.set_prompt_kind(&m, openshape_core_model::sketch::DistanceKind::Horizontal);
    let prompt = ed.prompt.as_mut().unwrap();
    assert_eq!(prompt.text, "30");
    prompt.text = "45".into();
    ed.commit_prompt(&mut m);
    assert!(ed.prompt.is_none());
    let s = solved(&m, f);
    let xs: Vec<f64> = s
        .entities
        .iter()
        .filter_map(|e| match e.geometry {
            Geometry::Point { at } if at.y > 50.0 * MM => Some(at.x),
            _ => None,
        })
        .collect();
    assert_eq!(xs.len(), 2);
    assert!(close((xs[0] - xs[1]).abs(), 45.0 * MM), "{xs:?}");
    assert_eq!(
        constraints(&m, f, |c| matches!(
            c,
            Constraint::Distance {
                kind: openshape_core_model::sketch::DistanceKind::Horizontal,
                ..
            }
        )),
        1
    );
}
