//! Helpers for driving sketch mode headlessly through a top-view camera,
//! shared by the sketch tool tests.

#![allow(dead_code)]

use egui::{Pos2, Rect};
use opendatum_app::profile_pick::{ProfilePick, SweepKind};
use opendatum_app::sketch_mode::{normal_view, ScreenMap, SketchEditor, Tool};
use opendatum_app::Model;
use opendatum_core_model::sketch::{Constraint, Geometry};
use opendatum_core_model::{DVec2, DVec3, FeatureId, PlaneRef, SketchId};
use opendatum_kernel::{Frame, Kernel};
use opendatum_render::{Camera, Projection};
use opendatum_sketch_solver::{ResolvedCurve, SolveStatus};

/// A millimetre in metres.
pub const MM: f64 = 1e-3;
/// Relative accuracy of geometry placed by clicking: screen positions are
/// `f32`, so clicked coordinates carry errors around 1e-8 m.
pub const CLICKED: f64 = 1e-6;

/// `actual` within [`CLICKED`] (relative) of `expected`.
pub fn close(actual: f64, expected: f64) -> bool {
    (actual - expected).abs() <= CLICKED * expected.abs()
}

/// A top view showing about ±133 mm × ±100 mm in an 800 × 600 viewport
/// (one screen point is a third of a millimetre).
pub fn view() -> ScreenMap {
    let (yaw, pitch) = normal_view(&Frame::top());
    let camera = Camera {
        target: DVec3::ZERO,
        distance: 0.1 / (opendatum_render::DEFAULT_FOV_Y_RAD as f64 * 0.5).tan(),
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
pub fn at(map: &ScreenMap, x: f64, y: f64) -> Pos2 {
    map.to_screen(DVec2::new(x, y) * MM).expect("on screen")
}

/// A model with one empty sketch on the top plane, its editor and the view.
pub fn setup() -> (Model, SketchEditor, ScreenMap) {
    let mut model = Model::new();
    let f = model.add_sketch(PlaneRef::Top);
    (model, SketchEditor::new(f), view())
}

/// Hovers and clicks at a plane point given in millimetres.
pub fn click(m: &mut Model, ed: &mut SketchEditor, map: &ScreenMap, x: f64, y: f64) {
    let p = at(map, x, y);
    ed.hover(m, map, Some(p));
    ed.click(m, map, p, false);
}

/// Hovers at a plane point given in millimetres.
pub fn hover(m: &Model, ed: &mut SketchEditor, map: &ScreenMap, x: f64, y: f64) {
    ed.hover(m, map, Some(at(map, x, y)));
}

/// Uses `tool` for the clicks `points` (millimetres), then returns to Select.
pub fn draw(
    m: &mut Model,
    ed: &mut SketchEditor,
    map: &ScreenMap,
    tool: Tool,
    points: &[(f64, f64)],
) {
    ed.set_tool(tool);
    for (x, y) in points {
        click(m, ed, map, *x, *y);
    }
    ed.set_tool(Tool::Select);
}

/// Number of entities of the sketch matching `pred`.
pub fn count<F: Fn(&Geometry<SketchId>) -> bool>(m: &Model, f: FeatureId, pred: F) -> usize {
    m.sketch_def(f)
        .unwrap()
        .sketch
        .entities
        .iter()
        .filter(|e| pred(&e.geometry))
        .count()
}

/// Number of construction curves of the sketch.
pub fn construction(m: &Model, f: FeatureId) -> usize {
    m.sketch_def(f)
        .unwrap()
        .sketch
        .entities
        .iter()
        .filter(|e| e.construction && e.geometry.is_curve())
        .count()
}

/// Number of constraints of the sketch matching `pred`.
pub fn constraints<F: Fn(&Constraint<SketchId>) -> bool>(
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

/// The solver status of a sketch.
pub fn status(m: &Model, f: FeatureId) -> SolveStatus {
    m.sketch_result(f).unwrap().solution.status.clone()
}

/// The degrees of freedom left in a sketch.
pub fn dof(m: &Model, f: FeatureId) -> usize {
    m.sketch_result(f).unwrap().solution.dof
}

/// The areas of a sketch's regions, largest first.
pub fn region_areas(m: &Model, f: FeatureId) -> Vec<f64> {
    let mut a: Vec<f64> = m
        .sketch_result(f)
        .unwrap()
        .regions
        .as_ref()
        .expect("regions")
        .arrangement
        .regions
        .iter()
        .map(|r| r.area)
        .collect();
    a.sort_by(|x, y| y.total_cmp(x));
    a
}

/// The solved curves of a sketch (non-construction unless `all`).
pub fn curves(m: &Model, f: FeatureId, all: bool) -> Vec<(SketchId, ResolvedCurve)> {
    m.sketch_result(f)
        .unwrap()
        .curves
        .iter()
        .filter(|(_, c, _)| all || !c)
        .map(|(id, _, c)| (*id, c.clone()))
        .collect()
}

/// The volume of the body a feature made.
pub fn volume(m: &Model, feature: FeatureId) -> f64 {
    let b = m
        .bodies()
        .bodies
        .iter()
        .find(|b| b.feature == feature)
        .expect("the feature made a body");
    m.kernel().mass_properties(&b.body, 1e-9).unwrap().volume
}

/// Extrudes the regions containing the given points (millimetres).
pub fn extrude(m: &mut Model, map: &ScreenMap, inside: &[(f64, f64)], distance: &str) -> FeatureId {
    let mut pick = ProfilePick::new(SweepKind::Extrude);
    for (x, y) in inside {
        pick.click(m, map, at(map, *x, *y));
    }
    pick.distance = distance.into();
    let feature = pick.build().expect("extrude builds");
    m.add_feature("Extrude 1", feature)
}
