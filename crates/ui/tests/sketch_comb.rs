//! Curvature combs on sketch curves: which curves take one, showing and
//! hiding them, what they show across a curvature-continuous (G2) joint
//! and across a merely tangent one, and the right-click menu, driven
//! offscreen. Joint curvatures are checked against the solved curves'
//! own end-derivative formulas (as in `sketch_curvature.rs`), not against
//! the sampler that draws the combs.

mod common;

use std::collections::BTreeSet;

use common::*;
use egui::Pos2;
use openshape_app::comb::{self, LONGEST_TOOTH};
use openshape_app::commands::Command;
use openshape_app::screenshots::AppDriver;
use openshape_app::sketch_mode::{ConstraintTool, ScreenMap, SketchEditor, Tool};
use openshape_app::Model;
use openshape_core_model::sketch::Geometry;
use openshape_core_model::{DVec2, FeatureId, SketchId};
use openshape_kernel::{CurvatureSample, Frame};
use openshape_render::RgbaImage;
use openshape_sketch_solver::{ResolvedCurve, SolveStatus};
use openshape_ui_kit::theme::Theme;

/// The id of the sketch's first entity matching `pred`.
fn id<F: Fn(&Geometry<SketchId>) -> bool>(m: &Model, f: FeatureId, pred: F) -> SketchId {
    m.sketch_def(f)
        .unwrap()
        .sketch
        .entities
        .iter()
        .find(|e| pred(&e.geometry))
        .map(|e| e.id)
        .expect("an entity")
}

/// A counter-clockwise arc about (0, 20) from (−20, 20) down to (0, 0),
/// turning left, and a cubic spline leaving (0, 0) up and to the right and
/// turning right; their ids (arc, spline).
fn arc_and_spline() -> (Model, SketchEditor, FeatureId, SketchId, SketchId) {
    let (mut m, mut ed, map) = setup();
    let f = ed.feature;
    ed.set_tool(Tool::CenterArc);
    click(&mut m, &mut ed, &map, 0.0, 20.0);
    click(&mut m, &mut ed, &map, -20.0, 20.0);
    hover(&m, &mut ed, &map, -14.0, 6.0);
    click(&mut m, &mut ed, &map, 0.0, 0.0);
    ed.set_tool(Tool::ControlSpline);
    for (x, y) in [(0.0, 0.0), (15.0, 6.0), (30.0, 2.0), (45.0, -4.0)] {
        click(&mut m, &mut ed, &map, x, y);
    }
    ed.finish(&mut m);
    ed.set_tool(Tool::Select);
    assert_eq!(status(&m, f), SolveStatus::Solved);
    let arc = id(&m, f, |g| matches!(g, Geometry::Arc { .. }));
    let spline = id(&m, f, |g| matches!(g, Geometry::ControlSpline { .. }));
    (m, ed, f, arc, spline)
}

/// The solved curve `id`.
fn resolved(m: &Model, f: FeatureId, id: SketchId) -> ResolvedCurve {
    curves(m, f, true)
        .into_iter()
        .find(|(k, _)| *k == id)
        .map(|(_, c)| c)
        .expect("a solved curve")
}

/// A spline's signed curvature at its start, from its control points:
/// `C′(0) = Q₀`, `C″(0) = (p − 1)(Q₁ − Q₀)/(u_{p+1} − u₂)` with
/// `Qᵢ = p(Pᵢ₊₁ − Pᵢ)/(u_{i+p+1} − u_{i+1})` (Piegl & Tiller eq. 3.8).
fn spline_start_curvature(c: &ResolvedCurve) -> f64 {
    let ResolvedCurve::Spline {
        degree,
        knots: u,
        points: p,
    } = c
    else {
        panic!("not a spline: {c:?}");
    };
    let d = *degree;
    let q = |i: usize| (p[i + 1] - p[i]) * (d as f64 / (u[i + d + 1] - u[i + 1]));
    let c1 = q(0);
    let c2 = (q(1) - q(0)) * ((d as f64 - 1.0) / (u[d + 1] - u[2]));
    c1.perp_dot(c2) / c1.length().powi(3)
}

/// The comb samples shown for `id`.
fn shown(m: &Model, ed: &SketchEditor, id: SketchId) -> Vec<CurvatureSample> {
    ed.display(m)
        .expect("a sketch display")
        .combs
        .into_iter()
        .find(|(k, _)| *k == id)
        .map(|(_, s)| s)
        .unwrap_or_default()
}

/// The arc's comb sample where it ends and the spline's where it starts:
/// the two sides of their joint.
fn joint(m: &Model, ed: &SketchEditor, arc: SketchId, spline: SketchId) -> [CurvatureSample; 2] {
    let a = shown(m, ed, arc);
    let s = shown(m, ed, spline);
    assert_eq!(a.len(), comb::TEETH + 1);
    assert!(s.len() > comb::TEETH, "{}", s.len());
    [*a.last().unwrap(), s[0]]
}

#[test]
fn combs_carry_on_across_a_curvature_continuous_joint() {
    let (mut m, mut ed, f, arc, spline) = arc_and_spline();
    ed.selection = vec![spline, arc];
    ed.constrain(&mut m, ConstraintTool::Curvature);
    assert_eq!(status(&m, f), SolveStatus::Solved, "{:?}", ed.message);
    ed.selection = vec![spline, arc];
    ed.toggle_combs(&m);
    assert_eq!(ed.combs, BTreeSet::from([arc, spline]));

    let [a, s] = joint(&m, &ed, arc, spline);
    assert!((a.point - s.point).length() < 1e-12, "{a:?} {s:?}");
    assert!(a.tangent.dot(s.tangent) > 1.0 - 1e-12, "{a:?} {s:?}");
    // The arc's curvature, 1/r, carried on by the spline: G2.
    let ResolvedCurve::Arc { radius, .. } = resolved(&m, f, arc) else {
        panic!("not an arc");
    };
    let k = spline_start_curvature(&resolved(&m, f, spline));
    assert!((a.curvature * radius - 1.0).abs() < 1e-12, "{a:?} {radius}");
    assert!((s.curvature / k - 1.0).abs() < 1e-9, "{s:?} {k}");
    assert!(
        (s.curvature / a.curvature - 1.0).abs() < 1e-9,
        "{a:?} {s:?}"
    );
    // So the teeth there are one tooth, on the outside of the turn.
    let scale = 1e-4;
    let (ta, ts) = (comb::comb(&[a], scale), comb::comb(&[s], scale));
    let (da, ds) = (
        ta.teeth[0][1] - ta.teeth[0][0],
        ts.teeth[0][1] - ts.teeth[0][0],
    );
    assert!((da - ds).length() < 1e-9 * da.length(), "{da} {ds}");
    let centre = a.point + a.tangent.perp() * radius;
    assert!(
        da.dot(centre - a.point) < 0.0,
        "the tooth points away from the centre"
    );
}

#[test]
fn combs_jump_where_curves_are_only_tangent() {
    let (mut m, mut ed, f, arc, spline) = arc_and_spline();
    ed.selection = vec![spline, arc];
    ed.constrain(&mut m, ConstraintTool::Tangent);
    assert_eq!(status(&m, f), SolveStatus::Solved, "{:?}", ed.message);
    ed.selection = vec![spline, arc];
    ed.toggle_combs(&m);
    let [a, s] = joint(&m, &ed, arc, spline);
    assert!((a.point - s.point).length() < 1e-12, "{a:?} {s:?}");
    assert!(a.tangent.dot(s.tangent) > 1.0 - 1e-9, "{a:?} {s:?}");
    // One direction, but the spline still bends the other way: the teeth
    // jump to the other side.
    let k = spline_start_curvature(&resolved(&m, f, spline));
    assert!((s.curvature / k - 1.0).abs() < 1e-9, "{s:?} {k}");
    assert!(a.curvature > 0.0 && s.curvature < 0.0, "{a:?} {s:?}");
}

#[test]
fn combs_are_shown_and_hidden_for_curves_that_bend() {
    let (mut m, mut ed, f, arc, spline) = arc_and_spline();
    let map = view();
    draw(
        &mut m,
        &mut ed,
        &map,
        Tool::Line,
        &[(-20.0, 20.0), (-20.0, 40.0)],
    );
    let line = id(&m, f, |g| matches!(g, Geometry::Line { .. }));
    let refused = |m: &Model, ed: &mut SketchEditor, sel: Vec<SketchId>| {
        ed.selection = sel;
        ed.message = None;
        assert!(ed.comb_candidates(m).is_empty());
        assert!(!ed.combs_shown(m));
        ed.toggle_combs(m);
        assert!(ed.combs.is_empty());
        let message = ed.message.clone().unwrap_or_default();
        assert!(message.contains("Select an arc"), "{message}");
    };
    // Nothing selected, or only a line or points: nothing bends.
    refused(&m, &mut ed, vec![]);
    refused(&m, &mut ed, vec![line]);
    let points: Vec<SketchId> = m
        .sketch_def(f)
        .unwrap()
        .sketch
        .entities
        .iter()
        .filter(|e| matches!(e.geometry, Geometry::Point { .. }))
        .map(|e| e.id)
        .collect();
    refused(&m, &mut ed, points);
    // A line and the arc: the arc's comb only.
    ed.selection = vec![line, arc];
    ed.toggle_combs(&m);
    assert_eq!(ed.combs, BTreeSet::from([arc]));
    assert!(ed.combs_shown(&m));
    assert_eq!(shown(&m, &ed, arc).len(), comb::TEETH + 1);
    assert!(shown(&m, &ed, spline).is_empty());
    // The arc and the spline, one not shown yet: both show.
    ed.selection = vec![arc, spline];
    assert!(!ed.combs_shown(&m));
    ed.toggle_combs(&m);
    assert_eq!(ed.combs, BTreeSet::from([arc, spline]));
    // A comb follows its curve through edits: after dragging the arc's
    // centre, it shows the arc as solved now.
    let before = shown(&m, &ed, arc)[0];
    let centre = at(&map, 0.0, 20.0);
    ed.selection.clear();
    ed.hover(&m, &map, Some(centre));
    ed.drag_start(&mut m, &map, centre);
    ed.drag_to(&m, &map, at(&map, 0.0, 23.0));
    ed.drag_end(&mut m, &map, at(&map, 0.0, 25.0));
    let after = shown(&m, &ed, arc)[0];
    let ResolvedCurve::Arc { center, radius, .. } = resolved(&m, f, arc) else {
        panic!("not an arc");
    };
    assert!(
        (after.point - before.point).length() > 1e-3,
        "{before:?} {after:?}"
    );
    assert!(((after.point - center).length() / radius - 1.0).abs() < 1e-12);
    assert!(
        (after.curvature * radius - 1.0).abs() < 1e-12,
        "{after:?} {radius}"
    );
    // Both shown: the same command hides both.
    ed.selection = vec![arc, spline];
    assert!(ed.combs_shown(&m));
    ed.toggle_combs(&m);
    assert!(ed.combs.is_empty());
    assert!(ed.display(&m).unwrap().combs.is_empty());
}

/// Screen positions of a plane point (millimetres) in the app's sketch view.
fn on_plane(d: &mut AppDriver, x: f64, y: f64) -> Pos2 {
    let rect = d.app.layout().graphics;
    let camera = d.app.camera_mut().clone();
    ScreenMap::new(&camera, Frame::top(), rect)
        .to_screen(DVec2::new(x, y) * MM)
        .expect("on screen")
}

/// The share of points around a circle `radius` points about `centre`
/// (those inside the graphics area) that have a pixel of the comb's
/// envelope colour within two pixels.
fn envelope_share(d: &AppDriver, img: &RgbaImage, centre: Pos2, radius: f32) -> f64 {
    let want = Theme::light().sketch.comb_envelope;
    let rect = d.app.layout().graphics.shrink(4.0);
    let (mut seen, mut hit) = (0, 0);
    for k in 0..72 {
        let a = k as f32 * 5f32.to_radians();
        let p = centre + egui::vec2(a.cos(), a.sin()) * radius;
        if !rect.contains(p) {
            continue;
        }
        seen += 1;
        let near = (-2..=2).any(|dx| {
            (-2..=2).any(|dy| {
                let q = img.pixel((p.x as i32 + dx) as u32, (p.y as i32 + dy) as u32);
                (q[0] as i32 - want.r() as i32).abs() < 40
                    && (q[1] as i32 - want.g() as i32).abs() < 40
                    && (q[2] as i32 - want.b() as i32).abs() < 40
            })
        });
        if near {
            hit += 1;
        }
    }
    assert!(seen >= 24, "only {seen} points of the envelope on screen");
    hit as f64 / seen as f64
}

#[test]
fn the_sketch_menu_shows_and_hides_a_comb() {
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
        &[(-60.0, -40.0), (-30.0, -40.0)],
    );
    let circle = id(&m, f, |g| matches!(g, Geometry::Circle { .. }));
    let line = id(&m, f, |g| matches!(g, Geometry::Line { .. }));
    let mut d = AppDriver::new([1200, 800], 1.0, m).expect("an adapter (Vulkan or GL)");
    d.app.enter_sketch(f);
    d.settle().unwrap();
    let centre = on_plane(&mut d, 0.0, 0.0);
    let r = (on_plane(&mut d, 20.0, 0.0) - centre).length();
    let envelope = r + LONGEST_TOOTH;
    let img = d.settle().unwrap();
    assert!(envelope_share(&d, &img, centre, envelope) < 0.1);

    // A line has no curvature: the menu offers the comb but cannot show it.
    let on_line = on_plane(&mut d, -45.0, -40.0);
    d.right_click(on_line).unwrap();
    d.settle().unwrap();
    assert_eq!(d.app.sketch_editor().unwrap().selection, vec![line]);
    assert!(d.menu_row("Show curvature comb").is_some());
    assert!(!d.app.enabled(Command::CurvatureComb));
    d.key(egui::Key::Escape, egui::Modifiers::NONE).unwrap();
    d.settle().unwrap();
    // Escape closes the menu, and only the menu.
    assert!(d.menu_row("Show curvature comb").is_none());
    assert_eq!(d.app.sketch_editor().unwrap().selection, vec![line]);

    // Right-clicking the circle selects it; the menu shows its comb: teeth
    // of the longest length all round, their tips on a circle that much
    // outside it.
    let on_circle = on_plane(&mut d, 0.0, -20.0);
    d.right_click(on_circle).unwrap();
    d.settle().unwrap();
    assert_eq!(d.app.sketch_editor().unwrap().selection, vec![circle]);
    assert!(d.app.enabled(Command::CurvatureComb));
    let show = d.menu_row("Show curvature comb").expect("the sketch menu");
    d.click(show.center()).unwrap();
    let img = d.settle().unwrap();
    assert!(
        d.menu_row("Show curvature comb").is_none(),
        "the menu closed"
    );
    assert_eq!(
        d.app.sketch_editor().unwrap().combs,
        BTreeSet::from([circle])
    );
    assert!(d.app.active(Command::CurvatureComb));
    let share = envelope_share(&d, &img, centre, envelope);
    assert!(share > 0.9, "envelope found at {share} of the points");

    // The menu now offers to hide it.
    d.right_click(on_circle).unwrap();
    d.settle().unwrap();
    assert!(d.menu_row("Show curvature comb").is_none());
    let hide = d.menu_row("Hide curvature comb").expect("the sketch menu");
    d.click(hide.center()).unwrap();
    let img = d.settle().unwrap();
    assert!(d.app.sketch_editor().unwrap().combs.is_empty());
    assert!(envelope_share(&d, &img, centre, envelope) < 0.1);
}
