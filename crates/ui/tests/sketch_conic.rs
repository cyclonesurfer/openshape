//! The conic tool, the ρ dimension and dragging ρ (design note 0009 §7–§8),
//! driven by clicks through a top-view camera. Points are checked against
//! Lee's implicit equation of the control triangle, and areas and volumes
//! against the closed form of a conic segment, r(w) times its triangle.

mod common;

use common::*;
use openshape_app::sketch_mode::{
    annotations, sample_resolved, ConstraintTool, SketchEditor, SnapTarget, Tool,
};
use openshape_app::sketch_ops::trim::{self, Preview};
use openshape_app::sketch_ops::Batch;
use openshape_app::Model;
use openshape_core_model::sketch::{Constraint, ConstraintEntry, Geometry, TangentKind};
use openshape_core_model::{DVec2, FeatureId, SketchEdit, SketchId};
use openshape_sketch_solver::conic::point as conic_point;
use openshape_sketch_solver::{ResolvedCurve, SolveStatus};

/// `a` within `rel` (relative) of `b`.
fn near(a: f64, b: f64, rel: f64) -> bool {
    (a - b).abs() <= rel * b.abs()
}

fn is_conic(g: &Geometry<SketchId>) -> bool {
    matches!(g, Geometry::Conic { .. })
}

/// The ids of the sketch's conics, in order.
fn conic_ids(m: &Model, f: FeatureId) -> Vec<SketchId> {
    m.sketch_def(f)
        .unwrap()
        .sketch
        .entities
        .iter()
        .filter(|e| is_conic(&e.geometry))
        .map(|e| e.id)
        .collect()
}

/// A solved conic: its control points and weight.
type Solved = ([DVec2; 3], f64);

fn conics(m: &Model, f: FeatureId) -> Vec<Solved> {
    curves(m, f, true)
        .into_iter()
        .filter_map(|(_, c)| match c {
            ResolvedCurve::Conic {
                start,
                control,
                end,
                weight,
            } => Some(([start, control, end], weight)),
            _ => None,
        })
        .collect()
}

/// The ρ stored on each conic of the solved sketch, in order.
fn rhos(m: &Model, f: FeatureId) -> Vec<f64> {
    m.sketch_result(f)
        .unwrap()
        .solution
        .sketch
        .entities
        .iter()
        .filter_map(|e| match e.geometry {
            Geometry::Conic { rho, .. } => Some(rho),
            _ => None,
        })
        .collect()
}

/// Lee's implicit equation `τ₁² − 4w²τ₀τ₂` in the barycentric coordinates
/// of the control triangle, scaled to be comparable across weights, and
/// the smallest coordinate (every one is non-negative on the arc).
fn implicit(c: &Solved, x: DVec2) -> (f64, f64) {
    let (p, w) = (c.0, c.1);
    let area = (p[1] - p[0]).perp_dot(p[2] - p[0]);
    let t0 = (p[1] - x).perp_dot(p[2] - x) / area;
    let t1 = (p[2] - x).perp_dot(p[0] - x) / area;
    let t2 = (p[0] - x).perp_dot(p[1] - x) / area;
    (
        (t1 * t1 - 4.0 * w * w * t0 * t2) / (1.0 + w * w),
        t0.min(t1).min(t2),
    )
}

/// arcosh w for w ≥ 1 as ln(1 + t + √(t(t + 2))), t = w − 1 (exact for
/// w ≤ 2), accurate to a few units in the last place near 1 whichever
/// toolchain builds the tests. `f64::acosh` is not: Rust 1.95's computes
/// ln(w + √(w − 1)√(w + 1)) and is 25 units in the last place off at
/// w = 1 + 1e-5 (1.97's is correctly rounded there), which the cancelling
/// segment formula below turns into a 4e-10 relative error.
fn acosh(w: f64) -> f64 {
    let t = w - 1.0;
    (t + (t * (t + 2.0)).sqrt()).ln_1p()
}

/// The area between a conic of weight `w` and its chord over its control
/// triangle's (design note 0009 §9).
fn ratio(w: f64) -> f64 {
    if w < 1.0 {
        let s = ((1.0 - w) * (1.0 + w)).sqrt();
        w * (w.acos() - w * s) / (s * s * s)
    } else if w == 1.0 {
        2.0 / 3.0
    } else {
        let s = ((w - 1.0) * (w + 1.0)).sqrt();
        w * (w * s - acosh(w)) / (s * s * s)
    }
}

/// The area of a solved conic's segment.
fn segment(c: &Solved) -> f64 {
    let p = c.0;
    ratio(c.1) * 0.5 * (p[1] - p[0]).perp_dot(p[2] - p[0]).abs()
}

/// The conic the tests draw: from (10, 5) to (70, 5) mm with its tangents
/// meeting at (30, 45) mm, its shoulder `height` mm above the chord; ρ is
/// `height` / 40 within the tool's range (clear of the origin, which fixes
/// a point clicked on it).
fn draw_conic(m: &mut Model, ed: &mut SketchEditor, height: f64) {
    let map = view();
    draw(
        m,
        ed,
        &map,
        Tool::Conic,
        &[(10.0, 5.0), (70.0, 5.0), (30.0, 45.0), (45.0, 5.0 + height)],
    );
}

/// The shoulder of the drawn conic for ρ: M + ρ(P₁ − M), in millimetres.
fn shoulder(rho: f64) -> (f64, f64) {
    (40.0 - 10.0 * rho, 5.0 + 40.0 * rho)
}

#[test]
fn four_clicks_make_a_conic_whose_shoulder_meets_the_pointer() {
    let (mut m, mut ed, map) = setup();
    let f = ed.feature;
    ed.set_tool(Tool::Conic);
    assert_eq!(ed.hint(), "Click the start point");
    click(&mut m, &mut ed, &map, 10.0, 5.0);
    assert_eq!(ed.hint(), "Click the end point");
    click(&mut m, &mut ed, &map, 70.0, 5.0);
    assert_eq!(ed.hint(), "Click where the tangents at the two ends meet");
    // Before the corner is chosen the band shows a parabola to the pointer.
    hover(&m, &mut ed, &map, 30.0, 45.0);
    let band = ed.rubber_band();
    assert_eq!(band.len(), 2);
    let parabola: Solved = (
        [
            DVec2::new(10.0, 5.0) * MM,
            DVec2::new(30.0, 45.0) * MM,
            DVec2::new(70.0, 5.0) * MM,
        ],
        1.0,
    );
    for p in &band[0] {
        assert!(implicit(&parabola, *p).0.abs() < 1e-5, "{p}");
    }
    click(&mut m, &mut ed, &map, 30.0, 45.0);
    assert!(ed.hint().starts_with("Click the shoulder"));
    // The curve, its control polygon and its median follow the pointer:
    // its highest point is level with it.
    hover(&m, &mut ed, &map, 45.0, 29.0);
    let band = ed.rubber_band();
    assert_eq!(band.len(), 3);
    let top = band[0].iter().map(|p| p.y).fold(f64::MIN, f64::max);
    assert!(near(top, 29.0 * MM, 1e-3), "{top}");
    click(&mut m, &mut ed, &map, 45.0, 29.0);
    assert!(ed.clicks.is_empty(), "the tool starts again");

    assert_eq!(count(&m, f, is_conic), 1);
    assert_eq!(count(&m, f, |g| matches!(g, Geometry::Point { .. })), 3);
    assert_eq!(constraints(&m, f, |_| true), 0);
    assert_eq!(status(&m, f), SolveStatus::Solved);
    assert_eq!(dof(&m, f), 7);
    let rho = rhos(&m, f)[0];
    assert!(close(rho, 0.6), "{rho}");
    let c = conics(&m, f).remove(0);
    assert!((c.0[1] - DVec2::new(30.0, 45.0) * MM).length() < 1e-8);
    assert!(close(c.1, 1.5), "w = ρ/(1 − ρ): {}", c.1);

    // Closed by its chord: one region of the segment's exact area, and an
    // extrusion of that times the depth.
    draw(
        &mut m,
        &mut ed,
        &map,
        Tool::Line,
        &[(70.0, 5.0), (10.0, 5.0)],
    );
    assert_eq!(status(&m, f), SolveStatus::Solved);
    let c = conics(&m, f).remove(0);
    let areas = region_areas(&m, f);
    assert_eq!(areas.len(), 1, "{areas:?}");
    assert!(near(areas[0], segment(&c), 1e-12), "{areas:?}");
    let x = extrude(&mut m, &map, &[(40.0, 15.0)], "8 mm");
    let v = volume(&m, x);
    assert!(near(v, segment(&c) * 8.0 * MM, 1e-9), "{v}");
}

#[test]
fn a_corner_in_line_with_the_ends_is_refused() {
    let (mut m, mut ed, map) = setup();
    let f = ed.feature;
    ed.set_tool(Tool::Conic);
    click(&mut m, &mut ed, &map, 10.0, 5.0);
    // The end on the start is not a second click.
    click(&mut m, &mut ed, &map, 10.0, 5.0);
    assert_eq!(ed.clicks.len(), 1);
    click(&mut m, &mut ed, &map, 70.0, 5.0);
    for (x, y) in [(40.0, 5.0), (100.0, 5.0), (-10.0, 5.0)] {
        click(&mut m, &mut ed, &map, x, y);
        let message = ed.message.clone().unwrap_or_default();
        assert!(message.contains("click off that line"), "{message}");
        assert_eq!(ed.clicks.len(), 2, "the clicks so far are kept");
    }
    // A corner beyond an end makes a lopsided but valid conic; a shoulder
    // click beyond the corner draws the tool's fullest.
    click(&mut m, &mut ed, &map, 90.0, 25.0);
    click(&mut m, &mut ed, &map, 50.0, 65.0);
    assert_eq!(count(&m, f, is_conic), 1);
    assert_eq!(status(&m, f), SolveStatus::Solved);
    assert_eq!(rhos(&m, f), vec![0.98]);
    // The drawn polyline lies on it.
    let (_, c) = curves(&m, f, false).remove(0);
    let solved = conics(&m, f).remove(0);
    for q in sample_resolved(&c) {
        let (g, inside) = implicit(&solved, q);
        assert!(g.abs() < 1e-12 && inside > -1e-12, "{q}: {g:e}");
    }
}

#[test]
fn the_rho_dimension_is_typed_as_a_plain_number() {
    let (mut m, mut ed, map) = setup();
    let f = ed.feature;
    draw_conic(&mut m, &mut ed, 24.0);
    let conic = conic_ids(&m, f)[0];
    // Selecting it at its shoulder and dimensioning it measures ρ.
    let (x, y) = shoulder(0.6);
    click(&mut m, &mut ed, &map, x, y);
    assert_eq!(ed.selection, vec![conic]);
    ed.dimension(&m, at(&map, x, y + 10.0));
    let prompt = ed.prompt.as_ref().expect("a prompt");
    assert!(matches!(prompt.constraint, Constraint::Rho { .. }));
    assert_eq!(prompt.text, "0.6");
    assert_eq!(prompt.unit(), "ρ");
    for bad in ["1.2", "0", "1", "-0.5", "abc", "0.4 mm"] {
        ed.prompt.as_mut().unwrap().text = bad.into();
        ed.commit_prompt(&mut m);
        let p = ed.prompt.as_ref().expect("the prompt stays open");
        assert!(
            p.error.as_deref().is_some_and(|e| e.contains("ρ")),
            "{bad}: {:?}",
            p.error
        );
    }
    let ops = m.session().log().len();
    ed.prompt.as_mut().unwrap().text = "0.3".into();
    ed.commit_prompt(&mut m);
    assert!(ed.prompt.is_none());
    assert_eq!(m.session().log().len(), ops + 1, "one undo step");
    assert_eq!(status(&m, f), SolveStatus::Solved);
    assert_eq!(dof(&m, f), 6);
    // The dimension's value exactly, in the document and the solution.
    assert_eq!(rhos(&m, f), vec![0.3]);
    let stored = &m.sketch_def(f).unwrap().sketch;
    assert!(matches!(
        stored.entity(conic).unwrap().geometry,
        Geometry::Conic { rho, .. } if rho == 0.3
    ));
    // Its label reads ρ beside the shoulder.
    let d = ed.display(&m).unwrap();
    let label = annotations(&d, &map)
        .into_iter()
        .find(|a| a.dimension)
        .expect("a label");
    assert_eq!(label.text, "ρ 0.3");
    let (x, y) = shoulder(0.3);
    assert!(label.rect.center().distance(at(&map, x, y)) < 60.0);
    // Its segment has the area of ρ = 0.3, w = 3/7.
    draw(
        &mut m,
        &mut ed,
        &map,
        Tool::Line,
        &[(70.0, 5.0), (10.0, 5.0)],
    );
    let areas = region_areas(&m, f);
    assert_eq!(areas.len(), 1);
    let c = conics(&m, f).remove(0);
    assert!(near(c.1, 3.0 / 7.0, 1e-15), "{}", c.1);
    assert!(near(areas[0], segment(&c), 1e-12), "{areas:?}");

    // Changing it: double-clicking the label opens the prompt on it.
    let dim = m.sketch_def(f).unwrap().sketch.constraints[0].id;
    ed.edit_dimension(&m, dim, at(&map, 0.0, 0.0));
    assert_eq!(ed.prompt.as_ref().unwrap().text, "0.3");
    ed.prompt.as_mut().unwrap().text = "0.85".into();
    ed.commit_prompt(&mut m);
    assert_eq!(rhos(&m, f), vec![0.85]);
    assert!(m.undo());
    assert_eq!(rhos(&m, f), vec![0.3]);
}

#[test]
fn dragging_the_curve_moves_its_shoulder_in_one_undo_step() {
    let (mut m, mut ed, map) = setup();
    let f = ed.feature;
    draw_conic(&mut m, &mut ed, 24.0);
    let before = conics(&m, f).remove(0);
    let rho0 = rhos(&m, f)[0];
    let ops = m.session().log().len();
    // Grabbed at its shoulder and pulled 8 mm further from the chord: ρ
    // grows by 8/40.
    let (x, y) = shoulder(0.6);
    ed.drag_start(&mut m, &map, at(&map, x, y));
    ed.drag_to(&m, &map, at(&map, x + 3.0, y + 4.0));
    ed.drag_to(&m, &map, at(&map, x + 5.0, y + 8.0));
    ed.drag_end(&mut m, &map, at(&map, x + 5.0, y + 8.0));
    assert_eq!(m.session().log().len(), ops + 1, "one undo step");
    assert_eq!(status(&m, f), SolveStatus::Solved);
    let rho = rhos(&m, f)[0];
    assert!((rho - (rho0 + 0.2)).abs() < 1e-6, "{rho}");
    // Its points stay put.
    let after = conics(&m, f).remove(0);
    for k in 0..3 {
        assert_eq!(after.0[k], before.0[k]);
    }
    assert!(m.undo());
    assert_eq!(rhos(&m, f), vec![rho0]);

    // Dragging its control point moves the point; a free ρ is kept.
    ed.drag_start(&mut m, &map, at(&map, 30.0, 45.0));
    ed.drag_to(&m, &map, at(&map, 40.0, 50.0));
    ed.drag_end(&mut m, &map, at(&map, 40.0, 50.0));
    let c = conics(&m, f).remove(0);
    assert!((c.0[1] - DVec2::new(40.0, 50.0) * MM).length() < 1e-6);
    assert_eq!(rhos(&m, f), vec![rho0]);

    // A dimensioned ρ does not follow the pointer.
    let conic = conic_ids(&m, f)[0];
    ed.selection = vec![conic];
    ed.dimension(&m, at(&map, 0.0, 0.0));
    ed.prompt.as_mut().unwrap().text = "0.5".into();
    ed.commit_prompt(&mut m);
    let ops = m.session().log().len();
    let c = conics(&m, f).remove(0);
    let mid = c.0[0].midpoint(c.0[2]);
    let s = (mid + (c.0[1] - mid) * 0.5) / MM;
    ed.drag_start(&mut m, &map, at(&map, s.x, s.y));
    ed.drag_to(&m, &map, at(&map, s.x, s.y + 6.0));
    ed.drag_end(&mut m, &map, at(&map, s.x, s.y + 6.0));
    assert_eq!(rhos(&m, f), vec![0.5]);
    assert!(m.session().log().len() <= ops + 1);
    assert_eq!(status(&m, f), SolveStatus::Solved);
}

#[test]
fn points_snap_onto_a_conic_exactly() {
    let (mut m, mut ed, map) = setup();
    let f = ed.feature;
    draw_conic(&mut m, &mut ed, 24.0);
    let conic = conic_ids(&m, f)[0];
    ed.set_tool(Tool::Point);
    // Near the curve between the drawn polyline's vertices.
    let c = conics(&m, f).remove(0);
    let q = openshape_sketch_solver::conic::point(c.0, c.1, 0.27) / MM;
    hover(&m, &mut ed, &map, q.x + 0.3, q.y - 0.2);
    let s = ed.cursor.expect("a snap");
    assert_eq!(s.target, SnapTarget::Curve(conic));
    let (g, inside) = implicit(&c, s.at);
    assert!(g.abs() < 1e-14 && inside >= 0.0, "{g:e}");
    click(&mut m, &mut ed, &map, q.x + 0.3, q.y - 0.2);
    assert_eq!(
        constraints(&m, f, |c| matches!(c, Constraint::PointOnCurve { .. })),
        1
    );
    assert_eq!(status(&m, f), SolveStatus::Solved);
    assert_eq!(dof(&m, f), 7 + 1);
    // The solved point is on the solved conic.
    let def = m.sketch_result(f).unwrap().solution.sketch.clone();
    let p = def
        .entities
        .iter()
        .rev()
        .find_map(|e| match e.geometry {
            Geometry::Point { at } => Some(at),
            _ => None,
        })
        .unwrap();
    let (g, inside) = implicit(&conics(&m, f).remove(0), p);
    assert!(g.abs() < 1e-12 && inside >= -1e-12, "{g:e}");
}

#[test]
fn curves_continue_tangentially_from_a_conics_ends() {
    // A line drawn on from the end along P₂ − P₁ is inferred tangent.
    let (mut m, mut ed, map) = setup();
    let f = ed.feature;
    draw_conic(&mut m, &mut ed, 24.0);
    // From (70, 5) along (40, −40): down to the right.
    draw(
        &mut m,
        &mut ed,
        &map,
        Tool::Line,
        &[(70.0, 5.0), (85.0, -10.2)],
    );
    assert_eq!(
        constraints(&m, f, |c| matches!(
            c,
            Constraint::Tangent {
                kind: TangentKind::AtPoint { .. },
                ..
            }
        )),
        1
    );
    assert_eq!(status(&m, f), SolveStatus::Solved);
    let c = conics(&m, f).remove(0);
    let (s, t) = curves(&m, f, false)
        .into_iter()
        .find_map(|(_, c)| match c {
            ResolvedCurve::Line { start, end } => Some((start, end)),
            _ => None,
        })
        .unwrap();
    assert!((s - c.0[2]).length() < 1e-12);
    let leg = (c.0[2] - c.0[1]).normalize();
    assert!((t - s).normalize().perp_dot(leg).abs() < 1e-12);

    // A tangent arc from its start: the arc's radius there is normal to
    // the first leg.
    draw(
        &mut m,
        &mut ed,
        &map,
        Tool::TangentArc,
        &[(10.0, 5.0), (-10.0, 15.0)],
    );
    assert_eq!(status(&m, f), SolveStatus::Solved);
    let c = conics(&m, f).remove(0);
    let (center, radius) = curves(&m, f, false)
        .into_iter()
        .find_map(|(_, c)| match c {
            ResolvedCurve::Arc { center, radius, .. } => Some((center, radius)),
            _ => None,
        })
        .unwrap();
    let p = c.0[0];
    assert!((((p - center).length() - radius) / radius).abs() < 1e-12);
    assert!(
        (p - center)
            .normalize()
            .dot((c.0[1] - c.0[0]).normalize())
            .abs()
            < 1e-9
    );

    // Tangency with a curve it does not share an end with is refused.
    draw(
        &mut m,
        &mut ed,
        &map,
        Tool::Circle,
        &[(40.0, 65.0), (40.0, 75.0)],
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
    let conic = conic_ids(&m, f)[0];
    ed.selection = vec![conic, circle];
    ed.constrain(&mut m, ConstraintTool::Tangent);
    let message = ed.message.clone().unwrap_or_default();
    assert!(
        message.contains("conic needs a shared end point"),
        "{message}"
    );
    for tool in [
        ConstraintTool::Equal,
        ConstraintTool::Concentric,
        ConstraintTool::Normal,
    ] {
        ed.selection = vec![conic, circle];
        ed.constrain(&mut m, tool);
        let message = ed.message.clone().unwrap_or_default();
        assert!(message.contains("needs"), "{tool:?}: {message}");
    }
}

#[test]
fn mirrored_and_patterned_conics_keep_their_rho() {
    let (mut m, mut ed, map) = setup();
    let f = ed.feature;
    draw_conic(&mut m, &mut ed, 24.0);
    ed.construction = true;
    draw(
        &mut m,
        &mut ed,
        &map,
        Tool::Line,
        &[(-10.0, -50.0), (-10.0, 50.0)],
    );
    ed.construction = false;
    let before = dof(&m, f);
    ed.set_tool(Tool::Mirror);
    click(&mut m, &mut ed, &map, -10.0, -40.0);
    let (x, y) = shoulder(0.6);
    click(&mut m, &mut ed, &map, x, y);
    assert_eq!(ed.message.as_deref(), Some("Mirrored 1 item"));
    assert_eq!(status(&m, f), SolveStatus::Solved);
    assert_eq!(dof(&m, f), before, "the image adds no freedom");
    let cs = conics(&m, f);
    assert_eq!(cs.len(), 2);
    let reflect = |p: DVec2| DVec2::new(-20.0 * MM - p.x, p.y);
    for k in 0..3 {
        assert!((cs[1].0[k] - reflect(cs[0].0[k])).length() < 1e-12);
    }
    let r = rhos(&m, f);
    assert_eq!(r[0], r[1]);
    // A ρ dimension on the source drives the image exactly.
    ed.set_tool(Tool::Select);
    ed.selection = vec![conic_ids(&m, f)[0]];
    ed.dimension(&m, at(&map, 0.0, 0.0));
    ed.prompt.as_mut().unwrap().text = "0.7".into();
    ed.commit_prompt(&mut m);
    assert_eq!(status(&m, f), SolveStatus::Solved);
    assert_eq!(rhos(&m, f), vec![0.7, 0.7]);

    // Three copies along a step, each with the source's ρ and area.
    let (mut m, mut ed, map) = setup();
    let f = ed.feature;
    draw_conic(&mut m, &mut ed, 24.0);
    draw(
        &mut m,
        &mut ed,
        &map,
        Tool::Line,
        &[(70.0, 5.0), (10.0, 5.0)],
    );
    ed.selection = m
        .sketch_def(f)
        .unwrap()
        .sketch
        .entities
        .iter()
        .filter(|e| e.geometry.is_curve())
        .map(|e| e.id)
        .collect();
    ed.set_tool(Tool::LinearPattern);
    click(&mut m, &mut ed, &map, 0.0, -10.0);
    click(&mut m, &mut ed, &map, 0.0, -60.0);
    assert_eq!(
        ed.message.as_deref(),
        Some("Patterned 2 items into 3 instances")
    );
    ed.prompt = None;
    assert_eq!(status(&m, f), SolveStatus::Solved);
    let cs = conics(&m, f);
    assert_eq!(cs.len(), 3);
    let r = rhos(&m, f);
    assert!(r.iter().all(|x| *x == r[0]), "{r:?}");
    let areas = region_areas(&m, f);
    assert_eq!(areas.len(), 3, "{areas:?}");
    for a in areas {
        assert!(near(a, segment(&cs[0]), 1e-12), "{a}");
    }
}

#[test]
fn editing_tools_that_cannot_handle_conics_say_so() {
    let (mut m, mut ed, map) = setup();
    let f = ed.feature;
    draw_conic(&mut m, &mut ed, 24.0);
    let before = m.sketch_def(f).unwrap().sketch.clone();
    // Offset, on the curve.
    let c = conics(&m, f).remove(0);
    let q = conic_point(c.0, c.1, 0.3) / MM;
    ed.set_tool(Tool::Offset);
    click(&mut m, &mut ed, &map, q.x, q.y);
    let message = ed.message.clone().unwrap_or_default();
    assert!(message.contains("its offset is not a conic"), "{message}");
    assert_eq!(m.sketch_def(f).unwrap().sketch, before);

    // A fillet where it meets a line.
    let (mut m, mut ed, map) = setup();
    let f = ed.feature;
    draw_conic(&mut m, &mut ed, 24.0);
    draw(
        &mut m,
        &mut ed,
        &map,
        Tool::Line,
        &[(70.0, 5.0), (70.0, -15.0)],
    );
    let before = m.sketch_def(f).unwrap().sketch.clone();
    ed.set_tool(Tool::SketchFillet);
    click(&mut m, &mut ed, &map, 70.0, 5.0);
    let message = ed.message.clone().unwrap_or_default();
    assert!(message.contains("a conic ends here"), "{message}");
    assert_eq!(m.sketch_def(f).unwrap().sketch, before);
}

#[test]
fn drawing_is_one_undo_step_and_deleting_takes_the_points() {
    let (mut m, mut ed, _) = setup();
    let f = ed.feature;
    let ops = m.session().log().len();
    draw_conic(&mut m, &mut ed, 24.0);
    assert_eq!(m.session().log().len(), ops + 1);
    ed.selection = conic_ids(&m, f);
    ed.delete_selection(&mut m);
    assert_eq!(m.sketch_def(f).unwrap().sketch.entities.len(), 0);
    assert!(m.undo());
    assert_eq!(count(&m, f, is_conic), 1);
    assert!(m.undo());
    assert_eq!(m.sketch_def(f).unwrap().sketch.entities.len(), 0);
}

/// The conic entity's (start, control, end) points.
fn conic_points(m: &Model, f: FeatureId, id: SketchId) -> [SketchId; 3] {
    match m.sketch_def(f).unwrap().sketch.entity(id).unwrap().geometry {
        Geometry::Conic {
            start,
            control,
            end,
            ..
        } => [start, control, end],
        ref g => panic!("{g:?}"),
    }
}

/// The solved position of point `id`.
fn solved_point(m: &Model, f: FeatureId, id: SketchId) -> DVec2 {
    m.sketch_result(f)
        .unwrap()
        .solution
        .sketch
        .point(id)
        .unwrap()
}

fn redundant(m: &Model, f: FeatureId) -> usize {
    m.sketch_result(f).unwrap().solution.redundant.len()
}

/// The ids of the sketch's lines, in order.
fn line_ids(m: &Model, f: FeatureId) -> Vec<SketchId> {
    m.sketch_def(f)
        .unwrap()
        .sketch
        .entities
        .iter()
        .filter(|e| matches!(e.geometry, Geometry::Line { .. }))
        .map(|e| e.id)
        .collect()
}

/// Every point of conic `c` (33 samples) is on the conic of `of`: Lee's
/// implicit equation of `of`'s control triangle, not `c`'s.
fn on_conic(c: &Solved, of: &Solved) -> bool {
    (0..=32).all(|i| {
        let x = conic_point(c.0, c.1, i as f64 / 32.0);
        implicit(of, x).0.abs() < 1e-9
    })
}

/// `a`, `b` and `c` lie on one line (relative to their spread).
fn in_line(a: DVec2, b: DVec2, c: DVec2) -> bool {
    (b - a).perp_dot(c - a).abs() <= 1e-9 * (b - a).length() * (c - a).length()
}

/// The first parameter of conic `c`, stepping from `t0` by `dt`, where
/// `f` of its point changes sign, refined by bisection.
fn first_crossing(c: &Solved, t0: f64, dt: f64, f: impl Fn(DVec2) -> f64) -> f64 {
    let g = |t: f64| f(conic_point(c.0, c.1, t));
    let (mut a, s0) = (t0, g(t0).signum());
    loop {
        let b = a + dt;
        assert!(b.abs() < 10.0, "no crossing");
        if g(b).signum() != s0 {
            let (mut lo, mut hi) = (a, b);
            for _ in 0..200 {
                let mid = 0.5 * (lo + hi);
                if g(mid).signum() == s0 {
                    lo = mid;
                } else {
                    hi = mid;
                }
            }
            return 0.5 * (lo + hi);
        }
        a = b;
    }
}

/// Uses `tool` once at a plane point (millimetres) and returns its message.
fn act(m: &mut Model, ed: &mut SketchEditor, tool: Tool, at_mm: DVec2) -> String {
    let map = view();
    ed.set_tool(tool);
    click(m, ed, &map, at_mm.x, at_mm.y);
    ed.set_tool(Tool::Select);
    ed.message.clone().unwrap_or_default()
}

/// Closes conic `c` with a line along its chord and checks the region
/// against the closed form of its segment.
fn check_segment(m: &mut Model, ed: &mut SketchEditor, c: &Solved) {
    let map = view();
    let f = ed.feature;
    let (a, b) = (c.0[0] / MM, c.0[2] / MM);
    draw(m, ed, &map, Tool::Line, &[(a.x, a.y), (b.x, b.y)]);
    assert_eq!(status(m, f), SolveStatus::Solved);
    let areas = region_areas(m, f);
    assert_eq!(areas.len(), 1, "{areas:?}");
    let c = conics(m, f)
        .into_iter()
        .find(|k| (k.0[0] - c.0[0]).length() < 1e-9)
        .unwrap();
    assert!(
        near(areas[0], segment(&c), 1e-11),
        "{areas:?} vs {}",
        segment(&c)
    );
}

#[test]
fn trimming_a_conic_keeps_an_exact_piece_of_it() {
    // Elliptic, parabolic and hyperbolic; the piece left of the line, and
    // the one right of it.
    for height in [16.0, 20.0, 24.0] {
        for left in [true, false] {
            let (mut m, mut ed, map) = setup();
            let f = ed.feature;
            draw_conic(&mut m, &mut ed, height);
            draw(
                &mut m,
                &mut ed,
                &map,
                Tool::Line,
                &[(34.0, -5.0), (34.0, 55.0)],
            );
            let line = line_ids(&m, f)[0];
            let conic = conic_ids(&m, f)[0];
            let [start, control, end] = conic_points(&m, f, conic);
            let old = conics(&m, f).remove(0);
            let before = dof(&m, f);
            let line_x = curves(&m, f, false)
                .into_iter()
                .find_map(|(id, c)| match c {
                    ResolvedCurve::Line { start, .. } if id == line => Some(start.x),
                    _ => None,
                })
                .unwrap();
            let tc = first_crossing(&old, 0.0, 1e-3, |p| p.x - line_x);
            let (moved, kept) = if left { (start, end) } else { (end, start) };
            let kept_at = solved_point(&m, f, kept);

            // The preview is the piece that goes, on the conic.
            let q = conic_point(old.0, old.1, if left { 0.1 } else { 0.9 }) / MM;
            ed.set_tool(Tool::Trim);
            hover(&m, &mut ed, &map, q.x, q.y);
            match ed.curve_edit.clone().unwrap().plan.unwrap().preview {
                Preview::Piece(ResolvedCurve::Conic {
                    start,
                    control,
                    end,
                    weight,
                }) => {
                    let gone = ([start, control, end], weight);
                    assert!(on_conic(&gone, &old));
                    let far = if left { old.0[0] } else { old.0[2] };
                    assert!((start - far).length().min((end - far).length()) < 1e-12);
                }
                p => panic!("{p:?}"),
            }
            let ops = m.session().log().len();
            act(&mut m, &mut ed, Tool::Trim, q);
            assert_eq!(m.session().log().len(), ops + 1, "one undo step");

            assert_eq!(conic_ids(&m, f), vec![conic], "the same conic, cut");
            assert_eq!(conic_points(&m, f, conic), [start, control, end]);
            let c = conics(&m, f).remove(0);
            let cut = conic_point(old.0, old.1, tc);
            assert!((solved_point(&m, f, moved) - cut).length() < 1e-9);
            assert!((solved_point(&m, f, kept) - kept_at).length() < 1e-12);
            let on = Constraint::PointOnCurve {
                point: moved,
                curve: line,
            };
            assert_eq!(constraints(&m, f, |k| *k == on), 1);
            // The kept end keeps its tangent: the new control point is on
            // the old one's leg from it.
            assert!(in_line(kept_at, c.0[1], old.0[1]), "{:?}", c.0);
            assert!(on_conic(&c, &old), "the piece is part of the conic");
            assert!(near(rhos(&m, f)[0], c.1 / (1.0 + c.1), 1e-12));
            assert_eq!(status(&m, f), SolveStatus::Solved);
            assert_eq!(redundant(&m, f), 0);
            assert_eq!(dof(&m, f), before - 1, "the cut end is held on the line");
            check_segment(&mut m, &mut ed, &c);
        }
    }
}

/// Adds `constraint` to the sketch as one edit.
fn constrain(m: &mut Model, f: FeatureId, constraint: Constraint<SketchId>) {
    m.edit_sketch(
        f,
        vec![SketchEdit::AddConstraint(ConstraintEntry {
            id: SketchId::new_in_batch(0),
            constraint,
        })],
    )
    .unwrap();
}

#[test]
fn trimming_the_middle_of_a_conic_leaves_two_conics() {
    let (mut m, mut ed, map) = setup();
    let f = ed.feature;
    draw_conic(&mut m, &mut ed, 16.0);
    let conic = conic_ids(&m, f)[0];
    let [start, control, end] = conic_points(&m, f, conic);
    constrain(
        &mut m,
        f,
        Constraint::Rho {
            curve: conic,
            value: 0.4,
        },
    );
    // A line on from its end, tangent there.
    draw(
        &mut m,
        &mut ed,
        &map,
        Tool::Line,
        &[(70.0, 5.0), (85.0, -10.0)],
    );
    let tangent_line = line_ids(&m, f)[0];
    for x in [25.0, 55.0] {
        draw(&mut m, &mut ed, &map, Tool::Line, &[(x, -5.0), (x, 55.0)]);
    }
    let old = conics(&m, f).remove(0);
    // Points held on the conic near its start, in the middle and near its
    // end.
    let mut b = Batch::default();
    for t in [0.05, 0.55, 0.95] {
        let p = b.point(conic_point(old.0, old.1, t));
        b.constraint(Constraint::PointOnCurve {
            point: p,
            curve: conic,
        });
    }
    m.edit_sketch(f, b.edits).unwrap();
    // Their ids in the document, along the conic.
    let mut held: Vec<SketchId> = m
        .sketch_def(f)
        .unwrap()
        .sketch
        .constraints
        .iter()
        .filter_map(|c| match c.constraint {
            Constraint::PointOnCurve { point, curve } if curve == conic => Some(point),
            _ => None,
        })
        .collect();
    held.sort_by_key(|p| (solved_point(&m, f, *p).x * 1e6) as i64);
    assert_eq!(status(&m, f), SolveStatus::Solved);
    assert_eq!(redundant(&m, f), 0);
    let before = dof(&m, f);

    let message = act(
        &mut m,
        &mut ed,
        Tool::Trim,
        conic_point(old.0, old.1, 0.45) / MM,
    );
    assert!(message.contains("removed 2 constraints"), "{message}");
    let ids = conic_ids(&m, f);
    assert_eq!(ids.len(), 2, "{message}");
    assert_eq!(ids[0], conic);
    let first = conic_points(&m, f, ids[0]);
    let second = conic_points(&m, f, ids[1]);
    assert_eq!(first, [start, control, end]);
    assert!((solved_point(&m, f, second[2]) - old.0[2]).length() < 1e-12);
    let cs = conics(&m, f);
    for c in &cs {
        assert!(on_conic(c, &old));
    }
    // Each keeps its outer end's tangent.
    assert!(in_line(old.0[0], cs[0].0[1], old.0[1]));
    assert!(in_line(old.0[2], cs[1].0[1], old.0[1]));
    // The ρ dimension went with the shape it measured; the middle point's
    // constraint went with the piece it was on; the end's tangency and the
    // point near the end moved to the new conic.
    assert_eq!(
        constraints(&m, f, |k| matches!(k, Constraint::Rho { .. })),
        0
    );
    let on = |p: SketchId, c: SketchId| Constraint::PointOnCurve { point: p, curve: c };
    assert_eq!(constraints(&m, f, |k| *k == on(held[0], conic)), 1);
    assert_eq!(
        constraints(&m, f, |k| matches!(k,
        Constraint::PointOnCurve { point, .. } if *point == held[1])),
        0
    );
    assert_eq!(constraints(&m, f, |k| *k == on(held[2], ids[1])), 1);
    assert_eq!(
        constraints(&m, f, |k| matches!(k,
            Constraint::Tangent { a, b, kind: TangentKind::AtPoint { .. } }
                if [*a, *b].contains(&tangent_line) && [*a, *b].contains(&ids[1]))),
        1
    );
    let line_start = match m
        .sketch_def(f)
        .unwrap()
        .sketch
        .entity(tangent_line)
        .unwrap()
        .geometry
    {
        Geometry::Line { start, .. } => start,
        _ => unreachable!(),
    };
    assert!((solved_point(&m, f, line_start) - solved_point(&m, f, second[2])).length() < 1e-12);
    assert_eq!(status(&m, f), SolveStatus::Solved);
    assert_eq!(redundant(&m, f), 0);
    // Per conic: 7 less its cut end held on a line; the dimension's 1
    // back, the middle point's 1 back. Before, one conic of 6.
    assert_eq!(dof(&m, f), before - 6 + 2 * 6 + 1);
    // The trim preview of the first conic's piece is gone, and the pieces
    // extrude as exact segments.
    check_segment(&mut m, &mut ed, &cs[1]);
}

#[test]
fn splitting_a_conic_leaves_two_conics_tangent_at_the_split() {
    for height in [16.0, 24.0] {
        let (mut m, mut ed, map) = setup();
        let f = ed.feature;
        draw_conic(&mut m, &mut ed, height);
        let conic = conic_ids(&m, f)[0];
        let [start, control, end] = conic_points(&m, f, conic);
        let old = conics(&m, f).remove(0);
        let before = dof(&m, f);
        let ops = m.session().log().len();
        let at_split = conic_point(old.0, old.1, 0.4);
        ed.set_tool(Tool::Split);
        hover(&m, &mut ed, &map, at_split.x / MM, at_split.y / MM);
        match ed.curve_edit.clone().unwrap().plan.unwrap().preview {
            Preview::Point(q) => assert!((q - at_split).length() < 1e-7),
            p => panic!("{p:?}"),
        }
        act(&mut m, &mut ed, Tool::Split, at_split / MM);
        assert_eq!(m.session().log().len(), ops + 1);
        let ids = conic_ids(&m, f);
        assert_eq!(ids.len(), 2, "{:?}", ed.message);
        assert_eq!(conic_points(&m, f, ids[0]), [start, control, end]);
        let second = conic_points(&m, f, ids[1]);
        let joint = Constraint::Coincident {
            a: end,
            b: second[0],
        };
        assert_eq!(constraints(&m, f, |k| *k == joint), 1);
        let smooth = Constraint::Tangent {
            a: ids[0],
            b: ids[1],
            kind: TangentKind::AtPoint { point: end },
        };
        assert_eq!(constraints(&m, f, |k| *k == smooth), 1);
        let cs = conics(&m, f);
        for c in &cs {
            assert!(on_conic(c, &old));
        }
        assert!((cs[0].0[2] - at_split).length() < 1e-7);
        assert_eq!(status(&m, f), SolveStatus::Solved);
        assert_eq!(redundant(&m, f), 0);
        // The second conic's control point and ρ, and its end, less the
        // tangency.
        assert_eq!(dof(&m, f), before + 4);
        // Together they bound the old segment.
        draw(
            &mut m,
            &mut ed,
            &map,
            Tool::Line,
            &[(10.0, 5.0), (70.0, 5.0)],
        );
        let areas = region_areas(&m, f);
        assert_eq!(areas.len(), 1, "{areas:?}");
        assert!(near(areas[0], segment(&old), 1e-11), "{areas:?}");
        // Dragging the joint keeps them tangent there.
        let j = solved_point(&m, f, end) / MM;
        ed.drag_start(&mut m, &map, at(&map, j.x, j.y));
        ed.drag_to(&m, &map, at(&map, j.x + 4.0, j.y + 3.0));
        ed.drag_end(&mut m, &map, at(&map, j.x + 4.0, j.y + 3.0));
        assert_eq!(status(&m, f), SolveStatus::Solved);
        let cs = conics(&m, f);
        assert!((cs[0].0[2] / MM - j).length() > 1.0, "it moved");
        let (u, v) = (cs[0].0[2] - cs[0].0[1], cs[1].0[1] - cs[1].0[0]);
        assert!(u.normalize().perp_dot(v.normalize()).abs() < 1e-9);
        assert!(u.dot(v) > 0.0);
    }
}

#[test]
fn extending_a_conic_follows_the_same_conic() {
    for height in [16.0, 24.0] {
        let (mut m, mut ed, map) = setup();
        let f = ed.feature;
        draw_conic(&mut m, &mut ed, height);
        let conic = conic_ids(&m, f)[0];
        let [start, control, end] = conic_points(&m, f, conic);
        // Below its end, and below its start.
        draw(
            &mut m,
            &mut ed,
            &map,
            Tool::Line,
            &[(40.0, -3.0), (110.0, -3.0)],
        );
        draw(
            &mut m,
            &mut ed,
            &map,
            Tool::Line,
            &[(-20.0, -2.0), (12.0, -2.0)],
        );
        let [below, left] = [line_ids(&m, f)[0], line_ids(&m, f)[1]];
        // `left` is the line below the start.
        let old = conics(&m, f).remove(0);
        let y = curves(&m, f, false)
            .into_iter()
            .find_map(|(id, c)| match c {
                ResolvedCurve::Line { start, .. } if id == below => Some(start.y),
                _ => None,
            })
            .unwrap();
        let y_left = curves(&m, f, false)
            .into_iter()
            .find_map(|(id, c)| match c {
                ResolvedCurve::Line { start, .. } if id == left => Some(start.y),
                _ => None,
            })
            .unwrap();
        let t_end = first_crossing(&old, 1.0, 1e-3, |p| p.y - y);
        let t_start = first_crossing(&old, 0.0, -1e-3, |p| p.y - y_left);
        let before = dof(&m, f);

        // At the end: on along the conic, down to the line below.
        let near_end = conic_point(old.0, old.1, 0.95) / MM;
        ed.set_tool(Tool::Extend);
        hover(&m, &mut ed, &map, near_end.x, near_end.y);
        match ed.curve_edit.clone().unwrap().plan.unwrap().preview {
            Preview::Piece(ResolvedCurve::Conic {
                start,
                control,
                end,
                weight,
            }) => {
                let added = ([start, control, end], weight);
                assert!(on_conic(&added, &old));
                assert!((start - old.0[2]).length() < 1e-12);
                assert!((end - conic_point(old.0, old.1, t_end)).length() < 1e-9);
            }
            p => panic!("{p:?}"),
        }
        act(&mut m, &mut ed, Tool::Extend, near_end);
        assert_eq!(conic_points(&m, f, conic), [start, control, end]);
        let c = conics(&m, f).remove(0);
        assert!((c.0[2] - conic_point(old.0, old.1, t_end)).length() < 1e-9);
        assert!((c.0[0] - old.0[0]).length() < 1e-12);
        assert!(on_conic(&c, &old));
        assert!(
            in_line(old.0[0], c.0[1], old.0[1]),
            "the start keeps its tangent"
        );
        let on = |p: SketchId, l: SketchId| Constraint::PointOnCurve { point: p, curve: l };
        assert_eq!(constraints(&m, f, |k| *k == on(end, below)), 1);
        assert_eq!(dof(&m, f), before - 1);

        // At the start: back along it to the line below it.
        let near_start = conic_point(c.0, c.1, 0.02) / MM;
        act(&mut m, &mut ed, Tool::Extend, near_start);
        let c = conics(&m, f).remove(0);
        assert!((c.0[0] - conic_point(old.0, old.1, t_start)).length() < 1e-9);
        assert!(on_conic(&c, &old));
        assert_eq!(constraints(&m, f, |k| *k == on(start, left)), 1);
        assert_eq!(status(&m, f), SolveStatus::Solved);
        assert_eq!(redundant(&m, f), 0);
        assert_eq!(dof(&m, f), before - 2);
        check_segment(&mut m, &mut ed, &c);
    }

    // Nothing beyond either end: refused.
    let (mut m, mut ed, _) = setup();
    draw_conic(&mut m, &mut ed, 16.0);
    let c = conics(&m, ed.feature).remove(0);
    let message = act(
        &mut m,
        &mut ed,
        Tool::Extend,
        conic_point(c.0, c.1, 0.9) / MM,
    );
    assert!(message.contains("Nothing to extend to"), "{message}");
}

#[test]
fn a_piece_too_short_to_keep_its_shape_is_refused() {
    // A long, flat conic (100 mm across, its corner 1 mm up; radius of
    // curvature 2.5 m at its start) crossed half a millimetre from its
    // start, where the piece's corner lies t² mm ≈ 0.025 µm off its chord.
    let mut m = Model::new();
    let f = m.add_sketch(openshape_core_model::PlaneRef::Top);
    let mut b = Batch::default();
    let p = [DVec2::ZERO, DVec2::new(50.0, 1.0), DVec2::new(100.0, 0.0)].map(|x| x * MM);
    let [start, control, end] = p.map(|x| b.point(x));
    b.entity(
        Geometry::Conic {
            start,
            control,
            end,
            rho: 0.5,
        },
        false,
    );
    b.line(DVec2::new(0.5, -5.0) * MM, DVec2::new(0.5, 5.0) * MM, false);
    m.edit_sketch(f, b.edits).unwrap();
    let conic = conic_ids(&m, f)[0];
    let r = m.sketch_result(f).unwrap();
    let curves: Vec<(SketchId, ResolvedCurve)> =
        r.curves.iter().map(|(id, _, c)| (*id, c.clone())).collect();
    let sketch = r.solution.sketch.clone();
    let cuts = trim::cuts(m.kernel(), &sketch, &curves, conic).unwrap();
    let short = "too short to keep its shape";

    // Keeping the short piece, or splitting one off, is refused.
    let message = trim::trim(&sketch, &cuts, conic_point(p, 1.0, 0.8)).unwrap_err();
    assert!(message.contains(short), "{message}");
    let message = trim::split(&sketch, &curves, conic, conic_point(p, 1.0, 0.004)).unwrap_err();
    assert!(message.contains(short), "{message}");

    // Trimming it away keeps the long piece.
    let plan = trim::trim(&sketch, &cuts, conic_point(p, 1.0, 0.002)).unwrap();
    m.edit_sketch_checked(f, plan.edits).unwrap();
    assert_eq!(status(&m, f), SolveStatus::Solved);
    let c = conics(&m, f).remove(0);
    assert!((c.0[0].x - 0.5 * MM).abs() < 1e-12, "{:?}", c.0);
    assert!(on_conic(&c, &(p, 1.0)));
}
