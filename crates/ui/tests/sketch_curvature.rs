//! The curvature (G2) constraint tool (design note 0010 §5), driven by
//! clicks through a top-view camera: joining a spline to an arc, the
//! refusals, replacing a tangency, and trim and split handing the
//! constraint on. Curvatures are checked from the solved control points by
//! the end-derivative formulas of Piegl & Tiller (eq. 3.8 applied twice),
//! independently of the solver's own evaluation.

mod common;

use common::*;
use opendatum_app::sketch_mode::{ConstraintTool, SketchEditor, Tool};
use opendatum_app::Model;
use opendatum_core_model::sketch::{Constraint, Geometry, TangentKind};
use opendatum_core_model::{DVec2, FeatureId, SketchId};
use opendatum_sketch_solver::{ResolvedCurve, SolveStatus};

/// The ids of the sketch's entities matching `pred`, in order.
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

fn spline_id(m: &Model, f: FeatureId) -> SketchId {
    ids(m, f, |g| matches!(g, Geometry::ControlSpline { .. }))[0]
}

fn arc_id(m: &Model, f: FeatureId) -> SketchId {
    ids(m, f, |g| matches!(g, Geometry::Arc { .. }))[0]
}

/// The solved curve `id`.
fn resolved(m: &Model, f: FeatureId, id: SketchId) -> ResolvedCurve {
    curves(m, f, true)
        .into_iter()
        .find(|(k, _)| *k == id)
        .map(|(_, c)| c)
        .expect("a solved curve")
}

/// A spline's start point, unit tangent and signed curvature there (for
/// travel away from the start): `C′(0) = Q₀`, `C″(0) = (p − 1)(Q₁ −
/// Q₀)/(u_{p+1} − u₂)` with `Qᵢ = p(Pᵢ₊₁ − Pᵢ)/(u_{i+p+1} − u_{i+1})`.
fn spline_start(c: &ResolvedCurve) -> (DVec2, DVec2, f64) {
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
    (p[0], c1.normalize(), c1.perp_dot(c2) / c1.length().powi(3))
}

/// An arc's end point, its counter-clockwise unit tangent there, its
/// curvature (1/r) and radius.
fn arc_end(c: &ResolvedCurve) -> (DVec2, DVec2, f64, f64) {
    let ResolvedCurve::Arc {
        center,
        radius,
        start_angle,
        sweep_angle,
    } = c
    else {
        panic!("not an arc: {c:?}");
    };
    let a = start_angle + sweep_angle;
    let radial = DVec2::new(a.cos(), a.sin());
    (
        *center + radial * *radius,
        radial.perp(),
        1.0 / radius,
        *radius,
    )
}

fn redundant(m: &Model, f: FeatureId) -> usize {
    m.sketch_result(f).unwrap().solution.redundant.len()
}

fn curvature_constraints(m: &Model, f: FeatureId) -> Vec<(SketchId, SketchId, SketchId)> {
    m.sketch_def(f)
        .unwrap()
        .sketch
        .constraints
        .iter()
        .filter_map(|c| match c.constraint {
            Constraint::Curvature { a, b, point } => Some((a, b, point)),
            _ => None,
        })
        .collect()
}

/// A counter-clockwise arc about (0, 20) from (−20, 20) down to (0, 0),
/// turning left, and a cubic spline leaving (0, 0) up and to the right and
/// turning right.
fn arc_and_spline() -> (Model, SketchEditor, FeatureId) {
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
    (m, ed, f)
}

#[test]
fn the_curvature_tool_joins_a_spline_to_an_arc_smoothly() {
    let (mut m, mut ed, f) = arc_and_spline();
    let (spline, arc) = (spline_id(&m, f), arc_id(&m, f));
    let (_, t0, k0) = spline_start(&resolved(&m, f, spline));
    let (_, ta0, ka0, r0) = arc_end(&resolved(&m, f, arc));
    // Not smooth to start with: the spline bends the other way.
    assert!(k0 * ka0 < 0.0, "{k0} {ka0}");
    assert!(t0.perp_dot(ta0).abs() > 1e-3);
    let before = dof(&m, f);

    ed.selection = vec![spline, arc];
    assert!(ed.can_constrain(&m, ConstraintTool::Curvature));
    ed.constrain(&mut m, ConstraintTool::Curvature);
    assert_eq!(curvature_constraints(&m, f).len(), 1, "{:?}", ed.message);
    assert_eq!(status(&m, f), SolveStatus::Solved);
    assert_eq!(redundant(&m, f), 0);
    assert_eq!(dof(&m, f), before - 2);

    let (p, t, k) = spline_start(&resolved(&m, f, spline));
    let (q, ta, ka, r) = arc_end(&resolved(&m, f, arc));
    assert!((p - q).length() < 1e-12, "{p} {q}");
    // One tangent, continuing the arc's direction of travel, and the
    // arc's curvature carried on: G2.
    assert!(t.dot(ta) > 1.0 - 1e-12, "{t} {ta}");
    assert!((k - ka).abs() <= 1e-9 * ka, "{k} vs {ka}");
    // The spline turned its curvature round; the arc kept its shape.
    assert!(k > 0.0);
    assert!((r / r0 - 1.0).abs() < 0.01, "radius {r0} became {r}");

    // One undo step takes it all back.
    assert!(m.undo());
    assert!(curvature_constraints(&m, f).is_empty());
    let (_, _, k_undone) = spline_start(&resolved(&m, f, spline));
    assert_eq!(k_undone, k0);
}

#[test]
fn the_curvature_tool_explains_what_it_cannot_join() {
    let (mut m, mut ed, f) = arc_and_spline();
    let map = view();
    let (spline, arc) = (spline_id(&m, f), arc_id(&m, f));
    let refused = |m: &mut Model, ed: &mut SketchEditor, sel: Vec<SketchId>, needle: &str| {
        ed.selection = sel;
        assert!(!ed.can_constrain(m, ConstraintTool::Curvature));
        let n = curvature_constraints(m, f).len();
        ed.constrain(m, ConstraintTool::Curvature);
        let message = ed.message.clone().unwrap_or_default();
        assert!(message.contains(needle), "{needle}: {message}");
        assert_eq!(curvature_constraints(m, f).len(), n);
    };
    // A line meeting the arc: no spline.
    draw(
        &mut m,
        &mut ed,
        &map,
        Tool::Line,
        &[(-20.0, 20.0), (-20.0, 40.0)],
    );
    let line = ids(&m, f, |g| matches!(g, Geometry::Line { .. }))[0];
    refused(&mut m, &mut ed, vec![line, arc], "needs a spline");
    // A circle has no end.
    draw(
        &mut m,
        &mut ed,
        &map,
        Tool::Circle,
        &[(60.0, 0.0), (70.0, 0.0)],
    );
    let circle = ids(&m, f, |g| matches!(g, Geometry::Circle { .. }))[0];
    refused(&mut m, &mut ed, vec![spline, circle], "no end to join");
    // Curves that do not meet.
    refused(&mut m, &mut ed, vec![spline, line], "share an end point");
    // A point is not a curve.
    let ends = ids(&m, f, |g| matches!(g, Geometry::Point { .. }));
    ed.selection = vec![spline, ends[0]];
    assert!(!ed.can_constrain(&m, ConstraintTool::Curvature));
}

#[test]
fn curvature_continuity_replaces_the_tangency_it_implies() {
    let (mut m, mut ed, f) = arc_and_spline();
    let (spline, arc) = (spline_id(&m, f), arc_id(&m, f));
    ed.selection = vec![spline, arc];
    ed.constrain(&mut m, ConstraintTool::Tangent);
    let tangencies = |m: &Model| {
        constraints(m, f, |c| {
            matches!(
                c,
                Constraint::Tangent {
                    kind: TangentKind::AtPoint { .. },
                    ..
                }
            )
        })
    };
    assert_eq!(tangencies(&m), 1);
    let before = dof(&m, f);

    ed.selection = vec![arc, spline];
    ed.constrain(&mut m, ConstraintTool::Curvature);
    let message = ed.message.clone().unwrap_or_default();
    assert!(
        message.starts_with("Replaced the tangency there"),
        "{message}"
    );
    assert_eq!(tangencies(&m), 0);
    assert_eq!(curvature_constraints(&m, f).len(), 1);
    assert_eq!(status(&m, f), SolveStatus::Solved);
    assert_eq!(redundant(&m, f), 0, "nothing left implied");
    assert_eq!(dof(&m, f), before - 1, "one freedom more than tangency");
}

#[test]
fn trim_and_split_hand_curvature_continuity_on_to_the_piece_at_the_joint() {
    let (mut m, mut ed, map) = setup();
    let f = ed.feature;
    draw(
        &mut m,
        &mut ed,
        &map,
        Tool::Line,
        &[(-40.0, 0.0), (0.0, 0.0)],
    );
    ed.set_tool(Tool::ControlSpline);
    for (x, y) in [(0.0, 0.0), (15.0, -6.0), (30.0, -2.0), (45.0, 4.0)] {
        click(&mut m, &mut ed, &map, x, y);
    }
    ed.finish(&mut m);
    for x in [-20.0, -10.0] {
        draw(&mut m, &mut ed, &map, Tool::Line, &[(x, -10.0), (x, 10.0)]);
    }
    let horizontal = |m: &Model| {
        let lines = ids(m, f, |g| matches!(g, Geometry::Line { .. }));
        lines
            .into_iter()
            .filter(|&id| match resolved(m, f, id) {
                ResolvedCurve::Line { start, end } => start.y.abs() < 1e-9 && end.y.abs() < 1e-9,
                _ => false,
            })
            .collect::<Vec<_>>()
    };
    let line = horizontal(&m)[0];
    let spline = spline_id(&m, f);
    ed.selection = vec![line, spline];
    ed.constrain(&mut m, ConstraintTool::Curvature);
    assert_eq!(status(&m, f), SolveStatus::Solved);
    let joint_line = |m: &Model| {
        let cs = curvature_constraints(m, f);
        assert_eq!(cs.len(), 1, "{cs:?}");
        let (a, b, _) = cs[0];
        if a == spline {
            b
        } else {
            a
        }
    };
    assert_eq!(joint_line(&m), line);

    // Trimming the middle piece leaves a new line at the joint, which
    // takes the constraint.
    ed.set_tool(Tool::Trim);
    click(&mut m, &mut ed, &map, -15.0, 0.0);
    ed.set_tool(Tool::Select);
    let near = joint_line(&m);
    assert_ne!(near, line);
    let ResolvedCurve::Line { start, end } = resolved(&m, f, near) else {
        panic!("not a line");
    };
    assert!((start.x.max(end.x)).abs() < 1e-9 && (start.x.min(end.x) + 10.0 * MM).abs() < 1e-9);
    assert_eq!(status(&m, f), SolveStatus::Solved);
    assert_eq!(redundant(&m, f), 0);

    // Splitting that line hands it on again.
    ed.set_tool(Tool::Split);
    click(&mut m, &mut ed, &map, -5.0, 0.0);
    ed.set_tool(Tool::Select);
    let piece = joint_line(&m);
    assert_ne!(piece, near);
    let ResolvedCurve::Line { start, end } = resolved(&m, f, piece) else {
        panic!("not a line");
    };
    assert!((start.x.min(end.x) + 5.0 * MM).abs() < 1e-9);
    assert_eq!(status(&m, f), SolveStatus::Solved);
    let (_, _, k) = spline_start(&resolved(&m, f, spline));
    assert!(
        k.abs() < 1e-9,
        "still curvature continuous with a line: {k}"
    );
}
