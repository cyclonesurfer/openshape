//! The ellipse and elliptical-arc tools (design note 0008 §7–§8), driven by
//! clicks through a top-view camera. Areas and volumes are checked against
//! the closed forms πab, (ab/2)(Δ − sin Δ) and πab·h on the solved axes.

mod common;

use common::*;
use openshape_app::sketch_mode::{ConstraintTool, SnapTarget, Tool};
use openshape_app::Model;
use openshape_core_model::sketch::{Constraint, Geometry, TangentKind};
use openshape_core_model::{DVec2, FeatureId, SketchId};
use openshape_sketch_solver::{ResolvedCurve, SolveStatus};
use std::f64::consts::{FRAC_PI_6, PI, TAU};

/// `a` within `rel` (relative) of `b`.
fn near(a: f64, b: f64, rel: f64) -> bool {
    (a - b).abs() <= rel * b.abs()
}

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

fn is_ellipse(g: &Geometry<SketchId>) -> bool {
    matches!(g, Geometry::Ellipse { .. } | Geometry::EllipticalArc { .. })
}

/// The solved ellipses and elliptical arcs: centre, semi-axes, and for an
/// arc its start anomaly and signed sweep.
type Solved = (DVec2, DVec2, DVec2, Option<(f64, f64)>);

fn ellipses(m: &Model, f: FeatureId) -> Vec<Solved> {
    curves(m, f, true)
        .into_iter()
        .filter_map(|(_, c)| match c {
            ResolvedCurve::Ellipse {
                center,
                major,
                minor,
            } => Some((center, major, minor, None)),
            ResolvedCurve::EllipticalArc {
                center,
                major,
                minor,
                start_param,
                sweep_param,
            } => Some((center, major, minor, Some((start_param, sweep_param)))),
            _ => None,
        })
        .collect()
}

/// E(t) on a solved ellipse.
fn point(e: &Solved, t: f64) -> DVec2 {
    e.0 + e.1 * t.cos() + e.2 * t.sin()
}

/// The ends of a solved elliptical arc.
fn ends(e: &Solved) -> (DVec2, DVec2) {
    let (t0, sweep) = e.3.expect("an arc");
    (point(e, t0), point(e, t0 + sweep))
}

/// Twice the area of the ellipse's axis rectangle's quarter: ab.
fn ab(e: &Solved) -> f64 {
    e.1.length() * e.2.length()
}

/// cos of the angle between the semi-axes.
fn skew(e: &Solved) -> f64 {
    e.1.dot(e.2) / (e.1.length() * e.2.length())
}

/// The point at anomaly `deg` of the nominal ellipse the tests draw: centre
/// (10, 5) mm, semi-axes (30, 0) and (0, `b`) mm.
fn nominal(b: f64, deg: f64) -> (f64, f64) {
    scaled(b, deg, 1.0)
}

/// The nominal ellipse's point at anomaly `deg`, `k` times as far from the
/// centre (the same anomaly).
fn scaled(b: f64, deg: f64, k: f64) -> (f64, f64) {
    let t = deg.to_radians();
    (10.0 + 30.0 * k * t.cos(), 5.0 + b * k * t.sin())
}

/// Draws the nominal ellipse with semi-minor `b` mm (negative: on the
/// right of the major axis).
fn draw_ellipse(m: &mut Model, ed: &mut openshape_app::sketch_mode::SketchEditor, b: f64) {
    let map = view();
    // The third click is off the minor axis: only its distance across the
    // major axis counts.
    draw(
        m,
        ed,
        &map,
        Tool::Ellipse,
        &[(10.0, 5.0), (40.0, 5.0), (17.0, 5.0 + b)],
    );
}

#[test]
fn three_clicks_make_an_ellipse_with_its_axes() {
    let (mut m, mut ed, map) = setup();
    let f = ed.feature;
    ed.set_tool(Tool::Ellipse);
    assert_eq!(ed.hint(), "Click the centre");
    click(&mut m, &mut ed, &map, 10.0, 5.0);
    assert_eq!(ed.hint(), "Click the end of the first axis");
    click(&mut m, &mut ed, &map, 40.0, 5.0);
    assert_eq!(ed.hint(), "Click to set the second axis");
    // The rubber band shows the ellipse and both semi-axes.
    hover(&m, &mut ed, &map, 17.0, 20.0);
    let band = ed.rubber_band();
    assert_eq!(band.len(), 2);
    let e0 = DVec2::new(10.0, 5.0) * MM;
    for p in &band[0] {
        let q = (*p - e0) / MM;
        let implicit = (q.x / 30.0).powi(2) + (q.y / 15.0).powi(2);
        assert!((implicit - 1.0).abs() < 1e-5, "{q:?}");
    }
    click(&mut m, &mut ed, &map, 17.0, 20.0);
    assert!(ed.clicks.is_empty(), "the tool starts again");

    assert_eq!(count(&m, f, |g| matches!(g, Geometry::Ellipse { .. })), 1);
    assert_eq!(count(&m, f, |g| matches!(g, Geometry::Point { .. })), 3);
    assert_eq!(construction(&m, f), 2, "two construction semi-axes");
    assert_eq!(constraints(&m, f, |_| true), 0);
    assert_eq!(status(&m, f), SolveStatus::Solved);
    // Centre 2, orientation 1, two radii.
    assert_eq!(dof(&m, f), 5);

    let e = ellipses(&m, f).remove(0);
    assert!((e.0 - e0).length() < 1e-8);
    assert!(close(e.1.length(), 30.0 * MM), "{:?}", e.1);
    assert!(close(e.2.length(), 15.0 * MM), "{:?}", e.2);
    assert!(skew(&e).abs() < 1e-12);
    assert!(e.1.perp_dot(e.2) > 0.0, "the second axis is on the left");

    // One region of area πab, exactly; extruded, a volume of πab·h.
    let areas = region_areas(&m, f);
    assert_eq!(areas.len(), 1, "{areas:?}");
    assert!(near(areas[0], PI * ab(&e), 1e-12), "{areas:?}");
    let x = extrude(&mut m, &map, &[(10.0, 0.0)], "12 mm");
    let v = volume(&m, x);
    assert!(near(v, PI * ab(&e) * 12.0 * MM, 1e-9), "{v}");
}

#[test]
fn a_second_axis_on_the_first_is_refused() {
    let (mut m, mut ed, map) = setup();
    let f = ed.feature;
    ed.set_tool(Tool::Ellipse);
    click(&mut m, &mut ed, &map, 10.0, 5.0);
    click(&mut m, &mut ed, &map, 40.0, 5.0);
    // On the line of the first axis: no second axis.
    click(&mut m, &mut ed, &map, 25.0, 5.0);
    assert_eq!(count(&m, f, is_ellipse), 0);
    assert_eq!(ed.clicks.len(), 2, "the clicks so far are kept");
    let message = ed.message.clone().unwrap_or_default();
    assert!(message.contains("second axis has no length"), "{message}");
    // A click off it then makes the ellipse.
    click(&mut m, &mut ed, &map, 25.0, -8.0);
    assert_eq!(count(&m, f, is_ellipse), 1);
    let e = ellipses(&m, f).remove(0);
    assert!(close(e.2.length(), 13.0 * MM));
    assert!(e.1.perp_dot(e.2) < 0.0, "the second axis is on the right");
}

/// Draws the arc of the nominal ellipse (semi-minor `b`) from anomaly
/// `from` to `to`, hovering through `via` on the way.
fn draw_arc(
    m: &mut Model,
    ed: &mut openshape_app::sketch_mode::SketchEditor,
    b: f64,
    from: f64,
    via: &[f64],
    to: f64,
) {
    let map = view();
    ed.set_tool(Tool::EllipticalArc);
    for (x, y) in [(10.0, 5.0), (40.0, 5.0), (17.0, 5.0 + b)] {
        click(m, ed, &map, x, y);
    }
    assert_eq!(ed.hint(), "Click the start point on the ellipse");
    // The start is taken at its eccentric anomaly: off the curve is fine.
    let (x, y) = scaled(b, from, 1.1);
    click(m, ed, &map, x, y);
    assert_eq!(ed.hint(), "Move around the ellipse and click the end point");
    for t in via {
        let (x, y) = scaled(b, *t, 0.9);
        hover(m, ed, &map, x, y);
    }
    let (x, y) = nominal(b, to);
    click(m, ed, &map, x, y);
    ed.set_tool(Tool::Select);
}

/// Closes an arc of the nominal ellipse with a line between its ends and
/// checks the region against the elliptic segment (ab/2)(Δ − sin Δ).
fn check_segment(m: &mut Model, ed: &mut openshape_app::sketch_mode::SketchEditor, b: f64) {
    let f = ed.feature;
    let e = ellipses(m, f).remove(0);
    let (s, t) = ends(&e);
    let map = view();
    draw(
        m,
        ed,
        &map,
        Tool::Line,
        &[(s.x / MM, s.y / MM), (t.x / MM, t.y / MM)],
    );
    assert_eq!(status(m, f), SolveStatus::Solved);
    let e = ellipses(m, f).remove(0);
    let sweep = e.3.unwrap().1.abs();
    let areas = region_areas(m, f);
    assert_eq!(areas.len(), 1, "{areas:?} (b = {b})");
    let expected = 0.5 * ab(&e) * (sweep - sweep.sin());
    assert!(near(areas[0], expected, 1e-11), "{areas:?} vs {expected}");
}

#[test]
fn an_elliptical_arc_runs_the_way_the_pointer_turned() {
    for (from, via, to) in [(30.0, vec![90.0], 150.0), (150.0, vec![90.0], 30.0)] {
        let (mut m, mut ed, _) = setup();
        let f = ed.feature;
        draw_arc(&mut m, &mut ed, 15.0, from, &via, to);
        assert_eq!(
            count(&m, f, |g| matches!(g, Geometry::EllipticalArc { .. })),
            1
        );
        assert_eq!(count(&m, f, |g| matches!(g, Geometry::Point { .. })), 5);
        assert_eq!(construction(&m, f), 2);
        assert_eq!(status(&m, f), SolveStatus::Solved);
        // The ellipse's 5, and the two ends' anomalies.
        assert_eq!(dof(&m, f), 7);
        // Either way it is the arc over the top, counter-clockwise from
        // 30° to 150°: drawn clockwise, its ends were swapped.
        let e = ellipses(&m, f).remove(0);
        let (t0, sweep) = e.3.unwrap();
        assert!((t0 - FRAC_PI_6).abs() < 1e-5, "{t0}");
        assert!((sweep - TAU / 3.0).abs() < 1e-5, "{sweep}");
        // The ends lie on the ellipse exactly (the solver holds them).
        let def = m.sketch_result(f).unwrap().solution.sketch.clone();
        let arc = ids(&m, f, is_ellipse)[0];
        let Geometry::EllipticalArc { start, end, .. } = def.entity(arc).unwrap().geometry else {
            panic!()
        };
        let (s, t) = ends(&e);
        assert!((def.point(start).unwrap() - s).length() < 1e-12);
        assert!((def.point(end).unwrap() - t).length() < 1e-12);
        check_segment(&mut m, &mut ed, 15.0);
    }
}

#[test]
fn arcs_of_left_handed_ellipses_run_counter_clockwise_too() {
    // The second axis on the right: increasing anomaly runs clockwise. The
    // pointer turns counter-clockwise over the top, from anomaly −30° to
    // −150° through −90°.
    let (mut m, mut ed, _) = setup();
    let f = ed.feature;
    draw_arc(
        &mut m,
        &mut ed,
        -10.0,
        -30.0,
        &[-60.0, -90.0, -120.0],
        -150.0,
    );
    assert_eq!(status(&m, f), SolveStatus::Solved);
    let e = ellipses(&m, f).remove(0);
    assert!(e.1.perp_dot(e.2) < 0.0);
    let (t0, sweep) = e.3.unwrap();
    assert!((t0 + FRAC_PI_6).abs() < 1e-5, "{t0}");
    assert!((sweep + TAU / 3.0).abs() < 1e-5, "{sweep}");
    // The arc's middle is the top of the ellipse.
    let mid = point(&e, t0 + sweep / 2.0);
    assert!((mid - DVec2::new(10.0, 15.0) * MM).length() < 1e-7, "{mid}");
    check_segment(&mut m, &mut ed, -10.0);

    // The long way round: the pointer turns clockwise through the bottom,
    // 240° of anomaly.
    let (mut m, mut ed, _) = setup();
    let f = ed.feature;
    draw_arc(
        &mut m,
        &mut ed,
        -10.0,
        -30.0,
        &[0.0, 60.0, 120.0, 180.0],
        210.0,
    );
    let e = ellipses(&m, f).remove(0);
    let (_, sweep) = e.3.unwrap();
    assert!((sweep.abs() - 2.0 * TAU / 3.0).abs() < 1e-5, "{sweep}");
    let mid = point(&e, e.3.unwrap().0 + sweep / 2.0);
    assert!((mid - DVec2::new(10.0, -5.0) * MM).length() < 1e-7, "{mid}");
    check_segment(&mut m, &mut ed, -10.0);
}

#[test]
fn points_snap_onto_an_ellipse_exactly() {
    let (mut m, mut ed, map) = setup();
    let f = ed.feature;
    draw_ellipse(&mut m, &mut ed, 15.0);
    let ellipse = ids(&m, f, is_ellipse)[0];
    ed.set_tool(Tool::Point);
    // Near the curve at about 50°, between the drawn polyline's vertices.
    let (x, y) = nominal(15.0, 50.3);
    hover(&m, &mut ed, &map, x + 0.3, y + 0.2);
    let s = ed.cursor.expect("a snap");
    assert_eq!(s.target, SnapTarget::Curve(ellipse));
    let e = ellipses(&m, f).remove(0);
    let q = s.at - e.0;
    let implicit =
        (q.dot(e.1) / e.1.length_squared()).powi(2) + (q.dot(e.2) / e.2.length_squared()).powi(2);
    assert!((implicit - 1.0).abs() < 1e-12, "{implicit}");
    click(&mut m, &mut ed, &map, x + 0.3, y + 0.2);
    assert_eq!(
        constraints(&m, f, |c| matches!(c, Constraint::PointOnCurve { .. })),
        1
    );
    assert_eq!(status(&m, f), SolveStatus::Solved);
    assert_eq!(dof(&m, f), 5 + 1);
}

#[test]
fn a_line_made_tangent_to_an_ellipse_touches_it() {
    let (mut m, mut ed, map) = setup();
    let f = ed.feature;
    draw_ellipse(&mut m, &mut ed, 15.0);
    // A sloping line above the ellipse, clear of it.
    draw(
        &mut m,
        &mut ed,
        &map,
        Tool::Line,
        &[(-40.0, 30.0), (60.0, 40.0)],
    );
    let line = ids(&m, f, |g| matches!(g, Geometry::Line { .. }))
        .into_iter()
        .last()
        .unwrap();
    let ellipse = ids(&m, f, is_ellipse)[0];
    ed.selection = vec![line, ellipse];
    ed.constrain(&mut m, ConstraintTool::Tangent);
    assert_eq!(
        constraints(&m, f, |c| matches!(
            c,
            Constraint::Tangent {
                kind: TangentKind::LineCircle { center_left: false },
                ..
            }
        )),
        1,
        "{:?}",
        ed.message
    );
    assert_eq!(status(&m, f), SolveStatus::Solved);
    // Independent check: the line's distance from the centre equals the
    // ellipse's support value in the line's normal direction.
    let e = ellipses(&m, f).remove(0);
    let (s, t) = curves(&m, f, false)
        .into_iter()
        .find_map(|(id, c)| match c {
            ResolvedCurve::Line { start, end } if id == line => Some((start, end)),
            _ => None,
        })
        .unwrap();
    let n = (t - s).perp().normalize();
    let h = DVec2::new(e.1.dot(n), e.2.dot(n)).length();
    let gap = n.dot(e.0 - s).abs() - h;
    assert!(gap.abs() < 1e-12, "{gap}");

    // A circle and an ellipse with no shared end are refused, by name.
    draw(
        &mut m,
        &mut ed,
        &map,
        Tool::Circle,
        &[(80.0, 0.0), (90.0, 0.0)],
    );
    let circle = ids(&m, f, |g| matches!(g, Geometry::Circle { .. }))[0];
    ed.selection = vec![ellipse, circle];
    ed.constrain(&mut m, ConstraintTool::Tangent);
    let message = ed.message.clone().unwrap_or_default();
    assert!(
        message.contains("ellipse can be tangent to a line"),
        "{message}"
    );
    for tool in [ConstraintTool::Equal, ConstraintTool::Concentric] {
        ed.selection = vec![ellipse, circle];
        ed.constrain(&mut m, tool);
        let message = ed.message.clone().unwrap_or_default();
        assert!(message.contains("needs"), "{message}");
    }
}

#[test]
fn curves_continue_tangentially_from_an_elliptical_arc() {
    let (mut m, mut ed, map) = setup();
    let f = ed.feature;
    draw_arc(&mut m, &mut ed, 15.0, 30.0, &[90.0], 150.0);
    let e = ellipses(&m, f).remove(0);
    let (_, end) = ends(&e);
    // A tangent arc from the arc's end, out to the left and down.
    draw(
        &mut m,
        &mut ed,
        &map,
        Tool::TangentArc,
        &[(end.x / MM, end.y / MM), (-30.0, -10.0)],
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
    let e = ellipses(&m, f).remove(0);
    let (t0, sweep) = e.3.unwrap();
    let t1 = t0 + sweep;
    let along = (e.2 * t1.cos() - e.1 * t1.sin()).normalize();
    let (center, radius) = curves(&m, f, false)
        .into_iter()
        .find_map(|(_, c)| match c {
            ResolvedCurve::Arc { center, radius, .. } => Some((center, radius)),
            _ => None,
        })
        .unwrap();
    let p = point(&e, t1);
    // The circle's radius at the joint is the ellipse's normal there.
    assert!((((p - center).length() - radius) / radius).abs() < 1e-12);
    assert!((p - center).normalize().dot(along).abs() < 1e-9);

    // A line drawn on along the arc's tangent is inferred tangent.
    let (mut m, mut ed, map) = setup();
    let f = ed.feature;
    draw_arc(&mut m, &mut ed, 15.0, 30.0, &[90.0], 150.0);
    let e = ellipses(&m, f).remove(0);
    let (start, _) = ends(&e);
    let t0 = e.3.unwrap().0;
    // Leaving the start, the curve runs backwards: −E′(t₀).
    let back = -(e.2 * t0.cos() - e.1 * t0.sin()).normalize();
    let to = start + back * 20.0 * MM;
    draw(
        &mut m,
        &mut ed,
        &map,
        Tool::Line,
        &[(start.x / MM, start.y / MM), (to.x / MM, to.y / MM)],
    );
    assert_eq!(
        constraints(&m, f, |c| matches!(c, Constraint::Tangent { .. })),
        1
    );
    assert_eq!(status(&m, f), SolveStatus::Solved);
}

#[test]
fn mirrored_arcs_have_the_other_handedness() {
    let (mut m, mut ed, map) = setup();
    let f = ed.feature;
    draw_arc(&mut m, &mut ed, 15.0, 30.0, &[90.0], 150.0);
    ed.construction = true;
    draw(
        &mut m,
        &mut ed,
        &map,
        Tool::Line,
        &[(60.0, -50.0), (60.0, 50.0)],
    );
    ed.construction = false;
    let before = dof(&m, f);
    ed.set_tool(Tool::Mirror);
    click(&mut m, &mut ed, &map, 60.0, 0.0);
    // The arc (its top is the minor axis's end point).
    let (x, y) = nominal(15.0, 60.0);
    click(&mut m, &mut ed, &map, x, y);
    assert_eq!(ed.message.as_deref(), Some("Mirrored 1 item"));
    assert_eq!(status(&m, f), SolveStatus::Solved);
    assert_eq!(dof(&m, f), before, "the image adds no freedom");
    let es = ellipses(&m, f);
    assert_eq!(es.len(), 2);
    let (src, img) = (&es[0], &es[1]);
    let x = 60.0 * MM;
    let reflect = |p: DVec2| DVec2::new(2.0 * x - p.x, p.y);
    assert!((img.0 - reflect(src.0)).length() < 1e-12);
    assert!(img.1.perp_dot(img.2) < 0.0, "a reflection flips handedness");
    // Still counter-clockwise: the image starts at the source's end.
    let ((s0, s1), (i0, i1)) = (ends(src), ends(img));
    assert!((i0 - reflect(s1)).length() < 1e-12, "{i0} {s1}");
    assert!((i1 - reflect(s0)).length() < 1e-12, "{i1} {s0}");
    assert!((img.3.unwrap().1 + src.3.unwrap().1).abs() < 1e-12);
    // Moving the source moves the image.
    let def = m.sketch_def(f).unwrap().sketch.clone();
    let arc = ids(&m, f, is_ellipse)[0];
    let Geometry::EllipticalArc { major, .. } = def.entity(arc).unwrap().geometry else {
        panic!()
    };
    ed.set_tool(Tool::Select);
    ed.drag_start(&mut m, &map, at(&map, 40.0, 5.0));
    ed.drag_to(&m, &map, at(&map, 45.0, 0.0));
    ed.drag_end(&mut m, &map, at(&map, 45.0, 0.0));
    assert_eq!(status(&m, f), SolveStatus::Solved);
    let p = m
        .sketch_result(f)
        .unwrap()
        .solution
        .sketch
        .point(major)
        .unwrap();
    assert!((p - DVec2::new(45.0, 0.0) * MM).length() < 1e-6, "{p}");
    let es = ellipses(&m, f);
    let (src, img) = (&es[0], &es[1]);
    assert!((img.0 - reflect(src.0)).length() < 1e-12);
    assert!((img.1 - DVec2::new(-src.1.x, src.1.y)).length() < 1e-12);
}

#[test]
fn patterned_ellipses_keep_their_shape() {
    let (mut m, mut ed, map) = setup();
    let f = ed.feature;
    draw_ellipse(&mut m, &mut ed, 15.0);
    let before = dof(&m, f);
    let ellipse = ids(&m, f, is_ellipse)[0];
    ed.selection = vec![ellipse];
    ed.set_tool(Tool::LinearPattern);
    click(&mut m, &mut ed, &map, 10.0, 5.0);
    click(&mut m, &mut ed, &map, 10.0, 40.0);
    assert_eq!(
        ed.message.as_deref(),
        Some("Patterned 1 item into 3 instances")
    );
    ed.prompt = None;
    assert_eq!(status(&m, f), SolveStatus::Solved);
    // The step's direction line adds its free end.
    assert_eq!(dof(&m, f), before + 2);
    let es = ellipses(&m, f);
    assert_eq!(es.len(), 3);
    for (k, e) in es.iter().enumerate() {
        let step = es[1].0 - es[0].0;
        assert!((e.0 - (es[0].0 + step * k as f64)).length() < 1e-12);
        assert!((e.1 - es[0].1).length() < 1e-12 && (e.2 - es[0].2).length() < 1e-12);
    }
    let areas = region_areas(&m, f);
    assert_eq!(areas.len(), 3, "{areas:?}");
    for a in areas {
        assert!(near(a, PI * ab(&es[0]), 1e-12), "{a}");
    }

    // A circular pattern turns the axes with the copies.
    let (mut m, mut ed, map) = setup();
    let f = ed.feature;
    draw_ellipse(&mut m, &mut ed, 15.0);
    ed.selection = ids(&m, f, is_ellipse);
    ed.pattern_count = 4;
    ed.set_tool(Tool::CircularPattern);
    click(&mut m, &mut ed, &map, -40.0, 5.0);
    ed.prompt = None;
    assert_eq!(status(&m, f), SolveStatus::Solved);
    let es = ellipses(&m, f);
    assert_eq!(es.len(), 4);
    let hub = DVec2::new(-40.0, 5.0) * MM;
    for (k, e) in es.iter().enumerate() {
        let turn = DVec2::from_angle(k as f64 * TAU / 4.0);
        assert!((e.0 - (hub + turn.rotate(es[0].0 - hub))).length() < 1e-9);
        assert!((e.1 - turn.rotate(es[0].1)).length() < 1e-9);
        assert!((e.2 - turn.rotate(es[0].2)).length() < 1e-9);
    }
}

#[test]
fn editing_tools_that_cannot_handle_ellipses_say_so() {
    let (mut m, mut ed, map) = setup();
    let f = ed.feature;
    draw_ellipse(&mut m, &mut ed, 15.0);
    // A line across it.
    draw(
        &mut m,
        &mut ed,
        &map,
        Tool::Line,
        &[(0.0, -30.0), (0.0, 40.0)],
    );
    let before = m.sketch_def(f).unwrap().sketch.clone();
    let (bx, by) = nominal(15.0, 270.0);
    for (tool, says) in [
        (Tool::Trim, "Trimming ellipses is not built yet"),
        (Tool::Extend, "Only lines and arcs can be extended"),
        (Tool::Split, "Only lines and arcs can be split"),
        (Tool::Offset, "its offset is not an ellipse"),
    ] {
        ed.set_tool(tool);
        click(&mut m, &mut ed, &map, bx, by);
        let message = ed.message.clone().unwrap_or_default();
        assert!(message.contains(says), "{tool:?}: {message}");
        assert_eq!(m.sketch_def(f).unwrap().sketch, before, "{tool:?}");
    }

    // A fillet where an elliptical arc meets a line.
    let (mut m, mut ed, map) = setup();
    let f = ed.feature;
    draw_arc(&mut m, &mut ed, 15.0, 30.0, &[90.0], 150.0);
    let e = ellipses(&m, f).remove(0);
    let (_, end) = ends(&e);
    draw(
        &mut m,
        &mut ed,
        &map,
        Tool::Line,
        &[(end.x / MM, end.y / MM), (end.x / MM, -20.0)],
    );
    let before = m.sketch_def(f).unwrap().sketch.clone();
    ed.set_tool(Tool::SketchFillet);
    click(&mut m, &mut ed, &map, end.x / MM, end.y / MM);
    let message = ed.message.clone().unwrap_or_default();
    assert!(message.contains("an elliptical arc ends here"), "{message}");
    assert_eq!(m.sketch_def(f).unwrap().sketch, before);
}

#[test]
fn dragging_an_axis_end_turns_the_ellipse_in_one_undo_step() {
    let (mut m, mut ed, map) = setup();
    let f = ed.feature;
    draw_ellipse(&mut m, &mut ed, 15.0);
    let before = ellipses(&m, f).remove(0);
    let ops = m.session().log().len();
    ed.set_tool(Tool::Select);
    ed.drag_start(&mut m, &map, at(&map, 40.0, 5.0));
    ed.drag_to(&m, &map, at(&map, 35.0, 20.0));
    ed.drag_to(&m, &map, at(&map, 30.0, 25.0));
    ed.drag_end(&mut m, &map, at(&map, 30.0, 25.0));
    assert_eq!(m.session().log().len(), ops + 1, "one undo step");
    assert_eq!(status(&m, f), SolveStatus::Solved);
    let after = ellipses(&m, f).remove(0);
    let end = after.0 + after.1;
    assert!((end - DVec2::new(30.0, 25.0) * MM).length() < 1e-6, "{end}");
    // The axes stay perpendicular and the handedness is kept.
    assert!(skew(&after).abs() < 1e-9, "{}", skew(&after));
    assert!(after.1.perp_dot(after.2) > 0.0);

    assert!(m.undo());
    let undone = ellipses(&m, f).remove(0);
    assert!((undone.1 - before.1).length() < 1e-12);

    // Dragging the curve itself moves the whole ellipse.
    let (x, y) = nominal(15.0, 60.0);
    ed.drag_start(&mut m, &map, at(&map, x, y));
    ed.drag_to(&m, &map, at(&map, x + 10.0, y - 5.0));
    ed.drag_end(&mut m, &map, at(&map, x + 10.0, y - 5.0));
    let moved = ellipses(&m, f).remove(0);
    assert!((moved.0 - before.0 - DVec2::new(10.0, -5.0) * MM).length() < 1e-6);
    assert!((moved.1 - before.1).length() < 1e-9);
}

#[test]
fn deleting_an_ellipse_takes_its_axes_and_points() {
    let (mut m, mut ed, _) = setup();
    let f = ed.feature;
    draw_ellipse(&mut m, &mut ed, 15.0);
    let (mut m2, mut ed2, _) = setup();
    draw_arc(&mut m2, &mut ed2, 15.0, 30.0, &[90.0], 150.0);
    for (m, ed) in [(&mut m, &mut ed), (&mut m2, &mut ed2)] {
        let f = ed.feature;
        ed.selection = ids(m, f, is_ellipse);
        ed.delete_selection(m);
        assert_eq!(m.sketch_def(f).unwrap().sketch.entities.len(), 0);
    }
    // An axis line on its own is just a line.
    draw_ellipse(&mut m, &mut ed, 15.0);
    let axis = ids(&m, f, |g| matches!(g, Geometry::Line { .. }))[0];
    ed.selection = vec![axis];
    ed.delete_selection(&mut m);
    assert_eq!(count(&m, f, is_ellipse), 1);
    assert_eq!(construction(&m, f), 1);
    assert_eq!(status(&m, f), SolveStatus::Solved);
}
