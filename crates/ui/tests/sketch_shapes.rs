//! Shape tools and line inference (design note 0007 §4–5), driven through
//! the editor exactly as the viewport drives it: entity and constraint
//! counts, degrees of freedom, no redundant constraint, and exact regions.

mod common;

use common::*;
use opendatum_app::sketch_mode::{Inference, Tool};
use opendatum_core_model::sketch::{Constraint, Geometry};
use opendatum_core_model::DVec2;
use opendatum_sketch_solver::{ResolvedCurve, SolveStatus};
use std::f64::consts::PI;

fn lines(m: &opendatum_app::Model, f: opendatum_core_model::FeatureId) -> Vec<(DVec2, DVec2)> {
    curves(m, f, false)
        .into_iter()
        .filter_map(|(_, c)| match c {
            ResolvedCurve::Line { start, end } => Some((start, end)),
            _ => None,
        })
        .collect()
}

#[test]
fn centre_rectangle_is_symmetric_about_its_centre() {
    let (mut m, mut ed, map) = setup();
    let f = ed.feature;
    draw(
        &mut m,
        &mut ed,
        &map,
        Tool::CenterRectangle,
        &[(10.0, 5.0), (40.0, 25.0)],
    );
    assert_eq!(count(&m, f, |g| matches!(g, Geometry::Line { .. })), 6);
    assert_eq!(construction(&m, f), 2);
    assert_eq!(
        constraints(&m, f, |c| matches!(c, Constraint::Midpoint { .. })),
        1
    );
    assert_eq!(status(&m, f), SolveStatus::Solved);
    // Centre, width and height.
    assert_eq!(dof(&m, f), 4);
    let a = region_areas(&m, f);
    assert_eq!(a.len(), 1);
    assert!(close(a[0], 60.0 * 40.0 * MM * MM), "{a:?}");
    // The far corner mirrors the clicked one through the centre.
    let far = DVec2::new(-20.0, -15.0) * MM;
    assert!(lines(&m, f)
        .iter()
        .any(|(s, e)| (*s - far).length() < 1e-9 || (*e - far).length() < 1e-9));
}

#[test]
fn aligned_rectangle_keeps_right_angles_at_any_angle() {
    let (mut m, mut ed, map) = setup();
    let f = ed.feature;
    let side = DVec2::new(30.0, 20.0);
    let n = DVec2::new(-side.y, side.x).normalize();
    let w = DVec2::new(35.0, 25.0) + n * 15.0 + side.normalize() * 4.0;
    draw(
        &mut m,
        &mut ed,
        &map,
        Tool::AlignedRectangle,
        &[(5.0, 5.0), (35.0, 25.0), (w.x, w.y)],
    );
    assert_eq!(count(&m, f, |g| matches!(g, Geometry::Line { .. })), 4);
    assert_eq!(
        constraints(&m, f, |c| matches!(c, Constraint::Perpendicular { .. })),
        1
    );
    assert_eq!(
        constraints(&m, f, |c| matches!(c, Constraint::Parallel { .. })),
        2
    );
    assert_eq!(status(&m, f), SolveStatus::Solved);
    // Position, angle, length and width.
    assert_eq!(dof(&m, f), 5);
    let a = region_areas(&m, f);
    assert!(close(a[0], side.length() * 15.0 * MM * MM), "{a:?}");
    let l = lines(&m, f);
    for i in 0..4 {
        let (d0, d1) = (l[i].1 - l[i].0, l[(i + 1) % 4].1 - l[(i + 1) % 4].0);
        assert!(d0.dot(d1).abs() < 1e-12 * d0.length() * d1.length() * 1e3);
    }
}

#[test]
fn slot_is_tangent_all_round_with_equal_ends() {
    let (mut m, mut ed, map) = setup();
    let f = ed.feature;
    draw(
        &mut m,
        &mut ed,
        &map,
        Tool::Slot,
        &[(-20.0, 10.0), (20.0, 10.0), (3.0, 22.0)],
    );
    assert_eq!(count(&m, f, |g| matches!(g, Geometry::Arc { .. })), 2);
    assert_eq!(count(&m, f, |g| matches!(g, Geometry::Line { .. })), 3);
    assert_eq!(construction(&m, f), 1);
    assert_eq!(
        constraints(&m, f, |c| matches!(c, Constraint::Tangent { .. })),
        4
    );
    assert_eq!(status(&m, f), SolveStatus::Solved);
    // Two centres and the radius.
    assert_eq!(dof(&m, f), 5);
    let r = 12.0;
    let a = region_areas(&m, f);
    assert_eq!(a.len(), 1);
    assert!(
        close(a[0], (2.0 * r * 40.0 + PI * r * r) * MM * MM),
        "{a:?}"
    );
    let e = extrude(&mut m, &map, &[(0.0, 10.0)], "5");
    let v = volume(&m, e);
    assert!(close(v, a[0] * 5.0 * MM), "{v}");
}

#[test]
fn polygons_are_regular_with_four_freedoms() {
    for sides in [3u32, 5, 6, 8] {
        let (mut m, mut ed, map) = setup();
        let f = ed.feature;
        ed.polygon_sides = sides;
        draw(
            &mut m,
            &mut ed,
            &map,
            Tool::Polygon,
            &[(10.0, 10.0), (40.0, 10.0)],
        );
        assert_eq!(
            count(&m, f, |g| matches!(g, Geometry::Line { .. })),
            sides as usize
        );
        assert_eq!(construction(&m, f), 1, "the circumcircle");
        assert_eq!(status(&m, f), SolveStatus::Solved, "{sides} sides");
        // Centre, size and rotation.
        assert_eq!(dof(&m, f), 4, "{sides} sides");
        let n = sides as f64;
        let expected = 0.5 * n * 30.0 * 30.0 * (2.0 * PI / n).sin() * MM * MM;
        let a = region_areas(&m, f);
        assert!(close(a[0], expected), "{sides} sides: {a:?} vs {expected}");
        let l = lines(&m, f);
        let len0 = (l[0].1 - l[0].0).length();
        for (s, e) in &l {
            assert!(((*e - *s).length() - len0).abs() < 1e-12);
        }
    }
}

#[test]
fn polygon_sides_are_limited() {
    let (mut m, mut ed, map) = setup();
    let f = ed.feature;
    ed.polygon_sides = 2;
    draw(
        &mut m,
        &mut ed,
        &map,
        Tool::Polygon,
        &[(10.0, 10.0), (40.0, 10.0)],
    );
    assert_eq!(count(&m, f, |_| true), 0);
    assert!(ed.message.as_deref().unwrap_or("").contains("3 to 64"));
}

#[test]
fn tangent_arcs_continue_lines_and_arcs() {
    let (mut m, mut ed, map) = setup();
    let f = ed.feature;
    ed.set_tool(Tool::Line);
    click(&mut m, &mut ed, &map, -40.0, 20.0);
    click(&mut m, &mut ed, &map, 0.0, 20.0);
    ed.escape();
    // From the line's end to (20, 40): a quarter circle about (0, 40).
    draw(
        &mut m,
        &mut ed,
        &map,
        Tool::TangentArc,
        &[(0.0, 20.0), (20.0, 40.0)],
    );
    let arcs: Vec<(DVec2, f64)> = curves(&m, f, false)
        .into_iter()
        .filter_map(|(_, c)| match c {
            ResolvedCurve::Arc { center, radius, .. } => Some((center, radius)),
            _ => None,
        })
        .collect();
    assert_eq!(arcs.len(), 1);
    assert!(
        (arcs[0].0 - DVec2::new(0.0, 40.0) * MM).length() < 1e-8,
        "{arcs:?}"
    );
    assert!(close(arcs[0].1, 20.0 * MM), "{arcs:?}");
    assert_eq!(
        constraints(&m, f, |c| matches!(c, Constraint::Tangent { .. })),
        1
    );
    assert_eq!(status(&m, f), SolveStatus::Solved);
    // Line: 4 − horizontal; arc: 5 − coincident 2 − tangent 1.
    assert_eq!(dof(&m, f), 5);

    // Continue the arc with another tangent arc turning the other way.
    draw(
        &mut m,
        &mut ed,
        &map,
        Tool::TangentArc,
        &[(20.0, 40.0), (40.0, 60.0)],
    );
    assert_eq!(
        constraints(&m, f, |c| matches!(c, Constraint::Tangent { .. })),
        2
    );
    assert_eq!(status(&m, f), SolveStatus::Solved);
    assert_eq!(dof(&m, f), 7);
    let arcs: Vec<(DVec2, f64)> = curves(&m, f, false)
        .into_iter()
        .filter_map(|(_, c)| match c {
            ResolvedCurve::Arc { center, radius, .. } => Some((center, radius)),
            _ => None,
        })
        .collect();
    // Leaving (20, 40) straight up and ending at (40, 60): turns right
    // about (40, 40).
    assert!(arcs
        .iter()
        .any(|(c, r)| (*c - DVec2::new(40.0, 40.0) * MM).length() < 1e-8 && close(*r, 20.0 * MM)));
}

#[test]
fn tangent_arc_needs_a_curve_end() {
    let (mut m, mut ed, map) = setup();
    let f = ed.feature;
    draw(
        &mut m,
        &mut ed,
        &map,
        Tool::TangentArc,
        &[(10.0, 10.0), (30.0, 30.0)],
    );
    assert_eq!(count(&m, f, |_| true), 0);
    assert!(ed
        .message
        .as_deref()
        .unwrap_or("")
        .contains("end of a line"));
}

/// Draws a line from `a` to `b` (millimetres) and ends the chain.
fn line(
    m: &mut opendatum_app::Model,
    ed: &mut opendatum_app::sketch_mode::SketchEditor,
    map: &opendatum_app::sketch_mode::ScreenMap,
    a: DVec2,
    b: DVec2,
) {
    ed.set_tool(Tool::Line);
    click(m, ed, map, a.x, a.y);
    click(m, ed, map, b.x, b.y);
    ed.escape();
}

#[test]
fn lines_infer_perpendicular_and_parallel_to_the_line_they_continue() {
    let (mut m, mut ed, map) = setup();
    let f = ed.feature;
    let a = DVec2::new(10.0, 10.0);
    let b = a + DVec2::from_angle(30f64.to_radians()) * 40.0;
    line(&mut m, &mut ed, &map, a, b);
    // From the end, 2° off perpendicular.
    let c = b + DVec2::from_angle(122f64.to_radians()) * 30.0;
    ed.set_tool(Tool::Line);
    click(&mut m, &mut ed, &map, b.x, b.y);
    hover(&m, &mut ed, &map, c.x, c.y);
    assert!(matches!(ed.inferred, Some(Inference::Perpendicular(_))));
    click(&mut m, &mut ed, &map, c.x, c.y);
    // 2.5° off parallel to the line just drawn (the chain continues from
    // its solved end).
    let c = ed.clicks[0].at / MM;
    let d = c + DVec2::from_angle(122.5f64.to_radians()) * 20.0;
    hover(&m, &mut ed, &map, d.x, d.y);
    assert!(matches!(ed.inferred, Some(Inference::Parallel(_))));
    click(&mut m, &mut ed, &map, d.x, d.y);
    ed.escape();
    assert_eq!(
        constraints(&m, f, |c| matches!(c, Constraint::Perpendicular { .. })),
        1
    );
    assert_eq!(
        constraints(&m, f, |c| matches!(c, Constraint::Parallel { .. })),
        1
    );
    assert_eq!(status(&m, f), SolveStatus::Solved);
    let l = lines(&m, f);
    let dir = |i: usize| (l[i].1 - l[i].0).normalize();
    assert!(dir(0).dot(dir(1)).abs() < 1e-12);
    assert!(dir(1).perp_dot(dir(2)).abs() < 1e-12);
}

#[test]
fn lines_infer_parallel_to_a_line_the_pointer_passed() {
    let (mut m, mut ed, map) = setup();
    let f = ed.feature;
    let a = DVec2::new(-60.0, -40.0);
    let b = a + DVec2::from_angle(30f64.to_radians()) * 40.0;
    line(&mut m, &mut ed, &map, a, b);
    ed.set_tool(Tool::Line);
    click(&mut m, &mut ed, &map, 20.0, -50.0);
    // Pass over the first line, then come back near its direction.
    let mid = (a + b) * 0.5;
    hover(&m, &mut ed, &map, mid.x, mid.y);
    let end = DVec2::new(20.0, -50.0) + DVec2::from_angle(31.5f64.to_radians()) * 40.0;
    click(&mut m, &mut ed, &map, end.x, end.y);
    ed.escape();
    assert_eq!(
        constraints(&m, f, |c| matches!(c, Constraint::Parallel { .. })),
        1
    );
    let l = lines(&m, f);
    let (d0, d1) = ((l[0].1 - l[0].0).normalize(), (l[1].1 - l[1].0).normalize());
    assert!(d0.perp_dot(d1).abs() < 1e-12);
    // Further off than 3°: nothing inferred.
    ed.set_tool(Tool::Line);
    click(&mut m, &mut ed, &map, 20.0, 40.0);
    hover(&m, &mut ed, &map, mid.x, mid.y);
    let off = DVec2::new(20.0, 40.0) + DVec2::from_angle(34f64.to_radians()) * 40.0;
    hover(&m, &mut ed, &map, off.x, off.y);
    assert_eq!(ed.inferred, None);
}

#[test]
fn lines_infer_tangency_leaving_an_arc() {
    let (mut m, mut ed, map) = setup();
    let f = ed.feature;
    // An arc from (10, 10) to (40, 30) through (20, 30): centre (30, 12.5),
    // so its tangent at (40, 30) is oblique (150°), far from the axes.
    draw(
        &mut m,
        &mut ed,
        &map,
        Tool::Arc,
        &[(10.0, 10.0), (40.0, 30.0), (20.0, 30.0)],
    );
    let arc = curves(&m, f, false)[0].1.clone();
    let (center, end) = match arc {
        ResolvedCurve::Arc {
            center,
            radius,
            start_angle,
            sweep_angle,
        } => {
            let p0 = center + DVec2::from_angle(start_angle) * radius;
            let p1 = center + DVec2::from_angle(start_angle + sweep_angle) * radius;
            // The clicked end at (40, 30).
            let e = if (p0 - DVec2::new(40.0, 30.0) * MM).length() < 1e-6 {
                p0
            } else {
                p1
            };
            (center, e)
        }
        _ => panic!("not an arc"),
    };
    let normal = (end - center).normalize();
    let tangent = DVec2::new(-normal.y, normal.x);
    let e = end / MM;
    let target = e + DVec2::from_angle(tangent.to_angle() + 2f64.to_radians()) * 25.0;
    ed.set_tool(Tool::Line);
    click(&mut m, &mut ed, &map, e.x, e.y);
    hover(&m, &mut ed, &map, target.x, target.y);
    assert!(
        matches!(ed.inferred, Some(Inference::Tangent(_))),
        "{:?}",
        ed.inferred
    );
    click(&mut m, &mut ed, &map, target.x, target.y);
    ed.escape();
    assert_eq!(
        constraints(&m, f, |c| matches!(c, Constraint::Tangent { .. })),
        1
    );
    assert_eq!(status(&m, f), SolveStatus::Solved);
}
