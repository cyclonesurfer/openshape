//! The control-point spline tool, driven by clicks through a top-view
//! camera. Areas are checked against an independent de Boor evaluation of
//! the solved control points, integrated exactly with Gauss–Legendre
//! quadrature (Green's theorem on each polynomial span).

mod common;

use common::*;
use opendatum_app::sketch_mode::{ConstraintTool, Tool};
use opendatum_core_model::sketch::Geometry;
use opendatum_core_model::{DVec2, FeatureId, SketchId};
use opendatum_sketch_solver::SolveStatus;

/// Clamped knots with uniform interior knots, written out here rather than
/// taken from the solver so the check is independent.
fn knots(n: usize, p: usize) -> Vec<f64> {
    let spans = n - p;
    let mut k = vec![0.0; p + 1];
    k.extend((1..spans).map(|i| i as f64 / spans as f64));
    k.extend(std::iter::repeat_n(1.0, p + 1));
    k
}

/// de Boor's algorithm: the point at `t` of the degree-`p` B-spline.
fn de_boor(pts: &[DVec2], k: &[f64], p: usize, t: f64) -> DVec2 {
    let n = pts.len();
    let mut s = p;
    while s + 1 < n && k[s + 1] <= t {
        s += 1;
    }
    let mut d: Vec<DVec2> = (0..=p).map(|j| pts[j + s - p]).collect();
    for r in 1..=p {
        for j in (r..=p).rev() {
            let i = j + s - p;
            let a = (t - k[i]) / (k[i + p + 1 - r] - k[i]);
            d[j] = d[j - 1] * (1.0 - a) + d[j] * a;
        }
    }
    d[p]
}

/// The signed area ½∮(x dy − y dx) swept by the spline (exact: on each
/// span the integrand is a polynomial of degree 2p − 1 ≤ 9).
fn spline_area(pts: &[DVec2], p: usize) -> f64 {
    const X: [f64; 5] = [
        0.0,
        -0.538_469_310_105_683_1,
        0.538_469_310_105_683_1,
        -0.906_179_845_938_664,
        0.906_179_845_938_664,
    ];
    const W: [f64; 5] = [
        0.568_888_888_888_888_9,
        0.478_628_670_499_366_5,
        0.478_628_670_499_366_5,
        0.236_926_885_056_189_1,
        0.236_926_885_056_189_1,
    ];
    let k = knots(pts.len(), p);
    // The derivative is a degree p − 1 spline on the inner knots.
    let q: Vec<DVec2> = (0..pts.len() - 1)
        .map(|i| (pts[i + 1] - pts[i]) * (p as f64 / (k[i + p + 1] - k[i + 1])))
        .collect();
    let dk = &k[1..k.len() - 1];
    let mut area = 0.0;
    for w in k.windows(2).filter(|w| w[1] > w[0]) {
        let (a, b) = (w[0], w[1]);
        for (x, wt) in X.iter().zip(W) {
            let t = 0.5 * (a + b) + 0.5 * (b - a) * x;
            let c = de_boor(pts, &k, p, t);
            let d = de_boor(&q, dk, p - 1, t);
            area += 0.5 * (b - a) * wt * 0.5 * (c.x * d.y - c.y * d.x);
        }
    }
    area
}

/// ½(x₀y₁ − x₁y₀): the line segment's share of the same integral.
fn segment_area(a: DVec2, b: DVec2) -> f64 {
    0.5 * (a.x * b.y - b.x * a.y)
}

/// The sketch's control-point splines: (control ids, degree).
fn splines(m: &opendatum_app::Model, f: FeatureId) -> Vec<(Vec<SketchId>, u32)> {
    m.sketch_def(f)
        .unwrap()
        .sketch
        .entities
        .iter()
        .filter_map(|e| match &e.geometry {
            Geometry::ControlSpline {
                control, degree, ..
            } => Some((control.clone(), *degree)),
            _ => None,
        })
        .collect()
}

/// The solved positions of `ids`.
fn solved(m: &opendatum_app::Model, f: FeatureId, ids: &[SketchId]) -> Vec<DVec2> {
    let sketch = &m.sketch_result(f).unwrap().solution.sketch;
    ids.iter().map(|id| sketch.point(*id).unwrap()).collect()
}

/// `a` within `rel` (relative) of `b`.
fn near(a: f64, b: f64, rel: f64) -> bool {
    (a - b).abs() <= rel * b.abs()
}

const ARCH: [(f64, f64); 5] = [
    (-30.0, 5.0),
    (-20.0, 30.0),
    (0.0, 40.0),
    (20.0, 30.0),
    (30.0, 5.0),
];

#[test]
fn an_open_spline_closed_by_a_line_bounds_the_exact_area() {
    let (mut m, mut ed, map) = setup();
    let f = ed.feature;
    ed.set_tool(Tool::ControlSpline);
    for (x, y) in ARCH {
        click(&mut m, &mut ed, &map, x, y);
    }
    ed.finish(&mut m);
    draw(
        &mut m,
        &mut ed,
        &map,
        Tool::Line,
        &[(30.0, 5.0), (-30.0, 5.0)],
    );
    assert_eq!(status(&m, f), SolveStatus::Solved);

    let s = splines(&m, f);
    assert_eq!(s.len(), 1);
    let (control, degree) = &s[0];
    assert_eq!((control.len(), *degree), (5, 3));
    // Five free control points; the line's ends are the curve's ends.
    assert_eq!(dof(&m, f), 10);

    let pts = solved(&m, f, control);
    let expected = (spline_area(&pts, 3) + segment_area(pts[4], pts[0])).abs();
    let areas = region_areas(&m, f);
    assert_eq!(areas.len(), 1, "{areas:?}");
    assert!(near(areas[0], expected, 1e-9), "{} vs {expected}", areas[0]);

    let e = extrude(&mut m, &map, &[(0.0, 20.0)], "10 mm");
    let v = volume(&m, e);
    assert!(
        near(v, expected * 10.0 * MM, 1e-7),
        "{v} vs {}",
        expected * 0.01
    );
}

#[test]
fn a_corner_closure_is_two_ends_made_coincident() {
    // Clicking the first point again closes a spline smoothly (ADR 0016,
    // tests/sketch_closed_spline.rs); a closure with a corner ends near
    // the start and joins the ends with a coincident constraint.
    let (mut m, mut ed, map) = setup();
    let f = ed.feature;
    ed.set_tool(Tool::ControlSpline);
    for (x, y) in [
        (-30.0, -20.0),
        (-30.0, 30.0),
        (30.0, 30.0),
        (30.0, -20.0),
        (-24.0, -20.0),
    ] {
        click(&mut m, &mut ed, &map, x, y);
    }
    ed.finish(&mut m);
    let s = splines(&m, f);
    assert_eq!(s.len(), 1);
    let (control, degree) = s[0].clone();
    assert_eq!((control.len(), degree), (5, 3));
    assert!(region_areas(&m, f).is_empty(), "open");
    ed.set_tool(Tool::Select);
    ed.selection = vec![control[0], control[4]];
    ed.constrain(&mut m, ConstraintTool::Coincident);
    assert_eq!(status(&m, f), SolveStatus::Solved);
    assert_eq!(dof(&m, f), 8);

    let pts = solved(&m, f, &control);
    assert!((pts[0] - pts[4]).length() < 1e-12);
    let expected = spline_area(&pts, 3).abs();
    let areas = region_areas(&m, f);
    assert_eq!(areas.len(), 1, "{areas:?}");
    assert!(near(areas[0], expected, 1e-9), "{} vs {expected}", areas[0]);

    // The closed curve alone bounds a solid.
    let e = extrude(&mut m, &map, &[(0.0, 10.0)], "5 mm");
    let v = volume(&m, e);
    assert!(
        near(v, expected * 5.0 * MM, 1e-7),
        "{v} vs {}",
        expected * 0.005
    );
}

#[test]
fn three_points_make_a_quadratic() {
    let (mut m, mut ed, map) = setup();
    let f = ed.feature;
    ed.set_tool(Tool::ControlSpline);
    for (x, y) in [(-20.0, 0.0), (0.0, 30.0), (20.0, 0.0)] {
        click(&mut m, &mut ed, &map, x, y);
    }
    ed.finish(&mut m);
    draw(
        &mut m,
        &mut ed,
        &map,
        Tool::Line,
        &[(20.0, 0.0), (-20.0, 0.0)],
    );
    let s = splines(&m, f);
    assert_eq!(s.len(), 1);
    let (control, degree) = &s[0];
    assert_eq!((control.len(), *degree), (3, 2));
    // A parabolic segment: two thirds of its control triangle.
    let pts = solved(&m, f, control);
    let triangle = 0.5 * (pts[2] - pts[0]).perp_dot(pts[1] - pts[0]).abs();
    let areas = region_areas(&m, f);
    assert_eq!(areas.len(), 1, "{areas:?}");
    assert!(near(areas[0], triangle * 2.0 / 3.0, 1e-9), "{areas:?}");
}

#[test]
fn two_points_are_not_a_spline() {
    let (mut m, mut ed, map) = setup();
    let f = ed.feature;
    ed.set_tool(Tool::ControlSpline);
    click(&mut m, &mut ed, &map, -20.0, 0.0);
    click(&mut m, &mut ed, &map, 20.0, 0.0);
    ed.finish(&mut m);
    assert!(splines(&m, f).is_empty());
    assert_eq!(count(&m, f, |g| matches!(g, Geometry::Point { .. })), 0);
    let message = ed.message.clone().unwrap_or_default();
    assert!(message.contains("at least 3 points"), "{message}");
}

#[test]
fn dragging_a_control_point_reshapes_the_curve_in_one_undo_step() {
    let (mut m, mut ed, map) = setup();
    let f = ed.feature;
    ed.set_tool(Tool::ControlSpline);
    for (x, y) in ARCH {
        click(&mut m, &mut ed, &map, x, y);
    }
    ed.finish(&mut m);
    ed.set_tool(Tool::Select);
    let (control, _) = splines(&m, f).remove(0);
    let before = solved(&m, f, &control);
    let area_before = spline_area(&before, 3);
    let ops = m.session().log().len();

    ed.drag_start(&mut m, &map, at(&map, 0.0, 40.0));
    ed.drag_to(&m, &map, at(&map, 0.0, 50.0));
    ed.drag_to(&m, &map, at(&map, 0.0, 60.0));
    ed.drag_end(&mut m, &map, at(&map, 0.0, 60.0));
    assert_eq!(m.session().log().len(), ops + 1, "one undo step");
    let after = solved(&m, f, &control);
    assert!(
        (after[2] - DVec2::new(0.0, 60.0) * MM).length() < 1e-6,
        "{after:?}"
    );
    for i in [0, 1, 3, 4] {
        assert!((after[i] - before[i]).length() < 1e-9, "point {i} moved");
    }
    // The curve follows its control point: the swept area grows.
    assert!(spline_area(&after, 3).abs() > area_before.abs() * 1.05);

    assert!(m.undo());
    let undone = solved(&m, f, &control);
    assert!((undone[2] - before[2]).length() < 1e-12);
}

#[test]
fn a_mirrored_spline_bounds_the_mirrored_area() {
    let (mut m, mut ed, map) = setup();
    let f = ed.feature;
    ed.set_tool(Tool::ControlSpline);
    for (x, y) in ARCH {
        click(&mut m, &mut ed, &map, x, y);
    }
    ed.finish(&mut m);
    draw(
        &mut m,
        &mut ed,
        &map,
        Tool::Line,
        &[(30.0, 5.0), (-30.0, 5.0)],
    );
    ed.construction = true;
    draw(
        &mut m,
        &mut ed,
        &map,
        Tool::Line,
        &[(-40.0, -50.0), (-40.0, 50.0)],
    );
    ed.construction = false;
    let before = dof(&m, f);

    // The arch's top (it is symmetric about x = 0) and its base line.
    let nominal: Vec<DVec2> = ARCH.iter().map(|(x, y)| DVec2::new(*x, *y)).collect();
    let top = de_boor(&nominal, &knots(5, 3), 3, 0.5);
    ed.set_tool(Tool::Mirror);
    click(&mut m, &mut ed, &map, -40.0, 0.0);
    click(&mut m, &mut ed, &map, top.x, top.y);
    click(&mut m, &mut ed, &map, 0.0, 5.0);
    assert_eq!(ed.message.as_deref(), Some("Mirrored 1 item"));

    let s = splines(&m, f);
    assert_eq!(s.len(), 2);
    assert_eq!(status(&m, f), SolveStatus::Solved);
    assert_eq!(dof(&m, f), before, "images add no freedom");
    let source = solved(&m, f, &s[0].0);
    let image = solved(&m, f, &s[1].0);
    let axis_x = -40.0 * MM;
    for (p, q) in source.iter().zip(&image) {
        let reflected = DVec2::new(2.0 * axis_x - p.x, p.y);
        assert!((reflected - *q).length() < 1e-9, "{p:?} → {q:?}");
    }
    let areas = region_areas(&m, f);
    assert_eq!(areas.len(), 2, "{areas:?}");
    assert!(near(areas[0], areas[1], 1e-9), "{areas:?}");
}
