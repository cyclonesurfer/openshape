//! Smooth closed splines (ADR 0016), driven by clicks through a top-view
//! camera: both spline tools close smoothly on their first point, the
//! closed curve is previewed, it bounds the exact area, it has no ends to
//! join, and mirror, pattern, copy and undo keep it closed.
//!
//! Shapes are checked against constructions independent of the solver's
//! clamped form: the uniform periodic B-spline on the wrapped control
//! polygon (Piegl & Tiller, *The NURBS Book*, §12.2), evaluated by de
//! Boor's algorithm on unclamped integer knots, and the periodic cubic
//! interpolant in second-derivative form (de Boor, *A Practical Guide to
//! Splines*, ch. IV). Areas are ½∮(x dy − y dx), integrated exactly by
//! Gauss–Legendre quadrature on each polynomial piece.

mod common;

use common::*;
use openshape_app::sketch_mode::{sample_resolved, ConstraintTool, SketchEditor, Tool};
use openshape_app::Model;
use openshape_core_model::sketch::{Constraint, Geometry};
use openshape_core_model::{DVec2, FeatureId, SketchId};
use openshape_kernel::Kernel;
use openshape_sketch_solver::{ResolvedCurve, SolveStatus};
use std::collections::BTreeSet;

/// Five-point Gauss–Legendre nodes and weights on [−1, 1].
const GL: [(f64, f64); 5] = [
    (0.0, 0.568_888_888_888_888_9),
    (-0.538_469_310_105_683_1, 0.478_628_670_499_366_5),
    (0.538_469_310_105_683_1, 0.478_628_670_499_366_5),
    (-0.906_179_845_938_664, 0.236_926_885_056_189_1),
    (0.906_179_845_938_664, 0.236_926_885_056_189_1),
];

/// ½∫(x y′ − y x′) over `[a, b]` of a curve given as point and derivative
/// (exact when the integrand is a polynomial of degree ≤ 9 there).
fn piece_area(a: f64, b: f64, f: impl Fn(f64) -> (DVec2, DVec2)) -> f64 {
    GL.iter()
        .map(|&(x, w)| {
            let (c, d) = f(0.5 * (a + b) + 0.5 * (b - a) * x);
            0.5 * (b - a) * w * 0.5 * c.perp_dot(d)
        })
        .sum()
}

/// de Boor's algorithm: the point at `t` of the degree-`p` B-spline on
/// `pts` and knots `k`, for `t` in `[k[p], k[pts.len()])`.
fn de_boor(pts: &[DVec2], k: &[f64], p: usize, t: f64) -> DVec2 {
    let mut s = p;
    while s + 1 < pts.len() && k[s + 1] <= t {
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

/// The point and derivative at `t` of the degree-`p` B-spline on `pts`
/// and knots `k` (the derivative is the degree p − 1 spline on the inner
/// knots, P&T eq. 3.8).
fn point_and_derivative(pts: &[DVec2], k: &[f64], p: usize, t: f64) -> (DVec2, DVec2) {
    let q: Vec<DVec2> = (0..pts.len() - 1)
        .map(|i| (pts[i + 1] - pts[i]) * (p as f64 / (k[i + p + 1] - k[i + 1])))
        .collect();
    (
        de_boor(pts, k, p, t),
        de_boor(&q, &k[1..k.len() - 1], p - 1, t),
    )
}

/// The uniform periodic B-spline of degree `p` on the closed polygon
/// `distinct` (no repeated last point): the wrapped polygon on the knots
/// 0, 1, …, n + 2p, over one period `[p, n + p]`. Returns the polygon and
/// the knots.
fn wrapped(distinct: &[DVec2], p: usize) -> (Vec<DVec2>, Vec<f64>) {
    let n = distinct.len();
    let pts: Vec<DVec2> = (0..n + p).map(|j| distinct[j % n]).collect();
    let knots: Vec<f64> = (0..=n + 2 * p).map(|i| i as f64).collect();
    (pts, knots)
}

/// The signed area inside the uniform periodic B-spline on `distinct`.
fn periodic_control_area(distinct: &[DVec2], p: usize) -> f64 {
    let (pts, k) = wrapped(distinct, p);
    (p..p + distinct.len())
        .map(|i| {
            piece_area(i as f64, i as f64 + 1.0, |t| {
                point_and_derivative(&pts, &k, p, t)
            })
        })
        .sum()
}

/// The periodic cubic interpolant of `y` at increasing `t` (`y` ends on its
/// first value again), in second-derivative form: the moments `M` from the
/// cyclic system `h_{i−1} M_{i−1} + 2(h_{i−1} + h_i) M_i + h_i M_{i+1} =
/// 6(δ_i − δ_{i−1})`, solved densely.
struct PeriodicCubic {
    t: Vec<f64>,
    y: Vec<f64>,
    m: Vec<f64>,
}

impl PeriodicCubic {
    fn new(t: &[f64], y: &[f64]) -> Self {
        let n = t.len() - 1;
        let h: Vec<f64> = (0..n).map(|i| t[i + 1] - t[i]).collect();
        let d: Vec<f64> = (0..n).map(|i| (y[i + 1] - y[i]) / h[i]).collect();
        let mut a = vec![vec![0.0; n + 1]; n];
        for i in 0..n {
            let prev = (i + n - 1) % n;
            a[i][prev] += h[prev];
            a[i][i] += 2.0 * (h[prev] + h[i]);
            a[i][(i + 1) % n] += h[i];
            a[i][n] = 6.0 * (d[i] - d[prev]);
        }
        for c in 0..n {
            let p = (c..n)
                .max_by(|&i, &j| a[i][c].abs().total_cmp(&a[j][c].abs()))
                .unwrap();
            a.swap(c, p);
            let pivot = a[c].clone();
            for row in a.iter_mut().skip(c + 1) {
                let f = row[c] / pivot[c];
                for (x, p) in row[c..].iter_mut().zip(&pivot[c..]) {
                    *x -= f * p;
                }
            }
        }
        let mut m = vec![0.0; n];
        for c in (0..n).rev() {
            let s: f64 = (c + 1..n).map(|k| a[c][k] * m[k]).sum();
            m[c] = (a[c][n] - s) / a[c][c];
        }
        Self {
            t: t.to_vec(),
            y: y.to_vec(),
            m,
        }
    }

    /// The value and derivative on piece `i` at `x`.
    fn eval(&self, i: usize, x: f64) -> (f64, f64) {
        let n = self.m.len();
        let h = self.t[i + 1] - self.t[i];
        let (m0, m1) = (self.m[i], self.m[(i + 1) % n]);
        let (a, b) = (self.t[i + 1] - x, x - self.t[i]);
        let (c0, c1) = (
            self.y[i] - m0 * h * h / 6.0,
            self.y[i + 1] - m1 * h * h / 6.0,
        );
        (
            m0 * a.powi(3) / (6.0 * h) + m1 * b.powi(3) / (6.0 * h) + c0 * a / h + c1 * b / h,
            -m0 * a * a / (2.0 * h) + m1 * b * b / (2.0 * h) + (c1 - c0) / h,
        )
    }
}

/// The closed through-point curve through `pts` (ending on the first
/// again) at `t`: its pieces' areas and a point evaluator.
struct ClosedThrough {
    x: PeriodicCubic,
    y: PeriodicCubic,
}

impl ClosedThrough {
    fn new(pts: &[DVec2], t: &[f64]) -> Self {
        let xs: Vec<f64> = pts.iter().map(|p| p.x).collect();
        let ys: Vec<f64> = pts.iter().map(|p| p.y).collect();
        Self {
            x: PeriodicCubic::new(t, &xs),
            y: PeriodicCubic::new(t, &ys),
        }
    }

    fn piece(&self, s: f64) -> usize {
        let t = &self.x.t;
        (0..t.len() - 1).rev().find(|&i| s >= t[i]).unwrap_or(0)
    }

    fn at(&self, s: f64) -> (DVec2, DVec2) {
        let i = self.piece(s);
        let (x, dx) = self.x.eval(i, s);
        let (y, dy) = self.y.eval(i, s);
        (DVec2::new(x, y), DVec2::new(dx, dy))
    }

    fn area(&self) -> f64 {
        self.x
            .t
            .windows(2)
            .enumerate()
            .map(|(i, w)| {
                piece_area(w[0], w[1], |s| {
                    let (x, dx) = self.x.eval(i, s);
                    let (y, dy) = self.y.eval(i, s);
                    (DVec2::new(x, y), DVec2::new(dx, dy))
                })
            })
            .sum()
    }
}

/// The sketch's splines: (id, points, periodic, degree or params).
#[derive(Debug, Clone)]
struct Spline {
    id: SketchId,
    points: Vec<SketchId>,
    periodic: bool,
    degree: Option<u32>,
    params: Vec<f64>,
}

fn splines(m: &Model, f: FeatureId) -> Vec<Spline> {
    m.sketch_def(f)
        .unwrap()
        .sketch
        .entities
        .iter()
        .filter_map(|e| match &e.geometry {
            Geometry::ControlSpline {
                control,
                degree,
                periodic,
            } => Some(Spline {
                id: e.id,
                points: control.clone(),
                periodic: *periodic,
                degree: Some(*degree),
                params: Vec::new(),
            }),
            Geometry::InterpolatedSpline {
                through,
                params,
                periodic,
            } => Some(Spline {
                id: e.id,
                points: through.clone(),
                periodic: *periodic,
                degree: None,
                params: params.clone(),
            }),
            _ => None,
        })
        .collect()
}

/// The solved positions of `ids`.
fn solved(m: &Model, f: FeatureId, ids: &[SketchId]) -> Vec<DVec2> {
    let sketch = &m.sketch_result(f).unwrap().solution.sketch;
    ids.iter().map(|id| sketch.point(*id).unwrap()).collect()
}

/// The solved clamped curve of spline `id`: (degree, knots, controls).
fn clamped(m: &Model, f: FeatureId, id: SketchId) -> (usize, Vec<f64>, Vec<DVec2>) {
    match curves(m, f, true).into_iter().find(|(k, _)| *k == id) {
        Some((
            _,
            ResolvedCurve::Spline {
                degree,
                knots,
                points,
            },
        )) => (degree, knots, points),
        other => panic!("{other:?}"),
    }
}

/// `a` within `rel` (relative) of `b`.
fn near(a: f64, b: f64, rel: f64) -> bool {
    (a - b).abs() <= rel * b.abs()
}

/// The first and last derivative of a clamped spline (P&T eq. 3.7): equal
/// for a curve with no seam.
fn end_derivatives(p: usize, k: &[f64], pts: &[DVec2]) -> (DVec2, DVec2) {
    let n = pts.len();
    (
        (pts[1] - pts[0]) * (p as f64 / (k[p + 1] - k[1])),
        (pts[n - 1] - pts[n - 2]) * (p as f64 / (k[n + p - 1] - k[n - 1])),
    )
}

/// The lengths of the seam and crease edges of the body `feature` made
/// (seams are where a surface continues smoothly; the display hides them).
fn seams_and_creases(m: &Model, feature: FeatureId) -> (f64, f64) {
    let b = m
        .bodies()
        .bodies
        .iter()
        .find(|b| b.feature == feature)
        .expect("a body");
    let opts = openshape_kernel::TessellationOptions {
        chord_tolerance: 1e-7,
        max_angle: 1f64.to_radians(),
    };
    let mesh = m.kernel().tessellate(&b.body, &opts).unwrap();
    let mut out = (0.0, 0.0);
    for (e, seam) in mesh.edges.iter().zip(&mesh.seams) {
        let l: f64 = e.windows(2).map(|w| (w[1] - w[0]).length()).sum();
        if *seam {
            out.0 += l;
        } else {
            out.1 += l;
        }
    }
    out
}

/// The length of a clamped spline, by a dense polyline.
fn length(k: &[f64], p: usize, ctrl: &[DVec2]) -> f64 {
    let n = 20_000;
    let pts: Vec<DVec2> = (0..=n)
        .map(|i| de_boor(ctrl, k, p, (i as f64 / n as f64).min(1.0 - 1e-15)))
        .collect();
    pts.windows(2).map(|w| (w[1] - w[0]).length()).sum()
}

/// Asserts the extruded closed spline's side is one smooth wall: its only
/// creases are the two rims, and the two lines where the kernel cuts the
/// curve in two are seams.
fn assert_smooth_wall(m: &Model, e: FeatureId, perimeter: f64, depth: f64) {
    let (seam, crease) = seams_and_creases(m, e);
    assert!(
        near(seam, 2.0 * depth, 1e-9),
        "seams {seam} vs {}",
        2.0 * depth
    );
    assert!(
        near(crease, 2.0 * perimeter, 1e-6),
        "creases {crease} vs {}",
        2.0 * perimeter
    );
}

const KITE: [(f64, f64); 4] = [(-30.0, -20.0), (-30.0, 30.0), (30.0, 30.0), (30.0, -20.0)];

/// Clicks `points` with `tool`, then clicks a pixel from the first point
/// again; returns the one spline made.
fn closed_with(m: &mut Model, ed: &mut SketchEditor, tool: Tool, points: &[(f64, f64)]) -> Spline {
    let map = view();
    let f = ed.feature;
    let before = splines(m, f).len();
    ed.set_tool(tool);
    for &(x, y) in points {
        click(m, ed, &map, x, y);
    }
    assert_eq!(splines(m, f).len(), before, "nothing until it closes");
    let (x, y) = points[0];
    click(m, ed, &map, x + 0.5, y + 0.3);
    let all = splines(m, f);
    assert_eq!(
        all.len(),
        before + 1,
        "the closing click finishes the curve"
    );
    all.last().unwrap().clone()
}

#[test]
fn the_control_spline_tool_closes_smoothly_on_its_first_point() {
    let (mut m, mut ed, map) = setup();
    let f = ed.feature;
    let s = closed_with(&mut m, &mut ed, Tool::ControlSpline, &KITE);
    assert!(s.periodic, "closing makes a smooth closed spline");
    assert_eq!((s.points.len(), s.degree), (5, Some(3)));
    assert_eq!(
        s.points.first(),
        s.points.last(),
        "closed on the same point"
    );
    assert_eq!(s.points.iter().collect::<BTreeSet<_>>().len(), 4);
    assert_eq!(status(&m, f), SolveStatus::Solved);
    assert_eq!(dof(&m, f), 8, "two freedoms per distinct control point");

    let pts = solved(&m, f, &s.points[..4]);
    assert!((pts[0] - DVec2::new(-30.0, -20.0) * MM).length() < 1e-9);
    // The solver's clamped curve is the wrapped uniform spline, point for
    // point: the clamped parameter s ∈ [0, 1] is u = 3 + 4s.
    let (p, k, ctrl) = clamped(&m, f, s.id);
    assert_eq!((p, ctrl.len()), (3, 7));
    let (wpts, wk) = wrapped(&pts, 3);
    for i in 0..=200 {
        let t = i as f64 / 200.0;
        let a = de_boor(&ctrl, &k, 3, t.min(1.0 - 1e-15));
        let b = de_boor(&wpts, &wk, 3, (3.0 + 4.0 * t).min(7.0 - 1e-14));
        assert!((a - b).length() < 1e-12, "at {t}: {a} vs {b}");
    }
    // No seam: it leaves its start as it arrives at its end, and it does
    // not pass through its first control point.
    assert_eq!(ctrl[0], ctrl[6], "it closes exactly");
    let (d0, d1) = end_derivatives(3, &k, &ctrl);
    assert!((d0 - d1).length() <= 1e-9 * d0.length(), "{d0} vs {d1}");
    assert!((ctrl[0] - pts[0]).length() > 5.0 * MM);

    let expected = periodic_control_area(&pts, 3).abs();
    let areas = region_areas(&m, f);
    assert_eq!(areas.len(), 1, "{areas:?}");
    assert!(near(areas[0], expected, 1e-9), "{} vs {expected}", areas[0]);
    // The closed curve alone bounds a solid.
    let e = extrude(&mut m, &map, &[(0.0, 5.0)], "5 mm");
    let v = volume(&m, e);
    assert!(
        near(v, expected * 5.0 * MM, 1e-7),
        "{v} vs {}",
        expected * 5.0 * MM
    );
    assert_smooth_wall(&m, e, length(&k, 3, &ctrl), 5.0 * MM);
}

#[test]
fn three_control_points_close_into_a_quadratic() {
    let (mut m, mut ed, _) = setup();
    let f = ed.feature;
    let tri = [(-30.0, -20.0), (30.0, -20.0), (0.0, 30.0)];
    let s = closed_with(&mut m, &mut ed, Tool::ControlSpline, &tri);
    assert!(s.periodic);
    assert_eq!((s.points.len(), s.degree), (4, Some(2)));
    assert_eq!(dof(&m, f), 6);
    let pts = solved(&m, f, &s.points[..3]);
    let (p, k, ctrl) = clamped(&m, f, s.id);
    let (d0, d1) = end_derivatives(p, &k, &ctrl);
    assert!((d0 - d1).length() <= 1e-9 * d0.length(), "{d0} vs {d1}");
    let expected = periodic_control_area(&pts, 2).abs();
    let areas = region_areas(&m, f);
    assert_eq!(areas.len(), 1, "{areas:?}");
    assert!(near(areas[0], expected, 1e-9), "{} vs {expected}", areas[0]);
}

#[test]
fn the_spline_tool_closes_smoothly_through_its_points() {
    let (mut m, mut ed, map) = setup();
    let f = ed.feature;
    let s = closed_with(&mut m, &mut ed, Tool::Spline, &KITE);
    assert!(s.periodic);
    assert_eq!(s.points.len(), 5);
    assert_eq!(s.points.first(), s.points.last());
    assert_eq!(dof(&m, f), 8);
    let pts = solved(&m, f, &s.points);
    // Chord-length parameters over the closed polygon, closing chord
    // included: sides of 50, 60, 50 and 60 mm.
    let want = [0.0, 50.0 / 220.0, 110.0 / 220.0, 160.0 / 220.0, 1.0];
    for (a, b) in s.params.iter().zip(want) {
        assert!((a - b).abs() < 1e-6, "{:?}", s.params);
    }
    assert_eq!(s.params[4], 1.0);

    // The solver's curve is the periodic cubic interpolant: through every
    // point at its parameter, the same shape everywhere, no seam.
    let (p, k, ctrl) = clamped(&m, f, s.id);
    assert_eq!((p, ctrl.len()), (3, 7));
    let exact = ClosedThrough::new(&pts, &s.params);
    for (q, &t) in pts.iter().zip(&s.params) {
        let c = de_boor(&ctrl, &k, 3, t.min(1.0 - 1e-15));
        assert!((c - *q).length() < 1e-12, "at {t}: {c} vs {q}");
    }
    for i in 0..=200 {
        let t = i as f64 / 200.0;
        let a = de_boor(&ctrl, &k, 3, t.min(1.0 - 1e-15));
        let (b, _) = exact.at(t);
        assert!((a - b).length() < 1e-12, "at {t}: {a} vs {b}");
    }
    let (d0, d1) = end_derivatives(3, &k, &ctrl);
    assert!((d0 - d1).length() <= 1e-9 * d0.length(), "{d0} vs {d1}");

    let expected = exact.area().abs();
    let areas = region_areas(&m, f);
    assert_eq!(areas.len(), 1, "{areas:?}");
    assert!(near(areas[0], expected, 1e-9), "{} vs {expected}", areas[0]);
    let e = extrude(&mut m, &map, &[(0.0, 5.0)], "4 mm");
    let v = volume(&m, e);
    assert!(
        near(v, expected * 4.0 * MM, 1e-7),
        "{v} vs {}",
        expected * 4.0 * MM
    );
    assert_smooth_wall(&m, e, length(&k, 3, &ctrl), 4.0 * MM);
}

#[test]
fn enter_still_finishes_an_open_spline() {
    let (mut m, mut ed, map) = setup();
    let f = ed.feature;
    for (tool, dx) in [(Tool::Spline, -60.0), (Tool::ControlSpline, 60.0)] {
        ed.set_tool(tool);
        for &(x, y) in &KITE {
            click(&mut m, &mut ed, &map, 0.5 * x + dx, 0.5 * y + 60.0);
        }
        ed.finish(&mut m);
    }
    let s = splines(&m, f);
    assert_eq!(s.len(), 2);
    for sp in &s {
        assert!(!sp.periodic);
        assert_eq!(sp.points.len(), 4);
        assert_ne!(sp.points.first(), sp.points.last());
    }
    assert!(region_areas(&m, f).is_empty());
}

#[test]
fn the_closed_curve_is_previewed_before_the_closing_click() {
    let (mut m, mut ed, map) = setup();
    let f = ed.feature;
    for tool in [Tool::Spline, Tool::ControlSpline] {
        ed.set_tool(tool);
        for &(x, y) in &KITE {
            click(&mut m, &mut ed, &map, x, y);
        }
        // Away from the first point: an open curve to the pointer.
        hover(&m, &mut ed, &map, 0.0, -40.0);
        let open = ed.rubber_band();
        let curve = open.last().unwrap();
        assert!((curve[0] - curve[curve.len() - 1]).length() > 10.0 * MM);
        // Over the first point: the smooth closed curve.
        hover(&m, &mut ed, &map, -29.5, -19.7);
        let shown = ed.rubber_band();
        let curve = shown.last().unwrap();
        let n = curve.len();
        assert!(n > 32, "{tool:?}: {n}");
        assert!((curve[0] - curve[n - 1]).length() < 1e-9, "{tool:?} closes");
        let (into, out) = (curve[n - 1] - curve[n - 2], curve[1] - curve[0]);
        assert!(
            into.normalize().dot(out.normalize()) > 0.99,
            "{tool:?}: no corner at the seam"
        );
        if tool == Tool::ControlSpline {
            // With the closed control polygon.
            assert_eq!(shown.len(), 2);
            assert_eq!(shown[0].len(), 5);
            assert_eq!(shown[0][0], shown[0][4]);
        } else {
            assert_eq!(shown.len(), 1);
            // Through the clicked points (to within the sampling's sag).
            for &(x, y) in &KITE {
                let q = DVec2::new(x, y) * MM;
                let d = curve
                    .windows(2)
                    .map(|w| {
                        let (a, b) = (w[0], w[1]);
                        let t = ((q - a).dot(b - a) / (b - a).length_squared()).clamp(0.0, 1.0);
                        (a + (b - a) * t - q).length()
                    })
                    .fold(f64::MAX, f64::min);
                assert!(d < 0.05 * MM, "{tool:?}: {d}");
            }
        }
        assert!(splines(&m, f).is_empty(), "previews make nothing");
        // The preview is the curve the closing click makes.
        click(&mut m, &mut ed, &map, -29.5, -19.7);
        let s = splines(&m, f);
        assert_eq!(s.len(), 1);
        let (degree, knots, points) = clamped(&m, f, s[0].id);
        let made = sample_resolved(&ResolvedCurve::Spline {
            degree,
            knots,
            points,
        });
        assert_eq!(made.len(), curve.len());
        for (a, b) in made.iter().zip(curve) {
            assert!((*a - *b).length() < 1e-12, "{tool:?}: {a} vs {b}");
        }
        assert!(m.undo());
        ed.set_tool(Tool::Select);
    }
}

#[test]
fn a_closed_spline_has_no_ends_to_join() {
    let (mut m, mut ed, map) = setup();
    let f = ed.feature;
    let s = closed_with(&mut m, &mut ed, Tool::Spline, &KITE);
    // A line from the spline's first point outwards.
    ed.set_tool(Tool::Line);
    click(&mut m, &mut ed, &map, -30.0, -20.0);
    click(&mut m, &mut ed, &map, -60.0, -20.0);
    ed.finish(&mut m);
    ed.set_tool(Tool::Select);
    let line = m
        .sketch_def(f)
        .unwrap()
        .sketch
        .entities
        .iter()
        .find(|e| matches!(e.geometry, Geometry::Line { .. }))
        .unwrap()
        .id;
    let n = m.sketch_def(f).unwrap().sketch.constraints.len();
    for tool in [ConstraintTool::Tangent, ConstraintTool::Curvature] {
        ed.selection = vec![s.id, line];
        assert!(!ed.can_constrain(&m, tool), "{tool:?}");
        ed.constrain(&mut m, tool);
        let message = ed.message.clone().unwrap_or_default();
        assert!(message.contains("has no end points"), "{tool:?}: {message}");
        assert_eq!(m.sketch_def(f).unwrap().sketch.constraints.len(), n);
    }
    // A point on it may still slide all the way round.
    let point = m.sketch_def(f).unwrap().sketch.entities.len();
    ed.set_tool(Tool::Point);
    click(&mut m, &mut ed, &map, 0.0, 0.0);
    ed.set_tool(Tool::Select);
    let p = m.sketch_def(f).unwrap().sketch.entities[point].id;
    ed.selection = vec![p, s.id];
    ed.constrain(&mut m, ConstraintTool::Coincident);
    assert_eq!(
        constraints(&m, f, |c| matches!(c, Constraint::PointOnCurve { .. })),
        1,
        "{:?}",
        ed.message
    );
    assert_eq!(status(&m, f), SolveStatus::Solved);
}

#[test]
fn mirror_pattern_and_copy_keep_a_spline_closed() {
    let (mut m, mut ed, map) = setup();
    let f = ed.feature;
    let ring = [(20.0, 10.0), (50.0, 10.0), (50.0, 40.0), (20.0, 40.0)];
    let s = closed_with(&mut m, &mut ed, Tool::Spline, &ring);
    let one = region_areas(&m, f)[0];
    // Mirrored about a vertical line at x = 0.
    ed.set_tool(Tool::Line);
    click(&mut m, &mut ed, &map, 0.0, 60.0);
    click(&mut m, &mut ed, &map, 0.0, 90.0);
    ed.finish(&mut m);
    ed.set_tool(Tool::Select);
    ed.selection = vec![s.id];
    ed.set_tool(Tool::Mirror);
    click(&mut m, &mut ed, &map, 0.0, 75.0);
    assert_eq!(ed.message.as_deref(), Some("Mirrored 1 item"));
    // Turned half round the origin by a circular pattern…
    ed.set_tool(Tool::Select);
    ed.selection = vec![s.id];
    ed.pattern_count = 2;
    ed.set_tool(Tool::CircularPattern);
    click(&mut m, &mut ed, &map, 0.0, 0.0);
    assert_eq!(
        ed.message.as_deref(),
        Some("Patterned 1 item into 2 instances")
    );
    // …and copied down.
    ed.set_tool(Tool::Select);
    ed.selection = vec![s.id];
    ed.set_tool(Tool::Move);
    ed.transform_copy = true;
    ed.transform_offset = DVec2::new(0.0, -80.0) * MM;
    ed.apply_typed_transform(&mut m);
    assert_eq!(ed.message.as_deref(), Some("Copied 1 item"));

    let all = splines(&m, f);
    assert_eq!(all.len(), 4);
    assert!(all.iter().all(|sp| sp.periodic), "{all:?}");
    assert!(all.iter().all(|sp| sp.points.first() == sp.points.last()));
    assert_eq!(status(&m, f), SolveStatus::Solved);
    let areas = region_areas(&m, f);
    assert_eq!(areas.len(), 4, "{areas:?}");
    assert!(
        areas.iter().all(|a| near(*a, one, 1e-9)),
        "{one}: {areas:?}"
    );
}

#[test]
fn one_undo_takes_a_closed_spline_back() {
    let (mut m, mut ed, _) = setup();
    let f = ed.feature;
    let before = m.sketch_def(f).unwrap().sketch.clone();
    closed_with(&mut m, &mut ed, Tool::ControlSpline, &KITE);
    assert!(m.undo());
    assert_eq!(m.sketch_def(f).unwrap().sketch, before);
    assert!(m.redo());
    let s = splines(&m, f);
    assert_eq!(s.len(), 1);
    assert!(s[0].periodic);
}
