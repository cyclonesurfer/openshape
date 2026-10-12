//! Smooth closed splines (ADR 0016): the resolved curves against
//! independent constructions, smoothness through the seam, degrees of
//! freedom, points sliding across the seam, refusals, mirror images, and
//! documents written before the flag.

mod common;

use common::{assert_satisfied, de_boor, dist_to_spline, pos, resolved, B};
use openshape_sketch_solver::{
    solve, Constraint, DVec2, Geometry, ResolvedCurve, Sketch, SolveError, SolveOptions,
    SolveStatus, TangentKind,
};

fn opts() -> SolveOptions<u32> {
    SolveOptions::default()
}

const POLY: [(f64, f64); 5] = [(0.0, 0.0), (3.0, -0.5), (4.0, 2.0), (1.5, 3.5), (-1.0, 2.0)];

/// The points of [`POLY`] and the reference list closed on the first.
fn closed_points(b: &mut B) -> Vec<u32> {
    let mut keys: Vec<u32> = POLY.iter().map(|&(x, y)| b.point(x, y)).collect();
    keys.push(keys[0]);
    keys
}

/// Chord-length parameters of a closed polygon (the tool's rule).
fn chord_params(pts: &[DVec2]) -> Vec<f64> {
    let mut acc = vec![0.0];
    for w in pts.windows(2) {
        acc.push(acc[acc.len() - 1] + (w[1] - w[0]).length());
    }
    let total = acc[acc.len() - 1];
    let mut t: Vec<f64> = acc.iter().map(|a| a / total).collect();
    *t.last_mut().unwrap() = 1.0;
    t
}

fn closed_through(b: &mut B) -> (u32, Vec<u32>, Vec<f64>) {
    let keys = closed_points(b);
    let mut pts: Vec<DVec2> = POLY.iter().map(|&(x, y)| DVec2::new(x, y)).collect();
    pts.push(pts[0]);
    let params = chord_params(&pts);
    let sp = b.entity(Geometry::InterpolatedSpline {
        through: keys.clone(),
        params: params.clone(),
        periodic: true,
    });
    (sp, keys, params)
}

fn closed_control(b: &mut B, degree: u32) -> (u32, Vec<u32>) {
    let keys = closed_points(b);
    let sp = b.entity(Geometry::ControlSpline {
        control: keys.clone(),
        degree,
        periodic: true,
    });
    (sp, keys)
}

fn spline_parts(c: ResolvedCurve) -> (usize, Vec<f64>, Vec<DVec2>) {
    match c {
        ResolvedCurve::Spline {
            degree,
            knots,
            points,
        } => (degree, knots, points),
        other => panic!("not a spline: {other:?}"),
    }
}

/// One coordinate of the periodic cubic spline interpolating `y` at `t`
/// (`t₀ = 0 < … < t_{n−1} < t_n = 1`, `y_n = y₀`) in its classical form:
/// second derivatives `M` from the cyclic system
/// `h_{i−1} M_{i−1} + 2(h_{i−1} + h_i) M_i + h_i M_{i+1} = 6(δ_i − δ_{i−1})`
/// (de Boor, *A Practical Guide to Splines*, ch. IV), solved densely, then
/// evaluated piece by piece. Independent of B-splines.
fn periodic_cubic(t: &[f64], y: &[f64], x: f64) -> f64 {
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
    // Gaussian elimination with partial pivoting.
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
    let i = (0..n).rev().find(|&i| x >= t[i]).unwrap_or(0);
    let (m0, m1) = (m[i], m[(i + 1) % n]);
    let (a_, b_) = (t[i + 1] - x, x - t[i]);
    m0 * a_.powi(3) / (6.0 * h[i])
        + m1 * b_.powi(3) / (6.0 * h[i])
        + (y[i] - m0 * h[i] * h[i] / 6.0) * a_ / h[i]
        + (y[i + 1] - m1 * h[i] * h[i] / 6.0) * b_ / h[i]
}

#[test]
fn a_closed_through_point_spline_is_the_periodic_cubic_interpolant() {
    let mut b = B::new();
    let (sp, keys, params) = closed_through(&mut b);
    let s = b.build();
    let (degree, knots, ctrl) = spline_parts(resolved(&s, sp));
    assert_eq!(degree, 3);
    assert_eq!(ctrl[0], ctrl[ctrl.len() - 1], "closed exactly");
    let ys: Vec<DVec2> = keys.iter().map(|&k| pos(&s, k)).collect();
    let xs: Vec<f64> = ys.iter().map(|p| p.x).collect();
    let yy: Vec<f64> = ys.iter().map(|p| p.y).collect();
    let mut worst = 0.0f64;
    for i in 0..=400 {
        let t = i as f64 / 400.0;
        let c = de_boor(degree, &knots, &ctrl, t);
        let e = DVec2::new(
            periodic_cubic(&params, &xs, t),
            periodic_cubic(&params, &yy, t),
        );
        worst = worst.max((c - e).length());
    }
    assert!(worst < 1e-12, "{worst}");
    for (t, q) in params.iter().zip(&ys) {
        assert!((de_boor(degree, &knots, &ctrl, *t) - *q).length() < 1e-13);
    }
}

#[test]
fn a_closed_control_spline_does_not_meet_its_control_points_and_has_no_seam() {
    for degree in [2, 3, 4] {
        let mut b = B::new();
        let (sp, keys) = closed_control(&mut b, degree);
        let s = b.build();
        let (p, knots, ctrl) = spline_parts(resolved(&s, sp));
        assert_eq!(
            (p, ctrl.len()),
            (degree as usize, POLY.len() + degree as usize)
        );
        assert_eq!(ctrl[0], ctrl[ctrl.len() - 1]);
        // Tangent continuity at the seam: the first and last legs of the
        // clamped polygon are parallel and point the same way.
        let (first, last) = (
            ctrl[1] - ctrl[0],
            ctrl[ctrl.len() - 1] - ctrl[ctrl.len() - 2],
        );
        assert!(first.perp_dot(last).abs() < 1e-12 * first.length() * last.length());
        assert!(first.dot(last) > 0.0);
        for &k in &keys[..POLY.len()] {
            assert!(dist_to_spline(&s, sp, pos(&s, k)) > 0.05, "degree {degree}");
        }
        let _ = knots;
    }
}

#[test]
fn closed_splines_have_two_freedoms_per_point() {
    let mut b = B::new();
    closed_through(&mut b);
    let s = solve(&b.build(), &opts());
    assert_eq!(s.status, SolveStatus::Solved);
    assert_eq!(s.dof, 2 * POLY.len());
    let mut b = B::new();
    closed_control(&mut b, 3);
    let s = solve(&b.build(), &opts());
    assert_eq!(s.dof, 2 * POLY.len());
}

#[test]
fn a_point_on_a_closed_spline_slides_across_the_seam() {
    for control in [false, true] {
        let mut b = B::new();
        let (sp, keys) = if control {
            closed_control(&mut b, 3)
        } else {
            let (sp, keys, _) = closed_through(&mut b);
            (sp, keys)
        };
        for &k in &keys[..POLY.len()] {
            b.fix(k);
        }
        // Start just after the seam, on the curve.
        let sk = b.build();
        let (deg, knots, ctrl) = spline_parts(resolved(&sk, sp));
        let start = de_boor(deg, &knots, &ctrl, 0.02);
        let p = b.point(start.x, start.y);
        b.c(Constraint::PointOnCurve {
            point: p,
            curve: sp,
        });
        let sk = b.build();
        // Drag it to just before the seam and back past it.
        for t in [0.97, 0.995, 0.01, 0.5] {
            let target = de_boor(deg, &knots, &ctrl, t) + DVec2::new(0.01, 0.01);
            let out = solve(
                &sk,
                &SolveOptions {
                    drag: vec![(p, target)],
                    ..SolveOptions::default()
                },
            );
            assert_eq!(out.status, SolveStatus::Solved, "control {control}, t {t}");
            assert_satisfied(&out.sketch);
            let got = pos(&out.sketch, p);
            let nearest = dist_to_spline(&out.sketch, sp, target);
            assert!(
                ((got - target).length() - nearest).abs() < 1e-9,
                "control {control}, t {t}: {} vs {nearest}",
                (got - target).length()
            );
        }
    }
}

fn invalid_with(s: &Sketch<u32>, text: &str) {
    match solve(s, &opts()).status {
        SolveStatus::Failed(SolveError::Invalid(m)) => assert!(m.contains(text), "{m}"),
        other => panic!("expected invalid ({text}), got {other:?}"),
    }
}

#[test]
fn a_closed_spline_has_no_end_to_join() {
    let mut b = B::new();
    let (sp, keys, _) = closed_through(&mut b);
    let (l, ls, _) = b.line(0.0, 0.0, -1.0, -1.0);
    b.coincident(ls, keys[0]);
    let mut tangent = b.clone();
    tangent.c(Constraint::Tangent {
        a: l,
        b: sp,
        kind: TangentKind::AtPoint { point: keys[0] },
    });
    invalid_with(
        &tangent.build(),
        "smooth closed spline and has no end points",
    );
    b.c(Constraint::Curvature {
        a: l,
        b: sp,
        point: keys[0],
    });
    invalid_with(&b.build(), "smooth closed spline and has no end points");
}

#[test]
fn malformed_closed_splines_are_invalid() {
    let mut b = B::new();
    let keys: Vec<u32> = POLY.iter().map(|&(x, y)| b.point(x, y)).collect();
    let mut open = b.clone();
    open.entity(Geometry::ControlSpline {
        control: keys.clone(),
        degree: 3,
        periodic: true,
    });
    invalid_with(&open.build(), "must end on its first point");
    let mut few = b.clone();
    few.entity(Geometry::ControlSpline {
        control: vec![keys[0], keys[1], keys[2], keys[0]],
        degree: 3,
        periodic: true,
    });
    invalid_with(&few.build(), "needs at least 4 distinct control points");
    let mut linear = b.clone();
    let mut closed = keys.clone();
    closed.push(keys[0]);
    linear.entity(Geometry::ControlSpline {
        control: closed.clone(),
        degree: 1,
        periodic: true,
    });
    invalid_with(&linear.build(), "degree 2 to 5");
    let mut two = b.clone();
    two.entity(Geometry::InterpolatedSpline {
        through: vec![keys[0], keys[1], keys[0]],
        params: vec![0.0, 0.5, 1.0],
        periodic: true,
    });
    invalid_with(&two.build(), "at least 3 distinct points");
}

#[test]
fn a_mirror_image_of_a_closed_spline_is_closed() {
    let mut b = B::new();
    let (sp, keys, params) = closed_through(&mut b);
    let (axis, ..) = b.line(6.0, -5.0, 6.0, 5.0);
    let mut image_keys: Vec<u32> = POLY.iter().map(|&(x, y)| b.point(12.0 - x, y)).collect();
    image_keys.push(image_keys[0]);
    let make = |b: &mut B, periodic: bool| {
        b.entity(Geometry::InterpolatedSpline {
            through: image_keys.clone(),
            params: params.clone(),
            periodic,
        })
    };
    let mut wrong = b.clone();
    let open_image = make(&mut wrong, false);
    wrong.c(Constraint::Mirror {
        axis,
        pairs: vec![(sp, open_image)],
    });
    invalid_with(&wrong.build(), "not splines of the same form");
    let image = make(&mut b, true);
    b.c(Constraint::Mirror {
        axis,
        pairs: vec![(sp, image)],
    });
    // Move the source: the image follows exactly.
    let out = solve(
        &b.build(),
        &SolveOptions {
            drag: vec![(keys[2], DVec2::new(4.5, 2.5))],
            ..SolveOptions::default()
        },
    );
    assert_eq!(out.status, SolveStatus::Solved);
    assert_satisfied(&out.sketch);
    let q = pos(&out.sketch, image_keys[2]);
    assert!((q - DVec2::new(12.0 - 4.5, 2.5)).length() < 1e-9, "{q}");
}

#[test]
fn documents_without_the_flag_read_as_open_and_stay_byte_identical() {
    let mut b = B::new();
    let keys = closed_points(&mut b);
    b.entity(Geometry::ControlSpline {
        control: keys.clone(),
        degree: 3,
        periodic: false,
    });
    let text = serde_json::to_string(&b.build()).unwrap();
    assert!(!text.contains("periodic"), "a false flag is not written");
    let back: Sketch<u32> = serde_json::from_str(&text).unwrap();
    assert!(matches!(
        back.entities.last().unwrap().geometry,
        Geometry::ControlSpline {
            periodic: false,
            ..
        }
    ));
    let mut b = B::new();
    closed_through(&mut b);
    let s = b.build();
    let text = serde_json::to_string(&s).unwrap();
    assert!(text.contains("\"periodic\":true"));
    let back: Sketch<u32> = serde_json::from_str(&text).unwrap();
    assert_eq!(back, s);
}
