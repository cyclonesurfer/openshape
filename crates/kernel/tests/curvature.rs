//! Curvature samples for curvature combs, checked against closed forms
//! computed here (circles, an ellipse, a parabola) and against finite
//! differences of the curve's own points (a spline whose curvature jumps).

use std::f64::consts::{PI, TAU};

use openshape_kernel::{CurvatureSample, DVec2, Kernel, KernelError, OsKernel, ProfileCurve};

fn k() -> OsKernel {
    OsKernel::default()
}

fn v(x: f64, y: f64) -> DVec2 {
    DVec2::new(x, y)
}

/// Lengths of the chords between consecutive samples.
fn chords(s: &[CurvatureSample]) -> Vec<f64> {
    s.windows(2)
        .map(|w| (w[1].point - w[0].point).length())
        .collect()
}

#[test]
fn a_circle_has_its_reciprocal_radius_everywhere() {
    let (c, r) = (v(0.02, -0.01), 0.005);
    let s = k()
        .curvature_samples(
            &ProfileCurve::Circle {
                center: c,
                radius: r,
            },
            40,
        )
        .unwrap();
    assert_eq!(s.len(), 41);
    for p in &s {
        assert!(((p.point - c).length() - r).abs() < 1e-15, "{p:?}");
        assert!((p.curvature * r - 1.0).abs() < 1e-12, "{p:?}");
        // Counter-clockwise: the tangent is the radius turned left.
        let e = (p.point - c) / r;
        assert!((p.tangent - e.perp()).length() < 1e-12, "{p:?}");
    }
    assert!((s[0].point - s[40].point).length() < 1e-15);
    for l in chords(&s) {
        let want = 2.0 * r * (PI / 40.0).sin();
        assert!((l / want - 1.0).abs() < 1e-9, "{l} {want}");
    }
}

#[test]
fn a_clockwise_arc_turns_right() {
    let r = 0.03;
    let s = k()
        .curvature_samples(
            &ProfileCurve::Arc {
                center: v(0.0, 0.0),
                radius: r,
                start_angle: 1.0,
                sweep_angle: -2.0,
            },
            8,
        )
        .unwrap();
    assert_eq!(s.len(), 9);
    for p in &s {
        assert!((p.curvature * r + 1.0).abs() < 1e-12, "{p:?}");
    }
    assert!((s[0].point - v(1f64.cos(), 1f64.sin()) * r).length() < 1e-15);
    assert!((s[8].point - v((-1f64).cos(), (-1f64).sin()) * r).length() < 1e-15);
}

#[test]
fn a_line_is_straight_and_evenly_cut() {
    let s = k()
        .curvature_samples(
            &ProfileCurve::Line {
                start: v(0.0, 0.0),
                end: v(0.03, 0.04),
            },
            5,
        )
        .unwrap();
    assert_eq!(s.len(), 6);
    for (i, p) in s.iter().enumerate() {
        assert_eq!(p.curvature, 0.0);
        assert!((p.tangent - v(0.6, 0.8)).length() < 1e-15);
        assert!((p.point - v(0.006, 0.008) * i as f64).length() < 1e-15);
    }
}

#[test]
fn an_ellipse_matches_the_closed_form() {
    // κ(θ) = ab / (a² sin²θ + b² cos²θ)^{3/2} at the point (a cos θ, b sin θ).
    let (a, b) = (0.05, 0.02);
    let c = v(0.01, 0.02);
    let e = ProfileCurve::elliptical_arc(c, v(a, 0.0), v(0.0, b), 0.3, TAU).unwrap();
    let s = k().curvature_samples(&e, 200).unwrap();
    assert_eq!(s.len(), 201, "no jumps at the knots of an exact ellipse");
    let mut worst: f64 = 0.0;
    for p in &s {
        let q = p.point - c;
        let th = (q.y / b).atan2(q.x / a);
        let (sn, cs) = th.sin_cos();
        let want = a * b / (a * a * sn * sn + b * b * cs * cs).powf(1.5);
        worst = worst.max((p.curvature / want - 1.0).abs());
        assert!(((q.x / a).powi(2) + (q.y / b).powi(2) - 1.0).abs() < 1e-12);
    }
    assert!(worst < 1e-11, "{worst}");
    // Evenly spaced by length: the arc lengths between samples (Simpson's
    // rule on the closed-form speed √(a² sin²θ + b² cos²θ)) agree to 1e-6.
    let anomaly = |p: DVec2| {
        let q = p - c;
        (q.y / b).atan2(q.x / a)
    };
    let mut lens = Vec::new();
    for w in s.windows(2) {
        let t0 = anomaly(w[0].point);
        let t1 = t0 + (anomaly(w[1].point) - t0).rem_euclid(TAU);
        let speed = |t: f64| (a * a * t.sin().powi(2) + b * b * t.cos().powi(2)).sqrt();
        let m = 64;
        let h = (t1 - t0) / m as f64;
        let sum: f64 = (0..=m)
            .map(|i| {
                let wgt = if i == 0 || i == m {
                    1.0
                } else if i % 2 == 1 {
                    4.0
                } else {
                    2.0
                };
                wgt * speed(t0 + h * i as f64)
            })
            .sum();
        lens.push(sum * h / 3.0);
    }
    let mean = lens.iter().sum::<f64>() / lens.len() as f64;
    for x in lens {
        assert!((x / mean - 1.0).abs() < 1e-6, "{x} {mean}");
    }
}

#[test]
fn a_parabola_matches_the_closed_form() {
    // The quadratic Bézier (-1, 1), (0, -1), (1, 1) is y = x² (scaled to
    // centimetres); κ = 2s / (1 + 4x²)^{3/2} in the unscaled units.
    let s0 = 0.01;
    let p =
        ProfileCurve::conic(v(-1.0, 1.0) * s0, v(0.0, -1.0) * s0, v(1.0, 1.0) * s0, 1.0).unwrap();
    let s = k().curvature_samples(&p, 64).unwrap();
    assert_eq!(s.len(), 65);
    for q in &s {
        let x = q.point.x / s0;
        assert!((q.point.y / s0 - x * x).abs() < 1e-12, "{q:?}");
        let want = 2.0 / (1.0 + 4.0 * x * x).powf(1.5) / s0;
        assert!((q.curvature / want - 1.0).abs() < 1e-11, "{q:?} {want}");
    }
}

#[test]
fn a_quadratic_spline_shows_both_sides_of_its_curvature_jumps() {
    // A uniform quadratic B-spline is only C1: its curvature jumps at the
    // interior knots 1 and 2.
    let sp = ProfileCurve::Spline {
        degree: 2,
        knots: vec![0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 3.0, 3.0],
        points: vec![
            v(0.0, 0.0),
            v(0.01, 0.0),
            v(0.02, 0.01),
            v(0.03, 0.0),
            v(0.04, 0.02),
        ],
        weights: None,
    };
    let kern = k();
    let s = kern.curvature_samples(&sp, 30).unwrap();
    assert_eq!(s.len(), 31 + 4);
    for w in s.windows(2) {
        assert!(w[0].t <= w[1].t);
    }
    // Curvature from differences of the curve's own points, with a stencil
    // inside one span (the span right of `t` when `right`): exact for the
    // quadratic polynomial pieces up to rounding.
    let fd = |t: f64, right: bool| {
        let span = if right { t.floor() } else { t.ceil() - 1.0 };
        let h = 0.2;
        let p = |u: f64| kern.curve_point(&sp, u).unwrap();
        let (d1, d2) = if t - 2.0 * h >= span {
            let (a, b, c) = (p(t - 2.0 * h), p(t - h), p(t));
            (
                (3.0 * c - 4.0 * b + a) / (2.0 * h),
                (c - 2.0 * b + a) / (h * h),
            )
        } else {
            let (a, b, c) = (p(t), p(t + h), p(t + 2.0 * h));
            (
                (-3.0 * a + 4.0 * b - c) / (2.0 * h),
                (a - 2.0 * b + c) / (h * h),
            )
        };
        d1.perp_dot(d2) / d1.length().powi(3)
    };
    for knot in [1.0, 2.0] {
        let at: Vec<&CurvatureSample> = s.iter().filter(|q| q.t == knot).collect();
        assert_eq!(at.len(), 2, "both sides of {knot}");
        let (left, right) = (fd(knot, false), fd(knot, true));
        assert!((at[0].curvature / left - 1.0).abs() < 1e-9, "{at:?} {left}");
        assert!(
            (at[1].curvature / right - 1.0).abs() < 1e-9,
            "{at:?} {right}"
        );
        assert!((left - right).abs() > 0.1 * left.abs().max(right.abs()));
    }
    for q in &s {
        let right = q.t < 3.0 && !(q.t == 1.0 || q.t == 2.0);
        let want = fd(q.t, right || q.t == 0.0);
        if q.t == 1.0 || q.t == 2.0 {
            continue;
        }
        assert!(
            (q.curvature - want).abs() < 1e-9 * want.abs().max(1.0),
            "{q:?} {want}"
        );
    }
}

#[test]
fn degenerate_curves_are_refused() {
    let e = k()
        .curvature_samples(
            &ProfileCurve::Line {
                start: v(0.0, 0.0),
                end: v(0.0, 0.0),
            },
            10,
        )
        .unwrap_err();
    assert!(matches!(e, KernelError::InvalidInput(_)), "{e:?}");
}
