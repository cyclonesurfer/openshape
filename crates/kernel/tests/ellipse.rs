//! Ellipses and elliptical arcs as exact rational quadratic profile curves
//! (design note 0008 §6): every point on the ellipse to rounding, ends and
//! direction as asked, and exact areas, extrusions and spheroids of
//! revolution against closed forms.

use std::f64::consts::{FRAC_PI_2, PI, TAU};

use openshape_kernel::{
    BodyMesh, DVec2, ExtrudeSpec, Frame, Kernel, OsKernel, Profile, ProfileCurve, RevolveSpec,
    SurfaceKind, TessellationOptions,
};
use proptest::prelude::*;

fn k() -> OsKernel {
    OsKernel::default()
}

/// `|q| − 1` with `p = c + [a b]·q`: zero exactly on the ellipse, for any
/// two conjugate semi-diameters.
fn off_ellipse(c: DVec2, a: DVec2, b: DVec2, p: DVec2) -> f64 {
    let w = p - c;
    let det = a.perp_dot(b);
    DVec2::new(w.perp_dot(b) / det, a.perp_dot(w) / det).length() - 1.0
}

fn at(c: DVec2, a: DVec2, b: DVec2, t: f64) -> DVec2 {
    c + a * t.cos() + b * t.sin()
}

/// The eccentric anomaly of a point of the ellipse.
fn anomaly(c: DVec2, a: DVec2, b: DVec2, p: DVec2) -> f64 {
    let w = p - c;
    let det = a.perp_dot(b);
    (a.perp_dot(w) / det).atan2(w.perp_dot(b) / det)
}

fn line(a: DVec2, b: DVec2) -> ProfileCurve {
    ProfileCurve::Line { start: a, end: b }
}

/// Checks an elliptical arc: every sample on the ellipse, the ends where
/// asked, and the anomaly moving monotonically by `sweep`.
fn check_arc(c: DVec2, a: DVec2, b: DVec2, start: f64, sweep: f64) {
    let curve = ProfileCurve::elliptical_arc(c, a, b, start, sweep).unwrap();
    let size = c.length() + a.length() + b.length();
    let pts = curve.sample(2000).unwrap();
    for p in &pts {
        let off = off_ellipse(c, a, b, *p);
        assert!(
            off.abs() < 1e-14 * size / a.length().min(b.length()),
            "off by {off:e}"
        );
    }
    assert!((pts[0] - at(c, a, b, start)).length() < 1e-15 * size);
    let last = *pts.last().unwrap();
    if sweep.abs() == TAU {
        assert_eq!(last, pts[0], "a full turn closes exactly");
    } else {
        assert!((last - at(c, a, b, start + sweep)).length() < 1e-14 * size);
    }
    // The anomaly runs from `start` by `sweep`, without turning back.
    let mut total = 0.0;
    let mut prev = anomaly(c, a, b, pts[0]);
    for p in &pts[1..] {
        let t = anomaly(c, a, b, *p);
        let step = (t - prev + PI).rem_euclid(TAU) - PI;
        assert!(step * sweep.signum() >= 0.0, "turned back");
        total += step;
        prev = t;
    }
    // Rounded coordinates fix the anomaly only to about ε·size/(shorter axis).
    let tol = 1e-14 * size / a.length().min(b.length());
    assert!((total - sweep).abs() < tol, "{total} vs {sweep}");
}

#[test]
fn every_point_is_on_the_ellipse_and_the_ends_are_exact() {
    let c = DVec2::new(0.3, -1.2);
    let frames = [
        (DVec2::new(2.0, 0.0), DVec2::new(0.0, 0.5)),
        (DVec2::new(0.0, 1.0), DVec2::new(3.0, 0.0)),
        // Rotated, and left-handed (the minor axis right of the major).
        (DVec2::new(1.2, 0.9), DVec2::new(0.3, -0.4)),
        // Conjugate semi-diameters that are not perpendicular.
        (DVec2::new(1.0, 0.2), DVec2::new(-0.3, 0.7)),
        // Very flat.
        (DVec2::new(1e-3, 0.0), DVec2::new(0.0, 1e-6)),
    ];
    for (a, b) in frames {
        for start in [0.0, -2.5, 1.1, FRAC_PI_2] {
            for sweep in [1e-3, 0.7, FRAC_PI_2, PI, 4.5, TAU] {
                check_arc(c, a, b, start, sweep);
                check_arc(c, a, b, start, -sweep);
            }
        }
    }
}

#[test]
fn segments_are_at_most_a_quarter_of_the_anomaly() {
    for (sweep, segments) in [
        (0.3, 1),
        (FRAC_PI_2, 1),
        (1.6, 2),
        (PI, 2),
        (4.0, 3),
        (TAU, 4),
    ] {
        let ProfileCurve::Spline {
            degree,
            knots,
            points,
            weights: Some(w),
        } = ProfileCurve::elliptical_arc(DVec2::ZERO, DVec2::X, DVec2::Y * 0.5, 0.0, sweep)
            .unwrap()
        else {
            panic!("not a rational spline")
        };
        assert_eq!(degree, 2);
        assert_eq!(points.len(), 2 * segments + 1);
        assert_eq!(knots.len(), points.len() + 3);
        let mid = (0.5 * sweep / segments as f64).cos();
        for (i, wi) in w.iter().enumerate() {
            assert_eq!(*wi, if i % 2 == 1 { mid } else { 1.0 });
        }
    }
}

#[test]
fn degenerate_ellipses_are_refused() {
    let c = DVec2::ZERO;
    let bad = [
        (DVec2::X, DVec2::X * 2.0, 0.0, 1.0),
        (DVec2::ZERO, DVec2::Y, 0.0, 1.0),
        (DVec2::X, DVec2::Y, 0.0, 0.0),
        (DVec2::X, DVec2::Y, 0.0, TAU + 1e-9),
        (DVec2::X, DVec2::Y, f64::NAN, 1.0),
        (DVec2::X, DVec2::Y, 0.0, f64::INFINITY),
        (DVec2::new(f64::NAN, 0.0), DVec2::Y, 0.0, 1.0),
    ];
    for (a, b, start, sweep) in bad {
        assert!(
            ProfileCurve::elliptical_arc(c, a, b, start, sweep).is_err(),
            "{a} {b} {start} {sweep}"
        );
    }
}

fn region_areas(curves: &[ProfileCurve]) -> Vec<f64> {
    let a = k().sketch_regions(curves).unwrap();
    let mut out: Vec<f64> = a.regions.iter().map(|r| r.area).collect();
    out.sort_by(f64::total_cmp);
    out
}

fn assert_rel(a: f64, b: f64, rel: f64) {
    assert!((a - b).abs() <= rel * b.abs(), "{a} vs {b}");
}

#[test]
fn a_full_ellipse_bounds_pi_a_b() {
    for (c, a, b) in [
        (
            DVec2::new(1.0, 2.0),
            DVec2::new(3.0, 0.0),
            DVec2::new(0.0, 1.0),
        ),
        (
            DVec2::new(-0.2, 0.1),
            DVec2::new(0.6, 0.8),
            DVec2::new(0.16, -0.12),
        ),
        (DVec2::ZERO, DVec2::new(1e-3, 0.0), DVec2::new(0.0, 2e-4)),
    ] {
        for sweep in [TAU, -TAU] {
            let e = ProfileCurve::elliptical_arc(c, a, b, 0.4, sweep).unwrap();
            let areas = region_areas(&[e]);
            assert_eq!(areas.len(), 1);
            assert_rel(areas[0], PI * a.length() * b.length(), 1e-12);
        }
    }
}

#[test]
fn sectors_and_segments_have_their_exact_areas() {
    let c = DVec2::new(0.5, -0.25);
    let (a, b) = (DVec2::new(1.6, 1.2), DVec2::new(-0.3, 0.4));
    let ab = a.length() * b.length();
    for t0 in [0.0, 2.0, -1.3] {
        for d in [0.7, 2.0, PI, 4.5] {
            let (p0, p1) = (at(c, a, b, t0), at(c, a, b, t0 + d));
            let arc = ProfileCurve::elliptical_arc(c, a, b, t0, d).unwrap();
            // The sector: the affine image of a circular sector of angle d.
            let areas = region_areas(&[arc.clone(), line(p1, c), line(c, p0)]);
            assert_eq!(areas.len(), 1);
            assert_rel(areas[0], 0.5 * ab * d, 1e-12);
            // The segment cut off by the chord.
            let areas = region_areas(&[arc, line(p1, p0)]);
            assert_eq!(areas.len(), 1);
            assert_rel(areas[0], 0.5 * ab * (d - d.sin()), 1e-12);
        }
    }
}

/// The perimeter by the trapezoidal rule, which converges geometrically for
/// a smooth periodic integrand (Trefethen & Weideman, SIAM Review 56, 2014).
fn perimeter(a: f64, b: f64) -> f64 {
    let n = 4096;
    let h = TAU / n as f64;
    (0..n)
        .map(|i| {
            let t = i as f64 * h;
            (a * a * t.sin().powi(2) + b * b * t.cos().powi(2)).sqrt()
        })
        .sum::<f64>()
        * h
}

fn assert_watertight(mesh: &BodyMesh) {
    let (_, tris) = mesh.welded();
    let mut directed = std::collections::HashMap::new();
    for t in &tris {
        for (p, q) in [(t[0], t[1]), (t[1], t[2]), (t[2], t[0])] {
            assert_ne!(p, q, "degenerate triangle");
            *directed.entry((p, q)).or_insert(0) += 1;
        }
    }
    for (&(p, q), &n) in &directed {
        assert_eq!(n, 1, "directed edge used {n} times");
        assert_eq!(directed.get(&(q, p)), Some(&1), "mesh not closed");
    }
}

fn check_body(body: &openshape_kernel::Body, volume: f64, area: f64) {
    let kern = k();
    let issues = kern.validate(body);
    assert!(issues.is_empty(), "invalid body: {issues:?}");
    let t = kern.topology(body);
    assert_eq!(t.genus, Some(0));
    assert_eq!(t.solids, 1);
    let m = kern.mass_properties(body, 1e-10).unwrap();
    assert_rel(m.volume, volume, 1e-9);
    assert_rel(m.area, area, 1e-9);
    for opts in [
        TessellationOptions::default(),
        TessellationOptions {
            chord_tolerance: 1e-3,
            max_angle: 0.5,
        },
    ] {
        assert_watertight(&kern.tessellate(body, &opts).unwrap());
    }
}

#[test]
fn an_extruded_ellipse_is_an_elliptic_cylinder() {
    let (a, b, h) = (0.8, 0.3, 1.5);
    let u = DVec2::from_angle(0.6);
    let e =
        ProfileCurve::elliptical_arc(DVec2::new(0.2, 0.1), u * a, u.perp() * b, 0.0, TAU).unwrap();
    let body = k()
        .extrude(
            &Profile {
                plane: Frame::WORLD,
                curves: vec![e],
            },
            &ExtrudeSpec {
                distance: h,
                symmetric: false,
            },
        )
        .unwrap();
    check_body(
        &body,
        PI * a * b * h,
        2.0 * PI * a * b + perimeter(a, b) * h,
    );
    assert!(k().face_kinds(&body).contains(&SurfaceKind::Extrusion));
}

/// The half ellipse with semi-axis `polar` along the y axis and `equatorial`
/// along x, closed by the axis, revolved a full turn about the y axis.
fn spheroid(equatorial: f64, polar: f64) -> openshape_kernel::Body {
    let arc = ProfileCurve::elliptical_arc(
        DVec2::ZERO,
        DVec2::new(equatorial, 0.0),
        DVec2::new(0.0, polar),
        -FRAC_PI_2,
        PI,
    )
    .unwrap();
    k().revolve(
        &Profile {
            plane: Frame::front(),
            curves: vec![arc, line(DVec2::new(0.0, polar), DVec2::new(0.0, -polar))],
        },
        &RevolveSpec {
            axis_origin: DVec2::ZERO,
            axis_dir: DVec2::Y,
            angle: TAU,
        },
    )
    .unwrap()
}

#[test]
fn revolved_half_ellipses_are_spheroids() {
    // Prolate (polar semi-axis c > equatorial a): S = 2πa²(1 + c·asin(e)/(a·e)),
    // e² = 1 − a²/c².
    let (a, c): (f64, f64) = (0.3, 0.5);
    let e = (1.0 - a * a / (c * c)).sqrt();
    let area = TAU * a * a * (1.0 + c * e.asin() / (a * e));
    let body = spheroid(a, c);
    check_body(&body, 4.0 / 3.0 * PI * a * a * c, area);
    assert!(k().face_kinds(&body).contains(&SurfaceKind::Revolution));
    // Oblate (polar c < equatorial a): S = 2πa² + π(c²/e)·ln((1 + e)/(1 − e)),
    // e² = 1 − c²/a².
    let (a, c): (f64, f64) = (0.5, 0.3);
    let e = (1.0 - c * c / (a * a)).sqrt();
    let area = TAU * a * a + PI * c * c / e * ((1.0 + e) / (1.0 - e)).ln();
    check_body(&spheroid(a, c), 4.0 / 3.0 * PI * a * a * c, area);
}

/// A full ellipse revolved about an axis that misses it: an elliptic torus,
/// by Pappus's theorems V = 2π·d·πab and S = 2π·d·perimeter, d the
/// centre's distance from the axis. The surface's profile is a closed
/// curve, so this also covers boundary points over its closure.
#[test]
fn a_revolved_ellipse_is_an_elliptic_torus() {
    let (a, b, d) = (0.2, 0.12, 0.7);
    for (major, minor) in [
        (DVec2::new(a, 0.0), DVec2::new(0.0, b)),
        (DVec2::new(0.0, a), DVec2::new(-b, 0.0)),
    ] {
        let e = ProfileCurve::elliptical_arc(DVec2::new(d, 0.1), major, minor, 0.3, TAU).unwrap();
        let body = k()
            .revolve(
                &Profile {
                    plane: Frame::front(),
                    curves: vec![e],
                },
                &RevolveSpec {
                    axis_origin: DVec2::ZERO,
                    axis_dir: DVec2::Y,
                    angle: TAU,
                },
            )
            .unwrap();
        let kern = k();
        assert!(kern.validate(&body).is_empty());
        assert_eq!(kern.topology(&body).genus, Some(1));
        let m = kern.mass_properties(&body, 1e-10).unwrap();
        assert_rel(m.volume, TAU * d * PI * a * b, 1e-9);
        assert_rel(m.area, TAU * d * perimeter(a, b), 1e-9);
        for opts in [
            TessellationOptions::default(),
            TessellationOptions {
                chord_tolerance: 1e-3,
                max_angle: 0.5,
            },
        ] {
            assert_watertight(&kern.tessellate(&body, &opts).unwrap());
        }
    }
}

proptest! {
    /// Extrusions of closed splines (an ellipse; a closed control-point
    /// spline) mesh watertight at any height. Their side faces are swept
    /// from a closed curve, whose closure point projects to either end of
    /// the curve's domain; before the tessellator put such points on the
    /// face's own end, some heights failed ("face boundary intersects
    /// itself in parameter space").
    #[test]
    fn closed_spline_extrusions_mesh_at_any_height(
        h in 0.01f64..10.0, cx in -2.0f64..2.0, rot in 0.0f64..TAU, start in -PI..PI,
    ) {
        let u = DVec2::from_angle(rot);
        let c = DVec2::new(cx, 0.3);
        let ellipse = ProfileCurve::elliptical_arc(c, u * 0.8, u.perp() * 0.3, start, TAU).unwrap();
        let ctrl: Vec<DVec2> = [(0.0, 0.0), (1.0, 0.2), (1.2, 1.0), (0.3, 1.1), (0.0, 0.0)]
            .iter()
            .map(|&(x, y)| c + DVec2::new(x, y).rotate(u))
            .collect();
        let spline = ProfileCurve::Spline {
            degree: 3,
            knots: vec![0.0, 0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0, 1.0],
            points: ctrl,
            weights: None,
        };
        for curve in [ellipse, spline] {
            let body = k()
                .extrude(
                    &Profile { plane: Frame::WORLD, curves: vec![curve] },
                    &ExtrudeSpec { distance: h, symmetric: false },
                )
                .unwrap();
            prop_assert!(k().validate(&body).is_empty());
            for opts in [
                TessellationOptions::default(),
                TessellationOptions { chord_tolerance: 1e-3, max_angle: 0.5 },
            ] {
                assert_watertight(&k().tessellate(&body, &opts).unwrap());
            }
        }
    }

    /// Random ellipses and ranges: every sample on the ellipse, the ends
    /// exact, the anomaly monotone.
    #[test]
    fn random_arcs_are_exact(
        cx in -5.0f64..5.0, cy in -5.0f64..5.0,
        la in 0.01f64..3.0, ratio in 0.02f64..1.0, rot in 0.0f64..TAU,
        flip in any::<bool>(), start in -PI..PI, sweep in 1e-4f64..TAU,
    ) {
        let u = DVec2::from_angle(rot);
        let b = if flip { -u.perp() } else { u.perp() } * la * ratio;
        check_arc(DVec2::new(cx, cy), u * la, b, start, sweep);
        check_arc(DVec2::new(cx, cy), u * la, b, start, -sweep);
    }
}
