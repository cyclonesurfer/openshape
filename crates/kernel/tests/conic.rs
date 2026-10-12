//! Conic arcs as exact rational quadratic profile curves (design note 0009
//! §6, §9): every point on the conic (Lee's implicit equation), the ends
//! exact, segment areas against their closed form and an independent
//! quadrature, a quarter ellipse, a paraboloid and a hyperboloid of
//! revolution against closed forms.

use std::f64::consts::PI;

use openshape_kernel::{
    BodyMesh, DVec2, ExtrudeSpec, Frame, Kernel, OsKernel, Profile, ProfileCurve, RevolveSpec,
    TessellationOptions,
};
use proptest::prelude::*;

fn k() -> OsKernel {
    OsKernel::default()
}

fn line(a: DVec2, b: DVec2) -> ProfileCurve {
    ProfileCurve::Line { start: a, end: b }
}

/// Lee's implicit equation in the barycentric coordinates of the control
/// triangle, `τ₁² − 4w²τ₀τ₂`, and the smallest coordinate (the arc is
/// where every one is non-negative).
fn implicit(p: [DVec2; 3], w: f64, x: DVec2) -> (f64, f64) {
    let area = (p[1] - p[0]).perp_dot(p[2] - p[0]);
    let t0 = (p[1] - x).perp_dot(p[2] - x) / area;
    let t1 = (p[2] - x).perp_dot(p[0] - x) / area;
    let t2 = (p[0] - x).perp_dot(p[1] - x) / area;
    (t1 * t1 - 4.0 * w * w * t0 * t2, t0.min(t1).min(t2))
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

/// The area between a conic arc and its chord as a fraction of its control
/// triangle's (design note 0009 §9), from the circular segment of
/// half-angle α with w = cos α, continued analytically to w ≥ 1. `1 − w²`
/// is formed as `(1 − w)(1 + w)`, whose first factor is exact, since near
/// w = 1 the numerator cancels to O(α³) and the rounding of `w²` would
/// dominate it.
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

/// `½∫₀¹ (C − P₀) × C′ dt` by composite 8-point Gauss–Legendre on 2000
/// panels, with `C = N/D` evaluated by the homogeneous de Casteljau
/// algorithm (P&T §1.4) and `C′ = (N′ − C D′)/D` from the hodograph of the
/// homogeneous control points (P&T §4.5): independent of the kernel and of
/// the closed form.
fn quadrature(p: [DVec2; 3], w: f64) -> f64 {
    // Homogeneous control points (wᵢPᵢ, wᵢ).
    let h = [(p[0], 1.0), (p[1] * w, w), (p[2], 1.0)];
    let c = |t: f64| {
        let lerp = |a: (DVec2, f64), b: (DVec2, f64)| {
            (a.0 * (1.0 - t) + b.0 * t, a.1 * (1.0 - t) + b.1 * t)
        };
        let (q0, q1) = (lerp(h[0], h[1]), lerp(h[1], h[2]));
        let r = lerp(q0, q1);
        // The hodograph of the homogeneous curve: 2((1 − t)Δ₀ + tΔ₁).
        let dn = ((h[1].0 - h[0].0) * (1.0 - t) + (h[2].0 - h[1].0) * t) * 2.0;
        let dd = ((h[1].1 - h[0].1) * (1.0 - t) + (h[2].1 - h[1].1) * t) * 2.0;
        let point = r.0 / r.1;
        (point, (dn - point * dd) / r.1)
    };
    const X: [f64; 4] = [
        0.183_434_642_495_649_8,
        0.525_532_409_916_329,
        0.796_666_477_413_626_7,
        0.960_289_856_497_536_3,
    ];
    const W: [f64; 4] = [
        0.362_683_783_378_362,
        0.313_706_645_877_887_3,
        0.222_381_034_453_374_5,
        0.101_228_536_290_376_3,
    ];
    let n = 2000;
    let mut sum = 0.0;
    for i in 0..n {
        let (a, b) = (i as f64 / n as f64, (i + 1) as f64 / n as f64);
        let (m, r) = (0.5 * (a + b), 0.5 * (b - a));
        for (x, wt) in X.iter().zip(W) {
            for t in [m - r * x, m + r * x] {
                let (point, d) = c(t);
                sum += wt * r * (point - p[0]).perp_dot(d);
            }
        }
    }
    0.5 * sum.abs()
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
fn every_point_is_on_the_conic_and_the_ends_are_exact() {
    let p = [
        DVec2::new(0.2, -0.1),
        DVec2::new(0.9, 1.4),
        DVec2::new(2.0, 0.3),
    ];
    for w in [1e-6, 0.05, 0.5, 1.0, 2.0, 40.0, 1e6] {
        let curve = ProfileCurve::conic(p[0], p[1], p[2], w).unwrap();
        let pts = curve.sample(2000).unwrap();
        assert_eq!(pts[0], p[0]);
        assert_eq!(*pts.last().unwrap(), p[2]);
        for x in &pts {
            let (f, inside) = implicit(p, w, *x);
            assert!(f.abs() < 1e-13 * (1.0 + w * w), "{w}: {f:e}");
            assert!(inside > -1e-14, "{w}: outside the triangle");
        }
    }
}

#[test]
fn degenerate_conics_are_refused() {
    let (a, b, c) = (DVec2::ZERO, DVec2::new(1.0, 1.0), DVec2::new(2.0, 0.0));
    for (p0, p1, p2, w) in [
        (a, b, c, 0.0),
        (a, b, c, -1.0),
        (a, b, c, f64::NAN),
        (a, b, c, f64::INFINITY),
        (a, DVec2::new(1.0, 0.0), c, 1.0),
        (a, b, a, 1.0),
        (a, a, c, 1.0),
        (DVec2::new(f64::NAN, 0.0), b, c, 1.0),
    ] {
        assert!(
            ProfileCurve::conic(p0, p1, p2, w).is_err(),
            "{p0} {p1} {p2} {w}"
        );
    }
}

#[test]
fn the_segment_area_formula_matches_a_quadrature() {
    let p = [
        DVec2::new(0.0, 0.0),
        DVec2::new(0.4, 1.0),
        DVec2::new(1.0, 0.2),
    ];
    let tri = 0.5 * (p[1] - p[0]).perp_dot(p[2] - p[0]).abs();
    for w in [0.01, 0.2, 0.5, 0.9, 0.999, 1.0, 1.001, 1.5, 4.0, 25.0] {
        assert_rel(quadrature(p, w), ratio(w) * tri, 1e-10);
    }
    // Continuity through the parabola: both branches follow the series
    // r = 2/3 − 2α²/15 + O(α⁴) about w = cos α = 1, with α² = −arcosh² w
    // beyond it.
    for (w, a2) in [
        (1.0 - 1e-5, (1.0 - 1e-5f64).acos().powi(2)),
        (1.0 + 1e-5, -acosh(1.0 + 1e-5).powi(2)),
    ] {
        assert_rel(ratio(w), 2.0 / 3.0 - 2.0 * a2 / 15.0, 1e-10);
    }
    let alpha: f64 = 0.7;
    let circle = alpha - alpha.sin() * alpha.cos();
    let triangle = alpha.sin().powi(3) / alpha.cos();
    assert_rel(ratio(alpha.cos()) * triangle, circle, 1e-14);
}

#[test]
fn segments_bound_their_exact_areas() {
    for (p, rel) in [
        (
            [
                DVec2::new(0.0, 0.0),
                DVec2::new(0.4, 1.0),
                DVec2::new(1.0, 0.2),
            ],
            1e-12,
        ),
        (
            [
                DVec2::new(3.0, -2.0),
                DVec2::new(1.0, 0.5),
                DVec2::new(-1.0, -1.0),
            ],
            1e-12,
        ),
        (
            [
                DVec2::new(0.0, 0.0),
                DVec2::new(2e-3, 1e-3),
                DVec2::new(1e-3, 3e-3),
            ],
            1e-12,
        ),
    ] {
        let tri = 0.5 * (p[1] - p[0]).perp_dot(p[2] - p[0]).abs();
        for w in [0.05, 0.3, 1.0 / 2f64.sqrt(), 1.0, 1.7, 6.0] {
            let arc = ProfileCurve::conic(p[0], p[1], p[2], w).unwrap();
            let areas = region_areas(&[arc, line(p[2], p[0])]);
            assert_eq!(areas.len(), 1);
            assert_rel(areas[0], ratio(w) * tri, rel);
        }
    }
}

#[test]
fn a_quarter_ellipse_as_a_conic_bounds_a_quarter_of_pi_a_b() {
    let (a, b) = (1.5, 0.6);
    let rho = 1.0 / (1.0 + 2f64.sqrt());
    let w = rho / (1.0 - rho);
    assert_rel(w, 0.5f64.sqrt(), 1e-15);
    let (pa, pb, o) = (DVec2::new(a, 0.0), DVec2::new(0.0, b), DVec2::ZERO);
    let arc = ProfileCurve::conic(pa, DVec2::new(a, b), pb, w).unwrap();
    // Every point on the ellipse x²/a² + y²/b² = 1.
    for q in arc.sample(500).unwrap() {
        assert!((q.x * q.x / (a * a) + q.y * q.y / (b * b) - 1.0).abs() < 1e-14);
    }
    let areas = region_areas(&[arc, line(pb, o), line(o, pa)]);
    assert_eq!(areas.len(), 1);
    assert_rel(areas[0], PI * a * b / 4.0, 1e-12);
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

fn revolve(curves: Vec<ProfileCurve>) -> openshape_kernel::Body {
    k().revolve(
        &Profile {
            plane: Frame::front(),
            curves,
        },
        &RevolveSpec {
            axis_origin: DVec2::ZERO,
            axis_dir: DVec2::Y,
            angle: 2.0 * PI,
        },
    )
    .unwrap()
}

#[test]
fn a_revolved_parabolic_arc_is_a_paraboloid() {
    // y = h·x²/r² from the vertex to (r, h): the tangents meet at (r/2, 0)
    // and the weight is 1 (ρ = ½). The bowl filled to its rim, revolved
    // about its axis: V = πr²h/2; its surface is the paraboloid,
    // S = πr/(6h²)·((r² + 4h²)^{3/2} − r³), and the lid πr².
    let (r, h): (f64, f64) = (0.4, 0.3);
    let arc =
        ProfileCurve::conic(DVec2::ZERO, DVec2::new(r / 2.0, 0.0), DVec2::new(r, h), 1.0).unwrap();
    let body = revolve(vec![
        arc,
        line(DVec2::new(r, h), DVec2::new(0.0, h)),
        line(DVec2::new(0.0, h), DVec2::ZERO),
    ]);
    let bowl = PI * r / (6.0 * h * h) * ((r * r + 4.0 * h * h).powf(1.5) - r.powi(3));
    check_body(&body, PI * r * r * h / 2.0, bowl + PI * r * r);
}

#[test]
fn a_revolved_hyperbolic_arc_is_a_hyperboloid() {
    // x² − y² = 1 for |y| ≤ Y = sinh s: from (cosh s, −Y) to (cosh s, Y),
    // the tangents meeting at (1/cosh s, 0), weight cosh s (the vertex is
    // the shoulder). Revolved about the y axis: the one-sheet hyperboloid
    // x² + z² − y² = 1 between its two disks, V = π(2Y + 2Y³/3) and
    // lateral area 2π∫√(1 + 2y²)dy.
    let s: f64 = 0.9;
    let (c, y) = (s.cosh(), s.sinh());
    let arc = ProfileCurve::conic(
        DVec2::new(c, -y),
        DVec2::new(1.0 / c, 0.0),
        DVec2::new(c, y),
        c,
    )
    .unwrap();
    for q in arc.sample(500).unwrap() {
        assert!((q.x * q.x - q.y * q.y - 1.0).abs() < 1e-14);
    }
    let body = revolve(vec![
        arc,
        line(DVec2::new(c, y), DVec2::new(0.0, y)),
        line(DVec2::new(0.0, y), DVec2::new(0.0, -y)),
        line(DVec2::new(0.0, -y), DVec2::new(c, -y)),
    ]);
    let r2 = 2f64.sqrt();
    let half = 0.5 * y * (1.0 + 2.0 * y * y).sqrt() + (r2 * y).asinh() / (2.0 * r2);
    let lateral = 2.0 * PI * 2.0 * half;
    check_body(
        &body,
        PI * (2.0 * y + 2.0 * y.powi(3) / 3.0),
        lateral + 2.0 * PI * c * c,
    );
}

proptest! {
    // Each case extrudes, integrates and tessellates a solid; 48 keep the
    // suite fast (2,000 were run once in release, design note 0009 §9).
    #![proptest_config(ProptestConfig::with_cases(48))]

    /// Random triangles and weights: the segment's region has the area of
    /// the closed form (affine invariance), and its extrusion is a valid,
    /// watertight solid of that area times the height.
    #[test]
    fn random_segments_have_their_exact_areas(
        x1 in -2.0f64..2.0, y1 in 0.3f64..3.0, x2 in 1.0f64..4.0, y2 in -1.0f64..1.0,
        u in -6.0f64..6.0, h in 0.05f64..2.0,
    ) {
        let p = [DVec2::ZERO, DVec2::new(x1, y1), DVec2::new(x2, y2)];
        let tri = 0.5 * (p[1] - p[0]).perp_dot(p[2] - p[0]).abs();
        prop_assume!(tri > 0.2);
        let w = u.exp();
        let arc = ProfileCurve::conic(p[0], p[1], p[2], w).unwrap();
        let curves = vec![arc, line(p[2], p[0])];
        let areas = region_areas(&curves);
        prop_assert_eq!(areas.len(), 1);
        prop_assert!((areas[0] - ratio(w) * tri).abs() <= 1e-11 * ratio(w) * tri);
        let body = k()
            .extrude(
                &Profile { plane: Frame::WORLD, curves },
                &ExtrudeSpec { distance: h, symmetric: false },
            )
            .unwrap();
        prop_assert!(k().validate(&body).is_empty());
        let m = k().mass_properties(&body, 1e-10).unwrap();
        prop_assert!((m.volume - ratio(w) * tri * h).abs() <= 1e-9 * ratio(w) * tri * h);
        assert_watertight(&k().tessellate(&body, &TessellationOptions::default()).unwrap());
    }
}
