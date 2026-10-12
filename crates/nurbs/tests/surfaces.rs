//! Surface evaluation, constructors and point inversion (acceptance test 7).

mod common;

use std::f64::consts::{FRAC_PI_2, PI, TAU};

use common::{arb_curve, arb_surface, close};
use glam::{DMat4, DQuat, DVec3};
use opendatum_nurbs::{NurbsCurve, NurbsSurface};
use proptest::prelude::*;

/// A parameter strictly inside the knot span selected by `which`, at
/// relative position `0.1 + 0.8 * rel`, plus the span width.
fn inside_span(knots: &[f64], degree: usize, which: f64, rel: f64) -> (f64, f64) {
    let mut breaks: Vec<f64> = knots[degree..knots.len() - degree].to_vec();
    breaks.dedup();
    let i = ((which * (breaks.len() - 1) as f64) as usize).min(breaks.len() - 2);
    let w = breaks[i + 1] - breaks[i];
    (breaks[i] + w * (0.1 + 0.8 * rel), w)
}

fn sphere(center: DVec3, radius: f64) -> NurbsSurface {
    let profile = NurbsCurve::circular_arc(center, -DVec3::Z, DVec3::X, radius, 0.0, PI).unwrap();
    NurbsSurface::revolution(&profile, center, DVec3::Z, TAU).unwrap()
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(256))]

    #[test]
    fn surface_derivatives_match_finite_differences(
        s in arb_surface(),
        wu in 0.0..1.0f64, ru in 0.0..1.0f64,
        wv in 0.0..1.0f64, rv in 0.0..1.0f64,
    ) {
        let (u, hu) = inside_span(s.knots_u().as_slice(), s.degree_u(), wu, ru);
        let (v, hv) = inside_span(s.knots_v().as_slice(), s.degree_v(), wv, rv);
        let (hu, hv) = (1e-6 * hu, 1e-6 * hv);
        let d = s.derivatives(u, v, 2);
        prop_assert_eq!(d.len(), 3);
        prop_assert_eq!(d[0].len(), 3);
        prop_assert_eq!(d[1].len(), 2);
        prop_assert_eq!(d[2].len(), 1);
        prop_assert!((d[0][0] - s.point(u, v)).length() <= 1e-12 * d[0][0].length().max(1.0));

        let fd_u = (s.point(u + hu, v) - s.point(u - hu, v)) / (2.0 * hu);
        let fd_v = (s.point(u, v + hv) - s.point(u, v - hv)) / (2.0 * hv);
        prop_assert!(close(d[1][0], fd_u, 1e-5), "S_u {:?} vs {:?}", d[1][0], fd_u);
        prop_assert!(close(d[0][1], fd_v, 1e-5), "S_v {:?} vs {:?}", d[0][1], fd_v);

        let du_p = s.derivatives(u + hu, v, 1);
        let du_m = s.derivatives(u - hu, v, 1);
        let dv_p = s.derivatives(u, v + hv, 1);
        let dv_m = s.derivatives(u, v - hv, 1);
        let fd_uu = (du_p[1][0] - du_m[1][0]) / (2.0 * hu);
        let fd_uv = (dv_p[1][0] - dv_m[1][0]) / (2.0 * hv);
        let fd_vu = (du_p[0][1] - du_m[0][1]) / (2.0 * hu);
        let fd_vv = (dv_p[0][1] - dv_m[0][1]) / (2.0 * hv);
        prop_assert!(close(d[2][0], fd_uu, 1e-5), "S_uu {:?} vs {:?}", d[2][0], fd_uu);
        prop_assert!(close(d[1][1], fd_uv, 1e-5), "S_uv {:?} vs {:?}", d[1][1], fd_uv);
        prop_assert!(close(d[1][1], fd_vu, 1e-5), "S_vu {:?} vs {:?}", d[1][1], fd_vu);
        prop_assert!(close(d[0][2], fd_vv, 1e-5), "S_vv {:?} vs {:?}", d[0][2], fd_vv);

        if let Some(n) = s.normal(u, v) {
            prop_assert!((n.length() - 1.0).abs() < 1e-12);
            prop_assert!(n.dot(d[1][0]).abs() <= 1e-9 * d[1][0].length());
            prop_assert!(n.dot(d[0][1]).abs() <= 1e-9 * d[0][1].length());
        }
    }

    #[test]
    fn surface_closest_point_round_trip(
        s in arb_surface(),
        ru in 0.0..=1.0f64,
        rv in 0.0..=1.0f64,
    ) {
        let ((ua, ub), (va, vb)) = s.domain();
        let (u0, v0) = (ua + (ub - ua) * ru, va + (vb - va) * rv);
        let p = s.point(u0, v0);
        let ((u, v), q) = s.closest_point(p);
        prop_assert!(u >= ua && u <= ub && v >= va && v <= vb);
        prop_assert!(q.distance(p) <= 1e-9 * p.length().max(1.0), "distance {}", q.distance(p));
        prop_assert!((s.point(u, v) - q).length() <= 1e-12 * q.length().max(1.0));
    }

    #[test]
    fn surface_closest_point_beats_grid_sampling(
        s in arb_surface(),
        x in -2.0..6.0f64,
        y in -2.0..6.0f64,
        z in -4.0..4.0f64,
    ) {
        let p = DVec3::new(x, y, z);
        let (_, q) = s.closest_point(p);
        let ((ua, ub), (va, vb)) = s.domain();
        let n = 120;
        let mut dense = f64::INFINITY;
        for i in 0..=n {
            for j in 0..=n {
                let u = ua + (ub - ua) * i as f64 / n as f64;
                let v = va + (vb - va) * j as f64 / n as f64;
                dense = dense.min(s.point(u, v).distance(p));
            }
        }
        prop_assert!(q.distance(p) <= dense + 1e-9, "closest {} > grid {}", q.distance(p), dense);
    }

    #[test]
    fn extrusion_is_curve_plus_offset(
        c in arb_curve(),
        dx in -5.0..5.0f64,
        dy in -5.0..5.0f64,
        dz in 0.1..5.0f64,
        rel in 0.0..=1.0f64,
        v in 0.0..=1.0f64,
    ) {
        let dir = DVec3::new(dx, dy, dz);
        let s = NurbsSurface::extrusion(&c, dir).unwrap();
        let (a, b) = c.domain();
        let t = a + (b - a) * rel;
        prop_assert_eq!(s.domain(), ((a, b), (0.0, 1.0)));
        let expected = c.point(t) + dir * v;
        prop_assert!((s.point(t, v) - expected).length() <= 1e-12 * expected.length().max(1.0));
        let d = s.derivatives(t, v, 1);
        prop_assert!((d[0][1] - dir).length() <= 1e-12 * dir.length().max(1.0));
        let ct = c.derivatives(t, 1)[1];
        prop_assert!((d[1][0] - ct).length() <= 1e-10 * ct.length().max(1.0));
    }

    #[test]
    fn iso_curves_and_transforms_match(
        s in arb_surface(),
        ru in 0.0..=1.0f64,
        rv in 0.0..=1.0f64,
        angle in -3.0..3.0f64,
    ) {
        let ((ua, ub), (va, vb)) = s.domain();
        let (u, v) = (ua + (ub - ua) * ru, va + (vb - va) * rv);
        let cu = s.iso_curve_u(u);
        let cv = s.iso_curve_v(v);
        prop_assert_eq!(cu.domain(), (va, vb));
        prop_assert_eq!(cv.domain(), (ua, ub));
        prop_assert!((cu.point(v) - s.point(u, v)).length() <= 1e-12 * 10.0);
        prop_assert!((cv.point(u) - s.point(u, v)).length() <= 1e-12 * 10.0);
        let m = DMat4::from_rotation_translation(
            DQuat::from_axis_angle(DVec3::new(0.3, 1.0, -0.2).normalize(), angle),
            DVec3::new(1.0, -2.0, 3.0),
        ) * DMat4::from_scale(DVec3::new(2.0, 0.5, 1.5));
        let moved = s.transformed(&m);
        let expected = m.transform_point3(s.point(u, v));
        prop_assert!((moved.point(u, v) - expected).length() <= 1e-11 * expected.length().max(1.0));
    }
}

#[test]
fn revolving_a_line_parallel_to_the_axis_gives_a_cylinder() {
    let origin = DVec3::new(1.0, 2.0, 3.0);
    let axis = DVec3::new(1.0, 1.0, 1.0);
    let a = axis.normalize();
    let radial = a.any_orthonormal_vector();
    let r = 2.5;
    let line =
        NurbsCurve::line(origin + radial * r - a * 1.0, origin + radial * r + a * 4.0).unwrap();
    for &angle in &[0.3, FRAC_PI_2, 2.0, PI, 4.0, TAU, -1.0, -TAU] {
        let s = NurbsSurface::revolution(&line, origin, axis, angle).unwrap();
        assert_eq!(s.degree_u(), 2);
        assert_eq!(s.degree_v(), 1);
        for i in 0..=30 {
            for j in 0..=10 {
                let (u, v) = (i as f64 / 30.0, j as f64 / 10.0);
                let p = s.point(u, v);
                let rel = p - origin;
                let dist_axis = (rel - a * rel.dot(a)).length();
                assert!(
                    (dist_axis - r).abs() <= 1e-12 * r,
                    "angle {angle}: {dist_axis}"
                );
                // The height along the axis is that of the profile point.
                let h = rel.dot(a);
                assert!((h - (-1.0 + 5.0 * v)).abs() < 1e-12);
            }
        }
        // S(0, v) is the profile; the swept angle is |angle|.
        assert!((s.point(0.0, 0.5) - line.point(0.5)).length() < 1e-14);
        let start = s.point(0.0, 0.0) - origin;
        let end = s.point(1.0, 0.0) - origin;
        let cos = (start - a * start.dot(a))
            .normalize()
            .dot((end - a * end.dot(a)).normalize());
        let expected = angle.abs().cos();
        assert!((cos - expected).abs() < 1e-12, "angle {angle}");
        // Rotation sense: positive angles turn counter-clockwise about the axis.
        let mid = s.point(0.01, 0.0) - origin;
        let turn = start.cross(mid).dot(a);
        assert!(turn * angle.signum() > 0.0);
    }
}

#[test]
fn full_revolution_closes_exactly() {
    let profile = NurbsCurve::interpolate(
        &[
            DVec3::new(1.0, 0.0, 0.0),
            DVec3::new(2.0, 0.0, 1.0),
            DVec3::new(1.5, 0.0, 2.0),
            DVec3::new(0.5, 0.0, 3.0),
        ],
        3,
    )
    .unwrap();
    let s = NurbsSurface::revolution(&profile, DVec3::ZERO, DVec3::Z, TAU).unwrap();
    for j in 0..=20 {
        let v = j as f64 / 20.0;
        assert_eq!(s.point(0.0, v), s.point(1.0, v));
        let p = profile.point(v);
        let r = (p.x * p.x + p.y * p.y).sqrt();
        for i in 0..=40 {
            let q = s.point(i as f64 / 40.0, v);
            assert!(((q.x * q.x + q.y * q.y).sqrt() - r).abs() < 1e-12);
            assert!((q.z - p.z).abs() < 1e-12);
        }
    }
}

#[test]
fn sphere_points_normals_and_closest_points() {
    let c = DVec3::new(0.5, -1.0, 2.0);
    let r = 3.0;
    let s = sphere(c, r);
    for i in 0..=24 {
        for j in 0..=24 {
            let (u, v) = (i as f64 / 24.0, j as f64 / 24.0);
            let p = s.point(u, v);
            assert!(((p - c).length() - r).abs() < 1e-12 * r);
            if let Some(n) = s.normal(u, v) {
                let radial = (p - c) / r;
                assert!((n.dot(radial).abs() - 1.0).abs() < 1e-10, "u {u} v {v}");
            } else {
                // Only the poles are degenerate.
                assert!(
                    j == 0 || j == 24,
                    "unexpected degenerate normal at {u}, {v}"
                );
            }
        }
    }
    // Closest points from outside and inside the sphere are radial.
    for k in 0..50 {
        let t = k as f64;
        let dir = DVec3::new((t * 0.7).cos(), (t * 1.3).sin(), (t * 0.4).cos() * 0.8).normalize();
        for &dist in &[0.2, 2.0, 2.999, 3.5, 40.0] {
            let p = c + dir * dist;
            let (_, q) = s.closest_point(p);
            let expected = c + dir * r;
            assert!(
                (q - expected).length() < 1e-9,
                "dir {dir:?} dist {dist}: {q:?}"
            );
        }
    }
}

#[test]
fn torus_closest_point() {
    let (major, minor) = (3.0, 1.0);
    let profile = NurbsCurve::circular_arc(
        DVec3::new(major, 0.0, 0.0),
        DVec3::X,
        DVec3::Z,
        minor,
        0.0,
        TAU,
    )
    .unwrap();
    let torus = NurbsSurface::revolution(&profile, DVec3::ZERO, DVec3::Z, TAU).unwrap();
    for k in 0..40 {
        let a = k as f64 * 0.37;
        let b = k as f64 * 0.91;
        let ring = DVec3::new(a.cos(), a.sin(), 0.0) * major;
        let off =
            (DVec3::new(a.cos(), a.sin(), 0.0) * b.cos() + DVec3::Z * b.sin()) * (minor * 1.7);
        let p = ring + off;
        let (_, q) = torus.closest_point(p);
        let expected = ring + off.normalize() * minor;
        assert!(
            (q - expected).length() < 1e-9,
            "k {k}: {q:?} vs {expected:?}"
        );
    }
}

#[test]
fn revolution_of_rational_profile_keeps_weights() {
    let profile =
        NurbsCurve::circular_arc(DVec3::new(2.0, 0.0, 0.0), DVec3::X, DVec3::Z, 0.5, 0.0, PI)
            .unwrap();
    let s = NurbsSurface::revolution(&profile, DVec3::ZERO, DVec3::Z, FRAC_PI_2).unwrap();
    assert!(s.is_rational());
    let (nu, nv) = s.control_net_size();
    assert_eq!((nu, nv), (3, profile.control_points().len()));
    let w = s.weights();
    let pw = profile.weights();
    let half = (FRAC_PI_2 / 2.0).cos();
    for j in 0..nv {
        assert!((w[0][j] - pw[j]).abs() < 1e-15);
        assert!((w[1][j] - pw[j] * half).abs() < 1e-15);
        assert!((w[2][j] - pw[j]).abs() < 1e-15);
    }
    // Every point lies on the torus with major radius 2 and minor radius 0.5.
    for i in 0..=10 {
        for j in 0..=10 {
            let p = s.point(i as f64 / 10.0, j as f64 / 10.0);
            let ring = DVec3::new(p.x, p.y, 0.0).normalize() * 2.0;
            assert!(((p - ring).length() - 0.5).abs() < 1e-12);
        }
    }
}

#[test]
fn bilinear_closest_point_and_boundaries() {
    let s = NurbsSurface::bilinear(
        DVec3::ZERO,
        DVec3::new(2.0, 0.0, 0.0),
        DVec3::new(0.0, 1.0, 0.0),
        DVec3::new(2.0, 1.0, 0.0),
    )
    .unwrap();
    let ((u, v), q) = s.closest_point(DVec3::new(0.5, 0.25, 3.0));
    assert!((u - 0.25).abs() < 1e-14 && (v - 0.25).abs() < 1e-14);
    assert!((q - DVec3::new(0.5, 0.25, 0.0)).length() < 1e-14);
    // Outside the patch: clamped to the boundary / corner.
    let ((u, v), q) = s.closest_point(DVec3::new(-1.0, 0.5, 1.0));
    assert_eq!(u, 0.0);
    assert!((v - 0.5).abs() < 1e-14);
    assert!((q - DVec3::new(0.0, 0.5, 0.0)).length() < 1e-14);
    let ((u, v), _) = s.closest_point(DVec3::new(5.0, 5.0, 0.0));
    assert_eq!((u, v), (1.0, 1.0));
}
