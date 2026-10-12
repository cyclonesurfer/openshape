//! Curve contacts for sketch tools (design note 0007 §6): where two curves
//! meet, in the arrangement's parameter conventions, against closed forms.

mod arrangement_util;

use std::f64::consts::{FRAC_PI_2, PI, TAU};

use arrangement_util::*;
use opendatum_kernel::arrangement::{arrange, curve_contacts};
use opendatum_kernel::{CurveContact, DVec2, Kernel, KernelError, ProfileCurve, Tolerance};
use proptest::prelude::*;

fn contacts(a: &ProfileCurve, b: &ProfileCurve) -> Vec<CurveContact> {
    let c = k().curve_contacts(a, b).expect("contacts");
    for x in &c {
        // Both parameters name the reported point, within the tolerance.
        let pa = k().curve_point(a, x.ta).unwrap();
        let pb = k().curve_point(b, x.tb).unwrap();
        assert!((pa - x.point).length() < 1e-15, "{x:?}");
        assert!((pa - pb).length() <= tol(), "{x:?}: {pa} vs {pb}");
    }
    for w in c.windows(2) {
        assert!(w[0].ta <= w[1].ta, "sorted along a");
        assert!((w[0].point - w[1].point).length() > tol(), "merged");
    }
    c
}

#[test]
fn crossing_lines_meet_once() {
    let c = contacts(
        &line(v(0.0, 0.0), v(4.0, 4.0)),
        &line(v(0.0, 4.0), v(4.0, 0.0)),
    );
    assert_eq!(c.len(), 1);
    assert_rel(c[0].ta, 0.5, 1e-15);
    assert_rel(c[0].tb, 0.5, 1e-15);
    assert!((c[0].point - v(2.0, 2.0)).length() < 1e-15);
}

#[test]
fn nearly_parallel_lines_far_from_the_origin_still_cross() {
    // Within the tolerance of each other over 1 km, crossing at a
    // 1e-10 rad angle near (1952, 1048): rounding moves the crossing's
    // parameters by about a millimetre along the lines, which must not
    // lose it (the arrangement then saw the lines swap sides between two
    // vertices with no crossing between them).
    let t = tol();
    let r = v(2000.0, 1000.0 + 1.1 * t);
    let chord = line(v(1000.0, 2000.0), r);
    let long = line(v(3000.0, 0.0), v(3.0 * t, 3000.0));
    let c = contacts(&chord, &long);
    assert_eq!(c.len(), 2, "{c:?}");
    assert!(c.iter().any(|x| (x.point - r).length() <= t), "{c:?}");
    let f = 0.1 / 2.1;
    let crossing = r + (v(1000.0, 2000.0) - r) * f;
    assert!(
        c.iter().any(|x| (x.point - crossing).length() < 1e-2),
        "{c:?}"
    );
}

#[test]
fn segments_that_miss_do_not_meet() {
    // The infinite lines cross at (5, 5), outside both segments.
    let c = contacts(
        &line(v(0.0, 0.0), v(4.0, 4.0)),
        &line(v(0.0, 10.0), v(4.0, 6.0)),
    );
    assert!(c.is_empty());
}

#[test]
fn line_through_a_circle_meets_it_twice() {
    let c = contacts(&line(v(-2.0, 0.5), v(2.0, 0.5)), &circle(v(0.0, 0.0), 1.0));
    assert_eq!(c.len(), 2);
    let x = 0.75f64.sqrt();
    assert!((c[0].point - v(-x, 0.5)).length() < 1e-15);
    assert!((c[1].point - v(x, 0.5)).length() < 1e-15);
    // Circle parameters are counter-clockwise angles from +x.
    assert_rel(c[0].tb, PI - (0.5f64).asin(), 1e-14);
    assert_rel(c[1].tb, (0.5f64).asin(), 1e-14);
}

#[test]
fn tangent_line_touches_once() {
    let c = contacts(&line(v(-2.0, 1.0), v(2.0, 1.0)), &circle(v(0.0, 0.0), 1.0));
    assert_eq!(c.len(), 1);
    assert!((c[0].point - v(0.0, 1.0)).length() <= tol());
    assert_rel(c[0].tb, FRAC_PI_2, 1e-6);
}

/// A line along the tangent at an end of a gently curved conic touches it
/// once, at that end. The curves stay within the tolerance of each other
/// for about 0.4 mm there (radius of curvature 0.87 m), and refinement used
/// to stop microns along it, reporting two more contacts.
#[test]
fn a_line_along_an_end_tangent_touches_once_at_the_end() {
    let (s, c, e) = (
        v(-0.008734990039930208, -0.03091914159280349),
        v(0.024128675730140697, 0.0),
        v(0.0, 0.0),
    );
    for w in [1.8832086295738748, 1.0, 0.4] {
        let conic = ProfileCurve::conic(s, c, e, w).unwrap();
        for (end, d) in [(s, (c - s).normalize()), (e, (e - c).normalize())] {
            for l in [
                line(end - d * 0.03, end + d * 0.03),
                line(end, end + d * 0.03),
                line(end - d * 0.03, end),
            ] {
                let x = contacts(&l, &conic);
                assert_eq!(x.len(), 1, "w {w}: {x:?}");
                assert!((x[0].point - end).length() < 1e-15, "w {w}: {x:?}");
            }
        }
    }
}

#[test]
fn arc_parameters_run_along_the_arc() {
    // A clockwise arc from angle π/2 sweeping −π (through +x): its
    // parameter is the angle travelled from its start.
    let a = arc(v(0.0, 0.0), 2.0, FRAC_PI_2, -PI);
    let c = contacts(&a, &line(v(0.0, 0.0), v(3.0, 0.0)));
    assert_eq!(c.len(), 1);
    assert_rel(c[0].ta, FRAC_PI_2, 1e-14);
    assert!((c[0].point - v(2.0, 0.0)).length() < 1e-15);
    // The line's own start is at the arc's centre, not on the arc.
    assert_rel(c[0].tb, 2.0 / 3.0, 1e-14);
}

#[test]
fn end_points_on_the_other_curve_count() {
    // A T junction: the stem ends on the bar.
    let c = contacts(
        &line(v(0.0, 0.0), v(4.0, 0.0)),
        &line(v(1.0, 0.0), v(1.0, 3.0)),
    );
    assert_eq!(c.len(), 1);
    assert_rel(c[0].ta, 0.25, 1e-15);
    assert_eq!(c[0].tb, 0.0);
}

#[test]
fn overlaps_are_reported_by_their_ends() {
    let c = contacts(
        &line(v(0.0, 0.0), v(4.0, 0.0)),
        &line(v(1.0, 0.0), v(6.0, 0.0)),
    );
    assert_eq!(c.len(), 2);
    assert_rel(c[0].ta, 0.25, 1e-15);
    assert_eq!(c[1].ta, 1.0);
    // Two arcs on one circle.
    let c = contacts(
        &arc(v(0.0, 0.0), 1.0, 0.0, PI),
        &arc(v(0.0, 0.0), 1.0, FRAC_PI_2, PI),
    );
    assert_eq!(c.len(), 2);
    assert_rel(c[0].ta, FRAC_PI_2, 1e-14);
    assert_rel(c[1].ta, PI, 1e-14);
}

#[test]
fn circles_crossing_and_touching() {
    let c = contacts(&circle(v(0.0, 0.0), 1.0), &circle(v(1.0, 0.0), 1.0));
    assert_eq!(c.len(), 2);
    assert_rel(c[0].ta, PI / 3.0, 1e-14);
    assert_rel(c[1].ta, 5.0 * PI / 3.0, 1e-14);
    // Externally tangent: one contact.
    let c = contacts(&circle(v(0.0, 0.0), 1.0), &circle(v(3.0, 0.0), 2.0));
    assert_eq!(c.len(), 1);
    assert!((c[0].point - v(1.0, 0.0)).length() <= tol());
}

#[test]
fn contacts_near_a_circle_seam_merge() {
    // A line through the circle's parameter seam (angle 0 = 2π).
    let c = contacts(&circle(v(0.0, 0.0), 1.0), &line(v(1.0, -1.0), v(1.0, 1.0)));
    assert_eq!(c.len(), 1, "{c:?}");
    assert!((c[0].point - v(1.0, 0.0)).length() <= tol());
}

#[test]
fn spline_contacts_match_the_polyline_closed_form() {
    // A degree-1 spline is its control polygon.
    let s = ProfileCurve::Spline {
        degree: 1,
        knots: vec![0.0, 0.0, 0.5, 1.0, 1.0],
        points: vec![v(0.0, 0.0), v(2.0, 2.0), v(4.0, 0.0)],
        weights: None,
    };
    let c = contacts(&s, &line(v(0.0, 1.0), v(4.0, 1.0)));
    assert_eq!(c.len(), 2);
    assert!((c[0].point - v(1.0, 1.0)).length() <= tol());
    assert!((c[1].point - v(3.0, 1.0)).length() <= tol());
    assert_rel(c[0].ta, 0.25, 1e-6);
    assert_rel(c[1].ta, 0.75, 1e-6);
}

#[test]
fn contacts_agree_with_the_arrangement() {
    // The same pair through arrange: every arrangement vertex on both
    // curves is a contact, and vice versa.
    let a =
        ProfileCurve::interpolate(&[v(0.0, 0.0), v(1.0, 2.0), v(3.0, -1.0), v(4.0, 1.0)]).unwrap();
    let b = circle(v(2.0, 0.5), 1.4);
    let c = contacts(&a, &b);
    let arr = arrange(&[a.clone(), b.clone()], &Tolerance::default()).unwrap();
    let both: Vec<DVec2> = arr
        .vertices
        .iter()
        .filter(|x| x.curves == vec![0, 1])
        .map(|x| x.point)
        .collect();
    assert_eq!(c.len(), both.len(), "{c:?} vs {both:?}");
    for x in &c {
        assert!(both.iter().any(|p| (*p - x.point).length() <= 2.0 * tol()));
    }
}

#[test]
fn degenerate_input_names_the_curve() {
    let e = curve_contacts(
        &line(v(0.0, 0.0), v(1.0, 0.0)),
        &circle(v(0.0, 0.0), 0.0),
        &Tolerance::default(),
    )
    .unwrap_err();
    match e {
        KernelError::InvalidInput(m) => assert!(m.contains("curve 1"), "{m}"),
        e => panic!("{e:?}"),
    }
}

proptest! {
    /// Random line–circle pairs against the quadratic's closed form.
    #[test]
    fn line_circle_contacts_match_the_closed_form(
        cx in -2.0..2.0f64, cy in -2.0..2.0f64, r in 0.1..3.0f64,
        ax in -4.0..4.0f64, ay in -4.0..4.0f64, angle in 0.0..TAU, len in 0.5..8.0f64,
    ) {
        let a = v(ax, ay);
        let d = DVec2::from_angle(angle);
        let b = a + d * len;
        let c = contacts(&line(a, b), &circle(v(cx, cy), r));
        // |a + d s − c|² = r², s in [0, len].
        let w = a - v(cx, cy);
        let (bq, cq) = (w.dot(d), w.length_squared() - r * r);
        let disc = bq * bq - cq;
        let mut expected: Vec<f64> = Vec::new();
        if disc > 1e-9 {
            for s in [-bq - disc.sqrt(), -bq + disc.sqrt()] {
                if (1e-6..=len - 1e-6).contains(&s) {
                    expected.push(s / len);
                }
            }
        }
        // Skip near-tangent and near-end cases, where tolerance decides.
        let near_edge = disc.abs() <= 1e-9
            || [-bq - disc.max(0.0).sqrt(), -bq + disc.max(0.0).sqrt()]
                .iter()
                .any(|s| s.abs() < 1e-6 || (s - len).abs() < 1e-6);
        prop_assume!(!near_edge);
        prop_assert_eq!(c.len(), expected.len());
        for (x, s) in c.iter().zip(&expected) {
            prop_assert!((x.ta - s).abs() < 1e-12, "{} vs {}", x.ta, s);
        }
    }
}
