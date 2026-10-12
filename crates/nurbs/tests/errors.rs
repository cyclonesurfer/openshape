//! Invalid input returns errors and never panics (acceptance test 8).

use glam::{DVec3, DVec4};
use opendatum_nurbs::{KnotVector, NurbsCurve, NurbsError, NurbsSurface};
use proptest::prelude::*;

fn pts(n: usize) -> Vec<DVec3> {
    (0..n)
        .map(|i| DVec3::new(i as f64, (i * i) as f64, 1.0))
        .collect()
}

#[test]
fn knot_vector_errors() {
    assert_eq!(
        KnotVector::new(vec![]),
        Err(NurbsError::TooFewKnots { found: 0 })
    );
    assert_eq!(
        KnotVector::new(vec![0.0, f64::INFINITY]),
        Err(NurbsError::NonFiniteKnot { index: 1 })
    );
    assert_eq!(
        KnotVector::new(vec![0.0, 1.0, 0.9, 2.0]),
        Err(NurbsError::DecreasingKnots { index: 2 })
    );
    assert_eq!(KnotVector::new(vec![1.0; 6]), Err(NurbsError::EmptyDomain));
    assert_eq!(
        KnotVector::clamped_uniform(0, 4),
        Err(NurbsError::ZeroDegree)
    );
}

#[test]
fn curve_constructor_errors() {
    let k3 = vec![0., 0., 0., 0., 1., 1., 1., 1.];
    assert_eq!(
        NurbsCurve::new(0, vec![0., 1.], pts(1), None),
        Err(NurbsError::ZeroDegree)
    );
    assert!(matches!(
        NurbsCurve::new(3, vec![0., 0., 0., 1., 1., 1.], pts(2), None),
        Err(NurbsError::TooFewControlPoints {
            degree: 3,
            required: 4,
            found: 2
        })
    ));
    assert!(matches!(
        NurbsCurve::new(3, k3.clone(), pts(5), None),
        Err(NurbsError::KnotCountMismatch {
            expected: 9,
            found: 8,
            ..
        })
    ));
    assert!(matches!(
        NurbsCurve::new(3, vec![0., 0., 0., 0., 1., 0.5, 1., 1., 1.], pts(5), None),
        Err(NurbsError::DecreasingKnots { .. })
    ));
    assert!(matches!(
        NurbsCurve::new(
            3,
            vec![0., 0., 0., 0., f64::NAN, 1., 1., 1., 1.],
            pts(5),
            None
        ),
        Err(NurbsError::NonFiniteKnot { index: 4 })
    ));
    // Interior knot of multiplicity > degree.
    assert!(matches!(
        NurbsCurve::new(2, vec![0., 0., 0., 0.5, 0.5, 0.5, 1., 1., 1.], pts(6), None),
        Err(NurbsError::KnotMultiplicityTooHigh {
            multiplicity: 3,
            max: 2,
            ..
        })
    ));
    // End knot of multiplicity > degree + 1.
    assert!(matches!(
        NurbsCurve::new(1, vec![0., 0., 0., 1., 1.], pts(3), None),
        Err(NurbsError::KnotMultiplicityTooHigh { .. })
    ));
    // Empty domain [u_p, u_{n+1}].
    assert!(matches!(
        NurbsCurve::new(2, vec![0., 1., 1., 1., 1., 2.], pts(3), None),
        Err(NurbsError::EmptyDomain) | Err(NurbsError::KnotMultiplicityTooHigh { .. })
    ));
    assert_eq!(
        NurbsCurve::new(3, k3.clone(), pts(4), Some(vec![1.0; 3])),
        Err(NurbsError::WeightCountMismatch {
            expected: 4,
            found: 3
        })
    );
    for bad in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(matches!(
            NurbsCurve::new(3, k3.clone(), pts(4), Some(vec![1.0, 1.0, bad, 1.0])),
            Err(NurbsError::InvalidWeight { index: 2, .. })
        ));
    }
    let mut bad_points = pts(4);
    bad_points[1].y = f64::NAN;
    assert_eq!(
        NurbsCurve::new(3, k3.clone(), bad_points, None),
        Err(NurbsError::NonFiniteControlPoint { index: 1 })
    );
    assert!(matches!(
        NurbsCurve::from_homogeneous(1, vec![0., 0., 1., 1.], vec![DVec4::W, DVec4::ZERO]),
        Err(NurbsError::InvalidWeight { index: 1, .. })
    ));
    assert!(matches!(
        NurbsCurve::from_homogeneous(
            1,
            vec![0., 0., 1., 1.],
            vec![DVec4::W, DVec4::new(1.0, 0.0, 0.0, 1e-320)]
        ),
        Err(NurbsError::NonFiniteControlPoint { index: 1 })
    ));
    assert!(NurbsCurve::line(DVec3::ZERO, DVec3::splat(f64::NAN)).is_err());
}

#[test]
fn surface_constructor_errors() {
    let k = vec![0., 0., 1., 1.];
    let net = vec![pts(2), pts(2)];
    assert!(NurbsSurface::new(1, 1, k.clone(), k.clone(), net.clone(), None).is_ok());
    assert_eq!(
        NurbsSurface::new(0, 1, k.clone(), k.clone(), net.clone(), None),
        Err(NurbsError::ZeroDegree)
    );
    assert!(matches!(
        NurbsSurface::new(1, 1, k.clone(), k.clone(), vec![pts(2), pts(3)], None),
        Err(NurbsError::RaggedControlNet {
            row: 1,
            expected: 2,
            found: 3
        })
    ));
    assert!(NurbsSurface::new(1, 1, k.clone(), k.clone(), vec![], None).is_err());
    assert!(NurbsSurface::new(1, 1, k.clone(), k.clone(), vec![vec![], vec![]], None).is_err());
    assert!(matches!(
        NurbsSurface::new(
            1,
            1,
            vec![0., 0., 0.5, 1., 1.],
            k.clone(),
            net.clone(),
            None
        ),
        Err(NurbsError::KnotCountMismatch { .. })
    ));
    assert!(matches!(
        NurbsSurface::new(1, 2, k.clone(), k.clone(), net.clone(), None),
        Err(NurbsError::TooFewControlPoints { .. })
    ));
    assert!(matches!(
        NurbsSurface::new(
            1,
            1,
            k.clone(),
            k.clone(),
            net.clone(),
            Some(vec![vec![1.0, 1.0]])
        ),
        Err(NurbsError::WeightCountMismatch { .. })
    ));
    assert!(matches!(
        NurbsSurface::new(
            1,
            1,
            k.clone(),
            k.clone(),
            net.clone(),
            Some(vec![vec![1.0, 1.0], vec![1.0]])
        ),
        Err(NurbsError::RaggedControlNet { row: 1, .. })
    ));
    assert!(matches!(
        NurbsSurface::new(
            1,
            1,
            k.clone(),
            k.clone(),
            net.clone(),
            Some(vec![vec![1.0, 1.0], vec![1.0, -2.0]])
        ),
        Err(NurbsError::InvalidWeight { index: 3, .. })
    ));
    let mut bad = net.clone();
    bad[1][0].z = f64::INFINITY;
    assert_eq!(
        NurbsSurface::new(1, 1, k.clone(), k.clone(), bad, None),
        Err(NurbsError::NonFiniteControlPoint { index: 2 })
    );
    assert!(matches!(
        NurbsSurface::from_homogeneous(
            1,
            1,
            k.clone(),
            k.clone(),
            vec![vec![DVec4::W, DVec4::W], vec![DVec4::W, -DVec4::W]]
        ),
        Err(NurbsError::InvalidWeight { index: 3, .. })
    ));
    assert!(NurbsSurface::bilinear(DVec3::ZERO, DVec3::X, DVec3::Y, DVec3::NAN).is_err());
}

#[test]
fn constructor_argument_errors() {
    let line = NurbsCurve::line(DVec3::X, DVec3::new(1.0, 0.0, 1.0)).unwrap();
    for (axis, angle) in [
        (DVec3::ZERO, 1.0),
        (DVec3::NAN, 1.0),
        (DVec3::Z, 0.0),
        (DVec3::Z, 7.0),
        (DVec3::Z, -7.0),
        (DVec3::Z, f64::NAN),
    ] {
        assert!(matches!(
            NurbsSurface::revolution(&line, DVec3::ZERO, axis, angle),
            Err(NurbsError::InvalidArgument(_))
        ));
    }
    assert!(NurbsSurface::revolution(&line, DVec3::NAN, DVec3::Z, 1.0).is_err());
    for dir in [DVec3::ZERO, DVec3::NAN, DVec3::INFINITY] {
        assert!(matches!(
            NurbsSurface::extrusion(&line, dir),
            Err(NurbsError::InvalidArgument(_))
        ));
    }
}

#[test]
fn errors_display() {
    let e = NurbsError::KnotCountMismatch {
        degree: 3,
        control_points: 5,
        expected: 9,
        found: 8,
    };
    assert_eq!(
        e.to_string(),
        "expected 9 knots for 5 control points of degree 3, found 8"
    );
}

fn arb_f64() -> impl Strategy<Value = f64> {
    prop_oneof![
        8 => -10.0..10.0f64,
        1 => Just(0.0),
        1 => Just(f64::NAN),
        1 => Just(f64::INFINITY),
        1 => Just(-1.0),
    ]
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(2048))]

    /// Arbitrary (mostly invalid) input never panics; anything accepted
    /// evaluates to finite values.
    #[test]
    fn arbitrary_curve_input_never_panics(
        degree in 0usize..5,
        mut knots in prop::collection::vec(arb_f64(), 0..14),
        sort in any::<bool>(),
        points in prop::collection::vec((arb_f64(), arb_f64(), arb_f64()), 0..9),
        weights in prop::option::of(prop::collection::vec(arb_f64(), 0..9)),
        t in arb_f64(),
    ) {
        if sort {
            knots.sort_by(f64::total_cmp);
        }
        let points: Vec<DVec3> = points.into_iter().map(|(x, y, z)| DVec3::new(x, y, z)).collect();
        if let Ok(c) = NurbsCurve::new(degree, knots, points, weights) {
            let (a, b) = c.domain();
            let tt = if t.is_finite() { t } else { 0.5 * (a + b) };
            prop_assert!(c.point(tt).is_finite());
            prop_assert!(c.derivatives(tt, 3).iter().all(|d| d.is_finite()));
            let (s, q) = c.closest_point(DVec3::new(1.0, 2.0, 3.0));
            prop_assert!(s >= a && s <= b && q.is_finite());
            prop_assert!(c.length().is_finite());
            let mid = 0.5 * (a + b);
            let (l, r) = c.split(mid);
            prop_assert!((l.point(mid) - r.point(mid)).length() <= 1e-9 * (1.0 + c.point(mid).length()));
            // Evaluation outside the domain or at NaN must not panic.
            let _ = c.point(f64::NAN);
            let _ = c.derivatives(f64::INFINITY, 2);
            let _ = c.closest_point(DVec3::NAN);
            let _ = c.insert_knot(f64::NAN, 2);
        }
    }

    #[test]
    fn arbitrary_surface_input_never_panics(
        du in 0usize..4,
        dv in 0usize..4,
        mut ku in prop::collection::vec(arb_f64(), 0..9),
        mut kv in prop::collection::vec(arb_f64(), 0..9),
        rows in prop::collection::vec(prop::collection::vec((arb_f64(), arb_f64(), arb_f64()), 0..5), 0..5),
        weights in prop::option::of(prop::collection::vec(prop::collection::vec(arb_f64(), 0..5), 0..5)),
    ) {
        ku.sort_by(f64::total_cmp);
        kv.sort_by(f64::total_cmp);
        let rows: Vec<Vec<DVec3>> = rows
            .into_iter()
            .map(|r| r.into_iter().map(|(x, y, z)| DVec3::new(x, y, z)).collect())
            .collect();
        if let Ok(s) = NurbsSurface::new(du, dv, ku, kv, rows, weights) {
            let ((ua, ub), (va, vb)) = s.domain();
            let (u, v) = (0.5 * (ua + ub), 0.5 * (va + vb));
            prop_assert!(s.point(u, v).is_finite());
            prop_assert!(s.derivatives(u, v, 2).iter().flatten().all(|d| d.is_finite()));
            let ((cu, cv), q) = s.closest_point(DVec3::new(1.0, 2.0, 3.0));
            prop_assert!(cu >= ua && cu <= ub && cv >= va && cv <= vb && q.is_finite());
            let _ = s.normal(u, v);
            let _ = s.point(f64::NAN, f64::INFINITY);
            let _ = s.closest_point(DVec3::NAN);
            let _ = s.iso_curve_u(f64::NAN);
        }
    }
}
