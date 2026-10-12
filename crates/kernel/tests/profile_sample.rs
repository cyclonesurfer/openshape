//! `ProfileCurve::sample` and the `Kernel::region_containing` delegation.

use std::f64::consts::{FRAC_PI_2, PI};

use opendatum_kernel::{BrepKernel, DVec2, Kernel, ProfileCurve};

fn close(a: DVec2, b: DVec2) -> bool {
    (a - b).length() < 1e-12
}

#[test]
fn samples_start_and_end_at_the_curve_ends() {
    let line = ProfileCurve::Line {
        start: DVec2::new(1.0, 2.0),
        end: DVec2::new(3.0, 6.0),
    };
    let p = line.sample(4).unwrap();
    assert_eq!(p.len(), 5);
    assert!(close(p[0], DVec2::new(1.0, 2.0)) && close(p[4], DVec2::new(3.0, 6.0)));
    assert!(close(p[2], DVec2::new(2.0, 4.0)));

    let arc = ProfileCurve::Arc {
        center: DVec2::new(1.0, 1.0),
        radius: 2.0,
        start_angle: 0.0,
        sweep_angle: FRAC_PI_2,
    };
    let p = arc.sample(8).unwrap();
    assert_eq!(p.len(), 9);
    assert!(close(p[0], DVec2::new(3.0, 1.0)) && close(p[8], DVec2::new(1.0, 3.0)));
    assert!(p
        .iter()
        .all(|q| ((*q - DVec2::ONE).length() - 2.0).abs() < 1e-12));

    let circle = ProfileCurve::Circle {
        center: DVec2::ZERO,
        radius: 0.5,
    };
    let p = circle.sample(16).unwrap();
    assert_eq!(p.len(), 17);
    assert!(close(p[0], p[16]), "a sampled circle is closed");
    assert!(close(p[8], DVec2::new(-0.5, 0.0)));

    // A quadratic Bézier as a spline: ends interpolated, middle at t = 1/2.
    let spline = ProfileCurve::Spline {
        degree: 2,
        knots: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
        points: vec![DVec2::ZERO, DVec2::new(1.0, 2.0), DVec2::new(2.0, 0.0)],
        weights: None,
    };
    let p = spline.sample(2).unwrap();
    assert!(close(p[0], DVec2::ZERO) && close(p[2], DVec2::new(2.0, 0.0)));
    assert!(close(p[1], DVec2::new(1.0, 1.0)));
    // Zero segments still gives both ends.
    assert_eq!(line.sample(0).unwrap().len(), 2);
}

#[test]
fn an_invalid_spline_is_an_error() {
    let bad = ProfileCurve::Spline {
        degree: 3,
        knots: vec![0.0, 1.0],
        points: vec![DVec2::ZERO, DVec2::X],
        weights: None,
    };
    assert!(bad.sample(4).is_err());
}

#[test]
fn region_containing_through_the_kernel_trait() {
    let k = BrepKernel::default();
    let curves = vec![
        ProfileCurve::Circle {
            center: DVec2::ZERO,
            radius: 1.0,
        },
        ProfileCurve::Circle {
            center: DVec2::ZERO,
            radius: 0.5,
        },
    ];
    let arr = k.sketch_regions(&curves).unwrap();
    assert_eq!(arr.regions.len(), 2);
    let ring = k
        .region_containing(&arr, &curves, DVec2::new(0.75, 0.0))
        .unwrap();
    let disc = k
        .region_containing(&arr, &curves, DVec2::new(0.1, 0.1))
        .unwrap();
    assert_ne!(ring, disc);
    let area = |i: usize| arr.regions[i].area;
    assert!((area(ring) - PI * 0.75).abs() < 1e-9);
    assert!((area(disc) - PI * 0.25).abs() < 1e-9);
    assert_eq!(
        k.region_containing(&arr, &curves, DVec2::new(2.0, 0.0)),
        None
    );
    // On a boundary: no region.
    assert_eq!(
        k.region_containing(&arr, &curves, DVec2::new(1.0, 0.0)),
        None
    );
}
