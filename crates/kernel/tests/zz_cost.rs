//! Scratch: the operations whose instructions callgrind counts.
#![allow(missing_docs)]
use std::f64::consts::TAU;

use opendatum_kernel::{BooleanOp, BrepKernel, DVec2, DVec3, ExtrudeSpec, Frame, Kernel, Profile};

fn k() -> BrepKernel {
    BrepKernel::default()
}

fn star(n: usize) -> Vec<DVec2> {
    (0..n)
        .map(|i| {
            let a = i as f64 * TAU / n as f64;
            let r = if i % 2 == 0 { 1.0 } else { 0.9 };
            DVec2::new(a.cos(), a.sin()) * r
        })
        .collect()
}

fn prism_of(corners: &[DVec2], frame: Frame) -> opendatum_kernel::Body {
    k().extrude(
        &Profile::polygon(frame, corners),
        &ExtrudeSpec {
            distance: 1.0,
            symmetric: false,
        },
    )
    .expect("a prism")
}

fn tilted() -> Frame {
    Frame::from_normal(
        DVec3::new(0.3, 0.1, 0.5),
        DVec3::new(0.2, 0.1, 1.0),
        DVec3::X,
    )
}

#[test]
fn zz_cost_star() {
    let (a, b) = (
        prism_of(&star(256), Frame::WORLD),
        prism_of(&star(256), tilted()),
    );
    k().boolean(&a, &b, BooleanOp::Union).expect("a union");
}

#[test]
fn zz_cost_pins() {
    let mut body = k()
        .make_box(&Frame::WORLD, DVec3::new(8.0, 8.0, 1.0))
        .expect("a box");
    for i in 0..8 {
        for j in 0..8 {
            let frame = Frame {
                origin: DVec3::new(i as f64 + 0.3, j as f64 + 0.3, -0.5),
                ..Frame::WORLD
            };
            let pin = k()
                .make_box(&frame, DVec3::new(0.4, 0.4, 2.0))
                .expect("a pin");
            body = k()
                .boolean(&body, &pin, BooleanOp::Subtract)
                .expect("a cut");
        }
    }
}

fn notched_block() -> opendatum_kernel::Body {
    let pts = [
        (-0.25, -1.0),
        (0.0, -1.0),
        (0.0, 0.0),
        (1.0, 0.0),
        (1.0, -1.0),
        (1.25, -1.0),
        (1.25, 1.0),
        (-0.25, 1.0),
    ]
    .map(|(x, y)| DVec2::new(x, y));
    let frame = Frame {
        origin: DVec3::new(0.0, 0.0, -1.0),
        ..Frame::WORLD
    };
    prism_of(&pts, frame)
}

fn zigzag(n: usize) -> opendatum_kernel::Body {
    let eps = k().tolerance().linear;
    let t = 10f64.to_radians();
    let frame = Frame {
        origin: DVec3::ZERO,
        x: DVec3::X,
        y: DVec3::new(0.0, t.cos(), t.sin()),
        z: DVec3::new(0.0, -t.sin(), t.cos()),
    };
    let mut pts = vec![DVec2::new(-0.5, -0.5), DVec2::new(0.0, -0.5)];
    pts.extend((0..=n).map(|i| {
        let w = if i % 2 == 0 { -4.0 } else { -2.0 };
        DVec2::new(i as f64 / n as f64, w * eps)
    }));
    pts.extend([(1.0, -0.5), (2.0, -0.5), (2.0, 0.5), (-0.5, 0.5)].map(|(x, y)| DVec2::new(x, y)));
    k().extrude(
        &Profile::polygon(frame, &pts),
        &ExtrudeSpec {
            distance: 0.5,
            symmetric: false,
        },
    )
    .expect("a prism")
}

#[test]
fn zz_cost_zigzag() {
    let n: usize = std::env::var("ZIG").map_or(256, |s| s.parse().unwrap());
    let (a, b) = (notched_block(), zigzag(n));
    let r = k()
        .boolean(&a, &b, BooleanOp::Intersect)
        .expect("two wedges");
    assert_eq!(k().topology(&r).solids, 2);
}

/// A comb: the zigzag crossing the axis, between w = ±1e-3 m, so the face
/// has n intervals on the line.
fn comb(n: usize) -> opendatum_kernel::Body {
    let t = 10f64.to_radians();
    let frame = Frame {
        origin: DVec3::ZERO,
        x: DVec3::X,
        y: DVec3::new(0.0, t.cos(), t.sin()),
        z: DVec3::new(0.0, -t.sin(), t.cos()),
    };
    let mut pts = vec![DVec2::new(-0.5, -0.5), DVec2::new(0.0, -0.5)];
    pts.extend((0..=n).map(|i| {
        let w = if i % 2 == 0 { -1e-3 } else { 1e-3 };
        DVec2::new(i as f64 / n as f64, w)
    }));
    pts.extend([(1.0, -0.5), (2.0, -0.5), (2.0, 0.5), (-0.5, 0.5)].map(|(x, y)| DVec2::new(x, y)));
    k().extrude(
        &Profile::polygon(frame, &pts),
        &ExtrudeSpec {
            distance: 0.5,
            symmetric: false,
        },
    )
    .expect("a prism")
}

#[test]
fn zz_cost_comb() {
    let n: usize = std::env::var("ZIG").map_or(256, |s| s.parse().unwrap());
    let a = k()
        .make_box(
            &Frame {
                origin: DVec3::new(-0.25, -1.0, -1.0),
                ..Frame::WORLD
            },
            DVec3::new(1.5, 2.0, 1.0),
        )
        .expect("a box");
    let b = comb(n);
    let t = std::time::Instant::now();
    let r = k().boolean(&a, &b, BooleanOp::Intersect);
    eprintln!(
        "COMB {n}: {:?} {:.2} ms",
        r.as_ref().map(|r| k().topology(r).solids),
        t.elapsed().as_secs_f64() * 1e3
    );
}
