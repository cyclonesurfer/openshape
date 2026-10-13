//! Time guards for the kernel's user-facing operations (docs/PERFORMANCE.md).
//!
//! Each operation is timed at a reference size (the best of a few runs) and
//! held to a limit about ten times its measured time, so ordinary noise and
//! slower machines pass while an accidental quadratic loop or a lost index
//! fails. The scaling checks time one operation at two sizes and bound the
//! ratio, which holds on any machine. Limits apply to release builds, which
//! CI runs (`cargo test --release -p opendatum-kernel --test timing`); a
//! debug build only reports the times.

use std::f64::consts::TAU;
use std::time::{Duration, Instant};

use opendatum_kernel::{
    Body, BooleanOp, BrepKernel, DVec2, DVec3, ExtrudeSpec, Frame, Kernel, Profile, ProfileCurve,
    RevolveSpec, TessellationOptions,
};

fn k() -> BrepKernel {
    BrepKernel::default()
}

/// The fastest of `runs` runs of `f`, and its last result.
fn best<T>(runs: usize, mut f: impl FnMut() -> T) -> (T, Duration) {
    let mut fastest = Duration::MAX;
    let mut out = None;
    for _ in 0..runs {
        let t0 = Instant::now();
        let r = f();
        fastest = fastest.min(t0.elapsed());
        out = Some(r);
    }
    (out.expect("at least one run"), fastest)
}

fn ms(t: Duration) -> f64 {
    t.as_secs_f64() * 1e3
}

/// Reports `t` and, in a release build, holds it to `limit_ms`.
fn within(label: &str, t: Duration, limit_ms: f64) {
    println!("{label}: {:.2} ms (limit {limit_ms} ms)", ms(t));
    if !cfg!(debug_assertions) {
        assert!(
            ms(t) <= limit_ms,
            "{label} took {:.2} ms, over its limit of {limit_ms} ms (docs/PERFORMANCE.md)",
            ms(t)
        );
    }
}

/// Times `small` and `large` alternately, `runs` times each, and reports
/// how the fastest time grows from the small input to the large one,
/// `factor` times its size. In a release build the growth is held to three
/// times `factor`: an n log n operation grows a little more than `factor`
/// (about 11 times for eightfold), a quadratic one `factor` squared (64).
/// Alternating the runs puts both sizes through the same machine noise.
fn scales<A, B>(
    label: &str,
    factor: f64,
    runs: usize,
    mut small: impl FnMut() -> A,
    mut large: impl FnMut() -> B,
) {
    let (mut s, mut l) = (Duration::MAX, Duration::MAX);
    for _ in 0..runs {
        let t0 = Instant::now();
        std::hint::black_box(small());
        s = s.min(t0.elapsed());
        let t0 = Instant::now();
        std::hint::black_box(large());
        l = l.min(t0.elapsed());
    }
    let ratio = l.as_secs_f64() / s.as_secs_f64().max(1e-9);
    let bound = 3.0 * factor;
    println!(
        "{label}: {:.2} ms → {:.2} ms at {factor} times the size, {ratio:.1} times the time \
         (bound {bound})",
        ms(s),
        ms(l)
    );
    if !cfg!(debug_assertions) {
        assert!(
            ratio <= bound,
            "{label}: time grew {ratio:.1} times for {factor} times the size, over {bound} \
             (docs/PERFORMANCE.md)"
        );
    }
}

/// A star with `n` corners alternating between radii 1 and 0.9.
fn star(n: usize) -> Vec<DVec2> {
    (0..n)
        .map(|i| {
            let a = i as f64 * TAU / n as f64;
            let r = if i % 2 == 0 { 1.0 } else { 0.9 };
            DVec2::new(a.cos(), a.sin()) * r
        })
        .collect()
}

/// A regular polygon with `n` corners on the unit circle.
fn regular(n: usize) -> Vec<DVec2> {
    (0..n)
        .map(|i| {
            let a = i as f64 * TAU / n as f64;
            DVec2::new(a.cos(), a.sin())
        })
        .collect()
}

fn prism_of(corners: &[DVec2], frame: Frame) -> Body {
    k().extrude(
        &Profile::polygon(frame, corners),
        &ExtrudeSpec {
            distance: 1.0,
            symmetric: false,
        },
    )
    .expect("a prism")
}

fn prism(n: usize, frame: Frame) -> Body {
    prism_of(&star(n), frame)
}

fn tilted() -> Frame {
    Frame::from_normal(
        DVec3::new(0.3, 0.1, 0.5),
        DVec3::new(0.2, 0.1, 1.0),
        DVec3::X,
    )
}

/// A plate with an `m` by `m` grid of round holes.
fn plate(m: usize) -> Body {
    let mut p = Profile::rectangle(Frame::WORLD, DVec2::ZERO, DVec2::splat(m as f64));
    for i in 0..m {
        for j in 0..m {
            p = p.with(Profile::circle(
                Frame::WORLD,
                DVec2::new(i as f64 + 0.5, j as f64 + 0.5),
                0.3,
            ));
        }
    }
    k().extrude(
        &p,
        &ExtrudeSpec {
            distance: 0.2,
            symmetric: false,
        },
    )
    .expect("a plate")
}

/// A ring with a rounded outside, revolved a full turn.
fn ring() -> Body {
    let p = Profile {
        plane: Frame::WORLD,
        curves: vec![
            ProfileCurve::Line {
                start: DVec2::new(0.5, 0.0),
                end: DVec2::new(1.0, 0.0),
            },
            ProfileCurve::Arc {
                center: DVec2::new(1.0, 0.25),
                radius: 0.25,
                start_angle: -TAU / 4.0,
                sweep_angle: TAU / 2.0,
            },
            ProfileCurve::Line {
                start: DVec2::new(1.0, 0.5),
                end: DVec2::new(0.5, 0.5),
            },
            ProfileCurve::Line {
                start: DVec2::new(0.5, 0.5),
                end: DVec2::new(0.5, 0.0),
            },
        ],
    };
    k().revolve(
        &p,
        &RevolveSpec {
            axis_origin: DVec2::ZERO,
            axis_dir: DVec2::Y,
            angle: TAU,
        },
    )
    .expect("a ring")
}

#[test]
fn solids_meshes_and_mass() {
    let opts = TessellationOptions::default();
    let (body, t) = best(3, || prism(1024, Frame::WORLD));
    within("extrude a 1024-sided prism", t, LIMIT.extrude_1024);
    let (_, t) = best(3, || k().tessellate(&body, &opts).expect("a mesh"));
    within("mesh the 1024-sided prism", t, LIMIT.mesh_1024);
    let (_, t) = best(3, || k().validate(&body));
    within("validate the 1024-sided prism", t, LIMIT.validate_1024);
    let (_, t) = best(2, || k().mass_properties(&body, 1e-9).expect("mass"));
    within(
        "mass properties of the 1024-sided prism (1e-9)",
        t,
        LIMIT.mass_1024,
    );

    let ring = ring();
    let (m, t) = best(3, || k().tessellate(&ring, &opts).expect("a mesh"));
    within(
        &format!("mesh a revolved ring ({} triangles)", m.triangle_count()),
        t,
        LIMIT.mesh_ring,
    );

    let (p, t) = best(2, plate10);
    within("extrude a plate with 100 holes", t, LIMIT.extrude_plate);
    let (_, t) = best(2, || k().tessellate(&p, &opts).expect("a mesh"));
    within("mesh the plate with 100 holes", t, LIMIT.mesh_plate);
}

fn plate10() -> Body {
    plate(10)
}

#[test]
fn booleans() {
    let (a, b) = (prism(256, Frame::WORLD), prism(256, tilted()));
    for op in [BooleanOp::Union, BooleanOp::Subtract] {
        let (_, t) = best(3, || k().boolean(&a, &b, op).expect("a boolean"));
        within(
            &format!("{op:?} of two 256-sided prisms"),
            t,
            LIMIT.boolean_256,
        );
    }
    // A box with an 8 by 8 pattern of pins cut one at a time: 64 booleans
    // on a body that grows to 325 faces.
    let (body, t) = best(1, || {
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
        body
    });
    assert_eq!(k().topology(&body).faces, 6 + 64 * 4);
    within("cut 64 holes into a box one at a time", t, LIMIT.pins_64);
}

#[test]
fn sketch_regions() {
    // A grid of 41 by 41 lines and 40 circles: 1,681 crossings among the
    // lines alone.
    let m = 40;
    let mut curves = Vec::new();
    for i in 0..=m {
        let x = i as f64;
        curves.push(ProfileCurve::Line {
            start: DVec2::new(x, -0.5),
            end: DVec2::new(x, m as f64 + 0.5),
        });
        curves.push(ProfileCurve::Line {
            start: DVec2::new(-0.5, x),
            end: DVec2::new(m as f64 + 0.5, x),
        });
    }
    for i in 0..m {
        curves.push(ProfileCurve::Circle {
            center: DVec2::new(i as f64 + 0.5, i as f64 + 0.5),
            radius: 0.7,
        });
    }
    let (_, t) = best(2, || k().sketch_regions(&curves).expect("regions"));
    within(
        &format!("sketch regions of {} curves", curves.len()),
        t,
        LIMIT.regions_122,
    );
}

#[test]
fn scaling() {
    let opts = TessellationOptions::default();
    scales(
        "extrude a prism, 256 to 2048 sides",
        8.0,
        5,
        || prism(256, Frame::WORLD),
        || prism(2048, Frame::WORLD),
    );

    // Regular prisms: the result grows in proportion to the input. (Two
    // star prisms, as in `booleans`, overlap in a band where every point of
    // one meets many of the other, so their result grows faster than their
    // size.)
    let small = (
        prism_of(&regular(128), Frame::WORLD),
        prism_of(&regular(128), tilted()),
    );
    let large = (
        prism_of(&regular(1024), Frame::WORLD),
        prism_of(&regular(1024), tilted()),
    );
    scales(
        "union of two regular prisms, 128 to 1024 sides",
        8.0,
        5,
        || {
            k().boolean(&small.0, &small.1, BooleanOp::Union)
                .expect("a union")
        },
        || {
            k().boolean(&large.0, &large.1, BooleanOp::Union)
                .expect("a union")
        },
    );

    let (p4, p11) = (plate(4), plate(11));
    scales(
        "mesh a plate, 16 to 121 holes",
        121.0 / 16.0,
        5,
        || k().tessellate(&p4, &opts).expect("a mesh"),
        || k().tessellate(&p11, &opts).expect("a mesh"),
    );
}

/// Release-build limits in milliseconds, about ten times the times measured
/// on the reference machine (docs/PERFORMANCE.md).
struct Limits {
    extrude_1024: f64,
    mesh_1024: f64,
    validate_1024: f64,
    mass_1024: f64,
    mesh_ring: f64,
    extrude_plate: f64,
    mesh_plate: f64,
    boolean_256: f64,
    pins_64: f64,
    regions_122: f64,
}

const LIMIT: Limits = Limits {
    extrude_1024: 120.0,
    mesh_1024: 100.0,
    validate_1024: 100.0,
    mass_1024: 700.0,
    mesh_ring: 3000.0,
    extrude_plate: 150.0,
    mesh_plate: 1500.0,
    boolean_256: 500.0,
    pins_64: 3000.0,
    regions_122: 1200.0,
};
