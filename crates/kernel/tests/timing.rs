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

/// A cylinder of radius `r` about `(x, y)`, from `z0` up `h`.
fn pin(x: f64, y: f64, z0: f64, r: f64, h: f64) -> Body {
    let frame = Frame {
        origin: DVec3::new(0.0, 0.0, z0),
        ..Frame::WORLD
    };
    k().extrude(
        &Profile::circle(frame, DVec2::new(x, y), r),
        &ExtrudeSpec {
            distance: h,
            symmetric: false,
        },
    )
    .expect("a pin")
}

/// Half the thickness of `plate(m)`'s plate and more: a box over its upper
/// half that cuts every hole across.
fn upper(m: usize) -> Body {
    k().make_box(
        &Frame {
            origin: DVec3::new(-0.5, -0.5, 0.1),
            ..Frame::WORLD
        },
        DVec3::new(m as f64 + 1.0, m as f64 + 1.0, 1.0),
    )
    .expect("a box")
}

#[test]
fn curved_booleans() {
    // A plate drilled with an 8 by 8 grid of holes, one boolean each: the
    // body grows to 134 faces, two half cylinders a hole (ADR 0008).
    let (body, t) = best(1, || {
        let mut body = k()
            .make_box(&Frame::WORLD, DVec3::new(8.0, 8.0, 0.2))
            .expect("a plate");
        for i in 0..8 {
            for j in 0..8 {
                let hole = pin(i as f64 + 0.5, j as f64 + 0.5, -0.5, 0.3, 1.2);
                body = k()
                    .boolean(&body, &hole, BooleanOp::Subtract)
                    .expect("a hole");
            }
        }
        body
    });
    assert_eq!(k().topology(&body).faces, 6 + 64 * 2);
    within(
        "drill 64 holes into a plate one at a time",
        t,
        LIMIT.drill_64,
    );

    // A plate with 100 holes and a box over its upper half, across them.
    let (p, top) = (plate(10), upper(10));
    for op in [BooleanOp::Intersect, BooleanOp::Subtract] {
        let (_, t) = best(3, || k().boolean(&p, &top, op).expect("a boolean"));
        within(
            &format!("{op:?} a plate with 100 holes and a box across them"),
            t,
            LIMIT.plate_box,
        );
    }
}

/// The vertical plane through `origin` holding the world's z axis as its
/// y axis: profiles in it revolve about the vertical.
fn upright(origin: DVec3) -> Frame {
    Frame {
        origin,
        x: DVec3::X,
        y: DVec3::Z,
        z: -DVec3::Y,
    }
}

/// A whole turn about the vertical.
fn turn() -> RevolveSpec {
    RevolveSpec {
        axis_origin: DVec2::ZERO,
        axis_dir: DVec2::Y,
        angle: TAU,
    }
}

/// A ball of radius `r` about `centre`.
fn ball(centre: DVec3, r: f64) -> Body {
    let p = Profile {
        plane: upright(centre),
        curves: vec![
            ProfileCurve::Arc {
                center: DVec2::ZERO,
                radius: r,
                start_angle: -TAU / 4.0,
                sweep_angle: TAU / 2.0,
            },
            ProfileCurve::Line {
                start: DVec2::new(0.0, r),
                end: DVec2::new(0.0, -r),
            },
        ],
    };
    k().revolve(&p, &turn()).expect("a ball")
}

/// The tool of a countersunk hole about the vertical through `(x, y)` in
/// a plate 0.4 thick: a shank of radius 0.15 and a 90° countersink 0.1
/// deep at the plate's top.
fn countersink(x: f64, y: f64) -> Body {
    let pts = [
        (0.0, -0.5),
        (0.15, -0.5),
        (0.15, 0.25),
        (0.35, 0.45),
        (0.0, 0.45),
    ]
    .map(|(r, h)| DVec2::new(r, h));
    k().revolve(
        &Profile::polygon(upright(DVec3::new(x, y, 0.0)), &pts),
        &turn(),
    )
    .expect("a countersink")
}

#[test]
fn spheres_and_cones() {
    let plate = || {
        k().make_box(&Frame::WORLD, DVec3::new(4.0, 4.0, 0.4))
            .expect("a plate")
    };
    // A plate dimpled by a 4 by 4 grid of balls, one boolean each.
    let (dimpled, t) = best(1, || {
        let mut body = plate();
        for i in 0..4 {
            for j in 0..4 {
                let b = ball(DVec3::new(i as f64 + 0.5, j as f64 + 0.5, 0.55), 0.3);
                body = k()
                    .boolean(&body, &b, BooleanOp::Subtract)
                    .expect("a dimple");
            }
        }
        body
    });
    assert_eq!(k().topology(&dimpled).faces, 6 + 16 * 2);
    within(
        "dimple a plate 16 times with balls, one boolean each",
        t,
        LIMIT.dimple_16,
    );
    // The same plate with 16 countersunk holes.
    let (sunk, t) = best(1, || {
        let mut body = plate();
        for i in 0..4 {
            for j in 0..4 {
                let tool = countersink(i as f64 + 0.5, j as f64 + 0.5);
                body = k()
                    .boolean(&body, &tool, BooleanOp::Subtract)
                    .expect("a countersunk hole");
            }
        }
        body
    });
    let topology = k().topology(&sunk);
    assert_eq!((topology.faces, topology.genus), (6 + 16 * 4, Some(16)));
    within(
        "cut 16 countersunk holes into a plate, one boolean each",
        t,
        LIMIT.sink_16,
    );
    let opts = TessellationOptions::default();
    for (name, body, [mesh, validate, mass]) in [
        ("dimpled", &dimpled, LIMIT.dimpled),
        ("countersunk", &sunk, LIMIT.countersunk),
    ] {
        let (m, t) = best(2, || k().tessellate(body, &opts).expect("a mesh"));
        within(
            &format!("mesh the {name} plate ({} triangles)", m.triangle_count()),
            t,
            mesh,
        );
        let (_, t) = best(2, || k().validate(body));
        within(&format!("validate the {name} plate"), t, validate);
        let (_, t) = best(2, || k().mass_properties(body, 1e-9).expect("mass"));
        within(
            &format!("mass properties of the {name} plate (1e-9)"),
            t,
            mass,
        );
    }
    let b = ball(DVec3::ZERO, 1.0);
    let (m, t) = best(2, || k().tessellate(&b, &opts).expect("a mesh"));
    within(
        &format!("mesh a unit ball ({} triangles)", m.triangle_count()),
        t,
        LIMIT.mesh_ball,
    );
}

/// The plane through `origin` whose y axis leans `lean` radians from the
/// vertical toward +x: profiles in it revolve about that line.
fn leaning(origin: DVec3, lean: f64) -> Frame {
    let (s, c) = lean.sin_cos();
    Frame {
        origin,
        x: DVec3::new(c, 0.0, -s),
        y: DVec3::new(s, 0.0, c),
        z: -DVec3::Y,
    }
}

/// The tool of a countersunk hole like [`countersink`]'s, its axis leaning
/// 10°, its cone reaching higher: the plate's top face meets the cone in
/// an ellipse.
fn tilted_countersink(x: f64, y: f64) -> Body {
    let pts = [
        (0.0, -0.5),
        (0.15, -0.5),
        (0.15, 0.25),
        (0.45, 0.55),
        (0.0, 0.55),
    ]
    .map(|(r, h)| DVec2::new(r, h));
    k().revolve(
        &Profile::polygon(leaning(DVec3::new(x, y, 0.0), 10f64.to_radians()), &pts),
        &turn(),
    )
    .expect("a tilted countersink")
}

/// A pin of 5° draft about the vertical through `(x, y)`, radius 0.3 at
/// its foot 0.2 up, inside a plate 0.4 thick, to 1.4 up.
fn drafted_pin(x: f64, y: f64) -> Body {
    let top = 0.3 - 1.2 * 5f64.to_radians().tan();
    let pts = [(0.0, 0.2), (0.3, 0.2), (top, 1.4), (0.0, 1.4)].map(|(r, h)| DVec2::new(r, h));
    k().revolve(
        &Profile::polygon(upright(DVec3::new(x, y, 0.0)), &pts),
        &turn(),
    )
    .expect("a drafted pin")
}

#[test]
fn planes_across_cones() {
    let plate = || {
        k().make_box(&Frame::WORLD, DVec3::new(4.0, 4.0, 0.4))
            .expect("a plate")
    };
    let grid = || (0..16).map(|n| ((n / 4) as f64 + 0.5, (n % 4) as f64 + 0.5));
    // Disjoint solids made one body, which a body with conic edges on
    // its curved faces cannot be an operand of yet (STATUS.md).
    let all = |solids: Vec<Body>| {
        solids
            .into_iter()
            .reduce(|a, b| k().boolean(&a, &b, BooleanOp::Union).expect("a union"))
            .expect("solids")
    };
    // 16 countersunk holes leaning 10°, in one boolean: ellipses where the
    // plate's faces cross the cones and the shanks.
    let tools = all(grid().map(|(x, y)| tilted_countersink(x, y)).collect());
    let (sunk, t) = best(3, || {
        k().boolean(&plate(), &tools, BooleanOp::Subtract)
            .expect("tilted countersunk holes")
    });
    let topology = k().topology(&sunk);
    assert_eq!((topology.faces, topology.genus), (6 + 16 * 4, Some(16)));
    within(
        "cut 16 countersunk holes leaning 10° into a plate in one boolean",
        t,
        LIMIT.tilted_sink_16,
    );
    // 16 pins of 5° draft on a plate, flattened on one side by boxes whose
    // faces run along their axes, in one boolean: hyperbolas.
    let pins = grid().fold(plate(), |body, (x, y)| {
        k().boolean(&body, &drafted_pin(x, y), BooleanOp::Union)
            .expect("a pin")
    });
    let boxes = all(grid()
        .map(|(x, y)| {
            k().make_box(
                &Frame {
                    origin: DVec3::new(x + 0.1, y - 0.4, 0.5),
                    ..Frame::WORLD
                },
                DVec3::new(0.3, 0.8, 1.5),
            )
            .expect("a box")
        })
        .collect());
    let (flat, t) = best(3, || {
        k().boolean(&pins, &boxes, BooleanOp::Subtract)
            .expect("flattened pins")
    });
    let topology = k().topology(&flat);
    assert_eq!((topology.solids, topology.genus), (1, Some(0)));
    within(
        "flatten 16 drafted pins along their axes in one boolean",
        t,
        LIMIT.flatten_16,
    );
    let opts = TessellationOptions::default();
    for (name, body, [mesh, validate, mass]) in [
        ("tilted countersunk", &sunk, LIMIT.tilted),
        ("flattened pins'", &flat, LIMIT.flattened),
    ] {
        let (m, t) = best(2, || k().tessellate(body, &opts).expect("a mesh"));
        within(
            &format!("mesh the {name} plate ({} triangles)", m.triangle_count()),
            t,
            mesh,
        );
        let (_, t) = best(2, || k().validate(body));
        within(&format!("validate the {name} plate"), t, validate);
        let (_, t) = best(2, || k().mass_properties(body, 1e-9).expect("mass"));
        within(
            &format!("mass properties of the {name} plate (1e-9)"),
            t,
            mass,
        );
    }
}

/// A solid turned from one rational cubic whose middle weights are far
/// apart (1, `w`, 0.2, 1): it bends nearly to a crease.
fn creased(w: f64) -> Body {
    let pts = [(0.2, 0.0), (1.2, 0.2), (1.2, 1.8), (0.2, 2.0)].map(|(x, y)| DVec2::new(x, y));
    k().revolve(
        &Profile {
            plane: Frame {
                origin: DVec3::ZERO,
                x: DVec3::X,
                y: DVec3::Z,
                z: -DVec3::Y,
            },
            curves: vec![
                ProfileCurve::Spline {
                    degree: 3,
                    knots: vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0],
                    points: pts.to_vec(),
                    weights: Some(vec![1.0, w, 0.2, 1.0]),
                },
                ProfileCurve::Line {
                    start: pts[3],
                    end: pts[0],
                },
            ],
        },
        &turn(),
    )
    .expect("a revolve")
}

/// A NURBS surface: a revolved rational cubic bending nearly to a crease,
/// built (validation meshes it), meshed at the defaults and finely, and
/// validated (found by review: its meshing took minutes).
#[test]
fn a_creased_rational_revolve() {
    let (body, t) = best(3, || creased(60.0));
    within("revolve a creased rational cubic", t, LIMIT.creased[0]);
    for (name, opts, limit) in [
        ("", TessellationOptions::default(), LIMIT.creased[1]),
        (
            " at 1e-3 and 5°",
            TessellationOptions {
                chord_tolerance: 1e-3,
                max_angle: 5f64.to_radians(),
            },
            LIMIT.creased[2],
        ),
    ] {
        let (m, t) = best(2, || k().tessellate(&body, &opts).expect("a mesh"));
        within(
            &format!("mesh it{name} ({} triangles)", m.triangle_count()),
            t,
            limit,
        );
    }
    let (_, t) = best(3, || k().validate(&body));
    within("validate it", t, LIMIT.creased[3]);
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
    let (top4, top11) = (upper(4), upper(11));
    scales(
        "cut a box across a plate's holes, 16 to 121 holes",
        121.0 / 16.0,
        5,
        || k().boolean(&p4, &top4, BooleanOp::Subtract).expect("a cut"),
        || {
            k().boolean(&p11, &top11, BooleanOp::Subtract)
                .expect("a cut")
        },
    );
    scales(
        "mesh a plate, 16 to 121 holes",
        121.0 / 16.0,
        5,
        || k().tessellate(&p4, &opts).expect("a mesh"),
        || k().tessellate(&p11, &opts).expect("a mesh"),
    );

    // A ball meshed eight times as finely: on a doubly curved face the
    // triangles grow as the tolerance shrinks.
    let b = ball(DVec3::ZERO, 1.0);
    let fine = |chord_tolerance: f64| TessellationOptions {
        chord_tolerance,
        ..opts
    };
    let (coarse, finer) = (fine(1e-3), fine(1.25e-4));
    let count =
        |o: &TessellationOptions| k().tessellate(&b, o).expect("a mesh").triangle_count() as f64;
    let (n, m) = (count(&coarse), count(&finer));
    scales(
        &format!("mesh a ball, {n} to {m} triangles"),
        m / n,
        5,
        || k().tessellate(&b, &coarse).expect("a mesh"),
        || k().tessellate(&b, &finer).expect("a mesh"),
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
    drill_64: f64,
    plate_box: f64,
    dimple_16: f64,
    sink_16: f64,
    /// Meshing, validation and mass properties of the dimpled plate.
    dimpled: [f64; 3],
    /// The same of the countersunk plate.
    countersunk: [f64; 3],
    tilted_sink_16: f64,
    flatten_16: f64,
    /// Meshing, validation and mass properties of the plate with tilted
    /// countersunk holes.
    tilted: [f64; 3],
    /// The same of the plate of flattened pins.
    flattened: [f64; 3],
    mesh_ball: f64,
    /// Revolving, meshing at the defaults and at 1e-3 and 5°, and
    /// validating the creased rational revolve.
    creased: [f64; 4],
    regions_122: f64,
}

const LIMIT: Limits = Limits {
    extrude_1024: 120.0,
    mesh_1024: 100.0,
    validate_1024: 100.0,
    mass_1024: 700.0,
    mesh_ring: 3000.0,
    extrude_plate: 150.0,
    mesh_plate: 2500.0,
    boolean_256: 500.0,
    pins_64: 3000.0,
    drill_64: 10000.0,
    plate_box: 600.0,
    dimple_16: 1200.0,
    sink_16: 2500.0,
    dimpled: [4000.0, 80.0, 150.0],
    countersunk: [300.0, 150.0, 250.0],
    tilted_sink_16: 450.0,
    flatten_16: 400.0,
    tilted: [600.0, 150.0, 200.0],
    flattened: [350.0, 100.0, 150.0],
    mesh_ball: 3500.0,
    creased: [150.0, 5000.0, 4000.0, 150.0],
    regions_122: 1200.0,
};
