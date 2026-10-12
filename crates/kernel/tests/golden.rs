//! Golden corpus: reference models with known volume, area, topology and
//! validity. Every model is also checked for a watertight, consistently
//! oriented tessellation.

use std::collections::HashMap;
use std::f64::consts::{PI, TAU};

use openshape_kernel::{
    BodyMesh, DVec2, DVec3, ExtrudeSpec, Frame, Kernel, OsKernel, Profile, ProfileCurve,
    RevolveSpec, SurfaceKind, TessellationOptions,
};
use proptest::prelude::*;

fn k() -> OsKernel {
    OsKernel::default()
}

/// Asserts the welded mesh is closed and consistently oriented: every
/// directed edge appears exactly once and its reverse exactly once.
fn assert_watertight(mesh: &BodyMesh) {
    let (pos, tris) = mesh.welded();
    let mut directed: HashMap<(u32, u32), u32> = HashMap::new();
    for t in &tris {
        for (a, b) in [(t[0], t[1]), (t[1], t[2]), (t[2], t[0])] {
            assert_ne!(a, b, "degenerate triangle");
            *directed.entry((a, b)).or_default() += 1;
        }
    }
    for (&(a, b), &n) in &directed {
        assert_eq!(n, 1, "directed edge {a}->{b} used {n} times");
        assert_eq!(
            directed.get(&(b, a)),
            Some(&1),
            "edge {a}->{b} at {:?} has no opposite (mesh not closed)",
            pos[a as usize]
        );
    }
}

/// Total polyline length of the seam edges and of the crease edges.
fn seam_and_crease_length(body: &openshape_kernel::Body) -> (f64, f64) {
    let opts = TessellationOptions {
        chord_tolerance: 1e-5,
        max_angle: 2f64.to_radians(),
    };
    let mesh = k().tessellate(body, &opts).unwrap();
    assert_eq!(mesh.seams.len(), mesh.edges.len());
    let mut out = (0.0, 0.0);
    for (e, seam) in mesh.edges.iter().zip(&mesh.seams) {
        let l: f64 = e.windows(2).map(|w| (w[1] - w[0]).length()).sum();
        if *seam {
            out.0 += l;
        } else {
            out.1 += l;
        }
    }
    out
}

fn assert_close(a: f64, b: f64, rel: f64, what: &str) {
    assert!(
        (a - b).abs() <= rel * b.abs().max(1e-12),
        "{what}: {a} vs {b}"
    );
}

struct Expect {
    volume: f64,
    area: f64,
    counts: Option<(usize, usize, usize)>, // V, E, F
    genus: usize,
    solids: usize,
    rel: f64,
}

fn check(body: &openshape_kernel::Body, e: Expect) {
    let kern = k();
    let issues = kern.validate(body);
    assert!(issues.is_empty(), "invalid body: {issues:?}");
    let t = kern.topology(body);
    if let Some((v, ed, f)) = e.counts {
        assert_eq!(
            (t.vertices, t.edges, t.faces),
            (v, ed, f),
            "topology counts"
        );
    }
    assert_eq!(t.genus, Some(e.genus), "genus");
    assert_eq!(t.solids, e.solids, "solids");
    let m = kern.mass_properties(body, 1e-4).unwrap();
    assert!(
        ((m.volume - e.volume) / e.volume).abs() < e.rel,
        "volume {} vs {}",
        m.volume,
        e.volume
    );
    assert!(
        ((m.area - e.area) / e.area).abs() < e.rel,
        "area {} vs {}",
        m.area,
        e.area
    );
    for opts in [
        TessellationOptions::default(),
        TessellationOptions {
            chord_tolerance: 1e-3,
            max_angle: 0.5,
        },
    ] {
        let mesh = kern.tessellate(body, &opts).unwrap();
        assert_watertight(&mesh);
    }
}

#[test]
fn unit_box() {
    let b = k()
        .make_box(&Frame::WORLD, DVec3::new(1.0, 2.0, 3.0))
        .unwrap();
    check(
        &b,
        Expect {
            volume: 6.0,
            area: 22.0,
            counts: Some((8, 12, 6)),
            genus: 0,
            solids: 1,
            rel: 1e-12,
        },
    );
    let m = k().mass_properties(&b, 1e-3).unwrap();
    assert!((m.centroid - DVec3::new(0.5, 1.0, 1.5)).length() < 1e-12);
    // Inertia of a cuboid about its centroid: m (b² + c²) / 12 etc.
    let mass = 6.0;
    let expect = [
        mass * (4.0 + 9.0) / 12.0,
        mass * (1.0 + 9.0) / 12.0,
        mass * (1.0 + 4.0) / 12.0,
    ];
    for (i, e) in expect.iter().enumerate() {
        assert!((m.inertia[i][i] - e).abs() < 1e-9, "I{i}{i}");
        for j in 0..3 {
            if i != j {
                assert!(m.inertia[i][j].abs() < 1e-9);
            }
        }
    }
    assert!(k().face_kinds(&b).iter().all(|k| *k == SurfaceKind::Plane));
    let mesh = k().tessellate(&b, &TessellationOptions::default()).unwrap();
    assert!(mesh.seams.iter().all(|s| !s), "a box has no seams");
}

#[test]
fn box_with_negative_size_and_offset_frame() {
    let f = Frame::from_normal(
        DVec3::new(1.0, 2.0, 3.0),
        DVec3::new(1.0, 1.0, 0.0),
        DVec3::Z,
    );
    let b = k().make_box(&f, DVec3::new(-0.1, 0.2, -0.3)).unwrap();
    check(
        &b,
        Expect {
            volume: 0.006,
            area: 2.0 * (0.02 + 0.03 + 0.06),
            counts: Some((8, 12, 6)),
            genus: 0,
            solids: 1,
            rel: 1e-10,
        },
    );
}

#[test]
fn degenerate_box_is_rejected() {
    assert!(k()
        .make_box(&Frame::WORLD, DVec3::new(1.0, 0.0, 1.0))
        .is_err());
    assert!(k()
        .make_box(&Frame::WORLD, DVec3::new(1.0, f64::NAN, 1.0))
        .is_err());
}

#[test]
fn cylinder_extrude() {
    let r = 0.5;
    let h = 2.0;
    let b = k()
        .extrude(
            &Profile::circle(Frame::WORLD, DVec2::new(1.0, -1.0), r),
            &ExtrudeSpec {
                distance: h,
                symmetric: false,
            },
        )
        .unwrap();
    check(
        &b,
        Expect {
            volume: PI * r * r * h,
            area: 2.0 * PI * r * r + 2.0 * PI * r * h,
            counts: Some((4, 6, 4)),
            genus: 0,
            solids: 1,
            rel: 1e-9,
        },
    );
    // The two lines where the side is split are seams; the circles are
    // creases.
    let (seam, crease) = seam_and_crease_length(&b);
    assert_close(seam, 2.0 * h, 1e-12, "seam length");
    assert_close(crease, 2.0 * TAU * r, 1e-4, "crease length");
    let kinds = k().face_kinds(&b);
    assert_eq!(
        kinds
            .iter()
            .filter(|k| **k == SurfaceKind::Cylinder)
            .count(),
        2
    );
    // Chord tolerance is honoured: every vertex of the side faces lies on the
    // cylinder and every triangle centroid is within tolerance of it.
    let tol = 1e-3;
    let mesh = k()
        .tessellate(
            &b,
            &TessellationOptions {
                chord_tolerance: tol,
                max_angle: 1.0,
            },
        )
        .unwrap();
    for (fi, f) in mesh.faces.iter().enumerate() {
        if kinds[fi] != SurfaceKind::Cylinder {
            continue;
        }
        for t in &f.triangles {
            let c = (f.positions[t[0] as usize]
                + f.positions[t[1] as usize]
                + f.positions[t[2] as usize])
                / 3.0;
            let d = (c.truncate() - glam::DVec2::new(1.0, -1.0)).length();
            assert!(
                d <= r + 1e-12 && d >= r - tol - 1e-12,
                "centroid off by {}",
                r - d
            );
        }
    }
}

#[test]
fn plate_with_hole_has_genus_one() {
    let p = Profile::rectangle(Frame::WORLD, DVec2::ZERO, DVec2::new(4.0, 3.0))
        .with(Profile::circle(Frame::WORLD, DVec2::new(2.0, 1.5), 0.5));
    let b = k()
        .extrude(
            &p,
            &ExtrudeSpec {
                distance: 0.5,
                symmetric: true,
            },
        )
        .unwrap();
    let area2d = 12.0 - PI * 0.25;
    check(
        &b,
        Expect {
            volume: area2d * 0.5,
            area: 2.0 * area2d + 14.0 * 0.5 + PI * 0.5,
            counts: Some((12, 18, 8)),
            genus: 1,
            solids: 1,
            rel: 1e-9,
        },
    );
    let (seam, crease) = seam_and_crease_length(&b);
    assert_close(seam, 2.0 * 0.5, 1e-12, "the hole's seams");
    assert_close(
        crease,
        2.0 * 14.0 + 4.0 * 0.5 + 2.0 * TAU * 0.5,
        1e-4,
        "creases",
    );
    let mesh = k().tessellate(&b, &TessellationOptions::default()).unwrap();
    let (lo, hi) = mesh.bounds().unwrap();
    assert!(
        (lo.z + 0.25).abs() < 1e-12 && (hi.z - 0.25).abs() < 1e-12,
        "symmetric"
    );
}

#[test]
fn downward_extrude_and_island_in_hole() {
    // Outer square with a square hole containing an island: two solids.
    let p = Profile::rectangle(Frame::WORLD, DVec2::splat(-3.0), DVec2::splat(3.0))
        .with(Profile::rectangle(
            Frame::WORLD,
            DVec2::splat(-2.0),
            DVec2::splat(2.0),
        ))
        .with(Profile::circle(Frame::WORLD, DVec2::ZERO, 1.0));
    let b = k()
        .extrude(
            &p,
            &ExtrudeSpec {
                distance: -1.0,
                symmetric: false,
            },
        )
        .unwrap();
    check(
        &b,
        Expect {
            volume: (36.0 - 16.0) + PI,
            area: 2.0 * (20.0 + PI) + 24.0 + 16.0 + 2.0 * PI,
            counts: None,
            genus: 1,
            solids: 2,
            rel: 1e-9,
        },
    );
}

#[test]
fn rounded_slot_from_lines_and_arcs() {
    // A slot: two lines and two half-circle arcs, given in mixed orientation.
    let p = Profile {
        plane: Frame::front(),
        curves: vec![
            ProfileCurve::Line {
                start: DVec2::new(0.0, 0.0),
                end: DVec2::new(2.0, 0.0),
            },
            ProfileCurve::Arc {
                center: DVec2::new(2.0, 0.5),
                radius: 0.5,
                start_angle: -PI / 2.0,
                sweep_angle: PI,
            },
            ProfileCurve::Line {
                start: DVec2::new(0.0, 1.0),
                end: DVec2::new(2.0, 1.0),
            },
            ProfileCurve::Arc {
                center: DVec2::new(0.0, 0.5),
                radius: 0.5,
                start_angle: -PI / 2.0,
                sweep_angle: -PI,
            },
        ],
    };
    let b = k()
        .extrude(
            &p,
            &ExtrudeSpec {
                distance: 0.3,
                symmetric: false,
            },
        )
        .unwrap();
    let a = 2.0 + PI * 0.25;
    check(
        &b,
        Expect {
            volume: a * 0.3,
            area: 2.0 * a + (4.0 + PI) * 0.3,
            counts: Some((8, 12, 6)),
            genus: 0,
            solids: 1,
            rel: 1e-9,
        },
    );
}

fn polygon_area(p: &[DVec2]) -> f64 {
    let n = p.len();
    (0..n).map(|i| p[i].perp_dot(p[(i + 1) % n])).sum::<f64>() * 0.5
}

/// Dense polyline of a profile spline, for reference integrals.
fn spline_samples(c: &ProfileCurve, n: usize) -> Vec<DVec2> {
    let ProfileCurve::Spline {
        degree,
        knots,
        points,
        weights,
    } = c
    else {
        panic!("not a spline")
    };
    let curve = openshape_nurbs::NurbsCurve::new(
        *degree,
        knots.clone(),
        points.iter().map(|p| p.extend(0.0)).collect(),
        weights.clone(),
    )
    .unwrap();
    let (a, b) = curve.domain();
    (0..=n)
        .map(|i| curve.point(a + (b - a) * i as f64 / n as f64).truncate())
        .collect()
}

#[test]
fn spline_profile_extrude() {
    let pts = [
        DVec2::new(0.0, 0.0),
        DVec2::new(1.0, 0.4),
        DVec2::new(2.0, -0.2),
        DVec2::new(3.0, 0.5),
        DVec2::new(4.0, 0.0),
    ];
    let spline = ProfileCurve::interpolate(&pts).unwrap();
    let p = Profile {
        plane: Frame::WORLD,
        curves: vec![
            spline.clone(),
            ProfileCurve::Line {
                start: DVec2::new(4.0, 0.0),
                end: DVec2::new(4.0, -2.0),
            },
            ProfileCurve::Line {
                start: DVec2::new(4.0, -2.0),
                end: DVec2::new(0.0, -2.0),
            },
            ProfileCurve::Line {
                start: DVec2::new(0.0, -2.0),
                end: DVec2::new(0.0, 0.0),
            },
        ],
    };
    let mut poly = spline_samples(&spline, 20000);
    poly.extend([DVec2::new(4.0, -2.0), DVec2::new(0.0, -2.0)]);
    let area = polygon_area(&poly).abs();
    let spline_len: f64 = spline_samples(&spline, 20000)
        .windows(2)
        .map(|w| (w[1] - w[0]).length())
        .sum();
    let h = 0.7;
    let b = k()
        .extrude(
            &p,
            &ExtrudeSpec {
                distance: h,
                symmetric: false,
            },
        )
        .unwrap();
    check(
        &b,
        Expect {
            volume: area * h,
            area: 2.0 * area + (spline_len + 2.0 + 4.0 + 2.0) * h,
            counts: Some((8, 12, 6)),
            genus: 0,
            solids: 1,
            rel: 1e-6,
        },
    );
    assert!(k().face_kinds(&b).contains(&SurfaceKind::Extrusion));
}

fn revolve_z(p: Profile, angle: f64) -> openshape_kernel::Body {
    k().revolve(
        &p,
        &RevolveSpec {
            axis_origin: DVec2::ZERO,
            axis_dir: DVec2::Y,
            angle,
        },
    )
    .unwrap()
}

#[test]
fn revolve_rectangle_on_axis_is_a_cylinder() {
    let (r, h) = (0.3, 0.8);
    let b = revolve_z(
        Profile::rectangle(Frame::front(), DVec2::ZERO, DVec2::new(r, h)),
        TAU,
    );
    check(
        &b,
        Expect {
            volume: PI * r * r * h,
            area: 2.0 * PI * r * r + 2.0 * PI * r * h,
            counts: Some((6, 10, 6)),
            genus: 0,
            solids: 1,
            rel: 1e-9,
        },
    );
    // Seams: the side's split lines and the caps' diameters; creases: the
    // two rims.
    let (seam, crease) = seam_and_crease_length(&b);
    assert_close(seam, 2.0 * h + 4.0 * r, 1e-12, "seams");
    assert_close(crease, 2.0 * TAU * r, 1e-4, "rims");
    let kinds = k().face_kinds(&b);
    assert_eq!(
        kinds
            .iter()
            .filter(|k| **k == SurfaceKind::Cylinder)
            .count(),
        2
    );
    assert_eq!(
        kinds.iter().filter(|k| **k == SurfaceKind::Plane).count(),
        4
    );
}

#[test]
fn revolve_circle_is_a_torus() {
    let (big, small) = (2.0, 0.5);
    let b = revolve_z(
        Profile::circle(Frame::front(), DVec2::new(big, 0.3), small),
        TAU,
    );
    check(
        &b,
        Expect {
            volume: 2.0 * PI * PI * big * small * small,
            area: 4.0 * PI * PI * big * small,
            counts: Some((4, 8, 4)),
            genus: 1,
            solids: 1,
            rel: 1e-9,
        },
    );
    assert!(k().face_kinds(&b).iter().all(|k| *k == SurfaceKind::Torus));
    let mesh = k().tessellate(&b, &TessellationOptions::default()).unwrap();
    assert!(mesh.seams.iter().all(|s| *s), "a torus has no creases");
}

#[test]
fn revolve_half_disc_is_a_sphere() {
    let r = 0.25;
    let p = Profile {
        plane: Frame::front(),
        curves: vec![
            ProfileCurve::Arc {
                center: DVec2::ZERO,
                radius: r,
                start_angle: -PI / 2.0,
                sweep_angle: PI,
            },
            ProfileCurve::Line {
                start: DVec2::new(0.0, r),
                end: DVec2::new(0.0, -r),
            },
        ],
    };
    let b = revolve_z(p, TAU);
    check(
        &b,
        Expect {
            volume: 4.0 / 3.0 * PI * r * r * r,
            area: 4.0 * PI * r * r,
            counts: Some((2, 2, 2)),
            genus: 0,
            solids: 1,
            rel: 1e-9,
        },
    );
    assert!(k().face_kinds(&b).iter().all(|k| *k == SurfaceKind::Sphere));
    let mesh = k().tessellate(&b, &TessellationOptions::default()).unwrap();
    assert!(mesh.seams.iter().all(|s| *s), "a sphere has no creases");
}

#[test]
fn revolve_triangle_is_a_cone() {
    let (r, h) = (0.4, 1.0);
    let p = Profile::polygon(
        Frame::front(),
        &[DVec2::ZERO, DVec2::new(r, 0.0), DVec2::new(0.0, h)],
    );
    let b = revolve_z(p, -TAU);
    let slant = (r * r + h * h).sqrt();
    check(
        &b,
        Expect {
            volume: PI * r * r * h / 3.0,
            area: PI * r * r + PI * r * slant,
            counts: Some((4, 6, 4)),
            genus: 0,
            solids: 1,
            rel: 1e-9,
        },
    );
    let kinds = k().face_kinds(&b);
    assert_eq!(kinds.iter().filter(|k| **k == SurfaceKind::Cone).count(), 2);
}

#[test]
fn partial_revolves() {
    // Off-axis rectangle, 90°: one piece, two caps.
    let (r1, r2, h) = (1.0, 1.5, 0.4);
    let rect = Profile::rectangle(Frame::front(), DVec2::new(r1, 0.0), DVec2::new(r2, h));
    let q = revolve_z(rect.clone(), PI / 2.0);
    let ann = 0.5 * (r2 * r2 - r1 * r1);
    check(
        &q,
        Expect {
            volume: PI / 2.0 * ann * h,
            area: 2.0 * (PI / 2.0 * ann) + PI / 2.0 * (r1 + r2) * h + 2.0 * (r2 - r1) * h,
            counts: Some((8, 12, 6)),
            genus: 0,
            solids: 1,
            rel: 1e-9,
        },
    );
    // 270° of a rectangle touching the axis: three pieces, the axis edge is
    // shared by both caps.
    let theta = 1.5 * PI;
    let (r, h) = (0.5, 0.2);
    let b = revolve_z(
        Profile::rectangle(Frame::front(), DVec2::ZERO, DVec2::new(r, h)),
        theta,
    );
    check(
        &b,
        Expect {
            volume: theta / 2.0 * r * r * h,
            area: theta * r * r + theta * r * h + 2.0 * r * h,
            counts: None,
            genus: 0,
            solids: 1,
            rel: 1e-9,
        },
    );
}

#[test]
fn revolve_spline_vase_matches_pappus() {
    let pts = [
        DVec2::new(0.3, 0.0),
        DVec2::new(0.5, 0.3),
        DVec2::new(0.25, 0.7),
        DVec2::new(0.35, 1.0),
    ];
    let spline = ProfileCurve::interpolate(&pts).unwrap();
    let p = Profile {
        plane: Frame::front(),
        curves: vec![
            spline.clone(),
            ProfileCurve::Line {
                start: DVec2::new(0.35, 1.0),
                end: DVec2::new(0.0, 1.0),
            },
            ProfileCurve::Line {
                start: DVec2::new(0.0, 1.0),
                end: DVec2::new(0.0, 0.0),
            },
            ProfileCurve::Line {
                start: DVec2::new(0.0, 0.0),
                end: DVec2::new(0.3, 0.0),
            },
        ],
    };
    // Pappus: V = 2π ∫∫ x dA over the profile region.
    let mut poly = spline_samples(&spline, 20000);
    poly.extend([DVec2::new(0.0, 1.0), DVec2::new(0.0, 0.0)]);
    let n = poly.len();
    let (mut a, mut mx) = (0.0, 0.0);
    for i in 0..n {
        let (p0, p1) = (poly[i], poly[(i + 1) % n]);
        let c = p0.perp_dot(p1);
        a += c * 0.5;
        mx += (p0.x + p1.x) * c / 6.0;
    }
    let volume = (2.0 * PI * mx).abs();
    let b = revolve_z(p, TAU);
    let kern = k();
    assert!(kern.validate(&b).is_empty());
    let m = kern.mass_properties(&b, 1e-4).unwrap();
    assert!(
        ((m.volume - volume) / volume).abs() < 1e-6,
        "{} vs {}",
        m.volume,
        volume
    );
    assert!(a.abs() > 0.0);
    assert_watertight(
        &kern
            .tessellate(&b, &TessellationOptions::default())
            .unwrap(),
    );
    assert!(kern.face_kinds(&b).contains(&SurfaceKind::Revolution));
}

#[test]
fn invalid_inputs_fail_cleanly() {
    let kern = k();
    let sq = Profile::rectangle(Frame::front(), DVec2::new(-1.0, 0.0), DVec2::new(1.0, 1.0));
    // Crosses the axis.
    assert!(kern
        .revolve(
            &sq,
            &RevolveSpec {
                axis_origin: DVec2::ZERO,
                axis_dir: DVec2::Y,
                angle: 1.0
            }
        )
        .is_err());
    // Bad angles and distances.
    let ok = Profile::rectangle(Frame::front(), DVec2::new(1.0, 0.0), DVec2::new(2.0, 1.0));
    for angle in [0.0, 7.0, f64::NAN] {
        assert!(kern
            .revolve(
                &ok,
                &RevolveSpec {
                    axis_origin: DVec2::ZERO,
                    axis_dir: DVec2::Y,
                    angle
                }
            )
            .is_err());
    }
    assert!(kern
        .extrude(
            &ok,
            &ExtrudeSpec {
                distance: 0.0,
                symmetric: false
            }
        )
        .is_err());
    // Circle tangent to the axis at a non-vertex point.
    let tangent = Profile::circle(Frame::front(), DVec2::new(1.0, 0.0), 1.0);
    assert!(kern
        .revolve(
            &tangent,
            &RevolveSpec {
                axis_origin: DVec2::ZERO,
                axis_dir: DVec2::Y,
                angle: TAU
            }
        )
        .is_err());
}

fn convex_polygon(n: usize, radii: &[f64], phase: f64) -> Vec<DVec2> {
    (0..n)
        .map(|i| {
            let a = phase + TAU * i as f64 / n as f64;
            DVec2::new(a.cos(), a.sin()) * radii[i % radii.len()]
        })
        .collect()
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(48))]

    /// Random star-shaped prisms: valid, watertight, V/E/F = 2n/3n/n+2 and
    /// volume = area × height exactly (planar faces only).
    #[test]
    fn random_prisms(n in 3usize..12, radii in prop::collection::vec(0.2f64..2.0, 12), phase in 0.0f64..TAU, h in -3.0f64..3.0, sym in any::<bool>()) {
        prop_assume!(h.abs() > 1e-3);
        let pts = convex_polygon(n, &radii, phase);
        let area = polygon_area(&pts).abs();
        let b = k().extrude(&Profile::polygon(Frame::WORLD, &pts), &ExtrudeSpec { distance: h, symmetric: sym }).unwrap();
        prop_assert!(k().validate(&b).is_empty());
        let t = k().topology(&b);
        prop_assert_eq!((t.vertices, t.edges, t.faces), (2 * n, 3 * n, n + 2));
        let m = k().mass_properties(&b, 1e-3).unwrap();
        prop_assert!((m.volume - area * h.abs()).abs() < 1e-9 * (1.0 + area * h.abs()));
        let mesh = k().tessellate(&b, &TessellationOptions::default()).unwrap();
        assert_watertight(&mesh);
        // Every edge of a convex prism is a crease.
        prop_assert!(mesh.seams.iter().all(|s| !s));
    }

    /// Random off-axis rectangles revolved by random angles match the
    /// analytic volume θ/2 (r2² − r1²) h.
    #[test]
    fn random_revolved_rings(r1 in 0.0f64..1.0, w in 0.05f64..1.0, h in 0.05f64..1.0, theta in 0.1f64..TAU) {
        let r2 = r1 + w;
        let b = revolve_z(Profile::rectangle(Frame::front(), DVec2::new(r1, 0.0), DVec2::new(r2, h)), theta);
        prop_assert!(k().validate(&b).is_empty());
        let m = k().mass_properties(&b, 1e-4).unwrap();
        let expect = theta / 2.0 * (r2 * r2 - r1 * r1) * h;
        prop_assert!(((m.volume - expect) / expect).abs() < 1e-9, "{} vs {}", m.volume, expect);
        assert_watertight(&k().tessellate(&b, &TessellationOptions::default()).unwrap());
        // The creases are the four rims and the two end rectangles (whose
        // edges on the axis is shared when r1 = 0); everything else is a
        // seam where a surface is split.
        let (_, crease) = seam_and_crease_length(&b);
        let ends = if r1 > 0.0 { 4.0 * (w + h) } else { 4.0 * w + 3.0 * h };
        let expect = theta * 2.0 * (r1 + r2) + ends;
        prop_assert!(((crease - expect) / expect).abs() < 1e-4, "creases {} vs {}", crease, expect);
    }
}

#[test]
fn cylinder_inertia_and_centroid() {
    let (r, h) = (0.2, 0.9);
    let b = k()
        .extrude(
            &Profile::circle(Frame::WORLD, DVec2::new(0.3, 0.1), r),
            &ExtrudeSpec {
                distance: h,
                symmetric: false,
            },
        )
        .unwrap();
    let m = k().mass_properties(&b, 1e-3).unwrap();
    let mass = PI * r * r * h;
    assert!((m.centroid - DVec3::new(0.3, 0.1, h / 2.0)).length() < 1e-12);
    let izz = mass * r * r / 2.0;
    let ixx = mass * (3.0 * r * r + h * h) / 12.0;
    assert!(
        (m.inertia[2][2] - izz).abs() < 1e-12 * mass,
        "{} vs {}",
        m.inertia[2][2],
        izz
    );
    assert!((m.inertia[0][0] - ixx).abs() < 1e-12 * mass);
    assert!((m.inertia[1][1] - ixx).abs() < 1e-12 * mass);
    assert!(m.inertia[0][1].abs() < 1e-12 * mass);
}

/// Regression: revolving about an axis away from the plane origin produced
/// parameter values of ~1e-60 that the triangulator rejected.
#[test]
fn revolve_about_offset_axis() {
    let x0 = -0.06;
    let (r, h) = (0.012, 0.05);
    let p = |x: f64, y: f64| DVec2::new(x0 + x, y);
    let spline = ProfileCurve::interpolate(&[
        p(r, 0.0),
        p(1.6 * r, 0.3 * h),
        p(0.8 * r, 0.7 * h),
        p(1.1 * r, h),
    ])
    .unwrap();
    let profile = Profile {
        plane: Frame::front(),
        curves: vec![
            spline,
            ProfileCurve::Line {
                start: p(1.1 * r, h),
                end: p(0.0, h),
            },
            ProfileCurve::Line {
                start: p(0.0, h),
                end: p(0.0, 0.0),
            },
            ProfileCurve::Line {
                start: p(0.0, 0.0),
                end: p(r, 0.0),
            },
        ],
    };
    let b = k()
        .revolve(
            &profile,
            &RevolveSpec {
                axis_origin: DVec2::new(x0, 0.0),
                axis_dir: DVec2::Y,
                angle: TAU,
            },
        )
        .unwrap();
    assert!(k().validate(&b).is_empty());
    assert_watertight(&k().tessellate(&b, &TessellationOptions::default()).unwrap());
    let m = k().mass_properties(&b, 1e-3).unwrap();
    assert!((m.centroid.x - x0).abs() < 1e-12 && m.centroid.y.abs() < 1e-12);
}
