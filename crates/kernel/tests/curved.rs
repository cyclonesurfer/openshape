//! Booleans of solids bounded by planes and cylinders (design note 0012
//! §10.7, item 2): golden cases with exact volumes, areas and topology.

use std::f64::consts::{FRAC_PI_2, FRAC_PI_4, PI, TAU};

use opendatum_kernel::{
    Body, BodyMesh, BooleanOp, BrepKernel, DVec2, DVec3, EdgeOrigin, ExtrudeSpec, FaceIx,
    FaceOrigin, Frame, Kernel, KernelError, Merged, Operand, Profile, ProfileCurve, RevolveSpec,
    TessellationOptions, VertexOrigin,
};

fn k() -> BrepKernel {
    BrepKernel::default()
}

fn cuboid(lo: DVec3, hi: DVec3) -> Body {
    k().make_box(
        &Frame {
            origin: lo,
            ..Frame::WORLD
        },
        hi - lo,
    )
    .expect("a box")
}

/// A cylinder of radius `r` on the plane `frame` about `centre`, `h` long
/// along the frame's normal.
fn cylinder(frame: Frame, centre: DVec2, r: f64, h: f64) -> Body {
    k().extrude(
        &Profile::circle(frame, centre, r),
        &ExtrudeSpec {
            distance: h,
            symmetric: false,
        },
    )
    .expect("a cylinder")
}

fn at_z(z: f64) -> Frame {
    Frame {
        origin: DVec3::new(0.0, 0.0, z),
        ..Frame::WORLD
    }
}

fn volume(b: &Body) -> f64 {
    k().mass_properties(b, 1e-12)
        .expect("mass properties")
        .volume
}

fn area(b: &Body) -> f64 {
    k().mass_properties(b, 1e-12).expect("mass properties").area
}

fn counts(b: &Body) -> (usize, usize, usize) {
    let t = k().topology(b);
    (t.faces, t.edges, t.vertices)
}

fn run(a: &Body, b: &Body, op: BooleanOp) -> Body {
    let r = k()
        .boolean(a, b, op)
        .unwrap_or_else(|e| panic!("{op:?}: {e}"));
    assert_eq!(k().validate(&r), Vec::<String>::new());
    origins_agree(a, b, &r, op).unwrap_or_else(|e| panic!("{op:?}: {e}"));
    r
}

/// Checks a result's origins against the operands' surfaces: a face on
/// the surface of each operand face it comes from, facing the same way
/// (the other way for the second operand of a difference); an edge on the
/// surfaces of the faces along an operand edge it is a piece of, or of the
/// two faces where the bodies meet; a vertex at an operand vertex. Within
/// three tolerances, the most a cluster moves a point.
fn origins_agree(a: &Body, b: &Body, r: &Body, op: BooleanOp) -> Result<(), String> {
    let mesh = |x: &Body| {
        k().tessellate(x, &mesh_options())
            .map_err(|e| e.to_string())
    };
    origins_agree_on((a, &mesh(a)?), (b, &mesh(b)?), (r, &mesh(r)?), op)
}

/// How far the meshes of the checks are from the surfaces: points a few
/// hundredths of the unit apart at most (the default mesh of a unit ball
/// has 52,720, which made the checks take seconds).
const CHORD: f64 = 1e-3;

/// The options of the checks' meshes.
fn mesh_options() -> TessellationOptions {
    TessellationOptions {
        chord_tolerance: CHORD,
        max_angle: 0.5,
    }
}

/// [`origins_agree`] on the bodies' meshes.
fn origins_agree_on(
    (a, ma): (&Body, &BodyMesh),
    (b, mb): (&Body, &BodyMesh),
    (r, mr): (&Body, &BodyMesh),
    op: BooleanOp,
) -> Result<(), String> {
    let near = 3.0 * k().tolerance().linear;
    let (ja, jb) = (k().adjacency(a), k().adjacency(b));
    let body = |o: Operand| match o {
        Operand::First => (a, ma, &ja),
        Operand::Second => (b, mb, &jb),
    };
    // How far `p` is from the surface of face `f` of `x`, and the outward
    // normal there; `None` where the surface has no tangent plane (a
    // cone's apex), which a face or an edge reaches at its corners only.
    let on = |x: &Body, f: FaceIx, p: DVec3| -> Option<(f64, DVec3)> {
        let t = k().tangent_plane(x, f, p)?;
        Some(((p - t.point).length(), t.normal))
    };
    // Mesh points of a face or an edge at an apex: at most its corners.
    const SINGULAR: usize = 2;
    let o = k().origins(r);
    let face = |i: usize, origin: &FaceOrigin| -> Result<(), String> {
        let FaceOrigin::Operand { operand, face } = *origin else {
            return Err(format!("face {i}: origin {origin:?}"));
        };
        let (x, _, _) = body(operand);
        let flip = op == BooleanOp::Subtract && operand == Operand::Second;
        let mut singular = 0;
        for p in &mr.faces[i].positions {
            let (Some((d, n)), Some((_, m))) = (on(x, face, *p), on(r, FaceIx(i as u32), *p))
            else {
                singular += 1;
                continue;
            };
            if d > near {
                return Err(format!(
                    "face {i} leaves the surface of {operand:?} face {} by {d:e}",
                    face.0
                ));
            }
            let facing = m.dot(n);
            if (facing < 0.0) != flip || facing.abs() < 1.0 - 1e-6 {
                return Err(format!(
                    "face {i} does not face as {operand:?} face {} does ({facing})",
                    face.0
                ));
            }
        }
        if singular > SINGULAR {
            return Err(format!("face {i}: no tangent plane at {singular} points"));
        }
        Ok(())
    };
    let edge = |i: usize, origin: &EdgeOrigin| -> Result<(), String> {
        let faces: Vec<(&Body, FaceIx)> = match *origin {
            EdgeOrigin::Operand { operand, edge } => {
                let (x, _, adj) = body(operand);
                adj.edge_faces[edge.0 as usize]
                    .iter()
                    .map(|&f| (x, f))
                    .collect()
            }
            EdgeOrigin::Meet { first, second } => vec![(a, first), (b, second)],
            _ => return Err(format!("edge {i}: origin {origin:?}")),
        };
        let mut singular = 0;
        for p in &mr.edges[i] {
            let mut apex = false;
            for &(x, f) in &faces {
                let Some((d, _)) = on(x, f, *p) else {
                    apex = true;
                    continue;
                };
                if d > near {
                    return Err(format!(
                        "edge {i} ({origin:?}) leaves the surface of face {} by {d:e}",
                        f.0
                    ));
                }
            }
            singular += usize::from(apex);
        }
        if singular > SINGULAR {
            return Err(format!("edge {i}: no tangent plane at {singular} points"));
        }
        Ok(())
    };
    let vertex = |i: usize, origin: &VertexOrigin| -> Result<(), String> {
        match *origin {
            VertexOrigin::Operand { operand, vertex } => {
                let (_, m, _) = body(operand);
                if (mr.vertices[i] - m.vertices[vertex.0 as usize]).length() > near {
                    return Err(format!("vertex {i} is off {operand:?} vertex {}", vertex.0));
                }
                Ok(())
            }
            VertexOrigin::Meet => Ok(()),
            _ => Err(format!("vertex {i}: origin {origin:?}")),
        }
    };
    for (i, f) in o.faces.iter().enumerate() {
        face(i, f)?;
    }
    for (i, e) in o.edges.iter().enumerate() {
        edge(i, e)?;
    }
    for (i, v) in o.vertices.iter().enumerate() {
        vertex(i, v)?;
    }
    for m in &o.merged {
        match m {
            Merged::Face(f, origin) => face(f.0 as usize, origin)?,
            Merged::Edge(e, origin) => edge(e.0 as usize, origin)?,
            Merged::Vertex(v, origin) => vertex(v.0 as usize, origin)?,
        }
    }
    Ok(())
}

fn close(x: f64, y: f64) -> bool {
    (x - y).abs() <= 1e-9 * (1.0 + y.abs())
}

#[test]
fn a_plate_with_a_round_hole() {
    let plate = cuboid(DVec3::new(-2.0, -2.0, 0.0), DVec3::new(2.0, 2.0, 1.0));
    let rod = cylinder(at_z(-1.0), DVec2::ZERO, 0.5, 3.0);
    let r = run(&plate, &rod, BooleanOp::Subtract);
    let hole = PI * 0.25;
    assert!(close(volume(&r), 16.0 - hole), "{}", volume(&r));
    let want = 2.0 * (16.0 - hole) + 4.0 * 4.0 + 2.0 * PI * 0.5;
    assert!(close(area(&r), want), "{} against {want}", area(&r));
    // The box's six faces and the hole's two halves (ADR 0008); its
    // twelve edges, two arcs at each end of the hole and two seams.
    assert_eq!(counts(&r), (8, 18, 12));
}

#[test]
fn a_boss_through_and_on_a_plate() {
    let plate = cuboid(DVec3::new(-2.0, -2.0, 0.0), DVec3::new(2.0, 2.0, 1.0));
    for (bottom, label) in [(0.5, "through the top"), (1.0, "standing on the top")] {
        let rod = cylinder(at_z(bottom), DVec2::ZERO, 0.5, 2.0 - bottom);
        let r = run(&plate, &rod, BooleanOp::Union);
        let disc = PI * 0.25;
        assert!(close(volume(&r), 16.0 + disc), "{label}: {}", volume(&r));
        let want = 48.0 + PI;
        assert!(
            close(area(&r), want),
            "{label}: {} against {want}",
            area(&r)
        );
        assert_eq!(counts(&r), (9, 18, 12), "{label}");
    }
}

#[test]
fn a_d_cut_shaft() {
    let shaft = cylinder(Frame::WORLD, DVec2::ZERO, 1.0, 3.0);
    let cut = cuboid(DVec3::new(0.6, -2.0, -1.0), DVec3::new(2.0, 2.0, 4.0));
    let r = run(&shaft, &cut, BooleanOp::Subtract);
    let segment = 0.6f64.acos() - 0.6 * 0.8;
    let section = PI - segment;
    assert!(close(volume(&r), 3.0 * section), "{}", volume(&r));
    let want = 2.0 * section + 3.0 * (2.0 * PI - 2.0 * 0.6f64.acos()) + 3.0 * 1.6;
    assert!(close(area(&r), want), "{} against {want}", area(&r));
    // Both halves of the side cut back to the flat; the ends a chord and
    // two arcs.
    assert_eq!(counts(&r), (5, 9, 6));
}

#[test]
fn two_parallel_cylinders() {
    let a = cylinder(Frame::WORLD, DVec2::ZERO, 1.0, 2.0);
    let b = cylinder(Frame::WORLD, DVec2::new(1.0, 0.0), 1.0, 2.0);
    let lens = 2.0 * 0.5f64.acos() - 0.5 * 3.0f64.sqrt();
    let arc = 2.0 * PI - 2.0 * 0.5f64.acos();
    for (op, section, rim) in [
        (BooleanOp::Union, 2.0 * PI - lens, 2.0 * arc),
        (BooleanOp::Intersect, lens, 2.0 * (2.0 * PI - arc)),
        (BooleanOp::Subtract, PI - lens, arc + (2.0 * PI - arc)),
    ] {
        let r = run(&a, &b, op);
        assert!(close(volume(&r), 2.0 * section), "{op:?}: {}", volume(&r));
        let want = 2.0 * section + 2.0 * rim;
        assert!(close(area(&r), want), "{op:?}: {} against {want}", area(&r));
        // Four pieces of the sides and the two merged ends.
        assert_eq!(counts(&r), (6, 12, 8), "{op:?}");
    }
}

#[test]
fn a_counterbore() {
    let plate = cuboid(DVec3::new(-2.0, -2.0, 0.0), DVec3::new(2.0, 2.0, 1.0));
    let hole = cylinder(at_z(-1.0), DVec2::ZERO, 0.5, 3.0);
    let bore = cylinder(at_z(0.6), DVec2::ZERO, 0.8, 1.4);
    let r = run(
        &run(&plate, &hole, BooleanOp::Subtract),
        &bore,
        BooleanOp::Subtract,
    );
    let want = 16.0 - PI * 0.25 - PI * (0.64 - 0.25) * 0.4;
    assert!(close(volume(&r), want), "{} against {want}", volume(&r));
    let want = 48.0 + 0.74 * PI;
    assert!(close(area(&r), want), "{} against {want}", area(&r));
    // The box, the hole's and the bore's halves and the bore's floor.
    assert_eq!(counts(&r), (11, 24, 16));
}

#[test]
fn voids_and_pieces_with_curved_faces() {
    // A box cut from inside a cylinder: a void of flat faces.
    let rod = cylinder(Frame::WORLD, DVec2::ZERO, 1.0, 3.0);
    let cube = cuboid(DVec3::new(-0.3, -0.3, 1.0), DVec3::new(0.3, 0.3, 2.0));
    let r = run(&rod, &cube, BooleanOp::Subtract);
    assert!(close(volume(&r), 3.0 * PI - 0.36), "{}", volume(&r));
    let want = 8.0 * PI + 3.12;
    assert!(close(area(&r), want), "{} against {want}", area(&r));
    let t = k().topology(&r);
    assert_eq!((t.solids, t.shells), (1, 2));
    assert_eq!(counts(&r), (10, 18, 12));

    // A cylinder cut from inside a box: a void with curved faces.
    let block = cuboid(DVec3::new(-2.0, -2.0, 0.0), DVec3::new(2.0, 2.0, 3.0));
    let pin = cylinder(at_z(1.0), DVec2::ZERO, 0.5, 1.0);
    let r = run(&block, &pin, BooleanOp::Subtract);
    assert!(close(volume(&r), 48.0 - 0.25 * PI), "{}", volume(&r));
    let want = 80.0 + 1.5 * PI;
    assert!(close(area(&r), want), "{} against {want}", area(&r));
    let t = k().topology(&r);
    assert_eq!((t.solids, t.shells), (1, 2));
    assert_eq!(counts(&r), (10, 18, 12));

    // A slab cut across a cylinder: two solids.
    let slab = cuboid(DVec3::new(-2.0, -2.0, 1.0), DVec3::new(2.0, 2.0, 2.0));
    let r = run(&rod, &slab, BooleanOp::Subtract);
    assert!(close(volume(&r), 2.0 * PI), "{}", volume(&r));
    let want = 8.0 * PI;
    assert!(close(area(&r), want), "{} against {want}", area(&r));
    let t = k().topology(&r);
    assert_eq!((t.solids, t.shells), (2, 2));
    assert_eq!(counts(&r), (8, 12, 8));
}

// ---------------------------------------------------------------------
// Random pairs: the volume identities of design note 0012 §7 on solids
// bounded by planes and cylinders.

/// SplitMix64: a small deterministic generator.
struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    fn unit(&mut self) -> f64 {
        (self.next() >> 11) as f64 / (1u64 << 53) as f64
    }

    fn range(&mut self, lo: f64, hi: f64) -> f64 {
        lo + (hi - lo) * self.unit()
    }

    fn chance(&mut self, p: f64) -> bool {
        self.unit() < p
    }

    /// A value on a grid of quarters (to make coincidences) or anywhere.
    fn coord(&mut self, lo: f64, hi: f64, snap: bool) -> f64 {
        let x = self.range(lo, hi);
        if snap {
            (x * 4.0).round() / 4.0
        } else {
            x
        }
    }
}

/// A frame at `origin` turned at random.
fn turned(r: &mut Rng, origin: DVec3) -> Frame {
    let n = DVec3::new(r.range(-1.0, 1.0), r.range(-1.0, 1.0), r.range(-1.0, 1.0));
    let n = if n.length() < 0.1 { DVec3::Z } else { n };
    let hint = DVec3::new(r.range(-1.0, 1.0), r.range(-1.0, 1.0), r.range(-1.0, 1.0));
    Frame::from_normal(origin, n, hint)
}

fn line(a: DVec2, b: DVec2) -> ProfileCurve {
    ProfileCurve::Line { start: a, end: b }
}

fn arc(center: DVec2, radius: f64, start_angle: f64, sweep_angle: f64) -> ProfileCurve {
    ProfileCurve::Arc {
        center,
        radius,
        start_angle,
        sweep_angle,
    }
}

/// A rectangle `w` by `h` about `c` with corners rounded to `rad` (a slot
/// where `rad` is half the lesser side).
fn rounded(c: DVec2, w: f64, h: f64, rad: f64) -> Vec<ProfileCurve> {
    let (x, y) = (0.5 * w - rad, 0.5 * h - rad);
    let mut out = Vec::new();
    let corners = [(x, -y), (x, y), (-x, y), (-x, -y)];
    for (k, (cx, cy)) in corners.into_iter().enumerate() {
        let start = (k as f64 - 1.0) * FRAC_PI_2;
        let centre = c + DVec2::new(cx, cy);
        let (from, _) = (centre + DVec2::from_angle(start) * rad, ());
        // The side before this corner.
        let (px, py) = corners[(k + 3) % 4];
        let prev = c + DVec2::new(px, py) + DVec2::from_angle(start) * rad;
        if (from - prev).length() > 1e-12 {
            out.push(line(prev, from));
        }
        out.push(arc(centre, rad, start, FRAC_PI_2));
    }
    out
}

/// A disc of radius `rad` about `c` cut off by a chord at `d` from its
/// centre, toward `angle`.
fn d_shape(c: DVec2, rad: f64, d: f64, angle: f64) -> Vec<ProfileCurve> {
    let half = (d / rad).acos();
    let a = DVec2::from_angle(angle + half) * rad;
    let b = DVec2::from_angle(angle - half) * rad;
    vec![
        arc(c, rad, angle + half, TAU - 2.0 * half),
        line(c + b, c + a),
    ]
}

/// A solid turned about the line through `origin` along `axis`'s normal:
/// the polygon `pts` of (radius, height) points revolved by `angle`. A full
/// turn makes faces of half a turn each, its caps planes recorded as
/// surfaces of revolution.
fn turned_solid(axis: Frame, origin: DVec3, pts: &[(f64, f64)], angle: f64) -> Body {
    let plane = Frame {
        origin,
        x: axis.x,
        y: axis.z,
        z: -axis.y,
    };
    let pts: Vec<DVec2> = pts.iter().map(|&(rho, h)| DVec2::new(rho, h)).collect();
    k().revolve(
        &Profile::polygon(plane, &pts),
        &RevolveSpec {
            axis_origin: DVec2::ZERO,
            axis_dir: DVec2::Y,
            angle,
        },
    )
    .expect("a revolve")
}

/// A random solid of the kinds step 2 handles: a solid turned about
/// `axis`'s normal (a rod, a tube or a stepped rod, turned a whole turn or
/// half of one), a box turned any way, or a prism along `axis`'s normal of
/// a profile with arcs (a circle, a slot or rounded rectangle, a D, a
/// polygon with a round hole, a disc with a polygonal hole).
fn random_curved(r: &mut Rng, axis: Frame) -> Body {
    let snap = r.chance(0.5);
    let centre = DVec3::new(
        r.coord(-0.5, 0.5, snap),
        r.coord(-0.5, 0.5, snap),
        r.coord(-0.5, 0.5, snap),
    );
    if r.chance(0.2) {
        let rad = r.coord(0.25, 1.0, snap).max(0.25);
        let (h0, h1) = (
            -r.coord(0.25, 0.75, snap).max(0.25),
            r.coord(0.25, 0.75, snap).max(0.25),
        );
        let pts = match r.next() % 3 {
            0 => vec![(0.0, h0), (rad, h0), (rad, h1), (0.0, h1)],
            1 => {
                let inner = rad * if snap { 0.5 } else { r.range(0.3, 0.7) };
                vec![(inner, h0), (rad, h0), (rad, h1), (inner, h1)]
            }
            _ => {
                let small = rad * if snap { 0.5 } else { r.range(0.3, 0.7) };
                let hm = 0.5 * (h0 + h1);
                vec![
                    (0.0, h0),
                    (rad, h0),
                    (rad, hm),
                    (small, hm),
                    (small, h1),
                    (0.0, h1),
                ]
            }
        };
        let angle = if r.chance(0.8) { TAU } else { PI };
        let origin = axis.origin + axis.x * centre.x + axis.y * centre.y + axis.z * centre.z;
        return turned_solid(axis, origin, &pts, angle);
    }
    if r.chance(0.3) {
        let size = DVec3::new(
            r.coord(0.25, 1.5, snap).max(0.25),
            r.coord(0.25, 1.5, snap).max(0.25),
            r.coord(0.25, 1.5, snap).max(0.25),
        );
        let frame = if r.chance(0.5) {
            Frame {
                origin: centre - size * 0.5,
                ..axis
            }
        } else {
            turned(r, centre)
        };
        return k().make_box(&frame, size).expect("a box");
    }
    // The profile's plane through `centre`, square to the pair's axis.
    let plane = Frame {
        origin: axis.origin + centre,
        ..axis
    };
    let c = DVec2::new(r.coord(-0.5, 0.5, snap), r.coord(-0.5, 0.5, snap));
    let rad = r.coord(0.25, 1.0, snap).max(0.25);
    let curves = match r.next() % 5 {
        0 => vec![ProfileCurve::Circle {
            center: c,
            radius: rad,
        }],
        1 => {
            let w = rad * 2.0 + r.coord(0.0, 1.0, snap);
            let h = rad * 2.0
                + if r.chance(0.5) {
                    0.0
                } else {
                    r.coord(0.25, 1.0, snap)
                };
            rounded(c, w, h, rad)
        }
        2 => {
            let d = rad * r.coord(0.25, 0.75, snap);
            let angle = if snap {
                (r.next() % 4) as f64 * FRAC_PI_2
            } else {
                r.range(0.0, TAU)
            };
            d_shape(c, rad, d, angle)
        }
        3 => {
            let n = 3 + (r.next() % 4) as usize;
            let big = rad + 0.6;
            let mut pts: Vec<DVec2> = (0..n)
                .map(|i| c + DVec2::from_angle(i as f64 * TAU / n as f64 + 0.2) * big * 1.6)
                .collect();
            pts.dedup();
            let mut out: Vec<ProfileCurve> =
                (0..n).map(|i| line(pts[i], pts[(i + 1) % n])).collect();
            out.push(ProfileCurve::Circle {
                center: c,
                radius: rad * 0.5,
            });
            out
        }
        _ => {
            let m = 3 + (r.next() % 3) as usize;
            let mut out = vec![ProfileCurve::Circle {
                center: c,
                radius: rad + 0.25,
            }];
            let pts: Vec<DVec2> = (0..m)
                .map(|i| c + DVec2::from_angle(i as f64 * TAU / m as f64 + 0.3) * rad * 0.5)
                .collect();
            out.extend((0..m).map(|i| line(pts[i], pts[(i + 1) % m])));
            out
        }
    };
    let distance = if snap {
        r.coord(0.25, 1.5, true).max(0.25)
    } else {
        r.range(0.3, 1.5)
    };
    k().extrude(
        &Profile { plane, curves },
        &ExtrudeSpec {
            distance,
            symmetric: r.chance(0.3),
        },
    )
    .expect("an extrude")
}

/// A random pair: prisms along one axis (the world's or a turned one), so
/// that their cylinders are parallel or coaxial, and boxes at any angle.
fn random_pair(r: &mut Rng) -> (Body, Body) {
    let axis = if r.chance(0.5) {
        Frame::WORLD
    } else {
        turned(r, DVec3::ZERO)
    };
    (random_curved(r, axis), random_curved(r, axis))
}

/// A random pair of prisms along different axes: their cylinders meet at
/// an angle (a quartic curve, refused) or not at all, and their planes
/// meet the other's cylinders in ellipses.
fn random_pair_at_an_angle(r: &mut Rng) -> (Body, Body) {
    let first = if r.chance(0.5) {
        Frame::WORLD
    } else {
        turned(r, DVec3::ZERO)
    };
    let second = turned(r, DVec3::ZERO);
    (random_curved(r, first), random_curved(r, second))
}

#[derive(Default)]
struct Tally {
    pairs: usize,
    checked: usize,
    overlapping: usize,
    refused: Vec<String>,
    bugs: Vec<String>,
}

/// A boolean's outcome: `None` for a refusal, `Some(None)` for an empty
/// result, else the result's volume and mesh.
type Outcome = Option<Option<(f64, Mesh)>>;

/// A boolean's outcome; refusals and kernel bugs are recorded.
fn measured(
    t: &mut Tally,
    (a, ma): (&Body, &Mesh),
    (b, mb): (&Body, &Mesh),
    op: BooleanOp,
    case: &str,
) -> Outcome {
    match k().boolean(a, b, op) {
        Ok(r) => {
            let issues = k().validate(&r);
            if !issues.is_empty() {
                t.bugs
                    .push(format!("{case} {op:?}: invalid result: {issues:?}"));
                return None;
            }
            let mr = match Mesh::of(&r) {
                Ok(m) => m,
                Err(e) => {
                    t.bugs
                        .push(format!("{case} {op:?}: meshing the result: {e}"));
                    return None;
                }
            };
            if let Err(e) = origins_agree_on((a, &ma.raw), (b, &mb.raw), (&r, &mr.raw), op) {
                t.bugs.push(format!("{case} {op:?}: origins: {e}"));
                return None;
            }
            match k().mass_properties(&r, 1e-12) {
                Ok(m) => Some(Some((m.volume, mr))),
                Err(e) => {
                    t.bugs.push(format!("{case} {op:?}: {e}"));
                    None
                }
            }
        }
        Err(KernelError::InvalidInput(m)) if m.contains("empty") => Some(None),
        Err(KernelError::InvalidInput(m) | KernelError::Unsupported(m)) => {
            t.refused.push(format!("{case} {op:?}: {m}"));
            None
        }
        Err(e) => {
            t.bugs.push(format!("{case} {op:?}: {e}"));
            None
        }
    }
}

/// A body's surface as triangles within [`CHORD`] of it, for deciding
/// whether points well away from it are inside: a check of the booleans
/// that shares none of their code (the volume identities hold for a result
/// that misses where the bodies meet altogether, as if they were apart).
/// The triangles are in a tree of boxes (a bounding volume hierarchy: T.
/// L. Kay and J. T. Kajiya, "Ray tracing complex scenes", SIGGRAPH 1986),
/// so a point walks down the tree instead of visiting every triangle.
struct Mesh {
    raw: BodyMesh,
    tris: Vec<[DVec3; 3]>,
    nodes: Vec<Node>,
}

/// A node of a mesh's tree: the box of the triangles `start..end`, and
/// its two children unless it is a leaf.
struct Node {
    lo: DVec3,
    hi: DVec3,
    start: usize,
    end: usize,
    children: Option<(usize, usize)>,
}

/// Directions of rays for [`Mesh::inside`], away from the axes and the
/// diagonals so that they rarely run along a face of a test body.
const RAYS: [DVec3; 3] = [
    DVec3::new(0.544_2, 0.327_1, 0.772_3),
    DVec3::new(-0.615_7, 0.741_1, 0.268_3),
    DVec3::new(0.291_7, -0.439_2, 0.849_7),
];

impl Mesh {
    fn of(b: &Body) -> Result<Mesh, String> {
        let raw = k()
            .tessellate(b, &mesh_options())
            .map_err(|e| e.to_string())?;
        let (pos, tris) = raw.welded();
        let mut tris: Vec<[DVec3; 3]> = tris.iter().map(|t| t.map(|i| pos[i as usize])).collect();
        let mut nodes = Vec::new();
        Self::build(&mut tris, 0, &mut nodes);
        Ok(Mesh { raw, tris, nodes })
    }

    /// Adds the node of `tris`, which start at `start` in the mesh, and
    /// below it nodes of their halves at the median of their centres along
    /// the box's longest side, down to four triangles; returns its index.
    fn build(tris: &mut [[DVec3; 3]], start: usize, nodes: &mut Vec<Node>) -> usize {
        let (lo, hi) = tris.iter().flatten().fold(
            (DVec3::splat(f64::MAX), DVec3::splat(f64::MIN)),
            |(lo, hi), q| (lo.min(*q), hi.max(*q)),
        );
        let i = nodes.len();
        nodes.push(Node {
            lo,
            hi,
            start,
            end: start + tris.len(),
            children: None,
        });
        if tris.len() > 4 {
            let size = hi - lo;
            let axis = if size.x >= size.y.max(size.z) {
                0
            } else if size.y >= size.z {
                1
            } else {
                2
            };
            let key = |t: &[DVec3; 3]| (t[0] + t[1] + t[2])[axis];
            let mid = tris.len() / 2;
            tris.select_nth_unstable_by(mid, |s, t| key(s).total_cmp(&key(t)));
            let (first, second) = tris.split_at_mut(mid);
            let l = Self::build(first, start, nodes);
            let r = Self::build(second, start + mid, nodes);
            nodes[i].children = Some((l, r));
        }
        i
    }

    /// Whether `p` is inside, if it is farther than `margin` from the
    /// mesh (so the mesh and the exact surface agree there): by the parity
    /// of the triangles a ray from `p` crosses, or where every ray passes
    /// too close to an edge or runs along a triangle, by the winding number.
    fn inside(&self, p: DVec3, margin: f64) -> Option<bool> {
        if self.near(p, margin + CHORD) {
            return None;
        }
        Some(
            match RAYS.iter().find_map(|d| self.crossings(p, d.normalize())) {
                Some(n) => n % 2 == 1,
                None => self.winding(p) > 0.5,
            },
        )
    }

    /// Whether a triangle comes within `r` of `p`.
    fn near(&self, p: DVec3, r: f64) -> bool {
        let mut stack = vec![0];
        while let Some(i) = stack.pop() {
            let n = &self.nodes[i];
            if (n.lo - p).max(p - n.hi).max(DVec3::ZERO).length() > r {
                continue;
            }
            match n.children {
                Some((a, b)) => stack.extend([a, b]),
                None => {
                    if self.tris[n.start..n.end]
                        .iter()
                        .any(|t| triangle_distance(p, t) <= r)
                    {
                        return true;
                    }
                }
            }
        }
        false
    }

    /// How many triangles the ray from `p` along the unit vector `d`
    /// crosses (T. Möller and B. Trumbore, "Fast, minimum storage
    /// ray-triangle intersection", J. Graphics Tools 2(1), 1997), or `None`
    /// if it runs within a millionth of a radian of a triangle's plane, or
    /// meets a triangle so near its boundary that rounding could put the
    /// crossing on either side: there it could be counted twice or not at
    /// all.
    fn crossings(&self, p: DVec3, d: DVec3) -> Option<usize> {
        let inv = d.recip();
        let reaches = |lo: DVec3, hi: DVec3| {
            let (t0, t1) = ((lo - p) * inv, (hi - p) * inv);
            let far = t0.max(t1).min_element();
            t0.min(t1).max_element() <= far && far >= 0.0
        };
        let mut count = 0;
        let mut stack = vec![0];
        while let Some(i) = stack.pop() {
            let n = &self.nodes[i];
            if !reaches(n.lo, n.hi) {
                continue;
            }
            if let Some((a, b)) = n.children {
                stack.extend([a, b]);
                continue;
            }
            for t in &self.tris[n.start..n.end] {
                let (e1, e2) = (t[1] - t[0], t[2] - t[0]);
                let h = d.cross(e2);
                let det = e1.dot(h);
                if det.abs() <= 1e-6 * e1.cross(e2).length() {
                    if reaches(t[0].min(t[1]).min(t[2]), t[0].max(t[1]).max(t[2])) {
                        return None;
                    }
                    continue;
                }
                let s = p - t[0];
                let q = s.cross(e1);
                if e2.dot(q) / det <= 0.0 {
                    continue;
                }
                let (u, v) = (s.dot(h) / det, d.dot(q) / det);
                let least = u.min(v).min(1.0 - u - v);
                // A generous bound on the rounding of the three
                // coordinates: of `s`, of the products over the sizes of
                // their factors, and of the division by `det`.
                let size = e1.length() + e2.length();
                let reach = s.length() + p.length() + t[0].length() + (u.abs() + v.abs()) * size;
                if least.abs() <= 64.0 * f64::EPSILON * reach * size / det.abs() {
                    return None;
                }
                count += usize::from(least > 0.0);
            }
        }
        Some(count)
    }

    /// The winding number of the mesh about `p`: the sum of the
    /// triangles' solid angles about it over 4π (A. van Oosterom and J.
    /// Strackee, "The solid angle of a plane triangle", IEEE Trans.
    /// Biomed. Eng. 30, 1983).
    fn winding(&self, p: DVec3) -> f64 {
        let mut total = 0.0;
        for t in &self.tris {
            let [a, b, c] = t.map(|q| q - p);
            let (la, lb, lc) = (a.length(), b.length(), c.length());
            let num = a.dot(b.cross(c));
            let den = la * lb * lc + a.dot(b) * lc + b.dot(c) * la + c.dot(a) * lb;
            total += 2.0 * num.atan2(den);
        }
        total / (4.0 * PI)
    }
}

/// The distance from `p` to the triangle `t`.
fn triangle_distance(p: DVec3, t: &[DVec3; 3]) -> f64 {
    let [a, b, c] = *t;
    let n = (b - a).cross(c - a);
    if n.length_squared() > 0.0 {
        // Inside the triangle's prism: the distance to its plane.
        let side = |u: DVec3, v: DVec3| (v - u).cross(p - u).dot(n);
        if side(a, b) >= 0.0 && side(b, c) >= 0.0 && side(c, a) >= 0.0 {
            return (p - a).dot(n).abs() / n.length();
        }
    }
    let seg = |u: DVec3, v: DVec3| {
        let d = v - u;
        let s = if d.length_squared() > 0.0 {
            ((p - u).dot(d) / d.length_squared()).clamp(0.0, 1.0)
        } else {
            0.0
        };
        (p - (u + d * s)).length()
    };
    seg(a, b).min(seg(b, c)).min(seg(c, a))
}

/// Points sampled where the bodies' boxes meet (or anywhere in their
/// boxes if they do not) must be inside each result exactly when the
/// operation's rule on the operands says so.
fn membership(
    t: &mut Tally,
    ma: &Mesh,
    mb: &Mesh,
    results: &[(BooleanOp, bool, &Outcome)],
    case: &str,
    rng: &mut Rng,
) {
    const MARGIN: f64 = 1e-3;
    let ((la, ha), (lb, hb)) = (
        (ma.nodes[0].lo, ma.nodes[0].hi),
        (mb.nodes[0].lo, mb.nodes[0].hi),
    );
    let (lo, hi) = if la.cmple(hb).all() && lb.cmple(ha).all() {
        (la.max(lb), ha.min(hb))
    } else {
        (la.min(lb), ha.max(hb))
    };
    let meshes: Vec<(BooleanOp, bool, Option<&Mesh>)> = results
        .iter()
        .filter_map(|&(op, swapped, r)| Some((op, swapped, r.as_ref()?.as_ref().map(|x| &x.1))))
        .collect();
    for _ in 0..24 {
        let p = DVec3::new(
            rng.range(lo.x, hi.x),
            rng.range(lo.y, hi.y),
            rng.range(lo.z, hi.z),
        );
        let (Some(ia), Some(ib)) = (ma.inside(p, MARGIN), mb.inside(p, MARGIN)) else {
            continue;
        };
        for (op, swapped, m) in &meshes {
            let (x, y) = if *swapped { (ib, ia) } else { (ia, ib) };
            let want = match op {
                BooleanOp::Union => x || y,
                BooleanOp::Intersect => x && y,
                BooleanOp::Subtract => x && !y,
            };
            let got = match m {
                None => Some(false),
                Some(m) => m.inside(p, MARGIN),
            };
            if got.is_some_and(|g| g != want) {
                t.bugs.push(format!(
                    "{case} {op:?}{}: the result is {} at {p:?}, which the operands put {}",
                    if *swapped { " (B − A)" } else { "" },
                    if want { "outside" } else { "inside" },
                    if want { "inside" } else { "outside" },
                ));
                return;
            }
        }
    }
}

fn identities(t: &mut Tally, a: &Body, b: &Body, case: &str, slack: f64, rng: &mut Rng) {
    t.pairs += 1;
    let (va, vb) = (volume(a), volume(b));
    let (ma, mb) = match (Mesh::of(a), Mesh::of(b)) {
        (Ok(x), Ok(y)) => (x, y),
        (Err(e), _) | (_, Err(e)) => {
            t.bugs.push(format!("{case}: meshing an operand: {e}"));
            return;
        }
    };
    let (a, b) = ((a, &ma), (b, &mb));
    let ru = measured(t, a, b, BooleanOp::Union, case);
    let ri = measured(t, a, b, BooleanOp::Intersect, case);
    let rab = measured(t, a, b, BooleanOp::Subtract, case);
    let rba = measured(t, b, a, BooleanOp::Subtract, case);
    let vol = |r: &Outcome| r.as_ref().map(|x| x.as_ref().map_or(0.0, |y| y.0));
    let (u, i, ab, ba) = (vol(&ru), vol(&ri), vol(&rab), vol(&rba));
    let tol = 1e-9 * va.max(vb) + slack;
    let mut check = |name: &str, lhs: f64, rhs: f64| {
        if (lhs - rhs).abs() > tol {
            t.bugs.push(format!(
                "{case}: {name}: {lhs} vs {rhs} (|A| {va}, |B| {vb})"
            ));
        }
    };
    if let Some(i) = i {
        if i > 0.0 {
            t.overlapping += 1;
        }
        if let Some(u) = u {
            check("|A∪B| + |A∩B| = |A| + |B|", u + i, va + vb);
        }
        if let Some(ab) = ab {
            check("|A−B| + |A∩B| = |A|", ab + i, va);
        }
        if let Some(ba) = ba {
            check("|B−A| + |A∩B| = |B|", ba + i, vb);
        }
        if u.is_some() && ab.is_some() && ba.is_some() {
            t.checked += 1;
        }
    }
    membership(
        t,
        &ma,
        &mb,
        &[
            (BooleanOp::Union, false, &ru),
            (BooleanOp::Intersect, false, &ri),
            (BooleanOp::Subtract, false, &rab),
            (BooleanOp::Subtract, true, &rba),
        ],
        case,
        rng,
    );
}

fn report(name: &str, t: &Tally) {
    println!(
        "{name}: {} pairs ({} overlapping), {} with every identity checked; {} refusals in {} \
         operations ({:.1} per thousand)",
        t.pairs,
        t.overlapping,
        t.checked,
        t.refused.len(),
        4 * t.pairs,
        1000.0 * t.refused.len() as f64 / (4 * t.pairs).max(1) as f64
    );
    let limit = if std::env::var_os("BOOLEAN_REFUSALS").is_some() {
        usize::MAX
    } else {
        12
    };
    for m in t.refused.iter().take(limit) {
        println!("  refused: {m}");
    }
    assert!(
        t.bugs.is_empty(),
        "{} failures:\n{}",
        t.bugs.len(),
        t.bugs.join("\n")
    );
}

/// A test's case count, times `BOOLEAN_SCALE` when set (for long runs).
fn scale(n: usize) -> usize {
    std::env::var("BOOLEAN_SCALE")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .map_or(n, |k| n * k)
}

#[test]
fn volume_identities_on_random_curved_pairs() {
    let mut t = Tally::default();
    let mut r = Rng(0x0DA7_0C11);
    let mut points = Rng(0x5A3F_0E11);
    for case in 0..scale(200) {
        let (a, b) = random_pair(&mut r);
        identities(&mut t, &a, &b, &format!("pair {case}"), 0.0, &mut points);
    }
    report("random solids of planes and cylinders", &t);
}

#[test]
fn volume_identities_on_random_pairs_at_an_angle() {
    let mut t = Tally::default();
    let mut r = Rng(0x0A76_1E55);
    let mut points = Rng(0x5A3F_0E12);
    for case in 0..scale(100) {
        let (a, b) = random_pair_at_an_angle(&mut r);
        identities(&mut t, &a, &b, &format!("pair {case}"), 0.0, &mut points);
    }
    report("random solids of planes and cylinders at an angle", &t);
}

/// A rod dipping into a pipe's wall from outside: the faces of the rod
/// have their boundaries wholly outside the pipe, yet meet it, in quartic
/// curves (two cylinders at an angle). Refused, never taken as apart.
#[test]
fn a_rod_dipping_into_a_pipe_wall_is_refused() {
    let pipe = cylinder(at_z(-5.0), DVec2::ZERO, 10.0, 10.0);
    for (y, z) in [(10.5, 0.0), (10.5, 3.0), (-10.5, -3.0)] {
        let f = Frame::from_normal(DVec3::new(-4.0, y, z), DVec3::X, DVec3::Z);
        let rod = cylinder(f, DVec2::ZERO, 1.0, 24.0);
        for op in [BooleanOp::Union, BooleanOp::Intersect, BooleanOp::Subtract] {
            match k().boolean(&pipe, &rod, op) {
                Err(KernelError::Unsupported(m)) => {
                    assert!(
                        m.contains("quartic curve (design note 0012 §10.3)"),
                        "{op:?}: {m}"
                    )
                }
                other => panic!("({y}, {z}) {op:?}: {:?}", other.map(|b| volume(&b))),
            }
        }
    }
}

/// Turned rods and tubes (their caps planes recorded as surfaces of
/// revolution, in two halves each) cut by boxes: a cap kept whole is one
/// face again, recorded as a plane, and the volumes are exact.
/// A vertex on the other body's surface whose face touches the surface
/// there from one side, while another edge at the vertex leaves to the
/// other side: the bodies cross, and the curve where they meet passes
/// through the vertex. The rule that refuses faces touching at a vertex
/// refused three of these as tangent contacts until it asked the vertex's
/// other edges and faces (found by diffing the fuzz's refusals): a box's
/// corner on a cylinder, two edges going in and one out; a cylinder
/// crossing a turned rod along its seam, where it touches one half of
/// the rod's cap at the seam's end and enters the other; and a cylinder
/// through the centre of a rod's cap, along the line that halves it. The
/// same corner on a box's face was done throughout (the face that
/// touches is apart from it). Each is done, with the volume identities
/// and the membership of sampled points.
#[test]
fn corners_and_seams_crossing_the_other_surface() {
    // Two of the box's edges at its corner go in, one goes out.
    let corner = |at: DVec3| {
        k().make_box(
            &Frame {
                origin: at,
                x: DVec3::new(-2.0, 1.0, -2.0) / 3.0,
                y: DVec3::new(1.0, -2.0, -2.0) / 3.0,
                z: DVec3::new(-2.0, -2.0, 1.0) / 3.0,
            },
            DVec3::splat(0.5),
        )
        .expect("a box")
    };
    let pairs = [
        (
            corner(DVec3::new(0.25, 0.0, 0.0)),
            cylinder(at_z(-1.0), DVec2::ZERO, 0.25, 2.0),
        ),
        (
            corner(DVec3::new(0.0, 0.0, 0.0)),
            cuboid(DVec3::splat(-1.0), DVec3::new(1.0, 1.0, 0.0)),
        ),
        (
            turned_solid(
                Frame::WORLD,
                DVec3::new(0.25, -0.5, 0.0),
                &[(0.0, -0.375), (0.125, -0.375), (0.125, 0.0), (0.0, 0.0)],
                TAU,
            ),
            cylinder(at_z(-1.0), DVec2::new(0.375, 0.0), 0.5, 2.0),
        ),
        (
            turned_solid(
                Frame::WORLD,
                DVec3::ZERO,
                &[(0.0, -0.5), (0.5, -0.5), (0.5, 0.0), (0.0, 0.0)],
                TAU,
            ),
            cylinder(at_z(-1.0), DVec2::new(0.0, 0.2), 0.2, 2.0),
        ),
    ];
    let mut t = Tally::default();
    let mut rng = Rng(0x5A3F_0E19);
    for (i, (a, b)) in pairs.iter().enumerate() {
        identities(&mut t, a, b, &format!("pair {i}"), 0.0, &mut rng);
    }
    assert!(t.refused.is_empty(), "{:#?}", t.refused);
    assert!(t.bugs.is_empty(), "{:#?}", t.bugs);
    assert_eq!(t.checked, pairs.len());
}

#[test]
fn turned_rods_and_tubes_cut_by_boxes() {
    let rod = turned_solid(
        Frame::WORLD,
        DVec3::ZERO,
        &[(0.0, 0.0), (1.0, 0.0), (1.0, 2.0), (0.0, 2.0)],
        TAU,
    );
    let tube = turned_solid(
        Frame::WORLD,
        DVec3::ZERO,
        &[(0.5, 0.0), (1.0, 0.0), (1.0, 2.0), (0.5, 2.0)],
        TAU,
    );
    let top = cuboid(DVec3::new(-2.0, -2.0, 1.5), DVec3::new(2.0, 2.0, 3.0));
    let side = cuboid(DVec3::new(0.5, -2.0, -1.0), DVec3::new(2.0, 2.0, 3.0));
    // The part of the unit disc beyond x = 0.5, twice over (z from 0 to 2).
    let segment = 2.0 * (PI / 3.0 - 3f64.sqrt() / 4.0);
    let cases = [
        (&rod, &top, [24.0 + 1.5 * PI, 0.5 * PI, 1.5 * PI], [9, 4, 4]),
        (
            &rod,
            &side,
            [24.0 + 2.0 * PI - segment, segment, 2.0 * PI - segment],
            [10, 5, 5],
        ),
        (
            &tube,
            &top,
            [24.0 + 1.125 * PI, 0.375 * PI, 1.125 * PI],
            [12, 6, 6],
        ),
    ];
    for (a, b, volumes, faces) in cases {
        for (k, op) in [BooleanOp::Union, BooleanOp::Intersect, BooleanOp::Subtract]
            .into_iter()
            .enumerate()
        {
            let r = run(a, b, op);
            assert!(
                close(volume(&r), volumes[k]),
                "{op:?}: {} vs {}",
                volume(&r),
                volumes[k]
            );
            assert_eq!(counts(&r).0, faces[k], "{op:?}");
        }
    }
}

/// A half-turn rod (a D-shaped half cylinder) cut down its side by a box:
/// its two flat sides lie on one plane and join into one face, and each
/// cap's two radii, both then between the cap and that face, used to join
/// into one diameter through the axis, where the cap's surface of
/// revolution is singular. The mass properties of such a result did not
/// converge (or, rarely, converged to a wrong volume); the centre stays a
/// vertex now, and the validator refuses an edge through a face's axis.
#[test]
fn a_half_turned_rod_keeps_the_centres_of_its_caps() {
    let rod = turned_solid(
        Frame::WORLD,
        DVec3::ZERO,
        &[(0.0, 0.0), (1.0, 0.0), (1.0, 2.0), (0.0, 2.0)],
        PI,
    );
    let side = cuboid(DVec3::new(0.5, -2.0, -1.0), DVec3::new(2.0, 2.0, 3.0));
    // Half of the part of the unit disc beyond x = 0.5, twice over.
    let piece = PI / 3.0 - 3f64.sqrt() / 4.0;
    let volumes = [24.0 + PI - piece, piece, PI - piece];
    // Faces, edges and vertices: each cap keeps its centre.
    let topology = [(10, 23, 16), (5, 9, 6), (5, 11, 8)];
    for (k, op) in [BooleanOp::Union, BooleanOp::Intersect, BooleanOp::Subtract]
        .into_iter()
        .enumerate()
    {
        let r = run(&rod, &side, op);
        assert!(
            close(volume(&r), volumes[k]),
            "{op:?}: {} vs {}",
            volume(&r),
            volumes[k]
        );
        assert_eq!(counts(&r), topology[k], "{op:?}");
    }
}

/// A prism of `curves` on `plane`, `distance` long.
fn prism(plane: Frame, curves: Vec<ProfileCurve>, distance: f64, symmetric: bool) -> Body {
    k().extrude(
        &Profile { plane, curves },
        &ExtrudeSpec {
            distance,
            symmetric,
        },
    )
    .expect("a prism")
}

/// Pairs a fuzz with turned solids found faulty (2026-10-13), with their
/// operands as drawn: every operation must succeed, valid, with the volume
/// identities and the membership of sampled points.
///
/// 1. A disc with a square hole and a stepped rod turned half a turn on
///    it: a cap's centre came out of its chart at 8e-49, which the
///    triangulation of a piece refuses as too small.
/// 2. A slot and a stepped rod at an angle, and 3. a slot and a box at an
///    angle: a curve's end moved 2.7e-8 to its cluster, more than a
///    sixteenth of the tolerance, and its end chord was halved without end
///    to narrow a band that halving cannot narrow.
#[test]
fn pairs_a_fuzz_of_turned_solids_found() {
    let pairs = [
        (
            prism(
                Frame {
                    origin: DVec3::new(0.0, 0.0, 0.0),
                    x: DVec3::new(1.0, 0.0, 0.0),
                    y: DVec3::new(0.0, 1.0, 0.0),
                    z: DVec3::new(0.0, 0.0, 1.0),
                },
                vec![
                    ProfileCurve::Circle {
                        center: DVec2::new(-0.5, -0.25),
                        radius: 1.0,
                    },
                    ProfileCurve::Line {
                        start: DVec2::new(-0.14174881657789773, -0.13917992250199768),
                        end: DVec2::new(-0.6108200774980024, 0.10825118342210227),
                    },
                    ProfileCurve::Line {
                        start: DVec2::new(-0.6108200774980024, 0.10825118342210227),
                        end: DVec2::new(-0.8582511834221023, -0.3608200774980022),
                    },
                    ProfileCurve::Line {
                        start: DVec2::new(-0.8582511834221023, -0.3608200774980022),
                        end: DVec2::new(-0.3891799225019978, -0.6082511834221023),
                    },
                    ProfileCurve::Line {
                        start: DVec2::new(-0.3891799225019978, -0.6082511834221023),
                        end: DVec2::new(-0.14174881657789773, -0.13917992250199768),
                    },
                ],
                0.5,
                false,
            ),
            turned_solid(
                Frame {
                    origin: DVec3::new(0.0, 0.0, 0.0),
                    x: DVec3::new(1.0, 0.0, 0.0),
                    y: DVec3::new(0.0, 1.0, 0.0),
                    z: DVec3::new(0.0, 0.0, 1.0),
                },
                DVec3::new(-0.25, 0.25, 0.5),
                &[
                    (0.0, -0.25),
                    (0.75, -0.25),
                    (0.75, 0.25),
                    (0.375, 0.25),
                    (0.375, 0.75),
                    (0.0, 0.75),
                ],
                PI,
            ),
        ),
        (
            prism(
                Frame {
                    origin: DVec3::new(
                        0.3198798658815176,
                        0.009509306982644183,
                        0.14267829601661042,
                    ),
                    x: DVec3::new(
                        -0.5095606575811561,
                        -0.30556614931921056,
                        0.8043489694378242,
                    ),
                    y: DVec3::new(0.42423748938859496, 0.7240690809675032, 0.5438258164009292),
                    z: DVec3::new(
                        -0.7485789796960115,
                        0.6183472280014727,
                        -0.23932450100266184,
                    ),
                },
                vec![
                    ProfileCurve::Line {
                        start: DVec2::new(-0.6128598662078211, -0.8688728997063812),
                        end: DVec2::new(-0.07551125986357149, -0.8688728997063812),
                    },
                    ProfileCurve::Arc {
                        center: DVec2::new(-0.07551125986357154, 0.06831019919853742),
                        radius: 0.9371830989049186,
                        start_angle: -FRAC_PI_2,
                        sweep_angle: FRAC_PI_2,
                    },
                    ProfileCurve::Line {
                        start: DVec2::new(0.861671839041347, 0.06831019919853742),
                        end: DVec2::new(0.861671839041347, 0.6849782878507099),
                    },
                    ProfileCurve::Arc {
                        center: DVec2::new(-0.07551125986357154, 0.6849782878507099),
                        radius: 0.9371830989049186,
                        start_angle: 0.0,
                        sweep_angle: FRAC_PI_2,
                    },
                    ProfileCurve::Line {
                        start: DVec2::new(-0.07551125986357149, 1.6221613867556286),
                        end: DVec2::new(-0.6128598662078211, 1.6221613867556286),
                    },
                    ProfileCurve::Arc {
                        center: DVec2::new(-0.6128598662078212, 0.6849782878507099),
                        radius: 0.9371830989049186,
                        start_angle: FRAC_PI_2,
                        sweep_angle: FRAC_PI_2,
                    },
                    ProfileCurve::Line {
                        start: DVec2::new(-1.5500429651127399, 0.68497828785071),
                        end: DVec2::new(-1.5500429651127399, 0.06831019919853754),
                    },
                    ProfileCurve::Arc {
                        center: DVec2::new(-0.6128598662078212, 0.06831019919853742),
                        radius: 0.9371830989049186,
                        start_angle: PI,
                        sweep_angle: FRAC_PI_2,
                    },
                ],
                0.8174169161354687,
                true,
            ),
            turned_solid(
                Frame {
                    origin: DVec3::new(0.0, 0.0, 0.0),
                    x: DVec3::new(
                        -0.5095606575811561,
                        -0.30556614931921056,
                        0.8043489694378242,
                    ),
                    y: DVec3::new(0.42423748938859496, 0.7240690809675032, 0.5438258164009292),
                    z: DVec3::new(
                        -0.7485789796960115,
                        0.6183472280014727,
                        -0.23932450100266184,
                    ),
                },
                DVec3::new(-0.3145349093192919, 0.07819526967056553, 0.1412561171087906),
                &[
                    (0.0, -0.25),
                    (0.75, -0.25),
                    (0.75, 0.125),
                    (0.375, 0.125),
                    (0.375, 0.5),
                    (0.0, 0.5),
                ],
                TAU,
            ),
        ),
        (
            prism(
                Frame {
                    origin: DVec3::new(
                        0.44989271050355284,
                        0.14026408013845004,
                        -0.4181833832883749,
                    ),
                    x: DVec3::new(1.0, 0.0, 0.0),
                    y: DVec3::new(0.0, 1.0, 0.0),
                    z: DVec3::new(0.0, 0.0, 1.0),
                },
                vec![
                    ProfileCurve::Line {
                        start: DVec2::new(-0.3980656106028941, -0.8620680373061178),
                        end: DVec2::new(0.41308136483509467, -0.8620680373061178),
                    },
                    ProfileCurve::Arc {
                        center: DVec2::new(0.4130813648350946, 0.05737362008177527),
                        radius: 0.9194416573878931,
                        start_angle: -FRAC_PI_2,
                        sweep_angle: FRAC_PI_2,
                    },
                    ProfileCurve::Arc {
                        center: DVec2::new(0.4130813648350946, 0.05737362008177527),
                        radius: 0.9194416573878931,
                        start_angle: 0.0,
                        sweep_angle: FRAC_PI_2,
                    },
                    ProfileCurve::Line {
                        start: DVec2::new(0.41308136483509467, 0.9768152774696683),
                        end: DVec2::new(-0.3980656106028941, 0.9768152774696683),
                    },
                    ProfileCurve::Arc {
                        center: DVec2::new(-0.39806561060289414, 0.05737362008177527),
                        radius: 0.9194416573878931,
                        start_angle: FRAC_PI_2,
                        sweep_angle: FRAC_PI_2,
                    },
                    ProfileCurve::Arc {
                        center: DVec2::new(-0.39806561060289414, 0.05737362008177527),
                        radius: 0.9194416573878931,
                        start_angle: PI,
                        sweep_angle: FRAC_PI_2,
                    },
                ],
                0.3403025177356964,
                false,
            ),
            k().make_box(
                &Frame {
                    origin: DVec3::new(
                        -0.06718975807541638,
                        -0.2819102258866296,
                        -0.6291372072111716,
                    ),
                    x: DVec3::new(0.09344293357510233, 0.6881420902924062, -0.7195337947121573),
                    y: DVec3::new(
                        -0.9949161505073416,
                        0.09179772557000873,
                        -0.041412933243448254,
                    ),
                    z: DVec3::new(0.03755358337804729, 0.7197455391651773, 0.6932215282485797),
                },
                DVec3::new(1.0409626070707354, 0.878961188058442, 1.1840879921925254),
            )
            .expect("a box"),
        ),
    ];
    let mut t = Tally::default();
    let mut points = Rng(0x5A3F_0E13);
    for (i, (a, b)) in pairs.iter().enumerate() {
        identities(&mut t, a, b, &format!("pair {}", i + 1), 0.0, &mut points);
    }
    assert!(t.refused.is_empty(), "{:?}", t.refused);
    assert!(t.bugs.is_empty(), "{:?}", t.bugs);
}

/// A ball of radius `r` about `centre`: a half disc turned about the
/// world's z axis, in two faces on the sphere.
fn ball(centre: DVec3, r: f64) -> Body {
    let plane = Frame {
        origin: centre,
        x: DVec3::X,
        y: DVec3::Z,
        z: -DVec3::Y,
    };
    k().revolve(
        &Profile {
            plane,
            curves: vec![
                arc(DVec2::ZERO, r, -FRAC_PI_2, PI),
                line(DVec2::new(0.0, r), DVec2::new(0.0, -r)),
            ],
        },
        &RevolveSpec {
            axis_origin: DVec2::ZERO,
            axis_dir: DVec2::Y,
            angle: TAU,
        },
    )
    .expect("a ball")
}

/// The volume of a cap of height `h` of a ball of radius `r`.
fn cap(r: f64, h: f64) -> f64 {
    PI * h * h * (3.0 * r - h) / 3.0
}

/// The volume of a cone of base radius `r` and height `h`.
fn cone_volume(r: f64, h: f64) -> f64 {
    PI * r * r * h / 3.0
}

/// Faces, edges and vertices.
type Topology = (usize, usize, usize);

/// A golden case: its name, its operands, and the volumes and topology of
/// A ∪ B, A ∩ B, A − B and B − A.
type Golden<'a> = (&'a str, &'a Body, &'a Body, [f64; 4], [Topology; 4]);

/// Balls and cones with planes, cylinders and each other (design note
/// 0012 §12.4, item 3a): every result valid with its origins on the
/// operands' surfaces, its volume exact to 1e-9, and its faces, edges and
/// vertices as first found, each count checked against Euler's formula
/// (a face with a hole adds a loop). Each case gives A ∪ B, A ∩ B, A − B
/// and B − A.
#[test]
fn balls_and_cones_with_plates_boxes_and_rods() {
    let plate = cuboid(DVec3::new(-2.0, -2.0, -1.0), DVec3::new(2.0, 2.0, 0.0));
    let ball_volume = |r: f64| 4.0 / 3.0 * PI * r.powi(3);
    // A ball resting half a radius deep in a plate.
    let resting = ball(DVec3::new(0.0, 0.0, 0.5), 1.0);
    let dip = cap(1.0, 0.5);
    // A ball cut by a box along its axis, off its centre.
    let unit = ball(DVec3::ZERO, 1.0);
    let side = cuboid(DVec3::new(0.3, -2.0, -2.0), DVec3::new(2.0, 2.0, 2.0));
    let slice = cap(1.0, 0.7);
    // A spherical pocket: a ball centred on the plate's top face.
    let pocket = ball(DVec3::ZERO, 0.8);
    let half = 0.5 * ball_volume(0.8);
    // Two balls across their axes, and two along their common axis.
    let beside = ball(DVec3::new(1.2, 0.0, 0.0), 1.0);
    let above = ball(DVec3::new(0.0, 0.0, 1.2), 1.0);
    let lens = 2.0 * cap(1.0, 0.4);
    // A countersink: a cone, apex down, 0.6 into the plate.
    let sink = turned_solid(
        Frame::WORLD,
        DVec3::new(0.0, 0.0, 0.2),
        &[(0.0, 0.0), (0.8, 0.0), (0.0, -0.8)],
        TAU,
    );
    let (sunk, sink_all) = (cone_volume(0.6, 0.6), cone_volume(0.8, 0.8));
    // A cone on a rod, coaxial, its base inside the rod.
    let rod = cylinder(at_z(0.0), DVec2::ZERO, 1.0, 1.0);
    let point = turned_solid(
        Frame::WORLD,
        DVec3::new(0.0, 0.0, 0.5),
        &[(0.0, 0.0), (0.9, 0.0), (0.0, 0.9)],
        TAU,
    );
    let (tip, point_all) = (cone_volume(0.4, 0.4), cone_volume(0.9, 0.9));
    // A pointed cone cut square to its axis half way up.
    let tall = turned_solid(
        Frame::WORLD,
        DVec3::ZERO,
        &[(0.0, 0.0), (1.0, 0.0), (0.0, 2.0)],
        TAU,
    );
    let top = cuboid(DVec3::new(-2.0, -2.0, 1.0), DVec3::new(2.0, 2.0, 3.0));
    let (tall_tip, tall_all) = (cone_volume(0.5, 1.0), cone_volume(1.0, 2.0));
    // A cone, apex down onto a ball from above, whose side crosses the
    // sphere twice close together: besides its cap below the cone's base,
    // the ball keeps a ring 0.011 thick outside the cone between the two
    // circles, which a coarse mesh of the result flattens to nothing.
    let (apex, slope, base) = (2.0, 0.57, -0.2);
    let crossing = turned_solid(
        Frame::WORLD,
        DVec3::ZERO,
        &[(0.0, base), (slope * (apex - base), base), (0.0, apex)],
        TAU,
    );
    // The sphere's radius squared less the cone's at height z,
    // (1 − z²) − m²(a − z)², is positive between the circles, its roots.
    let (qa, qb, qc) = (
        -(1.0 + slope * slope),
        2.0 * slope * slope * apex,
        1.0 - slope * slope * apex * apex,
    );
    let root = (qb * qb - 4.0 * qa * qc).sqrt();
    let (z1, z2) = ((-qb + root) / (2.0 * qa), (-qb - root) / (2.0 * qa));
    let integral = |z: f64| PI * (qa * z.powi(3) / 3.0 + 0.5 * qb * z * z + qc * z);
    let ring = (integral(z2) - integral(z1)).abs();
    let outside = cap(1.0, 1.0 + base) + ring;
    let (inside, crossing_all) = (
        ball_volume(1.0) - outside,
        cone_volume(slope * (apex - base), apex - base),
    );
    let cases: [Golden; 9] = [
        (
            "a ball resting in a plate",
            &resting,
            &plate,
            [
                ball_volume(1.0) + 16.0 - dip,
                dip,
                ball_volume(1.0) - dip,
                16.0 - dip,
            ],
            [(8, 16, 11), (3, 4, 3), (3, 4, 3), (8, 16, 11)],
        ),
        (
            "a ball cut by a box",
            &unit,
            &side,
            [
                ball_volume(1.0) + 27.2 - slice,
                slice,
                ball_volume(1.0) - slice,
                27.2 - slice,
            ],
            [(8, 17, 12), (3, 3, 2), (3, 5, 4), (8, 15, 10)],
        ),
        (
            "a spherical pocket",
            &plate,
            &pocket,
            [16.0 + half, half, 16.0 - half, half],
            [(8, 16, 11), (3, 4, 3), (8, 16, 11), (3, 4, 3)],
        ),
        (
            "two balls side by side",
            &unit,
            &beside,
            [
                2.0 * ball_volume(1.0) - lens,
                lens,
                ball_volume(1.0) - lens,
                ball_volume(1.0) - lens,
            ],
            [(4, 8, 6), (4, 4, 2), (4, 6, 4), (4, 6, 4)],
        ),
        (
            "two balls on one axis",
            &unit,
            &above,
            [
                2.0 * ball_volume(1.0) - lens,
                lens,
                ball_volume(1.0) - lens,
                ball_volume(1.0) - lens,
            ],
            [(4, 6, 4), (4, 6, 4), (4, 6, 4), (4, 6, 4)],
        ),
        (
            "a countersink",
            &plate,
            &sink,
            [16.0 + sink_all - sunk, sunk, 16.0 - sunk, sink_all - sunk],
            [(9, 18, 12), (3, 4, 3), (8, 16, 11), (4, 6, 4)],
        ),
        (
            "a cone on a rod",
            &rod,
            &point,
            [PI + tip, point_all - tip, PI - point_all + tip, tip],
            [(6, 10, 7), (4, 6, 4), (7, 12, 8), (3, 4, 3)],
        ),
        (
            "a pointed cone cut square",
            &tall,
            &top,
            [
                tall_all + 32.0 - tall_tip,
                tall_tip,
                tall_all - tall_tip,
                32.0 - tall_tip,
            ],
            [(9, 18, 12), (3, 4, 3), (4, 6, 4), (8, 16, 11)],
        ),
        (
            "a ball crossed twice by a cone",
            &unit,
            &crossing,
            [
                ball_volume(1.0) + crossing_all - inside,
                inside,
                outside,
                crossing_all - inside,
            ],
            [(9, 16, 10), (7, 12, 7), (7, 12, 7), (9, 16, 10)],
        ),
    ];
    for (name, a, b, volumes, topology) in cases {
        for (k, (x, y, op)) in [
            (a, b, BooleanOp::Union),
            (a, b, BooleanOp::Intersect),
            (a, b, BooleanOp::Subtract),
            (b, a, BooleanOp::Subtract),
        ]
        .into_iter()
        .enumerate()
        {
            let r = run(x, y, op);
            assert!(
                close(volume(&r), volumes[k]),
                "{name}, {k}: {} vs {}",
                volume(&r),
                volumes[k]
            );
            assert_eq!(counts(&r), topology[k], "{name}, {k}");
        }
    }
}

/// The area of the part of a disc of radius `r` beyond a chord `e` from
/// its centre.
fn segment(r: f64, e: f64) -> f64 {
    r * r * (e / r).acos() - e * (r * r - e * e).sqrt()
}

/// A box `across` square and `depth` deep standing on the plane through
/// `point` square to `normal`, centred on `point`, on the side `normal`
/// points to.
fn slab(point: DVec3, normal: DVec3, across: f64, depth: f64) -> Body {
    let f = Frame::from_normal(point, normal, DVec3::Y);
    k().make_box(
        &Frame {
            origin: point - (f.x + f.y) * (0.5 * across),
            ..f
        },
        DVec3::new(across, across, depth),
    )
    .expect("a box")
}

/// The volume of the part of a cone of apex (0, 0, 1) and base radius 1/2
/// on z = 0 beyond the plane through (p, 0, 0) square to `n`, a unit
/// normal in the xz plane pointing to the apex side, from first
/// principles: the plane meets the rulings in the xz plane at its
/// section's vertices, the near one inside the cone, and the base in a
/// chord across the section's axis 2w long, h from the near vertex. A
/// conic through a vertex and a chord square to its axis, with its other
/// vertex 2a away, has semi-axis b = w/√(u(2 ∓ u)), u = h/a, and its
/// segment has an area of ab(2 asin √(u/2) − (1 − u)√(u(2 − u))) for an
/// ellipse, ab((1 + u)√(u(2 + u)) − acosh(1 + u)) for a hyperbola, and
/// 4wh/3 for a parabola (Archimedes). The part is a cone over the base
/// and the section from the apex: (1/3) Σ hᵢAᵢ over its plane faces, hᵢ
/// the apex's distance below face i.
fn beyond(p: f64, n: DVec3) -> f64 {
    let (apex, at) = (DVec3::Z, DVec3::new(p, 0.0, 0.0));
    let hits = [0.5, -0.5].map(|s: f64| {
        let d = DVec3::new(s, 0.0, -1.0);
        let t = n.dot(at - apex) / n.dot(d);
        (t, apex + d * t)
    });
    let (near, far) = if (0.0..1.0).contains(&hits[0].0) {
        (hits[0], hits[1])
    } else {
        (hits[1], hits[0])
    };
    assert!((0.0..1.0).contains(&near.0), "no vertex in the cone");
    let (w, h) = ((0.25 - p * p).sqrt(), (at - near.1).length());
    let section = if !far.0.is_finite() || far.0.abs() > 1e12 {
        4.0 / 3.0 * w * h
    } else {
        let a = 0.5 * (far.1 - near.1).length();
        let u = h / a;
        if far.0 > 0.0 {
            let b = w / (u * (2.0 - u)).sqrt();
            a * b * (2.0 * (0.5 * u).sqrt().asin() - (1.0 - u) * (u * (2.0 - u)).sqrt())
        } else {
            let root = (u * (2.0 + u)).sqrt();
            let b = w / root;
            a * b * ((1.0 + u) * root - (u + root).ln_1p())
        }
    };
    let base = if n.x > 0.0 {
        0.25 * PI - segment(0.5, p)
    } else {
        segment(0.5, p)
    };
    (base + (at - apex).dot(n) * section) / 3.0
}

/// Planes across cones (design note 0012 §14, item 3b): each conic as
/// whole ellipses, half ellipses close to a parabola, a parabola and
/// hyperbolas, and two rulings where the plane passes through a cone's
/// apex outside the solid. Every result valid with its origins on the
/// operands' surfaces, its volume exact to 1e-9 against a derivation from
/// the rulings and chords (never from the kernel's own formulas), its
/// faces, edges and vertices as first found, and one solid of genus 0.
/// Each case gives A ∪ B, A ∩ B, A − B and B − A.
#[test]
fn planes_across_cones_in_ellipses_parabolas_and_hyperbolas() {
    // A pointed cone (apex (0, 0, 2), base radius 1) under a box tilted
    // 20°, through its axis half way up: a whole ellipse. Its major axis
    // runs between the rulings in the xz plane, (±t/2, 0, 2 − t), and its
    // minor axis is the chord across the circle of the cone at its centre.
    let tall = turned_solid(
        Frame::WORLD,
        DVec3::ZERO,
        &[(0.0, 0.0), (1.0, 0.0), (0.0, 2.0)],
        TAU,
    );
    let tilt = 20f64.to_radians();
    let n = DVec3::new(tilt.sin(), 0.0, tilt.cos());
    let mid = DVec3::new(0.0, 0.0, 1.0);
    let roof = slab(mid, n, 6.0, 3.0);
    let [e1, e2] = [0.5, -0.5].map(|s: f64| {
        let t = 1.0 / (1.0 - s * tilt.tan());
        DVec3::new(s * t, 0.0, 2.0 - t)
    });
    let (centre, a) = (0.5 * (e1 + e2), 0.5 * (e1 - e2).length());
    let b = ((0.5 * (2.0 - centre.z)).powi(2) - centre.x * centre.x).sqrt();
    let tip = (DVec3::new(0.0, 0.0, 2.0) - mid).dot(n) * PI * a * b / 3.0;
    let tall_all = cone_volume(1.0, 2.0);
    // A countersink turned 15° from upright, 45° each side of its axis,
    // apex down into a plate: the plate's top face meets it in a whole
    // ellipse. Its major axis runs between the rulings at 60° and −30°
    // from upright; across its centre m, a point m + y ŷ is on the cone
    // where its distance from the axis, √(|w|² + y² − s²) with w = m − apex
    // and s = w·axis, equals s.
    let plate = cuboid(DVec3::new(-2.0, -2.0, -1.0), DVec3::new(2.0, 2.0, 0.0));
    let lean = 15f64.to_radians();
    let up = DVec3::new(lean.sin(), 0.0, lean.cos());
    let sink = turned_solid(
        Frame::from_normal(DVec3::ZERO, up, DVec3::X),
        DVec3::new(0.0, 0.0, 0.4),
        &[(0.0, 0.0), (0.8, 0.0), (0.0, -0.8)],
        TAU,
    );
    let low = DVec3::new(0.0, 0.0, 0.4) - up * 0.8;
    let [x1, x2] = [lean + PI / 4.0, lean - PI / 4.0].map(|f| low.x - low.z * f.tan());
    let (m, a) = (DVec3::new(0.5 * (x1 + x2), 0.0, 0.0), 0.5 * (x1 - x2).abs());
    let (w, s) = (m - low, (m - low).dot(up));
    let b = (2.0 * s * s - w.length_squared()).sqrt();
    let sunk = -low.z * PI * a * b / 3.0;
    let sink_all = cone_volume(0.8, 0.8);
    // A cone (apex (0, 0, 1), base radius 1/2) cut by planes through the
    // base at x = 0.2: along a ruling (a parabola), and a thousandth of a
    // radian either side of one (an ellipse whose far end is 240 below,
    // and a hyperbola), each slab on the apex's side, the ellipse's also
    // a slab wide enough to hold all of it; and by a box from x = 0.1 (a
    // hyperbola, its area also ∫ 1 − 2√(0.01 + y²) dy).
    let spike = turned_solid(
        Frame::WORLD,
        DVec3::ZERO,
        &[(0.0, 0.0), (0.5, 0.0), (0.0, 1.0)],
        TAU,
    );
    let spike_all = cone_volume(0.5, 1.0);
    let ruling = 0.5f64.atan();
    let at = DVec3::new(0.2, 0.0, 0.0);
    let cut = |f: f64| DVec3::new(f.cos(), 0.0, f.sin());
    let along = slab(at, cut(ruling), 4.0, 3.0);
    let long = slab(at, cut(ruling + 1e-3), 4.0, 3.0);
    let wide = slab(at, cut(ruling + 1e-3), 1000.0, 3.0);
    let open = slab(at, cut(ruling - 1e-3), 4.0, 3.0);
    let [near_along, near_long, near_open] =
        [0.0, 1e-3, -1e-3].map(|f| spike_all - beyond(0.2, cut(ruling + f)));
    let side = cuboid(DVec3::new(0.1, -2.0, -1.0), DVec3::new(2.0, 2.0, 2.0));
    let y = 0.24f64.sqrt();
    let flank = (segment(0.5, 0.1) - 0.1 * (y - 0.02 * (y / 0.1).asinh())) / 3.0;
    assert!(close(flank, beyond(0.1, -DVec3::X)), "{flank}");
    // A frustum of 1° draft (radius 1 at z = 0, apex 1/tan 1° up) cut
    // 0.3 off its axis: a hyperbola whose vertices are far apart. Its
    // section has an area of ∫₀² 2√(ρ² − 0.09) dz, ρ = 1 − z tan 1°;
    // slicing across the axis checks the volume.
    let draft = 1f64.to_radians().tan();
    let top = 1.0 - 2.0 * draft;
    let tapered = turned_solid(
        Frame::WORLD,
        DVec3::ZERO,
        &[(0.0, 0.0), (1.0, 0.0), (top, 2.0), (0.0, 2.0)],
        TAU,
    );
    let tapered_all = 2.0 * PI * (1.0 + top + top * top) / 3.0;
    let (c, high) = (0.3, 1.0 / draft);
    let g = |u: f64| u * (u * u - c * c).sqrt() - c * c * (u / c).acosh();
    let face = (g(1.0) - g(top)) / draft;
    let shaved = (segment(1.0, c) * high + (2.0 - high) * segment(top, c) - c * face) / 3.0;
    let f = |u: f64| {
        u.powi(3) / 3.0 * (c / u).acos() - 2.0 * c / 3.0 * u * (u * u - c * c).sqrt()
            + c.powi(3) / 3.0 * (u / c).acosh()
    };
    assert!(close(shaved, (f(1.0) - f(top)) / draft), "{shaved}");
    let wall = cuboid(DVec3::new(c, -2.0, -1.0), DVec3::new(2.0, 2.0, 3.0));
    // A frustum (radius 1 at z = 0, 1/2 at z = 1) under a slab whose face
    // passes through the cone's apex at (0, 0, 2): two rulings. The part
    // beyond is a cone from the apex over a segment of the base, less one
    // over a segment of the top, a quarter its area: (2 − 1/4)/3 of it.
    let frustum = turned_solid(
        Frame::WORLD,
        DVec3::ZERO,
        &[(0.0, 0.0), (1.0, 0.0), (0.5, 1.0), (0.0, 1.0)],
        TAU,
    );
    let frustum_all = PI * 1.75 / 3.0;
    let wedge = 7.0 / 12.0 * segment(1.0, 0.4);
    let through = slab(
        DVec3::new(0.0, 0.0, 2.0),
        DVec3::new(1.0, 0.0, 0.2).normalize(),
        6.0,
        3.0,
    );
    let cases: [Golden; 9] = [
        (
            "a pointed cone cut at a slant",
            &tall,
            &roof,
            [tall_all + 108.0 - tip, tip, tall_all - tip, 108.0 - tip],
            [(9, 18, 12), (3, 4, 3), (4, 6, 4), (8, 16, 11)],
        ),
        (
            "a countersink at an angle",
            &plate,
            &sink,
            [16.0 + sink_all - sunk, sunk, 16.0 - sunk, sink_all - sunk],
            [(9, 18, 12), (3, 4, 3), (8, 16, 11), (4, 6, 4)],
        ),
        (
            "a cone cut along a ruling",
            &spike,
            &along,
            [
                spike_all + 48.0 - near_along,
                near_along,
                spike_all - near_along,
                48.0 - near_along,
            ],
            [(9, 18, 12), (4, 7, 5), (4, 6, 4), (9, 19, 13)],
        ),
        (
            "a cone cut in a long ellipse",
            &spike,
            &long,
            [
                spike_all + 48.0 - near_long,
                near_long,
                spike_all - near_long,
                48.0 - near_long,
            ],
            [(9, 18, 12), (4, 7, 5), (4, 6, 4), (9, 19, 13)],
        ),
        (
            "a cone cut in a long ellipse by a wide slab",
            &spike,
            &wide,
            [
                spike_all + 3e6 - near_long,
                near_long,
                spike_all - near_long,
                3e6 - near_long,
            ],
            [(9, 18, 12), (4, 7, 5), (4, 6, 4), (9, 19, 13)],
        ),
        (
            "a cone cut in a hyperbola close to a parabola",
            &spike,
            &open,
            [
                spike_all + 48.0 - near_open,
                near_open,
                spike_all - near_open,
                48.0 - near_open,
            ],
            [(9, 18, 12), (4, 7, 5), (4, 6, 4), (9, 19, 13)],
        ),
        (
            "a cone cut along its axis",
            &spike,
            &side,
            [
                spike_all + 22.8 - flank,
                flank,
                spike_all - flank,
                22.8 - flank,
            ],
            [(9, 19, 13), (4, 6, 4), (4, 7, 5), (9, 18, 12)],
        ),
        (
            "a slight taper cut along its axis",
            &tapered,
            &wall,
            [
                tapered_all + 6.8 * 4.0 - shaved,
                shaved,
                tapered_all - shaved,
                6.8 * 4.0 - shaved,
            ],
            [(10, 21, 14), (5, 9, 6), (5, 9, 6), (10, 21, 14)],
        ),
        (
            "a frustum cut through its apex",
            &frustum,
            &through,
            [
                frustum_all + 108.0 - wedge,
                wedge,
                frustum_all - wedge,
                108.0 - wedge,
            ],
            [(10, 21, 14), (5, 9, 6), (5, 9, 6), (10, 21, 14)],
        ),
    ];
    for (name, a, b, volumes, topology) in cases {
        for (i, (x, y, op)) in [
            (a, b, BooleanOp::Union),
            (a, b, BooleanOp::Intersect),
            (a, b, BooleanOp::Subtract),
            (b, a, BooleanOp::Subtract),
        ]
        .into_iter()
        .enumerate()
        {
            let r = run(x, y, op);
            // To 1e-9 and the rounding of a large volume (the wide slab's
            // 3e6): `close` would allow 3e-3 there, a hundredth of the cone.
            assert!(
                (volume(&r) - volumes[i]).abs() <= 1e-9 + 1e-12 * volumes[i].abs(),
                "{name}, {i}: {} vs {}",
                volume(&r),
                volumes[i]
            );
            assert_eq!(counts(&r), topology[i], "{name}, {i}");
            let t = k().topology(&r);
            assert_eq!(
                (t.solids, t.shells, t.genus),
                (1, 1, Some(0)),
                "{name}, {i}"
            );
        }
    }
}

/// A solid turned a whole turn about the line through `origin` along
/// `axis`'s normal from `curves`, (radius, height) curves in the half
/// plane of that line.
fn turned_curves(axis: Frame, origin: DVec3, curves: Vec<ProfileCurve>) -> Body {
    let plane = Frame {
        origin,
        x: axis.x,
        y: axis.z,
        z: -axis.y,
    };
    k().revolve(
        &Profile { plane, curves },
        &RevolveSpec {
            axis_origin: DVec2::ZERO,
            axis_dir: DVec2::Y,
            angle: TAU,
        },
    )
    .expect("a revolve")
}

/// A random solid turned about `axis`'s normal, its origin on that line:
/// a ball, a pointed cone (point up or down), a frustum, a rod with a
/// round end (a cylinder and a sphere meeting smoothly), or one of
/// `random_curved`'s turned rods, tubes and stepped rods. True with a
/// ball.
fn random_revolve(r: &mut Rng, axis: Frame) -> (Body, bool) {
    let snap = r.chance(0.5);
    let origin = axis.origin + axis.z * r.coord(-0.5, 0.5, snap);
    let rad = r.coord(0.25, 1.0, snap).max(0.25);
    let (h0, h1) = (
        -r.coord(0.25, 0.75, snap).max(0.25),
        r.coord(0.25, 0.75, snap).max(0.25),
    );
    let small = rad * if snap { 0.5 } else { r.range(0.3, 0.7) };
    let polygon = |pts: &[(f64, f64)]| -> Vec<ProfileCurve> {
        (0..pts.len())
            .map(|k| {
                let (a, b) = (pts[k], pts[(k + 1) % pts.len()]);
                line(DVec2::new(a.0, a.1), DVec2::new(b.0, b.1))
            })
            .collect()
    };
    let kind = r.next() % 8;
    let curves = match kind {
        0 | 1 => vec![
            arc(DVec2::ZERO, rad, -FRAC_PI_2, PI),
            line(DVec2::new(0.0, rad), DVec2::new(0.0, -rad)),
        ],
        2 => polygon(&[(0.0, h0), (rad, h0), (0.0, h1)]),
        3 => polygon(&[(0.0, h0), (rad, h1), (0.0, h1)]),
        4 => polygon(&[(0.0, h0), (rad, h0), (small, h1), (0.0, h1)]),
        5 => {
            let top = 0.5 * (h0 + h1);
            vec![
                line(DVec2::new(0.0, h0), DVec2::new(rad, h0)),
                line(DVec2::new(rad, h0), DVec2::new(rad, top)),
                arc(DVec2::new(0.0, top), rad, 0.0, FRAC_PI_2),
                line(DVec2::new(0.0, top + rad), DVec2::new(0.0, h0)),
            ]
        }
        6 => polygon(&[(small, h0), (rad, h0), (rad, h1), (small, h1)]),
        _ => polygon(&[
            (0.0, h0),
            (rad, h0),
            (rad, 0.5 * (h0 + h1)),
            (small, 0.5 * (h0 + h1)),
            (small, h1),
            (0.0, h1),
        ]),
    };
    (turned_curves(axis, origin, curves), kind < 2)
}

/// A random pair of solids turned about one axis, with faces on spheres,
/// cones, cylinders and planes square to the axis, or a ball and a box
/// turned any way or another ball anywhere: every curve where two of their
/// faces meet is a circle (design note 0012 §12.4, item 3a).
fn random_revolve_pair(r: &mut Rng) -> (Body, Body) {
    let axis = if r.chance(0.5) {
        Frame::WORLD
    } else {
        turned(r, DVec3::ZERO)
    };
    let (first, round) = random_revolve(r, axis);
    let second = if r.chance(0.5) {
        random_revolve(r, axis).0
    } else if round && r.chance(0.5) {
        if r.chance(0.5) {
            // A ball anywhere: two spheres meet in a circle about the line
            // through their centres.
            let snap = r.chance(0.5);
            let mut at = || r.coord(-1.0, 1.0, snap);
            let centre = axis.origin + DVec3::new(at(), at(), at());
            return (first, ball(centre, r.coord(0.25, 1.0, snap).max(0.25)));
        }
        let size = DVec3::new(r.range(0.25, 1.5), r.range(0.25, 1.5), r.range(0.25, 1.5));
        let corner =
            axis.origin + DVec3::new(r.range(-1.0, 0.5), r.range(-1.0, 0.5), r.range(-1.0, 0.5));
        k().make_box(&turned(r, corner), size).expect("a box")
    } else {
        // A slab square to the axis, wider than any of the solids.
        let snap = r.chance(0.5);
        let z0 = r.coord(-1.0, 0.75, snap);
        let thick = r.coord(0.25, 1.0, snap).max(0.25);
        let frame = Frame {
            origin: axis.origin - axis.x * 3.0 - axis.y * 3.0 + axis.z * z0,
            ..axis
        };
        k().make_box(&frame, DVec3::new(6.0, 6.0, thick))
            .expect("a slab")
    };
    (first, second)
}

#[test]
fn volume_identities_on_random_revolves() {
    let mut t = Tally::default();
    let mut r = Rng(0x0BA1_1C0E);
    let mut points = Rng(0x5A3F_0E14);
    for case in 0..scale(100) {
        let (a, b) = random_revolve_pair(&mut r);
        identities(&mut t, &a, &b, &format!("pair {case}"), 0.0, &mut points);
    }
    report("random balls and cones about one axis", &t);
}

/// A conical face of a turned solid: the cone through the circle of
/// radius `rho` about the solid's axis `u` at `point`, its radius growing
/// by `slope` per unit along `u`, between heights 0 and `h` from `point`.
struct Taper {
    point: DVec3,
    u: DVec3,
    rho: f64,
    slope: f64,
    h: f64,
}

/// A random solid with conical faces turned a whole turn about the line
/// through `axis`'s origin along its normal: a pointed cone either way
/// up, a frustum, a tube with a tapered bore, or a double cone; with its
/// cones.
fn random_tapered(r: &mut Rng, axis: Frame) -> (Body, Vec<Taper>) {
    let snap = r.chance(0.5);
    let rad = r.coord(0.25, 1.0, snap).max(0.25);
    let (h0, h1) = (
        -r.coord(0.25, 0.75, snap).max(0.25),
        r.coord(0.25, 0.75, snap).max(0.25),
    );
    let small = rad * if snap { 0.5 } else { r.range(0.2, 0.8) };
    let other = rad * if snap { 0.25 } else { r.range(0.1, 0.8) };
    let pts = match r.next() % 5 {
        0 => vec![(0.0, h0), (rad, h0), (0.0, h1)],
        1 => vec![(0.0, h0), (rad, h1), (0.0, h1)],
        2 => vec![(0.0, h0), (rad, h0), (small, h1), (0.0, h1)],
        3 => vec![(small, h0), (rad, h0), (rad, h1), (other, h1)],
        _ => vec![(0.0, h0), (rad, 0.5 * (h0 + h1)), (0.0, h1)],
    };
    let tapers = (0..pts.len())
        .filter_map(|k| {
            let (a, b): ((f64, f64), (f64, f64)) = (pts[k], pts[(k + 1) % pts.len()]);
            (a.0 != b.0 && a.1 != b.1).then(|| Taper {
                point: axis.origin + axis.z * a.1,
                u: axis.z,
                rho: a.0,
                slope: (b.0 - a.0) / (b.1 - a.1),
                h: b.1 - a.1,
            })
        })
        .collect();
    (turned_solid(axis, axis.origin, &pts, TAU), tapers)
}

/// A random tapered solid and a box across one of its cones (design note
/// 0012 §14, item 3b): a box turned any way (an ellipse or a hyperbola,
/// mostly); a slab whose face is parallel to the plane tangent along a
/// ruling (a parabola), or tilted from it by up to a tenth of a radian,
/// often by a thousandth to a billionth (an ellipse or a hyperbola close
/// to a parabola); or a slab whose face passes through the cone's apex or
/// within a few tolerances of it, crossing the cone (two rulings, a
/// hyperbola with its vertex by the apex, or a refusal where the apex is
/// on the solid), the apex in the slab's face or not; or a slab whose face
/// is within a ten-thousandth to a ten-millionth of a radian of tangent
/// along a ruling, through the apex or a few tolerances off it (found by
/// review: such planes made two solids a nanometre apart, and results
/// the kernel's own check refused).
fn random_conic_pair(r: &mut Rng) -> (Body, Body) {
    let axis = if r.chance(0.5) {
        Frame::WORLD
    } else {
        turned(r, DVec3::ZERO)
    };
    let (first, tapers) = random_tapered(r, axis);
    let taper = &tapers[r.next() as usize % tapers.len()];
    let u = taper.u;
    let phi = if r.chance(0.5) {
        FRAC_PI_2 * (r.next() % 4) as f64
    } else {
        r.range(0.0, TAU)
    };
    let e = axis.x * phi.cos() + axis.y * phi.sin();
    // A point of the face, the ruling there and the cone's normal.
    let s = r.range(0.2, 0.8) * taper.h;
    let rho = taper.rho + taper.slope * s;
    let on = taper.point + u * s + e * rho;
    let ruling = (u + e * taper.slope).normalize();
    let out = (e - u * taper.slope).normalize();
    let flip = if r.chance(0.5) { 1.0 } else { -1.0 };
    let second = match r.next() % 4 {
        0 => {
            let size = DVec3::new(r.range(0.25, 1.5), r.range(0.25, 1.5), r.range(0.25, 1.5));
            let corner = axis.origin
                + DVec3::new(r.range(-1.0, 0.5), r.range(-1.0, 0.5), r.range(-1.0, 0.5));
            k().make_box(&turned(r, corner), size).expect("a box")
        }
        1 => {
            let tilt = match r.next() % 5 {
                0 => 0.0,
                1 => r.range(-1e-3, 1e-3),
                2 => flip * [1e-6, 1e-9][r.next() as usize % 2],
                _ => r.range(-0.1, 0.1),
            };
            let n = out * tilt.cos() + ruling * tilt.sin();
            let point = on - out * (r.range(0.05, 0.6) * rho);
            slab(point, n * flip, r.range(1.0, 6.0), 3.0)
        }
        2 => {
            let apex = taper.point - u * (taper.rho / taper.slope);
            let psi = r.range(-0.9, 0.9) * taper.slope.abs().atan();
            let n = u.cross(e) * psi.cos() + u * psi.sin();
            let off = match r.next() % 4 {
                0 | 1 => 0.0,
                2 => r.range(-1.0, 1.0) * 6e-9,
                _ => r.range(-1.0, 1.0) * 1e-6,
            };
            let w = on - apex;
            let centre = apex + n * off + (w - n * w.dot(n)) * r.range(0.0, 1.2);
            slab(centre, n * flip, r.range(0.5, 5.0), 3.0)
        }
        _ => {
            let apex = taper.point - u * (taper.rho / taper.slope);
            let side = if r.chance(0.5) { 1.0 } else { -1.0 };
            let tilt: f64 = side * [1e-4, 1e-5, 1e-6, 1e-7][r.next() as usize % 4];
            let n = out * tilt.cos() + ruling * tilt.sin();
            let off = [0.0, 0.5, 2.0, 10.0][r.next() as usize % 4] * r.range(-1.0, 1.0).signum();
            slab(apex + n * (off * 1e-7), n * flip, r.range(1.0, 6.0), 3.0)
        }
    };
    (first, second)
}

/// A box beside a wide frustum whose axis leans, its face on y = 0 through
/// the apex: that plane meets the cone in two rulings, found from the
/// region where the faces can meet, whose centre lies below the apex's
/// height. The rulings were taken on the cone's other nappe, and every
/// operation refused (found by review).
#[test]
fn a_box_beside_a_leaning_frustum_in_a_plane_through_its_apex() {
    let (s, c) = 80f64.to_radians().sin_cos();
    let frustum = turned_solid(
        Frame::from_normal(DVec3::ZERO, DVec3::new(1.0, 0.0, 1.0).normalize(), DVec3::Y),
        DVec3::ZERO,
        &[
            (0.0, 10.0 * c),
            (10.0 * s, 10.0 * c),
            (11.0 * s, 11.0 * c),
            (0.0, 11.0 * c),
        ],
        TAU,
    );
    let block = cuboid(DVec3::new(2.0, -1.0, -7.0), DVec3::new(8.8, 0.0, -5.9));
    let mut t = Tally::default();
    identities(&mut t, &frustum, &block, "beside", 0.0, &mut Rng(0xB5DE));
    assert!(t.refused.is_empty(), "{:?}", t.refused);
    assert_eq!((t.checked, t.overlapping), (1, 1));
    report("a box beside a leaning frustum", &t);
}

/// Small tools across large cones, where the faces' boxes reach beyond
/// the region where both faces can be: a box 0.1 wide whose face cuts a
/// frustum 8 high in an ellipse whose other half (beyond the region,
/// past the ellipse's centre) the frustum's base and seams cross, and a
/// slab through a frustum's apex along one of its two rulings whose other
/// ruling the frustum's edges cross. Both were refused as meeting at too
/// shallow an angle: the tracks left out the curves' parts away from the
/// tool, which a face's boundary crossing must lie on.
#[test]
fn small_tools_across_large_cones() {
    let mut t = Tally::default();
    let mut rng = Rng(0x5A11_7001);
    // From z = 0 (radius 2) to z = 8 (radius 0.4), its apex at z = 10.
    let tall = turned_solid(
        Frame::WORLD,
        DVec3::ZERO,
        &[(0.0, 0.0), (2.0, 0.0), (0.4, 8.0), (0.0, 8.0)],
        TAU,
    );
    // Planes turned from the one tangent along the ruling at −x toward it.
    let ruling = DVec3::new(2.0, 0.0, 10.0).normalize();
    let out = DVec3::new(-10.0, 0.0, 2.0).normalize();
    for (tilt, size, side) in [(0.1, 0.1, 1.0), (0.05, 0.05, -1.0), (0.2, 0.05, -1.0)] {
        let n = out * f64::cos(tilt) + ruling * f64::sin(tilt);
        let tool = slab(DVec3::new(0.58, 0.0, 7.0), n * side, size, size);
        identities(&mut t, &tall, &tool, &format!("tilt {tilt}"), 0.0, &mut rng);
    }
    // From z = 0 (radius 1) to z = 1 (radius 0.5), its apex at z = 2; a
    // plane through the apex, its normal 0.3 rad from square to the axis,
    // crosses it in the rulings at θ with sin θ = 2 tan 0.3, the slab's
    // face 0.6 wide about the first of them halfway up.
    let short = turned_solid(
        Frame::WORLD,
        DVec3::ZERO,
        &[(0.0, 0.0), (1.0, 0.0), (0.5, 1.0), (0.0, 1.0)],
        TAU,
    );
    let n = DVec3::new(0.0, f64::cos(0.3), f64::sin(0.3));
    let theta = (2.0 * f64::tan(0.3)).asin();
    let at = DVec3::new(0.75 * theta.cos(), 0.75 * theta.sin(), 0.5);
    for side in [1.0, -1.0] {
        let tool = slab(at, n * side, 0.6, 3.0);
        identities(
            &mut t,
            &short,
            &tool,
            &format!("apex {side}"),
            0.0,
            &mut rng,
        );
    }
    assert!(t.refused.is_empty(), "{:?}", t.refused);
    assert_eq!((t.checked, t.overlapping), (5, 5));
    report("small tools across large cones", &t);
}

/// Sections of a frustum by boxes whose edge runs along its axis, or a
/// few tolerances off it, past the apex: the box's faces through the axis
/// meet the cone in two rulings from the apex, beyond the frustum, and the
/// box's edge crosses them there (a section view of a countersink, whose
/// cone's apex lies in the hole). Every operation was refused, as the
/// crossing lay on both rulings, and one panicked, until the rulings were
/// taken from where they come near the region where both faces can be
/// (design note 0012 §15).
#[test]
fn sections_of_a_frustum_by_boxes_along_its_axis() {
    let mut t = Tally::default();
    let mut rng = Rng(0xA9E_0FF);
    // From z = 0 (radius 1) to z = 1 (radius 0.5), its apex at z = 2.
    let frustum = turned_solid(
        Frame::WORLD,
        DVec3::ZERO,
        &[(0.0, 0.0), (1.0, 0.0), (0.5, 1.0), (0.0, 1.0)],
        TAU,
    );
    let whole = PI * 1.75 / 3.0;
    // The box from (x, y) on, past the apex or short of it (top 3 or
    // 1.5); a half section from x = −3. Its volume in the frustum, to
    // first order in the offsets, which move faces across the frustum's
    // section through the axis, 0.75 in area each side of the axis.
    for (x, y, top) in [
        (0.0, 0.0, 3.0),
        (0.0, 0.0, 1.5),
        (-3.0, 0.0, 3.0),
        (0.0, -5e-7, 3.0),
        (0.0, -8e-7, 3.0),
        (1e-6, 2e-6, 3.0),
        (-3.0, 4e-5, 3.0),
    ] {
        let block = cuboid(DVec3::new(x, y, -1.0), DVec3::new(3.0, 3.0, top));
        let inside = if x < -1.0 {
            0.5 * whole - 1.5 * y
        } else {
            0.25 * whole - 0.75 * (x + y)
        };
        let i = run(&frustum, &block, BooleanOp::Intersect);
        let d = run(&frustum, &block, BooleanOp::Subtract);
        let u = run(&frustum, &block, BooleanOp::Union);
        let case = format!("({x}, {y}) to {top}");
        assert!((volume(&i) - inside).abs() < 1e-9, "{case}: {}", volume(&i));
        assert!(
            (volume(&d) - (whole - inside)).abs() < 1e-9,
            "{case}: {}",
            volume(&d)
        );
        assert!(
            (volume(&u) - (volume(&block) + whole - inside)).abs() < 1e-9,
            "{case}: {}",
            volume(&u)
        );
        identities(&mut t, &frustum, &block, &case, 0.0, &mut rng);
    }
    assert!(t.refused.is_empty(), "{:?}", t.refused);
    report("sections of a frustum by boxes along its axis", &t);
}

/// Sections of a countersunk plate through the hole's axis: halves along
/// the cone's seams and square to them and turned 30°, and a quarter.
/// Each plane passes through the countersink's apex, in the hole below
/// the cone's face, where the region's ball (but not the common box of
/// the faces) reached it, and every half section was refused as a plane
/// through a pointed cone's apex until the apex rules asked whether faces
/// on both surfaces can reach it (design note 0012 §15).
#[test]
fn sections_through_a_countersinks_axis() {
    let plate = cuboid(DVec3::new(-2.0, -2.0, 0.0), DVec3::new(2.0, 2.0, 1.0));
    // A hole of radius 0.3 with a 90° countersink to radius 0.6 at the
    // top: its cone's apex at z = 0.4.
    let tool = turned_solid(
        Frame::WORLD,
        DVec3::ZERO,
        &[
            (0.0, -1.0),
            (0.3, -1.0),
            (0.3, 0.7),
            (0.6, 1.0),
            (0.6, 2.0),
            (0.0, 2.0),
        ],
        TAU,
    );
    let sunk = run(&plate, &tool, BooleanOp::Subtract);
    let whole = 16.0 - PI * 0.09 * 0.7 - PI * 0.1 * (0.09 + 0.18 + 0.36);
    assert!((volume(&sunk) - whole).abs() < 1e-9, "{}", volume(&sunk));
    let turned_half = |angle: f64| {
        let x = DVec3::new(angle.cos(), angle.sin(), 0.0);
        let f = Frame::from_normal(DVec3::new(0.0, 0.0, -3.0) - x * 3.0, DVec3::Z, x);
        k().make_box(&f, DVec3::new(6.0, 3.0, 6.0)).expect("a box")
    };
    let mut t = Tally::default();
    let mut rng = Rng(0xC5_0A15);
    for (name, block, part) in [
        ("along the seams", turned_half(0.0), 0.5),
        ("square to them", turned_half(FRAC_PI_2), 0.5),
        ("turned 30°", turned_half(PI / 6.0), 0.5),
        (
            "a quarter",
            cuboid(DVec3::new(0.0, 0.0, -3.0), DVec3::new(3.0, 3.0, 3.0)),
            0.25,
        ),
    ] {
        let i = run(&sunk, &block, BooleanOp::Intersect);
        let d = run(&sunk, &block, BooleanOp::Subtract);
        let u = run(&sunk, &block, BooleanOp::Union);
        assert!((volume(&i) - part * whole).abs() < 1e-9, "{name}");
        assert!((volume(&d) - (1.0 - part) * whole).abs() < 1e-9, "{name}");
        let expected = volume(&block) + (1.0 - part) * whole;
        assert!((volume(&u) - expected).abs() < 1e-9, "{name}");
        identities(&mut t, &sunk, &block, name, 0.0, &mut rng);
    }
    assert!(t.refused.is_empty(), "{:?}", t.refused);
    report("sections through a countersink's axis", &t);
}

/// Frustums leaning along (1, 1, 1), their apexes beyond their tops, cut
/// by planes through the apex: one square to x, at a slant to the axis,
/// and one along the axis. The faces' boxes, square to the world's axes,
/// hold the apex beside a leaning axis although no face comes near it, and
/// these were refused as planes through a pointed cone's apex (found by
/// review) until the apex rules asked how far above its apex the cone's
/// faces stay along its axis (design note 0012 §15). A plane through the
/// apex cuts a cone into cones over the pieces of its base, so the piece
/// of the frustum is (1 − s³) / 3 · H · (the base disc's segment), the
/// top's cone a copy of the whole scaled by s about the apex.
#[test]
fn leaning_frustums_cut_through_their_apexes() {
    let u = DVec3::ONE.normalize();
    let axis = Frame::from_normal(DVec3::ZERO, u, DVec3::X);
    let mut t = Tally::default();
    let mut rng = Rng(0x1EA_4A9E);
    for top in [0.5, 0.25] {
        // Radius 1 at the base to `top` half a unit along the axis.
        let frustum = turned_solid(
            axis,
            DVec3::ZERO,
            &[(0.0, 0.0), (1.0, 0.0), (top, 0.5), (0.0, 0.5)],
            TAU,
        );
        let whole = PI * 0.5 * (1.0 + top + top * top) / 3.0;
        assert!((volume(&frustum) - whole).abs() < 1e-12);
        let height = 0.5 / (1.0 - top);
        let apex = u * height;
        let s = (height - 0.5) / height;
        // The plane x = apex.x crosses the base's plane `d` from its centre.
        let d = apex.x / (1.0 - u.x * u.x).sqrt();
        let segment = d.acos() - d * (1.0 - d * d).sqrt();
        let beyond = (1.0 - s * s * s) / 3.0 * height * segment;
        let across = Frame::from_normal(apex, DVec3::new(1.0, -1.0, 0.0).normalize(), u);
        let half = k()
            .make_box(
                &Frame {
                    origin: apex - (across.x + across.y) * 3.0,
                    ..across
                },
                DVec3::splat(6.0),
            )
            .expect("a box");
        for (name, block, part) in [
            (
                "square to x",
                cuboid(DVec3::new(apex.x, -3.0, -3.0), DVec3::splat(3.0)),
                beyond,
            ),
            ("along the axis", half, 0.5 * whole),
        ] {
            let case = format!("top {top}, {name}");
            let i = run(&frustum, &block, BooleanOp::Intersect);
            let dif = run(&frustum, &block, BooleanOp::Subtract);
            let uni = run(&frustum, &block, BooleanOp::Union);
            assert!((volume(&i) - part).abs() < 1e-9, "{case}: {}", volume(&i));
            assert!((volume(&dif) - (whole - part)).abs() < 1e-9, "{case}");
            let expected = volume(&block) + whole - part;
            assert!((volume(&uni) - expected).abs() < 1e-9, "{case}");
            identities(&mut t, &frustum, &block, &case, 0.0, &mut rng);
        }
    }
    assert!(t.refused.is_empty(), "{:?}", t.refused);
    report("leaning frustums cut through their apexes", &t);
}

/// A slender pointed cone (2° half-angle) cut by planes along its axis a
/// few tenths of a millimetre from it: the hyperbola comes no nearer the
/// axis than the plane, thousands of tolerances, but passes within
/// 0.37 mm of the apex, where the cone is within 128 tolerances of its
/// axis. The rule refusing curves that turn by the apex took that as too
/// near (found by review); it now asks how near the curve comes to the
/// axis. The volumes are `beyond` on the cone scaled to a base radius of
/// one half (x and y by s, volumes by s²).
#[test]
fn planes_beside_a_slender_cones_apex() {
    let slope = 2f64.to_radians().tan();
    let cone = turned_solid(
        Frame::WORLD,
        DVec3::ZERO,
        &[(0.0, 0.0), (slope, 0.0), (0.0, 1.0)],
        TAU,
    );
    let whole = cone_volume(slope, 1.0);
    let s = 0.5 / slope;
    let mut t = Tally::default();
    let mut rng = Rng(0x5_1E4D);
    for d in [2e-4, 3e-4, 3.5e-4] {
        let block = cuboid(DVec3::new(d, -1.0, -1.0), DVec3::new(1.0, 1.0, 2.0));
        let past = beyond(d * s, DVec3::NEG_X) / (s * s);
        let i = run(&cone, &block, BooleanOp::Intersect);
        let c = run(&cone, &block, BooleanOp::Subtract);
        assert!(
            (volume(&i) - past).abs() < 1e-12,
            "{d}: {} {past}",
            volume(&i)
        );
        assert!(
            (volume(&c) - (whole - past)).abs() < 1e-12,
            "{d}: {} {}",
            volume(&c),
            whole - past
        );
        identities(&mut t, &cone, &block, &format!("{d}"), 0.0, &mut rng);
    }
    // A plane 0.1 tolerances nearer the axis than the rule's 128 is refused.
    let near = cuboid(DVec3::new(1.27e-5, -1.0, -1.0), DVec3::new(1.0, 1.0, 2.0));
    match k().boolean(&cone, &near, BooleanOp::Intersect) {
        Err(KernelError::Unsupported(m)) => {
            assert!(m.contains("near a cone's apex"), "{m}");
            assert!(m.contains("design note 0012 §15"), "{m}");
        }
        r => panic!("{r:?}"),
    }
    assert!(t.refused.is_empty(), "{:?}", t.refused);
    assert_eq!((t.checked, t.overlapping), (3, 3));
    report("planes beside a slender cone's apex", &t);
}

/// A box whose face crosses a pointed cone 16 µm below its apex, square
/// to the axis or tilted by 1e-5: a circle (or an ellipse) 8 µm round
/// the axis, where the validator holds a face's edges off its singular
/// point by 100 tolerances. Both are refused with a message (found by
/// review: the square one made an invalid body, caught only by the
/// kernel's own check of its result).
#[test]
fn a_box_across_a_cone_just_below_its_apex_is_refused() {
    let cone = turned_solid(
        Frame::WORLD,
        DVec3::ZERO,
        &[(0.0, 0.0), (0.5, 0.0), (0.0, 1.0)],
        TAU,
    );
    let mut t = Tally::default();
    let mut rng = Rng(0x0A9E_C1C1);
    let z = 1.0 - 1.6e-5;
    let square = cuboid(DVec3::new(-1.0, -1.0, z), DVec3::new(1.0, 1.0, 2.0));
    let tilted = slab(
        DVec3::new(0.0, 0.0, z),
        DVec3::new(1e-5, 0.0, 1.0).normalize(),
        2.0,
        1.0,
    );
    identities(&mut t, &cone, &square, "square", 0.0, &mut rng);
    identities(&mut t, &cone, &tilted, "tilted", 0.0, &mut rng);
    assert_eq!(t.refused.len(), 8, "{:?}", t.refused);
    assert!(
        t.refused.iter().all(|m| m.contains("axis")),
        "{:?}",
        t.refused
    );
    report("a box across a cone just below its apex", &t);
}

/// A box whose face runs a few tolerances off a turned rod's axis: its
/// edge across the rod's flat ends passes their centres, singular points
/// of the ends' surfaces of revolution, within the hundred tolerances the
/// validator holds edges off them. Where a result keeps such an edge it
/// is refused with a message naming them (it was an invalid result,
/// caught only by the kernel's own check); other results are exact, and
/// a face 200 tolerances off works.
#[test]
fn a_box_a_hair_off_a_turned_rods_axis() {
    let rod = turned_solid(
        Frame::WORLD,
        DVec3::ZERO,
        &[(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)],
        TAU,
    );
    let mut refused = 0;
    for y in [5e-7, 2e-6, 5e-6, 2e-5] {
        let block = cuboid(DVec3::new(-3.0, y, -3.0), DVec3::new(3.0, 3.0, 3.0));
        let part = segment(1.0, y);
        for (op, expected) in [
            (BooleanOp::Union, volume(&block) + PI - part),
            (BooleanOp::Intersect, part),
            (BooleanOp::Subtract, PI - part),
        ] {
            match k().boolean(&rod, &block, op) {
                Ok(r) if y < 1e-5 => {
                    assert_eq!(k().validate(&r), Vec::<String>::new());
                    assert!((volume(&r) - expected).abs() < 1e-9, "{y} {op:?}");
                }
                Err(KernelError::Unsupported(m)) if y < 1e-5 => {
                    assert!(m.contains("a turned face's centre"), "{y} {op:?}: {m}");
                    refused += 1;
                }
                _ => {
                    let r = run(&rod, &block, op);
                    assert!((volume(&r) - expected).abs() < 1e-9, "{y} {op:?}");
                }
            }
        }
    }
    assert!(refused >= 3, "{refused}");
}

/// A ball halved by a box whose face runs through its axis. At the angle
/// of the ball's seams (the plane of the half disc it was turned from,
/// either way round) the box's edges and the seams meet at the poles
/// within the tolerance, and the operation is refused as a
/// near-coincidence there; at other angles the halves are exact.
#[test]
fn a_ball_halved_through_its_axis() {
    let b = ball(DVec3::ZERO, 1.0);
    let half = 2.0 * PI / 3.0;
    for degrees in [0.0f64, 30.0, 45.0, 90.0, 135.0, 180.0] {
        let (s, c) = degrees.to_radians().sin_cos();
        let x = DVec3::new(c, s, 0.0);
        let y = DVec3::new(-s, c, 0.0);
        let frame = Frame {
            origin: -3.0 * x - 3.0 * DVec3::Z,
            x,
            y,
            z: DVec3::Z,
        };
        let block = k()
            .make_box(&frame, DVec3::new(6.0, 3.0, 6.0))
            .expect("a box");
        let seam = degrees % 180.0 == 0.0;
        for (op, expected) in [
            (BooleanOp::Union, volume(&block) + half),
            (BooleanOp::Intersect, half),
            (BooleanOp::Subtract, half),
        ] {
            if seam {
                match k().boolean(&b, &block, op) {
                    Err(KernelError::InvalidInput(m)) => {
                        assert!(m.contains("within the tolerance"), "{degrees} {op:?}: {m}");
                    }
                    Err(e) => panic!("{degrees} {op:?}: {e}"),
                    Ok(_) => panic!("{degrees} {op:?}: done (update STATUS.md)"),
                }
            } else {
                let r = run(&b, &block, op);
                assert!((volume(&r) - expected).abs() < 1e-9, "{degrees} {op:?}");
            }
        }
    }
}

#[test]
fn volume_identities_on_planes_across_cones() {
    let mut t = Tally::default();
    let mut r = Rng(0x0C0A_1C5E);
    let mut points = Rng(0x5A3F_0E15);
    for case in 0..scale(100) {
        let (a, b) = random_conic_pair(&mut r);
        identities(&mut t, &a, &b, &format!("pair {case}"), 0.0, &mut points);
    }
    report("random planes across cones", &t);
}

/// Pairs the fuzz of planes across cones found faulty at fifty times its
/// size (2026-10-15), with their operands as drawn: a pointed cone and a
/// slab through or by its apex. Meshing the results (and the kernel's own
/// check of them) refused "a triangle faces against the surface where it
/// is too small to split or flip": needles at the apex, their corners on
/// nearly one ruling, in a line within the model tolerance, whose side
/// rounding decides (design note 0002 §4). Every operation must succeed,
/// valid, with the volume identities and the membership of sampled points.
#[test]
fn pairs_a_fuzz_of_cones_found() {
    let axis = |x: [f64; 3], y: [f64; 3], z: [f64; 3]| Frame {
        origin: DVec3::ZERO,
        x: DVec3::from(x),
        y: DVec3::from(y),
        z: DVec3::from(z),
    };
    let pairs = [
        (
            turned_solid(
                axis(
                    [
                        -0.8048554945410071,
                        0.5870975376206801,
                        -0.08674165220287786,
                    ],
                    [
                        -0.08192847943355583,
                        -0.25467734399089315,
                        -0.9635492591017074,
                    ],
                    [
                        -0.5877885309912539,
                        -0.7684113037803896,
                        0.25307846818262764,
                    ],
                ),
                DVec3::ZERO,
                &[(0.0, -0.75), (0.75, -0.75), (0.0, 0.5)],
                TAU,
            ),
            slab(
                DVec3::new(
                    -0.3157350251667678,
                    0.011842351692448017,
                    -0.16814781541494764,
                ),
                DVec3::new(
                    -0.3381028400462773,
                    -0.5737269166800482,
                    -0.7460052912878341,
                ),
                0.8650562169587592,
                3.0,
            ),
        ),
        (
            turned_solid(
                axis(
                    [-0.796864608307113, -0.577118279788069, 0.17872125548465564],
                    [-0.320602570501996, 0.6546701098179141, 0.6845590106765911],
                    [-0.5120749826100316, 0.4882023539941819, -0.7067090446141717],
                ),
                DVec3::ZERO,
                &[
                    (0.0, -0.3543928295055908),
                    (0.8862456149308802, -0.3543928295055908),
                    (0.0, 0.38345553452980197),
                ],
                TAU,
            ),
            slab(
                DVec3::new(
                    -0.22598923947435284,
                    0.24327038384263963,
                    0.18347693192756537,
                ),
                DVec3::new(0.9001119687661934, 0.4356305875738727, 0.004943162336990964),
                0.8521832505725762,
                3.0,
            ),
        ),
        (
            turned_solid(
                axis(
                    [
                        0.9870518578515366,
                        -0.06995171628194835,
                        0.14434468227489314,
                    ],
                    [0.1561623788948448, 0.6245996115040233, -0.765172292184529],
                    [-0.03663251738186855, 0.7778059415421259, 0.6274360333707607],
                ),
                DVec3::ZERO,
                &[(0.0, -0.75), (0.25, -0.75), (0.0, 0.75)],
                TAU,
            ),
            slab(
                DVec3::new(-0.0658854127382031, 0.3498087349890746, 0.2556507156004063),
                DVec3::new(0.15457599309365386, 0.6551435845970321, -0.7395222416671876),
                0.5891312740501218,
                3.0,
            ),
        ),
    ];
    let mut t = Tally::default();
    let mut points = Rng(0x5A3F_0E17);
    for (i, (a, b)) in pairs.iter().enumerate() {
        identities(&mut t, a, b, &format!("pair {}", i + 1), 0.0, &mut points);
    }
    assert!(t.refused.is_empty(), "{:?}", t.refused);
    assert!(t.bugs.is_empty(), "{:?}", t.bugs);
    assert_eq!(t.checked, 3);
}

/// The membership check's rays against its winding numbers, at points
/// near and far from random solids: every ray that answers gives the
/// winding number's answer, and the winding numbers are whole (the meshes
/// are closed and face out).
#[test]
fn rays_and_winding_numbers_agree() {
    let mut r = Rng(0x0BA1_1C0F);
    let (mut points, mut rays, mut answered) = (0, 0, 0);
    for case in 0..16 {
        let axis = turned(&mut r, DVec3::ZERO);
        let body = if case % 2 == 0 {
            random_revolve(&mut r, axis).0
        } else {
            random_curved(&mut r, axis)
        };
        let m = Mesh::of(&body).expect("a mesh");
        let (lo, hi) = (m.nodes[0].lo - 0.2, m.nodes[0].hi + 0.2);
        for _ in 0..200 {
            let p = DVec3::new(
                r.range(lo.x, hi.x),
                r.range(lo.y, hi.y),
                r.range(lo.z, hi.z),
            );
            if m.near(p, 2.0 * CHORD) {
                continue;
            }
            points += 1;
            let w = m.winding(p);
            assert!(
                w.abs() < 1e-9 || (w - 1.0).abs() < 1e-9,
                "solid {case}: winding number {w} at {p:?}"
            );
            for d in RAYS {
                rays += 1;
                if let Some(n) = m.crossings(p, d.normalize()) {
                    answered += 1;
                    assert_eq!(n % 2 == 1, w > 0.5, "solid {case}: {n} crossings at {p:?}");
                }
            }
        }
    }
    assert!(points > 2000, "{points} points");
    assert!(
        answered + rays / 1000 >= rays,
        "{answered} of {rays} rays answered"
    );
}

/// Frustums cut by planes within a ten-thousandth of a radian of tangent
/// along a ruling, through the apex or a few tolerances off it, the top
/// a hundredth to a tenth of a percent of the height below the apex (found
/// by review, in a sweep of 8,000 such cases: 474 unions gave two solids a
/// few hundredths of a micrometre apart, where a rim lay within the
/// tolerance of the plane along the ruling, and 582 operations failed the
/// kernel's own check of its result, where the plane cut a sliver thinner
/// than the tolerance from the top, a rim touched the plane where a seam
/// meets it, or the plane's ellipse grazed the slab's edge metres from the
/// frustum). Each is now done or refused with a reason; a union of bodies
/// within the tolerance of each other is one solid or refused.
#[test]
fn planes_nearly_tangent_to_a_frustum_by_its_apex() {
    let eps = 1e-7;
    let axes = [
        Frame::WORLD,
        Frame::from_normal(DVec3::ZERO, DVec3::ONE.normalize(), DVec3::X),
    ];
    // (axis, slope, gap below the apex as a part of the height, azimuth,
    // tilt from tangent as a part of the cone's angle, offset in
    // tolerances): the cases the review found, then a sweep.
    let mut cases = vec![
        (0, 1.0, 1e-3, 0.0, 1e-6 / FRAC_PI_4, 0.0),
        (0, 1.0, 1e-3, 0.4, -1e-5, 0.0),
        (0, 1.0, 1e-3, 0.4, -1e-6, 0.0),
        (0, 1.0, 1e-4, 0.4, -1e-7, 0.5),
        (0, 2.0, 1e-3, 0.0, -1e-7, 2.0),
        (0, 2.0, 1e-3, PI, -1e-7, 2.0),
        (0, 10.0, 1e-3, 0.4, 1e-7, -10.0),
        (1, 1.0, 1e-3, 0.0, 1e-4, -2.0),
    ];
    for slope in [1.0, 10.0] {
        for phi in [0.0, 0.4] {
            for f in [-1e-7, 1e-7, 1e-5] {
                for d in [0.0, 2.0, -2.0, -10.0] {
                    cases.push((1, slope, 1e-3, phi, f, d));
                }
            }
        }
    }
    let mut t = Tally::default();
    let mut rng = Rng(0x7A_4E47);
    let mut lumps = Vec::new();
    for (i, &(ai, slope, gap, phi, f, d)) in cases.iter().enumerate() {
        let axis = axes[ai];
        let za = 1.0 / slope;
        let zt = za * (1.0 - gap);
        let rt = 1.0 - slope * zt;
        let frustum = turned_solid(
            axis,
            axis.origin,
            &[(0.0, 0.0), (1.0, 0.0), (rt, zt), (0.0, zt)],
            TAU,
        );
        let apex = axis.origin + axis.z * za;
        let e = axis.x * phi.cos() + axis.y * phi.sin();
        let psi = slope.atan() * (1.0 + f);
        let n = e * psi.cos() + axis.z * psi.sin();
        let q = apex + n * (d * eps);
        let block = slab(q, n, 12.0, 12.0);
        let case =
            format!("case {i}: slope {slope}, gap {gap}, azimuth {phi}, tilt {f:e}, offset {d}");
        // How far the frustum reaches toward the slab: from its rims.
        let reach = [(0.0, 1.0), (zt, rt)]
            .map(|(z, rho): (f64, f64)| {
                let w = n - axis.z * n.dot(axis.z);
                (axis.origin + axis.z * z - q).dot(n) + rho * w.length()
            })
            .into_iter()
            .fold(f64::MIN, f64::max);
        if let Ok(r) = k().boolean(&frustum, &block, BooleanOp::Union) {
            let solids = k().topology(&r).solids;
            if solids > 1 && reach >= -eps {
                lumps.push(format!("{case}: {solids} solids, {reach:e} apart"));
            }
        }
        if i == 0 {
            // A plane through the apex cutting a sliver thinner than the
            // tolerance from the top: refused.
            for op in [BooleanOp::Union, BooleanOp::Intersect, BooleanOp::Subtract] {
                let r = k().boolean(&frustum, &block, op);
                assert!(
                    matches!(
                        r,
                        Err(KernelError::Unsupported(_) | KernelError::InvalidInput(_))
                    ),
                    "{op:?}: {r:?}"
                );
            }
        }
        identities(&mut t, &frustum, &block, &case, 0.0, &mut rng);
    }
    assert!(lumps.is_empty(), "{lumps:#?}");
    assert!(t.checked > 0);
    report("planes nearly tangent to a frustum by its apex", &t);
}
