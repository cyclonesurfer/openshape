//! Booleans of solids bounded by planes and cylinders (design note 0012
//! §10.7, item 2): golden cases with exact volumes, areas and topology.

use std::f64::consts::{FRAC_PI_2, PI, TAU};

use opendatum_kernel::{
    Body, BooleanOp, BrepKernel, DVec2, DVec3, EdgeOrigin, ExtrudeSpec, FaceIx, FaceOrigin, Frame,
    Kernel, KernelError, Merged, Operand, Profile, ProfileCurve, RevolveSpec, TessellationOptions,
    VertexOrigin,
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
    let near = 3.0 * k().tolerance().linear;
    let opts = TessellationOptions::default();
    let mesh = |x: &Body| k().tessellate(x, &opts).map_err(|e| e.to_string());
    let (ma, mb, mr) = (mesh(a)?, mesh(b)?, mesh(r)?);
    let (ja, jb) = (k().adjacency(a), k().adjacency(b));
    let body = |o: Operand| match o {
        Operand::First => (a, &ma, &ja),
        Operand::Second => (b, &mb, &jb),
    };
    // How far `p` is from the surface of face `f` of `x`, and the outward
    // normal there.
    let on = |x: &Body, f: FaceIx, p: DVec3| -> Result<(f64, DVec3), String> {
        let t = k()
            .tangent_plane(x, f, p)
            .ok_or_else(|| format!("face {} has no tangent plane", f.0))?;
        Ok(((p - t.point).length(), t.normal))
    };
    let o = k().origins(r);
    let face = |i: usize, origin: &FaceOrigin| -> Result<(), String> {
        let FaceOrigin::Operand { operand, face } = *origin else {
            return Err(format!("face {i}: origin {origin:?}"));
        };
        let (x, _, _) = body(operand);
        let flip = op == BooleanOp::Subtract && operand == Operand::Second;
        for p in &mr.faces[i].positions {
            let (d, n) = on(x, face, *p)?;
            if d > near {
                return Err(format!(
                    "face {i} leaves the surface of {operand:?} face {} by {d:e}",
                    face.0
                ));
            }
            let (_, m) = on(r, FaceIx(i as u32), *p)?;
            let facing = m.dot(n);
            if (facing < 0.0) != flip || facing.abs() < 1.0 - 1e-6 {
                return Err(format!(
                    "face {i} does not face as {operand:?} face {} does ({facing})",
                    face.0
                ));
            }
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
        for p in &mr.edges[i] {
            for &(x, f) in &faces {
                let (d, _) = on(x, f, *p)?;
                if d > near {
                    return Err(format!(
                        "edge {i} ({origin:?}) leaves the surface of face {} by {d:e}",
                        f.0
                    ));
                }
            }
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
/// result, else the result and its volume.
type Outcome = Option<Option<(Body, f64)>>;

/// A boolean's outcome; refusals and kernel bugs are recorded.
fn measured(t: &mut Tally, a: &Body, b: &Body, op: BooleanOp, case: &str) -> Outcome {
    match k().boolean(a, b, op) {
        Ok(r) => {
            let issues = k().validate(&r);
            if !issues.is_empty() {
                t.bugs
                    .push(format!("{case} {op:?}: invalid result: {issues:?}"));
                return None;
            }
            if let Err(e) = origins_agree(a, b, &r, op) {
                t.bugs.push(format!("{case} {op:?}: origins: {e}"));
                return None;
            }
            match k().mass_properties(&r, 1e-12) {
                Ok(m) => Some(Some((r, m.volume))),
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

/// A body's surface as triangles within `chord` of it, for deciding
/// whether points well away from it are inside: a check of the booleans
/// that shares none of their code (the volume identities hold for a result
/// that misses where the bodies meet altogether, as if they were apart).
struct Mesh {
    tris: Vec<[DVec3; 3]>,
    chord: f64,
}

impl Mesh {
    fn of(b: &Body, chord: f64) -> Result<Mesh, String> {
        let opts = TessellationOptions {
            chord_tolerance: chord,
            max_angle: 0.5,
        };
        let (pos, tris) = k()
            .tessellate(b, &opts)
            .map_err(|e| e.to_string())?
            .welded();
        Ok(Mesh {
            tris: tris.iter().map(|t| t.map(|i| pos[i as usize])).collect(),
            chord,
        })
    }

    /// Whether `p` is inside, if it is farther than `margin` from the
    /// mesh (so the mesh and the exact surface agree there): by the winding
    /// number, the sum of the triangles' solid angles about `p` (A. van
    /// Oosterom and J. Strackee, "The solid angle of a plane triangle",
    /// IEEE Trans. Biomed. Eng. 30, 1983).
    fn inside(&self, p: DVec3, margin: f64) -> Option<bool> {
        let mut total = 0.0;
        for t in &self.tris {
            if triangle_distance(p, t) <= margin + self.chord {
                return None;
            }
            let [a, b, c] = t.map(|q| q - p);
            let (la, lb, lc) = (a.length(), b.length(), c.length());
            let num = a.dot(b.cross(c));
            let den = la * lb * lc + a.dot(b) * lc + b.dot(c) * la + c.dot(a) * lb;
            total += 2.0 * num.atan2(den);
        }
        Some(total / (4.0 * PI) > 0.5)
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
    a: &Body,
    b: &Body,
    results: &[(BooleanOp, bool, &Outcome)],
    case: &str,
    rng: &mut Rng,
) {
    const CHORD: f64 = 1e-3;
    const MARGIN: f64 = 1e-3;
    let mesh = |x: &Body| Mesh::of(x, CHORD);
    let (ma, mb) = match (mesh(a), mesh(b)) {
        (Ok(x), Ok(y)) => (x, y),
        (Err(e), _) | (_, Err(e)) => {
            t.bugs.push(format!("{case}: meshing an operand: {e}"));
            return;
        }
    };
    let bounds = |m: &Mesh| {
        m.tris.iter().flatten().fold(
            (DVec3::splat(f64::MAX), DVec3::splat(f64::MIN)),
            |(lo, hi), q| (lo.min(*q), hi.max(*q)),
        )
    };
    let ((la, ha), (lb, hb)) = (bounds(&ma), bounds(&mb));
    let (lo, hi) = if la.cmple(hb).all() && lb.cmple(ha).all() {
        (la.max(lb), ha.min(hb))
    } else {
        (la.min(lb), ha.max(hb))
    };
    let mut meshes = Vec::new();
    for &(op, swapped, r) in results {
        match r {
            None => {}
            Some(None) => meshes.push((op, swapped, None)),
            Some(Some((body, _))) => match mesh(body) {
                Ok(m) => meshes.push((op, swapped, Some(m))),
                Err(e) => t
                    .bugs
                    .push(format!("{case} {op:?}: meshing the result: {e}")),
            },
        }
    }
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
    let ru = measured(t, a, b, BooleanOp::Union, case);
    let ri = measured(t, a, b, BooleanOp::Intersect, case);
    let rab = measured(t, a, b, BooleanOp::Subtract, case);
    let rba = measured(t, b, a, BooleanOp::Subtract, case);
    let vol = |r: &Outcome| r.as_ref().map(|x| x.as_ref().map_or(0.0, |y| y.1));
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
        a,
        b,
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
                Err(KernelError::Unsupported(m)) => assert!(m.contains("quartic"), "{op:?}: {m}"),
                other => panic!("({y}, {z}) {op:?}: {:?}", other.map(|b| volume(&b))),
            }
        }
    }
}

/// Turned rods and tubes (their caps planes recorded as surfaces of
/// revolution, in two halves each) cut by boxes: a cap kept whole is one
/// face again, recorded as a plane, and the volumes are exact.
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
