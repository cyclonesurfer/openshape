//! Booleans of solids bounded by planes (design note 0012 §8 step 1, ADR
//! 0020): golden cases with exact volumes and topology, the volume
//! identities on random pairs, and motions by fractions of the tolerance.

use std::f64::consts::TAU;

mod self_contact;

use opendatum_kernel::{
    Adjacency, Body, BodyMesh, BooleanOp, BrepKernel, DVec2, DVec3, EdgeOrigin, ExtrudeSpec,
    FaceIx, FaceOrigin, Frame, Kernel, KernelError, Merged, Operand, Profile, RevolveSpec,
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

/// Checks a result's origins against the operands' geometry: a face on
/// an operand face's plane, facing the same way (the other way for the
/// second operand of a difference); an edge along an operand edge, or on
/// both faces' planes where the bodies meet; a vertex at an operand
/// vertex. Within three tolerances, the most a cluster moves a point.
fn origins_agree(a: &Body, b: &Body, r: &Body, op: BooleanOp) -> Result<(), String> {
    let near = 3.0 * k().tolerance().linear;
    let opts = TessellationOptions::default();
    let mesh = |x: &Body| k().tessellate(x, &opts).map_err(|e| e.to_string());
    let (ma, mb, mr) = (mesh(a)?, mesh(b)?, mesh(r)?);
    let body = |o: Operand| match o {
        Operand::First => (a, &ma),
        Operand::Second => (b, &mb),
    };
    let plane = |x: &Body, f: u32| {
        k().face_plane(x, FaceIx(f))
            .ok_or_else(|| format!("face {f} has no plane"))
    };
    let o = k().origins(r);
    // An edge from one operand edge is a piece of it; an edge joined from
    // pieces of several lies on each one's line and reaches it.
    let joined: Vec<usize> = o
        .merged
        .iter()
        .filter_map(|m| match m {
            Merged::Edge(e, _) => Some(e.0 as usize),
            _ => None,
        })
        .collect();
    let along = |i: usize, p: &[DVec3], q: &[DVec3]| -> bool {
        let (q0, q1) = (q[0], q[q.len() - 1]);
        let d = (q1 - q0).normalize();
        let len = (q1 - q0).length();
        let t = |x: DVec3| (x - q0).dot(d);
        let line = |x: &DVec3| (*x - q0).reject_from_normalized(d).length() <= near;
        let inside = |x: &DVec3| t(*x) >= -near && t(*x) <= len + near;
        if joined.contains(&i) {
            let (p0, p1) = (t(p[0]), t(p[p.len() - 1]));
            p.iter().all(line) && p0.max(p1) >= -near && p0.min(p1) <= len + near
        } else {
            p.iter().all(|x| line(x) && inside(x))
        }
    };
    let face = |i: usize, origin: &FaceOrigin| -> Result<(), String> {
        let FaceOrigin::Operand { operand, face } = *origin else {
            return Err(format!("face {i}: origin {origin:?}"));
        };
        let (x, _) = body(operand);
        let (pr, po) = (plane(r, i as u32)?, plane(x, face.0)?);
        let flip = op == BooleanOp::Subtract && operand == Operand::Second;
        let facing = pr.normal.dot(po.normal);
        if (facing < 0.0) != flip {
            return Err(format!(
                "face {i} faces the wrong way from {operand:?} face {}",
                face.0
            ));
        }
        if facing.abs() < 1.0 - 1e-9 {
            return Err(format!(
                "face {i} is not parallel to {operand:?} face {}",
                face.0
            ));
        }
        for p in &mr.faces[i].positions {
            if (*p - po.point).dot(po.normal).abs() > near {
                return Err(format!(
                    "face {i} leaves the plane of {operand:?} face {}",
                    face.0
                ));
            }
        }
        Ok(())
    };
    let edge = |i: usize, origin: &EdgeOrigin| -> Result<(), String> {
        let p = &mr.edges[i];
        match *origin {
            EdgeOrigin::Operand { operand, edge } => {
                let (_, m) = body(operand);
                if !along(i, p, &m.edges[edge.0 as usize]) {
                    return Err(format!("edge {i} is not along {operand:?} edge {}", edge.0));
                }
            }
            EdgeOrigin::Meet { first, second } => {
                for (x, f) in [(a, first), (b, second)] {
                    let q = plane(x, f.0)?;
                    if p.iter().any(|x| (*x - q.point).dot(q.normal).abs() > near) {
                        return Err(format!("edge {i} leaves the plane of face {}", f.0));
                    }
                }
            }
            _ => return Err(format!("edge {i}: origin {origin:?}")),
        }
        Ok(())
    };
    let vertex = |i: usize, origin: &VertexOrigin| -> Result<(), String> {
        match *origin {
            VertexOrigin::Operand { operand, vertex } => {
                let (_, m) = body(operand);
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

fn refused(a: &Body, b: &Body, op: BooleanOp) -> String {
    match k().boolean(a, b, op) {
        Err(KernelError::InvalidInput(m)) => m,
        other => panic!("{op:?}: expected a refusal, got {other:?}"),
    }
}

fn close(x: f64, y: f64) -> bool {
    (x - y).abs() <= 1e-12 * (1.0 + x.abs().max(y.abs()))
}

#[test]
fn overlapping_cubes() {
    let a = cuboid(DVec3::ZERO, DVec3::ONE);
    let b = cuboid(DVec3::splat(0.5), DVec3::splat(1.5));
    let u = run(&a, &b, BooleanOp::Union);
    assert!(close(volume(&u), 1.875), "{}", volume(&u));
    assert!(close(area(&u), 10.5), "{}", area(&u));
    assert_eq!(counts(&u), (12, 30, 20));
    let i = run(&a, &b, BooleanOp::Intersect);
    assert!(close(volume(&i), 0.125), "{}", volume(&i));
    assert_eq!(counts(&i), (6, 12, 8));
    let s = run(&a, &b, BooleanOp::Subtract);
    assert!(close(volume(&s), 0.875), "{}", volume(&s));
    assert!(close(area(&s), 6.0), "{}", area(&s));
    assert_eq!(counts(&s), (9, 21, 14));
}

#[test]
fn stacked_boxes_merge_into_one() {
    let a = cuboid(DVec3::ZERO, DVec3::ONE);
    let b = cuboid(DVec3::new(0.0, 0.0, 1.0), DVec3::new(1.0, 1.0, 2.0));
    let u = run(&a, &b, BooleanOp::Union);
    assert!(close(volume(&u), 2.0));
    assert_eq!(counts(&u), (6, 12, 8));
    assert!(refused(&a, &b, BooleanOp::Intersect).contains("empty"));
    let s = run(&a, &b, BooleanOp::Subtract);
    assert!(close(volume(&s), 1.0));
    assert_eq!(counts(&s), (6, 12, 8));
}

#[test]
fn a_box_against_part_of_a_side() {
    let a = cuboid(DVec3::ZERO, DVec3::ONE);
    let b = cuboid(DVec3::new(1.0, 0.25, 0.0), DVec3::new(2.0, 0.75, 1.0));
    let u = run(&a, &b, BooleanOp::Union);
    assert!(close(volume(&u), 1.5));
    // A T-shaped prism: eight sides and two caps.
    assert_eq!(counts(&u), (10, 24, 16));
    assert!(close(area(&u), 2.0 * 1.5 + 6.0));
    let s = run(&a, &b, BooleanOp::Subtract);
    assert!(close(volume(&s), 1.0));
    assert_eq!(counts(&s), (6, 12, 8));
}

#[test]
fn flush_overlaps() {
    let a = cuboid(DVec3::ZERO, DVec3::new(2.0, 1.0, 1.0));
    let b = cuboid(DVec3::new(1.0, 0.0, 0.0), DVec3::new(3.0, 1.0, 1.0));
    for (op, v) in [
        (BooleanOp::Union, 3.0),
        (BooleanOp::Intersect, 1.0),
        (BooleanOp::Subtract, 1.0),
    ] {
        let r = run(&a, &b, op);
        assert!(close(volume(&r), v), "{op:?}: {}", volume(&r));
        assert_eq!(counts(&r), (6, 12, 8), "{op:?}");
    }
}

#[test]
fn a_slot_and_a_hole() {
    let a = cuboid(DVec3::ZERO, DVec3::new(4.0, 2.0, 1.0));
    let slot = cuboid(DVec3::new(1.0, -1.0, 0.5), DVec3::new(3.0, 3.0, 2.0));
    let s = run(&a, &slot, BooleanOp::Subtract);
    assert!(close(volume(&s), 6.0));
    assert_eq!(counts(&s), (10, 24, 16));
    let plate = cuboid(DVec3::ZERO, DVec3::new(3.0, 3.0, 1.0));
    let pin = cuboid(DVec3::new(1.0, 1.0, -1.0), DVec3::new(2.0, 2.0, 2.0));
    let h = run(&plate, &pin, BooleanOp::Subtract);
    assert!(close(volume(&h), 8.0));
    assert_eq!(counts(&h), (10, 24, 16));
    assert_eq!(k().topology(&h).genus, Some(1));
    let u = run(&plate, &pin, BooleanOp::Union);
    assert!(close(volume(&u), 9.0 + 2.0));
    let i = run(&plate, &pin, BooleanOp::Intersect);
    assert!(close(volume(&i), 1.0));
    assert_eq!(counts(&i), (6, 12, 8));
}

#[test]
fn nested_bodies() {
    let a = cuboid(DVec3::ZERO, DVec3::splat(3.0));
    let b = cuboid(DVec3::ONE, DVec3::splat(2.0));
    let s = run(&a, &b, BooleanOp::Subtract);
    assert!(close(volume(&s), 26.0));
    let t = k().topology(&s);
    assert_eq!((t.solids, t.shells, t.faces), (1, 2, 12));
    let u = run(&a, &b, BooleanOp::Union);
    assert!(close(volume(&u), 27.0));
    assert_eq!(counts(&u), (6, 12, 8));
    let i = run(&a, &b, BooleanOp::Intersect);
    assert!(close(volume(&i), 1.0));
    assert!(refused(&b, &a, BooleanOp::Subtract).contains("empty"));
}

#[test]
fn disjoint_and_identical_bodies() {
    let a = cuboid(DVec3::ZERO, DVec3::ONE);
    let b = cuboid(DVec3::splat(2.0), DVec3::splat(3.0));
    let u = run(&a, &b, BooleanOp::Union);
    assert_eq!(k().topology(&u).solids, 2);
    assert!(close(volume(&u), 2.0));
    assert!(refused(&a, &b, BooleanOp::Intersect).contains("empty"));
    assert!(close(volume(&run(&a, &b, BooleanOp::Subtract)), 1.0));
    for op in [BooleanOp::Union, BooleanOp::Intersect] {
        let r = run(&a, &a, op);
        assert!(close(volume(&r), 1.0));
        assert_eq!(counts(&r), (6, 12, 8));
    }
    assert!(refused(&a, &a, BooleanOp::Subtract).contains("empty"));
}

/// `a o b` is `a` unchanged, each of its faces, edges and vertices named
/// once by its origin in `a`, with nothing merged: the contact with `b`
/// leaves no trace (ADR 0024 §3: `b` keeps no face at the shared entity,
/// so its entities there are not origins).
fn the_first_operand(a: &Body, b: &Body, op: BooleanOp) {
    let r = run(a, b, op);
    assert_eq!(counts(&r), counts(a), "{op:?}");
    assert!(
        close(volume(&r), volume(a)) && close(area(&r), area(a)),
        "{op:?}"
    );
    let o = k().origins(&r);
    let mut faces: Vec<u32> = (o.faces.iter())
        .map(|g| match *g {
            FaceOrigin::Operand {
                operand: Operand::First,
                face,
            } => face.0,
            ref other => panic!("{op:?}: a face from {other:?}"),
        })
        .collect();
    let mut edges: Vec<u32> = (o.edges.iter())
        .map(|g| match *g {
            EdgeOrigin::Operand {
                operand: Operand::First,
                edge,
            } => edge.0,
            ref other => panic!("{op:?}: an edge from {other:?}"),
        })
        .collect();
    let mut vertices: Vec<u32> = (o.vertices.iter())
        .map(|g| match *g {
            VertexOrigin::Operand {
                operand: Operand::First,
                vertex,
            } => vertex.0,
            ref other => panic!("{op:?}: a vertex from {other:?}"),
        })
        .collect();
    faces.sort_unstable();
    edges.sort_unstable();
    vertices.sort_unstable();
    let (f, e, v) = counts(a);
    assert_eq!(faces, (0..f as u32).collect::<Vec<_>>(), "{op:?}");
    assert_eq!(edges, (0..e as u32).collect::<Vec<_>>(), "{op:?}");
    assert_eq!(vertices, (0..v as u32).collect::<Vec<_>>(), "{op:?}");
    assert_eq!(o.merged, Vec::<Merged>::new(), "{op:?}");
}

#[test]
fn bodies_touching_along_an_edge_or_at_a_point() {
    let a = cuboid(DVec3::ZERO, DVec3::ONE);
    let edge = cuboid(DVec3::new(1.0, 1.0, 0.0), DVec3::new(2.0, 2.0, 1.0));
    let corner = cuboid(DVec3::ONE, DVec3::splat(2.0));
    let along = "the result would be non-manifold: the bodies touch along an edge near \
                 (1000.000000, 1000.000000, 500.000000) mm";
    let at = "the result would be non-manifold: the bodies touch at a point near \
              (1000.000000, 1000.000000, 1000.000000) mm";
    for (b, why) in [(&edge, along), (&corner, at)] {
        assert_eq!(refused(&a, b, BooleanOp::Union), why);
        assert_eq!(refused(b, &a, BooleanOp::Union), why);
        assert_eq!(
            refused(&a, b, BooleanOp::Intersect),
            "the result would be empty"
        );
        assert_eq!(
            refused(b, &a, BooleanOp::Intersect),
            "the result would be empty"
        );
        // Before ADR 0024 the shared edge and its ends also named the other
        // cube's edge and vertices, which bound no face of the result.
        the_first_operand(&a, b, BooleanOp::Subtract);
        the_first_operand(b, &a, BooleanOp::Subtract);
    }
}

#[test]
fn golden_cases_at_millimetre_and_kilometre_scales() {
    for s in [1e-3, 1e3] {
        let a = cuboid(DVec3::ZERO, DVec3::ONE * s);
        let b = cuboid(DVec3::splat(0.5 * s), DVec3::splat(1.5 * s));
        let u = run(&a, &b, BooleanOp::Union);
        assert!(close(volume(&u), 1.875 * s * s * s), "{s}: {}", volume(&u));
        assert!(close(area(&u), 10.5 * s * s), "{s}: {}", area(&u));
        assert_eq!(counts(&u), (12, 30, 20), "{s}");
        let plate = cuboid(DVec3::ZERO, DVec3::new(3.0, 3.0, 1.0) * s);
        let pin = cuboid(
            DVec3::new(1.0, 1.0, -1.0) * s,
            DVec3::new(2.0, 2.0, 2.0) * s,
        );
        let h = run(&plate, &pin, BooleanOp::Subtract);
        assert!(close(volume(&h), 8.0 * s * s * s), "{s}: {}", volume(&h));
        assert_eq!(counts(&h), (10, 24, 16), "{s}");
        assert_eq!(k().topology(&h).genus, Some(1), "{s}");
    }
}

#[test]
fn a_plane_through_a_cones_apex_is_not_supported_yet() {
    // A cone: a right triangle revolved about its upright side.
    let cone = k()
        .revolve(
            &Profile::polygon(
                Frame::WORLD,
                &[DVec2::ZERO, DVec2::new(0.5, 0.0), DVec2::new(0.0, 1.0)],
            ),
            &RevolveSpec {
                axis_origin: DVec2::ZERO,
                axis_dir: DVec2::Y,
                angle: TAU,
            },
        )
        .unwrap();
    // A box whose side runs through the apex along the axis: the plane
    // crosses the cone in two rulings from its apex, which lies in faces
    // of both bodies, so the result's edges would meet at the cone's
    // singular point (design note 0012 §15). Planes across the cone
    // elsewhere meet it in conics, and planes through the apex of a cone
    // whose faces stop short of it in rays (tests/curved.rs).
    let cut = k()
        .make_box(
            &Frame {
                origin: DVec3::new(0.0, -1.0, -1.0),
                ..Frame::WORLD
            },
            DVec3::new(1.0, 3.0, 2.0),
        )
        .unwrap();
    for op in [BooleanOp::Union, BooleanOp::Intersect, BooleanOp::Subtract] {
        match k().boolean(&cone, &cut, op) {
            Err(KernelError::Unsupported(m)) => {
                assert!(m.contains("a plane through a cone's apex"), "{m}");
                assert!(m.contains("design note 0012 §14"), "{m}");
            }
            other => panic!("{:?}", other.map(|_| ())),
        }
    }
}

#[test]
fn origins_cover_the_result_and_name_the_operands() {
    let a = cuboid(DVec3::ZERO, DVec3::ONE);
    let b = cuboid(DVec3::splat(0.5), DVec3::splat(1.5));
    let u = run(&a, &b, BooleanOp::Union);
    let o = k().origins(&u);
    assert_eq!(o.faces.len(), 12);
    let from = |op| {
        o.faces
            .iter()
            .filter(|f| matches!(f, FaceOrigin::Operand { operand, .. } if *operand == op))
            .count()
    };
    assert_eq!((from(Operand::First), from(Operand::Second)), (6, 6));
    // The six edges where the cubes meet, and the six new vertices.
    let meets = o
        .edges
        .iter()
        .filter(|e| matches!(e, EdgeOrigin::Meet { .. }))
        .count();
    assert_eq!(meets, 6);
    let new = o
        .vertices
        .iter()
        .filter(|v| matches!(v, VertexOrigin::Meet))
        .count();
    assert_eq!(new, 6);
    // Stacked boxes: the merged sides carry both operands' faces.
    let top = cuboid(DVec3::new(0.0, 0.0, 1.0), DVec3::new(1.0, 1.0, 2.0));
    let m = run(&a, &top, BooleanOp::Union);
    let o = k().origins(&m);
    let merged_faces = o
        .merged
        .iter()
        .filter(|x| matches!(x, Merged::Face(..)))
        .count();
    assert_eq!(merged_faces, 4, "{:?}", o.merged);
}

// ---------------------------------------------------------------------
// Random pairs: the volume identities of design note 0012 §7.

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

    /// A value on a coarse grid (to make coincidences) or anywhere.
    fn coord(&mut self, lo: f64, hi: f64, snap: bool) -> f64 {
        let x = self.range(lo, hi);
        if snap {
            (x * 4.0).round() / 4.0
        } else {
            x
        }
    }
}

fn random_frame(r: &mut Rng, aligned: bool, centre: DVec3) -> Frame {
    if aligned {
        return Frame {
            origin: centre,
            ..Frame::WORLD
        };
    }
    let n = DVec3::new(r.range(-1.0, 1.0), r.range(-1.0, 1.0), r.range(-1.0, 1.0));
    let n = if n.length() < 0.1 { DVec3::Z } else { n };
    let hint = DVec3::new(r.range(-1.0, 1.0), r.range(-1.0, 1.0), r.range(-1.0, 1.0));
    Frame::from_normal(centre, n, hint)
}

/// A random polyhedron: a box (often on a grid), or an extruded
/// star-shaped polygon, sometimes with a hole.
fn random_body(r: &mut Rng) -> Body {
    random_body_at(r, DVec3::ZERO)
}

/// [`random_body`] moved by `shift`.
fn random_body_at(r: &mut Rng, shift: DVec3) -> Body {
    let snap = r.chance(0.5);
    let aligned = snap || r.chance(0.3);
    if r.chance(0.4) {
        let lo = DVec3::new(
            r.coord(-1.0, 1.0, snap),
            r.coord(-1.0, 1.0, snap),
            r.coord(-1.0, 1.0, snap),
        );
        let size = DVec3::new(
            r.coord(0.25, 1.5, snap).max(0.25),
            r.coord(0.25, 1.5, snap).max(0.25),
            r.coord(0.25, 1.5, snap).max(0.25),
        );
        let frame = random_frame(r, aligned, lo + shift);
        return k().make_box(&frame, size).expect("a box");
    }
    let centre = DVec3::new(
        r.coord(-0.5, 0.5, snap),
        r.coord(-0.5, 0.5, snap),
        r.coord(-0.5, 0.5, snap),
    );
    let frame = random_frame(r, aligned, centre + shift);
    let n = 3 + (r.next() % 6) as usize;
    let mut angles: Vec<f64> = (0..n)
        .map(|i| (i as f64 + r.range(0.3, 0.7)) * TAU / n as f64)
        .collect();
    angles.sort_by(f64::total_cmp);
    let radii: Vec<f64> = (0..n).map(|_| r.range(0.5, 1.2)).collect();
    let pts: Vec<DVec2> = angles
        .iter()
        .zip(&radii)
        .map(|(a, rad)| DVec2::new(a.cos(), a.sin()) * *rad)
        .collect();
    let mut profile = Profile::polygon(frame, &pts);
    if r.chance(0.3) {
        // A hole well inside the polygon, which contains the disc of
        // radius min(r)·cos(g/2) for its largest angle step g ≤ 1.4·2π/n.
        let inner = 0.8 * 0.5 * (1.4 * std::f64::consts::PI / n as f64).cos();
        let m = 3 + (r.next() % 3) as usize;
        let hole: Vec<DVec2> = (0..m)
            .map(|i| {
                let a = i as f64 * TAU / m as f64 + 0.3;
                DVec2::new(a.cos(), a.sin()) * inner
            })
            .collect();
        profile = profile.with(Profile::polygon(frame, &hole));
    }
    let distance = if snap {
        r.coord(0.25, 1.5, true).max(0.25)
    } else {
        r.range(0.3, 1.5)
    };
    k().extrude(
        &profile,
        &ExtrudeSpec {
            distance,
            symmetric: r.chance(0.3),
        },
    )
    .expect("an extrude")
}

#[derive(Default)]
struct Tally {
    pairs: usize,
    checked: usize,
    /// Pairs whose intersection has volume.
    overlapping: usize,
    refused: Vec<String>,
    bugs: Vec<String>,
}

/// The volume of a boolean's result, `Some(0)` for an empty result, `None`
/// for a refusal (recorded); kernel bugs are recorded too, among them a
/// result that touches itself ([`self_contact::touches_itself`]). With it,
/// the outcome in short: the result's counts, "empty", the refusal's text
/// without its place, or "bug".
fn measured(t: &mut Tally, a: &Body, b: &Body, op: BooleanOp, case: &str) -> (Option<f64>, String) {
    let bug = |t: &mut Tally, m: String| {
        t.bugs.push(format!("{case} {op:?}: {m}"));
        (None, "bug".to_string())
    };
    match k().boolean(a, b, op) {
        Ok(r) => {
            let issues = k().validate(&r);
            if !issues.is_empty() {
                return bug(t, format!("invalid result: {issues:?}"));
            }
            if let Err(e) = origins_agree(a, b, &r, op) {
                return bug(t, format!("origins: {e}"));
            }
            match self_contact::touches_itself(&k(), &r, &TessellationOptions::default()) {
                Ok(None) => {}
                Ok(Some(m)) => return bug(t, format!("the result touches itself: {m}")),
                Err(e) => return bug(t, format!("mesh: {e}")),
            }
            match k().mass_properties(&r, 1e-12) {
                Ok(m) => (Some(m.volume), format!("{:?}", counts(&r))),
                Err(e) => bug(t, e.to_string()),
            }
        }
        Err(KernelError::InvalidInput(m)) if m.contains("empty") => (Some(0.0), "empty".into()),
        Err(KernelError::InvalidInput(m)) => {
            t.refused.push(format!("{case} {op:?}: {m}"));
            let short = m.split(" near ").next().unwrap_or(&m).to_string();
            (None, short)
        }
        Err(e) => bug(t, e.to_string()),
    }
}

/// Checks the identities on one pair; `slack` is the volume a coincidence
/// within the tolerance may move.
fn identities(t: &mut Tally, a: &Body, b: &Body, case: &str, slack: f64) {
    identities_with(t, a, b, case, slack);
}

/// [`identities`], returning the outcomes of A ∪ B, A ∩ B, A − B and
/// B − A in short ([`measured`]).
fn identities_with(t: &mut Tally, a: &Body, b: &Body, case: &str, slack: f64) -> [String; 4] {
    t.pairs += 1;
    let (va, vb) = (volume(a), volume(b));
    let (u, su) = measured(t, a, b, BooleanOp::Union, case);
    let (i, si) = measured(t, a, b, BooleanOp::Intersect, case);
    let (ab, sab) = measured(t, a, b, BooleanOp::Subtract, case);
    let (ba, sba) = measured(t, b, a, BooleanOp::Subtract, case);
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
    [su, si, sab, sba]
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
fn volume_identities_on_random_pairs() {
    let mut t = Tally::default();
    let mut r = Rng(0x0DA7_0012);
    let n = scale(300);
    for case in 0..n {
        let a = random_body(&mut r);
        let b = random_body(&mut r);
        identities(&mut t, &a, &b, &format!("pair {case}"), 0.0);
    }
    report("random polyhedra", &t);
    assert!(t.checked * 10 >= t.pairs * 8, "too few pairs fully checked");
}

/// The motions' fractions of the tolerance, case by case in turn.
const FRACTIONS: [f64; 7] = [0.25, 0.5, 0.9, 1.1, 1.5, 3.0, 10.0];

/// The motions, case by case: the seed of a random body and the direction
/// its copy moves in.
fn motions() -> impl Iterator<Item = (u64, DVec3)> {
    let mut r = Rng(0x7015_E000);
    (0usize..).map(move |case| {
        let seed = r.next();
        let dir = match case % 4 {
            0 => DVec3::X,
            1 => DVec3::Y,
            2 => DVec3::Z,
            _ => DVec3::new(r.range(-1.0, 1.0), r.range(-1.0, 1.0), r.range(-1.0, 1.0))
                .normalize_or(DVec3::X),
        };
        (seed, dir)
    })
}

/// Motion `case`: a random body, its copy moved by `FRACTIONS[case % 7]`
/// tolerances, and the slack the motion allows the volume identities.
fn motion(case: usize, (seed, dir): (u64, DVec3)) -> (Body, Body, f64) {
    let eps = k().tolerance().linear;
    let shift = eps * FRACTIONS[case % 7];
    let a = random_body(&mut Rng(seed));
    let b = random_body_at(&mut Rng(seed), dir * shift);
    let slack = 4.0 * area(&a) * (shift + eps);
    (a, b, slack)
}

/// Motions by 1.1 tolerances whose difference is a sliver between a face
/// and its moved copy, the two bounded by one loop of points (two faces,
/// four edges, four vertices; found by the motions' run a hundred times
/// longer): it is a solid, its volume what the identities give.
#[test]
fn slivers_between_a_face_and_its_moved_copy() {
    for (case, first) in [(836, true), (4539, false), (5883, false), (6751, false)] {
        let (a, b, _) = motion(case, motions().nth(case).expect("a motion"));
        let (x, y) = if first { (&a, &b) } else { (&b, &a) };
        let r = (k().boolean(x, y, BooleanOp::Subtract)).unwrap_or_else(|e| panic!("{case}: {e}"));
        assert_eq!(k().validate(&r), Vec::<String>::new(), "{case}");
        assert_eq!(counts(&r), (2, 4, 4), "{case}");
        let common = k()
            .boolean(x, y, BooleanOp::Intersect)
            .expect("the common part");
        let v = volume(x) - volume(&common);
        assert!(
            v > 0.0 && (volume(&r) - v).abs() < 1e-12,
            "{case}: {} {v}",
            volume(&r)
        );
    }
}

#[test]
fn motions_by_fractions_of_the_tolerance() {
    let mut t = Tally::default();
    let mut by_shift = [0usize; 7];
    for (case, m) in motions().take(scale(70)).enumerate() {
        let (a, b, slack) = motion(case, m);
        let before = t.refused.len();
        identities(&mut t, &a, &b, &format!("shift {case}"), slack);
        by_shift[case % 7] += t.refused.len() - before;
    }
    let ops = 4 * t.pairs / 7;
    let line: Vec<String> = FRACTIONS
        .iter()
        .zip(by_shift)
        .map(|(f, n)| format!("{f} tolerance: {n}/{ops}"))
        .collect();
    println!(
        "motions by the tolerance, refusals by shift: {}",
        line.join(", ")
    );
    report("motions by the tolerance", &t);
    // A copy moved by half the tolerance or less is the same body, and
    // every operation on the two is answered (none of 8,000 refused at a
    // hundred times this run).
    assert_eq!(
        by_shift[0] + by_shift[1],
        0,
        "copies moved by at most half a tolerance refused"
    );
}

/// How a body of the contact fuzz rests on a face.
#[derive(Clone, Copy, Debug, PartialEq)]
enum Rest {
    /// A box on one of its corners.
    Corner,
    /// A triangular prism on one of its edges.
    Edge,
}

/// The body of contact case `seed`: resting by `rest` on the plane through
/// `x` with unit normal `n`, on the side `side` (+1 the side `n` points
/// to), then moved by `delta` along `n`. A corner's three edges all rise
/// from the plane; an edge runs along `along` if given (an edge of the
/// face), or a random direction in the plane, and its prism's two faces
/// rise from it at 9° or more. The same seed gives the same body at every
/// `delta`.
fn resting(
    seed: u64,
    rest: Rest,
    (x, n): (DVec3, DVec3),
    side: f64,
    delta: f64,
    along: Option<DVec3>,
) -> Body {
    let mut r = Rng(seed);
    let (up, at) = (n * side, x + n * delta);
    match rest {
        Rest::Corner => {
            let frame = loop {
                let f = random_frame(&mut r, false, at);
                if [f.x, f.y, f.z].iter().all(|a| a.dot(up) > 0.15) {
                    break f;
                }
            };
            let size = DVec3::new(r.range(0.1, 0.8), r.range(0.1, 0.8), r.range(0.1, 0.8));
            k().make_box(&frame, size).expect("a box on its corner")
        }
        Rest::Edge => {
            let v = DVec3::new(r.range(-1.0, 1.0), r.range(-1.0, 1.0), r.range(-1.0, 1.0));
            let d = along.unwrap_or_else(|| {
                v.reject_from_normalized(n)
                    .normalize_or(n.any_orthonormal_vector())
            });
            let len = r.range(0.2, 1.6);
            let start = r.range(0.0, 1.0) * len;
            let frame = Frame {
                origin: at - d * start,
                x: n.cross(d),
                y: n,
                z: d,
            };
            let a1 = r.range(-0.6, 0.5);
            let a2 = r.range(a1 + 0.1, 0.6);
            let (h1, h2) = (r.range(0.1, 0.6), r.range(0.1, 0.6));
            let mut pts = vec![
                DVec2::ZERO,
                DVec2::new(a2, side * h2),
                DVec2::new(a1, side * h1),
            ];
            if (pts[1] - pts[0]).perp_dot(pts[2] - pts[0]) < 0.0 {
                pts.swap(1, 2);
            }
            k().extrude(
                &Profile::polygon(frame, &pts),
                &ExtrudeSpec {
                    distance: len,
                    symmetric: false,
                },
            )
            .expect("a prism on its edge")
        }
    }
}

/// A point of face `f` of a body: inside one of its triangles (true), or,
/// one case in four, on one of its edges or at a corner, with the edge's
/// direction.
fn on_face(
    r: &mut Rng,
    mesh: &BodyMesh,
    adj: &Adjacency,
    f: usize,
) -> (DVec3, Option<DVec3>, bool) {
    let edges: Vec<usize> = (0..adj.edge_faces.len())
        .filter(|&e| adj.edge_faces[e].contains(&FaceIx(f as u32)))
        .collect();
    if r.chance(0.25) && !edges.is_empty() {
        let line = &mesh.edges[edges[(r.next() % edges.len() as u64) as usize]];
        let (p, q) = (line[0], line[line.len() - 1]);
        let t = if r.chance(0.3) {
            0.0
        } else {
            r.range(0.2, 0.8)
        };
        return (p.lerp(q, t), Some((q - p).normalize()), false);
    }
    let face = &mesh.faces[f];
    let tri = face.triangles[(r.next() % face.triangles.len() as u64) as usize];
    let w = [r.range(0.05, 1.0), r.range(0.05, 1.0), r.range(0.05, 1.0)];
    let sum = w[0] + w[1] + w[2];
    let x = (0..3)
        .map(|i| face.positions[tri[i] as usize] * (w[i] / sum))
        .fold(DVec3::ZERO, |a, b| a + b);
    (x, None, true)
}

/// The motions of the contact fuzz, in tolerances along the touched face's
/// normal: across the band in quarters, and four tolerances either way.
const STEPS: [f64; 11] = [
    -4.0, -1.0, -0.75, -0.5, -0.25, 0.0, 0.25, 0.5, 0.75, 1.0, 4.0,
];

/// A case of the contact fuzz: the body A, the point `x` of its face with
/// normal `n`, whether `x` is inside the face, and the body resting there
/// at each motion.
struct ContactCase {
    a: Body,
    x: DVec3,
    n: DVec3,
    inside: bool,
    rest: Rest,
    side: f64,
    along: Option<DVec3>,
    seed: u64,
}

impl ContactCase {
    fn new(r: &mut Rng) -> ContactCase {
        let a = random_body(r);
        let mesh = k()
            .tessellate(&a, &TessellationOptions::default())
            .expect("a mesh");
        let adj = k().adjacency(&a);
        let f = (r.next() % mesh.faces.len() as u64) as usize;
        let n = k()
            .face_plane(&a, FaceIx(f as u32))
            .expect("a flat face")
            .normal;
        let (x, edge, inside) = on_face(r, &mesh, &adj, f);
        let rest = if r.chance(0.5) {
            Rest::Corner
        } else {
            Rest::Edge
        };
        let side = if r.chance(0.5) { 1.0 } else { -1.0 };
        let along = edge.filter(|_| r.chance(0.5));
        let seed = r.next();
        ContactCase {
            a,
            x,
            n,
            inside,
            rest,
            side,
            along,
            seed,
        }
    }

    /// The resting body moved by `q` tolerances along the face's normal.
    fn body(&self, q: f64) -> Body {
        let eps = k().tolerance().linear;
        let (c, x) = (self, (self.x, self.n));
        resting(c.seed, c.rest, x, c.side, q * eps, c.along)
    }
}

/// Case `n` of the contact fuzz.
fn contact_case(n: usize) -> ContactCase {
    let mut r = Rng(0xC047_AC75);
    let mut c = ContactCase::new(&mut r);
    for _ in 0..n {
        c = ContactCase::new(&mut r);
    }
    c
}

/// Edges lying on a face for a few tolerances (contacts 86, 2504 and 5245,
/// found at a hundred times the default run): B rests by an edge on a
/// corner or an edge of A's face from inside, and a piece of B's edge
/// three to four tolerances long, from where it passes the corner or meets
/// A's faces to where it crosses another face of A, lies on that face.
/// A − B kept the piece with B's faces on both sides, so the result
/// touched itself; such a piece is now checked against the faces at its
/// ends, and A − B is refused. Every operation passes the identities.
#[test]
fn edges_lying_on_a_face_for_a_few_tolerances() {
    let eps = k().tolerance().linear;
    let mut t = Tally::default();
    for (case, q) in [(86, -0.5), (2504, -4.0), (5245, -0.5)] {
        let c = contact_case(case);
        let b = c.body(q);
        match k().boolean(&c.a, &b, BooleanOp::Subtract) {
            Err(KernelError::InvalidInput(m)) => assert!(
                m.starts_with("the result would be non-manifold: the bodies touch along an edge"),
                "{case}: {m}"
            ),
            other => panic!("{case}: {other:?}"),
        }
        let slack = 4.0 * (area(&c.a) + area(&b)) * (q.abs() + 1.0) * eps;
        identities_with(&mut t, &c.a, &b, &format!("contact {case}"), slack);
    }
    report("edges lying on a face for a few tolerances", &t);
}

/// The contact fuzz (ADR 0024): a random polyhedron A and a body B resting
/// on a point of one of A's faces by a corner or an edge, from outside or
/// inside, swept across the tolerance band. Every result passes the
/// identities, the origins check and the oracle for results that touch
/// themselves; every operation is answered or refused, never a kernel
/// error. Where B rests from outside inside A's face, A ∪ B pinches and is
/// refused; from inside, A − B (design note 0012 §18.1; B is convex, so
/// A's face is outside B about the contact). Moved four tolerances clear,
/// neither is refused as non-manifold for that contact.
#[test]
fn contacts_swept_across_the_band() {
    let eps = k().tolerance().linear;
    let mut t = Tally::default();
    let mut by_step = [0usize; STEPS.len()];
    // Per step, operations whose outcome differs from the contact's own.
    let mut changed = [0usize; STEPS.len()];
    let mut pinch: Vec<String> = Vec::new();
    let mut r = Rng(0xC047_AC75);
    let cases = scale(60);
    let zero = STEPS.iter().position(|&q| q == 0.0).expect("a step of 0");
    for case in 0..cases {
        let c = ContactCase::new(&mut r);
        let (a, rest, side) = (&c.a, c.rest, c.side);
        let label = |q: f64| format!("contact {case} ({rest:?}, side {side}, {q} tolerance)");
        let slack = |b: &Body, q: f64| 4.0 * (area(a) + area(b)) * (q.abs() + 1.0) * eps;
        let got: Vec<[String; 4]> = (STEPS.iter().enumerate())
            .map(|(j, &q)| {
                let b = c.body(q);
                let before = t.refused.len();
                let got = identities_with(&mut t, a, &b, &label(q), slack(&b, q));
                by_step[j] += t.refused.len() - before;
                got
            })
            .collect();
        for (j, g) in got.iter().enumerate() {
            changed[j] += g.iter().zip(&got[zero]).filter(|(x, y)| x != y).count();
        }
        if !c.inside {
            continue;
        }
        // Where B rests inside A's face, the operation the table refuses.
        let (k_op, name) = if side > 0.0 {
            (0, "A ∪ B")
        } else {
            (2, "A − B")
        };
        if !got[zero][k_op].starts_with("the result would be non-manifold") {
            pinch.push(format!("{}: {name} gave {}", label(0.0), got[zero][k_op]));
        }
        let op = [BooleanOp::Union, BooleanOp::Intersect, BooleanOp::Subtract][k_op];
        if let Err(KernelError::InvalidInput(m)) = k().boolean(a, &c.body(4.0 * side), op) {
            if m.starts_with("the result would be non-manifold") {
                pinch.push(format!("{} moved clear: {name}: {m}", label(0.0)));
            }
        }
    }
    let ops = 4 * cases;
    let line: Vec<String> = STEPS
        .iter()
        .zip(by_step.iter().zip(changed))
        .map(|(q, (n, c))| format!("{q}: {n} refused, {c} changed"))
        .collect();
    println!(
        "contacts by step, of {ops} operations each: {}",
        line.join("; ")
    );
    for p in &pinch {
        println!("  pinch: {p}");
    }
    report("contacts swept across the band", &t);
    assert!(pinch.is_empty(), "{} pinches", pinch.len());
}

/// Debugging aid (temporary): one case of the contact fuzz, `CONTACT` and
/// `STEP` from the environment.
#[test]
fn zz_contact_case() {
    let Some(n) = std::env::var("CONTACT")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
    else {
        return;
    };
    let q: f64 = std::env::var("STEP")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(0.0);
    let c = contact_case(n);
    let b = c.body(q);
    println!(
        "x {:?} n {:?} inside {} rest {:?} side {} along {:?}",
        c.x, c.n, c.inside, c.rest, c.side, c.along
    );
    println!(
        "A {:?} vol {} B {:?} vol {}",
        counts(&c.a),
        volume(&c.a),
        counts(&b),
        volume(&b)
    );
    for (x, y, op) in [
        (&c.a, &b, BooleanOp::Union),
        (&c.a, &b, BooleanOp::Intersect),
        (&c.a, &b, BooleanOp::Subtract),
        (&b, &c.a, BooleanOp::Subtract),
    ] {
        match k().boolean(x, y, op) {
            Ok(res) => {
                if let Err(e) = origins_agree(x, y, &res, op) {
                    println!("{op:?} origins: {e}");
                    let m = k()
                        .tessellate(&res, &TessellationOptions::default())
                        .unwrap();
                    let (mx, my) = (
                        k().tessellate(x, &TessellationOptions::default()).unwrap(),
                        k().tessellate(y, &TessellationOptions::default()).unwrap(),
                    );
                    let og = k().origins(&res);
                    for (e, l) in m.edges.iter().enumerate() {
                        let src = match og.edges[e] {
                            EdgeOrigin::Operand {
                                operand: Operand::First,
                                edge,
                            } => Some(&mx.edges[edge.0 as usize]),
                            EdgeOrigin::Operand {
                                operand: Operand::Second,
                                edge,
                            } => Some(&my.edges[edge.0 as usize]),
                            _ => None,
                        };
                        println!(
                            "  edge {e} {:?} {:?} -> {:?} src {:?}",
                            og.edges[e],
                            l[0],
                            l[l.len() - 1],
                            src.map(|q| (q[0], q[q.len() - 1]))
                        );
                    }
                    for mm in &og.merged {
                        println!("  merged {mm:?}");
                    }
                }
                let o = self_contact::touches_itself(&k(), &res, &TessellationOptions::default());
                println!(
                    "{op:?}: {:?} vol {} oracle {:?}",
                    counts(&res),
                    volume(&res),
                    o
                );
                if o.as_ref().is_ok_and(|x| x.is_some()) && std::env::var_os("DUMP").is_some() {
                    let m = k()
                        .tessellate(&res, &TessellationOptions::default())
                        .unwrap();
                    let adj = k().adjacency(&res);
                    let og = k().origins(&res);
                    for (f, fo) in og.faces.iter().enumerate() {
                        println!(
                            "  face {f} {:?} {fo:?}",
                            k().face_plane(&res, FaceIx(f as u32))
                        );
                    }
                    for (e, l) in m.edges.iter().enumerate() {
                        println!(
                            "  edge {e} {:?} faces {:?} {:?} -> {:?}",
                            og.edges[e],
                            adj.edge_faces[e],
                            l[0],
                            l[l.len() - 1]
                        );
                    }
                    for (v, p) in m.vertices.iter().enumerate() {
                        println!(
                            "  vertex {v} {:?} {:?} edges {:?}",
                            og.vertices[v], p, adj.vertex_edges[v]
                        );
                    }
                }
            }
            Err(e) => println!("{op:?}: {e:?}"),
        }
    }
}

/// Copies of a body with a hole moved by three tolerances (motions 3449
/// and 6781, found at a hundred times the default run): where the holes'
/// edges are 2.4 tolerances apart, the strip of one body's face between
/// them is off the other's face but within the tolerance of a corner of
/// it, and was taken for a piece on that face, making an inconsistent
/// result. Every operation must be answered or refused, never
/// inconsistent (today all four are refused: the strip is within the
/// tolerance of the other body's surface).
#[test]
fn a_strip_between_moved_holes() {
    for case in [3449, 6781] {
        let m = motions().nth(case).expect("a motion");
        let (a, b, slack) = motion(case, m);
        let mut t = Tally::default();
        identities(&mut t, &a, &b, &format!("shift {case}"), slack);
        report(&format!("motion {case}"), &t);
    }
}

/// Two sides of a prism meeting at a crease of 5e-6 rad, which a test on
/// the cosine of the angle between their normals (1 − cos θ = 1.25e-11)
/// took for one plane, so that a boolean anywhere on the body merged them
/// into one face with a corner 5 µm off its plane. They stay two faces,
/// and their edge is a crease, not a seam.
#[test]
fn a_slight_crease_stays_a_crease() {
    let turn = 5e-6;
    let pts = [
        DVec2::new(0.0, 0.0),
        DVec2::new(1.0, 0.0),
        DVec2::new(2.0, turn),
        DVec2::new(2.0, 1.0),
        DVec2::new(0.0, 1.0),
    ];
    let a = k()
        .extrude(
            &Profile::polygon(Frame::WORLD, &pts),
            &ExtrudeSpec {
                distance: 1.0,
                symmetric: false,
            },
        )
        .expect("a prism");
    assert_eq!(counts(&a), (7, 15, 10));
    let mesh = k()
        .tessellate(&a, &TessellationOptions::default())
        .expect("a mesh");
    assert!(mesh.seams.iter().all(|s| !s), "the crease is not a seam");
    // A box through the far side, away from the crease.
    let b = cuboid(DVec3::new(1.5, 0.25, 0.25), DVec3::new(2.5, 0.75, 0.75));
    let u = run(&a, &b, BooleanOp::Union);
    assert_eq!(counts(&u), (12, 27, 18));
    let expect = (2.0 - 0.5 * turn) + 0.25 - 0.125;
    assert!(close(volume(&u), expect), "{} vs {expect}", volume(&u));
}
