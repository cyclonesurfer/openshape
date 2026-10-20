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

/// A prism on the plane through the x axis leaning 10° from z = 0, its
/// face's edge over x from 0 to 1 zigzagging in `n` segments between 1 mm
/// to either side of the axis, against a block under z = 0 from x = −0.25
/// to 1.25: the face crosses the block's top along the axis `n` times, so
/// the faces' regions have cycles of more sides than are walked one by
/// one (placed through trees, `planar::regions`). A ∩ B is the prism's two
/// ends under the top, each with half a tooth, and the `n / 2 − 1` teeth
/// between them, of volume tan 10° times the integral of the depth below
/// the axis over the profile: 2 × 0.25 × 0.5² / 2 for the ends beyond x =
/// 0 and 1, and w² / 12 for the teeth and the two half teeth (`n / 2`
/// triangles' worth, each of base 1 / n and depth w = 1 mm, w² / 6n).
#[test]
fn a_face_crossing_another_many_times_along_one_line() {
    let n = 64;
    let w = 1e-3;
    let t = 10f64.to_radians();
    let frame = Frame {
        origin: DVec3::ZERO,
        x: DVec3::X,
        y: DVec3::new(0.0, t.cos(), t.sin()),
        z: DVec3::new(0.0, -t.sin(), t.cos()),
    };
    let mut pts = vec![DVec2::new(-0.5, -0.5), DVec2::new(0.0, -0.5)];
    pts.extend((0..=n).map(|i| DVec2::new(i as f64 / n as f64, if i % 2 == 0 { -w } else { w })));
    pts.extend([(1.0, -0.5), (2.0, -0.5), (2.0, 0.5), (-0.5, 0.5)].map(|(x, y)| DVec2::new(x, y)));
    let comb = k()
        .extrude(
            &Profile::polygon(frame, &pts),
            &ExtrudeSpec {
                distance: 0.5,
                symmetric: false,
            },
        )
        .expect("a comb");
    let block = cuboid(DVec3::new(-0.25, -1.0, -1.0), DVec3::new(1.25, 1.0, 0.0));
    let i = run(&block, &comb, BooleanOp::Intersect);
    assert_eq!(k().topology(&i).solids, n / 2 + 1);
    let expected = t.tan() * (0.0625 + w * w / 12.0);
    assert!(
        (volume(&i) - expected).abs() <= 1e-9 * expected,
        "{} {expected}",
        volume(&i)
    );
    let (va, vb) = (volume(&block), volume(&comb));
    let u = run(&block, &comb, BooleanOp::Union);
    assert!((volume(&u) - (va + vb - expected)).abs() <= 1e-9 * (va + vb));
    for (x, y, v) in [(&block, &comb, va), (&comb, &block, vb)] {
        let d = run(x, y, BooleanOp::Subtract);
        assert!((volume(&d) - (v - expected)).abs() <= 1e-9 * v);
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
/// 2.98 to 3.60 tolerances long, from where it passes the corner or meets
/// A's faces to where it crosses another face of A, lies on that face.
/// A − B kept the piece with B's faces on both sides, so the result
/// touched itself; the middle of such a piece is now checked against a
/// face with corners at both its ends, and A − B is refused; the other
/// three are answered and pass the identities and the checks of
/// [`measured`].
#[test]
fn edges_lying_on_a_face_for_a_few_tolerances() {
    let along = "the result would be non-manifold: the bodies touch along an edge";
    pinned_contacts(&[
        (
            86,
            -0.5,
            ["(20, 45, 29)", "(7, 15, 10)", along, "(13, 21, 14)"],
        ),
        (
            2504,
            -4.0,
            ["(15, 34, 22)", "(7, 15, 10)", along, "(9, 15, 10)"],
        ),
        (
            5245,
            -0.5,
            ["(19, 42, 27)", "(8, 18, 12)", along, "(14, 24, 16)"],
        ),
    ]);
}

/// Contacts at the boundary of a face, which the checks of design note 0012
/// §21 leave to the settling of touched points and the corners within
/// their reach (found by the contact fuzz at a hundred times the default
/// run). B rests by a corner (contacts 548 and 2237 at 1 tolerance) or an
/// edge (86 at 1) on A's face, with a corner of a face of B across the line
/// where that face's plane meets A's face's, both its edges crossing the
/// line (a flap), and the clustering puts the corner on the boundary of A's
/// face. Main answered these, A − B touching along an edge or at a point;
/// the flap's strip now leaves out only its points within the tolerance of
/// a member of the corner's cluster, the same in either order of the
/// operands (§21, rule E), A's face holds the flap beyond those, and each is
/// refused in all four operations: the cost of that rule (with discs of
/// three tolerances about the cluster's place, main's outcomes). Or B
/// rests by an edge passing a corner of A's face (4917 at −0.5 tolerance
/// and 2383 at −0.75), the edge lying on the face's plane as the whole
/// interval of B's face on the line, which no interval of A's face meets,
/// so that B's face lies within the tolerance of A's over a sliver at the
/// corner, A's flap there, which its own rule decides and the run's strip
/// leaves out (§21, rule D; tested as a part of the run's strip, each was
/// refused in all four operations). Their outcomes are main's, and the
/// answered results pass the identities, the checks of [`measured`] and
/// the faces-crossing check.
#[test]
fn contacts_at_the_boundary_of_a_face_across_the_line_where_the_faces_cross() {
    let point = "the result would be non-manifold: the bodies touch at a point";
    let band = "a face of one body runs within the tolerance of the other's surface without lying \
                on one of its faces,";
    let cases = [
        (548, 1.0, [band; 4]),
        (2237, 1.0, [band; 4]),
        (86, 1.0, [band; 4]),
        (4917, -0.5, [point, "empty", "(6, 12, 8)", "(5, 9, 6)"]),
        (2383, -0.75, [point, "empty", "(6, 12, 8)", "(5, 9, 6)"]),
    ];
    pinned_contacts(&cases);
    for (case, q, _) in cases {
        let c = contact_case(case);
        let b = c.body(q);
        for (x, y, op) in [
            (&c.a, &b, BooleanOp::Union),
            (&c.a, &b, BooleanOp::Intersect),
            (&c.a, &b, BooleanOp::Subtract),
            (&b, &c.a, BooleanOp::Subtract),
        ] {
            if let Ok(r) = k().boolean(x, y, op) {
                let crossing = self_contact::faces_cross(&k(), &r, 1e-8).expect("a mesh");
                assert_eq!(crossing, None, "contact {case} ({q}) {op:?}");
            }
        }
    }
}

/// Contacts at the boundary of a face held beyond the reach of its places
/// (design note 0012 §21), found by the contact fuzz at a hundred times the
/// default run. A flap held 3.2 tolerances from its corner's place
/// (contact 3129 at 0.5 tolerance), where main returned A − B with the
/// flap's face and A's crossing 1.36 tolerances inside both, as the first
/// build did, which left the whole flap out; and a flap held 3.3
/// tolerances from its corner's place (5385 at 1), where main refused
/// three operations as a point within the tolerance of an edge and B − A
/// with the strip's text: each is refused with the strip's text in all
/// four. And B resting by an edge that passes a corner of A's face as the
/// whole interval of its face on the line, 3.1 and 3.9 tolerances from the
/// corner (3613 and 3492 at −0.75) and, the edge 0.55° from the face's
/// edge, 54 and 105 tolerances along it (4917 at −0.75 and −1): what A's
/// face holds of the run's strip there is A's flap at the corner, which
/// its own rule decides and the run's strip leaves out (§21, rule D), and
/// the outcomes are main's, A ∪ B touching at a point (tested as a part of
/// the run's strip, the flap refused every operation).
#[test]
fn contacts_at_the_boundary_of_a_face_held_beyond_the_reach_of_its_places() {
    let band = "a face of one body runs within the tolerance of the other's surface without lying \
                on one of its faces,";
    let point = "the result would be non-manifold: the bodies touch at a point";
    pinned_contacts(&[
        (3129, 0.5, [band; 4]),
        (5385, 1.0, [band; 4]),
        (3613, -0.75, [point, "empty", "(5, 9, 6)", "(5, 9, 6)"]),
        (3492, -0.75, [point, "empty", "(7, 15, 10)", "(5, 9, 6)"]),
        (4917, -0.75, [point, "empty", "(6, 12, 8)", "(5, 9, 6)"]),
        (4917, -1.0, [point, "empty", "(6, 12, 8)", "(5, 9, 6)"]),
    ]);
}

/// Sharp corners that the check of a result's places lets off (design note
/// 0012 §19, `check_apart`), each found by disabling its rules over longer
/// runs of the contact fuzz: a piece of an edge shorter than four
/// tolerances between two corners of a face, its middle outside the face
/// but within the tolerance of its boundary (contacts 481 at 0.5 tolerance
/// and 2896 at 0.75); one shorter than four tolerances leaving a face's
/// corner (84 at 1, and golden 26); and a longer one leaving a corner whose
/// far end is more than the tolerance off the face's plane (22 at −1, 281
/// at 1, 1335 and 5358 at −4, crossing the face's plane at 14° and 3.6°;
/// and 2714 at −1, at the edge of the band, an edge a tolerance under the
/// face whose far end the imprint calls off it, 1.0000000002 tolerances,
/// and whose other end it places on the face's corner). With its rule
/// disabled, each is refused as touching along an edge; the results are
/// valid and pass the oracle, which reported the last three until it let
/// off a flat face that a straight edge leaves at a corner.
#[test]
fn sharp_corners_the_check_of_places_lets_off() {
    let along = "the result would be non-manifold: the bodies touch along an edge";
    pinned_contacts(&[
        (
            84,
            1.0,
            [
                "(15, 38, 25)",
                "(8, 18, 12)",
                "(12, 28, 18)",
                "(11, 27, 18)",
            ],
        ),
        (
            481,
            0.5,
            ["(18, 42, 27)", "(8, 18, 12)", along, "(15, 28, 19)"],
        ),
        (
            2896,
            0.75,
            [
                "(18, 49, 32)",
                "(8, 18, 12)",
                "(16, 42, 27)",
                "(10, 24, 16)",
            ],
        ),
        (
            22,
            -1.0,
            ["(14, 31, 20)", "(5, 9, 6)", "(9, 22, 14)", "(10, 18, 12)"],
        ),
        (
            281,
            1.0,
            [
                "(23, 64, 42)",
                "(14, 36, 24)",
                "(22, 61, 40)",
                "(15, 39, 26)",
            ],
        ),
        (
            1335,
            -4.0,
            ["(21, 54, 35)", "(4, 6, 4)", "(16, 42, 28)", "(9, 18, 11)"],
        ),
        (
            2714,
            -1.0,
            ["(12, 27, 17)", "(4, 6, 4)", "(10, 21, 13)", "(6, 12, 8)"],
        ),
        (
            5358,
            -4.0,
            ["(15, 39, 26)", "(5, 9, 6)", "(12, 30, 20)", "(8, 18, 12)"],
        ),
    ]);
}

/// Cases of the contact fuzz at a hundred times its default run where a
/// corner of a result is beyond the tolerance of a face it is not on, off
/// the face's boundary: 1.08 tolerances (contact 3827 at −4) and 1.31
/// (369 at −1), so the result keeps clear of the face. The check of a
/// result's places first took a point within the tolerance of a face's
/// plane and, in the plane, of its loops for on the face, up to √2
/// tolerances from it off a corner, and refused these; it now takes the
/// two distances together. The results are valid and pass the oracle.
#[test]
fn corners_just_beyond_the_tolerance_of_a_face() {
    let along = "the result would be non-manifold: the bodies touch along an edge";
    pinned_contacts(&[
        (
            3827,
            -4.0,
            ["(20, 54, 36)", "(6, 12, 8)", "(18, 48, 32)", "(8, 18, 12)"],
        ),
        (
            369,
            -1.0,
            ["(14, 36, 24)", "(6, 12, 8)", along, "(6, 12, 8)"],
        ),
    ]);
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

/// Cases of the contact fuzz at a hundred times its length where a bridge
/// of a face carried no edge of the other body, a kernel error before:
/// the curves where the face meets two faces of the other body merged
/// into one over a sliver of the face (contact 4996), or at a corner of
/// the other body within the tolerance of the face (2152, 5387); an edge
/// lying on the face and leaving it at a slant, resting on it from
/// outside (1203 at 0.75, the union pinching) or with its faces crossing
/// the face beside it (1203 at −0.75 and 661, refused as a point the
/// tolerance cannot place); the tree of a corner's bridges (161, refused
/// by the bridge's test of a piece of the other body's edge leaving its
/// end); 477, whose bridge carries the curves of the two faces of an edge
/// lying on the face, refused where the edge leaves it (B − A as edges
/// passing within the tolerance of each other); and 149, whose merged
/// curves run along an edge of the other body lying on the face, refused
/// where the faces meet. 149, 661 and 1203 at −0.75 are now
/// refused by the check of the faces' meeting's ends
/// (`imprint::plane_meeting`) before a bridge forms, with the same text;
/// 161 and 5387 by the check of a flap of a face across the line where it
/// crosses the other's (design note 0012 §21: each a flap whose layer lies
/// in its own body and outside the other's, its touch unpinched, which its
/// own rule now tests, 161's held beyond the tolerance of the members of
/// its corner's cluster), where the face runs within the tolerance of the
/// other's surface, in every operation.
#[test]
fn contacts_whose_bridges_carry_no_edge() {
    let leaves = "an edge of one body lies on a face of the other and leaves it at a point the \
                  tolerance cannot place";
    let sliver = "the bodies leave a sliver thinner than the tolerance";
    let band = "a face of one body runs within the tolerance of the other's surface without lying \
                on one of its faces,";
    pinned_contacts(&[
        (149, -1.0, [leaves; 4]),
        (4996, -4.0, [sliver; 4]),
        (661, 0.75, [leaves; 4]),
        (
            1203,
            0.75,
            [
                "the result would be non-manifold: the bodies touch along an edge",
                "empty",
                "(6, 12, 8)",
                "(5, 9, 6)",
            ],
        ),
        (1203, -0.75, [leaves; 4]),
        (161, 1.0, [band; 4]),
        (
            477,
            -1.0,
            [
                leaves,
                leaves,
                leaves,
                "edges where the bodies meet pass within the tolerance of each other",
            ],
        ),
        (2152, -1.0, [sliver, sliver, sliver, leaves]),
        (5387, -4.0, [band; 4]),
    ]);
}

/// Cases of the contact fuzz at a hundred times its length where the
/// faces at an edge lying on a face left it to both sides, or a contact's
/// pieces disagreed on its side, kernel errors before: a box's corner by
/// an edge of the other body that crosses one of the box's faces at a
/// slant a few tolerances from the corner (681, 1583, 1666, 4183, 4906),
/// refused as an edge running within the tolerance of a face over a
/// stretch (1583's B − A as a sliver); a prism's face within the tolerance of the other body's over
/// a band (2080), refused there in every operation, where one order of the
/// operands reached the edge first; a prism's edge parallel to an edge of
/// the face and about the tolerance from it, that edge lying on the
/// prism's faces (828, as golden 24); and an edge lying on a face that
/// leaves the body through a neighbour at an acute angle before it reaches
/// their edge (3003, as golden 25).
#[test]
fn contacts_whose_lying_edges_leave_the_face_to_both_sides() {
    let stretch = "an edge of one body runs within the tolerance of a face of the other over a \
                   stretch without lying on it,";
    let band = "a face of one body runs within the tolerance of the other's surface without lying \
                on one of its faces,";
    let leaves = "an edge of one body lies on a face of the other and leaves it at a point the \
                  tolerance cannot place";
    let sliver = "the bodies leave a sliver thinner than the tolerance";
    pinned_contacts(&[
        (681, -1.0, [stretch; 4]),
        (1583, 1.0, [stretch, stretch, stretch, sliver]),
        (1666, -1.0, [stretch; 4]),
        (4183, -0.75, [stretch; 4]),
        (4906, -1.0, [stretch; 4]),
        (2080, 4.0, [band; 4]),
        (828, -1.0, [leaves; 4]),
        (3003, -1.0, [leaves; 4]),
    ]);
}

/// Cases of the contact fuzz where a box's corner, resting inside a face
/// of A, is pushed four tolerances through it, a metre or less from the
/// world's origin, so that one result holds the box's tip as a tetrahedron
/// a few tolerances across: A ∩ B where the box rests from outside
/// (contacts 2 and 42), and B − A where it rests from inside (1, beside a
/// larger piece of the box, and 15). Taken from the world's origin, the
/// volume of that shell took the wrong sign, in the booleans (whose faces,
/// projected onto their planes, fail to close by up to the tolerance) and
/// in the validator (by rounding), and the operation failed with an
/// invalid body, a void inside no solid (1, 2) or an outer shell turned
/// inside out (15, 42); both now take their volumes from a corner of the
/// shell.
#[test]
fn corners_dipping_four_tolerances_through_a_face() {
    pinned_contacts(&[
        (
            1,
            4.0,
            [
                "(15, 33, 22)",
                "(7, 15, 10)",
                "(11, 27, 18)",
                "(11, 21, 14)",
            ],
        ),
        (
            2,
            -4.0,
            ["(14, 33, 22)", "(4, 6, 4)", "(11, 24, 16)", "(7, 15, 10)"],
        ),
        (
            15,
            4.0,
            ["(9, 18, 12)", "(7, 15, 10)", "(12, 27, 18)", "(4, 6, 4)"],
        ),
        (
            42,
            -4.0,
            ["(12, 27, 18)", "(4, 6, 4)", "(9, 18, 12)", "(7, 15, 10)"],
        ),
    ]);
}

/// Asserts the outcomes of A ∪ B, A ∩ B, A − B and B − A for cases of the
/// contact fuzz, each with its motion, and that every operation passes
/// the identities and the checks of [`measured`].
fn pinned_contacts(cases: &[(usize, f64, [&str; 4])]) {
    let eps = k().tolerance().linear;
    let mut t = Tally::default();
    for &(case, q, want) in cases {
        let c = contact_case(case);
        let b = c.body(q);
        let label = format!("contact {case} ({q} tolerance)");
        let slack = 4.0 * (area(&c.a) + area(&b)) * (q.abs() + 1.0) * eps;
        let got = identities_with(&mut t, &c.a, &b, &label, slack);
        assert_eq!(got, want, "{label}");
    }
    assert_eq!(t.bugs, Vec::<String>::new());
}

/// A prism 0.5 thick whose face lies on the plane through the x axis
/// tilted `t` degrees toward +y, with the profile `pts` (u along x, w up
/// the plane); `below`: the prism lies below its face.
fn tilted_prism(t: f64, pts: &[DVec2], below: bool) -> Body {
    let t = t.to_radians();
    let (y, z) = (
        DVec3::new(0.0, t.cos(), t.sin()),
        DVec3::new(0.0, -t.sin(), t.cos()),
    );
    let (frame, pts) = if below {
        let frame = Frame {
            origin: DVec3::ZERO,
            x: -DVec3::X,
            y,
            z: -z,
        };
        (
            frame,
            pts.iter().rev().map(|p| DVec2::new(-p.x, p.y)).collect(),
        )
    } else {
        let frame = Frame {
            origin: DVec3::ZERO,
            x: DVec3::X,
            y,
            z,
        };
        (frame, pts.to_vec())
    };
    k().extrude(
        &Profile::polygon(frame, &pts),
        &ExtrudeSpec {
            distance: 0.5,
            symmetric: false,
        },
    )
    .expect("a tilted prism")
}

/// The distance from `p` to the segment from `a` to `b`.
fn segment_distance(p: DVec2, a: DVec2, b: DVec2) -> f64 {
    let d = b - a;
    let t = ((p - a).dot(d) / d.length_squared()).clamp(0.0, 1.0);
    (a + d * t - p).length()
}

/// A random case of [`notches_and_spikes_across_the_line`]: a profile on
/// the tilted plane and its angle, with a label. Its corners, one or two,
/// lie within a few tolerances of the plane z = 0, across the axis or short
/// of it: the tip of a notch from the face's far side or near side, the
/// tip of a spike rising to the axis or falling to it, or two notches'
/// tips (a W) a few tolerances to a few hundred apart. Every vertex of the
/// profile is more than 10.5 tolerances from every edge it is not on (the
/// profile's own merge distance).
fn notch_case(r: &mut Rng, eps: f64) -> (Vec<DVec2>, f64, String) {
    let log = |r: &mut Rng, lo: f64, hi: f64| lo * (hi / lo).powf(r.unit());
    loop {
        let t = if r.chance(0.5) {
            [0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 45.0][(r.next() % 7) as usize]
        } else {
            log(r, 0.5, 45.0)
        };
        let st = t.to_radians().sin();
        // How far a tip lies across the plane z = 0, in tolerances, and so
        // across the axis in the face's plane.
        let depth = |r: &mut Rng| {
            if r.chance(0.2) {
                [-0.5, 0.0, 0.5, 1.0, 1.5][(r.next() % 5) as usize]
            } else {
                r.range(-1.0, 2.0)
            }
        };
        let (z, z2) = (depth(r), depth(r));
        let (d, d2) = (z * eps / st, z2 * eps / st);
        let x = match r.next() % 4 {
            0 => 1.25 + r.range(-10.0, 10.0) * eps,
            1 => -0.25 + r.range(-10.0, 10.0) * eps,
            _ => r.range(0.0, 1.0),
        };
        let p = DVec2::new;
        let kind = r.next() % 5;
        let (a, b) = (log(r, 0.005, 0.4), log(r, 0.005, 0.4));
        let (label, pts) = match kind {
            0 => (
                format!("far notch at {x:.9} {a:.4}/{b:.4} depth {z:+.3}"),
                vec![
                    p(-1.2, -0.5),
                    p(2.2, -0.5),
                    p(2.2, 0.5),
                    p(x + a, 0.5),
                    p(x, -d),
                    p(x - b, 0.5),
                    p(-1.2, 0.5),
                ],
            ),
            1 => (
                format!("near notch at {x:.9} {a:.4}/{b:.4} depth {z:+.3}"),
                vec![
                    p(-1.2, -0.5),
                    p(x - b, -0.5),
                    p(x, d),
                    p(x + a, -0.5),
                    p(2.2, -0.5),
                    p(2.2, 0.5),
                    p(-1.2, 0.5),
                ],
            ),
            2 => (
                format!("spike up at {x:.9} {a:.4}/{b:.4} depth {z:+.3}"),
                vec![
                    p(-1.2, -0.5),
                    p(2.2, -0.5),
                    p(2.2, -0.4),
                    p(x + a, -0.4),
                    p(x, d),
                    p(x - b, -0.4),
                    p(-1.2, -0.4),
                ],
            ),
            3 => (
                format!("spike down at {x:.9} {a:.4}/{b:.4} depth {z:+.3}"),
                vec![
                    p(-1.2, 0.4),
                    p(x - b, 0.4),
                    p(x, -d),
                    p(x + a, 0.4),
                    p(2.2, 0.4),
                    p(2.2, 0.5),
                    p(-1.2, 0.5),
                ],
            ),
            _ => {
                let (g1, g2) = (log(r, 0.3, 150.0) * eps, log(r, 0.3, 150.0) * eps);
                let m = match r.next() % 3 {
                    0 => 0.5,
                    1 => log(r, 1e-4, 0.1),
                    _ => r.range(-0.5, 5.0) * eps / st,
                };
                (
                    format!(
                        "W at {x:.9} tips {:.3}/{:.3} tolerances apart, tooth {m:.3e}, mouths \
                         {a:.4}, depths {z:+.3}/{z2:+.3}",
                        g1 / eps,
                        g2 / eps
                    ),
                    vec![
                        p(-1.2, -0.5),
                        p(2.2, -0.5),
                        p(2.2, 0.5),
                        p(x + g2 + a, 0.5),
                        p(x + g2, -d2),
                        p(x, m),
                        p(x - g1, -d),
                        p(x - g1 - a, 0.5),
                        p(-1.2, 0.5),
                    ],
                )
            }
        };
        // Each vertex more than the profile's merge distance from every
        // edge it is not on (the profile refuses curves touching).
        let n = pts.len();
        let apart = (0..n).all(|i| {
            (0..n)
                .filter(|&j| i != j && i != (j + 1) % n)
                .all(|j| segment_distance(pts[i], pts[j], pts[(j + 1) % n]) > 10.5 * eps)
        });
        if apart {
            return (pts, t, label);
        }
    }
}

/// Random corners of a tilted prism's face within a few tolerances of a
/// box's top z = 0 (x from −0.25 to 1.25), across the line where the face
/// crosses it or short of it, beside the box's ends or over its top, the
/// prism above or below its face ([`notch_case`]; design note 0012 §21):
/// every operation is answered or refused, never a kernel error, and the
/// results pass the identities, the origins check and the oracle. Before
/// §21 the notch tips more than a tolerance across the line failed with
/// the box first, and spikes there returned faces crossing each other.
#[test]
fn notches_and_spikes_across_the_line() {
    let eps = k().tolerance().linear;
    let a = cuboid(DVec3::new(-0.25, -1.0, -1.0), DVec3::new(1.25, 1.0, 0.0));
    let mut t = Tally::default();
    let mut r = Rng(0x0DA7_0021);
    for case in 0..scale(150) {
        let (pts, angle, label) = notch_case(&mut r, eps);
        let below = r.chance(0.3);
        let b = tilted_prism(angle, &pts, below);
        let label = format!(
            "notch {case} ({label} at {angle:.3}°{})",
            if below { ", below its face" } else { "" }
        );
        identities(&mut t, &a, &b, &label, 4.0 * (area(&a) + area(&b)) * eps);
    }
    report("notches and spikes across the line", &t);
}

/// A notch tip of the notch fuzz ([`notches_and_spikes_across_the_line`],
/// case 1376, found at ten times its default run) whose foot on the box's
/// top lies within ρ = ε (1 + 2⁻¹⁰) of the top's edge: a notch from the
/// face's near side at 45°, its tip a tolerance across the line and 0.63
/// tolerances inside the box's end x = −0.25, the prism below its face.
/// Rule 3 of design note 0012 §21 tests V's foot as the end of V's edges,
/// at ρ from the top's edges as the strip search tests its points, and
/// passes it; the split refuses it as main did, as a point within the
/// tolerance of an edge it is not on. Testing the foot to a rounding margin
/// of the edges refused it at the tip, with the strip's text.
#[test]
fn a_notch_tip_whose_foot_lies_within_rho_of_the_tops_edge() {
    let eps = k().tolerance().linear;
    let a = cuboid(DVec3::new(-0.25, -1.0, -1.0), DVec3::new(1.25, 1.0, 0.0));
    let edge = "a point where the bodies meet lies within the tolerance of an edge it is not on,";
    let mut r = Rng(0x0DA7_0021);
    let mut t = Tally::default();
    for case in 0..=1376 {
        let (pts, angle, label) = notch_case(&mut r, eps);
        let below = r.chance(0.3);
        if case == 1376 {
            assert!(
                below && label.starts_with("near notch at -0.249999937"),
                "{label}"
            );
            let b = tilted_prism(angle, &pts, below);
            let slack = 4.0 * (area(&a) + area(&b)) * eps;
            assert_eq!(identities_with(&mut t, &a, &b, &label, slack), [edge; 4]);
        }
    }
    assert_eq!(t.bugs, Vec::<String>::new());
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
