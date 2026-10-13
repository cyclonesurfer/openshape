//! Booleans of solids bounded by planes (design note 0012 §8 step 1, ADR
//! 0020): golden cases with exact volumes and topology, the volume
//! identities on random pairs, and motions by fractions of the tolerance.

use std::f64::consts::TAU;

use opendatum_kernel::{
    Body, BooleanOp, BrepKernel, DVec2, DVec3, EdgeOrigin, ExtrudeSpec, FaceIx, FaceOrigin, Frame,
    Kernel, KernelError, Merged, Operand, Profile, RevolveSpec, TessellationOptions, VertexOrigin,
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

#[test]
fn bodies_touching_along_an_edge_or_at_a_point() {
    let a = cuboid(DVec3::ZERO, DVec3::ONE);
    let edge = cuboid(DVec3::new(1.0, 1.0, 0.0), DVec3::new(2.0, 2.0, 1.0));
    assert!(refused(&a, &edge, BooleanOp::Union).contains("non-manifold"));
    assert!(refused(&a, &edge, BooleanOp::Intersect).contains("empty"));
    assert!(close(volume(&run(&a, &edge, BooleanOp::Subtract)), 1.0));
    let corner = cuboid(DVec3::ONE, DVec3::splat(2.0));
    assert!(refused(&a, &corner, BooleanOp::Union).contains("non-manifold"));
    assert!(close(volume(&run(&a, &corner, BooleanOp::Subtract)), 1.0));
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
fn planes_across_a_cone_off_its_axis_are_not_supported_yet() {
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
    // Planes parallel to the axis (a hyperbola), across it at a slant (an
    // ellipse) and through the apex (two lines).
    let cuts = [
        k().make_box(
            &Frame {
                origin: DVec3::new(0.1, -1.0, -1.0),
                ..Frame::WORLD
            },
            DVec3::new(1.0, 3.0, 2.0),
        ),
        k().make_box(
            &Frame::from_normal(
                DVec3::new(-1.0, 0.5, -1.0),
                DVec3::new(0.3, 1.0, 0.0),
                DVec3::Z,
            ),
            DVec3::new(2.0, 2.0, 2.0),
        ),
        k().make_box(
            &Frame {
                origin: DVec3::new(0.0, -1.0, -1.0),
                ..Frame::WORLD
            },
            DVec3::new(1.0, 3.0, 2.0),
        ),
    ];
    for cut in cuts {
        match k().boolean(&cone, &cut.unwrap(), BooleanOp::Union) {
            Err(KernelError::Unsupported(m)) => assert!(m.contains("conic"), "{m}"),
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
/// for a refusal (recorded); kernel bugs are recorded too.
fn measured(t: &mut Tally, a: &Body, b: &Body, op: BooleanOp, case: &str) -> Option<f64> {
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
                Ok(m) => Some(m.volume),
                Err(e) => {
                    t.bugs.push(format!("{case} {op:?}: {e}"));
                    None
                }
            }
        }
        Err(KernelError::InvalidInput(m)) if m.contains("empty") => Some(0.0),
        Err(KernelError::InvalidInput(m)) => {
            t.refused.push(format!("{case} {op:?}: {m}"));
            None
        }
        Err(e) => {
            t.bugs.push(format!("{case} {op:?}: {e}"));
            None
        }
    }
}

/// Checks the identities on one pair; `slack` is the volume a coincidence
/// within the tolerance may move.
fn identities(t: &mut Tally, a: &Body, b: &Body, case: &str, slack: f64) {
    t.pairs += 1;
    let (va, vb) = (volume(a), volume(b));
    let u = measured(t, a, b, BooleanOp::Union, case);
    let i = measured(t, a, b, BooleanOp::Intersect, case);
    let ab = measured(t, a, b, BooleanOp::Subtract, case);
    let ba = measured(t, b, a, BooleanOp::Subtract, case);
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
