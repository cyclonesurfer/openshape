//! Bodies that touch without crossing (design note 0012 §18, ADR 0024):
//! an edge or a corner of one body resting on a face of the other, from
//! outside or inside, on its own or with crossings in the same face.
//!
//! Each golden states the outcome of all six operations (both unions,
//! both intersections, both differences) with exact texts and places for
//! refusals, and counts (faces, edges, vertices, solids), volume and area
//! for results; every result must validate, have origins that agree with
//! the operands, and not touch itself (the oracle in `self_contact`).
//! Where a golden checks that a contact which does not make the result
//! non-manifold leaves no trace, the result is the operand it equals,
//! entity for entity with that operand's origins (`is_operand`), or the
//! same as with the touching body moved clear of the other
//! (`as_if_apart`).

use opendatum_kernel::{
    Body, BooleanOp, BrepKernel, CurveKind, DVec2, DVec3, EdgeOrigin, ExtrudeSpec, FaceIx,
    FaceOrigin, Frame, Kernel, KernelError, Merged, Operand, Origins, Profile, TessellationOptions,
    Tolerance, VertexOrigin,
};

mod self_contact;

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

fn prism(frame: Frame, pts: &[DVec2], h: f64) -> Body {
    k().extrude(
        &Profile::polygon(frame, pts),
        &ExtrudeSpec {
            distance: h,
            symmetric: false,
        },
    )
    .expect("a prism")
}

fn op(a: &Body, b: &Body, o: BooleanOp) -> Body {
    k().boolean(a, b, o)
        .unwrap_or_else(|e| panic!("building an operand: {o:?}: {e}"))
}

/// A profile's frame in the xz plane at `origin`, extruded toward −y.
fn xz(origin: DVec3) -> Frame {
    Frame {
        origin,
        x: DVec3::X,
        y: DVec3::Z,
        z: -DVec3::Y,
    }
}

/// A profile's frame in the yz plane at `origin`, extruded toward +x.
fn yz(origin: DVec3) -> Frame {
    Frame {
        origin,
        x: DVec3::Y,
        y: DVec3::Z,
        z: DVec3::X,
    }
}

/// The triangle of three points, scaled by `s`.
fn tri(s: f64, p: [(f64, f64); 3]) -> [DVec2; 3] {
    p.map(|(x, y)| DVec2::new(x, y) * s)
}

/// A cube of side `side` with a corner at `at` and its three edges there
/// rising at the same angle (`down`: falling).
fn corner_cube(at: DVec3, side: f64, down: bool) -> Body {
    let s3 = 3f64.sqrt();
    let x = DVec3::new((2.0f64 / 3.0).sqrt(), 0.0, 1.0 / s3);
    let y = DVec3::new(-1.0 / 6f64.sqrt(), 1.0 / 2f64.sqrt(), 1.0 / s3);
    let z = DVec3::new(-1.0 / 6f64.sqrt(), -1.0 / 2f64.sqrt(), 1.0 / s3);
    let frame = if down {
        Frame {
            origin: at,
            x: -x,
            y: -z,
            z: -y,
        }
    } else {
        Frame {
            origin: at,
            x,
            y,
            z,
        }
    };
    k().make_box(&frame, DVec3::splat(side)).expect("a cube")
}

/// The plate P: `cuboid((−2, −2, −1), (2, 2, 0))`, scaled by `s`.
fn plate(s: f64) -> Body {
    cuboid(
        DVec3::new(-2.0, -2.0, -1.0) * s,
        DVec3::new(2.0, 2.0, 0.0) * s,
    )
}

/// The wedge W, scaled by `s` and raised by `z`: its edge x = 0.3, z = 0,
/// |y| ≤ 0.5, its faces rising to z = 0.5 at x = −0.2 and x = 0.8.
fn wedge(s: f64, z: f64) -> Body {
    prism(
        xz(DVec3::new(0.0, 0.5 * s, z)),
        &tri(s, [(0.3, 0.0), (0.8, 0.5), (-0.2, 0.5)]),
        s,
    )
}

/// The cube C of side 0.5 on its corner at (0.3, 0.2, 0), scaled by `s`
/// and raised by `z`.
fn cube_on_corner(s: f64, z: f64) -> Body {
    corner_cube(DVec3::new(0.3 * s, 0.2 * s, z), 0.5 * s, false)
}

/// H: the plate with a square hole through (−1.2, −1.2)–(−0.8, −0.8).
fn holed_plate(s: f64) -> Body {
    let hole = cuboid(
        DVec3::new(-1.2, -1.2, -2.0) * s,
        DVec3::new(-0.8, -0.8, 1.0) * s,
    );
    op(&plate(s), &hole, BooleanOp::Subtract)
}

/// K: a knife whose edge y = −1, z = 0 runs from x = −2.5 beyond the
/// plate's edge to x = −1, over the hole of [`holed_plate`].
fn knife(s: f64, z: f64) -> Body {
    prism(
        yz(DVec3::new(-2.5 * s, 0.0, z)),
        &tri(s, [(-1.0, 0.0), (-0.8, 0.3), (-1.2, 0.3)]),
        1.5 * s,
    )
}

fn volume(b: &Body) -> f64 {
    k().mass_properties(b, 1e-12)
        .expect("mass properties")
        .volume
}

fn area(b: &Body) -> f64 {
    k().mass_properties(b, 1e-12).expect("mass properties").area
}

fn counts(b: &Body) -> (usize, usize, usize, usize) {
    let t = k().topology(b);
    (t.faces, t.edges, t.vertices, t.solids)
}

/// The outcome of `a op b`: a refusal's text, "empty", or the result's
/// counts, volume over `s`³ and area over `s`²; a result must validate and
/// have origins that agree with the operands.
fn outcome(a: &Body, b: &Body, o: BooleanOp, s: f64) -> String {
    match k().boolean(a, b, o) {
        Ok(r) => {
            assert_eq!(k().validate(&r), Vec::<String>::new(), "{o:?}");
            origins_agree(a, b, &r, o).unwrap_or_else(|e| panic!("{o:?}: {e}"));
            alone(&r, o);
            format!(
                "{:?} V={:.9} A={:.9}",
                counts(&r),
                volume(&r) / (s * s * s),
                area(&r) / (s * s)
            )
        }
        Err(KernelError::InvalidInput(m)) if m == "the result would be empty" => "empty".into(),
        Err(KernelError::InvalidInput(m)) => format!("refused: {m}"),
        Err(e) => panic!("{o:?}: {e:?}"),
    }
}

/// The outcomes of A ∪ B, B ∪ A, A ∩ B, B ∩ A, A − B and B − A.
fn six(a: &Body, b: &Body, s: f64) -> [String; 6] {
    use BooleanOp::*;
    [
        outcome(a, b, Union, s),
        outcome(b, a, Union, s),
        outcome(a, b, Intersect, s),
        outcome(b, a, Intersect, s),
        outcome(a, b, Subtract, s),
        outcome(b, a, Subtract, s),
    ]
}

fn edge_refusal(at: &str) -> String {
    format!("refused: the result would be non-manifold: the bodies touch along an edge near {at}")
}

fn point_refusal(at: &str) -> String {
    format!("refused: the result would be non-manifold: the bodies touch at a point near {at}")
}

fn band_refusal(at: &str) -> String {
    format!(
        "refused: a face of one body runs within the tolerance of the other's surface without \
         lying on one of its faces, near {at}"
    )
}

fn run(a: &Body, b: &Body, o: BooleanOp) -> Body {
    let r = k()
        .boolean(a, b, o)
        .unwrap_or_else(|e| panic!("{o:?}: {e}"));
    assert_eq!(k().validate(&r), Vec::<String>::new(), "{o:?}");
    origins_agree(a, b, &r, o).unwrap_or_else(|e| panic!("{o:?}: {e}"));
    alone(&r, o);
    r
}

/// Asserts that a result does not touch itself
/// ([`self_contact::touches_itself`]).
fn alone(r: &Body, o: BooleanOp) {
    let opts = TessellationOptions {
        chord_tolerance: 1e-4,
        max_angle: 0.5,
    };
    let found =
        self_contact::touches_itself(&k(), r, &opts).unwrap_or_else(|e| panic!("{o:?}: {e}"));
    assert_eq!(found, None, "{o:?}: the result touches itself");
}

/// `a o b` leaves no trace of the contact: it is operand `side` (`a` for
/// the first), with each of its faces, edges and vertices named once by
/// its origin, as if the other were apart.
fn is_operand(a: &Body, b: &Body, o: BooleanOp, side: Operand) {
    let r = run(a, b, o);
    let x = match side {
        Operand::First => a,
        Operand::Second => b,
    };
    assert_eq!(counts(&r), counts(x), "{o:?}");
    assert!(
        close(volume(&r), volume(x)) && close(area(&r), area(x)),
        "{o:?}"
    );
    let got = k().origins(&r);
    let (f, e, v, _) = counts(x);
    let mut faces: Vec<u32> = got
        .faces
        .iter()
        .map(|g| match *g {
            FaceOrigin::Operand { operand, face } if operand == side => face.0,
            ref other => panic!("{o:?}: a face from {other:?}"),
        })
        .collect();
    let mut edges: Vec<u32> = got
        .edges
        .iter()
        .map(|g| match *g {
            EdgeOrigin::Operand { operand, edge } if operand == side => edge.0,
            ref other => panic!("{o:?}: an edge from {other:?}"),
        })
        .collect();
    let mut vertices: Vec<u32> = got
        .vertices
        .iter()
        .map(|g| match *g {
            VertexOrigin::Operand { operand, vertex } if operand == side => vertex.0,
            ref other => panic!("{o:?}: a vertex from {other:?}"),
        })
        .collect();
    faces.sort_unstable();
    edges.sort_unstable();
    vertices.sort_unstable();
    assert_eq!(faces, (0..f as u32).collect::<Vec<_>>(), "{o:?}");
    assert_eq!(edges, (0..e as u32).collect::<Vec<_>>(), "{o:?}");
    assert_eq!(vertices, (0..v as u32).collect::<Vec<_>>(), "{o:?}");
    assert_eq!(got.merged, Vec::<Merged>::new(), "{o:?}");
}

/// The origins of a result as sorted lists, entity numbers aside.
fn origin_lists(o: &Origins) -> (Vec<String>, Vec<String>, Vec<String>, Vec<String>) {
    let list = |v: Vec<String>| {
        let mut v = v;
        v.sort();
        v
    };
    (
        list(o.faces.iter().map(|x| format!("{x:?}")).collect()),
        list(o.edges.iter().map(|x| format!("{x:?}")).collect()),
        list(o.vertices.iter().map(|x| format!("{x:?}")).collect()),
        list(
            o.merged
                .iter()
                .map(|m| match m {
                    Merged::Face(_, x) => format!("face {x:?}"),
                    Merged::Edge(_, x) => format!("edge {x:?}"),
                    Merged::Vertex(_, x) => format!("vertex {x:?}"),
                })
                .collect(),
        ),
    )
}

/// `r` and `q` are the same result: counts, volume, area and origins.
fn same(r: &Body, q: &Body, what: &str) {
    assert_eq!(counts(r), counts(q), "{what}");
    assert!(
        close(volume(r), volume(q)) && close(area(r), area(q)),
        "{what}"
    );
    assert_eq!(
        origin_lists(k().origins(r)),
        origin_lists(k().origins(q)),
        "{what}"
    );
}

/// `a o b` leaves no trace of the contact: it is the same as `apart o b`
/// (`apart` being `a` moved clear), in counts, volume, area and origins.
fn as_if_apart(a: &Body, apart: &Body, b: &Body, o: BooleanOp) {
    same(&run(a, b, o), &run(apart, b, o), &format!("{o:?}"));
}

/// As [`as_if_apart`], the touching body being the second operand.
fn as_if_apart_second(a: &Body, b: &Body, apart: &Body, o: BooleanOp) {
    same(&run(a, b, o), &run(a, apart, o), &format!("{o:?}"));
}

fn close(x: f64, y: f64) -> bool {
    (x - y).abs() <= 1e-12 * (1.0 + x.abs().max(y.abs()))
}

/// Checks a result's origins against the operands' surfaces, as the
/// booleans' tests do: a face on the surface of the operand face it comes
/// from, facing the same way (the other way for the second operand of a
/// difference); an edge on the surfaces of the faces along the operand edge
/// it is a piece of, and along that edge if it is straight, or on both
/// faces' surfaces where the bodies meet; a vertex at an operand vertex;
/// within three tolerances.
fn origins_agree(a: &Body, b: &Body, r: &Body, o: BooleanOp) -> Result<(), String> {
    let near = 3.0 * k().tolerance().linear;
    // Mesh points a few hundredths of the unit apart at most on curved
    // faces; straight edges are meshed by their ends.
    let opts = TessellationOptions {
        chord_tolerance: 1e-3,
        max_angle: 0.5,
    };
    let mesh = |x: &Body| k().tessellate(x, &opts).map_err(|e| e.to_string());
    let (ma, mb, mr) = (mesh(a)?, mesh(b)?, mesh(r)?);
    let (ja, jb) = (k().adjacency(a), k().adjacency(b));
    let (ka, kb) = (k().edge_kinds(a), k().edge_kinds(b));
    let body = |x: Operand| match x {
        Operand::First => (a, &ma, &ja, &ka),
        Operand::Second => (b, &mb, &jb, &kb),
    };
    // How far `p` is from the surface of face `f` of `x`, and the outward
    // normal there (no input here has a surface without a tangent plane).
    let on = |x: &Body, f: FaceIx, p: DVec3| -> Result<(f64, DVec3), String> {
        let t = k()
            .tangent_plane(x, f, p)
            .ok_or_else(|| format!("face {} has no tangent plane near {p}", f.0))?;
        Ok(((p - t.point).length(), t.normal))
    };
    let origins = k().origins(r);
    let joined: Vec<usize> = origins
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
        let (x, _, _, _) = body(operand);
        let flip = o == BooleanOp::Subtract && operand == Operand::Second;
        for p in &mr.faces[i].positions {
            let ((d, n), (_, m)) = (on(x, face, *p)?, on(r, FaceIx(i as u32), *p)?);
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
        Ok(())
    };
    let edge = |i: usize, origin: &EdgeOrigin| -> Result<(), String> {
        let p = &mr.edges[i];
        let faces: Vec<(&Body, FaceIx)> = match *origin {
            EdgeOrigin::Operand { operand, edge } => {
                let (x, m, adj, kinds) = body(operand);
                let e = edge.0 as usize;
                if kinds[e] == CurveKind::Line && !along(i, p, &m.edges[e]) {
                    return Err(format!("edge {i} is not along {operand:?} edge {e}"));
                }
                adj.edge_faces[e].iter().map(|&f| (x, f)).collect()
            }
            EdgeOrigin::Meet { first, second } => vec![(a, first), (b, second)],
            _ => return Err(format!("edge {i}: origin {origin:?}")),
        };
        for x in p {
            for &(y, f) in &faces {
                let (d, _) = on(y, f, *x)?;
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
                let (_, m, _, _) = body(operand);
                if (mr.vertices[i] - m.vertices[vertex.0 as usize]).length() > near {
                    return Err(format!("vertex {i} is off {operand:?} vertex {}", vertex.0));
                }
                Ok(())
            }
            VertexOrigin::Meet => Ok(()),
            _ => Err(format!("vertex {i}: origin {origin:?}")),
        }
    };
    for (i, f) in origins.faces.iter().enumerate() {
        face(i, f)?;
    }
    for (i, e) in origins.edges.iter().enumerate() {
        edge(i, e)?;
    }
    for (i, v) in origins.vertices.iter().enumerate() {
        vertex(i, v)?;
    }
    for m in &origins.merged {
        match m {
            Merged::Face(f, origin) => face(f.0 as usize, origin)?,
            Merged::Edge(e, origin) => edge(e.0 as usize, origin)?,
            Merged::Vertex(v, origin) => vertex(v.0 as usize, origin)?,
        }
    }
    Ok(())
}

// ---------------------------------------------------------------- probes

const RW: &str = "a face of one body runs within the tolerance of the other's surface";

/// A placement: a rotation by `phi` radians about z, then a shift.
#[derive(Clone, Copy)]
struct W {
    phi: f64,
    t: DVec3,
}

impl W {
    fn id() -> W {
        W {
            phi: 0.0,
            t: DVec3::ZERO,
        }
    }
    fn rot(&self, v: DVec3) -> DVec3 {
        let (s, c) = self.phi.sin_cos();
        DVec3::new(c * v.x - s * v.y, s * v.x + c * v.y, v.z)
    }
    fn frame(&self, f: Frame) -> Frame {
        Frame {
            origin: self.rot(f.origin) + self.t,
            x: self.rot(f.x),
            y: self.rot(f.y),
            z: self.rot(f.z),
        }
    }
    fn cuboid(&self, lo: DVec3, hi: DVec3) -> Body {
        k().make_box(
            &self.frame(Frame {
                origin: lo,
                ..Frame::WORLD
            }),
            hi - lo,
        )
        .expect("a box")
    }
    fn prism(&self, frame: Frame, pts: &[DVec2], h: f64) -> Body {
        prism(self.frame(frame), pts, h)
    }
    /// [`tilted`] placed.
    fn tilted(&self, t: f64, pts: &[(f64, f64, bool)]) -> Body {
        let eps = k().tolerance().linear;
        let t = t.to_radians();
        let frame = Frame {
            origin: DVec3::ZERO,
            x: DVec3::X,
            y: DVec3::new(0.0, t.cos(), t.sin()),
            z: DVec3::new(0.0, -t.sin(), t.cos()),
        };
        let p: Vec<DVec2> = (pts.iter())
            .map(|&(u, w, tol)| DVec2::new(u, if tol { w * eps } else { w }))
            .collect();
        self.prism(frame, &p, 0.5)
    }
    /// A block from `lo` to `hi` in x, y, from z = −1 to 0, with an xy
    /// profile instead (a prism of height 1 under z = 0).
    fn block(&self, pts: &[(f64, f64)]) -> Body {
        let p: Vec<DVec2> = pts.iter().map(|&(x, y)| DVec2::new(x, y)).collect();
        self.prism(
            Frame {
                origin: DVec3::new(0.0, 0.0, -1.0),
                ..Frame::WORLD
            },
            &p,
            1.0,
        )
    }
    fn beside_t(&self, t: f64, w: f64) -> Body {
        self.tilted(
            t,
            &[
                (-0.5, -0.5, false),
                (0.0, -w, true),
                (1.0, -w, true),
                (1.5, -0.5, false),
                (2.0, -0.5, false),
                (2.0, 0.5, false),
                (-0.5, 0.5, false),
            ],
        )
    }
}

fn rep(a: &Body, b: &Body, o: BooleanOp) -> String {
    let r = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        match k().boolean(a, b, o) {
            Ok(r) => {
                let v = k().validate(&r);
                let opts = TessellationOptions {
                    chord_tolerance: 1e-4,
                    max_angle: 0.5,
                };
                let oracle = self_contact::touches_itself(&k(), &r, &opts);
                let org = origins_agree(a, b, &r, o);
                format!(
                    "{:?} V={:.12} A={:.9} valid={} origins={} oracle={:?}",
                    counts(&r),
                    volume(&r),
                    area(&r),
                    if v.is_empty() {
                        "ok".to_string()
                    } else {
                        format!("{v:?}")
                    },
                    if org.is_ok() {
                        "ok".to_string()
                    } else {
                        "BAD".to_string()
                    },
                    oracle
                )
            }
            Err(KernelError::InvalidInput(m)) if m == "the result would be empty" => "empty".into(),
            Err(e) => format!("refused: {e}"),
        }
    }));
    r.unwrap_or_else(|_| "PANICKED".into())
}

fn show(name: &str, a: &Body, b: &Body) {
    use BooleanOp::*;
    eprintln!("=== {name}");
    for (tag, o, x, y) in [
        ("A+B", Union, a, b),
        ("B+A", Union, b, a),
        ("A*B", Intersect, a, b),
        ("B*A", Intersect, b, a),
        ("A-B", Subtract, a, b),
        ("B-A", Subtract, b, a),
    ] {
        eprintln!("  {tag} {}", rep(x, y, o));
    }
}

fn holed_prism(w: &W, t: f64, x0: f64, x1: f64, w0: f64, w1: f64) -> Body {
    let full = w.tilted(
        t,
        &[
            (-0.5, -0.5, false),
            (2.0, -0.5, false),
            (2.0, 0.5, false),
            (-0.5, 0.5, false),
        ],
    );
    let tr = t.to_radians();
    let n = DVec3::new(0.0, -tr.sin(), tr.cos());
    let frame = Frame {
        origin: -n * 0.3,
        x: DVec3::X,
        y: DVec3::new(0.0, tr.cos(), tr.sin()),
        z: n,
    };
    let cutter = w.prism(
        frame,
        &[
            DVec2::new(x0, w0),
            DVec2::new(x1, w0),
            DVec2::new(x1, w1),
            DVec2::new(x0, w1),
        ],
        1.1,
    );
    op(&full, &cutter, BooleanOp::Subtract)
}

#[test]
fn zz_overcheck() {
    let eps = k().tolerance().linear;
    let w = W::id();
    let narrow = w.cuboid(DVec3::new(0.2, -1.0, -1.0), DVec3::new(0.8, 1.0, 0.0));
    let wide = w.cuboid(DVec3::new(-1.25, -1.0, -1.0), DVec3::new(1.25, 1.0, 0.0));
    let far = w.cuboid(DVec3::new(-0.25, -1.0, -1.0), DVec3::new(1.9, 1.0, 0.0));
    let boxed = w.cuboid(DVec3::new(-0.25, -1.0, -1.0), DVec3::new(1.25, 1.0, 0.0));
    show(
        "ov58 hole -5e..0.3",
        &narrow,
        &holed_prism(&w, 10.0, 0.0, 1.0, -5.0 * eps, 0.3),
    );
    show(
        "ov58 hole -0.3..+5e",
        &narrow,
        &holed_prism(&w, 10.0, 0.0, 1.0, -0.3, 5.0 * eps),
    );
    show(
        "ov58 hole -0.3..0.01",
        &narrow,
        &holed_prism(&w, 10.0, 0.0, 1.0, -0.3, 0.01),
    );
    show(
        "ov58 t3 hole -15e..0.3",
        &narrow,
        &holed_prism(&w, 3.0, 0.0, 1.0, -15.0 * eps, 0.3),
    );
    show(
        "ov58 t3 hole -5e..0.3",
        &narrow,
        &holed_prism(&w, 3.0, 0.0, 1.0, -5.0 * eps, 0.3),
    );
    let u = w.tilted(
        10.0,
        &[
            (-0.5, -0.5, false),
            (2.0, -0.5, false),
            (2.0, 0.5, false),
            (1.0, 0.5, false),
            (1.0, -5.0, true),
            (0.0, -5.0, true),
            (0.0, 0.5, false),
            (-0.5, 0.5, false),
        ],
    );
    show("ov58 U-notch floor -5e", &narrow, &u);
    let tab = w.tilted(
        10.0,
        &[
            (-0.5, -0.5, false),
            (2.0, -0.5, false),
            (2.0, 0.5, false),
            (1.5, 0.5, false),
            (1.5, -5.0, true),
            (-0.5, -5.0, true),
        ],
    );
    show("ov58 tab top -5e", &narrow, &tab);
    let above = w.tilted(
        10.0,
        &[
            (-0.5, 5.0, true),
            (2.0, 5.0, true),
            (2.0, 0.5, false),
            (-0.5, 0.5, false),
        ],
    );
    show("ov58 F above, bottom +5e", &narrow, &above);
    let tab2 = w.tilted(
        10.0,
        &[
            (-0.5, 5.0, true),
            (1.5, 5.0, true),
            (1.5, -0.5, false),
            (2.0, -0.5, false),
            (2.0, 0.5, false),
            (-0.5, 0.5, false),
        ],
    );
    show("ov58 F above with tab, bottom +5e", &narrow, &tab2);
    let step = w.tilted(
        10.0,
        &[
            (-1.0, -0.5, false),
            (0.0, -0.5, false),
            (0.0, 2.0, true),
            (1.0, 2.0, true),
            (1.5, 0.5, false),
            (-1.0, 0.5, false),
        ],
    );
    show("zza_step2", &wide, &step);
    let fc = w.tilted(
        10.0,
        &[
            (0.0, 0.5, false),
            (0.0, 3.0, true),
            (0.5, 3.0, true),
            (1.5, -0.01, false),
            (2.0, -0.01, false),
            (2.0, 0.5, false),
        ],
    );
    show("farcross", &far, &fc);
    let m8 = w.tilted(
        10.0,
        &[
            (-0.5, 0.15, false),
            (0.0, 0.15, false),
            (0.5, -3.0, true),
            (0.6, -0.5, false),
            (2.0, -0.5, false),
            (2.0, 0.5, false),
            (-0.5, 0.5, false),
        ],
    );
    show("m8 reflex corner crossing", &boxed, &m8);
    let m8b = w.tilted(
        10.0,
        &[
            (-0.5, -0.5, false),
            (0.4, -0.5, false),
            (0.5, -3.0, true),
            (1.0, 0.15, false),
            (2.0, 0.15, false),
            (2.0, 0.5, false),
            (-0.5, 0.5, false),
        ],
    );
    show("m8 reflex corner crossing, mirrored", &boxed, &m8b);
}

fn zig(n: usize) -> Body {
    let w = W::id();
    let mut pts = vec![(-0.5, -0.5, false), (0.0, -0.5, false)];
    for i in 0..=n {
        let u = i as f64 / n as f64;
        pts.push((u, if i % 2 == 0 { -4.0 } else { -2.0 }, true));
    }
    pts.extend([
        (1.0, -0.5, false),
        (2.0, -0.5, false),
        (2.0, 0.5, false),
        (-0.5, 0.5, false),
    ]);
    w.tilted(10.0, &pts)
}

fn notch_block() -> Body {
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
    prism(frame, &pts, 1.0)
}

#[test]
fn zz_zigtime() {
    let a = notch_block();
    let box_ = W::id().cuboid(DVec3::new(-0.25, -1.0, -1.0), DVec3::new(1.25, 1.0, 0.0));
    for n in [16usize, 64, 256, 1024] {
        let b = zig(n);
        for (name, x) in [("notch", &a), ("box", &box_)] {
            let t0 = std::time::Instant::now();
            let r = rep(x, &b, BooleanOp::Union);
            eprintln!(
                "ZIG n={n} {name}: {:.2} ms {}",
                t0.elapsed().as_secs_f64() * 1e3,
                &r[..r.len().min(160)]
            );
        }
    }
}

#[test]
fn zz_convex() {
    let w = W::id();
    let boxed = w.cuboid(DVec3::new(-0.25, -1.0, -1.0), DVec3::new(1.25, 1.0, 0.0));
    let tri = w.tilted(
        10.0,
        &[(0.0, 0.5, false), (0.5, -5.0, true), (1.0, 0.5, false)],
    );
    show("convex tri tip -5", &boxed, &tri);
    let wide = w.tilted(
        10.0,
        &[
            (0.5, -5.0, true),
            (1.5, 0.2, false),
            (2.0, 0.5, false),
            (-0.5, 0.5, false),
            (-0.5, 0.2, false),
        ],
    );
    show("convex wide tip -5", &boxed, &wide);
    let wide3 = w.tilted(
        10.0,
        &[
            (0.5, -3.0, true),
            (1.5, 0.2, false),
            (2.0, 0.5, false),
            (-0.5, 0.5, false),
            (-0.5, 0.2, false),
        ],
    );
    show("convex wide tip -3", &boxed, &wide3);
}

#[test]
fn zz_zigops() {
    let a = notch_block();
    for n in [256usize, 1024] {
        let b = zig(n);
        for (name, x, y, o) in [
            ("A*B", &a, &b, BooleanOp::Intersect),
            ("A-B", &a, &b, BooleanOp::Subtract),
            ("B-A", &b, &a, BooleanOp::Subtract),
        ] {
            let t0 = std::time::Instant::now();
            let r = rep(x, y, o);
            eprintln!(
                "ZOP n={n} {name}: {:.2} ms {}",
                t0.elapsed().as_secs_f64() * 1e3,
                &r[..r.len().min(160)]
            );
        }
    }
}
