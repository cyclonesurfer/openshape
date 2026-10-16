//! Bodies that touch without crossing (design note 0012 §18, ADR 0024):
//! an edge or a corner of one body resting on a face of the other, from
//! outside or inside, on its own or with crossings in the same face.
//!
//! Each golden states the outcome of all six operations (both unions,
//! both intersections, both differences) with exact texts and places for
//! refusals, and counts (faces, edges, vertices, solids), volume and area
//! for results. A contact that does not make the result non-manifold
//! leaves no trace: the result is the operand it equals, entity for entity
//! with that operand's origins, or the same as with the touching body
//! moved clear of the other.

use opendatum_kernel::{
    Body, BooleanOp, BrepKernel, CurveKind, DVec2, DVec3, EdgeOrigin, ExtrudeSpec, FaceIx,
    FaceOrigin, Frame, Kernel, KernelError, Merged, Operand, Origins, Profile, TessellationOptions,
    VertexOrigin,
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
        origin_lists(&k().origins(r)),
        origin_lists(&k().origins(q)),
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

const P: &str = "(6, 12, 8, 1) V=16.000000000 A=48.000000000";
const W: &str = "(5, 9, 6, 1) V=0.250000000 A=2.914213562";
const C: &str = "(6, 12, 8, 1) V=0.125000000 A=1.500000000";

/// Golden 1: W resting on P by its edge (from outside). The union pinches
/// along the edge; the differences are the operands, as if apart (on main
/// the union returned two solids touching along the edge, and validated).
#[test]
fn a_wedge_resting_on_a_plate_by_its_edge() {
    let (p, w) = (plate(1.0), wedge(1.0, 0.0));
    let along = edge_refusal("(300.000000, 0.000000, 0.000000) mm");
    assert_eq!(
        six(&p, &w, 1.0),
        [along.as_str(), along.as_str(), "empty", "empty", P, W]
    );
    is_operand(&p, &w, BooleanOp::Subtract, Operand::First);
    is_operand(&w, &p, BooleanOp::Subtract, Operand::First);
}

/// Golden 2: C resting on P by its corner. The union pinches at the
/// corner; the differences are the operands.
#[test]
fn a_cube_resting_on_a_plate_by_its_corner() {
    let (p, c) = (plate(1.0), cube_on_corner(1.0, 0.0));
    let at = point_refusal("(300.000000, 200.000000, 0.000000) mm");
    assert_eq!(
        six(&p, &c, 1.0),
        [at.as_str(), at.as_str(), "empty", "empty", P, C]
    );
    is_operand(&p, &c, BooleanOp::Subtract, Operand::First);
    is_operand(&c, &p, BooleanOp::Subtract, Operand::First);
}

/// Golden 3: goldens 1 and 2 on P with a round hole elsewhere, which makes
/// its top face not flat: the answers do not change (on main they were
/// refused as unsupported tangent contacts).
#[test]
fn an_unrelated_round_hole_changes_nothing() {
    let rod = k()
        .extrude(
            &Profile::circle(
                Frame {
                    origin: DVec3::new(0.0, 0.0, -2.0),
                    ..Frame::WORLD
                },
                DVec2::new(-1.5, -1.5),
                0.2,
            ),
            &ExtrudeSpec {
                distance: 3.0,
                symmetric: false,
            },
        )
        .expect("a rod");
    let p = op(&plate(1.0), &rod, BooleanOp::Subtract);
    let holed = "(8, 18, 12, 1) V=15.874336294 A=49.005309649";
    let (w, c) = (wedge(1.0, 0.0), cube_on_corner(1.0, 0.0));
    let along = edge_refusal("(300.000000, 0.000000, 0.000000) mm");
    assert_eq!(
        six(&p, &w, 1.0),
        [along.as_str(), along.as_str(), "empty", "empty", holed, W]
    );
    let at = point_refusal("(300.000000, 200.000000, 0.000000) mm");
    assert_eq!(
        six(&p, &c, 1.0),
        [at.as_str(), at.as_str(), "empty", "empty", holed, C]
    );
    is_operand(&p, &w, BooleanOp::Subtract, Operand::First);
    is_operand(&p, &c, BooleanOp::Subtract, Operand::First);
}

/// Golden 4: W's edge through the centre of P's top face, where the point
/// classifying the face's piece would fall without the contact as an
/// obstacle (main refused every operation there).
#[test]
fn a_contact_through_the_middle_of_a_face() {
    let p = plate(1.0);
    let w = prism(
        xz(DVec3::new(0.0, 0.5, 0.0)),
        &tri(1.0, [(0.0, 0.0), (0.5, 0.5), (-0.5, 0.5)]),
        1.0,
    );
    let along = edge_refusal("(0.000000, 0.000000, 0.000000) mm");
    assert_eq!(
        six(&p, &w, 1.0),
        [along.as_str(), along.as_str(), "empty", "empty", P, W]
    );
    is_operand(&p, &w, BooleanOp::Subtract, Operand::First);
}

/// Golden 5: W inverted inside P, its edge on P's top face from inside
/// (nested). The union is P, the intersection W, and P − W pinches (on
/// main it returned a body whose cavity touched the top face along the
/// edge, and validated).
#[test]
fn a_wedge_inside_a_plate_touching_its_face() {
    let p = plate(1.0);
    let w = prism(
        xz(DVec3::new(0.0, 0.5, 0.0)),
        &tri(1.0, [(0.3, 0.0), (0.8, -0.5), (-0.2, -0.5)]),
        1.0,
    );
    let along = edge_refusal("(300.000000, 0.000000, 0.000000) mm");
    assert_eq!(six(&p, &w, 1.0), [P, P, W, W, along.as_str(), "empty"]);
    is_operand(&p, &w, BooleanOp::Union, Operand::First);
    is_operand(&w, &p, BooleanOp::Union, Operand::Second);
    is_operand(&p, &w, BooleanOp::Intersect, Operand::Second);
    is_operand(&w, &p, BooleanOp::Intersect, Operand::First);
}

/// Golden 6: a cube nested at its corner against the top face of a box
/// (main returned the difference touching itself at the corner).
#[test]
fn a_cube_inside_a_box_touching_its_face_by_a_corner() {
    let q = cuboid(DVec3::new(-2.0, -2.0, -2.0), DVec3::new(2.0, 2.0, 0.0));
    let c = corner_cube(DVec3::new(0.3, 0.2, 0.0), 0.5, true);
    let box_q = "(6, 12, 8, 1) V=32.000000000 A=64.000000000";
    let at = point_refusal("(300.000000, 200.000000, 0.000000) mm");
    assert_eq!(six(&q, &c, 1.0), [box_q, box_q, C, C, at.as_str(), "empty"]);
    is_operand(&q, &c, BooleanOp::Union, Operand::First);
    is_operand(&q, &c, BooleanOp::Intersect, Operand::Second);
}

/// Golden 7: the root of a groove in B resting on the top face of A, A
/// inside B. Only B − A pinches.
#[test]
fn a_grooves_root_on_a_face() {
    let a = cuboid(DVec3::new(-0.5, -0.5, -1.0), DVec3::new(0.5, 0.5, 0.0));
    let groove = prism(
        xz(DVec3::new(0.0, 2.5, 0.0)),
        &tri(1.0, [(-0.5, 1.0), (0.5, 1.0), (0.0, 0.0)]),
        5.0,
    );
    let b = op(
        &cuboid(DVec3::new(-2.0, -2.0, -2.0), DVec3::new(2.0, 2.0, 1.0)),
        &groove,
        BooleanOp::Subtract,
    );
    let grooved = "(9, 21, 14, 1) V=46.000000000 A=83.944271910";
    let cube = "(6, 12, 8, 1) V=1.000000000 A=6.000000000";
    let along = edge_refusal("(0.000000, 0.000000, 0.000000) mm");
    assert_eq!(
        six(&a, &b, 1.0),
        [grooved, grooved, cube, cube, "empty", along.as_str()]
    );
    is_operand(&a, &b, BooleanOp::Union, Operand::Second);
    is_operand(&a, &b, BooleanOp::Intersect, Operand::First);
}

/// Golden 8: the root of a notch in N lying in P's top face, N crossing
/// P: the root is a chord of P's face, and the intersection pinches along
/// it (as on main).
#[test]
fn a_notchs_root_across_a_face() {
    let p = plate(1.0);
    let notch = prism(
        xz(DVec3::new(0.0, 2.5, 0.0)),
        &tri(1.0, [(-0.25, -0.5), (0.25, -0.5), (0.0, 0.0)]),
        5.0,
    );
    let n = op(
        &cuboid(DVec3::new(-1.0, -1.0, -0.5), DVec3::new(1.0, 1.0, 0.5)),
        &notch,
        BooleanOp::Subtract,
    );
    let along = edge_refusal("(0.000000, 0.000000, 0.000000) mm");
    assert_eq!(
        six(&p, &n, 1.0),
        [
            "(11, 24, 16, 1) V=18.000000000 A=52.000000000",
            "(11, 24, 16, 1) V=18.000000000 A=52.000000000",
            along.as_str(),
            along.as_str(),
            "(16, 35, 22, 1) V=14.250000000 A=52.986067977",
            "(6, 12, 8, 1) V=2.000000000 A=12.000000000",
        ]
    );
}

/// Golden 9: a knife from beyond the plate's edge to over a hole, its
/// edge on the face between the outer loop and the hole's. The difference
/// is the holed plate (main kept the edge as a slit: 10 faces, 27 edges
/// and 18 vertices, and validated).
#[test]
fn a_knife_between_a_faces_edge_and_a_hole_leaves_no_slit() {
    let (h, kn) = (holed_plate(1.0), knife(1.0, 0.0));
    let holed = "(10, 24, 16, 1) V=15.840000000 A=49.280000000";
    let kf = "(5, 9, 6, 1) V=0.090000000 A=1.801665383";
    let along = edge_refusal("(-1600.000000, -1000.000000, 0.000000) mm");
    assert_eq!(
        six(&h, &kn, 1.0),
        [along.as_str(), along.as_str(), "empty", "empty", holed, kf]
    );
    is_operand(&h, &kn, BooleanOp::Subtract, Operand::First);
    as_if_apart(&kn, &knife(1.0, 1e-3), &h, BooleanOp::Subtract);
}

/// Golden 10: contacts on P's boundary: W's edge running out past P's
/// edge (a bridge from the face's boundary), C's corner on the middle of
/// P's top edge and on P's corner, and a wedge's edge along P's top edge.
#[test]
fn contacts_on_a_faces_boundary() {
    let p = plate(1.0);
    let out = prism(
        xz(DVec3::new(0.0, 2.5, 0.0)),
        &tri(1.0, [(0.3, 0.0), (0.8, 0.5), (-0.2, 0.5)]),
        1.0,
    );
    let along = edge_refusal("(300.000000, 1750.000000, 0.000000) mm");
    assert_eq!(
        six(&p, &out, 1.0),
        [along.as_str(), along.as_str(), "empty", "empty", P, W]
    );
    for (corner, at) in [
        (
            DVec3::new(0.0, 2.0, 0.0),
            "(0.000000, 2000.000000, 0.000000) mm",
        ),
        (
            DVec3::new(2.0, 2.0, 0.0),
            "(2000.000000, 2000.000000, 0.000000) mm",
        ),
    ] {
        let c = corner_cube(corner, 0.5, false);
        let at = point_refusal(at);
        assert_eq!(
            six(&p, &c, 1.0),
            [at.as_str(), at.as_str(), "empty", "empty", P, C]
        );
        is_operand(&p, &c, BooleanOp::Subtract, Operand::First);
        is_operand(&c, &p, BooleanOp::Subtract, Operand::First);
    }
    let on_edge = prism(
        yz(DVec3::new(-0.5, 0.0, 0.0)),
        &tri(1.0, [(2.0, 0.0), (2.5, 0.5), (1.5, 0.5)]),
        1.0,
    );
    let along = edge_refusal("(0.000000, 2000.000000, 0.000000) mm");
    assert_eq!(
        six(&p, &on_edge, 1.0),
        [along.as_str(), along.as_str(), "empty", "empty", P, W]
    );
    for w in [&out, &on_edge] {
        // P's edge keeps P's edge alone as its origin, not the wedge's.
        is_operand(&p, w, BooleanOp::Subtract, Operand::First);
        is_operand(w, &p, BooleanOp::Subtract, Operand::First);
    }
}

/// Golden 12: a wedge inside a box with its edge along the box's top front
/// edge, its faces going into the box: the union is the box, and the box's
/// edge keeps only the box's edge as its origin though the wedge's is the
/// first operand's (before ADR 0024 the wedge's edge was the origin).
#[test]
fn a_nested_edge_on_an_edge() {
    let q = cuboid(DVec3::new(-2.0, -2.0, -2.0), DVec3::new(2.0, 2.0, 0.0));
    let a = prism(
        yz(DVec3::new(-0.5, 0.0, 0.0)),
        &tri(1.0, [(-2.0, 0.0), (-1.0, -0.5), (-1.5, -1.0)]),
        1.0,
    );
    let box_q = "(6, 12, 8, 1) V=32.000000000 A=64.000000000";
    let wedge_a = "(5, 9, 6, 1) V=0.375000000 A=3.693174759";
    let along = edge_refusal("(0.000000, -2000.000000, 0.000000) mm");
    assert_eq!(
        six(&a, &q, 1.0),
        [box_q, box_q, wedge_a, wedge_a, "empty", along.as_str()]
    );
    is_operand(&a, &q, BooleanOp::Union, Operand::Second);
    is_operand(&q, &a, BooleanOp::Union, Operand::First);
    is_operand(&a, &q, BooleanOp::Intersect, Operand::First);
}

/// Golden 13, review's counterexamples to designs that were not built.
/// CE1: a block with a through V-groove crossing P's top face, the root
/// inside the face; the block minus P is two prisms touching along the
/// root, refused as on main.
#[test]
fn a_grooved_block_crossing_a_face() {
    let p = plate(1.0);
    let block = op(
        &cuboid(DVec3::new(0.2, -0.5, -0.5), DVec3::new(1.2, 0.5, 0.5)),
        &prism(
            xz(DVec3::new(0.0, 1.0, 0.0)),
            &tri(1.0, [(0.7, 0.0), (1.0, 0.6), (0.4, 0.6)]),
            2.0,
        ),
        BooleanOp::Subtract,
    );
    let union = "(16, 35, 22, 1) V=16.375000000 A=50.368033989";
    let lower = "(6, 12, 8, 1) V=0.500000000 A=4.000000000";
    let along = edge_refusal("(700.000000, 0.000000, 0.000000) mm");
    assert_eq!(
        six(&block, &p, 1.0),
        [
            union,
            union,
            lower,
            lower,
            along.as_str(),
            "(11, 24, 16, 1) V=15.500000000 A=50.000000000",
        ]
    );
}

/// CE2: W hanging inside P from its top face, its end flush with P's
/// side: P − W pinches along the edge (main returned it).
#[test]
fn a_wedge_hanging_inside_a_plate_flush_with_its_side() {
    let p = plate(1.0);
    let w = prism(
        xz(DVec3::new(0.0, 2.0, 0.0)),
        &tri(1.0, [(0.3, 0.0), (0.8, -0.5), (-0.2, -0.5)]),
        1.0,
    );
    let along = edge_refusal("(300.000000, 1500.000000, 0.000000) mm");
    assert_eq!(six(&p, &w, 1.0), [P, P, W, W, along.as_str(), "empty"]);
}

/// CE3: an L-shaped post standing on P with a knife from the foot of its
/// re-entrant edge lying on P: the post's foot lies on P, the knife is a
/// bridge from that footprint, and the union pinches along the knife.
#[test]
fn a_knife_from_a_posts_foot() {
    let p = plate(1.0);
    let post = prism(
        Frame::WORLD,
        &[
            DVec2::new(-1.0, -1.0),
            DVec2::new(1.0, -1.0),
            DVec2::new(1.0, 0.0),
            DVec2::new(0.0, 0.0),
            DVec2::new(0.0, 1.0),
            DVec2::new(-1.0, 1.0),
        ],
        1.0,
    );
    let d = DVec3::new(1.0, 1.0, 0.0).normalize();
    let blade = prism(
        Frame {
            origin: DVec3::new(-0.2, -0.2, 0.0),
            x: DVec3::new(-1.0, 1.0, 0.0).normalize(),
            y: DVec3::Z,
            z: d,
        },
        &tri(1.0, [(0.0, 0.0), (0.3, 0.5), (-0.3, 0.5)]),
        1.2,
    );
    let a = op(&post, &blade, BooleanOp::Union);
    let along = edge_refusal("(324.264069, 324.264069, 0.000000) mm");
    let post_knife = "(12, 28, 18, 1) V=3.122573593 A=15.292813786";
    assert_eq!(
        six(&p, &a, 1.0),
        [
            along.as_str(),
            along.as_str(),
            "empty",
            "empty",
            P,
            post_knife
        ]
    );
    is_operand(&p, &a, BooleanOp::Subtract, Operand::First);
    is_operand(&a, &p, BooleanOp::Subtract, Operand::First);
}

/// A post with a knife whose cross-section meets the post's side at the
/// middle of its bottom edge: the union's side face has a loop through one
/// vertex twice, which the classification of a later boolean must take
/// (main refused it as "two vertices at one point of the face's plane").
#[test]
fn an_operand_face_whose_loop_passes_a_vertex_twice() {
    let p = plate(1.0);
    let post = cuboid(DVec3::new(-0.5, -0.5, 0.0), DVec3::new(0.5, 0.5, 1.0));
    let blade = prism(
        yz(DVec3::new(0.3, 0.0, 0.0)),
        &tri(1.0, [(0.0, 0.0), (0.3, 0.5), (-0.3, 0.5)]),
        1.2,
    );
    let a = op(&post, &blade, BooleanOp::Union);
    let along = edge_refusal("(1000.000000, 0.000000, 0.000000) mm");
    let post_knife = "(10, 22, 14, 1) V=1.150000000 A=7.766190379";
    assert_eq!(
        six(&p, &a, 1.0),
        [
            along.as_str(),
            along.as_str(),
            "empty",
            "empty",
            P,
            post_knife
        ]
    );
    is_operand(&p, &a, BooleanOp::Subtract, Operand::First);
    is_operand(&a, &p, BooleanOp::Subtract, Operand::First);
}

/// CE3′: a knife between two posts through P, lying on P's top face from
/// one hole to the other: P minus the posts and knife is P with two clean
/// holes, as with the knife raised clear.
#[test]
fn a_knife_between_two_posts() {
    let p = plate(1.0);
    let posts = |z: f64| {
        let one = cuboid(DVec3::new(-1.5, -0.25, -2.0), DVec3::new(-1.0, 0.25, 1.0));
        let two = cuboid(DVec3::new(1.0, -0.25, -2.0), DVec3::new(1.5, 0.25, 1.0));
        let blade = prism(
            yz(DVec3::new(-1.2, 0.0, z)),
            &tri(1.0, [(0.0, 0.0), (0.2, 0.5), (-0.2, 0.5)]),
            2.4,
        );
        op(&op(&one, &blade, BooleanOp::Union), &two, BooleanOp::Union)
    };
    let a = posts(0.0);
    let along = edge_refusal("(0.000000, 0.000000, 0.000000) mm");
    let inside = "(12, 24, 16, 2) V=0.500000000 A=5.000000000";
    assert_eq!(
        six(&p, &a, 1.0),
        [
            along.as_str(),
            along.as_str(),
            inside,
            inside,
            "(14, 36, 24, 1) V=15.500000000 A=51.000000000",
            "(27, 59, 38, 3) V=1.200000000 A=12.754065923",
        ]
    );
    // The knife's own apexes fall on the posts' cut edges in A − P, so
    // only P − A has the raised knife's topology.
    as_if_apart_second(&p, &a, &posts(0.1), BooleanOp::Subtract);
}

/// Golden 14: the reflex corner of an L-prism's foot at the centre of a
/// tilted slab's face: the foot crosses the slab's face through the corner
/// while the corner's edges all rise from it, an ordinary crossing that
/// main answered correctly in all four operations; unchanged.
#[test]
fn a_reflex_corner_crossing_a_face() {
    let n = DVec3::ONE.normalize();
    let x = DVec3::new(1.0, -1.0, 0.0).normalize();
    let y = n.cross(x);
    let slab = k()
        .make_box(
            &Frame {
                origin: -2.0 * x - 2.0 * y - n,
                x,
                y,
                z: n,
            },
            DVec3::new(4.0, 4.0, 1.0),
        )
        .expect("a slab");
    let l = prism(
        Frame::WORLD,
        &[
            DVec2::new(-1.0, -1.0),
            DVec2::new(1.0, -1.0),
            DVec2::new(1.0, 0.0),
            DVec2::new(0.0, 0.0),
            DVec2::new(0.0, 1.0),
            DVec2::new(-1.0, 1.0),
        ],
        1.0,
    );
    let union = "(18, 38, 24, 1) V=17.836539648 A=53.947441117";
    let common = "(6, 12, 8, 1) V=1.163460352 A=8.052558883";
    assert_eq!(
        six(&l, &slab, 1.0),
        [
            union,
            union,
            common,
            common,
            "(14, 26, 16, 2) V=1.836539648 A=11.267949192",
            "(10, 24, 16, 1) V=14.836539648 A=50.732050808",
        ]
    );
}

/// Golden 15: W's edge across an inverted wedge's edge, the two touching
/// at one point: the union pinches there, the differences are the
/// operands.
#[test]
fn two_edges_crossing_at_a_point() {
    let w = wedge(1.0, 0.0);
    let under = prism(
        yz(DVec3::new(-0.2, 0.0, 0.0)),
        &tri(1.0, [(0.0, 0.0), (0.5, -0.5), (-0.5, -0.5)]),
        1.0,
    );
    let at = point_refusal("(300.000000, 0.000000, 0.000000) mm");
    assert_eq!(
        six(&w, &under, 1.0),
        [at.as_str(), at.as_str(), "empty", "empty", W, W]
    );
    is_operand(&w, &under, BooleanOp::Subtract, Operand::First);
    is_operand(&under, &w, BooleanOp::Subtract, Operand::First);
}

/// Golden 16: a crossing and a touch in one face: a post through P
/// carrying a beam with a wedge hanging from it onto P's top face. The
/// union pinches along the wedge's edge; the rest is as with the beam and
/// the wedge raised clear.
#[test]
fn a_crossing_and_a_touch_in_one_face() {
    let p = plate(1.0);
    let carried = |z: f64| {
        let post = cuboid(DVec3::new(-1.5, -0.25, -2.0), DVec3::new(-1.0, 0.25, 1.5));
        let beam = cuboid(
            DVec3::new(-1.5, -0.25, 1.0 + z),
            DVec3::new(1.0, 0.25, 1.2 + z),
        );
        let hung = prism(
            xz(DVec3::new(0.0, 0.25, z)),
            &tri(1.0, [(0.3, 0.0), (0.8, 1.05), (-0.2, 1.05)]),
            0.5,
        );
        op(&op(&post, &beam, BooleanOp::Union), &hung, BooleanOp::Union)
    };
    let (a, raised) = (carried(0.0), carried(0.1));
    let got = six(&p, &a, 1.0);
    let along = edge_refusal("(300.000000, 0.000000, 0.000000) mm");
    assert_eq!([&got[0], &got[1]], [&along, &along]);
    let apart = six(&p, &raised, 1.0);
    assert_eq!(got[2..], apart[2..]);
    for o in [BooleanOp::Intersect, BooleanOp::Subtract] {
        as_if_apart(&a, &raised, &p, o);
        as_if_apart_second(&p, &a, &raised, o);
    }
}

/// Golden 17: many contacts in one face: five wedges, and a 4 × 4 grid of
/// cubes on their corners (main refused the wedges, a test point falling
/// on a contact, and returned the grid's union touching itself).
#[test]
fn many_contacts_in_one_face() {
    let p = plate(1.0);
    let small = |x0: f64| {
        prism(
            xz(DVec3::new(0.0, 0.5, 0.0)),
            &tri(1.0, [(x0, 0.0), (x0 + 0.2, 0.3), (x0 - 0.2, 0.3)]),
            1.0,
        )
    };
    let mut five = small(-1.6);
    for x0 in [-0.8, 0.0, 0.8, 1.6] {
        five = op(&five, &small(x0), BooleanOp::Union);
    }
    let along = edge_refusal("(-1600.000000, 0.000000, 0.000000) mm");
    let wedges = "(25, 45, 30, 5) V=0.300000000 A=6.205551275";
    assert_eq!(
        six(&p, &five, 1.0),
        [along.as_str(), along.as_str(), "empty", "empty", P, wedges]
    );
    is_operand(&p, &five, BooleanOp::Subtract, Operand::First);
    is_operand(&five, &p, BooleanOp::Subtract, Operand::First);
    let mut grid: Option<Body> = None;
    for i in 0..4 {
        for j in 0..4 {
            let c = corner_cube(
                DVec3::new(-1.5 + i as f64, -1.5 + j as f64, 0.0),
                0.25,
                false,
            );
            grid = Some(match grid {
                None => c,
                Some(g) => op(&g, &c, BooleanOp::Union),
            });
        }
    }
    let grid = grid.expect("a grid");
    let at = point_refusal("(-1500.000000, -1500.000000, 0.000000) mm");
    let cubes = "(96, 192, 128, 16) V=0.250000000 A=6.000000000";
    assert_eq!(
        six(&p, &grid, 1.0),
        [at.as_str(), at.as_str(), "empty", "empty", P, cubes]
    );
    is_operand(&p, &grid, BooleanOp::Subtract, Operand::First);
    is_operand(&grid, &p, BooleanOp::Subtract, Operand::First);
}

/// Golden 18: goldens 1, 2 and 9 at millimetre and kilometre scale.
#[test]
fn contacts_at_millimetre_and_kilometre_scales() {
    for (s, at) in [
        (
            1e-3,
            [
                "0.300000, 0.000000",
                "0.300000, 0.200000",
                "-1.600000, -1.000000",
            ],
        ),
        (
            1e3,
            [
                "300000.000000, 0.000000",
                "300000.000000, 200000.000000",
                "-1600000.000000, -1000000.000000",
            ],
        ),
    ] {
        let place = |xy: &str| format!("({xy}, 0.000000) mm");
        let p = plate(s);
        let (w, c) = (wedge(s, 0.0), cube_on_corner(s, 0.0));
        let along = edge_refusal(&place(at[0]));
        assert_eq!(
            six(&p, &w, s),
            [along.as_str(), along.as_str(), "empty", "empty", P, W],
            "{s}"
        );
        let point = point_refusal(&place(at[1]));
        assert_eq!(
            six(&p, &c, s),
            [point.as_str(), point.as_str(), "empty", "empty", P, C],
            "{s}"
        );
        let (h, kn) = (holed_plate(s), knife(s, 0.0));
        let along = edge_refusal(&place(at[2]));
        assert_eq!(
            six(&h, &kn, s),
            [
                along.as_str(),
                along.as_str(),
                "empty",
                "empty",
                "(10, 24, 16, 1) V=15.840000000 A=49.280000000",
                "(5, 9, 6, 1) V=0.090000000 A=1.801665383",
            ],
            "{s}"
        );
        for (a, b) in [(&p, &w), (&p, &c), (&h, &kn)] {
            is_operand(a, b, BooleanOp::Subtract, Operand::First);
            is_operand(b, a, BooleanOp::Subtract, Operand::First);
        }
    }
}

/// Golden 19: a wedge whose edge pierces the top face of an L-prism on
/// its edge from the inner corner, four tolerances from the corner, with a
/// face across the top at a slight angle to that edge, passing within the
/// tolerance of the corner. The intersection's edge across the top joins
/// the four tolerances of the prism's edge to the long meeting beyond the
/// corner; it is where the faces meet, not a piece of the prism's edge
/// (on main it was named a piece of the edge).
#[test]
fn a_wedge_leaving_an_edge_near_a_corner() {
    let eps = k().tolerance().linear;
    let l = prism(
        Frame::WORLD,
        &[
            DVec2::new(0.0, 0.0),
            DVec2::new(1.0, 0.0),
            DVec2::new(1.0, 1.0),
            DVec2::new(2.0, 1.0),
            DVec2::new(2.0, 2.0),
            DVec2::new(0.0, 2.0),
        ],
        1.0,
    );
    // The wedge's edge runs down into the prism's arm along y through a
    // point four tolerances along the edge from the corner (1, 1, 1); one
    // face holds the line across the top through it toward (-1, 0.1).
    let at = DVec3::new(1.0 + 4.0 * eps, 1.0, 1.0);
    let x = DVec3::new(-1.0, 0.1, 0.0).normalize();
    let z = DVec3::new(0.0, 1.0, -1.0)
        .reject_from_normalized(x)
        .normalize();
    let y = z.cross(x);
    let side = y.dot(DVec3::new(1.0, 1.5, 1.0) - at).signum();
    let (c, s) = (60f64.to_radians().cos(), 60f64.to_radians().sin());
    let wedge = prism(
        Frame {
            origin: at - 3.0 * z,
            x,
            y,
            z,
        },
        &[
            DVec2::ZERO,
            DVec2::new(3.0, 0.0),
            DVec2::new(3.0 * c, 3.0 * s * side),
        ],
        6.0,
    );
    let union = "(15, 37, 24, 1) V=25.964366372 A=71.757842260";
    let common = "(6, 12, 8, 1) V=0.418319530 A=4.036386340";
    assert_eq!(
        six(&l, &wedge, 1.0),
        [
            union,
            union,
            common,
            common,
            "(12, 28, 18, 1) V=2.581680470 A=13.923572009",
            "(9, 21, 14, 1) V=22.964366372 A=61.870656619",
        ]
    );
    // The edges across the top from the point: where the top meets each
    // of the wedge's faces at the edge, with no other origin.
    let r = run(&l, &wedge, BooleanOp::Intersect);
    let mesh = k()
        .tessellate(&r, &TessellationOptions::default())
        .expect("a mesh");
    let origins = k().origins(&r);
    let across: Vec<EdgeOrigin> = (mesh.edges.iter().zip(&origins.edges))
        .filter(|(p, _)| {
            p.iter().all(|q| (q.z - 1.0).abs() < 1e-9)
                && p.iter().any(|q| (*q - at).length() < 1e-9)
        })
        .map(|(_, &o)| o)
        .collect();
    let meet = |second| EdgeOrigin::Meet {
        first: FaceIx(7),
        second: FaceIx(second),
    };
    assert_eq!(across, [meet(2), meet(0)]);
    assert_eq!(origins.merged, []);
}

/// Golden 20: a prism inside the unit cube with its edge along the top
/// face, its ends out through the sides x = 0 and x = 1, and both its
/// faces at the edge falling to the same side, at slopes of 1 in 2 and 3
/// in 4 (as contact 469 of the contact fuzz, at a hundred times its default
/// run). Raised 0.75 tolerances, the edge lies on the top face within the
/// tolerance, but the prism's faces cross the top's plane 1.5 and 1
/// tolerances beside it, on their side, where the cube's top edges cross
/// them: 1.25 tolerances along those edges from where the prism's edge
/// leaves the top, so the two points are not one. The top's meeting with
/// the prism's faces ran between those crossings, a copy of the prism's
/// edge on its faces beside it; the intersection and the prism minus the
/// cube were inconsistent ("an edge with no origin"), and the other two
/// refused as touching along the copy. Every operation is refused at the
/// crossing. Lowered 0.75 tolerances, the faces cross the top's plane
/// beside the edge off them, where the meeting lies off the prism's faces
/// and splits the top into pieces of one class: the outcomes are the exact
/// contact's.
#[test]
fn an_edge_on_a_face_leaving_it_where_its_faces_cross_beside_it() {
    let eps = k().tolerance().linear;
    let cube = cuboid(DVec3::ZERO, DVec3::ONE);
    let rod = |z: f64| {
        prism(
            yz(DVec3::new(-0.5, 0.0, z)),
            &[
                DVec2::new(0.5, 1.0),
                DVec2::new(0.1, 0.8),
                DVec2::new(0.1, 0.7),
            ],
            2.0,
        )
    };
    // Refused where the first operand's faces meet the other's first: the
    // cube's side x = 0, or the prism's at x = 1.
    let at = |x: &str| {
        format!(
            "refused: an edge of one body lies on a face of the other and leaves it at a point \
             the tolerance cannot place near ({x}, 499.999875, 1000.000000) mm"
        )
    };
    let (c, r) = (at("0.000000"), at("1000.000000"));
    assert_eq!(
        six(&cube, &rod(0.75 * eps), 1.0),
        [c.clone(), r.clone(), c.clone(), r.clone(), c, r]
    );
    let lowered = six(&cube, &rod(-0.75 * eps), 1.0);
    let union = "(14, 32, 20, 1) V=1.020000000 A=7.047213549";
    let common = "(5, 9, 6, 1) V=0.020000000 A=1.087213615";
    assert_eq!(
        lowered,
        [
            union.into(),
            union.into(),
            common.into(),
            common.into(),
            edge_refusal("(500.000000, 500.000000, 999.999944) mm"),
            "(10, 18, 12, 2) V=0.020000000 A=1.127213605".into(),
        ]
    );
    // The exact contact's outcomes: the same counts and refusals, at places
    // within the tolerance of its own.
    let head = |s: &String| match s.split_once(" near ") {
        Some((text, _)) => text.to_string(),
        None => s.split(" V=").next().unwrap_or_default().to_string(),
    };
    let exact = six(&cube, &rod(0.0), 1.0);
    assert_eq!(lowered.each_ref().map(head), exact.each_ref().map(head));
}

/// Golden 21: golden 20 from outside, as contact 323 of the fuzz was: the
/// prism lies across the floor of a groove with its ends in the groove's
/// walls, its faces rising to the same side at slopes of 1 in 2 and 3 in
/// 4. Lowered 0.75 tolerances, the faces cross the floor's plane beside
/// the edge on them and every operation is refused at the crossing (on the
/// code before, the intersection and the prism minus the block were
/// inconsistent, and the other two refused as touching along the copy of
/// the edge); raised 0.75 tolerances, the outcomes are the exact
/// contact's.
#[test]
fn an_edge_across_a_grooves_floor_leaving_it_where_its_faces_cross_beside_it() {
    let eps = k().tolerance().linear;
    let block = op(
        &cuboid(DVec3::ZERO, DVec3::new(3.0, 1.0, 1.0)),
        &cuboid(DVec3::new(1.0, -1.0, 0.5), DVec3::new(2.0, 2.0, 2.0)),
        BooleanOp::Subtract,
    );
    let rod = |z: f64| {
        prism(
            yz(DVec3::new(0.5, 0.0, z)),
            &[
                DVec2::new(0.5, 0.5),
                DVec2::new(0.1, 0.8),
                DVec2::new(0.1, 0.7),
            ],
            2.0,
        )
    };
    // Refused where the first operand's faces meet the other's first: the
    // block's wall at x = 2, or the prism's face at x = 1.
    let at = |x: &str| {
        format!(
            "refused: an edge of one body lies on a face of the other and leaves it at a point \
             the tolerance cannot place near ({x}, 499.999875, 500.000000) mm"
        )
    };
    let (b, r) = (at("2000.000000"), at("1000.000000"));
    assert_eq!(
        six(&block, &rod(-0.75 * eps), 1.0),
        [b.clone(), r.clone(), b.clone(), r.clone(), b, r]
    );
    let raised = six(&block, &rod(0.75 * eps), 1.0);
    let pinch = edge_refusal("(1500.000000, 500.000000, 500.000056) mm");
    let ends = "(10, 18, 12, 2) V=0.020000000 A=1.127213605";
    assert_eq!(
        raised,
        [
            pinch.clone(),
            pinch,
            ends.into(),
            ends.into(),
            "(18, 44, 28, 1) V=2.480000000 A=15.047213549".into(),
            "(5, 9, 6, 1) V=0.020000000 A=1.087213615".into(),
        ]
    );
    let head = |s: &String| match s.split_once(" near ") {
        Some((text, _)) => text.to_string(),
        None => s.split(" V=").next().unwrap_or_default().to_string(),
    };
    let exact = six(&block, &rod(0.0), 1.0);
    assert_eq!(raised.each_ref().map(head), exact.each_ref().map(head));
}

/// A blade: a thin prism whose edge runs at 45° across P's top edge
/// x = 2, through (2, 0, `z`), from 0.71 inside P's top face to as far
/// outside; its two faces leave the edge to the same side at 12° and 14°
/// to the top face, rising (`down`: falling).
fn blade(z: f64, down: bool) -> Body {
    let d = DVec3::new(1.0, 1.0, 0.0).normalize();
    let frame = Frame {
        origin: DVec3::new(2.0, 0.0, z) - d,
        x: DVec3::new(-d.y, d.x, 0.0),
        y: DVec3::Z,
        z: d,
    };
    let at = |deg: f64| {
        let t = if down { -deg } else { deg }.to_radians();
        DVec2::new(t.cos(), t.sin()) * 0.5
    };
    prism(frame, &[DVec2::ZERO, at(12.0), at(14.0)], 2.0)
}

/// Golden 22: an edge lying on a face within the tolerance and crossing
/// the face's boundary at a slant, its faces at a shallow angle to the
/// face. The curves where the face meets them run from the edge's end
/// inside the face to the boundary, ending where the boundary crosses
/// their planes, a few tolerances from the edge's own crossing: a bridge
/// carrying no edge, which is the edge's contact the tolerance moved
/// (found by the contact fuzz, as a kernel error). Resting on P (0.75
/// tolerances above), the union pinches along the edge; nested in P
/// (0.75 below, falling), P − blade pinches (on main it returned a body
/// touching itself); 0.75 above but falling, the faces cross the top
/// face a few tolerances from the edge, and every operation is refused
/// (on main P − blade touched itself).
#[test]
fn an_edge_lying_on_a_face_leaving_it_at_a_slant() {
    let eps = k().tolerance().linear;
    let p = plate(1.0);
    let rest = blade(0.75 * eps, false);
    let b = "(5, 9, 6, 1) V=0.008724874 A=2.043629687";
    let along = edge_refusal("(1646.446609, -353.553622, 0.000037) mm");
    assert_eq!(
        six(&p, &rest, 1.0),
        [along.as_str(), along.as_str(), "empty", "empty", P, b]
    );
    is_operand(&p, &rest, BooleanOp::Subtract, Operand::First);
    is_operand(&rest, &p, BooleanOp::Subtract, Operand::First);
    let nested = blade(-0.75 * eps, true);
    let common = "(5, 9, 6, 1) V=0.005779097 A=1.280040940";
    let union = "(10, 22, 14, 1) V=16.002945777 A=48.763588643";
    let along = edge_refusal("(1646.446609, -353.553622, -0.000037) mm");
    assert_eq!(
        six(&p, &nested, 1.0),
        [
            union,
            union,
            common,
            common,
            along.as_str(),
            "(5, 9, 6, 1) V=0.002945777 A=0.775927591",
        ]
    );
    let crossing = blade(0.75 * eps, true);
    let leaves = "refused: an edge of one body lies on a face of the other and leaves it at a \
                  point the tolerance cannot place near (2000.000000, 0.000462, 0.000000) mm";
    assert_eq!(six(&p, &crossing, 1.0), [leaves; 6]);
}
