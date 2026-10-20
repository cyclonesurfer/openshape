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

use std::f64::consts::{FRAC_PI_2, PI};

use opendatum_kernel::{
    Body, BooleanOp, BrepKernel, CurveKind, DVec2, DVec3, EdgeOrigin, ExtrudeSpec, FaceIx,
    FaceOrigin, Frame, Kernel, KernelError, Merged, Operand, Origins, Profile, ProfileCurve,
    TessellationOptions, Tolerance, VertexOrigin,
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
/// counts, volume over `s`³ and area over `s`²; a result must validate,
/// have origins that agree with the operands, and not touch itself, as
/// measured on a mesh within `chord` of it ([`alone_meshed`]).
fn outcome_meshed(a: &Body, b: &Body, o: BooleanOp, s: f64, chord: f64) -> String {
    match k().boolean(a, b, o) {
        Ok(r) => {
            assert_eq!(k().validate(&r), Vec::<String>::new(), "{o:?}");
            origins_agree(a, b, &r, o).unwrap_or_else(|e| panic!("{o:?}: {e}"));
            alone_meshed(&r, o, chord);
            format!(
                "{:?} V={:.9} A={:.9}",
                counts(&r),
                volume(&r) / (s * s * s),
                area(&r) / (s * s)
            )
        }
        Err(KernelError::InvalidInput(m)) if m == "the result would be empty" => "empty".into(),
        Err(KernelError::InvalidInput(m)) => format!("refused: {m}"),
        Err(KernelError::Unsupported(m)) => format!("unsupported: {m}"),
        Err(e) => panic!("{o:?}: {e:?}"),
    }
}

/// The outcomes of A ∪ B, B ∪ A, A ∩ B, B ∩ A, A − B and B − A.
fn six(a: &Body, b: &Body, s: f64) -> [String; 6] {
    six_meshed(a, b, s, CHORD)
}

/// [`six`] with [`outcome_meshed`].
fn six_meshed(a: &Body, b: &Body, s: f64, chord: f64) -> [String; 6] {
    use BooleanOp::*;
    [
        outcome_meshed(a, b, Union, s, chord),
        outcome_meshed(b, a, Union, s, chord),
        outcome_meshed(a, b, Intersect, s, chord),
        outcome_meshed(b, a, Intersect, s, chord),
        outcome_meshed(a, b, Subtract, s, chord),
        outcome_meshed(b, a, Subtract, s, chord),
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

/// How far the meshes the self-contact oracle measures may be from a
/// result's curves and surfaces, unless a test says otherwise.
const CHORD: f64 = 1e-4;

/// Asserts that a result does not touch itself
/// ([`self_contact::touches_itself`]).
fn alone(r: &Body, o: BooleanOp) {
    alone_meshed(r, o, CHORD);
}

/// [`alone`], measured on a mesh within `chord` of the result. The oracle
/// takes a point within `chord` and half the tolerance of a face with
/// curved edges or a curved face for a contact if it is on the face's
/// surface, so a finer mesh tells a point just outside such a face from
/// one on it.
fn alone_meshed(r: &Body, o: BooleanOp, chord: f64) {
    let opts = TessellationOptions {
        chord_tolerance: chord,
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

/// Golden 3: goldens 1 and 2 on P with a round hole elsewhere, which gives
/// its top face a curved side: the answers do not change but for the hole
/// (on main they were refused as unsupported tangent contacts).
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
/// edge on its faces beside it. On main the prism minus the cube came back
/// as one solid, its two ends joined over the top by a sliver 0.75
/// tolerances thick, and validated; the code before this fix found that
/// the copy had no origin. Every operation is refused at the crossing. Lowered 0.75 tolerances, the faces cross the top's plane
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
/// the edge on them and every operation is refused at the crossing (on
/// main the intersection came back as one solid, the prism's two pieces in
/// the walls joined along the floor by a sliver 0.75 tolerances thick, and
/// validated; on the code before this fix the intersection and the prism
/// minus the block were inconsistent); raised 0.75 tolerances, the
/// outcomes are the exact contact's.
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
/// tolerances above), the union pinches along the edge (main refused it
/// as touching at a point at the edge's far end); nested in P (0.75
/// below, falling), P − blade pinches (main returned a body that
/// validates with the blade's edge 0.75 tolerances under P's top face,
/// which it does not bound); 0.75 above but falling, the faces cross the
/// top face a few tolerances from the edge, and every operation is
/// refused (main returned P − blade with the edge 0.375 tolerances from
/// the top face).
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

/// Golden 23: a cube of side 0.3 on its corner on the plane of the unit
/// box's front face, 0.75 tolerances below its top face, so within the
/// tolerance of the box's edge between them. The cube's bottom face holds
/// the corner, tilted 18° from level; the box's edge crosses it at 14.5°
/// three tolerances from the corner, 1.5 tolerances from the cube's edge
/// in that face. The piece of the box's edge from the corner to the
/// crossing is within the tolerance of the bottom face, the box's faces
/// at it leaving that face to both sides: the edge crosses the face at a
/// slant where the tolerance does not decide, and every operation is
/// refused there (before, a kernel error: "the faces at a lying edge leave
/// the face it lies on to both sides"; found by the contact fuzz).
#[test]
fn a_corner_by_an_edge_that_crosses_its_face_at_a_slant() {
    let eps = k().tolerance().linear;
    let unit = cuboid(DVec3::ZERO, DVec3::ONE);
    let h = 0.75 * eps;
    let corner = DVec3::new(0.5, 0.0, 1.0 - h);
    // The cube's edges: `up` normal to its bottom face; `along` in it, 30°
    // from the direction from the corner to where the box's edge crosses
    // the face.
    let up = DVec3::new(-0.25, 0.187, 0.95).normalize();
    let crossing = DVec3::new(0.5 - up.z * h / up.x, 0.0, 1.0);
    let u = (crossing - corner).normalize();
    let (c, s) = (30f64.to_radians().cos(), 30f64.to_radians().sin());
    let along = u * c + up.cross(u) * s;
    let cube = k()
        .make_box(
            &Frame {
                origin: corner,
                x: along.cross(up),
                y: along,
                z: up,
            },
            DVec3::splat(0.3),
        )
        .expect("a cube");
    let refused = "refused: an edge of one body runs within the tolerance of a face of the other \
                   over a stretch without lying on it, near (500.000143, 0.000000, 999.999963) mm";
    assert_eq!(six(&unit, &cube, 1.0), [(); 6].map(|_| refused.to_string()));
}

/// Golden 24: a prism's edge lying on the unit box's front face, parallel
/// to the box's top front edge and 1.1 tolerances below it, its faces
/// rising on both sides of the front face (as contact 828 of the contact
/// fuzz). The box's top front edge then lies within the tolerance of both
/// of the prism's faces (0.57 and 0.41 tolerances from them), and the
/// box's top face meets each along that edge, from where the prism's ends
/// cross the top face's plane, points the tolerance puts at the prism's
/// corners, 1.1 tolerances below the edge: the box's edge leaves the
/// prism's faces where the tolerance cannot place the point, and every
/// operation is refused at one of the prism's ends (before, a kernel
/// error: the prism's faces at its edge left the front face to both
/// sides).
#[test]
fn an_edge_lying_on_a_face_by_its_edge_beyond_the_tolerance() {
    let eps = k().tolerance().linear;
    let unit = cuboid(DVec3::ZERO, DVec3::ONE);
    let blade = prism(
        yz(DVec3::new(0.25, 0.0, 1.0 - 1.1 * eps)),
        &[DVec2::ZERO, DVec2::new(0.3, 0.5), DVec2::new(-0.2, 0.5)],
        0.5,
    );
    let at = |x: &str| {
        format!(
            "refused: an edge of one body lies on a face of the other and leaves it at a point \
             the tolerance cannot place near ({x}, 0.000000, 999.999890) mm"
        )
    };
    let (far, near) = (at("750.000000"), at("250.000000"));
    assert_eq!(
        six(&unit, &blade, 1.0),
        [
            far.clone(),
            near.clone(),
            far.clone(),
            near.clone(),
            far,
            near
        ]
    );
}

/// A prism whose top face (z = 0, x ≤ 0, |y| ≤ 1) meets its slanted face
/// at `degrees` along the edge x = z = 0.
fn acute(degrees: f64) -> Body {
    let t = degrees.to_radians().tan();
    prism(
        xz(DVec3::new(0.0, 1.0, 0.0)),
        &[
            DVec2::new(0.0, 0.0),
            DVec2::new(-2.0, -2.0 * t),
            DVec2::new(-2.0, 0.0),
        ],
        2.0,
    )
}

/// A knife whose edge runs along x from x = −0.5 to 0.5, `depth`
/// tolerances under the plane z = 0, its faces falling from the edge.
fn knife_under(depth: f64) -> Body {
    let eps = k().tolerance().linear;
    prism(
        yz(DVec3::new(-0.5, 0.0, -depth * eps)),
        &[
            DVec2::new(0.0, 0.0),
            DVec2::new(-0.3, -0.5),
            DVec2::new(0.3, -0.5),
        ],
        1.0,
    )
}

/// Golden 25: an edge lying on a face, within the tolerance under it,
/// that leaves the body through a neighbour meeting the face at an acute
/// angle before it reaches their edge (contact 3003 of the contact fuzz at
/// a hundred times its default run, a tolerance inside). Nine tenths of a
/// tolerance down, the knife's edge leaves the prism through the slanted
/// face at x = −1.56 tolerances, within the tolerance of the top face but
/// 1.8 tolerances from its edge, and passes the edge itself 1.56
/// tolerances on: the tolerance cannot place where it leaves the top face,
/// and all six operations are refused there (they were kernel errors, the
/// knife's pieces past the slanted face disagreeing with the side the
/// contact has on the top face). A quarter of a tolerance down, the two
/// points are one and the operations are answered, the difference of the
/// prism pinching along the edge as with the knife resting on the face.
/// A cube's corner within the tolerance of both faces of a 20° edge,
/// three tolerances from it, is refused the same way.
#[test]
fn an_edge_leaving_a_face_through_an_acute_neighbour() {
    let a = acute(30.0);
    let leaves = "refused: an edge of one body lies on a face of the other and leaves it at a \
                  point the tolerance cannot place near (-0.000156, 0.000000, -0.000090) mm";
    assert_eq!(six(&a, &knife_under(0.9), 1.0), [leaves; 6]);
    let union = "(10, 22, 14, 1) V=2.451067746 A=14.985469630";
    let common = "(4, 6, 4, 1) V=0.008333331 A=0.318325044";
    assert_eq!(
        six(&a, &knife_under(0.25), 1.0),
        [
            union,
            union,
            common,
            common,
            &edge_refusal("(-250.000005, 0.000000, -0.000017) mm"),
            "(6, 12, 8, 1) V=0.141666669 A=1.947865326",
        ]
    );
    let eps = k().tolerance().linear;
    let z = -0.9 * eps;
    let at = DVec3::new(z / 20f64.to_radians().tan() - 0.5 * eps, 0.0, z);
    let two = "refused: a corner of one body lies within the tolerance of two faces of the other \
               away from their edge near (-0.000297, 0.000000, -0.000090) mm";
    assert_eq!(
        six(&acute(20.0), &corner_cube(at, 0.1, true), 1.0),
        [two; 6]
    );
}

/// A thin wedge whose edge runs level from `start` along the unit
/// horizontal `d` for `len`, its two faces leaving the edge at angles `t1`
/// < `t2` (degrees) from the horizontal `Z × d`, toward +z (as [`blade`]).
fn level_wedge(start: DVec3, d: DVec3, len: f64, (t1, t2): (f64, f64)) -> Body {
    let frame = Frame {
        origin: start,
        x: DVec3::Z.cross(d),
        y: DVec3::Z,
        z: d,
    };
    let leg = |t: f64| DVec2::new(t.to_radians().cos(), t.to_radians().sin()) * 0.5;
    prism(frame, &[DVec2::ZERO, leg(t1), leg(t2)], len)
}

/// A wedge over the plate's corner (2, 2, 0): its edge, `z` tolerances
/// up, crosses the plate's edge x = 2 `delta` tolerances from the corner,
/// running at `alpha` degrees from y.
fn corner_blade(delta: f64, alpha: f64, z: f64, faces: (f64, f64)) -> Body {
    let eps = k().tolerance().linear;
    let a = alpha.to_radians();
    let d = DVec3::new(a.sin(), a.cos(), 0.0);
    let at = DVec3::new(2.0, 2.0 - delta * eps, z * eps);
    level_wedge(at - 0.71 * d, d, 1.42, faces)
}

/// Golden 26: a thin wedge's edge lying on the plate's top and leaving it
/// across the edge x = 2 three tolerances from the corner, its faces
/// falling into the plate at 12° and 14° (found by review). The curves
/// where the top meets the blade's faces ended at the top's corner, and
/// the point of the blade's edge nearest the corner lay beyond the top,
/// past where the edge crosses its boundary: the contact placed along the
/// edge to that point took in a piece outside the top, and every
/// operation was a kernel error (a contact's pieces disagreeing on their
/// side). The contact now stops where the edge leaves the top. Exact and
/// 0.75 tolerances down, the outcomes are the same, the plate minus the
/// wedge pinching along the edge (on main, exact, it came back with the
/// wedge's edge lying on the top face, which does not bound it).
#[test]
fn a_blade_leaving_a_face_near_its_corner() {
    let p = plate(1.0);
    let faces = (-168.0, -166.0);
    for (z, area, at) in [
        (0.0, ["49.145490820", "1.157848476"], "0.000000"),
        (-0.75, ["49.145490714", "1.157848482"], "-0.000037"),
    ] {
        let union = format!("(11, 25, 16, 1) V=16.005463499 A={}", area[0]);
        let common = "(4, 6, 4, 1) V=0.000731161 A=0.308016471".to_string();
        assert_eq!(
            six(&p, &corner_blade(3.0, 30.0, z, faces), 1.0),
            [
                union.clone(),
                union,
                common.clone(),
                common,
                edge_refusal(&format!("(1822.500000, 1692.560832, {at}) mm")),
                format!("(6, 12, 8, 1) V=0.005463499 A={}", area[1]),
            ]
        );
    }
}

/// Golden 27: a thin wedge crossing the plate's side x = 2 eight tolerances
/// from its corner, at 60°, its faces falling at 20° and 23° (found by
/// review). The curves where the side meets the blade's two faces leave
/// the edge's crossing at a slant and end within the tolerance of each
/// other on the corner's vertical edge, so one segment carries both: the
/// blade's section by the side's plane is thinner than the tolerance
/// there. The intersection and the blade minus the plate had the blade's
/// two faces meet along that segment, an edge of neither body, a kernel
/// error (main answered, the side's piece dropped and the edge named
/// after the side, which does not bound it). Refused as a sliver in every
/// operation, exact, 0.75 tolerances down, and 1.5 down, where the edge
/// is clear of the top.
#[test]
fn a_blade_whose_section_by_a_face_is_a_sliver() {
    let p = plate(1.0);
    for (z, at) in [
        (0.0, "-0.000137"),
        (-0.75, "-0.000202"),
        (-1.5, "-0.000287"),
    ] {
        let sliver = format!(
            "refused: the bodies leave a sliver thinner than the tolerance near (2000.000000, \
             1999.999600, {at}) mm"
        );
        let b = corner_blade(8.0, 60.0, z, (-23.0, -20.0));
        assert_eq!(six(&p, &b, 1.0), [(); 6].map(|_| sliver.clone()), "{z}");
    }
}

/// A prism of thickness 0.5 whose bottom face lies on the plane through
/// the x axis tilted `t` degrees up toward +y, with the profile `pts`:
/// (u along x, w across, whether w is in tolerances).
fn tilted(t: f64, pts: &[(f64, f64, bool)]) -> Body {
    let eps = k().tolerance().linear;
    let p: Vec<DVec2> = (pts.iter())
        .map(|&(u, w, tol)| DVec2::new(u, if tol { w * eps } else { w }))
        .collect();
    prism(tilted_frame(t), &p, 0.5)
}

/// [`tilted`]'s frame: its x along the x axis, its y up the plane through
/// the x axis tilted `t` degrees toward +y.
fn tilted_frame(t: f64) -> Frame {
    let t = t.to_radians();
    Frame {
        origin: DVec3::ZERO,
        x: DVec3::X,
        y: DVec3::new(0.0, t.cos(), t.sin()),
        z: DVec3::new(0.0, -t.sin(), t.cos()),
    }
}

/// Golden 28: a prism tilted 10° about the x axis whose bottom face's
/// edge along x bends in the middle, its three vertices 5 tolerances
/// to either side of the axis, so within the tolerance of the plane
/// z = 0, on a box whose top's edge x = 0.45 crosses the axis (found by
/// review). The check of a lying run's ends measured the side against the
/// line from its first vertex to its last: with the bend away from the
/// face (reflex), the end at x = 0.45 lies 4 tolerances off the face
/// beside the edge there, but the line put it 5 inside, and every
/// operation was refused; with the bend into the face (convex), the end
/// lies 4 tolerances inside, which the line put outside, and other checks
/// refused it with other texts. Measured against the run's edge nearest
/// the end, the reflex bend gives the exact contact's outcomes and the
/// convex one is refused at the end.
#[test]
fn a_bent_edge_lying_on_a_face() {
    let a = cuboid(DVec3::new(0.45, -1.0, -1.0), DVec3::new(2.0, 1.0, 0.0));
    let bent = |w: f64| {
        tilted(
            10.0,
            &[
                (0.0, -w, true),
                (0.5, w, true),
                (1.0, -w, true),
                (1.0, 0.5, false),
                (0.0, 0.5, false),
            ],
        )
    };
    let touch = edge_refusal("(750.000000, 0.000000, 0.000000) mm");
    assert_eq!(
        six(&a, &bent(5.0), 1.0),
        [
            touch.clone(),
            touch,
            "empty".into(),
            "empty".into(),
            "(6, 12, 8, 1) V=3.100000000 A=13.300000000".into(),
            "(7, 15, 10, 1) V=0.250000000 A=2.500000500".into(),
        ]
    );
    let leaves =
        "refused: an edge of one body lies on a face of the other and leaves it at a point \
                  the tolerance cannot place near (450.000000, 0.000000, 0.000000) mm";
    assert_eq!(six(&a, &bent(-5.0), 1.0), [leaves; 6]);
}

/// Golden 29: the tilted prism's edge straight, 5 tolerances off the axis
/// on the far side from its face, and ending at a reflex corner, where its
/// face's boundary turns away and the axis runs on into the face (found
/// by review). The check of a lying run's ends looked only at a meeting
/// that is the run itself, so it passed this one, which goes on past the
/// run; the intersection then held a fin 0.55 long and under a tolerance
/// thick between the box's top and the prism's face (as on main). As
/// golden 20, every
/// operation is refused where the box's edge crosses the axis, 5
/// tolerances inside the face beside the run.
#[test]
fn an_edge_lying_on_a_face_ending_at_a_reflex_corner() {
    let a = cuboid(DVec3::new(0.45, -1.0, -1.0), DVec3::new(1.8, 1.0, 0.0));
    let b = tilted(
        10.0,
        &[
            (0.0, -5.0, true),
            (1.0, -5.0, true),
            (1.5, -0.5, false),
            (2.0, -0.5, false),
            (2.0, 0.5, false),
            (0.0, 0.5, false),
        ],
    );
    let leaves =
        "refused: an edge of one body lies on a face of the other and leaves it at a point \
                  the tolerance cannot place near (450.000000, 0.000000, 0.000000) mm";
    assert_eq!(six(&a, &b, 1.0), [leaves; 6]);
}

/// A wedge from x = 0.2 to 0.8 whose face rises at 10° from the x axis
/// toward +y, its edge on that face's plane `off` tolerances short of the
/// axis.
fn rising_wedge(off: f64) -> Body {
    let eps = k().tolerance().linear;
    let t = 10f64.to_radians();
    let up = DVec2::new(t.cos(), t.sin());
    let at = up * (-off * eps);
    prism(
        yz(DVec3::new(0.2, 0.0, 0.0)),
        &[
            at,
            at + up * 0.5,
            at + DVec2::new(-0.1, 1.0).normalize() * 0.5,
        ],
        0.6,
    )
}

/// Golden 30: an edge of each body lying on the other's face, the faces
/// at 10° (found by review): the box's top edge y = 4 tolerances on the
/// wedge's face, and the wedge's edge on the top, 4 tolerances to the
/// other side of the axis. The ends of the meeting are the runs' own
/// vertices and the end check passes them; the faces lie within the
/// tolerance of each other over the band between the edges, which
/// classifying their pieces refuses in every operation. With the wedge's
/// edge 4 tolerances to the box's side instead, the runs do not overlap
/// and the outcomes are the exact contact's along the edge.
#[test]
fn edges_of_both_bodies_lying_on_the_others_face() {
    let eps = k().tolerance().linear;
    let a = cuboid(DVec3::new(0.0, -1.0, -1.0), DVec3::new(1.0, 4.0 * eps, 0.0));
    let band = "refused: a face of one body runs within the tolerance of the other's surface \
                without lying on one of its faces, near (500.000000, -0.000001, -0.000000) mm";
    assert_eq!(six(&a, &rising_wedge(4.0), 1.0), [band; 6]);
    let touch = edge_refusal("(500.000000, 0.000394, 0.000069) mm");
    assert_eq!(
        six(&a, &rising_wedge(-4.0), 1.0),
        [
            touch.clone(),
            touch,
            "empty".into(),
            "empty".into(),
            "(6, 12, 8, 1) V=1.000000400 A=6.000001600".into(),
            "(5, 9, 6, 1) V=0.074789923 A=1.257389214".into(),
        ]
    );
}

/// Golden 31: a wedge whose edge passes half a tolerance under the corner
/// (1, 1, 1) of the unit cube into it, rising at about 9.6° to cross the
/// top face three tolerances from the corner (contact 86 of the contact
/// fuzz, from its triage). The corner is placed on the edge, so the edge's
/// piece inside the cube runs from the top face's corner to the corner of
/// the hole the wedge makes in it and lies on the face; A − B keeps the
/// top face and the wedge's faces along it, so it would touch itself
/// there and is refused (on main it was returned). With a round hole
/// through the cube's sides away from the corner, which gives the cube,
/// and so A − B and the unions, curved faces, A − B is refused at the same
/// place and the intersections and B − A are the same (the unions differ
/// by the hole): the check of the result's places first ran only on
/// results with no curved face.
#[test]
fn a_wedge_under_a_corner_crossing_a_face_beside_it() {
    let eps = k().tolerance().linear;
    let a = cuboid(DVec3::ZERO, DVec3::ONE);
    let d = DVec3::new(-1.0, -1.0, (2.0f64 / 35.0).sqrt()).normalize();
    let y = (DVec3::Z - d * d.z).normalize();
    let frame = Frame {
        origin: DVec3::ONE - DVec3::Z * (0.5 * eps) - d * 0.3,
        x: y.cross(d),
        y,
        z: d,
    };
    let b = prism(
        frame,
        &tri(1.0, [(0.0, 0.0), (0.3, -0.5), (-0.3, -0.5)]),
        0.8,
    );
    let rod = k()
        .extrude(
            &Profile::circle(yz(DVec3::new(-1.0, 0.0, 0.0)), DVec2::new(0.3, 0.3), 0.1),
            &ExtrudeSpec {
                distance: 3.0,
                symmetric: false,
            },
        )
        .expect("a rod");
    let holed = op(&a, &rod, BooleanOp::Subtract);
    let pinch = edge_refusal("(999.999895, 999.999895, 1000.000000) mm");
    let common = "(7, 14, 9, 1) V=0.067524482 A=1.051687610";
    // B − A: the wedge beyond the corner, and its tip through the top.
    let lobes = "(10, 17, 11, 2) V=0.052475518 A=1.142975557";
    for (x, union) in [
        (&a, "(13, 28, 18, 1) V=1.052475518 A=6.661264706"),
        (&holed, "(15, 34, 22, 1) V=1.021059592 A=7.226751384"),
    ] {
        assert_eq!(
            six(x, &b, 1.0),
            [union, union, common, common, &pinch, lobes]
        );
    }
}

/// The tilted prism of the known fault (design note 0012 §19, §20): its
/// face's edge from x = 0 to 1 lies `w` tolerances off the axis, away from
/// the face, between two reflex corners past which the face's boundary
/// turns away and the axis runs on into the face.
fn beside(w: f64) -> Body {
    beside_at(10.0, w)
}

/// [`beside`]'s prism tilted `t` degrees.
fn beside_at(t: f64, w: f64) -> Body {
    tilted(
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

/// Golden 32: golden 29's prism with its edge between two reflex corners,
/// so that the line where its face crosses the box's top passes beside
/// the edge with both its ends beyond it (§19's known fault). The strip of
/// the face between the edge and the line lies within the tolerance of
/// the top, and A ∩ B was one solid: the wedges beyond the corners joined
/// by a fin under a tolerance thick along the edge. With the edge 2 and 5
/// tolerances off the axis, the top holds the strip, and every operation
/// is refused (design note 0012 §20), as with the box's top ending on the
/// axis or 4 tolerances past it (found by review: that edge of the top
/// then lies on the face and is the meeting), or ending on the axis at x
/// = 0.5, beside the edge (A ∩ B a wedge and half the fin; found by
/// probing), and with golden 29's prism against its box ending on the axis
/// (found by review). Half a tolerance off the axis, edge and line are one
/// place, the strip is empty, and B rests on the box by that edge: the
/// unions pinch, and A ∩ B is the two wedges.
#[test]
fn an_edge_lying_on_a_face_beside_the_line_where_the_faces_cross() {
    let eps = k().tolerance().linear;
    let boxed = |y: f64| cuboid(DVec3::new(-0.25, -1.0, -1.0), DVec3::new(1.25, y, 0.0));
    let touch = edge_refusal("(500.000000, -0.000049, -0.000009) mm");
    assert_eq!(
        six(&boxed(1.0), &beside(0.5), 1.0),
        [
            touch.clone(),
            touch,
            "(8, 12, 8, 2) V=0.000918370 A=0.152569852".into(),
            "(8, 12, 8, 2) V=0.000918370 A=0.152569852".into(),
            "(10, 24, 16, 1) V=2.999081630 A=13.003600651".into(),
            "(13, 31, 20, 1) V=0.874081667 A=7.203506245".into(),
        ]
    );
    let at = band_refusal("(500.000000, -0.000197, -0.000035) mm");
    assert_eq!(six(&boxed(1.0), &beside(2.0), 1.0), [at.as_str(); 6]);
    let at = band_refusal("(500.000000, -0.000492, -0.000087) mm");
    for y in [1.0, 0.0, 4.0 * eps] {
        assert_eq!(six(&boxed(y), &beside(5.0), 1.0), [at.as_str(); 6]);
    }
    let a = cuboid(DVec3::new(-0.25, -1.0, -1.0), DVec3::new(0.5, 0.0, 0.0));
    let at = band_refusal("(249.999950, -0.000492, -0.000087) mm");
    assert_eq!(six(&a, &beside(5.0), 1.0), [at.as_str(); 6]);
    let a = cuboid(DVec3::new(0.45, -1.0, -1.0), DVec3::new(1.8, 0.0, 0.0));
    let b = tilted(
        10.0,
        &[
            (0.0, -5.0, true),
            (1.0, -5.0, true),
            (1.5, -0.5, false),
            (2.0, -0.5, false),
            (2.0, 0.5, false),
            (0.0, 0.5, false),
        ],
    );
    let at = band_refusal("(725.000050, -0.000492, -0.000087) mm");
    assert_eq!(six(&a, &b, 1.0), [at.as_str(); 6]);
}

/// A block from x = −0.25 to 1.25 and y = −1 to 1, its top at z = 0,
/// with a notch through it from x = −0.1 to 1.1 whose floor is at `y`.
fn notched(y: f64) -> Body {
    let pts = [
        (-0.25, -1.0),
        (-0.1, -1.0),
        (-0.1, y),
        (1.1, y),
        (1.1, -1.0),
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

/// Golden 33: golden 32's prism, its edge 5 tolerances off the axis, on
/// bodies whose top does or does not hold the strip between the edge and
/// the line (design note 0012 §20). A box from y = 0 to 1, its edge on the
/// axis and its material away from the prism's edge; a slab with a hole
/// whose back edge is on the axis; and a block with a notch whose floor is
/// on it (found by review): the prism's edge lies over no face of A, and B
/// rests on that edge of A. The unions pinch; the intersections are empty
/// and the differences the operands, but for the block, whose prongs
/// beside the notch hold the wedges. The hole's front edge on the axis,
/// its material on the prism's edge's side, and the notch's floor 2
/// tolerances off the axis toward the edge (found by probing), hold the
/// strip: every operation is refused. The notch's floor is itself an edge
/// of A lying on the prism's face 2 tolerances beside the axis, whose own
/// strip the face holds: with A first it is found there, with B first at
/// the prism's edge's corner.
#[test]
fn an_edge_lying_beside_the_line_where_the_faces_cross_over_a_gap_in_the_face() {
    let eps = k().tolerance().linear;
    let b = beside(5.0);
    let mirror = cuboid(DVec3::new(-0.25, 0.0, -1.0), DVec3::new(1.25, 1.0, 0.0));
    let slab = cuboid(DVec3::new(-1.0, -1.0, -1.0), DVec3::new(3.0, 1.0, 0.0));
    let holed = |y0: f64, y1: f64| {
        let hole = cuboid(DVec3::new(-0.7, y0, -2.0), DVec3::new(2.2, y1, 1.0));
        op(&slab, &hole, BooleanOp::Subtract)
    };
    let prism_b = "(9, 21, 14, 1) V=0.875000375 A=7.207107928";
    let touch = edge_refusal("(500.000000, 0.000000, 0.000000) mm");
    for (a, a_alone) in [
        (&mirror, "(6, 12, 8, 1) V=1.500000000 A=8.000000000"),
        (
            &holed(-0.6, 0.0),
            "(10, 24, 16, 1) V=6.260000000 A=31.520000000",
        ),
    ] {
        assert_eq!(
            six(a, &b, 1.0),
            [
                touch.clone(),
                touch.clone(),
                "empty".into(),
                "empty".into(),
                a_alone.into(),
                prism_b.into(),
            ]
        );
    }
    let touch = edge_refusal("(-49.999750, 0.000000, 0.000000) mm");
    let wedges = "(10, 22, 16, 2) V=0.000859597 A=0.131685433";
    assert_eq!(
        six(&notched(0.0), &b, 1.0),
        [
            touch.clone(),
            touch,
            wedges.into(),
            wedges.into(),
            "(16, 40, 26, 1) V=1.799140403 A=12.599497987".into(),
            "(15, 43, 30, 1) V=0.874140778 A=7.207609909".into(),
        ]
    );
    let at = band_refusal("(500.000000, -0.000492, -0.000087) mm");
    assert_eq!(six(&holed(0.0, 0.6), &b, 1.0), [at.as_str(); 6]);
    let (on_a, on_b) = (
        band_refusal("(500.000000, -0.000200, 0.000000) mm"),
        band_refusal("(0.000000, -0.000099, -0.000017) mm"),
    );
    assert_eq!(
        six(&notched(-2.0 * eps), &b, 1.0),
        [&on_a, &on_b, &on_a, &on_b, &on_a, &on_b].map(|s| s.as_str())
    );
}

/// Golden 32's prism tilted `t` degrees, its edge `d` tolerances off the
/// axis, with steps from `r` tolerances off it to the edge's corners,
/// `delta` tolerances along it.
fn stepped(t: f64, d: f64, r: f64, delta: f64) -> Body {
    let eps = k().tolerance().linear;
    tilted(
        t,
        &[
            (-0.5, -0.5, false),
            (-delta * eps, -r, true),
            (0.0, -d, true),
            (1.0, -d, true),
            (1.0 + delta * eps, -r, true),
            (1.5, -0.5, false),
            (2.0, -0.5, false),
            (2.0, 0.5, false),
            (-0.5, 0.5, false),
        ],
    )
}

/// Golden 32's box with a triangular hole through it under each corner of
/// the prism's lying edge (x = 0 and x = 1): its apex at (x, `apex`
/// tolerances), its base from x − `half` to x + `half` at y = `base`.
fn nicked(apex: f64, half: f64, base: f64) -> Body {
    let eps = k().tolerance().linear;
    let mut boxed = cuboid(DVec3::new(-0.25, -1.0, -1.0), DVec3::new(1.25, 1.0, 0.0));
    for x in [0.0, 1.0] {
        let frame = Frame {
            origin: DVec3::new(0.0, 0.0, -2.0),
            ..Frame::WORLD
        };
        let tri =
            [(x - half, base), (x, apex * eps), (x + half, base)].map(|(x, y)| DVec2::new(x, y));
        boxed = op(&boxed, &prism(frame, &tri, 3.0), BooleanOp::Subtract);
    }
    boxed
}

/// Golden 34: golden 32's prism where the box's top holds only part of the
/// strip between the edge and the line, or the line does not reach the
/// edge (found by the second review, design note 0012 §20, where points
/// measured at the edge's corners missed them): the edge's corners
/// replaced by short steps nearly square to the axis (3°, the edge 18
/// tolerances off it, the steps from 1.5; and 1°, 50 and 12); triangular
/// nicks in the box under the edge's corners, their apexes 2 to 4
/// tolerances off the axis toward the edge; and a box whose top's edge
/// crosses the axis at a slant at x = 1.2, beyond the edge, 6 tolerances
/// to either side of it at its ends, its material on either side. Main
/// answered the steps, the nicks and the slanted edge with the material on
/// the edge's side with A ∩ B one solid joined along the edge; every
/// operation is refused. A nick's apex is a reflex corner of the box's top
/// lying on the prism's face, and the face holds the top between it and
/// the axis, so with the box first the refusal is found there, at the
/// neck. With the edge half a tolerance off the axis, the nicked box gives
/// the contact.
#[test]
fn an_edge_lying_beside_the_line_where_the_faces_cross_held_in_part() {
    let eps = k().tolerance().linear;
    let a = cuboid(DVec3::new(-0.25, -1.0, -1.0), DVec3::new(1.25, 1.0, 0.0));
    for (b, at) in [
        (
            stepped(3.0, 18.0, 1.5, 0.5),
            "(1000.000025, -0.000974, -0.000051) mm",
        ),
        (
            stepped(3.0, 18.0, 1.5, -0.5),
            "(500.000000, -0.001798, -0.000094) mm",
        ),
        (
            stepped(1.0, 50.0, 12.0, 0.5),
            "(1000.000025, -0.003100, -0.000054) mm",
        ),
    ] {
        let at = band_refusal(at);
        assert_eq!(six(&a, &b, 1.0), [at.as_str(); 6]);
    }
    let on_b = band_refusal("(500.000000, -0.000492, -0.000087) mm");
    for (apex, half, base, at) in [
        (-2.0, 0.2, -0.1, "(0.000000, -0.000200, 0.000000) mm"),
        (-3.0, 0.1, -0.02, "(0.000000, -0.000300, 0.000000) mm"),
        (-4.0, 0.2, -0.3, "(0.000000, -0.000400, 0.000000) mm"),
    ] {
        let on_a = band_refusal(at);
        assert_eq!(
            six(&nicked(apex, half, base), &beside(5.0), 1.0),
            [&on_a, &on_b, &on_a, &on_b, &on_a, &on_b].map(|s| s.as_str())
        );
    }
    let touch = edge_refusal("(0.000117, -0.000049, -0.000009) mm");
    assert_eq!(
        six(&nicked(-4.0, 0.2, -0.3), &beside(0.5), 1.0),
        [
            touch.clone(),
            touch,
            "(8, 14, 10, 2) V=0.000918370 A=0.152569842".into(),
            "(8, 14, 10, 2) V=0.000918370 A=0.152569842".into(),
            "(16, 44, 30, 1) V=2.879081790 A=15.005820142".into(),
            "(13, 31, 20, 1) V=0.874081667 A=7.203506245".into(),
        ]
    );
    let frame = Frame {
        origin: DVec3::new(0.0, 0.0, -1.0),
        ..Frame::WORLD
    };
    let (lo, hi) = ((-1.05, -6.0 * eps), (3.45, 6.0 * eps));
    for (pts, at) in [
        (
            [(-1.05, -1.0), (3.45, -1.0), hi, lo],
            "(500.000000, -0.000492, -0.000087) mm",
        ),
        (
            [lo, hi, (3.45, 1.0), (-1.05, 1.0)],
            "(0.000000, -0.000159, -0.000028) mm",
        ),
    ] {
        let a = prism(frame, &pts.map(|(x, y)| DVec2::new(x, y)), 1.0);
        let at = band_refusal(at);
        assert_eq!(six(&a, &beside(5.0), 1.0), [at.as_str(); 6]);
    }
}

/// Golden 35: golden 32's prism, its edge 5 tolerances off the axis, over a
/// gap of the other body's top whose back edge is on the axis: a notch and
/// a hole whose walls are flush with the edge's corners, a hole flush with
/// one of them, and a pocket whose walls stand 0.1 beyond the corners and
/// whose far lip lies 0.3 tolerances beyond the edge (found by the second
/// review: measured at the edge's corners, the first rule took the walls
/// for the top holding the strip). The strip
/// lies over the gap, so nothing of A bounds it, and each is answered as
/// the contact: the unions pinch, A ∩ B is the two wedges, and the
/// differences are the operands less the wedges.
#[test]
fn an_edge_lying_beside_the_line_over_a_gap_flush_with_its_corners() {
    let eps = k().tolerance().linear;
    let b = beside(5.0);
    let touch = edge_refusal("(0.000250, 0.000000, 0.000000) mm");
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
    assert_eq!(
        six(&prism(frame, &pts, 1.0), &b, 1.0),
        [
            touch.clone(),
            touch.clone(),
            "(8, 16, 12, 2) V=0.000918373 A=0.152570232".into(),
            "(8, 16, 12, 2) V=0.000918373 A=0.152570232".into(),
            "(16, 38, 24, 1) V=1.999081627 A=13.003600551".into(),
            "(13, 35, 24, 1) V=0.874082002 A=7.203507268".into(),
        ]
    );
    let slab = cuboid(DVec3::new(-1.0, -1.0, -1.0), DVec3::new(3.0, 1.0, 0.0));
    let holed = |x1: f64| {
        let hole = cuboid(DVec3::new(0.0, -0.6, -2.0), DVec3::new(x1, 0.0, 1.0));
        op(&slab, &hole, BooleanOp::Subtract)
    };
    assert_eq!(
        six(&holed(1.0), &b, 1.0),
        [
            touch.clone(),
            touch.clone(),
            "(9, 17, 12, 2) V=0.018367401 A=1.158218263".into(),
            "(9, 17, 12, 2) V=0.018367401 A=1.158218263".into(),
            "(17, 41, 26, 1) V=7.381632599 A=30.142791107".into(),
            "(11, 27, 18, 1) V=0.856632974 A=7.064316769".into(),
        ]
    );
    assert_eq!(
        six(&holed(1.5), &b, 1.0),
        [
            touch.clone(),
            touch,
            "(9, 17, 12, 2) V=0.014693919 A=0.897160054".into(),
            "(9, 17, 12, 2) V=0.014693919 A=0.897160054".into(),
            "(16, 39, 25, 1) V=7.085306081 A=30.491508100".into(),
            "(12, 30, 20, 1) V=0.860306456 A=7.115599780".into(),
        ]
    );
    let lip = -5.0 * eps * 10f64.to_radians().cos() - 0.3 * eps;
    let mut a = cuboid(DVec3::new(-0.25, -1.0, -1.0), DVec3::new(1.25, 1.0, -0.3));
    for (lo, hi) in [
        ((-0.25, -1.0), (1.25, lip)),
        ((-0.25, 0.0), (1.25, 1.0)),
        ((-0.25, -0.5), (-0.1, 0.5)),
        ((1.1, -0.5), (1.25, 0.5)),
    ] {
        let part = cuboid(DVec3::new(lo.0, lo.1, -0.3), DVec3::new(hi.0, hi.1, 0.0));
        a = op(&a, &part, BooleanOp::Union);
    }
    assert_eq!(
        six(&a, &b, 1.0),
        [
            edge_refusal("(-49.999750, 0.000000, 0.000000) mm"),
            edge_refusal("(500.000000, -0.000492, -0.000087) mm"),
            "(8, 18, 14, 2) V=0.000918373 A=0.152570029".into(),
            "(8, 18, 14, 2) V=0.000918373 A=0.152570029".into(),
            "(16, 38, 24, 1) V=2.999081439 A=13.723600889".into(),
            "(13, 37, 26, 1) V=0.874082002 A=7.203507268".into(),
        ]
    );
}

/// Golden 36: an edge lying on a face beside the line where the faces
/// cross, with an edge of the other body lying on the face between them
/// (found by the third review, design note 0012 §20): boxes from y = −3
/// and −2 tolerances against golden 32's prism, its edge 5 tolerances off
/// the axis; a box from −20 against the prism at 1°, its edge 40 off; and
/// golden 33's notched block, its floor 5 off, against golden 30's wedge,
/// its edge 3 and 2 off. The other edge is the meeting the imprint makes
/// (§18.2), and the strip between the face's edge and that meeting lies
/// over no part of the other face, so each is answered as on main: the
/// unions pinch along the meeting, the intersections are empty and the
/// differences are the operands. Swept on to the line, the strips were
/// held, and every operation was refused.
#[test]
fn an_edge_lying_beside_the_line_with_the_meeting_between() {
    let eps = k().tolerance().linear;
    let boxed = |y: f64| cuboid(DVec3::new(-0.25, y, -1.0), DVec3::new(1.25, 1.0, 0.0));
    let rests = |at: &str, a: &str, b: &str| {
        let touch = edge_refusal(at);
        [
            touch.clone(),
            touch,
            "empty".into(),
            "empty".into(),
            a.into(),
            b.into(),
        ]
    };
    let prism_b = "(9, 21, 14, 1) V=0.875000375 A=7.207107928";
    for (d, at, a) in [
        (
            3.0,
            "(500.000000, -0.000300, 0.000000) mm",
            "(6, 12, 8, 1) V=1.500000450 A=8.000001500",
        ),
        (
            2.0,
            "(500.000000, -0.000200, 0.000000) mm",
            "(6, 12, 8, 1) V=1.500000300 A=8.000001000",
        ),
    ] {
        assert_eq!(
            six(&boxed(-d * eps), &beside(5.0), 1.0),
            rests(at, a, prism_b)
        );
    }
    assert_eq!(
        six(&boxed(-20.0 * eps), &beside_at(1.0, 40.0), 1.0),
        rests(
            "(500.000000, -0.002000, 0.000000) mm",
            "(6, 12, 8, 1) V=1.500003000 A=8.000010000",
            "(9, 21, 14, 1) V=0.875003000 A=7.207115953",
        )
    );
    for (off, at) in [
        (3.0, "(500.000000, -0.000295, -0.000052) mm"),
        (2.0, "(500.000000, -0.000197, -0.000035) mm"),
    ] {
        assert_eq!(
            six(&notched(-5.0 * eps), &rising_wedge(off), 1.0),
            rests(
                at,
                "(10, 24, 16, 1) V=1.800000600 A=12.600000200",
                "(5, 9, 6, 1) V=0.074789923 A=1.257389214",
            )
        );
    }
}

/// A plate on golden 32's plane at `t` degrees, its face from x = −0.5 to
/// 2 and from 0.5 to either side of the axis, with a hole through it from
/// x = 0 to 1 and from `w0` to `w1` tolerances off the axis.
fn tilted_holed(t: f64, w0: f64, w1: f64) -> Body {
    holed(tilted(t, &PLATE), t, w0, w1)
}

/// The profile of [`tilted_holed`]'s plate.
const PLATE: [(f64, f64, bool); 4] = [
    (-0.5, -0.5, false),
    (2.0, -0.5, false),
    (2.0, 0.5, false),
    (-0.5, 0.5, false),
];

/// `plate`, a prism on golden 32's plane at `t` degrees half a unit thick
/// to either side, with a hole through it from x = 0 to 1 and from `w0`
/// to `w1` tolerances off the axis.
fn holed(plate: Body, t: f64, w0: f64, w1: f64) -> Body {
    let eps = k().tolerance().linear;
    let t = t.to_radians();
    let n = DVec3::new(0.0, -t.sin(), t.cos());
    let frame = Frame {
        origin: -n * 0.8,
        x: DVec3::X,
        y: DVec3::new(0.0, t.cos(), t.sin()),
        z: n,
    };
    let hole = [(0.0, w0), (1.0, w0), (1.0, w1), (0.0, w1)].map(|(x, y)| DVec2::new(x, y * eps));
    op(&plate, &prism(frame, &hole, 1.6), BooleanOp::Subtract)
}

/// [`tilted`]'s prism turned over: its face from the same profile on the
/// same plane, its body below the face.
fn tilted_below(t: f64, pts: &[(f64, f64, bool)]) -> Body {
    let eps = k().tolerance().linear;
    let t = t.to_radians();
    let frame = Frame {
        origin: DVec3::ZERO,
        x: -DVec3::X,
        y: DVec3::new(0.0, t.cos(), t.sin()),
        z: DVec3::new(0.0, t.sin(), -t.cos()),
    };
    let mut p: Vec<DVec2> = (pts.iter())
        .map(|&(u, w, tol)| DVec2::new(-u, if tol { w * eps } else { w }))
        .collect();
    p.reverse();
    prism(frame, &p, 0.5)
}

/// Golden 37: boundaries of a face lying on the other's plane beside the
/// line where they cross that are not runs of its outer loop (found by the
/// third review, design note 0012 §20): a hole in a plate at 1° from 50 to
/// 20 tolerances off the axis, its loop within 0.9 tolerances of the box's
/// top throughout; and a reflex corner of golden 32's prism's face 5
/// tolerances off the axis at the apex of a wide notch and of a narrow
/// one, its edges leaving the top's tolerance. The box's top holds the
/// face between each and the line, and main answered with A ∩ B one
/// solid: a flap under a tolerance thick from the hole to the line, or a
/// neck from the corner to the line joining the wedges beside it. Every
/// operation is refused, as with the hole reaching 200 tolerances off the
/// axis, its near edge a run. With the hole's near edge half a tolerance
/// off the axis there is no flap, and with a hole in the box under the
/// corner to a tenth of a tolerance off the axis the top does not hold the
/// neck: each is answered as on main.
#[test]
fn a_hole_or_a_corner_lying_beside_the_line_where_the_faces_cross() {
    let eps = k().tolerance().linear;
    let a = cuboid(DVec3::new(-0.25, -1.0, -1.0), DVec3::new(1.25, 1.0, 0.0));
    let at = band_refusal("(500.000000, -0.002000, -0.000035) mm");
    for w0 in [-50.0, -200.0] {
        assert_eq!(
            six(&a, &tilted_holed(1.0, w0, -20.0), 1.0),
            [at.as_str(); 6]
        );
    }
    let touch = edge_refusal("(500.000000, -0.000050, -0.000001) mm");
    let flap = "(5, 15, 12, 1) V=0.003272825 A=1.517559231";
    assert_eq!(
        six(&a, &tilted_holed(1.0, -50.0, -0.5), 1.0),
        [
            touch.clone(),
            touch,
            flap.into(),
            flap.into(),
            "(9, 25, 18, 1) V=2.996727175 A=13.008613286".into(),
            "(13, 36, 24, 1) V=1.246724700 A=9.491381721".into(),
        ]
    );
    let notch = |half: f64| {
        let mut pts = vec![(-0.5, -0.5, false)];
        if half < 1.0 {
            pts.push((0.5 - half, -0.5, false));
        }
        pts.extend([
            (0.5, -5.0, true),
            (0.5 + half, -0.5, false),
            (2.0, -0.5, false),
            (2.0, 0.5, false),
            (-0.5, 0.5, false),
        ]);
        tilted(10.0, &pts)
    };
    let at = band_refusal("(500.000000, -0.000492, -0.000087) mm");
    for half in [1.0, 0.1] {
        assert_eq!(six(&a, &notch(half), 1.0), [at.as_str(); 6]);
    }
    let hole = cuboid(
        DVec3::new(0.45, -0.3, -2.0),
        DVec3::new(0.55, -0.1 * eps, 1.0),
    );
    let touch = edge_refusal("(474.999510, -0.000010, 0.000000) mm");
    let wedges = "(10, 18, 12, 2) V=0.006197171 A=0.644425638";
    assert_eq!(
        six(&op(&a, &hole, BooleanOp::Subtract), &notch(1.0), 1.0),
        [
            touch.clone(),
            touch,
            wedges.into(),
            wedges.into(),
            "(15, 40, 26, 1) V=2.963802830 A=13.765973465".into(),
            "(14, 36, 24, 1) V=0.993803079 A=7.592061282".into(),
        ]
    );
}

/// Golden 38: golden 32's run and golden 37's corner 5 tolerances off the
/// axis, and golden 37's hole 20 to 50 off it in its plate at 1°, to
/// either side of the axis, on the prism or plate as in those goldens and
/// turned over, its body below its face (design note 0012 §20). The layer between the face and the
/// box's top over a strip lies in the prism's body on one side of the axis
/// only: below it for the prism above its face (A ∩ B holds it, golden
/// 37), above it for the prism below (B − A holds it). There every
/// operation is refused. On the other side the layer lies outside the
/// prism, in the box or in neither body, and the corner and the hole are
/// answered as on main: where the layer is in neither (the prism above,
/// the strip above the axis) the unions pinch along the strip's edge and
/// are refused as touching, and where it is in the box (the prism below,
/// the strip below the axis) A − B is. A run there is refused as well:
/// main failed every operation on it ("the faces at a lying edge leave the
/// face it lies on to both sides").
#[test]
fn a_run_a_hole_or_a_corner_beside_the_line_on_either_side() {
    let a = cuboid(DVec3::new(-0.25, -1.0, -1.0), DVec3::new(1.25, 1.0, 0.0));
    let body = |above: bool, pts: &[(f64, f64, bool)]| {
        if above {
            tilted(10.0, pts)
        } else {
            tilted_below(10.0, pts)
        }
    };
    let run = |w: f64| {
        let mut pts = vec![
            (-0.5, -0.5, false),
            (0.0, -w.abs(), true),
            (1.0, -w.abs(), true),
            (1.5, -0.5, false),
            (2.0, -0.5, false),
            (2.0, 0.5, false),
            (-0.5, 0.5, false),
        ];
        if w > 0.0 {
            pts = pts.iter().map(|&(u, w, t)| (u, -w, t)).rev().collect();
        }
        pts
    };
    let apex = |w: f64| {
        let mut pts = vec![
            (-0.5, -0.5, false),
            (0.5, -w.abs(), true),
            (1.5, -0.5, false),
            (2.0, -0.5, false),
            (2.0, 0.5, false),
            (-0.5, 0.5, false),
        ];
        if w > 0.0 {
            pts = pts.iter().map(|&(u, w, t)| (u, -w, t)).rev().collect();
        }
        pts
    };
    let hole = |above: bool, w: f64| {
        let (w0, w1) = if w < 0.0 {
            (-50.0, -20.0)
        } else {
            (20.0, 50.0)
        };
        let plate = if above {
            tilted(1.0, &PLATE)
        } else {
            tilted_below(1.0, &PLATE)
        };
        holed(plate, 1.0, w0, w1)
    };
    let refused = |at: &str| -> [String; 6] { std::array::from_fn(|_| band_refusal(at)) };
    for above in [true, false] {
        for w in [-5.0, 5.0] {
            let at = if w < 0.0 {
                "(500.000000, -0.000492, -0.000087) mm"
            } else {
                "(500.000000, 0.000492, 0.000087) mm"
            };
            assert_eq!(six(&a, &body(above, &run(w)), 1.0), refused(at));
        }
    }
    // The prism above: the layer below the axis is in both bodies.
    let at = "(500.000000, -0.000492, -0.000087) mm";
    assert_eq!(six(&a, &body(true, &apex(-5.0)), 1.0), refused(at));
    let at = "(500.000000, -0.002000, -0.000035) mm";
    assert_eq!(six(&a, &hole(true, -5.0), 1.0), refused(at));
    // The prism below: the layer above the axis is in the prism alone.
    let at = "(500.000000, 0.000492, 0.000087) mm";
    assert_eq!(six(&a, &body(false, &apex(5.0)), 1.0), refused(at));
    let at = "(500.000000, 0.002000, 0.000035) mm";
    assert_eq!(six(&a, &hole(false, 5.0), 1.0), refused(at));
    // The prism above, the layer above the axis in neither body.
    let touch = point_refusal("(500.000000, 0.000492, 0.000087) mm");
    let tip = "(5, 9, 6, 1) V=0.033061309 A=1.687896940";
    assert_eq!(
        six(&a, &body(true, &apex(5.0)), 1.0),
        [
            touch.clone(),
            touch,
            tip.into(),
            tip.into(),
            "(9, 21, 14, 1) V=2.966938691 A=13.076593531".into(),
            "(11, 27, 18, 1) V=0.966938941 A=7.541441234".into(),
        ]
    );
    let touch = edge_refusal("(500.000000, 0.002000, 0.000035) mm");
    let tip = "(5, 9, 6, 1) V=0.003272825 A=1.517569311";
    assert_eq!(
        six(&a, &hole(true, 5.0), 1.0),
        [
            touch.clone(),
            touch,
            tip.into(),
            tip.into(),
            "(9, 21, 14, 1) V=2.996727175 A=13.008613286".into(),
            "(13, 33, 22, 1) V=1.246725675 A=9.491383714".into(),
        ]
    );
    // The prism below, the layer below the axis in the box alone.
    let both = "(16, 42, 28, 1) V=3.517436325 A=16.256765401";
    let common = "(8, 18, 12, 1) V=0.482563925 A=4.361269365";
    assert_eq!(
        six(&a, &body(false, &apex(-5.0)), 1.0),
        [
            both.into(),
            both.into(),
            common.into(),
            common.into(),
            point_refusal("(500.000000, -0.000492, -0.000087) mm"),
            "(12, 30, 20, 1) V=0.517436325 A=6.441742078".into(),
        ]
    );
    let both = "(15, 41, 28, 1) V=3.503272825 A=16.017340775";
    let common = "(11, 27, 18, 1) V=0.746725675 A=6.482656181";
    assert_eq!(
        six(&a, &hole(false, -5.0), 1.0),
        [
            both.into(),
            both.into(),
            common.into(),
            common.into(),
            edge_refusal("(500.000000, -0.004999, -0.000087) mm"),
            "(11, 27, 18, 1) V=0.503272825 A=6.508841778".into(),
        ]
    );
}

/// Golden 39: edges lying on the box's top beside the line where the faces
/// cross, with the face beyond them from the line (found by the audit of
/// the third rule, design note 0012 §20): in golden 32's prism at 10°, a
/// hole from 5 tolerances below the axis, a U-notch and a tab down to 5
/// below it, against a box from x = 0.2 to 0.8 under them; an edge 2
/// tolerances above the axis, the face above it, against a box from −1.25
/// to 1.25; and an edge 3 above it whose face crosses the axis farther on,
/// against a box to x = 1.9. The sweep of each edge toward the line lies
/// outside the face, so there is no strip, and each is answered as on
/// main.
#[test]
fn an_edge_lying_beside_the_line_with_the_face_beyond_it() {
    let eps = k().tolerance().linear;
    let narrow = cuboid(DVec3::new(0.2, -1.0, -1.0), DVec3::new(0.8, 1.0, 0.0));
    let common = "(5, 9, 6, 1) V=0.013224524 A=0.701607182";
    let less = "(9, 21, 14, 1) V=1.186775476 A=7.604188344";
    let full = tilted(10.0, &PLATE);
    let notch = tilted(
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
    let tab = tilted(
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
    for (b, both, rest) in [
        (
            holed(full, 10.0, -5.0, 0.3 / eps),
            "(17, 46, 30, 1) V=2.286775226 A=16.098392220",
            "(13, 34, 22, 1) V=1.086775226 A=9.195811058",
        ),
        (
            notch,
            "(18, 46, 30, 1) V=2.186775226 A=14.898392220",
            "(14, 34, 22, 1) V=0.986775226 A=7.995811058",
        ),
        (
            tab,
            "(16, 40, 26, 1) V=1.936774976 A=13.398390699",
            "(12, 28, 18, 1) V=0.736774976 A=6.495809537",
        ),
    ] {
        assert_eq!(
            six(&narrow, &b, 1.0),
            [both, both, common, common, less, rest].map(String::from)
        );
    }
    let wide = cuboid(DVec3::new(-1.25, -1.0, -1.0), DVec3::new(1.25, 1.0, 0.0));
    let step = tilted(
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
    let touch = edge_refusal("(500.000000, 0.000197, 0.000035) mm");
    let common = "(5, 9, 6, 1) V=0.022040873 A=1.139958749";
    assert_eq!(
        six(&wide, &step, 1.0),
        [
            touch.clone(),
            touch,
            common.into(),
            common.into(),
            "(10, 21, 14, 1) V=4.977959127 A=19.124531940".into(),
            "(9, 20, 13, 1) V=0.790459002 A=6.479020980".into(),
        ]
    );
    let long = cuboid(DVec3::new(-0.25, -1.0, -1.0), DVec3::new(1.9, 1.0, 0.0));
    let far = tilted(
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
    let touch = edge_refusal("(250.000000, 0.000295, 0.000052) mm");
    let common = "(5, 9, 6, 1) V=0.000006465 A=0.019734762";
    assert_eq!(
        six(&long, &far, 1.0),
        [
            touch.clone(),
            touch,
            common.into(),
            common.into(),
            "(9, 21, 14, 1) V=4.299993535 A=16.901439332".into(),
            "(10, 23, 15, 1) V=0.504993385 A=4.523584919".into(),
        ]
    );
}

/// Golden 40: corners of golden 32's prism's face within the tolerance of
/// the box's top beside the line, under it, that are not apexes (design
/// note 0012 §20). A convex corner 5 and 3 tolerances across the line,
/// its edges crossing the line, the tip of a triangle and of a wider face:
/// the face's part across the line tapers to the corner, and B rests on
/// the top by it (the unions touch at a point, A ∩ B is empty). A convex
/// corner 5 and 3 tolerances off the line pointing at it, its edges
/// running away from it under the top, the tip of a triangle and of a
/// spike on a wider face: the face lies beyond the corner from the line,
/// and A ∩ B is the wedge of B under the top, ending at the corner in a
/// tip under the tolerance thick. A reflex corner 3 tolerances off the
/// line, one of its edges crossing the line 10 tolerances short of the
/// corner and the other running away below it (and mirrored): the region
/// between the corner and the line is the end of the face's part across
/// the line, bounded by that edge. Each is answered as on main; taken for
/// an apex (the second without the rule's test that the corner is reflex,
/// the third without its test that both edges run away from the plane,
/// the first without both), each is refused in every operation.
#[test]
fn a_corner_beside_the_line_that_is_no_apex() {
    let a = cuboid(DVec3::new(-0.25, -1.0, -1.0), DVec3::new(1.25, 1.0, 0.0));
    let wide = |w: f64| {
        [
            (0.5, w, true),
            (1.5, 0.2, false),
            (2.0, 0.5, false),
            (-0.5, 0.5, false),
            (-0.5, 0.2, false),
        ]
    };
    let tri = tilted(
        10.0,
        &[(0.0, 0.5, false), (0.5, -5.0, true), (1.0, 0.5, false)],
    );
    for (b, at, rest) in [
        (
            tri,
            "(500.000000, -0.000492, -0.000087) mm",
            "(5, 9, 6, 1) V=0.125000125 A=1.707107635",
        ),
        (
            tilted(10.0, &wide(-5.0)),
            "(500.000000, -0.000492, -0.000087) mm",
            "(7, 15, 10, 1) V=0.437500250 A=4.461352596",
        ),
        (
            tilted(10.0, &wide(-3.0)),
            "(500.000000, -0.000295, -0.000052) mm",
            "(7, 15, 10, 1) V=0.437500150 A=4.461352156",
        ),
    ] {
        let touch = point_refusal(at);
        assert_eq!(
            six(&a, &b, 1.0),
            [
                touch.clone(),
                touch,
                "empty".into(),
                "empty".into(),
                "(6, 12, 8, 1) V=3.000000000 A=13.000000000".into(),
                rest.into(),
            ]
        );
    }
    let spike = |w: f64| {
        [
            (-0.5, -0.5, false),
            (0.0, -0.5, false),
            (0.5, w, true),
            (1.0, -0.5, false),
            (2.0, -0.5, false),
            (2.0, -0.6, false),
            (-0.5, -0.6, false),
        ]
    };
    let triangle = [(0.0, -0.5, false), (0.5, -5.0, true), (1.0, -0.5, false)];
    for (pts, both, common, less, rest) in [
        (
            &triangle[..],
            "(10, 21, 14, 1) V=3.110305967 A=14.052745310",
            "(5, 8, 5, 1) V=0.014693908 A=0.654360618",
            "(10, 20, 13, 1) V=2.985306092 A=13.146647804",
            "(5, 9, 6, 1) V=0.110305967 A=1.560458123",
        ),
        (
            &spike(-5.0)[..],
            "(17, 45, 30, 1) V=3.220758991 A=15.716422813",
            "(9, 20, 13, 1) V=0.029240884 A=1.090683115",
            "(13, 32, 21, 1) V=2.970759116 A=13.239550382",
            "(13, 33, 22, 1) V=0.220758991 A=3.567555546",
        ),
        (
            &spike(-3.0)[..],
            "(17, 45, 30, 1) V=3.220759038 A=15.716422943",
            "(9, 20, 13, 1) V=0.029240887 A=1.090683326",
            "(13, 32, 21, 1) V=2.970759113 A=13.239550396",
            "(13, 33, 22, 1) V=0.220759038 A=3.567555873",
        ),
    ] {
        assert_eq!(
            six(&a, &tilted(10.0, pts), 1.0),
            [both, both, common, common, less, rest].map(String::from)
        );
    }
    let left = [
        (-0.5, 0.15, false),
        (0.0, 0.15, false),
        (0.5, -3.0, true),
        (0.6, -0.5, false),
        (2.0, -0.5, false),
        (2.0, 0.5, false),
        (-0.5, 0.5, false),
    ];
    let right = [
        (-0.5, -0.5, false),
        (0.4, -0.5, false),
        (0.5, -3.0, true),
        (1.0, 0.15, false),
        (2.0, 0.15, false),
        (2.0, 0.5, false),
        (-0.5, 0.5, false),
    ];
    let common = "(5, 9, 6, 1) V=0.015061263 A=0.807223617";
    let less = "(9, 21, 14, 1) V=2.984938737 A=13.052343435";
    for (pts, both, rest) in [
        (
            left,
            "(15, 38, 25, 1) V=3.916188782 A=19.308735098",
            "(11, 26, 17, 1) V=0.916188782 A=7.063615279",
        ),
        (
            right,
            "(15, 38, 25, 1) V=3.753688782 A=18.658735098",
            "(11, 26, 17, 1) V=0.753688782 A=6.413615279",
        ),
    ] {
        assert_eq!(
            six(&a, &tilted(10.0, &pts), 1.0),
            [both, both, common, common, less, rest].map(String::from)
        );
    }
}

/// A body extruded 3 up from z = −2 by a circle in the xy plane.
fn rod(c: DVec2, r: f64) -> Body {
    let frame = Frame {
        origin: DVec3::new(0.0, 0.0, -2.0),
        ..Frame::WORLD
    };
    k().extrude(
        &Profile::circle(frame, c, r),
        &ExtrudeSpec {
            distance: 3.0,
            symmetric: false,
        },
    )
    .expect("a rod")
}

/// Golden 41: golden 32's prism, its edge 5 tolerances off the axis,
/// against faces with curved edges (design note 0012 §20). A round bar
/// whose top holds the strip, as golden 32's box does: every operation is
/// refused. A slab whose round hole lies under the strip, which nothing
/// then holds: answered as the contact, A ∩ B the two wedges. (Placed by
/// their outlines in their planes, which take an arc by its chord, the
/// bar's top held no point and the hole every point the slab's top did:
/// the bar's A ∩ B was the wedges joined by the fin, and the hole was
/// refused; found by review.) The prism's face turning away from the axis
/// past the edge's end along a quarter circle instead of a line, which
/// leaves the strip: refused as golden 32's box is. A round hole in the
/// box's top whose rim passes 0.9 tolerances from the strip's far side:
/// with the box first, its strip check refuses it as not supported; with
/// the prism first the rim is refused first as a tangent contact; 1.5
/// tolerances off, it is refused as golden 32's box is. A round hole of
/// the prism's face whose rim dips 2 tolerances across the axis into the
/// strip, over a gap in the slab's top under the whole strip: refused as
/// not supported in every operation; dipping half a tolerance across, it
/// stays within the band of the line, outside the strip, and is answered
/// as the plain prism is over the gap.
#[test]
fn a_face_with_curved_edges_beside_the_line_where_the_faces_cross() {
    let eps = k().tolerance().linear;
    let b = beside(5.0);
    let at = band_refusal("(500.000000, -0.000492, -0.000087) mm");
    let bar = {
        let frame = Frame {
            origin: DVec3::new(0.0, 0.0, -1.0),
            ..Frame::WORLD
        };
        k().extrude(
            &Profile::circle(frame, DVec2::new(0.5, 0.0), 0.8),
            &ExtrudeSpec {
                distance: 1.0,
                symmetric: false,
            },
        )
        .expect("a bar")
    };
    assert_eq!(six(&bar, &b, 1.0), [at.as_str(); 6]);
    let slab = cuboid(DVec3::new(-1.0, -1.0, -1.0), DVec3::new(3.0, 1.0, 0.0));
    let holed = op(&slab, &rod(DVec2::new(0.5, 0.1), 0.75), BooleanOp::Subtract);
    let both = "(17, 53, 36, 1) V=7.090039217 A=35.326561698";
    let wedges = "(11, 23, 16, 2) V=0.017815290 A=1.058643475";
    assert_eq!(
        six_meshed(&holed, &b, 1.0, 1e-7),
        [
            both,
            both,
            wedges,
            wedges,
            "(15, 41, 28, 1) V=6.215038842 A=29.305379246",
            "(13, 35, 24, 1) V=0.857185085 A=7.079825927",
        ]
        .map(String::from)
    );
    // The face's boundary past the edge's end at x = 1 a quarter circle.
    let w = 5.0 * eps;
    let line = |a: (f64, f64), b: (f64, f64)| ProfileCurve::Line {
        start: DVec2::new(a.0, a.1),
        end: DVec2::new(b.0, b.1),
    };
    let curves = vec![
        line((-0.5, -0.5), (0.0, -w)),
        line((0.0, -w), (1.0, -w)),
        ProfileCurve::Arc {
            center: DVec2::new(1.5, -w),
            radius: 0.5,
            start_angle: PI,
            sweep_angle: FRAC_PI_2,
        },
        line((1.5, -0.5 - w), (2.0, -0.5 - w)),
        line((2.0, -0.5 - w), (2.0, 0.5)),
        line((2.0, 0.5), (-0.5, 0.5)),
        line((-0.5, 0.5), (-0.5, -0.5)),
    ];
    let rounded = k()
        .extrude(
            &Profile {
                plane: tilted_frame(10.0),
                curves,
            },
            &ExtrudeSpec {
                distance: 0.5,
                symmetric: false,
            },
        )
        .expect("a prism");
    let boxed = cuboid(DVec3::new(-0.25, -1.0, -1.0), DVec3::new(1.25, 1.0, 0.0));
    assert_eq!(six(&boxed, &rounded, 1.0), [at.as_str(); 6]);
    // The box with a round hole whose rim passes `d` tolerances beyond the
    // strip's far side, at x = 0.5.
    let far = w * 10f64.to_radians().cos();
    let rim = |d: f64| {
        let r = 0.2;
        op(
            &boxed,
            &rod(DVec2::new(0.5, -far - d * eps - r), r),
            BooleanOp::Subtract,
        )
    };
    let curved = |at: &str| {
        format!(
            "unsupported: booleans where a face of one body runs within the tolerance of the \
             other's surface beside a curved edge (design note 0012 §20) near {at}"
        )
    };
    let (beside_rim, tangent) = (
        curved("(1000.000000, -0.000492, -0.000087) mm"),
        "unsupported: booleans where an edge of one body touches a face of the other without \
         crossing it (a tangent contact; design note 0011 §4.4) near (500.000000, -0.000492, \
         -0.000087) mm",
    );
    assert_eq!(
        six(&rim(0.9), &b, 1.0),
        [
            &beside_rim,
            tangent,
            &beside_rim,
            tangent,
            &beside_rim,
            tangent
        ]
        .map(String::from)
    );
    assert_eq!(six(&rim(1.5), &b, 1.0), [at.as_str(); 6]);
    // The prism's face with a round hole whose rim dips `dip` tolerances
    // across the axis, over a gap in a slab's top under the whole strip.
    let gap = cuboid(DVec3::new(-0.2, -0.6, -2.0), DVec3::new(1.2, 0.6, 1.0));
    let gapped = op(&slab, &gap, BooleanOp::Subtract);
    let pierced = |dip: f64| {
        let r = 0.1;
        let frame = tilted_frame(10.0);
        let pts = [
            (-0.5, -0.5),
            (0.0, -w),
            (1.0, -w),
            (1.5, -0.5),
            (2.0, -0.5),
            (2.0, 0.5),
            (-0.5, 0.5),
        ]
        .map(|(x, y)| DVec2::new(x, y));
        k().extrude(
            &Profile::polygon(frame, &pts).with(Profile::circle(
                frame,
                DVec2::new(0.5, r - dip * eps),
                r,
            )),
            &ExtrudeSpec {
                distance: 0.5,
                symmetric: false,
            },
        )
        .expect("a pierced prism")
    };
    let (both, wedges, less) = (
        "(19, 57, 38, 1) V=7.177103182 A=35.972428476",
        "(11, 23, 16, 2) V=0.017897193 A=1.074679434",
        "(17, 45, 30, 1) V=6.302102807 A=29.966380585",
    );
    assert_eq!(
        six(&gapped, &b, 1.0),
        [
            both,
            both,
            wedges,
            wedges,
            less,
            "(13, 35, 24, 1) V=0.857103182 A=7.080727325",
        ]
        .map(String::from)
    );
    assert_eq!(six(&gapped, &pierced(2.0), 1.0), [beside_rim.as_str(); 6]);
    let both = "(21, 63, 42, 1) V=7.161395218 A=36.223755888";
    assert_eq!(
        six(&gapped, &pierced(0.5), 1.0),
        [
            both,
            both,
            wedges,
            wedges,
            less,
            "(15, 41, 28, 1) V=0.841395218 A=7.332054737",
        ]
        .map(String::from)
    );
}

/// A prism tilted `t` degrees as [`tilted`] whose face has a V-notch from
/// its far side, its tip at `x` along the axis and `v` tolerances across
/// it, its mouth 0.2 wide at w = 0.5.
fn notch_at(t: f64, x: f64, v: f64) -> Body {
    tilted(
        t,
        &[
            (-0.5, -0.5, false),
            (2.0, -0.5, false),
            (2.0, 0.5, false),
            (x + 0.1, 0.5, false),
            (x, v, true),
            (x - 0.1, 0.5, false),
            (-0.5, 0.5, false),
        ],
    )
}

/// Golden 42: golden 32's prism with a V-notch from its face's far side
/// whose tip V lies across the line where the face crosses the box's top,
/// within the tolerance of the top but more than a tolerance from the line
/// (design note 0012 §21). V is on the top by its sign, so its edges get
/// no crossings and the face's interval on the line passes over the
/// notch's mouth; but the meeting's split does not take V, which is too
/// far from it, and the two disagree. On main, with the box first, every
/// operation failed (a contact's pieces disagreeing on their side), and
/// with the prism first each was refused as edges passing within the
/// tolerance of each other (the notch mirrored and the prism below its
/// face were refused so in every operation, in both orders). Each is now
/// refused at V in every operation:
/// at 10° with V 2, 3 and 5 tolerances across, at 1° with V 5, 30 and 50;
/// V the tip of a triangular hole of the face; the notch's walls leaning,
/// one meeting the line at 120°; V 2 and 4 tolerances inside the box's end
/// x = 1.25; and two notches (a W), over the box's top or over the two
/// faces of a box grooved between them, reported at the tip with the least
/// coordinates in both orders. The notch mirrored, from the face's near
/// side with V above the top, and the prism below its face, V 3
/// tolerances across at 10°, are no isolated touch of the top at V (the
/// prism, and its edge at V out of the face, lie on the other side of the
/// top from V's neighbours), so rule 3 leaves them to the split, which
/// refuses V's edges' crossings of the line as passing within the
/// tolerance of each other in every operation, in both orders, as on main.
#[test]
fn a_notch_tip_across_the_line_where_the_faces_cross() {
    let eps = k().tolerance().linear;
    let a = cuboid(DVec3::new(-0.25, -1.0, -1.0), DVec3::new(1.25, 1.0, 0.0));
    let refused = |b: &Body, at: &str| {
        assert_eq!(six(&a, b, 1.0), [(); 6].map(|_| band_refusal(at)), "{at}");
    };
    for (t, v, at) in [
        (10.0, -2.0, "(500.000000, -0.000197, -0.000035) mm"),
        (10.0, -3.0, "(500.000000, -0.000295, -0.000052) mm"),
        (10.0, -5.0, "(500.000000, -0.000492, -0.000087) mm"),
        (1.0, -5.0, "(500.000000, -0.000500, -0.000009) mm"),
        (1.0, -30.0, "(500.000000, -0.003000, -0.000052) mm"),
        (1.0, -50.0, "(500.000000, -0.004999, -0.000087) mm"),
    ] {
        refused(&notch_at(t, 0.5, v), at);
    }
    let at = "(500.000000, -0.000295, -0.000052) mm";
    let mirrored = tilted(
        10.0,
        &[
            (-0.5, -0.5, false),
            (0.4, -0.5, false),
            (0.5, 3.0, true),
            (0.6, -0.5, false),
            (2.0, -0.5, false),
            (2.0, 0.5, false),
            (-0.5, 0.5, false),
        ],
    );
    let split = |b: &Body, at: &str| {
        let m = format!(
            "refused: edges where the bodies meet pass within the tolerance of each other near {at}"
        );
        assert_eq!(six(&a, b, 1.0), [(); 6].map(|_| m.clone()), "{at}");
    };
    split(&mirrored, "(500.000000, 0.000295, 0.000052) mm");
    let plate = tilted(
        10.0,
        &[
            (-0.5, -0.5, false),
            (2.0, -0.5, false),
            (2.0, 0.5, false),
            (-0.5, 0.5, false),
        ],
    );
    let f = tilted_frame(10.0);
    let hole = prism(
        Frame {
            origin: f.origin - f.z * 0.25,
            ..f
        },
        &[
            DVec2::new(0.4, 0.3),
            DVec2::new(0.5, -3.0 * eps),
            DVec2::new(0.6, 0.3),
        ],
        1.0,
    );
    refused(&op(&plate, &hole, BooleanOp::Subtract), at);
    let leaning = tilted(
        10.0,
        &[
            (-0.5, -0.5, false),
            (2.0, -0.5, false),
            (2.0, 0.5, false),
            (0.5 - 0.5 / 3f64.sqrt(), 0.5, false),
            (0.5, -3.0, true),
            (-0.2, 0.5, false),
            (-0.5, 0.5, false),
        ],
    );
    refused(&leaning, at);
    let below = [
        (-0.5, -0.5, false),
        (2.0, -0.5, false),
        (2.0, 0.5, false),
        (0.6, 0.5, false),
        (0.5, -3.0, true),
        (0.4, 0.5, false),
        (-0.5, 0.5, false),
    ];
    split(&tilted_below(10.0, &below), at);
    refused(
        &notch_at(10.0, 1.25 - 2.0 * eps, -3.0),
        "(1249.999800, -0.000295, -0.000052) mm",
    );
    refused(
        &notch_at(10.0, 1.25 - 4.0 * eps, -3.0),
        "(1249.999600, -0.000295, -0.000052) mm",
    );
    let w = tilted(
        10.0,
        &[
            (-0.5, -0.5, false),
            (2.0, -0.5, false),
            (2.0, 0.5, false),
            (0.8, 0.5, false),
            (0.7, -3.0, true),
            (0.5, 0.3, false),
            (0.3, -3.0, true),
            (0.2, 0.5, false),
            (-0.5, 0.5, false),
        ],
    );
    let at = band_refusal("(300.000000, -0.000295, -0.000052) mm");
    let slot = cuboid(DVec3::new(0.45, -2.0, -0.5), DVec3::new(0.55, 2.0, 1.0));
    for a in [a.clone(), op(&a, &slot, BooleanOp::Subtract)] {
        assert_eq!(six(&a, &w, 1.0), [at.as_str(); 6]);
    }
}

/// Golden 43: golden 42's notch where the rule leaves main's outcomes
/// (design note 0012 §21). With V half a tolerance and one tolerance
/// across the line, at 10° and at 1°, the meeting's split takes V and
/// bends through it: answered. With V 5 tolerances beyond the box's end,
/// the top does not hold V's foot: answered. With V half a tolerance
/// inside or beyond the end, V is a point of the top's boundary, and the
/// notch's edge is refused as leaving the top at a point the tolerance
/// cannot place, as on main. With the box's top ending along the line
/// between it and V's foot, at 1° with V 50 tolerances across, the top's
/// edge lies on the face and the top does not hold the foot: the contact's
/// outcomes. With the prism below its face and V half a tolerance across:
/// answered, the differences touching at a point.
#[test]
fn a_notch_tip_across_the_line_that_is_answered() {
    let eps = k().tolerance().linear;
    let a = cuboid(DVec3::new(-0.25, -1.0, -1.0), DVec3::new(1.25, 1.0, 0.0));
    let common = "(5, 9, 6, 1) V=0.033061309 A=1.687896940";
    let less = "(9, 21, 14, 1) V=2.966938691 A=13.076593531";
    for (b, union, rest) in [
        (
            notch_at(10.0, 0.5, -0.5),
            "(17, 43, 28, 1) V=4.191938689 A=20.122005125",
            "(13, 31, 20, 1) V=1.191938689 A=8.733308460",
        ),
        (
            notch_at(10.0, 0.5, -1.0),
            "(17, 43, 28, 1) V=4.191938686 A=20.122005239",
            "(13, 31, 20, 1) V=1.191938686 A=8.733308499",
        ),
        (
            notch_at(10.0, 1.25 + 5.0 * eps, -3.0),
            "(16, 42, 28, 1) V=4.191938676 A=20.122005246",
            "(12, 30, 20, 1) V=1.191938676 A=8.733308654",
        ),
    ] {
        assert_eq!(
            six(&a, &b, 1.0),
            [union, union, common, common, less, rest].map(String::from)
        );
    }
    let common = "(5, 9, 6, 1) V=0.003272825 A=1.517569311";
    let less = "(9, 21, 14, 1) V=2.996727175 A=13.008613286";
    for (v, union, rest) in [
        (
            -0.5,
            "(17, 43, 28, 1) V=4.221727173 A=20.292332754",
            "(13, 31, 20, 1) V=1.221727173 A=8.801288704",
        ),
        (
            -1.0,
            "(17, 43, 28, 1) V=4.221727170 A=20.292332868",
            "(13, 31, 20, 1) V=1.221727170 A=8.801288743",
        ),
    ] {
        assert_eq!(
            six(&a, &notch_at(1.0, 0.5, v), 1.0),
            [union, union, common, common, less, rest].map(String::from)
        );
    }
    for (x, at) in [
        (1.25 - 0.5 * eps, "(1250.000001, -0.000026, -0.000001) mm"),
        (1.25 + 0.5 * eps, "(1249.999999, -0.000026, -0.000001) mm"),
    ] {
        let leaves = format!(
            "refused: an edge of one body lies on a face of the other and leaves it at a point \
             the tolerance cannot place near {at}"
        );
        assert_eq!(
            six(&a, &notch_at(10.0, x, -3.0), 1.0),
            [(); 6].map(|_| leaves.clone())
        );
    }
    for (y, at, less) in [
        (
            -10.0,
            "(124.999600, -0.001000, -0.000004) mm",
            "(6, 12, 8, 1) V=1.500001500 A=8.000005000",
        ),
        (
            -25.0,
            "(124.999750, -0.002500, -0.000011) mm",
            "(6, 12, 8, 1) V=1.500003750 A=8.000012500",
        ),
        (
            -40.0,
            "(124.999900, -0.004000, -0.000017) mm",
            "(6, 12, 8, 1) V=1.500006000 A=8.000020000",
        ),
    ] {
        let short = cuboid(DVec3::new(-0.25, y * eps, -1.0), DVec3::new(1.25, 1.0, 0.0));
        let touch = edge_refusal(at);
        assert_eq!(
            six(&short, &notch_at(1.0, 0.5, -50.0), 1.0),
            [
                touch.clone(),
                touch,
                "empty".into(),
                "empty".into(),
                less.into(),
                "(9, 21, 14, 1) V=1.224999750 A=8.809905854".into(),
            ]
        );
    }
    let below = [
        (-0.5, -0.5, false),
        (2.0, -0.5, false),
        (2.0, 0.5, false),
        (0.6, 0.5, false),
        (0.5, -0.5, true),
        (0.4, 0.5, false),
        (-0.5, 0.5, false),
    ];
    let union = "(19, 48, 31, 1) V=3.530122526 A=16.192850395";
    let common = "(11, 25, 16, 1) V=0.694877472 A=5.617051595";
    let touch = point_refusal("(500.000000, -0.000049, -0.000009) mm");
    assert_eq!(
        six(&a, &tilted_below(10.0, &below), 1.0),
        [union, union, common, common, &touch, &touch]
    );
}

/// A prism tilted `t` degrees as [`tilted`] (`below`: [`tilted_below`])
/// whose face has a spike 0.2 wide at its root, its tip at x = 0.5 and `d`
/// tolerances across the axis: rising from a strip w ≤ −0.2 (`up`), or
/// falling from a strip w ≥ 0.2.
fn spiked(t: f64, d: f64, up: bool, below: bool) -> Body {
    let pts = if up {
        [
            (-0.5, -0.5, false),
            (2.0, -0.5, false),
            (2.0, -0.2, false),
            (0.6, -0.2, false),
            (0.5, d, true),
            (0.4, -0.2, false),
            (-0.5, -0.2, false),
        ]
    } else {
        [
            (-0.5, 0.2, false),
            (0.4, 0.2, false),
            (0.5, d, true),
            (0.6, 0.2, false),
            (2.0, 0.2, false),
            (2.0, 0.5, false),
            (-0.5, 0.5, false),
        ]
    };
    if below {
        tilted_below(t, &pts)
    } else {
        tilted(t, &pts)
    }
}

/// Golden 44: the tip V of a spike of the prism's face across the line
/// where the face crosses the box's top, within the tolerance of the top
/// (design note 0012 §21). The face's interval on the line at the spike is
/// a single point, so no meeting is made there: the flap of the face
/// across the line, between V's edges, stays in the face, and so does the
/// piece of the top under it. Where the layer between them lies outside
/// the prism (V above the top with the prism above its face, or V below
/// the top with the prism below its face), A ∩ B or B − A kept both, and
/// their faces crossed along the line, the result validating; now, with V
/// 3 tolerances across at 10° and 30 at 1°, where the flap holds points
/// more than a tolerance from both of V's edges over the top, every
/// operation is refused. With V 1.5 tolerances across at 10°, every point
/// of the flap is within a tolerance of V's edges: answered. Where the
/// layer lies in the prism, with V above the top and the prism below its
/// face, the operation that would keep both touches at V and is refused,
/// and the others are answered, as on main (with V below the top and the
/// prism above: golden 40).
#[test]
fn a_spike_tip_across_the_line_where_the_faces_cross() {
    let a = cuboid(DVec3::new(-0.25, -1.0, -1.0), DVec3::new(1.25, 1.0, 0.0));
    let at = band_refusal("(500.000000, 0.000000, 0.000000) mm");
    for (t, d) in [(10.0, 3.0), (1.0, 30.0)] {
        assert_eq!(six(&a, &spiked(t, d, true, false), 1.0), [at.as_str(); 6]);
        assert_eq!(six(&a, &spiked(t, -d, false, true), 1.0), [at.as_str(); 6]);
    }
    let union = "(17, 45, 30, 1) V=3.356758303 A=16.293351921";
    let common = "(9, 20, 13, 1) V=0.028241705 A=1.170255041";
    assert_eq!(
        six(&a, &spiked(10.0, 1.5, true, false), 1.0),
        [
            union,
            union,
            common,
            common,
            "(13, 32, 21, 1) V=2.971758295 A=13.141696664",
            "(13, 33, 22, 1) V=0.356758303 A=4.321910298",
        ]
    );
    let union = "(17, 44, 29, 1) V=3.178241705 A=14.815753996";
    let common = "(9, 21, 14, 1) V=0.206758303 A=2.647852966";
    assert_eq!(
        six(&a, &spiked(10.0, -1.5, false, true), 1.0),
        [
            union,
            union,
            common,
            common,
            "(13, 33, 22, 1) V=2.793241697 A=14.167409253",
            "(13, 32, 21, 1) V=0.178241705 A=3.296197709",
        ]
    );
    let union = "(16, 36, 24, 1) V=3.150000000 A=14.600000000";
    let common = "(9, 21, 14, 1) V=0.235000015 A=2.863607126";
    assert_eq!(
        six(&a, &spiked(10.0, 3.0, true, true), 1.0),
        [
            union.into(),
            union.into(),
            common.into(),
            common.into(),
            point_refusal("(500.000000, 0.000295, 0.000052) mm"),
            "(12, 24, 16, 2) V=0.150000000 A=2.200000000".into(),
        ]
    );
}

/// Golden 45: a V-notch from the far side of golden 32's prism's face,
/// tilted 33°, its tip 1.7 tolerances under the box's top and 0.8
/// tolerances beyond the box's side x = −0.25 (or x = 1.25), its edges at
/// 78.7° to the line (found by the notch-tip probes, design note 0012
/// §21). The notch's root, the prism's edge from the tip square to its
/// face, lies on the box's side as the whole interval of a wall of the
/// notch on the line where their planes meet, but no part of the side's
/// interval runs along it: the wall leans 11.3° from the side and lies within
/// the tolerance of it over a triangle between the root, the box's top
/// edge and the notch's edge, which neither face's imprint cut out, and
/// every operation failed in both orders as kept pieces that do not close.
/// A whole interval is the meeting only where a stretch runs along it, so
/// the root's strip is tested, the side holds it, and every operation is
/// refused. With the tip 1.2 tolerances beyond the side, the root no
/// longer lies on it: answered.
#[test]
fn a_notch_wall_within_the_tolerance_of_a_side() {
    let eps = k().tolerance().linear;
    let a = cuboid(DVec3::new(-0.25, -1.0, -1.0), DVec3::new(1.25, 1.0, 0.0));
    let tip = -1.7 / 33f64.to_radians().sin();
    let notch = |x: f64| {
        tilted(
            33.0,
            &[
                (-1.2, -0.5, false),
                (2.2, -0.5, false),
                (2.2, 0.5, false),
                (x + 0.1, 0.5, false),
                (x, tip, true),
                (x - 0.1, 0.5, false),
                (-1.2, 0.5, false),
            ],
        )
    };
    for (x, at) in [
        (-0.25 - 0.8 * eps, "(-250.000080, -0.000284, -0.000135) mm"),
        (1.25 + 0.8 * eps, "(1250.000080, -0.000284, -0.000135) mm"),
    ] {
        assert_eq!(
            six(&a, &notch(x), 1.0),
            [(); 6].map(|_| band_refusal(at)),
            "{x}"
        );
    }
    let union = "(17, 43, 28, 1) V=4.553236061 A=22.216222097";
    let common = "(5, 9, 6, 1) V=0.121763924 A=2.293680094";
    for x in [-0.25 - 1.2 * eps, 1.25 + 1.2 * eps] {
        assert_eq!(
            six(&a, &notch(x), 1.0),
            [
                union,
                union,
                common,
                common,
                "(9, 21, 14, 1) V=2.878236076 A=13.180431349",
                "(13, 31, 20, 1) V=1.553236061 A=11.329470841",
            ],
            "{x}"
        );
    }
}

/// A prism from z = −1 up to the plane z = 0, its top the polygon `pts`.
fn slab(pts: &[DVec2]) -> Body {
    let frame = Frame {
        origin: DVec3::new(0.0, 0.0, -1.0),
        ..Frame::WORLD
    };
    prism(frame, pts, 1.0)
}

/// Golden 46: golden 44's rising spike at 2°, its tip V 28 tolerances
/// across the line and 0.98 above the top, over a slab whose top's edge
/// runs at 45° to the line through V's foot, the top holding all of the
/// flap on its side of that edge (found by review, design note 0012 §21).
/// V's cluster lies on the top's edge, which once let the whole flap off
/// as a contact at the top's boundary: A ∩ B validated with the flap and
/// the top crossing 12.5 tolerances inside both, as on main. Only the
/// flap's points within the tolerance of the place of a member of V's
/// cluster are the contact's (the members' discs; the first build took
/// those within three tolerances of the cluster's place); the top holds
/// the rest, and every operation is refused.
#[test]
fn a_spike_tip_on_the_edge_of_a_face_holding_its_flap() {
    let eps = k().tolerance().linear;
    let (t, d) = (2.0f64, 28.0);
    let y0 = d * eps * t.to_radians().cos();
    let x = |y: f64| 0.5 + d * eps - y * (d * eps) / y0;
    let a = slab(&[
        DVec2::new(-3.0, -1.0),
        DVec2::new(x(-1.0), -1.0),
        DVec2::new(x(1.0), 1.0),
        DVec2::new(-3.0, 1.0),
    ]);
    let at = band_refusal("(500.000000, 0.000000, 0.000000) mm");
    assert_eq!(
        six(&a, &spiked(t, d, true, false), 1.0),
        [(); 6].map(|_| at.clone())
    );
}

/// Golden 47: a run that is its face's whole interval on the line where
/// the face crosses a slab's top, and passes a corner of the top (found by
/// review, design note 0012 §21). B, tilted `t` degrees, has an edge R at
/// `h` tolerances across the line from x = 0 to 1, within the tolerance of
/// the top, its face crossing the top only beyond x = 5. The top lies 20
/// tolerances beyond R, but for a tongue entering the strip between R and
/// the line from x > 1, its tip `w2` tolerances from the line at x = 0.7,
/// and a bump whose apex lies on R at x = 0.3. No meeting covers R, and the
/// tongue holds R's strip more than a tolerance inside: every operation
/// is refused, with the bump as without it. Main, and this rule while it
/// let a run passing a corner of the top off whole, failed A ∩ B, B ∩ A
/// and B − A as kept pieces that do not close.
#[test]
fn a_run_no_meeting_covers_passing_a_corner_of_the_face_holding_it() {
    let eps = k().tolerance().linear;
    for (t, h, w2, at) in [
        (1.0f64, 55.0, 30.0, "(1000.000000, -0.004875, -0.000085) mm"),
        (0.5, 110.0, 40.0, "(1000.000000, -0.007750, -0.000068) mm"),
    ] {
        let yr = -h * eps * t.to_radians().cos();
        let ym = yr - 20.0 * eps;
        let b = tilted(
            t,
            &[
                (0.0, -h, true),
                (1.0, -h, true),
                (1.0, 0.3, false),
                (5.0, 0.3, false),
                (5.0, -0.5, false),
                (6.0, -0.5, false),
                (6.0, 0.5, false),
                (0.0, 0.5, false),
            ],
        );
        for bump in [false, true] {
            let mut top = vec![
                DVec2::new(-0.5, -1.0),
                DVec2::new(4.0, -1.0),
                DVec2::new(4.0, 0.5),
                DVec2::new(3.5, 0.5),
                DVec2::new(3.5, ym),
                DVec2::new(1.6, ym),
                DVec2::new(0.7, -w2 * eps),
                DVec2::new(1.3, ym),
            ];
            if bump {
                top.extend([
                    DVec2::new(0.301, ym),
                    DVec2::new(0.3, yr),
                    DVec2::new(0.299, ym),
                ]);
            }
            top.push(DVec2::new(-0.5, ym));
            let refused = band_refusal(at);
            assert_eq!(
                six(&slab(&top), &b, 1.0),
                [(); 6].map(|_| refused.clone()),
                "{t} {bump}"
            );
        }
    }
}

/// [`tilted`]'s prism from a profile in metres (u along x, w up the
/// plane), above its face or, `below`, under it ([`tilted_below`]).
fn on_plane(t: f64, pts: &[DVec2], below: bool) -> Body {
    let p: Vec<(f64, f64, bool)> = pts.iter().map(|q| (q.x, q.y, false)).collect();
    if below {
        tilted_below(t, &p)
    } else {
        tilted(t, &p)
    }
}

/// The profile of a face from u = −`l` to `r` and w = −0.5 to 0.5 with a
/// V-notch from its far side, its tip at (`x`, `w`) and its mouth from
/// x − `b` to x + `a` at w = 0.5; `near`: from its near side, its mouth at
/// w = −0.5.
fn notch_profile((l, r): (f64, f64), x: f64, w: f64, (a, b): (f64, f64), near: bool) -> Vec<DVec2> {
    let p = DVec2::new;
    if near {
        vec![
            p(-l, -0.5),
            p(x - b, -0.5),
            p(x, w),
            p(x + a, -0.5),
            p(r, -0.5),
            p(r, 0.5),
            p(-l, 0.5),
        ]
    } else {
        vec![
            p(-l, -0.5),
            p(r, -0.5),
            p(r, 0.5),
            p(x + a, 0.5),
            p(x, w),
            p(x - b, 0.5),
            p(-l, 0.5),
        ]
    }
}

/// The profile of a face from u = −`l` to `r` with a spike, its tip at
/// (`x`, `w`) and its root from x − `b` to x + `a`: rising from the strip
/// from w = −0.5 to −`root` (`up`), or falling from the strip from w =
/// `root` to 0.5.
fn spiked_profile(
    (l, r): (f64, f64),
    x: f64,
    w: f64,
    (a, b): (f64, f64),
    root: f64,
    up: bool,
) -> Vec<DVec2> {
    let p = DVec2::new;
    if up {
        vec![
            p(-l, -0.5),
            p(r, -0.5),
            p(r, -root),
            p(x + a, -root),
            p(x, w),
            p(x - b, -root),
            p(-l, -root),
        ]
    } else {
        vec![
            p(-l, root),
            p(x - b, root),
            p(x, w),
            p(x + a, root),
            p(r, root),
            p(r, 0.5),
            p(-l, 0.5),
        ]
    }
}

/// The convex polygon `poly` less what lies beyond the line through `o`
/// square to `n`, keeping n · (q − o) ≤ 0 (Sutherland and Hodgman 1974).
fn clip(poly: &[DVec2], o: DVec2, n: DVec2) -> Vec<DVec2> {
    let f = |q: DVec2| n.dot(q - o);
    let mut out = Vec::new();
    for i in 0..poly.len() {
        let (a, b) = (poly[i], poly[(i + 1) % poly.len()]);
        let (fa, fb) = (f(a), f(b));
        if fa <= 0.0 {
            out.push(a);
        }
        if (fa < 0.0 && fb > 0.0) || (fa > 0.0 && fb < 0.0) {
            out.push(a + (b - a) * (fa / (fa - fb)));
        }
    }
    out
}

/// The box's top: x from −0.25 to 1.25 and y from −1 to 1.
fn top_rect() -> Vec<DVec2> {
    [(-0.25, -1.0), (1.25, -1.0), (1.25, 1.0), (-0.25, 1.0)]
        .map(|(x, y)| DVec2::new(x, y))
        .to_vec()
}

/// `pts` turning counter-clockwise.
fn ccw(mut pts: Vec<DVec2>) -> Vec<DVec2> {
    let n = pts.len();
    if (0..n)
        .map(|i| pts[i].perp_dot(pts[(i + 1) % n]))
        .sum::<f64>()
        < 0.0
    {
        pts.reverse();
    }
    pts
}

/// The box from z = −1 to 0 with x from −0.25 to 1.25 and y from −1 to 1.
fn the_box() -> Body {
    cuboid(DVec3::new(-0.25, -1.0, -1.0), DVec3::new(1.25, 1.0, 0.0))
}

/// `a` less a column from z = −2 to 1 whose section is the polygon `hole`.
fn pierced_by(a: &Body, hole: Vec<DVec2>) -> Body {
    let frame = Frame {
        origin: DVec3::new(0.0, 0.0, -2.0),
        ..Frame::WORLD
    };
    op(a, &prism(frame, &ccw(hole), 3.0), BooleanOp::Subtract)
}

/// An oblique prism (the soundness review's): the profile `pts` of a face
/// on golden 32's plane at `t` degrees, the prism above its face or,
/// `below`, under it, swept along D, which leans `alpha` degrees from the
/// face's inward normal toward the direction `phi` degrees from x in the
/// face's plane, and cut to a slab 0.43 thick behind the face: its edge
/// out of the face at each corner runs along D.
fn oblique(t: f64, pts: &[DVec2], below: bool, alpha: f64, phi: f64) -> Body {
    let thick = 0.43;
    let tf = tilted_frame(t);
    let inward = if below { -tf.z } else { tf.z };
    let (al, ph) = (alpha.to_radians(), phi.to_radians());
    let d = (inward * al.cos() + (tf.x * ph.cos() + tf.y * ph.sin()) * al.sin()).normalize();
    let slab = k()
        .make_box(
            &Frame {
                origin: tf.x * -6.0 + tf.y * -6.0 + if below { tf.z * -thick } else { DVec3::ZERO },
                ..tf
            },
            DVec3::new(12.0, 12.0, thick),
        )
        .expect("a slab");
    let len = thick / al.cos() + 1.0;
    let seed = if d.cross(DVec3::X).length() > 0.3 {
        DVec3::X
    } else {
        DVec3::Y
    };
    let fx = seed.reject_from_normalized(d).normalize();
    let fy = d.cross(fx);
    let o = -d * len;
    let q: Vec<DVec2> = (pts.iter())
        .map(|p| {
            let x = tf.x * p.x + tf.y * p.y;
            let x = x - d * (x - o).dot(d);
            DVec2::new((x - o).dot(fx), (x - o).dot(fy))
        })
        .collect();
    let frame = Frame {
        origin: o,
        x: fx,
        y: fy,
        z: d,
    };
    op(
        &prism(frame, &ccw(q), 2.0 * len),
        &slab,
        BooleanOp::Intersect,
    )
}

/// The refusal of edges where the bodies meet passing within the tolerance
/// of each other, near `at`.
fn crossing_refusal(at: &str) -> String {
    format!(
        "refused: edges where the bodies meet pass within the tolerance of each other near {at}"
    )
}

/// The refusal of a point where the bodies meet within the tolerance of an
/// edge it is not on, near `at`.
fn near_edge_refusal(at: &str) -> String {
    format!(
        "refused: a point where the bodies meet lies within the tolerance of an edge it is not \
         on, near {at}"
    )
}

/// Six times `s`.
fn all(s: String) -> [String; 6] {
    [(); 6].map(|_| s.clone())
}

/// Asserts the six outcomes ([`six`]) of the case `what`.
fn outcomes(what: &str, a: &Body, b: &Body, want: [String; 6]) {
    assert_eq!(six(a, b, 1.0), want, "{what}");
}

/// Answered outcomes: A ∪ B and B ∪ A as `union`, A ∩ B and B ∩ A as
/// `common`, then A − B and B − A.
fn answered(union: &str, common: &str, a_less_b: String, b_less_a: &str) -> [String; 6] {
    [
        union.into(),
        union.into(),
        common.into(),
        common.into(),
        a_less_b,
        b_less_a.into(),
    ]
}

/// The outcomes of a touch the result would keep as a non-manifold point
/// or edge (`refusal`): A ∪ B and B ∪ A refused, A ∩ B and B ∩ A empty,
/// then A − B and B − A.
fn touching(refusal: String, a_less_b: &str, b_less_a: &str) -> [String; 6] {
    [
        refusal.clone(),
        refusal,
        "empty".into(),
        "empty".into(),
        a_less_b.into(),
        b_less_a.into(),
    ]
}

/// Golden 48: notch tips of golden 32's prism's face across the line L
/// where the face crosses the box's top, within the tolerance of the top,
/// where the top holds part of rule 3's test set more than ρ from its
/// edges (design note 0012 §21, rule 3; found by the synthesis families).
/// Tilted 3° (1°), the tip 15 (30) tolerances across L in the face's
/// plane, under the top, over a box pierced under it by a triangle whose
/// apex lies half (0.49) of the way from the tip's foot to L: the foot
/// lies in the hole, and the tip's edges leave it over the top. Tilted 1°:
/// the tip 17.2 tolerances across L under a top cut along a line from 6
/// tolerances along L from the tip to a quarter of the way to its foot and
/// pierced likewise (the hole's apex 0.3 of the way); and the tip 34.4
/// (17.2) tolerances across L, from the near side above the top (from the
/// far side under it), under a top cut by a line from 6 tolerances before
/// the tip along L to half way to its foot and by a line at 70° (140°) to
/// L beyond the foot. Main failed every operation of each but the 140°
/// case, where it refused A ∪ B and B ∪ A, found A ∩ B and B ∩ A empty and
/// answered A − B and B − A, though the faces lie within the tolerance of
/// each other over the tip's edges' feet, which neither the meeting nor
/// classification separates; every operation is now refused at the tip. Without rule 3, without its test
/// of V's edges (V's foot alone), or passing a notch tip with no stretch
/// in front of it as an apex is passed, an operation fails as on main.
#[test]
fn a_notch_tip_is_refused_where_the_top_holds_its_edges() {
    let eps = k().tolerance().linear;
    let p = DVec2::new;
    for (what, t, v, k_apex, at) in [
        (
            "3°",
            3.0f64,
            -15.0,
            0.5,
            "(500.000000, -0.001498, -0.000079) mm",
        ),
        (
            "1°",
            1.0,
            -30.0,
            0.49,
            "(500.000000, -0.003000, -0.000052) mm",
        ),
    ] {
        let b = on_plane(
            t,
            &notch_profile((0.5, 2.0), 0.5, v * eps, (0.1, 0.1), false),
            false,
        );
        let yv = v * eps * t.to_radians().cos();
        let (s, apex) = (yv.signum(), yv * (1.0 - k_apex));
        let hole = vec![p(0.5 - 0.2, s * 0.5), p(0.5 + 0.2, s * 0.5), p(0.5, apex)];
        outcomes(
            what,
            &pierced_by(&the_box(), hole),
            &b,
            all(band_refusal(at)),
        );
    }
    let t = 1.0f64;
    {
        let h = 0.3 / t.to_radians().sin();
        let b = on_plane(
            t,
            &notch_profile((0.5, 2.0), 0.5, -h * eps, (0.1, 0.1), false),
            false,
        );
        let yv = -h * eps * t.to_radians().cos();
        let (sb, foot) = (yv.signum(), p(0.5, yv));
        let o1 = p(0.5 + 6.0 * eps, 0.0);
        let mut n1 = (p(0.5, yv * 0.25) - o1).normalize().perp();
        if n1.dot(foot - o1) > 0.0 {
            n1 = -n1;
        }
        let top = slab(&ccw(clip(&top_rect(), o1, n1)));
        let apex = yv * (1.0 - 0.3);
        let hole = vec![p(0.5 - 0.2, sb * 0.5), p(0.5 + 0.2, sb * 0.5), p(0.5, apex)];
        let at = "(500.000000, -0.001719, -0.000030) mm";
        outcomes(
            "cut short",
            &pierced_by(&top, hole),
            &b,
            all(band_refusal(at)),
        );
    }
    for (what, z, near, phi, at) in [
        (
            "70°",
            0.6f64,
            true,
            70.0f64,
            "(500.000000, 0.003437, 0.000060) mm",
        ),
        (
            "140°",
            0.3,
            false,
            140.0,
            "(500.000000, -0.001719, -0.000030) mm",
        ),
    ] {
        let h = z / t.to_radians().sin();
        let len = (1.0 - z) / t.to_radians().tan();
        let w = if near { h } else { -h };
        let b = on_plane(
            t,
            &notch_profile((0.5, 2.0), 0.5, w * eps, (0.1, 0.1), near),
            near,
        );
        let yv = w * eps * t.to_radians().cos();
        let (sb, foot) = (yv.signum(), p(0.5, yv));
        let o1 = p(0.5 - 6.0 * eps, 0.0);
        let mut n1 = (p(0.5, yv * 0.5) - o1).normalize().perp();
        if n1.dot(foot - o1) > 0.0 {
            n1 = -n1;
        }
        let o2 = p(0.5, yv + sb * 0.15 * len * eps);
        let mut n2 = p(phi.to_radians().cos(), phi.to_radians().sin()).perp();
        if n2.dot(foot - o2) > 0.0 {
            n2 = -n2;
        }
        let a = slab(&ccw(clip(&clip(&top_rect(), o1, n1), o2, n2)));
        outcomes(what, &a, &b, all(band_refusal(at)));
    }
}

/// Golden 49: apexes of golden 32's prism's face across L within the
/// tolerance of the top: reflex corners both of whose edges run away from
/// the top (design note 0012 §21, rule 3 for apexes; the synthesis
/// families and the soundness review's apex family). A notch from the far
/// side, its apex 0.3 tolerances above the top at 1°, and from the near
/// side, its apex 0.3 tolerances under the top with its prism below its
/// face at 3° and 10°, each over a top ending 3 tolerances beyond the
/// apex's foot; and a notch at x = 0.631703 with a lopsided mouth from the
/// far side at 1°, its apex 0.05 (0.118) tolerances above a top ending
/// 28.65 tolerances across L. Main failed every operation of each, a
/// contact's pieces disagreeing on their side; every operation is now
/// refused at the apex. Without rule 3's test of apexes, they fail as on
/// main.
#[test]
fn an_apex_across_the_line_is_refused_where_the_top_holds_its_foot() {
    let eps = k().tolerance().linear;
    let p = DVec2::new;
    for (t, near, at) in [
        (1.0f64, false, "(500.000000, 0.001719, 0.000030) mm"),
        (3.0, true, "(500.000000, -0.000572, -0.000030) mm"),
        (10.0, true, "(500.000000, -0.000170, -0.000030) mm"),
    ] {
        let h = 0.3 * eps / t.to_radians().sin();
        let w = if near { -h } else { h };
        let b = on_plane(
            t,
            &notch_profile((1.2, 2.2), 0.5, w, (0.1, 0.1), near),
            near,
        );
        let y = w * t.to_radians().cos();
        let s = y.signum();
        let a = slab(&ccw(clip(
            &top_rect(),
            p(0.5, y + s * 3.0 * eps),
            p(0.0, s),
        )));
        outcomes(&format!("{t}°"), &a, &b, all(band_refusal(at)));
    }
    for (z, at) in [
        (0.05f64, "(631.703000, 0.000286, 0.000005) mm"),
        (0.118, "(631.703000, 0.000676, 0.000012) mm"),
    ] {
        let t = 1.0f64;
        let st = t.to_radians().sin();
        let (x, h) = (0.631703, z * eps / st);
        let b = on_plane(
            t,
            &notch_profile((1.2, 2.2), x, h, (0.0129, 0.018), false),
            false,
        );
        let e = 0.5 / st;
        let a = if z < 0.1 {
            cuboid(
                DVec3::new(-0.25, -1.0, -1.0),
                DVec3::new(1.25, e * eps, 0.0),
            )
        } else {
            slab(&ccw(clip(&top_rect(), p(x, e * eps), p(0.0, 1.0))))
        };
        outcomes(&format!("lopsided {z}"), &a, &b, all(band_refusal(at)));
    }
}

/// Golden 50: notch tips across L within the tolerance of the top that
/// rule 3 does not refuse, every outcome as on main (design note 0012
/// §21; found by the synthesis families and the notch fuzz). The tip from
/// the near side 0.3 tolerances above the top at 1°, the top ending 1.05
/// tolerances beyond the tip's foot: every operation answered. At 0.5°,
/// the top's edge 2 tolerances across L on the side away from the tip: the
/// tip is no isolated touch of the top (rule 3's isolation), and a touch
/// along the top's edge refuses the unions. Two notches the fuzz found at 1° and
/// 0.5°, their tips within the tolerance of a stretch of the meeting: the
/// split refuses them, at the tip or (its stretch meeting the top's cut
/// edge) at the box's end, about 740 mm away. Without rule 3's isolation
/// the tip over the top's edge is refused in every operation with the
/// strip's text, and with V's foot tested a rounding margin from the top's
/// edges rather than at ρ so is the fuzz's notch at 1°.
#[test]
fn a_notch_tip_left_to_the_meeting_or_the_split() {
    let eps = k().tolerance().linear;
    let p = DVec2::new;
    {
        let t = 1.0f64;
        let h = 0.3 * eps / t.to_radians().sin();
        let b = on_plane(
            t,
            &notch_profile((1.2, 2.2), 0.5, h, (0.1, 0.1), true),
            false,
        );
        let foot = h * t.to_radians().cos();
        let a = cuboid(
            DVec3::new(-0.25, -1.0, -1.0),
            DVec3::new(1.25, foot + 1.05 * eps, 0.0),
        );
        let want = [
            "(15, 39, 26, 1) V=3.172020743 A=18.089640880",
            "(15, 39, 26, 1) V=3.172020743 A=18.089640880",
            "(8, 17, 11, 1) V=0.002981906 A=1.420271531",
            "(8, 17, 11, 1) V=0.002981906 A=1.420271531",
            "(11, 26, 17, 1) V=1.497020829 A=8.011334699",
            "(12, 30, 20, 1) V=1.672018008 A=11.498577713",
        ];
        outcomes("beside", &a, &b, want.map(String::from));
    }
    {
        let t = 0.5f64;
        let h = 0.3 * eps / t.to_radians().sin();
        let b = on_plane(
            t,
            &notch_profile((1.2, 2.2), 0.5, h, (0.1, 0.1), true),
            false,
        );
        let a = cuboid(
            DVec3::new(-0.25, -2.0 * eps, -1.0),
            DVec3::new(1.25, 1.0, 0.0),
        );
        let want = touching(
            edge_refusal("(124.999636, -0.000200, -0.000000) mm"),
            "(6, 12, 8, 1) V=1.500000300 A=8.000001000",
            "(9, 21, 14, 1) V=1.674999828 A=11.509904635",
        );
        outcomes("over the edge", &a, &b, want);
    }
    let b = on_plane(
        1.0,
        &[
            p(-1.2, -0.5),
            p(2.2, -0.5),
            p(2.2, 0.5),
            p(0.5585068932436936, 0.5),
            p(0.47350998224935675, -5.414409039095355e-6),
            p(0.296696907214841, 0.5),
            p(-1.2, 0.5),
        ],
        true,
    );
    let a = cuboid(
        DVec3::new(-0.25, -1.0, -1.0),
        DVec3::new(1.25, -5.353927718941095e-6, 0.0),
    );
    let at = "(473.509982, -0.005414, -0.000094) mm";
    let (near, cross) = (near_edge_refusal(at), crossing_refusal(at));
    let want = [0, 1, 2, 3, 4, 5].map(|i| {
        if i % 2 == 0 {
            near.clone()
        } else {
            cross.clone()
        }
    });
    outcomes("fuzz at 1°", &a, &b, want);
    let b = on_plane(
        0.5,
        &[
            p(-1.2, -0.5),
            p(0.29894098938327757, -0.5),
            p(0.4504511852793128, 5.729650674006515e-6),
            p(0.4559764177328263, -0.5),
            p(2.2, -0.5),
            p(2.2, 0.5),
            p(-1.2, 0.5),
        ],
        true,
    );
    let o = p(0.4504514275595844, 5.550401201270718e-6);
    let n = p(-0.31746232165223376, -0.9482708865778668);
    let a = slab(&ccw(clip(&top_rect(), o, n)));
    let at = "(-250.000000, 234.502823, 0.000000) mm";
    outcomes("fuzz at 0.5°", &a, &b, all(crossing_refusal(at)));
}

/// Golden 51: corners of oblique prisms' faces across L within the
/// tolerance of the box's top, whose edges out of the face at the corner
/// lean far from the face's normal (the soundness review's oblique family;
/// design note 0012 §21, rule 3's isolation). Rule 3 lets a corner off as
/// no isolated touch only where both of its tests say so: the face's
/// normals (its body behind the corner lying across the corner's way to L)
/// and the pair's touch record (the body's edges at the corner leaving P to
/// one side only). A notch tip from the near side 0.7 tolerances above the
/// top at 30°, its prism leaning 75° so its edge out of the face at the tip
/// leaves the tip to its neighbours' side, a touch the record holds, though
/// the normals do not call it isolated: main failed A ∪ B, A ∩ B and A − B
/// as a contact's pieces that disagree on their side. A notch tip from the
/// far side 0.3 tolerances under the top at 10°, its prism below its face
/// leaning 85°: main failed the same three. An apex of a notch 0.3
/// tolerances above a top ending 3 tolerances beyond its foot at 5°, its
/// prism leaning 87.5° so its edge out of the face leaves the apex to P's
/// other side, which the record does not count a touch while the normals
/// call it isolated: main failed every operation as kept pieces that do
/// not close. Every operation is now refused at the corner. With the
/// normals' test alone the near tip fails as on main, and with the
/// record's alone the apex does.
#[test]
fn an_oblique_corner_is_refused_unless_both_tests_call_it_no_isolated_touch() {
    let eps = k().tolerance().linear;
    let p = DVec2::new;
    {
        let t = 30.0f64;
        let h = 0.7 * eps / t.to_radians().sin();
        let prof = notch_profile((0.5, 2.0), 0.5, h, (0.35, 0.35), true);
        let b = oblique(t, &prof, false, 75.0, -90.0);
        let at = "(500.000000, 0.000121, 0.000070) mm";
        outcomes("near tip", &the_box(), &b, all(band_refusal(at)));
    }
    {
        let t = 10.0f64;
        let h = 0.3 * eps / t.to_radians().sin();
        let prof = notch_profile((0.5, 2.0), 0.5, -h / eps * eps, (0.35, 0.35), false);
        let b = oblique(t, &prof, true, 85.0, 90.0);
        let at = "(500.000000, -0.000170, -0.000030) mm";
        outcomes("far tip", &the_box(), &b, all(band_refusal(at)));
    }
    {
        let t = 5.0f64;
        let h = 0.3 * eps / t.to_radians().sin();
        let prof = notch_profile((1.2, 2.2), 0.5, h, (0.1, 0.1), false);
        let b = oblique(t, &prof, false, 87.5, -90.0);
        let foot = p(0.5, h * t.to_radians().cos());
        let a = slab(&ccw(clip(
            &top_rect(),
            p(foot.x, foot.y + 3.0 * eps),
            p(0.0, 1.0),
        )));
        let at = "(500.000000, 0.000343, 0.000030) mm";
        outcomes("apex", &a, &b, all(band_refusal(at)));
    }
}

/// Golden 52: spike tips V across L within the tolerance of the top, whose
/// flap (the face across L between V's edges) the flap rule's table
/// decides (design note 0012 §21; found by the synthesis and oblique
/// families). The flap's own side is tested where the layer between it
/// and the top lies outside the face's body, where the touch record does
/// not hold V, or where the layer lies outside the top's body and the
/// touch at V is unpinched, and each of these is refused: a rising spike
/// 50 tolerances across at 1°, its prism below its face, under a top cut
/// along a line at 20° through V's foot (the layer in the face's body and
/// outside the top's, V on the top's edge; main crossed the flap and the
/// top 22.4 tolerances inside in A − B); spikes 0.9 tolerances from the
/// top at 10° and 3°, their prisms leaning 85° and 88.5° so that the record
/// does not hold V (main crossed them 4.5 tolerances inside in A − B, and
/// answered every operation with the flap crossing the top 14.9 tolerances
/// inside). Where the record holds V and the layer lies in both bodies, or
/// outside the top's with the touch pinched and the top's edges out of the
/// tolerance's reach of the flap, the touch at V decides it, every outcome
/// as on main: a spike 0.3 tolerances above the top at 0.5° with a
/// lopsided root, its prism below its face, the same 0.9 tolerances from the top at 10° unleaning,
/// and 0.5 tolerances at 10° under a top ending 3 tolerances beyond V's
/// foot (A − B touching at V, refused); a falling spike 3 tolerances
/// across L at 10° under a top cut along L 0.75 of the way to V's foot, and the walls of a notch
/// whose tip lies 0.3 tolerances above the top at 30° under a top ending
/// at a corner 0.4 tolerances beyond its foot (the union touching at V,
/// refused, or every operation answered). Letting the touch decide a flap
/// the record does not hold answers the 85° spike as main did; testing
/// only flaps outside the face's body answers the 50-tolerance spike as
/// main did; testing every flap, or a pinched one's side whole, refuses
/// the lopsided spike; and testing an unpinched flap in both bodies
/// refuses the falling spike with the strip's text.
#[test]
fn a_spike_flap_is_tested_unless_the_touch_at_its_tip_decides_it() {
    let eps = k().tolerance().linear;
    let p = DVec2::new;
    let zero = all(band_refusal("(500.000000, 0.000000, 0.000000) mm"));
    {
        let (t, d) = (1.0f64, 50.0);
        let b = on_plane(
            t,
            &spiked_profile((0.5, 2.0), 0.5, d * eps, (0.1, 0.1), 0.2, true),
            true,
        );
        let yv = d * eps * t.to_radians().cos();
        let n = p(20f64.to_radians().cos(), 20f64.to_radians().sin()).perp();
        let a = slab(&ccw(clip(&top_rect(), p(0.5, yv), n)));
        outcomes("outside the top's body, unpinched", &a, &b, zero.clone());
    }
    let up_10 = |alpha: f64, phi: f64| {
        let h = 0.9 * eps / 10f64.to_radians().sin();
        let prof = spiked_profile((0.5, 2.0), 0.5, h / eps * eps, (0.35, 0.35), 0.2, true);
        (h, oblique(10.0, &prof, true, alpha, phi))
    };
    {
        let (h, b) = up_10(85.0, 90.0);
        let foot = p(0.5, h * 10f64.to_radians().cos());
        let a = slab(&ccw(clip(
            &top_rect(),
            p(foot.x, foot.y + 3.0 * eps),
            p(0.0, 1.0),
        )));
        outcomes("leaning 85°", &a, &b, zero.clone());
    }
    {
        let t = 3.0f64;
        let h = 0.9 * eps / t.to_radians().sin();
        let prof = spiked_profile((0.5, 2.0), 0.5, -h / eps * eps, (0.35, 0.35), 0.2, false);
        let b = oblique(t, &prof, false, 88.5, -90.0);
        let at = "(500.000000, -0.000000, 0.000000) mm";
        outcomes("leaning 88.5°", &the_box(), &b, all(band_refusal(at)));
    }
    {
        let t = 0.5f64;
        let h = 0.3 * eps / t.to_radians().sin();
        let b = on_plane(
            t,
            &spiked_profile((1.2, 2.2), 0.5, h, (0.12, 0.02), 0.4, true),
            true,
        );
        let want = answered(
            "(16, 36, 24, 1) V=3.095000000 A=15.280000000",
            "(9, 21, 14, 1) V=0.089000120 A=2.295059819",
            point_refusal("(500.000000, 0.003438, 0.000030) mm"),
            "(12, 24, 16, 2) V=0.095000000 A=2.480000000",
        );
        outcomes("lopsided", &the_box(), &b, want);
    }
    {
        let (_, b) = up_10(0.0, 0.0);
        let want = answered(
            "(16, 36, 24, 1) V=3.129000000 A=14.460000000",
            "(9, 21, 14, 1) V=0.223600078 A=2.633677667",
            point_refusal("(500.000000, 0.000510, 0.000090) mm"),
            "(12, 24, 16, 2) V=0.129000000 A=1.976000000",
        );
        outcomes("unleaning", &the_box(), &b, want);
    }
    {
        let t = 10.0f64;
        let d = 0.5 / t.to_radians().sin();
        let b = on_plane(
            t,
            &spiked_profile((0.5, 2.0), 0.5, d * eps, (0.5, 0.5), 0.2, true),
            true,
        );
        let yv = d * eps * t.to_radians().cos();
        let a = slab(&ccw(clip(&top_rect(), p(0.5, yv + 3.0 * eps), p(0.0, 1.0))));
        let want = answered(
            "(19, 42, 28, 1) V=1.653239504 A=9.626224014",
            "(10, 24, 16, 1) V=0.271761443 A=2.912295780",
            point_refusal("(500.000000, 0.000284, 0.000050) mm"),
            "(16, 30, 20, 3) V=0.153238629 A=2.450028113",
        );
        outcomes("short top", &a, &b, want);
    }
    {
        let (t, d) = (10.0f64, 3.0);
        let b = on_plane(
            t,
            &spiked_profile((0.5, 2.0), 0.5, -d * eps, (0.1, 0.1), 0.2, false),
            false,
        );
        let yv = -d * eps * t.to_radians().cos();
        let a = slab(&ccw(clip(&top_rect(), p(0.5, 0.75 * yv), p(0.0, -1.0))));
        let want = touching(
            point_refusal("(500.000000, -0.000295, -0.000052) mm"),
            "(6, 12, 8, 1) V=1.500000332 A=8.000001108",
            "(9, 21, 14, 1) V=0.385000015 A=4.463607126",
        );
        outcomes("falling", &a, &b, want);
    }
    {
        let t = 30.0f64;
        let h = 0.3 * eps / t.to_radians().sin();
        let b = on_plane(
            t,
            &notch_profile((0.5, 2.0), 0.5, h, (0.1, 0.1), true),
            true,
        );
        let foot = p(0.5, h * t.to_radians().cos());
        let o = foot + p(1.0, 1.0) * (0.4 * eps / 2f64.sqrt());
        let a = slab(&ccw(clip(
            &clip(&top_rect(), o, p(0.0, 1.0)),
            o,
            p(1.0, 0.0),
        )));
        let want = [
            "(15, 36, 23, 1) V=1.852737752 A=12.217545001",
            "(15, 36, 23, 1) V=1.852737752 A=12.217544999",
            "(6, 12, 8, 1) V=0.122262334 A=1.592357374",
            "(6, 12, 8, 1) V=0.122262334 A=1.592357367",
            "(10, 22, 14, 1) V=0.627737755 A=5.404003235",
            "(11, 26, 17, 1) V=1.102737663 A=8.405899169",
        ];
        outcomes("corner", &a, &b, want.map(String::from));
    }
}

/// Golden 53: spike tips V across L within the tolerance of the top, the
/// layer between the flap and the top in the face's body and outside the
/// top's, the record holding V, and the top holding V's foot (a pinched
/// flap; design note 0012 §21, the flap rule's crossing test; found by the
/// synthesis and oblique families). The touch at V decides the flap only
/// where the top's boundary stays out of the tolerance's reach of the
/// flap's strip; here an edge of the top passes within it, each prism
/// below its face: at 1°, V 0.3 tolerances above the top, the top cut
/// along a line at 150° 1.5 (3) tolerances beyond V's foot; at 3°, V 12
/// tolerances across, the top cut at 110° half way from V's foot to L; at
/// 30°, V 0.9 tolerances above the top, the top's edge along L a tolerance
/// short of V's foot, the top on V's side of it. Main failed every
/// operation, a contact's pieces disagreeing on their side; every
/// operation is now refused at V. Without the crossing test (the flap left
/// to the touch), or with the flap's own side tested instead, the first
/// fails as on main.
#[test]
fn a_pinched_flap_is_refused_where_the_tops_edge_comes_within_its_strip() {
    let eps = k().tolerance().linear;
    let p = DVec2::new;
    let spike = |t: f64, w: f64, (a, b): (f64, f64)| {
        on_plane(
            t,
            &spiked_profile((0.5, 2.0), 0.5, w, (a, b), 0.2, true),
            true,
        )
    };
    let at = all(band_refusal("(500.000000, 0.001719, 0.000030) mm"));
    for dl in [-1.5f64, -3.0] {
        let t = 1.0f64;
        let h = 0.3 * eps / t.to_radians().sin();
        let b = spike(t, h / eps * eps, (0.3, 0.3));
        let foot = p(0.5, h * t.to_radians().cos());
        let mut n = p(150f64.to_radians().cos(), 150f64.to_radians().sin()).perp();
        if n.dot(p(foot.x, 0.0) - foot) < 0.0 {
            n = -n;
        }
        let a = slab(&ccw(clip(&top_rect(), foot - n * (dl * eps), n)));
        outcomes(&format!("150° {dl}"), &a, &b, at.clone());
    }
    {
        let (t, d) = (3.0f64, 12.0);
        let b = spike(t, d * eps, (0.02, 0.3));
        let yv = d * eps * t.to_radians().cos();
        let n = p(110f64.to_radians().cos(), 110f64.to_radians().sin()).perp();
        let a = slab(&ccw(clip(&top_rect(), p(0.5, 0.5 * yv), n)));
        let at = "(500.000000, 0.001198, 0.000063) mm";
        outcomes("110°", &a, &b, all(band_refusal(at)));
    }
    {
        let t = 30.0f64;
        let d = 0.9 / t.to_radians().sin();
        let b = spike(t, d * eps, (0.2, 0.2));
        let yv = d * eps * t.to_radians().cos();
        let a = slab(&ccw(clip(&top_rect(), p(0.5, yv - eps), p(0.0, -1.0))));
        let at = "(500.000000, 0.000156, 0.000090) mm";
        outcomes("short", &a, &b, all(band_refusal(at)));
    }
}

/// Golden 54: spike tips whose flap lies over the strip of a run of the
/// other face (design note 0012 §21, the exclusion; found by the
/// synthesis families and the fuzz). The run's strip leaves out the flap's
/// triangle, which the flap's own rule decides by the touch at its tip: a
/// falling spike 0.3 tolerances under the top at 0.5°, under a slab whose
/// top has a notch from x = 0.3 to 0.7 ending 1.05 tolerances beyond the
/// tip's foot; and a notch the fuzz found at 2°, under a box turned 12.7°
/// whose corner lies at its tip. Main refused every operation, the run's
/// strip held over the flap; now the touch at the tip refuses the unions,
/// A ∩ B and B ∩ A are empty, and A − B and B − A are answered. Without the
/// exclusion, or with only runs that are their face's whole interval
/// leaving flaps out, the slab's case is refused in every operation as on
/// main.
#[test]
fn a_strip_leaves_out_the_other_faces_flap() {
    let eps = k().tolerance().linear;
    let p = DVec2::new;
    {
        let t = 0.5f64;
        let h = 0.3 * eps / t.to_radians().sin();
        let b = on_plane(
            t,
            &spiked_profile((1.2, 2.2), 0.5, -h, (0.12, 0.02), 0.4, false),
            false,
        );
        let y = -h * t.to_radians().cos() - 1.05 * eps;
        let a = slab(&ccw(vec![
            p(-0.25, 1.0),
            p(1.25, 1.0),
            p(1.25, -0.9),
            p(0.7, -0.9),
            p(0.7, y),
            p(0.3, y),
            p(0.3, -0.9),
            p(-0.25, -0.9),
        ]));
        let want = touching(
            point_refusal("(500.000000, -0.003438, -0.000030) mm"),
            "(10, 24, 16, 1) V=2.490001417 A=13.579995749",
            "(9, 21, 14, 1) V=0.184000120 A=4.575059819",
        );
        outcomes("slab", &a, &b, want);
    }
    let b = on_plane(
        2.0,
        &[
            p(-1.2, -0.5),
            p(0.5359088166745397, -0.5),
            p(0.5804842176658804, -2.578833751305944e-6),
            p(0.6103051633978451, -0.5),
            p(2.2, -0.5),
            p(2.2, 0.5),
            p(-1.2, 0.5),
        ],
        false,
    );
    let frame = Frame {
        origin: DVec3::new(0.5804844468118274, -2.612368254420224e-6, -1.0),
        x: DVec3::new(0.9757129766956367, 0.21905293220530012, 0.0),
        y: DVec3::new(-0.21905293220530012, 0.9757129766956367, 0.0),
        z: DVec3::Z,
    };
    let a = k().make_box(&frame, DVec3::ONE).expect("a box");
    let want = touching(
        point_refusal("(580.484447, -0.002612, 0.000000) mm"),
        "(6, 12, 8, 1) V=1.000000000 A=6.000000000",
        "(9, 21, 14, 1) V=1.690700505 A=11.627037040",
    );
    outcomes("fuzz", &a, &b, want);
}

/// Golden 55: spike tips V across L within the tolerance of a slab's top
/// whose cluster lies on the top's boundary, so the touch at V excuses
/// the flap's points about it (design note 0012 §21, the member discs).
/// The excused region is the union of discs of radius ε about the places
/// of V's cluster's members, which the settling of touched points decides,
/// the same in either order of the operands. At 3°, V 0.1 tolerances from
/// the top, the top ending at a corner 0.6 tolerances beyond V's foot, the
/// spike rising (falling); at 20° (5°), V 2.6 (1.6) tolerances across,
/// the top ending along L at V's foot (half a tolerance short of it); at
/// 10°, V 3 tolerances across, the top cut along a line at 20° through V's
/// foot. Main answered A ∩ B and B ∩ A, or A − B, or B − A, with the flap
/// crossing the top up to 2 tolerances inside; every operation is now
/// refused at V. With a disc of three tolerances about the cluster's
/// place, the rising spike at 3° is refused in one order of the operands
/// and answered in the other, as main did.
#[test]
fn a_flap_is_excused_only_within_the_tolerance_of_its_tips_members() {
    let eps = k().tolerance().linear;
    let p = DVec2::new;
    let spike = |t: f64, d: f64, wa: f64, up: bool, below: bool| {
        on_plane(
            t,
            &spiked_profile((0.5, 2.0), 0.5, d * eps, (wa, wa), 0.2, up),
            below,
        )
    };
    for (up, at) in [
        (true, "(499.999846, -0.000000, 0.000000) mm"),
        (false, "(499.999846, 0.000000, 0.000000) mm"),
    ] {
        let t = 3.0f64;
        let s = if up { 1.0 } else { -1.0 };
        let d = s * 0.1 / t.to_radians().sin();
        let b = spike(t, d, 0.5, up, !up);
        let o = p(0.5, d * eps * t.to_radians().cos() + s * 0.6 * eps);
        let a = slab(&ccw(clip(&clip(&top_rect(), o, p(0.0, s)), o, p(1.0, 0.0))));
        outcomes(&format!("corner {up}"), &a, &b, all(band_refusal(at)));
    }
    let zero = all(band_refusal("(500.000000, 0.000000, 0.000000) mm"));
    for (t, d, wa, below, dl) in [
        (20.0f64, 2.6, 0.25, true, 0.0),
        (5.0, 1.6, 0.2, false, -0.5),
    ] {
        let b = spike(t, d, wa, true, below);
        let yv = d * eps * t.to_radians().cos();
        let a = slab(&ccw(clip(&top_rect(), p(0.5, yv + dl * eps), p(0.0, 1.0))));
        outcomes(&format!("along L at {t}°"), &a, &b, zero.clone());
    }
    {
        let (t, d) = (10.0f64, 3.0);
        let b = spike(t, d, 0.1, true, false);
        let yv = d * eps * t.to_radians().cos();
        let n = p(20f64.to_radians().cos(), 20f64.to_radians().sin()).perp();
        let a = slab(&ccw(clip(&top_rect(), p(0.5, yv), n)));
        outcomes("at 20°", &a, &b, zero);
    }
}

/// A kernel a hundred times finer than the default, which builds bodies
/// that the default tolerance takes for touching themselves: the oracle's
/// own tests.
fn fine() -> BrepKernel {
    BrepKernel::with_tolerance(Tolerance {
        linear: 1e-9,
        ..Tolerance::default()
    })
}

/// The oracle reports a body whose wall is thinner than half the
/// tolerance, a C whose slot ends `w` from its outer side (built by a
/// finer kernel; the default refuses the profile), and passes it at
/// three fifths of the tolerance.
#[test]
fn the_oracle_finds_a_wall_thinner_than_half_the_tolerance() {
    let eps = k().tolerance().linear;
    let c = |w: f64| {
        let pts = [
            (0.0, 0.0),
            (1.0, 0.0),
            (1.0, 0.45),
            (w, 0.45),
            (w, 0.55),
            (1.0, 0.55),
            (1.0, 1.0),
            (0.0, 1.0),
        ]
        .map(|(x, y)| DVec2::new(x, y));
        let body = fine()
            .extrude(
                &Profile::polygon(Frame::WORLD, &pts),
                &ExtrudeSpec {
                    distance: 1.0,
                    symmetric: false,
                },
            )
            .expect("a C");
        self_contact::touches_itself(&k(), &body, &TessellationOptions::default()).expect("a mesh")
    };
    let thin = c(0.3 * eps).expect("a report");
    assert!(thin.contains("comes within 3.000e-8 m of face"), "{thin}");
    assert_eq!(c(0.6 * eps), None);
}

/// A body's diamond ridge along x (A) and the ridge of a second diamond
/// along y three tenths of the tolerance above it (B), crossing at a
/// quarter of each, made one body by a finer kernel. Neither ridge's ends
/// nor middle comes near the other body, so only measuring the edges whole
/// finds the crossing; the oracle sampled them before.
#[test]
fn the_oracle_measures_straight_edges_whole() {
    let eps = k().tolerance().linear;
    let diamond = |c: DVec2| {
        [(0.0, -0.5), (0.5, 0.0), (0.0, 0.5), (-0.5, 0.0)].map(|(x, y)| c + DVec2::new(x, y))
    };
    let solid = |frame: Frame, pts: &[DVec2]| {
        fine()
            .extrude(
                &Profile::polygon(frame, pts),
                &ExtrudeSpec {
                    distance: 1.0,
                    symmetric: false,
                },
            )
            .expect("a prism")
    };
    // A's ridge: (0..1, 0.5, 0.5). B's: (0.25, 0.25..1.25, 0.5 + 0.3 ε).
    let a = solid(yz(DVec3::ZERO), &diamond(DVec2::new(0.5, 0.0)));
    let b = solid(
        xz(DVec3::new(0.0, 1.25, 0.0)),
        &diamond(DVec2::new(0.25, 1.0 + 0.3 * eps)),
    );
    let both = fine()
        .boolean(&a, &b, BooleanOp::Union)
        .expect("two solids");
    assert_eq!(counts(&both).3, 2);
    let report = self_contact::touches_itself(&k(), &both, &TessellationOptions::default())
        .expect("a mesh")
        .expect("a report");
    assert!(report.contains("at (250.0"), "{report}");
}

/// The oracle lets off the acute tip of a notch split by a meeting (design
/// note 0012 §21, the notch-tip probes): B's bottom face, on the plane
/// tilted 40° about the x axis, has two notches 0.1 wide at their mouths
/// from its far side, their tips 2 tolerances before and 10 after x = 0.5,
/// half a tolerance above and below the top z = 0 of the box A, and
/// between them a spike up to the far side at x = 0.5. A ∪ B keeps B's
/// bottom face above z = 0 as three faces, the spike and one either side,
/// which meet only at the tips: there each notch's edge along the spike
/// leaves the face across the notch 11.3° from its edge there, within half
/// the tolerance of it for 2.6 tolerances from the tip. The mesh of the
/// right one has a point 1.6 tolerances from the tip and 0.31 from that
/// face, which the oracle reported before it let such corners off.
#[test]
fn the_oracle_lets_off_a_notch_tip_split_by_a_meeting() {
    let eps = k().tolerance().linear;
    let t = 40.0_f64;
    // The distance up the bottom face of a point `z` tolerances above z = 0.
    let up = |z: f64| z * eps / t.to_radians().sin();
    let (left, right) = (0.5 - 2.0 * eps, 0.5 + 10.0 * eps);
    let pts = [
        (-0.5, -0.5),
        (2.0, -0.5),
        (2.0, 0.5),
        (right + 0.1, 0.5),
        (right, up(-0.5)),
        (0.5, 0.5),
        (left, up(0.5)),
        (left - 0.1, 0.5),
        (-0.5, 0.5),
    ]
    .map(|(x, y)| DVec2::new(x, y));
    let a = cuboid(DVec3::new(-0.25, -1.0, -1.0), DVec3::new(1.25, 1.0, 0.0));
    let b = prism(tilted_frame(t), &pts, 0.5);
    let r = run(&a, &b, BooleanOp::Union);
    let opts = TessellationOptions::default();
    assert_eq!(self_contact::touches_itself(&k(), &r, &opts), Ok(None));
    // The closeness is in the mesh the oracle measures: a point of an edge
    // more than the tolerance from its ends, within half the tolerance of a
    // face it does not bound, about the right tip.
    let tip = tilted_frame(t).y * up(-0.5) + DVec3::new(right, 0.0, 0.0);
    let mesh = k().tessellate(&r, &opts).expect("a mesh");
    let adj = k().adjacency(&r);
    let mut close = Vec::new();
    for (e, line) in mesh.edges.iter().enumerate() {
        let (p0, p1) = (line[0], line[line.len() - 1]);
        let points = (line.iter().copied()).chain(line.windows(2).map(|s| 0.5 * (s[0] + s[1])));
        for p in points.filter(|&p| (p - p0).length() > eps && (p - p1).length() > eps) {
            for (f, face) in mesh.faces.iter().enumerate() {
                if adj.edge_faces[e].iter().any(|x| x.0 as usize == f) {
                    continue;
                }
                let plane = k().face_plane(&r, FaceIx(f as u32));
                let onto = |q: DVec3| match &plane {
                    Some(pl) => q - pl.normal * (q - pl.point).dot(pl.normal),
                    None => q,
                };
                let d = (face.triangles.iter())
                    .map(|t| {
                        let [a, b, c] = t.map(|i| onto(face.positions[i as usize]));
                        self_contact::triangle_distance(p, a, b, c)
                    })
                    .fold(f64::INFINITY, f64::min);
                if d <= 0.5 * eps {
                    close.push(((p - tip).length() / eps, d / eps));
                }
            }
        }
    }
    assert!(
        !close.is_empty(),
        "the mesh has no point the corner lets off"
    );
    for (from_tip, d) in close {
        assert!(
            from_tip > 1.0 && from_tip < 2.6 && d > 0.3,
            "{from_tip} {d}"
        );
    }
}

/// The oracle reports an edge running along an edge of a face from a
/// corner the two share: B's bottom face as above, with one notch from its
/// far side whose tip is on the x axis at x = 0.5, its right edge to
/// (0.6, 0.5) and its left edge running beside the right one to a point
/// 0.3 of the way along it and 0.6 tolerances to the left of it, where the
/// notch opens to (0.3, 0.5). A finer kernel's A ∪ B splits B's bottom
/// face at the tip, and the left edge is within half the tolerance of the
/// face across the notch for half its length, not only about the corner:
/// its middle, 0.3 tolerances from that face's edge, is reported.
#[test]
fn the_oracle_finds_an_edge_running_along_a_faces_edge_from_their_corner() {
    let eps = k().tolerance().linear;
    let (tip, mouth) = (DVec2::new(0.5, 0.0), DVec2::new(0.6, 0.5));
    let along = mouth - tip;
    let side = along.perp().normalize();
    let pts = [
        DVec2::new(-0.5, -0.5),
        DVec2::new(2.0, -0.5),
        DVec2::new(2.0, 0.5),
        mouth,
        tip,
        tip + 0.3 * along + 0.6 * eps * side,
        DVec2::new(0.3, 0.5),
        DVec2::new(-0.5, 0.5),
    ];
    let spec = ExtrudeSpec {
        distance: 0.5,
        symmetric: false,
    };
    let b = (fine().extrude(&Profile::polygon(tilted_frame(40.0), &pts), &spec)).expect("B");
    let a = (fine().make_box(
        &Frame {
            origin: DVec3::new(-0.25, -1.0, -1.0),
            ..Frame::WORLD
        },
        DVec3::new(1.5, 2.0, 1.0),
    ))
    .expect("A");
    let r = fine().boolean(&a, &b, BooleanOp::Union).expect("A ∪ B");
    assert_eq!(fine().validate(&r), Vec::<String>::new());
    let report = self_contact::touches_itself(&k(), &r, &TessellationOptions::default())
        .expect("a mesh")
        .expect("a report");
    assert!(
        report.contains("comes within 3.000e-8 m of face"),
        "{report}"
    );
}

/// The check for faces crossing ([`self_contact::faces_cross`]) passes
/// solids whose flat faces meet only along their edges: a box, an L-shaped
/// prism (whose bottom's plane meets the inner wall's along a line through
/// the bottom only), boxes stacked into a step, the plate with a hole, and
/// a cube standing on its corner.
#[test]
fn faces_cross_passes_faces_meeting_at_their_edges() {
    let l_prism = prism(
        Frame::WORLD,
        &[
            (0.0, 0.0),
            (2.0, 0.0),
            (2.0, 0.5),
            (1.0, 0.5),
            (1.0, 1.0),
            (0.0, 1.0),
        ]
        .map(|(x, y)| DVec2::new(x, y)),
        1.0,
    );
    let step = op(
        &cuboid(DVec3::ZERO, DVec3::new(2.0, 1.0, 0.5)),
        &cuboid(DVec3::ZERO, DVec3::new(1.0, 1.0, 1.0)),
        BooleanOp::Union,
    );
    for (name, body) in [
        ("box", cuboid(DVec3::ZERO, DVec3::ONE)),
        ("L", l_prism),
        ("step", step),
        ("holed plate", holed_plate(1.0)),
        ("cube on its corner", cube_on_corner(1.0, 0.0)),
    ] {
        let found = self_contact::faces_cross(&k(), &body, 1e-8).expect("a mesh");
        assert_eq!(found, None, "{name}");
    }
}

/// Two flat faces given as loops: the square |x|, |y| ≤ 1 of z = 0 and
/// the rectangle |x| ≤ w, |z| ≤ h of y = 0, crossing along the x axis.
fn crossing_faces(w: f64, h: f64) -> (self_contact::FlatFace, self_contact::FlatFace) {
    let loop_of = |pts: &[DVec3]| -> Vec<(DVec3, DVec3)> {
        (0..pts.len())
            .map(|i| (pts[i], pts[(i + 1) % pts.len()]))
            .collect()
    };
    let v = DVec3::new;
    let square = self_contact::FlatFace {
        point: DVec3::ZERO,
        normal: DVec3::Z,
        segments: loop_of(&[
            v(-1.0, -1.0, 0.0),
            v(1.0, -1.0, 0.0),
            v(1.0, 1.0, 0.0),
            v(-1.0, 1.0, 0.0),
        ]),
    };
    let wall = self_contact::FlatFace {
        point: DVec3::ZERO,
        normal: -DVec3::Y,
        segments: loop_of(&[v(-w, 0.0, -h), v(w, 0.0, -h), v(w, 0.0, h), v(-w, 0.0, h)]),
    };
    (square, wall)
}

/// [`self_contact::pair_crossing`] measures a crossing exactly: a wall
/// through a square, crossing it along the x axis, is as deep as its
/// half-height or half-width, whichever is less, at its middle; it is
/// reported a millionth of the tolerance deeper than the tolerance and not
/// a millionth shallower; faces sharing their edge along the line do not
/// cross; an open loop is an error.
#[test]
fn pair_crossing_measures_exactly() {
    let eps = k().tolerance().linear;
    let (square, wall) = crossing_faces(0.5, 0.25);
    let c = self_contact::pair_crossing(&square, &wall, eps)
        .expect("closed loops")
        .expect("a crossing");
    assert!((c.depth - 0.25).abs() < 1e-6 * eps, "{c:?}");
    assert!(c.point.length() < 1e-6 * eps, "{c:?}");
    assert!(
        (c.angle - std::f64::consts::FRAC_PI_2).abs() < 1e-15,
        "{c:?}"
    );
    for (h, found) in [(1.000001, true), (0.999999, false)] {
        let (square, wall) = crossing_faces(0.5, h * eps);
        let c = self_contact::pair_crossing(&square, &wall, eps).expect("closed loops");
        assert_eq!(c.is_some(), found, "{h}: {c:?}");
        if let Some(c) = c {
            assert!((c.depth - h * eps).abs() < 1e-6 * eps, "{c:?}");
        }
    }
    // The wall's edge on the square's: the half y ≥ 0 of the square and the
    // wall's half z ≤ 0 share the x axis.
    let (mut square, mut wall) = crossing_faces(0.5, 0.25);
    let v = DVec3::new;
    square.segments = vec![
        (v(-1.0, 0.0, 0.0), v(1.0, 0.0, 0.0)),
        (v(1.0, 0.0, 0.0), v(1.0, 1.0, 0.0)),
        (v(1.0, 1.0, 0.0), v(-1.0, 1.0, 0.0)),
        (v(-1.0, 1.0, 0.0), v(-1.0, 0.0, 0.0)),
    ];
    wall.segments = vec![
        (v(-1.0, 0.0, 0.0), v(-1.0, 0.0, -1.0)),
        (v(-1.0, 0.0, -1.0), v(1.0, 0.0, -1.0)),
        (v(1.0, 0.0, -1.0), v(1.0, 0.0, 0.0)),
        (v(1.0, 0.0, 0.0), v(-1.0, 0.0, 0.0)),
    ];
    assert_eq!(self_contact::pair_crossing(&square, &wall, eps), Ok(None));
    wall.segments.pop();
    wall.segments.pop();
    assert!(self_contact::pair_crossing(&square, &wall, eps).is_err());
}

/// The check for faces crossing ([`self_contact::faces_cross`]) on a
/// result whose faces fold under the tolerance (golden 39's edge 3
/// tolerances above the line, A ∩ B in both orders): B's corner 0.521
/// tolerances above A's top moves the planes' common line off the edge
/// the faces share, and the check finds them crossing 2.95 tolerances
/// deep, which is less than the corner's distance from the line
/// (0.521 / sin 10° = 3.0 tolerances), the bound its documentation gives
/// a valid fold. And a spike of the notch probe's random sweep (case 1165
/// of its sweep 6) whose A ∩ B, on main, validated and passed the
/// self-contact oracle while the spike's face and the piece of A's top
/// over it crossed 11.1 tolerances inside both, found by this check: B, a
/// prism on the plane through the x axis tilted 1°, its bottom face a
/// strip with a spike up to (−0.05, 57.3 ε), a tolerance above A's top.
/// The boolean is now refused at the spike's tip (design note 0012 §21).
#[test]
fn faces_cross_measures_a_fold_and_finds_no_spike() {
    let eps = k().tolerance().linear;
    let long = cuboid(DVec3::new(-0.25, -1.0, -1.0), DVec3::new(1.9, 1.0, 0.0));
    let far = tilted(
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
    // The corner, 3 tolerances across the line in B's plane, is 3 sin 10°
    // off A's top, and so 3 tolerances from the planes' common line.
    let bound = 3.0;
    for (x, y) in [(&long, &far), (&far, &long)] {
        let r = k().boolean(x, y, BooleanOp::Intersect).expect("a result");
        let c = self_contact::faces_crossing(&k(), &r, 1e-8)
            .expect("a mesh")
            .expect("a crossing");
        assert!((2.9..bound).contains(&(c.at.depth / eps)), "{c}");
        assert!((c.at.angle.to_degrees() - 10.0).abs() < 1e-6, "{c}");
        assert_eq!(c.shared, self_contact::Shared::Edge, "{c}");
        println!("{c}");
    }
    let a = cuboid(DVec3::new(-0.25, -1.0, -1.0), DVec3::new(1.25, 1.0, 0.0));
    let pts = [
        (-0.5, -0.5),
        (2.0, -0.5),
        (2.0, -0.4),
        (0.02914404710520814, -0.4),
        (-0.05, 5.729868849855019e-6),
        (-0.12914404710520816, -0.4),
        (-0.5, -0.4),
    ]
    .map(|(u, w)| DVec2::new(u, w));
    let b = prism(tilted_frame(1.0), &pts, 0.5);
    for (x, y) in [(&a, &b), (&b, &a)] {
        let e = k()
            .boolean(x, y, BooleanOp::Intersect)
            .expect_err("refused");
        assert!(
            e.to_string()
                .contains("runs within the tolerance of the other's surface"),
            "{e}"
        );
    }
}
