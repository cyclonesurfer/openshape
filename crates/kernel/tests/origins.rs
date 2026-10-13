//! Origin tables (ADR 0018, design note 0013 §9): complete, unique,
//! consistent with the body's topology, and right by an independent
//! geometric check, for boxes, extrudes and revolves.

use std::collections::{BTreeSet, HashSet};
use std::f64::consts::{FRAC_PI_2, PI, TAU};

use opendatum_kernel::{
    Body, BrepKernel, DVec2, DVec3, EdgeOrigin, ExtrudeSpec, FaceIx, FaceOrigin, Frame, Kernel,
    Origins, PieceEnd, Profile, ProfileCurve, ProfilePiece, ProfileVertex, RevolveSpec,
    SurfaceKind, TessellationOptions, VertexOrigin,
};

fn k() -> BrepKernel {
    BrepKernel::default()
}

#[derive(Clone, Copy)]
enum Sweep {
    Extrude(ExtrudeSpec),
    Revolve(RevolveSpec),
}

/// Pieces the kernel splits an arc of this sweep into (at most π each).
fn pieces(sweep: f64) -> usize {
    ((sweep.abs() / PI) - 1e-9).ceil().max(1.0) as usize
}

/// The point at fraction `f` along a piece of a line, arc or circle.
fn piece_point(c: &ProfileCurve, piece: u32, f: f64) -> DVec2 {
    let on_arc = |center: DVec2, radius: f64, start: f64, sweep: f64| {
        let a = start + sweep * (piece as f64 + f) / pieces(sweep) as f64;
        center + DVec2::from_angle(a) * radius
    };
    match *c {
        ProfileCurve::Line { start, end } => start + (end - start) * f,
        ProfileCurve::Arc {
            center,
            radius,
            start_angle,
            sweep_angle,
        } => on_arc(center, radius, start_angle, sweep_angle),
        ProfileCurve::Circle { center, radius } => on_arc(center, radius, 0.0, TAU),
        ProfileCurve::Spline { .. } => unreachable!("test profiles have no splines"),
    }
}

fn mid(p: &Profile, piece: ProfilePiece) -> DVec2 {
    piece_point(&p.curves[piece.curve as usize], piece.piece, 0.5)
}

fn end_point(p: &Profile, e: PieceEnd) -> DVec2 {
    let f = if e.end { 1.0 } else { 0.0 };
    piece_point(&p.curves[e.piece.curve as usize], e.piece.piece, f)
}

/// A profile point carried to sweep position `t` (fractional between
/// positions).
fn place(p: &Profile, sweep: Sweep, positions: u32, q: DVec2, t: f64) -> DVec3 {
    let x = p.plane.point2(q);
    match sweep {
        Sweep::Extrude(s) => {
            let base = if s.symmetric { -0.5 } else { 0.0 };
            x + p.plane.z * s.distance * (base + t)
        }
        Sweep::Revolve(s) => {
            let o = p.plane.point2(s.axis_origin);
            let axis = p.plane.vector2(s.axis_dir).normalize();
            let a = s.angle * t / (positions - 1) as f64;
            o + glam::DQuat::from_axis_angle(axis, a) * (x - o)
        }
    }
}

fn tri_distance(p: DVec3, a: DVec3, b: DVec3, c: DVec3) -> f64 {
    // Closest point on a triangle (Ericson, Real-Time Collision Detection
    // 5.1.5), by its barycentric regions.
    let (ab, ac, ap) = (b - a, c - a, p - a);
    let (d1, d2) = (ab.dot(ap), ac.dot(ap));
    if d1 <= 0.0 && d2 <= 0.0 {
        return (p - a).length();
    }
    let bp = p - b;
    let (d3, d4) = (ab.dot(bp), ac.dot(bp));
    if d3 >= 0.0 && d4 <= d3 {
        return (p - b).length();
    }
    let vc = d1 * d4 - d3 * d2;
    if vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0 {
        return (p - (a + ab * (d1 / (d1 - d3)))).length();
    }
    let cp = p - c;
    let (d5, d6) = (ab.dot(cp), ac.dot(cp));
    if d6 >= 0.0 && d5 <= d6 {
        return (p - c).length();
    }
    let vb = d5 * d2 - d1 * d6;
    if vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0 {
        return (p - (a + ac * (d2 / (d2 - d6)))).length();
    }
    let va = d3 * d6 - d5 * d4;
    if va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0 {
        let w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        return (p - (b + (c - b) * w)).length();
    }
    let denom = 1.0 / (va + vb + vc);
    (p - (a + ab * (vb * denom) + ac * (vc * denom))).length()
}

fn segment_distance(p: DVec3, a: DVec3, b: DVec3) -> f64 {
    let ab = b - a;
    let t = ((p - a).dot(ab) / ab.length_squared()).clamp(0.0, 1.0);
    (p - (a + ab * t)).length()
}

/// Index and distance of the nearest of `items` to `p` by `dist`.
fn nearest<T>(items: &[T], dist: impl Fn(&T) -> f64) -> (usize, f64) {
    items
        .iter()
        .enumerate()
        .map(|(i, x)| (i, dist(x)))
        .min_by(|a, b| a.1.total_cmp(&b.1))
        .unwrap()
}

/// The pieces of the loop holding `piece`, in order, and its position.
fn loop_of(o: &Origins, piece: ProfilePiece) -> (&[ProfilePiece], usize) {
    for region in &o.regions {
        for lp in region {
            if let Some(i) = lp.iter().position(|&p| p == piece) {
                return (lp, i);
            }
        }
    }
    panic!("{piece:?} is in no region");
}

/// The pieces next to `piece` in its loop.
fn neighbours(o: &Origins, piece: ProfilePiece) -> [ProfilePiece; 2] {
    let (lp, i) = loop_of(o, piece);
    let n = lp.len();
    [lp[(i + n - 1) % n], lp[(i + 1) % n]]
}

fn region_of(o: &Origins, piece: ProfilePiece) -> u32 {
    o.regions
        .iter()
        .position(|r| r.iter().flatten().any(|&p| p == piece))
        .expect("piece in a region") as u32
}

/// Every check of ADR 0018's consequences on one body.
/// Every planar face has a plane (and no other face has one) that its
/// mesh lies on, with the normal of its outward-wound triangles.
fn check_planes(b: &Body) {
    let kern = k();
    let mesh = kern.tessellate(b, &TessellationOptions::default()).unwrap();
    let kinds = kern.face_kinds(b);
    for (i, fm) in mesh.faces.iter().enumerate() {
        let plane = kern.face_plane(b, FaceIx(i as u32));
        let Some(plane) = plane else {
            assert_ne!(kinds[i], SurfaceKind::Plane, "face {i} has no plane");
            continue;
        };
        assert_eq!(kinds[i], SurfaceKind::Plane, "face {i}");
        assert!((plane.normal.length() - 1.0).abs() < 1e-12);
        let mut area = DVec3::ZERO;
        for t in &fm.triangles {
            let [a, b, c] = t.map(|j| fm.positions[j as usize]);
            area += (b - a).cross(c - a);
        }
        for q in &fm.positions {
            let off = (*q - plane.point).dot(plane.normal).abs();
            assert!(off < 1e-9, "face {i}: a mesh point is {off} off its plane");
        }
        let area_n = area.normalize();
        assert!(
            (area_n - plane.normal).length() < 1e-9,
            "face {i}: outward {area_n} vs plane normal {}",
            plane.normal
        );
    }
    assert_eq!(kern.face_plane(b, FaceIx(mesh.faces.len() as u32)), None);
}

fn check(p: &Profile, sweep: Sweep, b: &Body) {
    check_planes(b);
    let kern = k();
    let o = kern.origins(b);
    let t = kern.topology(b);
    let full = matches!(sweep, Sweep::Revolve(s) if s.angle.abs() >= TAU - 1e-9);
    let last = o.positions - 1;

    // Complete and unique.
    assert_eq!(
        (o.faces.len(), o.edges.len(), o.vertices.len()),
        (t.faces, t.edges, t.vertices),
        "one origin per entity"
    );
    assert_eq!(o.faces.iter().collect::<HashSet<_>>().len(), o.faces.len());
    assert_eq!(o.edges.iter().collect::<HashSet<_>>().len(), o.edges.len());
    assert_eq!(
        o.vertices.iter().collect::<HashSet<_>>().len(),
        o.vertices.len()
    );
    let mut all_pieces: Vec<ProfilePiece> = o.regions.iter().flatten().flatten().copied().collect();
    let n_pieces = all_pieces.len();
    all_pieces.sort();
    all_pieces.dedup();
    assert_eq!(all_pieces.len(), n_pieces, "each piece in one loop once");

    // Consistent with topology.
    let adj = kern.adjacency(b);
    let face_set = |fs: &[FaceIx]| -> BTreeSet<FaceOrigin> {
        fs.iter().map(|f| o.faces[f.0 as usize]).collect()
    };
    let span_wrap = |s: i64| -> Option<u32> {
        let k = last as i64;
        if full {
            Some(s.rem_euclid(k) as u32)
        } else {
            (0..k).contains(&s).then_some(s as u32)
        }
    };
    let side_spans = |piece: ProfilePiece| -> BTreeSet<u32> {
        o.faces
            .iter()
            .filter_map(|f| match *f {
                FaceOrigin::Side { piece: q, span } if q == piece => Some(span),
                _ => None,
            })
            .collect()
    };
    for (e, origin) in o.edges.iter().enumerate() {
        let got = face_set(&adj.edge_faces[e]);
        let want: BTreeSet<FaceOrigin> = match *origin {
            EdgeOrigin::Profile { piece, at } => {
                let region = region_of(o, piece);
                if side_spans(piece).is_empty() {
                    // A piece on a revolve's axis bounds both caps.
                    [0, last]
                        .map(|at| FaceOrigin::Cap { at, region })
                        .into_iter()
                        .collect()
                } else {
                    let mut w: BTreeSet<FaceOrigin> = [at as i64 - 1, at as i64]
                        .into_iter()
                        .filter_map(span_wrap)
                        .map(|span| FaceOrigin::Side { piece, span })
                        .collect();
                    if !full && (at == 0 || at == last) {
                        w.insert(FaceOrigin::Cap { at, region });
                    }
                    w
                }
            }
            EdgeOrigin::Swept { vertex, span } => vertex
                .ends
                .map(|e| FaceOrigin::Side {
                    piece: e.piece,
                    span,
                })
                .into_iter()
                .collect(),
            _ => panic!("unexpected origin {origin:?}"),
        };
        assert_eq!(got, want, "faces of edge {e} ({origin:?})");
    }
    for (v, origin) in o.vertices.iter().enumerate() {
        let VertexOrigin::Profile { vertex, at } = *origin else {
            panic!("unexpected origin {origin:?}")
        };
        let pieces = vertex.ends.map(|e| e.piece);
        // The two piece ends that meet there, sorted, of neighbouring
        // pieces.
        let x = end_point(p, vertex.ends[0]);
        let mut meeting: Vec<PieceEnd> = all_pieces
            .iter()
            .flat_map(|&piece| [false, true].map(|end| PieceEnd { piece, end }))
            .filter(|&e| (end_point(p, e) - x).length() < 1e-9)
            .collect();
        meeting.sort();
        assert_eq!(meeting, vertex.ends, "vertex {v} ({origin:?})");
        assert!(
            neighbours(o, pieces[0]).contains(&pieces[1]),
            "vertex {v} ({origin:?})"
        );
        let edges = &adj.vertex_edges[v];
        assert!(edges.len() >= 2, "vertex {v} {origin:?}");
        // A vertex on a revolve's axis is one vertex for every position: it
        // sweeps no edge.
        let on_axis = edges
            .iter()
            .all(|e| matches!(o.edges[e.0 as usize], EdgeOrigin::Profile { .. }));
        assert!(!on_axis || at == 0, "vertex {v} ({origin:?}) on the axis");
        for e in edges {
            let ok = match o.edges[e.0 as usize] {
                EdgeOrigin::Profile { piece, at: a } => {
                    pieces.contains(&piece) && (a == at || on_axis || side_spans(piece).is_empty())
                }
                EdgeOrigin::Swept { vertex: w, span } => {
                    w == vertex
                        && [at as i64 - 1, at as i64]
                            .into_iter()
                            .filter_map(span_wrap)
                            .any(|s| s == span)
                }
                other => panic!("unexpected origin {other:?}"),
            };
            assert!(
                ok,
                "vertex {v} ({origin:?}) on edge {e:?} ({:?})",
                o.edges[e.0 as usize]
            );
        }
    }

    // Right by geometry: a side face is the one nearest its piece swept to
    // mid-span, a profile edge the one nearest its piece at its position, a
    // swept edge runs between its vertex's two positions, and a cap lies in
    // the plane of its position.
    let opts = TessellationOptions {
        chord_tolerance: 1e-3,
        max_angle: 0.3,
    };
    let mesh = kern.tessellate(b, &opts).unwrap();
    let close = 2e-3;
    for (f, origin) in o.faces.iter().enumerate() {
        match *origin {
            FaceOrigin::Side { piece, span } => {
                let q = place(p, sweep, o.positions, mid(p, piece), span as f64 + 0.5);
                let (i, d) = nearest(&mesh.faces, |m| {
                    m.triangles
                        .iter()
                        .map(|t| {
                            let [a, b, c] = t.map(|i| m.positions[i as usize]);
                            tri_distance(q, a, b, c)
                        })
                        .fold(f64::INFINITY, f64::min)
                });
                assert!(
                    i == f && d < close,
                    "face {f} ({origin:?}): nearest {i} at {d}"
                );
            }
            FaceOrigin::Cap { at, .. } => {
                let o3 = place(p, sweep, o.positions, DVec2::ZERO, at as f64);
                let ex = place(p, sweep, o.positions, DVec2::X, at as f64) - o3;
                let ey = place(p, sweep, o.positions, DVec2::Y, at as f64) - o3;
                let n = ex.cross(ey).normalize();
                for q in &mesh.faces[f].positions {
                    assert!((*q - o3).dot(n).abs() < 1e-9, "cap {f} ({origin:?})");
                }
            }
            _ => panic!("unexpected origin {origin:?}"),
        }
    }
    for (e, origin) in o.edges.iter().enumerate() {
        match *origin {
            EdgeOrigin::Profile { piece, at } => {
                let q = place(p, sweep, o.positions, mid(p, piece), at as f64);
                let (i, d) = nearest(&mesh.edges, |l| {
                    l.windows(2)
                        .map(|w| segment_distance(q, w[0], w[1]))
                        .fold(f64::INFINITY, f64::min)
                });
                assert!(
                    i == e && d < close,
                    "edge {e} ({origin:?}): nearest {i} at {d}"
                );
            }
            EdgeOrigin::Swept { vertex, span } => {
                let x = end_point(p, vertex.ends[0]);
                let a = place(p, sweep, o.positions, x, span as f64);
                let b = place(p, sweep, o.positions, x, span as f64 + 1.0);
                let l = &mesh.edges[e];
                let (s, t) = (l[0], *l.last().unwrap());
                let ok = ((s - a).length() < 1e-9 && (t - b).length() < 1e-9)
                    || ((s - b).length() < 1e-9 && (t - a).length() < 1e-9);
                assert!(ok, "edge {e} ({origin:?}) runs {s} to {t}, not {a} to {b}");
            }
            _ => panic!("unexpected origin {origin:?}"),
        }
    }
}

fn extrude(p: &Profile, distance: f64, symmetric: bool) {
    let spec = ExtrudeSpec {
        distance,
        symmetric,
    };
    let b = k().extrude(p, &spec).unwrap();
    assert_eq!(k().origins(&b).positions, 2);
    check(p, Sweep::Extrude(spec), &b);
}

fn revolve(p: &Profile, axis_origin: DVec2, axis_dir: DVec2, angle: f64) {
    let spec = RevolveSpec {
        axis_origin,
        axis_dir,
        angle,
    };
    let b = k()
        .revolve(p, &spec)
        .unwrap_or_else(|e| panic!("angle {angle}: {e}"));
    check(p, Sweep::Revolve(spec), &b);
}

fn tilted() -> Frame {
    let z = DVec3::new(0.3, -0.4, 0.866).normalize();
    let x = DVec3::Y.cross(z).normalize();
    Frame {
        origin: DVec3::new(0.5, -1.0, 2.0),
        x,
        y: z.cross(x),
        z,
    }
}

/// A slot (two lines, two half-circle arcs) with a round hole, and a
/// separate rectangle: two regions, a hole, arcs of exactly π and a circle
/// split in two.
fn slot_and_square(plane: Frame) -> Profile {
    let (r, l) = (1.0, 3.0);
    Profile {
        plane,
        curves: vec![
            ProfileCurve::Line {
                start: DVec2::new(0.0, -r),
                end: DVec2::new(l, -r),
            },
            ProfileCurve::Arc {
                center: DVec2::new(l, 0.0),
                radius: r,
                start_angle: -FRAC_PI_2,
                sweep_angle: PI,
            },
            // Drawn against the loop's direction.
            ProfileCurve::Line {
                start: DVec2::new(0.0, r),
                end: DVec2::new(l, r),
            },
            ProfileCurve::Arc {
                center: DVec2::ZERO,
                radius: r,
                start_angle: FRAC_PI_2,
                sweep_angle: PI,
            },
            ProfileCurve::Circle {
                center: DVec2::new(1.5, 0.0),
                radius: 0.5,
            },
        ],
    }
    .with(Profile::rectangle(
        plane,
        DVec2::new(6.0, -1.0),
        DVec2::new(8.0, 1.0),
    ))
}

#[test]
fn boxes_have_complete_origins() {
    for (frame, size) in [
        (Frame::WORLD, DVec3::new(1.0, 2.0, 3.0)),
        (tilted(), DVec3::new(-1.0, 0.5, -2.0)),
    ] {
        let b = k().make_box(&frame, size).unwrap();
        let p = Profile::rectangle(
            frame,
            DVec2::ZERO.min(size.truncate()),
            DVec2::ZERO.max(size.truncate()),
        );
        let spec = ExtrudeSpec {
            distance: size.z,
            symmetric: false,
        };
        check(&p, Sweep::Extrude(spec), &b);
        let o = k().origins(&b);
        assert_eq!(o.regions.len(), 1);
        assert_eq!(o.regions[0].len(), 1, "one loop");
        assert_eq!(o.regions[0][0].len(), 4, "four sides");
    }
}

#[test]
fn extrudes_have_complete_origins() {
    for plane in [Frame::WORLD, tilted()] {
        let p = slot_and_square(plane);
        extrude(&p, 0.75, false);
        extrude(&p, -0.75, true);
        let b = k()
            .extrude(
                &p,
                &ExtrudeSpec {
                    distance: 1.0,
                    symmetric: false,
                },
            )
            .unwrap();
        let o = k().origins(&b);
        assert_eq!(o.regions.len(), 2);
        // The slot: outer loop of 4 pieces, the hole of 2 (a circle).
        let slot = o
            .regions
            .iter()
            .find(|r| r.len() == 2)
            .expect("the slot has a hole");
        assert_eq!((slot[0].len(), slot[1].len()), (4, 2));
        assert!(slot[1].iter().all(|q| q.curve == 4));
    }
    // An arc of more than π is two pieces, each with its own side face.
    let p = Profile {
        plane: Frame::WORLD,
        curves: vec![
            ProfileCurve::Arc {
                center: DVec2::ZERO,
                radius: 1.0,
                start_angle: -0.25 * PI,
                sweep_angle: 1.5 * PI,
            },
            ProfileCurve::Line {
                start: DVec2::from_angle(1.25 * PI),
                end: DVec2::ZERO,
            },
            ProfileCurve::Line {
                start: DVec2::ZERO,
                end: DVec2::from_angle(-0.25 * PI),
            },
        ],
    };
    extrude(&p, 0.5, false);
    let b = k()
        .extrude(
            &p,
            &ExtrudeSpec {
                distance: 0.5,
                symmetric: false,
            },
        )
        .unwrap();
    let sides = k()
        .origins(&b)
        .faces
        .iter()
        .filter(|f| matches!(f, FaceOrigin::Side { piece, .. } if piece.curve == 0))
        .count();
    assert_eq!(sides, 2);
}

#[test]
fn revolves_have_complete_origins() {
    // A rectangle off the axis: partial revolves of one and two pieces,
    // and a full revolve (a tube).
    let rect = Profile::rectangle(Frame::WORLD, DVec2::new(1.0, 0.0), DVec2::new(2.0, 1.0));
    for angle in [0.5 * PI, 1.5 * PI, -1.2 * PI, TAU] {
        revolve(&rect, DVec2::ZERO, DVec2::Y, angle);
    }
    // A half disc on the axis: its diameter lies on the axis (one edge
    // bounding both caps when partial; no edge when full), its ends are
    // vertices on the axis.
    let half = Profile {
        plane: Frame::WORLD,
        curves: vec![
            ProfileCurve::Arc {
                center: DVec2::ZERO,
                radius: 1.0,
                start_angle: -FRAC_PI_2,
                sweep_angle: PI,
            },
            ProfileCurve::Line {
                start: DVec2::new(0.0, 1.0),
                end: DVec2::new(0.0, -1.0),
            },
        ],
    };
    for angle in [0.75 * PI, TAU] {
        revolve(&half, DVec2::ZERO, DVec2::Y, angle);
    }
    // A slot and a square on a tilted plane about an axis below them.
    let p = slot_and_square(tilted());
    for angle in [0.4 * PI, TAU] {
        revolve(&p, DVec2::new(0.0, -3.0), DVec2::X, angle);
    }
}

#[test]
fn origins_follow_the_curves_not_their_order() {
    // The same slot with its curves listed in another order: every origin
    // is the same up to the renumbering of curves.
    let p = slot_and_square(Frame::WORLD);
    let perm = [3usize, 0, 4, 1, 2, 6, 5, 8, 7];
    let q = Profile {
        plane: p.plane,
        curves: perm.iter().map(|&i| p.curves[i].clone()).collect(),
    };
    let spec = ExtrudeSpec {
        distance: 1.0,
        symmetric: false,
    };
    let (a, b) = (
        k().extrude(&p, &spec).unwrap(),
        k().extrude(&q, &spec).unwrap(),
    );
    let renumber = |x: ProfilePiece| ProfilePiece {
        curve: perm[x.curve as usize] as u32,
        piece: x.piece,
    };
    let sides = |b: &Body, map: &dyn Fn(ProfilePiece) -> ProfilePiece| -> BTreeSet<ProfilePiece> {
        k().origins(b)
            .faces
            .iter()
            .filter_map(|f| match *f {
                FaceOrigin::Side { piece, .. } => Some(map(piece)),
                _ => None,
            })
            .collect()
    };
    assert_eq!(sides(&a, &|x| x), sides(&b, &renumber));
    let vertices =
        |b: &Body, map: &dyn Fn(ProfilePiece) -> ProfilePiece| -> BTreeSet<(ProfileVertex, u32)> {
            k().origins(b)
                .vertices
                .iter()
                .map(|v| match *v {
                    VertexOrigin::Profile { vertex, at } => {
                        let mut ends = vertex.ends.map(|e| PieceEnd {
                            piece: map(e.piece),
                            end: e.end,
                        });
                        ends.sort();
                        (ProfileVertex { ends }, at)
                    }
                    _ => unreachable!(),
                })
                .collect()
        };
    assert_eq!(vertices(&a, &|x| x), vertices(&b, &renumber));
}
