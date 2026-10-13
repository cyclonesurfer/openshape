//! Booleans of solids: design note 0012 and ADR 0020.
//!
//! The stages of design note 0012 §4, for polyhedra (§9) and, as step 2
//! lands (§10), for solids bounded by planes and cylinders:
//!
//! 1. **Points** (`points`): every vertex of either body, every point where
//!    an edge of one crosses a face's surface of the other, and every point
//!    where two edges touch, clustered within the tolerance into the
//!    candidate vertices of the result. Whether a point lies on an edge is
//!    decided once per point and edge, and every later stage reads that.
//! 2. **Imprint** (this module): each pair of faces that meet adds the
//!    segments along which they meet to both faces: for faces on two
//!    planes, the stretch of their common line inside both (closed) faces,
//!    found from the signs of the vertices against the other plane; for
//!    faces on one plane, the pieces of each face's edges inside the other.
//! 3. **Split** (`planar`): each face's edges and imprinted segments form a
//!    planar graph on known points, whose bounded regions inside the face
//!    are its pieces.
//! 4. **Classify** (`classify`): each piece is on a coincident face of the
//!    other body (with the same or the opposite outward normal), or inside
//!    or outside it by exact ray parity at a point well inside the piece.
//! 5. **Select, sew, merge** (`assemble`): pieces are kept by the table of
//!    §4.4, sewn along their segments, coplanar neighbours and collinear
//!    edges joined, shells and solids rebuilt, and the origins of ADR 0020
//!    recorded. The kernel's validator then checks the result.
//!
//! Segments are known by their index rather than by their ends, since two
//! curved segments (two halves of a circle) can join the same two points.
//!
//! Anything the stages cannot decide within the tolerance is refused with
//! a message and a location rather than guessed.

use std::f64::consts::{PI, TAU};

use glam::{DVec2, DVec3};
use rustc_hash::{FxHashMap, FxHashSet};

use crate::bvh::Bvh;
use crate::origin::{Operand, Origins};
use crate::topo::Brep;
use crate::{BooleanOp, KernelError, Tolerance};

mod assemble;
mod bernstein;
mod chart;
mod classify;
mod imprint;
mod meet;
mod planar;
mod points;
mod quadric;
mod shape;
#[cfg(test)]
mod testing;

use chart::Chart;
use classify::Class;
use imprint::Relation;
use meet::Track;
use points::Points;
use shape::Shape;

/// Work limits (docs/PERFORMANCE.md): pairs of faces that meet, pairs of
/// edges near each other, candidate points, and segments of the result's
/// edges. Every stage takes time in proportion to these (times a
/// logarithm), so the limits bound a boolean's time and memory whatever its
/// input. The defaults are far above any part a person models (a result of
/// a million faces); a boolean past one is refused rather than left to run
/// for minutes.
#[derive(Clone, Copy, Debug)]
struct Limits {
    face_pairs: usize,
    edge_pairs: usize,
    points: usize,
    segments: usize,
}

const LIMITS: Limits = Limits {
    face_pairs: 2_000_000,
    edge_pairs: 8_000_000,
    points: 4_000_000,
    segments: 8_000_000,
};

fn too_large(what: &str, limit: usize) -> KernelError {
    KernelError::TooLarge(format!(
        "the boolean would need more than {limit} {what}; it is refused rather than left \
         to run for minutes"
    ))
}

/// The two operands, indexed by [`Operand`].
struct Pair<'a> {
    first: Shape<'a>,
    second: Shape<'a>,
}

impl<'a> Pair<'a> {
    fn get(&self, op: Operand) -> &Shape<'a> {
        match op {
            Operand::First => &self.first,
            Operand::Second => &self.second,
        }
    }
}

fn other(op: Operand) -> Operand {
    match op {
        Operand::First => Operand::Second,
        Operand::Second => Operand::First,
    }
}

fn boxes_overlap(a: (DVec3, DVec3), b: (DVec3, DVec3), eps: f64) -> bool {
    a.0.cmple(b.1 + DVec3::splat(eps)).all() && b.0.cmple(a.1 + DVec3::splat(eps)).all()
}

/// A box grown by this finds, in a tree, every box that `boxes_overlap`
/// accepts at `eps`: twice the tolerance, so the rounding of the two tests
/// (far below the tolerance wherever coordinates are below `eps · 2⁵¹`,
/// hundreds of kilometres at the default) cannot drop one.
fn pad(eps: f64) -> DVec3 {
    DVec3::splat(2.0 * eps)
}

/// Two clusters, the lesser first: the ends of a segment.
type Key = (usize, usize);

fn key(a: usize, b: usize) -> Key {
    (a.min(b), a.max(b))
}

/// Where a segment came from.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
enum Source {
    /// A piece of an operand's edge.
    Edge(Operand, usize),
    /// Where face `.0` of the first operand meets face `.1` of the second.
    Meet(usize, usize),
}

/// The course of a segment between its two clusters.
#[derive(Clone, Copy, Debug, PartialEq)]
enum Course {
    /// The straight line between the clusters' places.
    Straight,
    /// Along curved edge `e` of operand `op`, from its parameter `ta` at
    /// the segment's `a` to `tb` at its `b`.
    Edge {
        op: Operand,
        e: usize,
        ta: f64,
        tb: f64,
    },
    /// Along closed track `id`, from parameter `ta` at `a` to `tb` at `b`,
    /// less than a turn apart.
    Track { id: usize, ta: f64, tb: f64 },
}

impl Course {
    /// The course run the other way.
    fn reversed(self) -> Course {
        match self {
            Course::Straight => Course::Straight,
            Course::Edge { op, e, ta, tb } => Course::Edge {
                op,
                e,
                ta: tb,
                tb: ta,
            },
            Course::Track { id, ta, tb } => Course::Track { id, ta: tb, tb: ta },
        }
    }
}

/// A segment of the result's candidate edges.
struct Seg {
    /// Its ends (clusters), the lesser first.
    a: usize,
    b: usize,
    course: Course,
    sources: Vec<Source>,
    /// For a curved segment, the point of its course halfway along its
    /// parameter.
    mid: Option<DVec3>,
}

/// A segment in one face's graph: material on the left of `a → b`
/// (`Some(true)`), of `b → a` (`Some(false)`), or on both sides (`None`,
/// an imprint inside the face).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct FaceSeg {
    seg: usize,
    side: Option<bool>,
}

/// A piece of an operand face, as cycles of segments (outer first),
/// counter-clockwise about the face's outward normal for the outer loop,
/// each segment with whether it is run from its `a` to its `b`.
struct Piece {
    operand: Operand,
    face: usize,
    loops: Vec<Vec<(usize, bool)>>,
    class: Class,
    /// For a piece on a coincident face of the other operand, that face.
    partner: Option<usize>,
}

/// A run of segments as the clusters it starts at, in order.
fn corners(segs: &Segments, l: &[(usize, bool)]) -> Vec<usize> {
    l.iter().map(|&(s, fw)| segs.start(s, fw)).collect()
}

pub(crate) fn boolean(
    a: &Brep,
    b: &Brep,
    op: BooleanOp,
    tol: &Tolerance,
) -> Result<(Brep, Origins), KernelError> {
    boolean_within(a, b, op, tol, &LIMITS)
}

fn boolean_within(
    a: &Brep,
    b: &Brep,
    op: BooleanOp,
    tol: &Tolerance,
    limits: &Limits,
) -> Result<(Brep, Origins), KernelError> {
    let pair = Pair {
        first: shape::shape(a, "first")?,
        second: shape::shape(b, "second")?,
    };
    let eps = tol.linear;
    let outlines = [
        imprint::outlines(&pair.first)?,
        imprint::outlines(&pair.second)?,
    ];
    // Each face's chart: its plane's frame, or its outline's.
    let charts: [Vec<Chart>; 2] = [Operand::First, Operand::Second].map(|op| {
        pair.get(op)
            .faces
            .iter()
            .zip(&outlines[slot(op)])
            .map(|(f, o)| o.as_ref().map_or(Chart::Plane(f.frame), |o| o.chart))
            .collect()
    });
    let mut contacts = imprint::Contacts::new(&pair);
    let mut pts = Points::new(&pair, eps);

    // Face pairs that can meet, and how.
    let mut crossing: Vec<(usize, usize)> = Vec::new();
    let mut coplanar: Vec<(usize, usize)> = Vec::new();
    let tree = Bvh::new(pair.second.faces.iter().map(|f| (f.lo, f.hi)).collect());
    let mut near = Vec::new();
    for (fa, a) in pair.first.faces.iter().enumerate() {
        tree.overlapping(a.lo - pad(eps), a.hi + pad(eps), &mut near);
        for &fb in &near {
            let b = &pair.second.faces[fb];
            if !boxes_overlap((a.lo, a.hi), (b.lo, b.hi), eps) {
                continue;
            }
            if crossing.len() + coplanar.len() >= limits.face_pairs {
                return Err(too_large("pairs of faces that meet", limits.face_pairs));
            }
            match imprint::relation(&pair, &mut pts, &mut contacts, fa, fb)? {
                Relation::Apart => {}
                Relation::Coplanar => coplanar.push((fa, fb)),
                Relation::Crossing => {
                    pts.add_crossings(&pair, Operand::First, fa, b.surface)?;
                    pts.add_crossings(&pair, Operand::Second, fb, a.surface)?;
                    if pts.raw_count() > limits.points {
                        return Err(too_large("candidate points", limits.points));
                    }
                    crossing.push((fa, fb));
                }
            }
        }
    }
    pts.add_edge_contacts(&pair, limits.edge_pairs)
        .map_err(|limit| too_large("pairs of edges near each other", limit))?;
    if pts.raw_count() > limits.points {
        return Err(too_large("candidate points", limits.points));
    }
    pts.cluster()?;
    pts.incidences(&pair)?;

    // Each face's segments: its own edges' pieces, then imprints.
    let mut segs = Segments::default();
    let mut face_segs: [Vec<Vec<FaceSeg>>; 2] = [
        vec![Vec::new(); pair.first.faces.len()],
        vec![Vec::new(); pair.second.faces.len()],
    ];
    // Per operand, per edge, its pieces from its start to its end.
    let mut edge_segs: [Vec<Vec<usize>>; 2] = [
        vec![Vec::new(); pair.first.edges.len()],
        vec![Vec::new(); pair.second.edges.len()],
    ];
    for op in [Operand::First, Operand::Second] {
        let shape = pair.get(op);
        for (fi, f) in shape.faces.iter().enumerate() {
            for l in &f.loops {
                for &(e, forward) in l {
                    if edge_segs[slot(op)][e].is_empty() {
                        let along = pts.on_edge(op, e);
                        edge_segs[slot(op)][e] = if shape.straight[e] {
                            along
                                .windows(2)
                                .map(|w| {
                                    segs.add(w[0], w[1], Course::Straight, Source::Edge(op, e))
                                })
                                .collect()
                        } else {
                            let ts = pts.on_edge_params(op, e);
                            let (curve, _, _) = shape.curve(e);
                            (0..along.len() - 1)
                                .map(|k| {
                                    let (ta, tb) = (ts[k], ts[k + 1]);
                                    segs.add_curved(
                                        along[k],
                                        along[k + 1],
                                        Course::Edge { op, e, ta, tb },
                                        curve.point(0.5 * (ta + tb)),
                                        Source::Edge(op, e),
                                        eps,
                                    )
                                })
                                .collect()
                        };
                    }
                    let along = pts.on_edge(op, e);
                    for (k, &s) in edge_segs[slot(op)][e].iter().enumerate() {
                        let from = if forward { along[k] } else { along[k + 1] };
                        face_segs[slot(op)][fi].push(FaceSeg {
                            seg: s,
                            side: Some(from == segs.list[s].a),
                        });
                    }
                }
            }
        }
    }
    // Closed tracks wholly in both faces of a pair, no edge reaching them.
    let mut whole: Vec<(usize, usize, usize)> = Vec::new();
    for &(fa, fb) in &crossing {
        let (a, b) = (&pair.first.faces[fa], &pair.second.faces[fb]);
        let planes = matches!(
            pair.first.surfaces[a.surface],
            quadric::Quadric::Plane { .. }
        ) && matches!(
            pair.second.surfaces[b.surface],
            quadric::Quadric::Plane { .. }
        );
        let mut imprint = |s: usize| {
            face_segs[0][fa].push(FaceSeg { seg: s, side: None });
            face_segs[1][fb].push(FaceSeg { seg: s, side: None });
        };
        if planes {
            for (lo, hi) in imprint::plane_meeting(&pair, &pts, &outlines, fa, fb)? {
                for w in pts.split_segment(lo, hi)?.windows(2) {
                    imprint(segs.add(w[0], w[1], Course::Straight, Source::Meet(fa, fb)));
                }
            }
            continue;
        }
        let (parts, full) = imprint::stretches(&pair, &pts, &outlines, &contacts, fa, fb)?;
        whole.extend(full.into_iter().map(|id| (fa, fb, id)));
        for st in parts {
            let track = contacts.tracks[st.track];
            if !track.closed() {
                for w in pts.split_segment(st.lo.0, st.hi.0)?.windows(2) {
                    imprint(segs.add(w[0], w[1], Course::Straight, Source::Meet(fa, fb)));
                }
                continue;
            }
            for w in pts.on_track(&track, st.lo, st.hi)?.windows(2) {
                let (ta, tb) = (w[0].1, w[1].1);
                imprint(segs.add_curved(
                    w[0].0,
                    w[1].0,
                    Course::Track {
                        id: st.track,
                        ta,
                        tb,
                    },
                    track.point(0.5 * (ta + tb)),
                    Source::Meet(fa, fb),
                    eps,
                ));
            }
        }
    }
    // Each such track is split at two seams, at its parameters 0 and π,
    // into halves with ends: new points, clear of every other (no edge of
    // either face comes within the tolerance of the track, which the split
    // of each face checks).
    let mut seams: Vec<(usize, usize)> = Vec::new();
    for &(_, _, id) in &whole {
        if !seams.iter().any(|&(t, _)| t == id) {
            seams.push((id, 0));
        }
    }
    let at: Vec<DVec3> = seams
        .iter()
        .flat_map(|&(id, _)| {
            [
                contacts.tracks[id].point(0.0),
                contacts.tracks[id].point(PI),
            ]
        })
        .collect();
    let first = pts.add_seams(&at);
    for (k, s) in seams.iter_mut().enumerate() {
        s.1 = first + 2 * k;
    }
    for &(fa, fb, id) in &whole {
        let track = contacts.tracks[id];
        let c = seams
            .iter()
            .find(|&&(t, _)| t == id)
            .map(|&(_, c)| c)
            .expect("every whole track has seams");
        for (lo, hi) in [((c, 0.0), (c + 1, PI)), ((c + 1, PI), (c, TAU))] {
            for w in pts.on_track(&track, lo, hi)?.windows(2) {
                let (ta, tb) = (w[0].1, w[1].1);
                let s = segs.add_curved(
                    w[0].0,
                    w[1].0,
                    Course::Track { id, ta, tb },
                    track.point(0.5 * (ta + tb)),
                    Source::Meet(fa, fb),
                    eps,
                );
                face_segs[0][fa].push(FaceSeg { seg: s, side: None });
                face_segs[1][fb].push(FaceSeg { seg: s, side: None });
            }
        }
    }
    for &(fa, fb) in &coplanar {
        for (op, f, g) in [(Operand::First, fa, fb), (Operand::Second, fb, fa)] {
            // The pieces of g's edges in f (closed).
            let shape = pair.get(op);
            let face = &shape.faces[f];
            let gs = pair.get(other(op));
            for l in &gs.faces[g].loops {
                for &(e, _) in l {
                    for &s in &edge_segs[slot(other(op))][e] {
                        let seg = &segs.list[s];
                        let m = seg
                            .mid
                            .unwrap_or_else(|| (pts.rep[seg.a] + pts.rep[seg.b]) * 0.5);
                        let outline = outlines[slot(op)][f].as_ref();
                        if imprint::locate(shape, outline, face, m, eps) != Some(false) {
                            face_segs[slot(op)][f].push(FaceSeg { seg: s, side: None });
                        }
                    }
                }
            }
        }
    }

    if segs.list.len() > limits.segments {
        return Err(too_large("edge segments", limits.segments));
    }

    // Split every face into pieces and classify them.
    let mut partners: [FxHashMap<usize, Vec<usize>>; 2] =
        [FxHashMap::default(), FxHashMap::default()];
    for &(fa, fb) in &coplanar {
        partners[0].entry(fa).or_default().push(fb);
        partners[1].entry(fb).or_default().push(fa);
    }
    let boundaries = [
        classify::Boundary::new(&pair.first, &outlines[0])?,
        classify::Boundary::new(&pair.second, &outlines[1])?,
    ];

    let boxes = [pair.first.extent(), pair.second.extent()];
    let mut pieces: Vec<Piece> = Vec::new();
    for op in [Operand::First, Operand::Second] {
        let shape = pair.get(op);
        let them = pair.get(other(op));
        for (fi, f) in shape.faces.iter().enumerate() {
            // A face farther than the tolerance from the other body's box
            // is untouched by it (no point where the bodies meet is on its
            // edges, and no face of the other body is on its surface) and
            // wholly outside it: its one piece, its own loops, needs no
            // point inside it and no test.
            let clear = !boxes_overlap((f.lo, f.hi), boxes[slot(other(op))], eps);
            let chart = charts[slot(op)][fi];
            let list = face_segs[slot(op)][fi].iter().map(|fs| (fs.seg, fs.side));
            let edges = planar_edges(&pair, &contacts.tracks, &segs, &chart, list)?;
            let pos = |c: usize| chart.map(pts.rep[c]);
            let reach = |q: DVec2| chart.reach(q, eps);
            let regions = planar::regions(&edges, &pos, &reach, eps).map_err(|p| match p {
                planar::Problem::Close(c) => KernelError::InvalidInput(format!(
                    "a point where the bodies meet lies within the tolerance of an edge \
                     it is not on, near {}",
                    show(pts.rep[c])
                )),
                planar::Problem::Cross(c) => KernelError::InvalidInput(format!(
                    "edges where the bodies meet pass within the tolerance of each other \
                     near {}",
                    show(pts.rep[c])
                )),
                planar::Problem::Inconsistent(m, c) => KernelError::InvalidResult(format!(
                    "boolean: face {fi} of the {} body: {m} near {}",
                    name(op),
                    show(pts.rep[c])
                )),
            })?;
            for region in regions {
                if clear {
                    pieces.push(Piece {
                        operand: op,
                        face: fi,
                        loops: region.loops,
                        class: Class::Outside,
                        partner: None,
                    });
                    continue;
                }
                let (inner, clearance) = interior(&region, &edges, &chart, eps)
                    .map_err(|m| KernelError::InvalidResult(format!("boolean: face {fi}: {m}")))?;
                let at = chart.unmap(inner);
                if clearance <= eps {
                    return Err(KernelError::InvalidInput(format!(
                        "the bodies leave a sliver thinner than the tolerance near {}",
                        show(at)
                    )));
                }
                let mut class = None;
                let mut partner = None;
                // On a coinciding face of the other body, or off it. Every
                // edge of that face inside this one bounds pieces, so a
                // piece is wholly on it or off it, and the point, more than
                // the tolerance inside the piece, is on the same side of
                // each of the face's edges as the piece: whether it is in
                // the face is decided by parity alone. A band of the
                // tolerance about the face's edges would take a piece off
                // the face, but within the tolerance of a corner of it, for
                // one on it.
                for &g in partners[slot(op)].get(&fi).into_iter().flatten() {
                    let gf = &them.faces[g];
                    let on = match &outlines[slot(other(op))][g] {
                        None => imprint::in_face(them, gf, at),
                        Some(o) => o.locate(o.chart.map(at), 0.0, eps / 16.0) == Some(true),
                    };
                    if on {
                        let same = if f.flat && gf.flat {
                            gf.normal().dot(f.normal()) > 0.0
                        } else {
                            outward(shape, f, at).dot(outward(them, gf, at)) > 0.0
                        };
                        class = Some(if same {
                            Class::OnSame
                        } else {
                            Class::OnOpposite
                        });
                        partner = Some(g);
                        break;
                    }
                }
                let class = match class {
                    Some(c) => c,
                    None => boundaries[slot(other(op))].classify(at, eps)?,
                };
                pieces.push(Piece {
                    operand: op,
                    face: fi,
                    loops: region.loops,
                    class,
                    partner,
                });
            }
        }
    }
    let coplanar: FxHashSet<(usize, usize)> = coplanar.into_iter().collect();
    let geometry = assemble::Geometry {
        pair: &pair,
        tracks: &contacts.tracks,
        charts: &charts,
    };
    assemble::assemble(&geometry, &pts, segs, pieces, &coplanar, op, tol)
}

/// The outward normal of face `f` of `shape` at `p` (on its surface).
fn outward(shape: &Shape, f: &shape::Face, p: DVec3) -> DVec3 {
    let g = shape.surfaces[f.surface].gradient(p);
    if f.outward {
        g
    } else {
        -g
    }
}

/// Points on a piece's boundary above which [`interior`] probes the piece
/// ([`classify::probe_point`]) before triangulating it.
const PROBE_ABOVE: usize = 64;

/// A point inside a region of a face and its distance in space from the
/// region's boundary at least: inside the polygon of the region's chords,
/// as far from its sides as can be found, less the chords' bands, as a
/// distance in space ([`Chart::space`]).
/// While that is within the tolerance, chords are halved, until their
/// sags are within a sixteenth of it.
fn interior(
    region: &planar::Region,
    edges: &[planar::Edge],
    chart: &Chart,
    eps: f64,
) -> Result<(DVec2, f64), String> {
    let mut refined: Option<planar::Region> = None;
    loop {
        let r = refined.as_ref().unwrap_or(region);
        let polygons = r.polygons();
        let band = r.band();
        // A large piece is probed first: a triangulation of it costs far
        // more, most of all when its points lie on circles.
        if polygons.iter().map(Vec::len).sum::<usize>() > PROBE_ABOVE {
            if let Some((inner, clearance)) = classify::probe_point(&polygons) {
                let space = chart.space(inner, clearance - band);
                if space > eps {
                    return Ok((inner, space));
                }
            }
        }
        let (inner, clearance) = classify::interior_point(&polygons)?;
        let space = chart.space(inner, clearance - band);
        let sag = r.sag();
        if space > eps || sag <= eps / 16.0 {
            return Ok((inner, space));
        }
        refined = Some(r.refined(edges, sag / 4.0));
    }
}

/// The segments `list` (with their sides) of a face whose chart is
/// `chart`, as edges of its planar graph: a curved one with its course in
/// the chart, unless that is straight there.
fn planar_edges<'a>(
    pair: &'a Pair,
    tracks: &[Track],
    segs: &Segments,
    chart: &Chart,
    list: impl Iterator<Item = (usize, Option<bool>)>,
) -> Result<Vec<planar::Edge<'a>>, KernelError> {
    let traced = |of: chart::TraceOf<'a>, ta: f64, tb: f64| {
        chart::trace(chart, of)?
            .map(|t| {
                let ts = t.division(ta, tb)?;
                Ok((Box::new(t) as Box<dyn planar::Bend + 'a>, ts))
            })
            .transpose()
    };
    list.map(|(s, side)| {
        let seg = &segs.list[s];
        let bend = match seg.course {
            Course::Straight => None,
            Course::Edge { op, e, ta, tb } => {
                traced(chart::TraceOf::Curve(pair.get(op).curve(e).0), ta, tb)?
            }
            Course::Track { id, ta, tb } => traced(chart::TraceOf::Track(tracks[id]), ta, tb)?,
        };
        Ok(planar::Edge {
            seg: s,
            a: seg.a,
            b: seg.b,
            side,
            bend,
        })
    })
    .collect()
}

fn slot(op: Operand) -> usize {
    match op {
        Operand::First => 0,
        Operand::Second => 1,
    }
}

fn name(op: Operand) -> &'static str {
    match op {
        Operand::First => "first",
        Operand::Second => "second",
    }
}

fn show(p: DVec3) -> String {
    format!("({:.6}, {:.6}, {:.6}) mm", p.x * 1e3, p.y * 1e3, p.z * 1e3)
}

/// The segments of the result's candidate edges, with their sources.
#[derive(Default)]
struct Segments {
    list: Vec<Seg>,
    /// The straight segment between two clusters, if there is one.
    straight: FxHashMap<Key, usize>,
    /// The curved segments between two clusters.
    curved: FxHashMap<Key, Vec<usize>>,
}

impl Segments {
    /// The segment from `from` to `to` along `course` from `source`: a new
    /// one, or the one already there (which gains the source).
    fn add(&mut self, from: usize, to: usize, course: Course, s: Source) -> usize {
        debug_assert_eq!(
            course,
            Course::Straight,
            "a curved segment needs its midpoint"
        );
        let k = key(from, to);
        let i = *self.straight.entry(k).or_insert_with(|| {
            self.list.push(Seg {
                a: k.0,
                b: k.1,
                course,
                sources: Vec::new(),
                mid: None,
            });
            self.list.len() - 1
        });
        self.source(i, s);
        i
    }

    /// The curved segment from `from` to `to` along `course` (given from
    /// `from`), with `mid` the point halfway along it, from source `s`: a
    /// new one, or one already between the two clusters whose midpoint is
    /// within the tolerance `eps` of `mid` (two pieces of one course, an
    /// operand's edge and the curve where faces meet along it), which
    /// gains the source.
    fn add_curved(
        &mut self,
        from: usize,
        to: usize,
        course: Course,
        mid: DVec3,
        s: Source,
        eps: f64,
    ) -> usize {
        let k = key(from, to);
        let found = self.curved.get(&k).and_then(|l| {
            l.iter()
                .copied()
                .find(|&i| self.list[i].mid.is_some_and(|m| (m - mid).length() <= eps))
        });
        let i = found.unwrap_or_else(|| {
            self.list.push(Seg {
                a: k.0,
                b: k.1,
                course: if from <= to {
                    course
                } else {
                    course.reversed()
                },
                sources: Vec::new(),
                mid: Some(mid),
            });
            let i = self.list.len() - 1;
            self.curved.entry(k).or_default().push(i);
            i
        });
        self.source(i, s);
        i
    }

    fn source(&mut self, i: usize, s: Source) {
        let list = &mut self.list[i].sources;
        if !list.contains(&s) {
            list.push(s);
        }
    }

    /// The cluster segment `s` starts at when run forward (from `a` to `b`)
    /// or not.
    fn start(&self, s: usize, forward: bool) -> usize {
        let seg = &self.list[s];
        if forward {
            seg.a
        } else {
            seg.b
        }
    }

    /// The cluster segment `s` ends at when run forward or not.
    fn end(&self, s: usize, forward: bool) -> usize {
        self.start(s, !forward)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::math::Frame;
    use crate::{BrepKernel, Kernel};

    fn cube(at: DVec3, size: f64) -> Brep {
        let frame = Frame {
            origin: at,
            ..Frame::WORLD
        };
        let body = BrepKernel::default()
            .make_box(&frame, DVec3::splat(size))
            .expect("a box");
        (*body.brep).clone()
    }

    /// A thin shell that a coarse mesh flattens to nothing takes its
    /// orientation from its mass properties: the ring a ball keeps outside
    /// a cone whose side crosses its sphere twice close together. (Found
    /// by the fuzz of balls and cones at 50 times its size: validation
    /// refused such rings as inside out, by the sign of rounding.)
    #[test]
    fn a_thin_ring_is_measured_exactly() {
        let tol = Tolerance::default();
        let upright = |z: f64| Frame {
            origin: DVec3::new(0.0, 0.0, z),
            x: DVec3::X,
            y: DVec3::Z,
            z: -DVec3::Y,
        };
        let (apex, slope, base) = (2.0, 0.57, -0.2);
        let ball = testing::ball(upright(0.0), 1.0);
        let cone = testing::cone(upright(base), slope * (apex - base), apex - base);
        let (r, _) = boolean(&ball, &cone, BooleanOp::Subtract, &tol).expect("a difference");
        assert_eq!(crate::check::validate(&r, &tol), Vec::<String>::new());
        // The ring's volume: the integral of π((1 − z²) − m²(a − z)²)
        // between its roots, the circles where the cone crosses the sphere.
        let (qa, qb, qc) = (
            -(1.0 + slope * slope),
            2.0 * slope * slope * apex,
            1.0 - slope * slope * apex * apex,
        );
        let root = (qb * qb - 4.0 * qa * qc).sqrt();
        let (z1, z2) = ((-qb + root) / (2.0 * qa), (-qb - root) / (2.0 * qa));
        let integral = |z: f64| PI * (qa * z.powi(3) / 3.0 + 0.5 * qb * z * z + qc * z);
        let ring = (integral(z2) - integral(z1)).abs();
        let coarse = crate::check::shell_volumes(&r, &tol).expect("coarse volumes");
        let signed = crate::check::shell_signed_volumes(&r, &tol).expect("volumes");
        let thin: Vec<usize> = (0..coarse.len())
            .filter(|&i| coarse[i].0.abs() <= coarse[i].1)
            .collect();
        assert_eq!(thin.len(), 1, "{coarse:?}");
        assert!(
            (signed[thin[0]] - ring).abs() <= 1e-6 * ring,
            "{} vs {ring}",
            signed[thin[0]]
        );
    }

    #[test]
    fn work_limits_refuse_rather_than_run_on() {
        let tol = Tolerance::default();
        // Flush top and bottom faces, so edges of the two cross in those planes.
        let (a, b) = (cube(DVec3::ZERO, 1.0), cube(DVec3::new(0.5, 0.5, 0.0), 1.0));
        assert!(boolean_within(&a, &b, BooleanOp::Union, &tol, &LIMITS).is_ok());
        let tight = [
            Limits {
                face_pairs: 3,
                ..LIMITS
            },
            Limits {
                edge_pairs: 0,
                ..LIMITS
            },
            Limits {
                points: 20,
                ..LIMITS
            },
            Limits {
                segments: 10,
                ..LIMITS
            },
        ];
        for limits in tight {
            match boolean_within(&a, &b, BooleanOp::Union, &tol, &limits) {
                Err(KernelError::TooLarge(m)) => assert!(m.contains("refused"), "{m}"),
                other => panic!("{limits:?}: {:?}", other.map(|_| ())),
            }
        }
    }
}
