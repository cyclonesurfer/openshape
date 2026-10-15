//! Imprint (design note 0012 §4.2 and §10.3): how two faces of the
//! operands lie, and the stretches of the curves where crossing faces meet
//! that lie in both (closed) faces.
//!
//! Faces on two planes are step 1's: their common line, and the stretches
//! of it in each face from the signs of the face's boundary against the
//! other plane. Faces with a curved surface meet along the curves (tracks)
//! of [`meet`], found once per pair of surfaces. A track enters and leaves
//! a face only where the face's boundary crosses the other surface, so the
//! cuts along it are those crossings and the boundary's stations on that
//! surface (decided once, in `points`), taken by step 1's two half-open
//! rules. On an open track (a line or a conic, whose parameter runs along
//! it one way) the cuts of one rule pair up in order; on a closed track (a
//! circle or an ellipse) each cut says whether the track enters or
//! leaves the face there, from the side of the other surface the boundary
//! goes to and the way the track runs across the face. A closed track that
//! no cut reaches lies wholly in the face or wholly out of it, decided at
//! its points.

use std::cmp::Ordering;
use std::collections::BTreeMap;
use std::f64::consts::{FRAC_PI_2, TAU};

use glam::DVec3;
use rustc_hash::FxHashMap;

use super::chart::Outline;
use super::meet::{self, Meet, Track};
use super::points::{self, Points};
use super::quadric::{self, Quadric};
use super::shape::{Face, Shape};
use super::{name, other, show, slot, Pair};
use crate::origin::Operand;
use crate::KernelError;

/// How two faces of different operands lie.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum Relation {
    /// One face lies wholly on one side of the other's surface.
    Apart,
    /// On one surface (either orientation).
    Coplanar,
    /// The faces' surfaces cross, and each face reaches the other's
    /// surface.
    Crossing,
}

/// The outlines of an operand's faces that are not flat; step 1's exact
/// tests serve the flat ones.
pub(super) fn outlines<'a>(shape: &'a Shape) -> Result<Vec<Option<Outline<'a>>>, KernelError> {
    shape
        .faces
        .iter()
        .map(|f| {
            if f.flat {
                Ok(None)
            } else {
                Outline::new(shape, f).map(Some)
            }
        })
        .collect()
}

/// Whether `p` (on the surface of face `f` of `shape`, whose outline is
/// `outline` if it is not flat) lies in the face: `Some(true)` inside,
/// `Some(false)` outside, `None` within the tolerance of its boundary.
pub(super) fn locate(
    shape: &Shape,
    outline: Option<&Outline>,
    f: &Face,
    p: DVec3,
    eps: f64,
) -> Option<bool> {
    match outline {
        Some(o) => o.locate(o.chart.map(p), eps, eps / 16.0),
        None => {
            for l in &f.loops {
                for &(e, _) in l {
                    let (a, b) = shape.edge_points(e);
                    if points::segment_distance(p, a, b) <= eps {
                        return None;
                    }
                }
            }
            Some(in_face(shape, f, p))
        }
    }
}

/// True if `p` (on the plane of the flat face `face`) lies inside it, by
/// the parity of the face's edges crossed by a ray from it; for a point
/// off its edges.
pub(super) fn in_face(shape: &Shape, face: &Face, p: DVec3) -> bool {
    let q = face.to2(p);
    let mut inside = false;
    for vs in &face.verts {
        for k in 0..vs.len() {
            let a = face.to2(shape.points[vs[k]]);
            let b = face.to2(shape.points[vs[(k + 1) % vs.len()]]);
            if (a.y > q.y) != (b.y > q.y) {
                let x = a.x + (q.y - a.y) / (b.y - a.y) * (b.x - a.x);
                if x > q.x {
                    inside = !inside;
                }
            }
        }
    }
    inside
}

/// How a surface of the first operand and a surface of the second meet.
#[derive(Clone, Debug)]
enum Contact {
    Apart,
    Same,
    Touch,
    Unsupported(&'static str),
    /// Along these tracks (indices into [`Contacts::tracks`]).
    Tracks(Vec<usize>),
}

/// How the surfaces of the two operands meet, decided once per pair of
/// surfaces (as step 1 decides each vertex's side of each plane once), over
/// the part of space where faces on both can be.
pub(super) struct Contacts {
    /// The curves where surfaces meet.
    pub tracks: Vec<Track>,
    meets: FxHashMap<(usize, usize), Contact>,
    /// Per operand, per surface, the box of its faces.
    boxes: [Vec<(DVec3, DVec3)>; 2],
    /// Per operand, per surface, for a cone the least height of its faces
    /// above its apex along its axis ([`meet::Region::lift`]); infinite for
    /// other surfaces.
    lifts: [Vec<f64>; 2],
}

impl Contacts {
    pub fn new(pair: &Pair) -> Contacts {
        let boxes = [&pair.first, &pair.second].map(|s| {
            let mut b = vec![(DVec3::splat(f64::MAX), DVec3::splat(f64::MIN)); s.surfaces.len()];
            for f in &s.faces {
                let x = &mut b[f.surface];
                *x = (x.0.min(f.lo), x.1.max(f.hi));
            }
            b
        });
        let lifts = [&pair.first, &pair.second].map(lifts);
        Contacts {
            tracks: Vec::new(),
            meets: FxHashMap::default(),
            boxes,
            lifts,
        }
    }

    /// How surface `sa` of the first operand and `sb` of the second meet,
    /// within the common part of the boxes of their faces (grown by the
    /// tolerance): any point where faces on the two meet is there.
    fn get(&mut self, pair: &Pair, sa: usize, sb: usize, eps: f64) -> &Contact {
        if !self.meets.contains_key(&(sa, sb)) {
            let (a, b) = (self.boxes[0][sa], self.boxes[1][sb]);
            let lo = a.0.max(b.0) - DVec3::splat(eps);
            let hi = a.1.min(b.1) + DVec3::splat(eps);
            let c = if lo.cmpgt(hi).any() {
                Contact::Apart
            } else {
                let region = meet::Region {
                    lift: self.lifts[0][sa].min(self.lifts[1][sb]),
                    ..meet::Region::of_box(lo, hi)
                };
                let (qa, qb) = (&pair.first.surfaces[sa], &pair.second.surfaces[sb]);
                match meet::meet(qa, qb, &region, eps) {
                    Meet::Apart => Contact::Apart,
                    Meet::Same => Contact::Same,
                    Meet::Touch => Contact::Touch,
                    Meet::Unsupported(m) => Contact::Unsupported(m),
                    Meet::Curves(list) => {
                        let first = self.tracks.len();
                        self.tracks.extend(list);
                        Contact::Tracks((first..self.tracks.len()).collect())
                    }
                }
            };
            self.meets.insert((sa, sb), c);
        }
        &self.meets[&(sa, sb)]
    }

    /// The tracks along which surface `sa` of the first operand and `sb`
    /// of the second cross, once decided.
    fn tracks_of(&self, sa: usize, sb: usize) -> &[usize] {
        match self.meets.get(&(sa, sb)) {
            Some(Contact::Tracks(ids)) => ids,
            _ => &[],
        }
    }
}

/// Per surface of `shape`, for a cone the least height of its faces above
/// its apex along its axis, from the control points of their edges' Bézier
/// pieces: a face on a cone lies in the convex hull of its boundary (each
/// of its rulings ends on it), and that in the hull of those points.
/// Infinite for other surfaces; minus infinity for a cone with an edge
/// that has no such pieces.
fn lifts(shape: &Shape) -> Vec<f64> {
    let mut out = vec![f64::INFINITY; shape.surfaces.len()];
    for f in &shape.faces {
        let Quadric::Cone {
            point,
            axis,
            radius,
            slope,
        } = shape.surfaces[f.surface]
        else {
            continue;
        };
        let apex = point - axis * (radius / slope);
        for &(e, _) in f.loops.iter().flatten() {
            let (curve, t0, t1) = shape.curve(e);
            let low = quadric::hull_points(curve, t0, t1).map_or(f64::NEG_INFINITY, |pieces| {
                pieces
                    .iter()
                    .flatten()
                    .map(|&p| (p - apex).dot(axis))
                    .fold(f64::INFINITY, f64::min)
            });
            out[f.surface] = out[f.surface].min(low);
        }
    }
    out
}

fn is_plane(q: &Quadric) -> bool {
    matches!(q, Quadric::Plane { .. })
}

fn is_sphere(q: &Quadric) -> bool {
    matches!(q, Quadric::Sphere { .. })
}

/// Which sides of `surface` (of the other operand) the boundary of face
/// `f` of `op` reaches: above it, below it, on it (within the tolerance).
/// The coedges near the surface by their stations (a straight edge against
/// a plane by its start vertex, as in step 1: its end is the next
/// coedge's start), the others by their boxes.
fn sides(
    pair: &Pair,
    pts: &mut Points,
    op: Operand,
    f: usize,
    surface: usize,
) -> Result<(bool, bool, bool), KernelError> {
    let shape = pair.get(op);
    let face = &shape.faces[f];
    let q = &pair.get(other(op)).surfaces[surface];
    let mut near = Vec::new();
    let (mut pos, mut neg) = face.near_surface(q, pts.eps, &mut near);
    let mut zero = false;
    let mut note = |s: i8| match s {
        1 => pos = true,
        -1 => neg = true,
        _ => zero = true,
    };
    for k in near {
        let (li, i) = face.coedge_at(k);
        let e = face.loops[li][i].0;
        if shape.straight[e] && is_plane(q) {
            note(pts.sign(op, face.verts[li][i], surface));
        } else {
            for st in pts.stations(pair, op, e, surface)? {
                note(st.sign);
            }
        }
    }
    Ok((pos, neg, zero))
}

/// How faces `fa` of the first operand and `fb` of the second lie.
pub(super) fn relation(
    pair: &Pair,
    pts: &mut Points,
    contacts: &mut Contacts,
    fa: usize,
    fb: usize,
) -> Result<Relation, KernelError> {
    let a = &pair.first.faces[fa];
    let b = &pair.second.faces[fb];
    let (qa, qb) = (
        &pair.first.surfaces[a.surface],
        &pair.second.surfaces[b.surface],
    );
    let (ap, an, az) = sides(pair, pts, Operand::First, fa, b.surface)?;
    let (bp, bn, bz) = sides(pair, pts, Operand::Second, fb, a.surface)?;
    let corner = || pair.first.points[a.verts[0][0]];
    if is_plane(qa) && is_plane(qb) {
        let a_on = az && !ap && !an;
        let b_on = bz && !bp && !bn;
        if a_on && b_on {
            return Ok(Relation::Coplanar);
        }
        if a_on || b_on {
            return Err(KernelError::InvalidInput(format!(
                "two faces of the bodies nearly coincide but are not on one plane, \
                 within the tolerance, near {}",
                show(corner())
            )));
        }
        let apart = |p: bool, n: bool, z: bool| !z && (p != n);
        if apart(ap, an, az) || apart(bp, bn, bz) {
            return Ok(Relation::Apart);
        }
        return Ok(Relation::Crossing);
    }
    // A boundary strictly on one side of the other surface keeps the face
    // there if the face cannot reach across without its boundary doing so.
    // A face on a plane, a cylinder or a cone lies in the convex hull of its
    // boundary (through each of its points runs a line of the surface that
    // leaves the face through its boundary both ways), so a boundary on one
    // side of a plane, or inside a cylinder, cone or sphere (whose insides
    // are convex), keeps the face there. A face on a sphere bulges out of
    // its boundary's hull, and a face outside a curved surface can still
    // hold a closed curve of their meeting (a plate around a hole), so
    // those are decided by how the surfaces meet.
    let hull = |own: &Quadric, theirs: &Quadric, p: bool, n: bool, z: bool| {
        !z && p != n && !is_sphere(own) && (is_plane(theirs) || n)
    };
    if hull(qa, qb, ap, an, az) || hull(qb, qa, bp, bn, bz) {
        return Ok(Relation::Apart);
    }
    // A face on a cylinder or a cone (at most half a turn, ADR 0008) holds
    // no closed curve of a meeting found in closed form with a curved
    // surface: those are lines along its axis, which run out through its
    // boundary, and circles about its axis, which go around it. Curves of
    // higher degree (two cylinders at an angle) can close inside the face,
    // so an unsupported meeting is refused even then.
    let wraps = |own: &Quadric, p: bool, n: bool, z: bool| {
        !z && p != n && matches!(own, Quadric::Cylinder { .. } | Quadric::Cone { .. })
    };
    let wrapped = wraps(qa, ap, an, az) || wraps(qb, bp, bn, bz);
    match contacts.get(pair, a.surface, b.surface, pts.eps) {
        Contact::Apart => Ok(Relation::Apart),
        Contact::Same => Ok(Relation::Coplanar),
        Contact::Touch | Contact::Tracks(_) if wrapped => Ok(Relation::Apart),
        Contact::Touch => Err(KernelError::Unsupported(format!(
            "booleans of bodies whose faces touch without crossing (a tangent contact; \
             design note 0011 §4.4) near {}",
            show(corner())
        ))),
        Contact::Unsupported(m) => Err(KernelError::Unsupported(format!(
            "booleans of bodies whose faces meet where {m}, near {}",
            show(corner())
        ))),
        Contact::Tracks(_) => Ok(Relation::Crossing),
    }
}

/// A place where a face's boundary crosses a surface of the other operand
/// under one of the half-open rules: its cluster, and whether the boundary
/// goes there to the surface's positive side (where its distance grows).
#[derive(Clone, Copy, Debug)]
struct Cut {
    cluster: usize,
    up: bool,
}

fn missing(op: Operand, e: usize) -> KernelError {
    KernelError::InvalidResult(format!(
        "boolean: no crossing recorded for edge {e} of the {} body",
        name(op)
    ))
}

/// The cuts of face `f` of `op` against `surface` of the other operand,
/// counting points within the tolerance of it as on its positive side
/// (`zero_positive`) or its negative side, over the coedges `near` it
/// ([`Face::near_surface`]; the others stay on one side), in loop order.
fn cuts(
    pair: &Pair,
    pts: &Points,
    op: Operand,
    f: usize,
    surface: usize,
    near: &[usize],
    zero_positive: bool,
) -> Result<Vec<Cut>, KernelError> {
    let shape = pair.get(op);
    let face = &shape.faces[f];
    let plane = is_plane(&pair.get(other(op)).surfaces[surface]);
    let side = |s: i8| if zero_positive { s >= 0 } else { s > 0 };
    let mut out = Vec::new();
    for &c in near {
        let (li, k) = face.coedge_at(c);
        let (e, forward) = face.loops[li][k];
        if shape.straight[e] && plane {
            let vs = &face.verts[li];
            let (va, vb) = (vs[k], vs[(k + 1) % vs.len()]);
            let (sa, sb) = (pts.sign(op, va, surface), pts.sign(op, vb, surface));
            if side(sa) == side(sb) {
                continue;
            }
            let cluster = if sa == 0 {
                pts.vertex_cluster(op, va)
            } else if sb == 0 {
                pts.vertex_cluster(op, vb)
            } else {
                pts.crossing_cluster(op, e, surface, 0)
                    .ok_or_else(|| missing(op, e))?
            };
            out.push(Cut {
                cluster,
                up: side(sb),
            });
            continue;
        }
        let st = pts
            .cached_stations(op, e, surface)
            .ok_or_else(|| missing(op, e))?;
        let n = st.len();
        for w in 0..n - 1 {
            // Stations `i` then `j` along the coedge.
            let (i, j) = if forward {
                (w, w + 1)
            } else {
                (n - 1 - w, n - 2 - w)
            };
            let (sa, sb) = (st[i].sign, st[j].sign);
            if side(sa) == side(sb) {
                continue;
            }
            let cluster = if sa == 0 {
                pts.station_cluster(pair, op, e, surface, i)
            } else if sb == 0 {
                pts.station_cluster(pair, op, e, surface, j)
            } else {
                pts.crossing_cluster(op, e, surface, i.min(j))
            };
            out.push(Cut {
                cluster: cluster.ok_or_else(|| missing(op, e))?,
                up: side(sb),
            });
        }
    }
    Ok(out)
}

/// Where a loop of a face meets a surface: a vertex of the face, or a
/// station (by its index along the edge) of an edge.
#[derive(Clone, Copy)]
enum Mark {
    Vertex(usize),
    Station(usize, usize),
}

/// The clusters where the boundary of face `f` of `op` touches `surface`
/// of the other operand from one side: a station or a vertex on it, or a
/// run of them, between two strictly on one side along a loop (tangent
/// contacts, at a vertex too: a rim touching a plane where a seam meets
/// it, found by review). A run holding a whole edge is that edge lying on
/// the surface (a rim on a plate), not a touch; nor is a run through a
/// vertex that another edge leaves to the other side, or whose faces all
/// lie on one surface (smooth there), where the operand crosses the
/// surface and the curve where they meet passes through the vertex (a
/// box's corner on a cylinder or a plane, a cylinder crossing a rod along
/// its seam, a circle through a rod's cap's centre along the line that
/// halves the cap; found refused as touches by diffing the fuzz's refusals
/// against those of the code before this rule).
fn touches(
    pair: &Pair,
    pts: &Points,
    op: Operand,
    f: usize,
    surface: usize,
    near: &[usize],
) -> Result<Vec<usize>, KernelError> {
    let shape = pair.get(op);
    let face = &shape.faces[f];
    let plane = is_plane(&pair.get(other(op)).surfaces[surface]);
    let mut out = Vec::new();
    // Along a run of coedges: each station's side, where it is, and
    // whether its edge lies wholly on the surface.
    let mut seq: Vec<(i8, Mark, bool)> = Vec::new();
    // Per loop, its coedges near the surface (by position in the loop).
    let mut k = 0;
    while k < near.len() {
        let li = face.coedge_at(near[k]).0;
        let mut idx = Vec::new();
        while k < near.len() && face.coedge_at(near[k]).0 == li {
            idx.push(face.coedge_at(near[k]).1);
            k += 1;
        }
        let len = face.loops[li].len();
        let whole = idx.len() == len;
        // Runs of consecutive coedges, from one after a coedge not near
        // (strictly on one side, so a run starts and ends off the surface
        // unless it is the whole loop).
        if !whole {
            let at = (0..idx.len())
                .find(|&i| idx[(i + idx.len() - 1) % idx.len()] != (idx[i] + len - 1) % len)
                .unwrap_or(0);
            idx.rotate_left(at);
        }
        let mut i = 0;
        while i < idx.len() {
            let mut j = i;
            while j + 1 < idx.len() && idx[j + 1] == (idx[j] + 1) % len {
                j += 1;
            }
            seq.clear();
            for (n, &c) in idx[i..=j].iter().enumerate() {
                let (e, forward) = face.loops[li][c];
                let (va, vb) = shape.edges[e];
                let mut marks: Vec<(i8, Mark)> = match pts.cached_stations(op, e, surface) {
                    Some(st) if !(shape.straight[e] && plane) => {
                        let last = st.len() - 1;
                        st.iter()
                            .enumerate()
                            .map(|(s, x)| {
                                let m = match s {
                                    0 => Mark::Vertex(va),
                                    s if s == last => Mark::Vertex(vb),
                                    s => Mark::Station(e, s),
                                };
                                (x.sign, m)
                            })
                            .collect()
                    }
                    _ => vec![
                        (pts.sign(op, va, surface), Mark::Vertex(va)),
                        (pts.sign(op, vb, surface), Mark::Vertex(vb)),
                    ],
                };
                if !forward {
                    marks.reverse();
                }
                let on = marks.iter().all(|x| x.0 == 0);
                if whole || n < j - i {
                    // The end is the next coedge's start.
                    marks.pop();
                }
                seq.extend(marks.into_iter().map(|(sign, m)| (sign, m, on)));
            }
            if whole {
                // Around the loop from a point off the surface back to it.
                let Some(at) = seq.iter().position(|x| x.0 != 0) else {
                    break;
                };
                seq.rotate_left(at);
                seq.push(seq[0]);
            }
            let mut a = 1;
            while a + 1 < seq.len() {
                if seq[a].0 != 0 || seq[a - 1].0 == 0 {
                    a += 1;
                    continue;
                }
                let mut b = a;
                while b + 1 < seq.len() && seq[b + 1].0 == 0 {
                    b += 1;
                }
                let lying = seq[a..=b].iter().any(|x| x.2);
                let mut crossed = false;
                for x in &seq[a..=b] {
                    if let Mark::Vertex(v) = x.1 {
                        crossed |= smooth_at(shape, v)
                            || leaves_to(pair, pts, op, v, surface, -seq[a - 1].0)?;
                    }
                }
                if b + 1 < seq.len() && seq[a - 1].0 == seq[b + 1].0 && !lying && !crossed {
                    let cluster = match seq[a].1 {
                        Mark::Vertex(v) => Some(pts.vertex_cluster(op, v)),
                        Mark::Station(e, s) => pts.station_cluster(pair, op, e, surface, s),
                    };
                    let e = match seq[a].1 {
                        Mark::Station(e, _) => e,
                        Mark::Vertex(_) => face.loops[li][idx[i]].0,
                    };
                    out.push(cluster.ok_or_else(|| missing(op, e))?);
                }
                a = b + 1;
            }
            i = j + 1;
        }
    }
    Ok(out)
}

/// Whether vertex `v` of `shape` lies inside a smooth piece of its
/// surface: all its faces on one surface, not a cone (whose apex is not
/// smooth). Such a piece crosses any other surface through `v` (one
/// touching it there is refused by the surfaces' contact instead).
fn smooth_at(shape: &Shape, v: usize) -> bool {
    shape
        .one_surface(v)
        .is_some_and(|s| !matches!(shape.surfaces[s], Quadric::Cone { .. }))
}

/// Whether an edge of `op` at its vertex `v`, which lies on `surface` of
/// the other operand, leaves it strictly on side `side` (by its first
/// station past `v`).
fn leaves_to(
    pair: &Pair,
    pts: &Points,
    op: Operand,
    v: usize,
    surface: usize,
    side: i8,
) -> Result<bool, KernelError> {
    let shape = pair.get(op);
    for &e in shape.edges_at(v) {
        if pts.leaving(pair, op, e, surface, shape.edges[e].0 == v)? == side {
            return Ok(true);
        }
    }
    Ok(false)
}

/// Refuses an edge of either face touching the other face's surface at a
/// point of the other face (closed): a tangent contact.
fn refuse_touches(
    pair: &Pair,
    pts: &Points,
    outlines: &[Vec<Option<Outline>>; 2],
    fa: usize,
    fb: usize,
) -> Result<(), KernelError> {
    let eps = pts.eps;
    for (op, f, g) in [(Operand::First, fa, fb), (Operand::Second, fb, fa)] {
        let (mine, theirs) = (pair.get(op), pair.get(other(op)));
        let gf = &theirs.faces[g];
        let q = &theirs.surfaces[gf.surface];
        let mut near = Vec::new();
        mine.faces[f].near_surface(q, eps, &mut near);
        for c in touches(pair, pts, op, f, gf.surface, &near)? {
            let p = pts.rep[c];
            let outline = outlines[slot(other(op))][g].as_ref();
            if locate(theirs, outline, gf, p, eps) != Some(false) {
                return Err(KernelError::Unsupported(format!(
                    "booleans where an edge of one body touches a face of the other without \
                     crossing it (a tangent contact; design note 0011 §4.4) near {}",
                    show(p)
                )));
            }
        }
    }
    Ok(())
}

/// A closed interval of the common line of two face planes, between two
/// clusters, with their positions along the line.
#[derive(Clone, Copy, Debug)]
struct Interval {
    lo: (f64, usize),
    hi: (f64, usize),
}

/// The stretches of the common line of two crossing faces on planes that
/// lie in both (closed) faces, as pairs of clusters (step 1's meeting).
pub(super) fn plane_meeting(
    pair: &Pair,
    pts: &Points,
    outlines: &[Vec<Option<Outline>>; 2],
    fa: usize,
    fb: usize,
) -> Result<Vec<(usize, usize)>, KernelError> {
    let a = &pair.first.faces[fa];
    let b = &pair.second.faces[fb];
    let dir = a.normal().cross(b.normal());
    if dir.length() == 0.0 {
        return Ok(Vec::new());
    }
    let dir = dir.normalize();
    let eps = pts.eps;
    let sa = face_line(pair, pts, Operand::First, fa, b.surface, dir)?;
    let sb = face_line(pair, pts, Operand::Second, fb, a.surface, dir)?;
    let mut out = Vec::new();
    for x in &sa {
        for y in &sb {
            let lo = later(x.lo, y.lo, eps, pts)?;
            let hi = earlier(x.hi, y.hi, eps, pts)?;
            if lo.1 != hi.1 && lo.0 < hi.0 {
                out.push((lo.1, hi.1));
            }
        }
    }
    if !(a.flat && b.flat) {
        refuse_touches(pair, pts, outlines, fa, fb)?;
    }
    Ok(out)
}

fn compare(
    x: (f64, usize),
    y: (f64, usize),
    eps: f64,
    pts: &Points,
) -> Result<Ordering, KernelError> {
    if x.1 == y.1 {
        return Ok(Ordering::Equal);
    }
    if (x.0 - y.0).abs() <= eps {
        return Err(close(pts.rep[x.1]));
    }
    Ok(x.0.total_cmp(&y.0))
}

fn close(p: DVec3) -> KernelError {
    KernelError::InvalidInput(format!(
        "the bodies meet in points closer than the tolerance near {}",
        show(p)
    ))
}

fn later(
    x: (f64, usize),
    y: (f64, usize),
    eps: f64,
    pts: &Points,
) -> Result<(f64, usize), KernelError> {
    Ok(match compare(x, y, eps, pts)? {
        Ordering::Less => y,
        _ => x,
    })
}

fn earlier(
    x: (f64, usize),
    y: (f64, usize),
    eps: f64,
    pts: &Points,
) -> Result<(f64, usize), KernelError> {
    Ok(match compare(x, y, eps, pts)? {
        Ordering::Greater => y,
        _ => x,
    })
}

/// The closed face `f` of `op` cut by the plane `surface` of the other
/// operand: the intervals of the planes' common line (direction `dir`)
/// inside the face.
///
/// Two half-open rules, counting a point on the plane as on its positive
/// side and then as on its negative side, give the interior together with
/// the boundary edges lying in the plane from either side; their union is
/// the closed face's intersection with the line.
fn face_line(
    pair: &Pair,
    pts: &Points,
    op: Operand,
    f: usize,
    surface: usize,
    dir: DVec3,
) -> Result<Vec<Interval>, KernelError> {
    let face = &pair.get(op).faces[f];
    let eps = pts.eps;
    let mut spans: Vec<Interval> = Vec::new();
    // Only coedges near the plane can cross it under either rule.
    let mut near = Vec::new();
    face.near_surface(&pair.get(other(op)).surfaces[surface], eps, &mut near);
    for zero_positive in [true, false] {
        let mut cuts: Vec<(f64, usize)> = cuts(pair, pts, op, f, surface, &near, zero_positive)?
            .into_iter()
            .map(|c| (pts.rep[c.cluster].dot(dir), c.cluster))
            .collect();
        cuts.sort_by(|x, y| x.0.total_cmp(&y.0).then(x.1.cmp(&y.1)));
        for w in cuts.windows(2) {
            if w[0].1 != w[1].1 && (w[1].0 - w[0].0) <= eps {
                return Err(close(pts.rep[w[0].1]));
            }
        }
        if cuts.len() % 2 != 0 {
            return Err(KernelError::InvalidResult(format!(
                "boolean: face {f} of the {} body crosses a plane an odd number of times",
                name(op)
            )));
        }
        for c in cuts.chunks(2) {
            spans.push(Interval { lo: c[0], hi: c[1] });
        }
    }
    // Union of the two rules' intervals.
    spans.sort_by(|x, y| x.lo.0.total_cmp(&y.lo.0));
    let mut out: Vec<Interval> = Vec::new();
    for s in spans {
        if let Some(last) = out.last_mut() {
            if s.lo.1 == last.hi.1 || s.lo.0 <= last.hi.0 {
                if s.hi.0 > last.hi.0 {
                    last.hi = s.hi;
                }
                continue;
            }
        }
        out.push(s);
    }
    Ok(out)
}

/// A stretch of a track inside two crossing faces: from cluster `lo.0` at
/// parameter `lo.1` to cluster `hi.0` at `hi.1`, the parameter increasing
/// (by less than a turn on a closed track).
#[derive(Clone, Copy, Debug)]
pub(super) struct Stretch {
    pub track: usize,
    pub lo: (usize, f64),
    pub hi: (usize, f64),
}

/// A point on a track: its cluster and its parameter along the track.
type End = (usize, f64);

/// The cluster of a span's end where an open track ends inside a face, at
/// a parameter of −∞ or +∞ (its start or its end): no point of the
/// result.
const OPEN: usize = usize::MAX;

/// The part of a track in a face under one rule.
#[derive(Clone, Copy, Debug)]
enum Span {
    /// The whole closed track.
    Full,
    /// From cluster `lo.0` at parameter `lo.1` to `hi.0` at `hi.1`, along
    /// the track's direction (on an open track `lo.1 < hi.1`, either end
    /// [`OPEN`] where the track ends in the face; on a closed track the
    /// parameters as the track gives them, the span running from `lo` the
    /// way the parameter grows).
    Part { lo: End, hi: End },
}

/// How far, in tolerances, a cut may lie from the track it is on: a point
/// within the tolerance of both surfaces lies within the tolerance over
/// the sine of their angle of the curve where they meet, and the result's
/// edges along tracks end at these points (the validator allows a hundred
/// tolerances).
const ON_TRACK: f64 = 16.0;

/// The track (by position in `tracks`) a cut at `p` is on, and its
/// parameter there: the nearest, within [`ON_TRACK`] tolerances, and no
/// other that near. `None` for a cut on a ray's line well before where it
/// starts beyond the region (§14), which no face on the other surface
/// reaches; a cut within half [`meet::CLIP`] tolerances of that start is
/// refused, as the face's side of the start is found from the face (where
/// the cut could be on either side of it).
fn assign(tracks: &[&Track], p: DVec3, eps: f64) -> Result<Option<(usize, f64)>, KernelError> {
    let mut before = false;
    for t in tracks {
        if let Track::Ray {
            origin,
            dir,
            start,
            from_apex: false,
        } = **t
        {
            let w = p - origin;
            let s = w.dot(dir);
            if s <= start + 0.5 * meet::CLIP * eps && (w - dir * s).length() <= ON_TRACK * eps {
                if s >= start - 0.5 * meet::CLIP * eps {
                    return Err(through_apex(p));
                }
                before = true;
            }
        }
    }
    let mut best: Option<(usize, f64, f64)> = None;
    let mut second = f64::INFINITY;
    for (k, t) in tracks.iter().enumerate() {
        if t.ends().0 == f64::INFINITY {
            // A ray that stays away from the region.
            continue;
        }
        let (u, d) = t.nearest(p);
        match best {
            Some(b) if d >= b.2 => second = second.min(d),
            _ => {
                if let Some(b) = best {
                    second = second.min(b.2);
                }
                best = Some((k, u, d));
            }
        }
    }
    let (k, u, d) = best.ok_or_else(|| {
        KernelError::InvalidResult(format!(
            "boolean: faces cross along no curve near {}",
            show(p)
        ))
    })?;
    if before {
        // Unless it is on a track where it matters, too.
        if d <= ON_TRACK * eps {
            return Err(nearly_touch(p));
        }
        return Ok(None);
    }
    if d > ON_TRACK * eps {
        return Err(KernelError::InvalidInput(format!(
            "faces of the bodies meet at too shallow an angle near {}",
            show(p)
        )));
    }
    if second <= ON_TRACK * eps {
        // Rays meet at their apex: a pointed cone's, on the plane but
        // outside the faces of it that meet the cone (§14).
        if let Track::Ray { start, .. } = *tracks[k] {
            if (tracks[k].point(start) - p).length() <= ON_TRACK * eps {
                return Err(through_apex(p));
            }
        }
        return Err(nearly_touch(p));
    }
    Ok(Some((k, u)))
}

fn nearly_touch(p: DVec3) -> KernelError {
    KernelError::Unsupported(format!(
        "booleans of bodies whose faces nearly touch (design note 0011 §4.4) near {}",
        show(p)
    ))
}

/// The refusal of a plane through a cone's apex where a face's boundary
/// reaches it, or crosses a ray from it or runs along one at the ray's
/// start (§14).
fn through_apex(p: DVec3) -> KernelError {
    KernelError::Unsupported(format!(
        "booleans of bodies whose faces meet where a plane through a cone's apex meets it in \
         two rulings (design note 0012 §14), near {}",
        show(p)
    ))
}

/// Per rule (points on the surface counted positive, then negative), per
/// track of `ids`, the spans of the tracks in face `f` of `op`, which
/// meets the surface of face `g` of the other operand along them.
#[allow(clippy::too_many_arguments)]
fn face_spans(
    pair: &Pair,
    pts: &Points,
    outlines: &[Vec<Option<Outline>>; 2],
    contacts: &Contacts,
    op: Operand,
    f: usize,
    g: usize,
    ids: &[usize],
) -> Result<[Vec<Vec<Span>>; 2], KernelError> {
    let shape = pair.get(op);
    let face = &shape.faces[f];
    let gf = &pair.get(other(op)).faces[g];
    let surface = gf.surface;
    let theirs = &pair.get(other(op)).surfaces[surface];
    let eps = pts.eps;
    let region = (face.lo.max(gf.lo) - eps, face.hi.min(gf.hi) + eps);
    let mut near = Vec::new();
    face.near_surface(theirs, eps, &mut near);
    let tracks: Vec<&Track> = ids.iter().map(|&i| &contacts.tracks[i]).collect();
    let mut out = [vec![Vec::new(); ids.len()], vec![Vec::new(); ids.len()]];
    for (r, zero_positive) in [true, false].into_iter().enumerate() {
        let mut on: Vec<Vec<(f64, usize, bool)>> = vec![Vec::new(); ids.len()];
        for c in cuts(pair, pts, op, f, surface, &near, zero_positive)? {
            if let Some((k, t)) = assign(&tracks, pts.rep[c.cluster], eps)? {
                on[k].push((t, c.cluster, c.up));
            }
        }
        for (k, list) in on.into_iter().enumerate() {
            out[r][k] = spans_on(pair, pts, outlines, op, f, theirs, tracks[k], region, list)?;
        }
    }
    Ok(out)
}

/// The spans of `track` in face `f` of `op` under one rule, from the cuts
/// on it: (parameter, cluster, up), `region` the box (`lo`, `hi`) where
/// faces of both operands can be.
#[allow(clippy::too_many_arguments)]
fn spans_on(
    pair: &Pair,
    pts: &Points,
    outlines: &[Vec<Option<Outline>>; 2],
    op: Operand,
    f: usize,
    theirs: &Quadric,
    track: &Track,
    region: (DVec3, DVec3),
    list: Vec<(f64, usize, bool)>,
) -> Result<Vec<Span>, KernelError> {
    let shape = pair.get(op);
    let face = &shape.faces[f];
    let own = &shape.surfaces[face.surface];
    let eps = pts.eps;
    let inconsistent = |p: DVec3| {
        KernelError::InvalidResult(format!(
            "boolean: face {f} of the {} body meets the other body's surface inconsistently \
             near {}",
            name(op),
            show(p)
        ))
    };
    // Which way the track runs across the face: with the face's outward
    // normal `n` and the other surface's gradient `g`, the track runs along
    // `n × g` or against it, the same way all along it (the surfaces cross
    // without touching), and the boundary crossing to the positive side of
    // the other surface enters the face exactly when it runs along it.
    // Decided where the surfaces cross most steeply among the cuts.
    let mut way = 0.0f64;
    for &(t, _, _) in &list {
        let p = track.point(t);
        let n = own.gradient(p) * if face.outward { 1.0 } else { -1.0 };
        let v = track
            .deriv(t)
            .normalize_or_zero()
            .dot(n.cross(theirs.gradient(p)));
        if v.abs() > way.abs() {
            way = v;
        }
    }
    if !list.is_empty() && way.abs() < 1e-6 {
        return Err(KernelError::Unsupported(format!(
            "booleans of bodies whose faces meet at a grazing angle (design note 0011 §4.4) \
             near {}",
            show(pts.rep[list[0].1])
        )));
    }
    // The cuts at one cluster together: the boundary passing through a
    // point of the track, or touching it there (no net change).
    let mut groups: BTreeMap<usize, (f64, i32)> = BTreeMap::new();
    for (t, c, up) in list {
        let g = groups.entry(c).or_insert((t, 0));
        g.1 += if up == (way > 0.0) { 1 } else { -1 };
    }
    let mut ends: Vec<(f64, usize, i32)> = groups
        .into_iter()
        .filter(|(_, (_, n))| *n != 0)
        .map(|(c, (t, n))| (t, c, n))
        .collect();
    if let Some(e) = ends.iter().find(|e| e.2.abs() > 1) {
        return Err(inconsistent(pts.rep[e.1]));
    }
    ends.sort_by(|x, y| x.0.total_cmp(&y.0).then(x.1.cmp(&y.1)));
    if !track.closed() {
        // An end of an open track inside the face (a ray from a cone's
        // apex or from where it comes near the region, or the half of an
        // ellipse ending on its minor axis): the face holds the track from
        // or to there. No edge of the result ends there, as the end lies
        // beyond the region where both faces can be (§14), so it is open,
        // the other face bounding the part in both.
        let (t0, t1) = track.ends();
        let outline = outlines[slot(op)][f].as_ref();
        let inside = |t: f64| {
            if !t.is_finite() {
                return Ok(false);
            }
            let p = track.point(t);
            locate(shape, outline, face, p, eps).ok_or_else(|| {
                if matches!(*track, Track::Ray { .. }) {
                    // The apex or the ray's start beyond the region, on
                    // the face's boundary (a pointed cone's seam along
                    // the ray).
                    return through_apex(p);
                }
                KernelError::Unsupported(format!(
                    "booleans where a curve where faces meet ends on an edge (design note \
                     0012 §14) near {}",
                    show(p)
                ))
            })
        };
        if inside(t0)? {
            ends.insert(0, (f64::NEG_INFINITY, OPEN, 1));
        }
        if inside(t1)? {
            ends.push((f64::INFINITY, OPEN, -1));
        }
    }
    if ends.is_empty() {
        if !track.closed() {
            return Ok(Vec::new());
        }
        // Wholly in the face or wholly out of it; points of the track
        // within the tolerance of the boundary are in the closed face.
        // Decided where faces of the other body can be: far from there
        // the track can come within the tolerance of the boundary without
        // a cut, which no stretch depends on (found by review: a plane
        // nearly tangent to a cone along a ruling, their ellipse metres
        // long, grazing the edge of a plate far from the cone's faces).
        // Where every point tried there is on the boundary (a track along
        // a rim), by the quarters of the track.
        let outline = outlines[slot(op)][f].as_ref();
        let (lo, hi) = region;
        let quarters = (0..4).map(|k| track.point(k as f64 * FRAC_PI_2));
        let near = quarters
            .clone()
            .filter(|p| p.cmpge(lo).all() && p.cmple(hi).all())
            .chain([track.point(track.nearest((lo + hi) * 0.5).0)]);
        let sides = |probes: &mut dyn Iterator<Item = DVec3>| {
            let (mut inside, mut outside) = (false, false);
            for p in probes {
                match locate(shape, outline, face, p, eps) {
                    Some(true) => inside = true,
                    Some(false) => outside = true,
                    None => {}
                }
            }
            (inside, outside)
        };
        let (mut inside, mut outside) = sides(&mut near.into_iter());
        if !inside && !outside {
            (inside, outside) = sides(&mut quarters.into_iter());
        }
        return match (inside, outside) {
            (true, true) => Err(inconsistent(track.point(0.0))),
            (false, true) => Ok(Vec::new()),
            _ => Ok(vec![Span::Full]),
        };
    }
    let n = ends.len();
    let pairs = if track.closed() { n } else { n - 1 };
    for k in 0..pairs {
        let (x, y) = (ends[k], ends[(k + 1) % n]);
        if x.1 != OPEN && y.1 != OPEN && (track.point(x.0) - track.point(y.0)).length() <= eps {
            return Err(close(pts.rep[y.1]));
        }
    }
    // Entries and exits alternate, from an entry (on a closed track, any).
    let start = if track.closed() {
        ends.iter().position(|e| e.2 == 1).unwrap_or(0)
    } else {
        0
    };
    ends.rotate_left(start);
    if n % 2 != 0
        || ends
            .iter()
            .enumerate()
            .any(|(k, e)| e.2 != if k % 2 == 0 { 1 } else { -1 })
    {
        // Near a cut, not the open end of a track inside the face.
        let near = match ends.iter().find(|e| e.1 != OPEN) {
            Some(e) => pts.rep[e.1],
            None => track.point(track.ends().0.max(0.0)),
        };
        return Err(inconsistent(near));
    }
    Ok(ends
        .chunks(2)
        .map(|c| Span::Part {
            lo: (c[0].1, c[0].0),
            hi: (c[1].1, c[1].0),
        })
        .collect())
}

/// The stretches of `track` in both of two spans, not both the whole
/// track.
fn both(track: &Track, x: Span, y: Span, pts: &Points) -> Result<Vec<(End, End)>, KernelError> {
    let eps = pts.eps;
    // A part of a closed track as parameters increasing from its start.
    let unrolled = |lo: End, hi: End| (lo, (hi.0, lo.1 + (hi.1 - lo.1).rem_euclid(TAU)));
    match (x, y) {
        (Span::Full, Span::Full) => unreachable!("a whole track is imprinted whole"),
        (Span::Full, Span::Part { lo, hi }) | (Span::Part { lo, hi }, Span::Full) => {
            Ok(vec![unrolled(lo, hi)])
        }
        (Span::Part { lo: xl, hi: xh }, Span::Part { lo: yl, hi: yh }) => {
            if !track.closed() {
                let lo = later((xl.1, xl.0), (yl.1, yl.0), eps, pts)?;
                let hi = earlier((xh.1, xh.0), (yh.1, yh.0), eps, pts)?;
                if lo.0 < hi.0 && (lo.1 == OPEN || hi.1 == OPEN) {
                    return Err(KernelError::InvalidResult(
                        "boolean: a curve where faces meet ends inside both of them".into(),
                    ));
                }
                return Ok(if lo.1 != hi.1 && lo.0 < hi.0 {
                    vec![((lo.1, lo.0), (hi.1, hi.0))]
                } else {
                    Vec::new()
                });
            }
            // Positions measured from the start of `x`, which runs to
            // `lx`; `y` runs from `ya` to `yb`, and is met again a turn
            // back.
            let base = xl.1;
            let lx = (xh.1 - base).rem_euclid(TAU);
            let ya = (yl.1 - base).rem_euclid(TAU);
            let yb = ya + (yh.1 - yl.1).rem_euclid(TAU);
            // The greater (or lesser) of two positions of clusters: equal
            // for one cluster (its positions differ by whole turns, and the
            // nearer copy is meant), refused for two within the tolerance.
            let pick = |greater: bool, p: (f64, usize), q: (f64, usize)| {
                if p.1 == q.1 && (p.0 - q.0).abs() < 1.0 {
                    return Ok(p);
                }
                if p.1 != q.1 && (track.point(base + p.0) - track.point(base + q.0)).length() <= eps
                {
                    return Err(close(pts.rep[q.1]));
                }
                Ok(if (p.0 < q.0) == greater { q } else { p })
            };
            let mut out = Vec::new();
            for shift in [0.0, -TAU] {
                let s = pick(true, (0.0, xl.0), (ya + shift, yl.0))?;
                let e = pick(false, (lx, xh.0), (yb + shift, yh.0))?;
                if s.1 != e.1 && s.0 < e.0 {
                    out.push(((s.1, base + s.0), (e.1, base + e.0)));
                }
            }
            Ok(out)
        }
    }
}

/// The stretches of the tracks where crossing faces `fa` (of the first
/// operand) and `fb` (of the second) meet, at least one of them on a
/// curved surface, that lie in both (closed) faces; and the closed tracks
/// that lie wholly in both, with no edge of either face reaching them (a
/// ball's face around a hole in a plate).
pub(super) fn stretches(
    pair: &Pair,
    pts: &Points,
    outlines: &[Vec<Option<Outline>>; 2],
    contacts: &Contacts,
    fa: usize,
    fb: usize,
) -> Result<(Vec<Stretch>, Vec<usize>), KernelError> {
    let a = &pair.first.faces[fa];
    let b = &pair.second.faces[fb];
    // An edge touching the other face first: its cuts would not pair up
    // along the tracks.
    refuse_touches(pair, pts, outlines, fa, fb)?;
    let ids = contacts.tracks_of(a.surface, b.surface);
    let fs = face_spans(pair, pts, outlines, contacts, Operand::First, fa, fb, ids)?;
    let gs = face_spans(pair, pts, outlines, contacts, Operand::Second, fb, fa, ids)?;
    let mut out = Vec::new();
    let mut whole = Vec::new();
    for (k, &id) in ids.iter().enumerate() {
        let track = &contacts.tracks[id];
        let full = |spans: &[Vec<Vec<Span>>; 2]| {
            let mut all = spans.iter().flat_map(|rule| &rule[k]).peekable();
            all.peek().is_some() && all.all(|s| matches!(s, Span::Full))
        };
        if full(&fs) && full(&gs) {
            whole.push(id);
            continue;
        }
        for x in fs.iter().flat_map(|rule| &rule[k]) {
            for y in gs.iter().flat_map(|rule| &rule[k]) {
                if let (Span::Full, Span::Full) = (x, y) {
                    // Whole in both faces under one rule, not under the
                    // other: an edge touches the track.
                    return Err(KernelError::Unsupported(format!(
                        "booleans where an edge touches a closed curve where faces meet \
                         (design note 0011 §4.4), near {}",
                        show(track.point(0.0))
                    )));
                }
                for (lo, hi) in both(track, *x, *y, pts)? {
                    out.push(Stretch { track: id, lo, hi });
                }
            }
        }
    }
    Ok((out, whole))
}
