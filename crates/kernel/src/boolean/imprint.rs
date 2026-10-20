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

use std::cell::OnceCell;
use std::cmp::Ordering;
use std::collections::BTreeMap;
use std::f64::consts::{FRAC_PI_2, TAU};

use glam::{DVec2, DVec3};
use rustc_hash::FxHashMap;

use super::chart::Outline;
use super::flat_outline::{FlatOutline, FlatOutlines};
use super::meet::{self, Meet, Track};
use super::points::{self, Points};
use super::quadric::{self, Quadric};
use super::shape::{Face, Shape};
use super::{name, other, show, slot, Pair};
use crate::bvh::Bvh;
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
    /// surface. `zero` says, per face (the first operand's first), whether
    /// its boundary has a mark within the tolerance of the other's surface:
    /// only then can it touch the other face without crossing it.
    Crossing { zero: [bool; 2] },
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
        return Ok(Relation::Crossing { zero: [az, bz] });
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
        Contact::Tracks(_) => Ok(Relation::Crossing { zero: [az, bz] }),
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
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub(super) enum Mark {
    Vertex(usize),
    Station(usize, usize),
}

/// A face's boundary touching a surface of the other operand from one
/// side ([`touches`]): the cluster, the mark of the boundary there, and the
/// side of the surface the boundary lies on next to it (by the surface's
/// signed distance).
#[derive(Clone, Copy, Debug)]
struct Touch {
    cluster: usize,
    mark: Mark,
    sign: i8,
}

/// A point of a face's boundary touching a face of the other operand
/// ([`touches_in`]): the touched face (`op`, `face`), the cluster, the
/// touching mark (of the other operand's boundary), and the side of the
/// touched face the touching boundary lies on, +1 outside the touched body
/// and −1 inside it.
#[derive(Clone, Copy, Debug)]
pub(super) struct Touched {
    pub op: Operand,
    pub face: usize,
    pub cluster: usize,
    pub mark: Mark,
    pub side: i8,
}

/// The places where the boundary of face `f` of `op` touches `surface` of
/// the other operand from one side: a station or a vertex on it, or a run
/// of them, between two strictly on one side along a loop (tangent
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
///
/// The marks of a run that touches, other than those inside an edge lying
/// on the surface, must all be in one cluster, or, in a run with lying
/// edges, in the clusters of those edges' ends: otherwise an edge runs
/// within the tolerance of the surface over a stretch, which the
/// tolerance cannot place, and the boolean is refused (design note 0012
/// §18.2). On a plane no straight edge can: two of its vertices within the
/// tolerance make it lie on the plane.
fn touches(
    pair: &Pair,
    pts: &Points,
    op: Operand,
    f: usize,
    surface: usize,
    near: &[usize],
) -> Result<Vec<Touch>, KernelError> {
    let shape = pair.get(op);
    let face = &shape.faces[f];
    let plane = is_plane(&pair.get(other(op)).surfaces[surface]);
    let mut out = Vec::new();
    // Along a run of coedges: each mark's side, where it is, whether its
    // edge lies wholly on the surface, and its edge.
    let mut seq: Vec<(i8, Mark, bool, usize)> = Vec::new();
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
                seq.extend(marks.into_iter().map(|(sign, m)| (sign, m, on, e)));
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
                if b + 1 < seq.len() && seq[a - 1].0 == seq[b + 1].0 && !crossed {
                    let cluster = |m: Mark| match m {
                        Mark::Vertex(v) => Some(pts.vertex_cluster(op, v)),
                        Mark::Station(e, s) => pts.station_cluster(pair, op, e, surface, s),
                    };
                    // The clusters the run may hold: its lying edges' ends,
                    // or else its first mark's.
                    let allowed: Vec<usize> = if lying {
                        seq[a..=b]
                            .iter()
                            .filter(|x| x.2)
                            .flat_map(|x| {
                                let (va, vb) = shape.edges[x.3];
                                [pts.vertex_cluster(op, va), pts.vertex_cluster(op, vb)]
                            })
                            .collect()
                    } else {
                        let first = cluster(seq[a].1).ok_or_else(|| match seq[a].1 {
                            Mark::Station(e, _) => missing(op, e),
                            Mark::Vertex(_) => missing(op, seq[a].3),
                        })?;
                        vec![first]
                    };
                    for x in &seq[a..=b] {
                        // A station inside a lying edge is on the surface
                        // with it.
                        if x.2 && matches!(x.1, Mark::Station(..)) {
                            continue;
                        }
                        if !cluster(x.1).is_some_and(|c| allowed.contains(&c)) {
                            let at = match x.1 {
                                Mark::Vertex(v) => shape.points[v],
                                Mark::Station(e, s) => {
                                    let st = pts
                                        .cached_stations(op, e, surface)
                                        .ok_or_else(|| missing(op, e))?;
                                    shape.curve(e).0.point(st[s].t)
                                }
                            };
                            return Err(KernelError::InvalidInput(format!(
                                "an edge of one body runs within the tolerance of a face of the \
                                 other over a stretch without lying on it, near {}",
                                show(at)
                            )));
                        }
                    }
                    if !lying {
                        out.push(Touch {
                            cluster: allowed[0],
                            mark: seq[a].1,
                            sign: seq[a - 1].0,
                        });
                    }
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

/// Records where the boundary of either face touches the other face from
/// one side (design note 0012 §18.2): `zero` says, per face (the first's
/// first), whether its boundary has a mark on the other's surface at all.
/// A face with an outline keeps only points in it (closed); on a flat face
/// the split places each point exactly.
pub(super) fn touches_in(
    pair: &Pair,
    pts: &Points,
    outlines: &[Vec<Option<Outline>>; 2],
    fa: usize,
    fb: usize,
    zero: [bool; 2],
    out: &mut Vec<Touched>,
) -> Result<(), KernelError> {
    let eps = pts.eps;
    for (op, f, g) in [(Operand::First, fa, fb), (Operand::Second, fb, fa)] {
        if !zero[slot(op)] {
            continue;
        }
        let (mine, theirs) = (pair.get(op), pair.get(other(op)));
        let gf = &theirs.faces[g];
        let q = &theirs.surfaces[gf.surface];
        let mut near = Vec::new();
        mine.faces[f].near_surface(q, eps, &mut near);
        for t in touches(pair, pts, op, f, gf.surface, &near)? {
            let outline = outlines[slot(other(op))][g].as_ref();
            if outline.is_some()
                && locate(theirs, outline, gf, pts.rep[t.cluster], eps) == Some(false)
            {
                continue;
            }
            out.push(Touched {
                op: other(op),
                face: g,
                cluster: t.cluster,
                mark: t.mark,
                side: if gf.outward { t.sign } else { -t.sign },
            });
        }
    }
    Ok(())
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
        for t in touches(pair, pts, op, f, gf.surface, &near)? {
            let p = pts.rep[t.cluster];
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
///
/// Where one face's interval is a run of its edges lying on the other's
/// plane, and the other's is not, the stretch is that run (design note
/// 0012 §18.2: a lying edge is imprinted as the meeting, and merges with
/// the edge's own segments), so an end of it that the other face's
/// boundary gives must be a point of the run. A run within the tolerance
/// of a plane can be a few tolerances from the planes' exact common line
/// where they meet at a slant, and the other face's boundary crosses that
/// line, not the run: the run then leaves the other face at one point and
/// the faces' meeting ends at another, more than the tolerance from the
/// run. Beside the run on its face's side, or within the tolerance of its
/// line, the stretch would bound a sliver of that face between it and the
/// run that no face of the other body bounds (found by the contact fuzz:
/// the edge of a prism resting within the tolerance on a face and crossing
/// its boundary, its faces rising at a slant); that is refused as input
/// ([`unplaced_leaving`]), as the face's planar split could not place it
/// either. Beyond that off its face, the stretch lies off that face, whose
/// split drops it, and it is kept: it splits the other face into pieces of
/// one class, which the result joins again, or, where it ends inside that
/// face, it is a bridge there carrying no edge, which
/// [`super::lying_pieces`] places on the run as the contact (golden 22). The same holds where the run
/// is only a part of the face's interval, the line going on into the face
/// past a reflex corner at an end of the run or the run being an edge of a
/// hole, for an end of the stretch beside the run (found by review: the
/// stretch there bounded a fin of the two faces thinner than the
/// tolerance). The side is taken from the run's edge nearest the end,
/// since a run of several edges can bend within the tolerance of the plane
/// (found by review). Where both intervals are runs (an edge of each face
/// on the other's plane, the planes meeting at a slight angle), the ends
/// are the runs' own vertices, not crossings of a plane, and are not
/// checked here: the faces then lie within the tolerance of each other
/// over a band, which classifying their pieces refuses. After those
/// checks, the strips between the line and the chains of either face's
/// boundary lying on the other's plane are checked ([`held_strips`]), for
/// each face whose `zero` ([`Relation::Crossing`]) says its boundary has a
/// mark within the tolerance of the other's plane: a face with none has no
/// chain. The checks of one boolean share `strips`, and place points in a
/// face with curved edges by its chart's outline in `outlines`.
pub(super) fn plane_meeting(
    pair: &Pair,
    pts: &Points,
    outlines: &[Vec<Option<Outline>>; 2],
    strips: &mut StripSearch,
    fa: usize,
    fb: usize,
    zero: [bool; 2],
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
    // Each face's runs lying on the other's plane, found once for the pair
    // when a stretch first needs them.
    let mut found: Option<[Runs; 2]> = None;
    let mut out = Vec::new();
    for x in &sa {
        // The intervals of `b` within ε of `x` (each face's intervals are
        // disjoint and in order along the line): no other can meet it.
        let first = sb.partition_point(|y| y.hi.0 < x.lo.0 - eps);
        for y in sb[first..].iter().take_while(|y| y.lo.0 <= x.hi.0 + eps) {
            let lo = later(x.lo, y.lo, eps, pts)?;
            let hi = earlier(x.hi, y.hi, eps, pts)?;
            if lo.1 != hi.1 && lo.0 < hi.0 {
                let spans = [(Operand::First, fa, x), (Operand::Second, fb, y)];
                let runs = found.get_or_insert_with(|| {
                    [
                        (Operand::First, fa, b.surface),
                        (Operand::Second, fb, a.surface),
                    ]
                    .map(|(op, f, surface)| {
                        Runs::new(lying_runs(pair, pts, op, f, surface, &[]), dir)
                    })
                });
                for (k, &(op, f, own)) in spans.iter().enumerate() {
                    // Where the other span is a run as well, the ends it
                    // gives are its run's own vertices.
                    if runs[1 - k].exact(spans[1 - k].2) {
                        continue;
                    }
                    let theirs = [lo.1, hi.1]
                        .into_iter()
                        .filter(|&c| c != own.lo.1 && c != own.hi.1);
                    for c in theirs {
                        let p = pts.rep[c];
                        for run in runs[k].beside(p.dot(dir), eps) {
                            let (r0, r1) = run.stretch;
                            if p.dot(dir) < r0 - eps
                                || p.dot(dir) > r1 + eps
                                || run.holds(c, pts, op)
                            {
                                continue;
                            }
                            let run = &run.run;
                            // The side of the run's face the end is on: the
                            // direction into the face across the run's
                            // edge nearest it is the face's outward normal
                            // crossed with the edge's direction as its
                            // loop runs it.
                            let (a, b) = run.nearest(p);
                            let face = &pair.get(op).faces[f];
                            let into = face.normal().cross(b - a).normalize();
                            if (p - a).dot(into) > -eps {
                                return Err(unplaced_leaving(p));
                            }
                        }
                    }
                }
                out.push((lo.1, hi.1));
            }
        }
    }
    held_strips(
        pair,
        pts,
        outlines,
        strips,
        dir,
        [
            (Operand::First, fa, &sa, zero[0]),
            (Operand::Second, fb, &sb, zero[1]),
        ],
        &out,
    )?;
    Ok(out)
}

/// How far past the tolerance ε the strip search reaches, as a fraction of
/// ε ([`held_strip`]): it finds every region of a strip that the other
/// face holds deeper than ρ = ε (1 + `DEEP`), and refuses only where it
/// holds one deeper than ε, so the difference is far below the tolerance
/// and far above rounding.
const DEEP: f64 = 1.0 / 1024.0;

/// Steps of the strip searches of one boolean, at most, by default
/// ([`held_strip`]; `Limits::strip_steps`): a stretch, an edge or a cut
/// taken, a box tested, a clip of a segment or of a curve's piece against
/// a piece of a strip, and a side, a loop or a chord of a face's outline
/// met in placing a point are each one step, and building a flat face's
/// outline is one for each of its sides.
pub(super) const STRIP_STEPS: usize = 1 << 22;

/// What [`STRIP_STEPS`] counts, for the refusal.
pub(super) const STRIP_WORK: &str = "steps to test where a face runs within the tolerance of \
                                     the other body's surface (design note 0012 §20)";

/// The steps the strip searches may still take, of `limit`.
struct Budget {
    left: usize,
    limit: usize,
}

impl Budget {
    fn new(limit: usize) -> Budget {
        Budget { left: limit, limit }
    }

    fn spend(&mut self, n: usize) -> Result<(), KernelError> {
        self.left =
            (self.left.checked_sub(n)).ok_or_else(|| super::too_large(STRIP_WORK, self.limit))?;
        Ok(())
    }
}

/// What the strip checks of one boolean share: the flat faces' outlines,
/// each built when a strip first needs it, and the steps left, so the
/// searches of all its strips together are bounded.
pub(super) struct StripSearch {
    flat: FlatOutlines,
    budget: Budget,
}

impl StripSearch {
    /// No outline built and no step taken yet, of `steps` at most.
    pub(super) fn new(steps: usize) -> StripSearch {
        StripSearch {
            flat: FlatOutlines::new(),
            budget: Budget::new(steps),
        }
    }
}

/// The check of the strips of a pair of crossing faces on planes (design
/// note 0012 §20): for each face, given with its operand and its intervals
/// on the planes' common line L (direction `dir`) and whether its boundary
/// has a mark within the tolerance of the other's plane, each chain of its
/// boundary lying on the other's plane ([`strip_runs`]; none without such
/// a mark), against the pair's `stretches` beside it ([`held_strip`]). A
/// run that is one of the face's intervals whole is tested only over the
/// parts of its span along L that no stretch covers (each stretch grown by
/// the tolerance): where none does, the run is no meeting, and nothing
/// separates its strip from the other face (§21: a notch's wall within the
/// tolerance of a box's side over a triangle beside the notch's root, the
/// root lying on the side as the wall's whole interval with no stretch,
/// failed as kept pieces that do not close).
///
/// Then the notch tips of either face across L ([`strip_runs`], §21): a
/// vertex V of the face within the tolerance of the other's plane P, at a
/// reflex corner, both its edges crossing L. The tolerance puts V on P, so
/// its edges meet P only at V, and the face, whose interior wraps around
/// V's other side, crosses P through V: the faces' meeting must pass
/// through V. Where no stretch does, as an end or as a cluster its split
/// takes ([`Points::on_split`]; a tip within the tolerance of a stretch
/// bends it and is answered), V's cluster is no point of the other face's
/// boundary (as the clustering decided it), and V's foot on P is not
/// outside that face, the boolean is refused before any face is split:
/// the meeting along L would leave V touching P from the side its edges
/// leave to while the face lies on both, and a meeting bent through V
/// would leave a strip of the face between it and L lying within the
/// tolerance of the other face, which holds it (§20's fin). The tip with
/// the least coordinates is reported, so the two orders of the operands
/// report one point.
fn held_strips(
    pair: &Pair,
    pts: &Points,
    outlines: &[Vec<Option<Outline>>; 2],
    strips: &mut StripSearch,
    dir: DVec3,
    sides: [(Operand, usize, &[Interval], bool); 2],
    stretches: &[(usize, usize)],
) -> Result<(), KernelError> {
    let eps = pts.eps;
    // The stretches along L, sorted when a chain or a tip first needs them:
    // a chain takes those over its span along L (and ε more), whose
    // meetings are the only ones that can cut its strip.
    let mut along: Option<Spans> = None;
    let mut tips = Vec::new();
    for (k, &(op, f, ivs, zero)) in sides.iter().enumerate() {
        if !zero {
            continue;
        }
        let (them, g, _, _) = sides[1 - k];
        let surface = pair.get(them).faces[g].surface;
        let (chains, corners) = strip_runs(pair, pts, op, f, surface, ivs);
        tips.extend(corners.into_iter().map(|v| (op, v, them, g)));
        for chain in chains {
            // A chain whose strip cannot come within the search's reach of
            // g has no point g holds, and is passed before any step: the
            // steps of a face's chains would otherwise grow with the
            // coplanar faces of the other body (found by review).
            if !strip_near(pair, [(op, f), (them, g)], dir, &chain.points, eps) {
                continue;
            }
            // Where the clustering put g's boundary on the chain: a flap's
            // corner that is a point of it, and the corners of g on a run
            // that is one of f's intervals whole. The settling of touched
            // points and the corners decide the contact there, as for a
            // notch tip, so the strip is tested less its points within
            // their reach (§21).
            let mut spots = Vec::new();
            if let Some(v) = chain.corner {
                strips.budget.spend(1)?;
                let c = pts.vertex_cluster(op, v);
                if on_boundary(pair, pts, (them, g), c) {
                    spots.push(pts.rep[c]);
                }
            }
            if chain.kind == Kind::Whole {
                corners_on(
                    pair,
                    pts,
                    op,
                    &chain.edges,
                    (them, g),
                    &mut spots,
                    &mut strips.budget,
                )?;
            }
            let along = sorted(&mut along, pts, dir, stretches, &mut strips.budget)?;
            let (lo, hi) = (chain.points.iter())
                .fold((f64::INFINITY, f64::NEG_INFINITY), |(l, h), p| {
                    (l.min(p.dot(dir)), h.max(p.dot(dir)))
                });
            let (over, passed) = along.meeting(lo - eps, hi + eps);
            strips.budget.spend(1 + passed)?;
            // In the pair's order, which the strip's own sort keeps for
            // stretches starting at one place.
            let mut beside: Vec<usize> = over.collect();
            beside.sort_unstable();
            let beside: Vec<(usize, usize)> = beside.into_iter().map(|i| stretches[i]).collect();
            // A run that is one of the face's intervals whole is the meeting
            // only where a stretch of the pair runs along it (§18.2); there
            // what lies between it and L is left to classification, and its
            // strip is tested over the rest of its span alone.
            let gaps = if chain.kind == Kind::Whole {
                let spans = beside.iter().map(|&(c, d)| {
                    let (s, t) = (pts.rep[c].dot(dir), pts.rep[d].dot(dir));
                    (s.min(t), s.max(t))
                });
                let gaps = uncovered((lo, hi), spans, eps);
                if gaps.is_empty() {
                    continue;
                }
                Some(gaps)
            } else {
                None
            };
            held_strip(
                pair,
                pts,
                outlines,
                strips,
                dir,
                [(op, f), (them, g)],
                &chain,
                &beside,
                gaps.as_deref(),
                &spots,
            )?;
        }
    }
    let mut least: Option<DVec3> = None;
    for (op, v, them, g) in tips {
        let p = pair.get(op).points[v];
        // A tip whose foot on P lies outside g's box, beyond the margins,
        // is placed outside g ([`tip_refused`]'s third test), and is
        // passed before any step, as a chain far from g is.
        let gf = &pair.get(them).faces[g];
        let foot = p - gf.normal() * (p - gf.frame.origin).dot(gf.normal());
        let pad = DVec3::splat(4.0 * eps);
        if !(foot.cmpge(gf.lo - pad).all() && foot.cmple(gf.hi + pad).all()) {
            continue;
        }
        let along = sorted(&mut along, pts, dir, stretches, &mut strips.budget)?;
        if tip_refused(
            pair,
            pts,
            outlines,
            strips,
            along,
            stretches,
            dir,
            (op, v),
            (them, g),
        )? && least.is_none_or(|q| {
            (p.x.total_cmp(&q.x))
                .then(p.y.total_cmp(&q.y))
                .then(p.z.total_cmp(&q.z))
                .is_lt()
        }) {
            least = Some(p);
        }
    }
    match least {
        Some(p) => Err(KernelError::InvalidInput(format!(
            "{}, near {}",
            super::RUNS_WITHIN,
            show(p)
        ))),
        None => Ok(()),
    }
}

/// Whether the strip of a chain of face `f` of `op` with `points`
/// ([`held_strip`]) may come within the search's reach of face `g` of
/// `them`: the strip lies in `f`'s plane between the points and their feet
/// on L (along `f`'s plane, square to `dir`), so in the box of those, and
/// its points lie no farther from P than the chain's farthest; the search
/// looks within twice ρ of their feet on P. True also where the planes are
/// too nearly parallel to find the feet, for [`held_strip`] to decide.
fn strip_near(
    pair: &Pair,
    [(op, f), (them, g)]: [(Operand, usize); 2],
    dir: DVec3,
    points: &[DVec3],
    eps: f64,
) -> bool {
    let (face, gf) = (&pair.get(op).faces[f], &pair.get(them).faces[g]);
    let ng = gf.normal();
    // Moving along `across` in f's plane changes the distance from P by
    // `sin` a unit.
    let across = face.normal().cross(dir);
    let sin = across.dot(ng);
    if sin.abs() <= 1e-10 {
        return true;
    }
    let (mut lo, mut hi) = (DVec3::splat(f64::INFINITY), DVec3::splat(f64::NEG_INFINITY));
    let mut off = 0.0f64;
    for &p in points {
        let d = (p - gf.frame.origin).dot(ng);
        let foot = p - across * (d / sin);
        (lo, hi) = (lo.min(p).min(foot), hi.max(p).max(foot));
        off = off.max(d.abs());
    }
    // Twice ρ, and as much again for rounding.
    let pad = DVec3::splat(off + 4.0 * eps * (1.0 + DEEP));
    (lo - pad).cmple(gf.hi).all() && gf.lo.cmple(hi + pad).all()
}

/// The pair's `stretches` sorted by their spans along `dir` ([`Spans`]),
/// kept in `along` once made (a step for each).
fn sorted<'s>(
    along: &'s mut Option<Spans>,
    pts: &Points,
    dir: DVec3,
    stretches: &[(usize, usize)],
    budget: &mut Budget,
) -> Result<&'s Spans, KernelError> {
    if along.is_none() {
        budget.spend(stretches.len())?;
    }
    Ok(along.get_or_insert_with(|| {
        Spans::new(stretches.iter().map(|&(c, d)| {
            let (s, t) = (pts.rep[c].dot(dir), pts.rep[d].dot(dir));
            (s.min(t), s.max(t))
        }))
    }))
}

/// Whether the notch tip at vertex `v` of `op` ([`held_strips`]) is
/// refused against face `g` of `them`: no stretch of the pair (`along`,
/// sorted from `stretches`) whose span along L comes within the tolerance
/// of the tip's takes the tip's cluster, the cluster holds no vertex of
/// `g` and lies on no edge of `g`, and the tip's foot on `g`'s plane is
/// not outside `g` (to the rounding margin of the strip search).
#[allow(clippy::too_many_arguments)]
fn tip_refused(
    pair: &Pair,
    pts: &Points,
    outlines: &[Vec<Option<Outline>>; 2],
    strips: &mut StripSearch,
    along: &Spans,
    stretches: &[(usize, usize)],
    dir: DVec3,
    (op, v): (Operand, usize),
    (them, g): (Operand, usize),
) -> Result<bool, KernelError> {
    let eps = pts.eps;
    let (flat, budget) = (&strips.flat, &mut strips.budget);
    let p = pair.get(op).points[v];
    let c = pts.vertex_cluster(op, v);
    // The stretches whose spans come within the tolerance of the cluster's
    // place, which the split measures (its representative may lie up to
    // twice the tolerance from V).
    let t = pts.rep[c].dot(dir);
    let (mut over, passed) = along.meeting(t - eps, t + eps);
    budget.spend(1 + passed)?;
    if over.any(|j| pts.on_split(stretches[j].0, stretches[j].1, c)) {
        return Ok(false);
    }
    // A point of g's boundary: the settling of touched points and the
    // corners decide it.
    budget.spend(1)?;
    if on_boundary(pair, pts, (them, g), c) {
        return Ok(false);
    }
    let theirs = pair.get(them);
    let gf = &theirs.faces[g];
    let ng = gf.normal();
    let foot = p - ng * (p - gf.frame.origin).dot(ng);
    let against = match &outlines[slot(them)][g] {
        Some(o) => Against::Chart(o),
        None => {
            let (o, fresh) = flat.get(pair, them, g);
            if fresh {
                budget.spend(o.len())?;
            }
            Against::Flat(o)
        }
    };
    Ok(place(&against, gf, foot, eps * DEEP, budget)? != Some(false))
}

/// Whether cluster `c` is a point of face `g` of `them`'s boundary, as the
/// clustering decided it: it holds a corner of `g` ([`holds_corner`]), or
/// lies on an edge bounding `g` between its ends.
fn on_boundary(pair: &Pair, pts: &Points, (them, g): (Operand, usize), c: usize) -> bool {
    let theirs = pair.get(them);
    holds_corner(pair, pts, (them, g), c)
        || (pts.inside_edge(them, c))
            .is_some_and(|e| theirs.edge_faces(e).iter().any(|&(h, _)| h == g))
}

/// Whether cluster `c` holds a corner of face `g` of `them`: a vertex of an
/// edge bounding `g`.
fn holds_corner(pair: &Pair, pts: &Points, (them, g): (Operand, usize), c: usize) -> bool {
    let theirs = pair.get(them);
    let vertex = match them {
        Operand::First => pts.vertices_of(c).0,
        Operand::Second => pts.vertices_of(c).1,
    };
    vertex.is_some_and(|u| {
        (theirs.edges_at(u).iter()).any(|&e| theirs.edge_faces(e).iter().any(|&(h, _)| h == g))
    })
}

/// The places of the clusters on `edges` of `op` (straight), at their
/// ends or between, that hold a corner of face `g` of `them`
/// ([`holds_corner`]), each once, added to `spots`. Such a cluster's place
/// lies within twice the tolerance of the corner, so of `g`'s box, and its
/// foot on the edge's line within twice the tolerance more (an end's place
/// within twice it of the end): only the clusters whose feet lie in `g`'s
/// box grown by that are looked at, found by their order along the edge
/// ([`Points::on_edge`]). A step for each edge and each cluster looked at
/// (found by review: looking at every cluster on a run in each pair took
/// steps growing with the square of the coplanar faces whose corners it
/// passes).
fn corners_on(
    pair: &Pair,
    pts: &Points,
    op: Operand,
    edges: &[usize],
    (them, g): (Operand, usize),
    spots: &mut Vec<DVec3>,
    budget: &mut Budget,
) -> Result<(), KernelError> {
    let (shape, gf) = (pair.get(op), &pair.get(them).faces[g]);
    // Four times the tolerance, and one more for rounding.
    let pad = DVec3::splat(5.0 * pts.eps);
    let (lo, hi) = (gf.lo - pad, gf.hi + pad);
    let mut held = Vec::new();
    for &e in edges {
        budget.spend(1)?;
        let (a, b) = shape.edge_points(e);
        let d = b - a;
        let Some((t0, t1)) = clip_line(a, d, (lo, hi)) else {
            continue;
        };
        let on = pts.on_edge(op, e);
        let len2 = d.length_squared();
        let at = |c: usize| (pts.rep[c] - a).dot(d) / len2;
        let first = on.partition_point(|&c| at(c) < t0);
        let end = first + on[first..].partition_point(|&c| at(c) <= t1);
        budget.spend(end - first)?;
        held.extend(
            on[first..end]
                .iter()
                .copied()
                .filter(|&c| holds_corner(pair, pts, (them, g), c)),
        );
    }
    held.sort_unstable();
    held.dedup();
    spots.extend(held.into_iter().map(|c| pts.rep[c]));
    Ok(())
}

/// The parameters `t` of the points `a + t d` of a line within the closed
/// box `lo`–`hi`, as an interval, if any (`d` not zero).
fn clip_line(a: DVec3, d: DVec3, (lo, hi): (DVec3, DVec3)) -> Option<(f64, f64)> {
    let (mut t0, mut t1) = (f64::NEG_INFINITY, f64::INFINITY);
    for i in 0..3 {
        if d[i] == 0.0 {
            if a[i] < lo[i] || a[i] > hi[i] {
                return None;
            }
        } else {
            let (u, v) = ((lo[i] - a[i]) / d[i], (hi[i] - a[i]) / d[i]);
            (t0, t1) = (t0.max(u.min(v)), t1.min(u.max(v)));
        }
    }
    (t0 <= t1).then_some((t0, t1))
}

/// A chain of a face's boundary lying on the other operand's plane, whose
/// strip [`held_strip`] tests ([`strip_runs`]).
struct Chain {
    /// The places of its points in order: one for a vertex alone, a closed
    /// loop's first again at its end, and for a flap the crossing of the
    /// corner's edge before it with the line, the corner, and the crossing
    /// of its edge after it.
    points: Vec<DVec3>,
    kind: Kind,
    /// A run's edges, in loop order (none for a loop or a vertex).
    edges: Vec<usize>,
    /// A flap's corner, the vertex across the line.
    corner: Option<usize>,
}

/// What a chain is ([`strip_runs`]), which decides where its strip is
/// tested ([`held_strip`]).
#[derive(Clone, Copy, Debug, PartialEq)]
enum Kind {
    /// A run of edges lying on the plane: on both sides of the planes'
    /// common line.
    Run,
    /// A run that is one of the face's intervals on the line whole: on
    /// both sides, over the parts of its span that no stretch covers.
    Whole,
    /// A loop lying on the plane whole (a hole) or an apex: on the side of
    /// the line where the layer between the faces lies in the face's own
    /// body.
    Inner,
    /// The flap of a corner across the line at a convex corner (design
    /// note 0012 §21): the triangle between the corner's edges and the
    /// line, less its points within the tolerance of the edges, where the
    /// layer between it and the other face lies outside the face's body.
    Flap,
}

/// The parts of the span `lo`..`hi` along the planes' common line that no
/// span of `spans`, each grown by `eps` at both ends (its end clusters'
/// reach), covers, in order along the line.
fn uncovered(
    (lo, hi): (f64, f64),
    spans: impl Iterator<Item = (f64, f64)>,
    eps: f64,
) -> Vec<(f64, f64)> {
    let mut spans: Vec<(f64, f64)> = spans.map(|(s, t)| (s - eps, t + eps)).collect();
    spans.sort_by(|p, q| p.0.total_cmp(&q.0));
    let mut gaps = Vec::new();
    let mut at = lo;
    for (s, t) in spans {
        if at >= hi {
            break;
        }
        if s > at {
            gaps.push((at, s.min(hi)));
        }
        at = at.max(t);
    }
    if at < hi {
        gaps.push((at, hi));
    }
    gaps
}

/// The chains of face `f` of `op`'s boundary lying on plane `surface` of
/// the other operand whose strips [`held_strip`] tests, and its notch tips
/// across the planes' common line L that [`held_strips`] tests:
///
/// - each run of its straight edges lying on the plane ([`lying_runs`]),
///   on either side, marked whole where it is one of its intervals `ivs`
///   on L whole: §18.2 takes such a run itself as the meeting where a
///   stretch of the pair runs along it, so no line is imprinted beside it
///   there and what lies between it and the other face is left to
///   classification, and [`held_strips`] tests its strip only over the
///   rest of its span;
/// - each loop lying on the plane whole, closed: a hole of the face (an
///   outer loop lies there whole only where the whole face does, which
///   [`relation`] takes for coplanar or refuses before this);
/// - the vertices within the tolerance of the plane whose two edges are
///   straight, sorted in one pass (design note 0012 §21), with `d` a
///   vertex's distance from the plane: an apex where `d` is not 0 and both
///   edges run away from the plane (the other end of each farther from it
///   than the vertex, on its side), at a reflex corner of the face, where
///   the boundary touches the plane without crossing it and the face lies
///   on both sides of the vertex's way to L (a neck); and a corner across
///   L where `d` is not 0 and both neighbours lie beyond the tolerance on
///   the plane's other side, so that both edges cross L and, by the
///   signs, the boundary meets the plane at the vertex alone: a notch tip
///   at a reflex corner, the face crossing the plane through it, or at a
///   convex corner the flap of the face across L between its edges. The
///   two are exclusive: neighbours beyond the tolerance on the vertex's
///   side are farther than it from the plane, so for them the first is
///   exactly `d` on their side and the second `d` on the other.
///
/// Only the coedges near the plane are looked at.
fn strip_runs(
    pair: &Pair,
    pts: &Points,
    op: Operand,
    f: usize,
    surface: usize,
    ivs: &[Interval],
) -> (Vec<Chain>, Vec<usize>) {
    let shape = pair.get(op);
    let face = &shape.faces[f];
    let plane = &pair.get(other(op)).surfaces[surface];
    let lies = |li: usize, i: usize| {
        let (e, _) = face.loops[li][i];
        let vs = &face.verts[li];
        shape.straight[e]
            && pts.sign(op, vs[i], surface) == 0
            && pts.sign(op, vs[(i + 1) % vs.len()], surface) == 0
    };
    let mut runs: Vec<Chain> = (lying_runs(pair, pts, op, f, surface, ivs).into_iter())
        .map(|run| Chain {
            points: run.points,
            kind: if run.exact { Kind::Whole } else { Kind::Run },
            edges: run.edges,
            corner: None,
        })
        .collect();
    let mut tips = Vec::new();
    let mut near = Vec::new();
    face.near_surface(plane, pts.eps, &mut near);
    let mut seen = vec![false; face.loops.len()];
    for k in near {
        let (li, i) = face.coedge_at(k);
        let vs = &face.verts[li];
        let n = vs.len();
        let before = (i + n - 1) % n;
        if !seen[li] && lies(li, i) {
            seen[li] = true;
            if (0..n).all(|j| lies(li, j)) {
                let mut points: Vec<DVec3> = vs.iter().map(|&v| shape.points[v]).collect();
                points.push(points[0]);
                runs.push(Chain {
                    points,
                    kind: Kind::Inner,
                    edges: Vec::new(),
                    corner: None,
                });
            }
        }
        // The vertex the coedge starts at.
        let (v, a, b) = (vs[i], vs[before], vs[(i + 1) % n]);
        let p = shape.points[v];
        let d = plane.distance(p);
        if pts.sign(op, v, surface) != 0
            || d == 0.0
            || !shape.straight[face.loops[li][before].0]
            || !shape.straight[face.loops[li][i].0]
        {
            continue;
        }
        let away = |u: usize| {
            let du = plane.distance(shape.points[u]);
            du.signum() == d.signum() && du.abs() > d.abs()
        };
        let beyond = -d.signum() as i8;
        let across = pts.sign(op, a, surface) == beyond && pts.sign(op, b, surface) == beyond;
        // Reflex: the loop turns clockwise about the face's outward normal
        // (holes run clockwise, so this holds for their corners too).
        let turn = (p - shape.points[a]).cross(shape.points[b] - p);
        let turn = turn.dot(face.normal());
        if away(a) && away(b) {
            if turn < 0.0 {
                runs.push(Chain {
                    points: vec![p],
                    kind: Kind::Inner,
                    edges: Vec::new(),
                    corner: None,
                });
            }
        } else if across && turn < 0.0 {
            tips.push(v);
        } else if across && turn > 0.0 {
            // Where each edge crosses the plane: on L.
            let cross = |u: usize| {
                let q = shape.points[u];
                p + (q - p) * (d / (d - plane.distance(q)))
            };
            runs.push(Chain {
                points: vec![cross(a), p, cross(b)],
                kind: Kind::Flap,
                edges: Vec::new(),
                corner: Some(v),
            });
        }
    }
    (runs, tips)
}

/// The check of the strip of `run`, a chain of face `f` of `op`'s boundary
/// lying on the plane P of face `g` of `them` ([`strip_runs`]), which
/// crosses `f`'s plane along the line L (direction `dir`; design note 0012
/// §20). The strip is the region of `f`'s plane swept from the chain to
/// the meeting the imprint makes, on the side of each of the chain's
/// segments where `f` lies toward L ([`strip_pieces`]; a lone vertex
/// sweeps a segment, [`point_pieces`]), cut to `f`: each point of the
/// chain to its foot on L, and where one of `stretches` (the pair's, those
/// over the chain's span along L) lies between the chain and L, no farther
/// than that stretch
/// ([`beyond_meeting`]); less its points within the tolerance of L or of
/// the stretch; and, for a run that is one of `f`'s intervals whole, only
/// over `gaps`, the parts of its span along L (as positions along `dir`)
/// that no stretch covers. A flap's is the triangle between the corner's
/// edges and L, less its points within the tolerance of either edge
/// ([`flap_piece`]; its points on L are kept, since no meeting lies
/// there), cut by the stretches likewise. Its points within three times
/// the tolerance of `spots` (places where the clustering put `g`'s
/// boundary on the chain, whose contact the settling of touched points
/// decides) are left out. Every point of it lies within
/// the tolerance of P (the chain does, L and each stretch are on P, and
/// the distance from P is linear along each segment). Where `g` holds a
/// point of it more than the tolerance inside its boundary (the point's
/// foot on P), the faces lie flush over a region that neither the
/// meeting's imprint nor classification separates, or cross there, and
/// the result would hold a fin thinner than the tolerance between them,
/// so the boolean is refused there.
///
/// The search runs in P ([`held_point`]), each point of `f`'s plane taken
/// to its foot there, which is linear and keeps straight edges straight.
/// It does not take curved edges: one of `g` that may come within the
/// search's reach of the strip, or one of `f` that may enter it, refuses
/// the boolean as not supported. A face with curved edges farther away
/// places its points by its chart's outline in `outlines` ([`place`]).
#[allow(clippy::too_many_arguments)]
fn held_strip(
    pair: &Pair,
    pts: &Points,
    outlines: &[Vec<Option<Outline>>; 2],
    strips: &mut StripSearch,
    dir: DVec3,
    [(op, f), (them, g)]: [(Operand, usize); 2],
    chain: &Chain,
    stretches: &[(usize, usize)],
    gaps: Option<&[(f64, f64)]>,
    spots: &[DVec3],
) -> Result<(), KernelError> {
    let eps = pts.eps;
    let rho = eps * (1.0 + DEEP);
    let (shape, theirs) = (pair.get(op), pair.get(them));
    let (face, gf) = (&shape.faces[f], &theirs.faces[g]);
    let (nf, ng) = (face.normal(), gf.normal());
    // In f's plane, `dir` along L and `across` square to it; P leans from
    // `across` by the angle between the planes, along `up`.
    let across = nf.cross(dir);
    let sin = across.dot(ng);
    if sin.abs() <= 1e-10 {
        return Ok(());
    }
    let lean = across - ng * sin;
    let cos = lean.length();
    // Faces this nearly square to each other have no strip: a run within
    // the tolerance of P lies at most ε / sin from L, so its strip is at
    // most ε (1 / sin − 1), under 10⁻¹² ε, wide.
    if cos <= 1e-6 {
        return Ok(());
    }
    let up = lean / cos;
    // The margin for rounding of the search, in P: ρ / 1024 in f's plane.
    let tau = rho * DEEP * cos;
    // The origin: on L, square across from the chain's first point.
    let a = chain.points[0];
    let a = a - nf * (a - face.frame.origin).dot(nf);
    let o = a - across * ((a - gf.frame.origin).dot(ng) / sin);
    // Coordinates in P along `dir` and `up`: a point of f's plane at its
    // foot, and a point of P itself; and back.
    let foot = |p: DVec3| DVec2::new((p - o).dot(dir), (p - o).dot(across) * cos);
    let in_p = |p: DVec3| DVec2::new((p - o).dot(dir), (p - o).dot(up));
    let on_f = |v: DVec2| o + dir * v.x + across * (v.y / cos);
    let on_g = |v: DVec2| o + dir * v.x + up * v.y;
    let (flat, budget) = (&strips.flat, &mut strips.budget);
    let feet: Vec<DVec2> = chain.points.iter().map(|&p| foot(p)).collect();
    // A point of f's plane `t` across L lies `t sin` off P along `ng`, so
    // the layer between it and P runs from it along −sign(t sin) `ng`, and
    // lies in f's body (behind f, against `nf`) where sign(t sin) nf·ng is
    // positive: on the side of L whose sign is `inner`.
    let inner = sin.signum() * nf.dot(ng).signum();
    // A flap, across L from the rest of f, is tested on its side where the
    // layer lies outside f's body (§21): on the other, the one operation
    // that would keep it and the piece of g beside it has them touch at
    // the corner from f's side, and is refused as touching at a point.
    let sides: &[f64] = match chain.kind {
        Kind::Run | Kind::Whole => &[1.0, -1.0],
        Kind::Inner | Kind::Flap if (inner > 0.0) == (chain.kind == Kind::Inner) => &[1.0],
        Kind::Inner | Kind::Flap => &[-1.0],
    };
    let lines: Vec<(DVec2, DVec2)> = (stretches.iter())
        .map(|&(c, d)| (foot(pts.rep[c]), foot(pts.rep[d])))
        .collect();
    let mut pieces = strip_of(&feet, chain.kind, sides, lines, eps, cos, budget)?;
    if let Some(gaps) = gaps {
        // A whole run's strip over its span's parts that no stretch
        // covers, their positions along L taken to the coordinates in P:
        // each piece meets the gaps found by a search, and is cut to each.
        let base = o.dot(dir);
        let mut kept = Vec::new();
        for piece in pieces {
            let (x0, x1) = piece.span();
            let first = gaps.partition_point(|g| g.1 - base <= x0);
            for &(u0, u1) in gaps[first..].iter().take_while(|g| g.0 - base < x1) {
                budget.spend(1)?;
                let part =
                    (piece.cut(DVec2::X, u0 - base)).and_then(|p| p.cut(-DVec2::X, base - u1));
                kept.extend(part);
            }
        }
        pieces = kept;
    }
    if pieces.is_empty() {
        return Ok(());
    }
    let strip = Strip::new(pieces, tau);
    let (lo, hi) = strip.bounds;
    let corners = [lo, DVec2::new(lo.x, hi.y), hi, DVec2::new(hi.x, lo.y)];
    let box_of = |map: &dyn Fn(DVec2) -> DVec3, pad: f64| {
        let (l, h) = (corners.iter()).fold(
            (DVec3::splat(f64::INFINITY), DVec3::splat(f64::NEG_INFINITY)),
            |(l, h), &v| (l.min(map(v)), h.max(map(v))),
        );
        (l - DVec3::splat(pad), h + DVec3::splat(pad))
    };
    // How far the strip's points are from P, at most (at its corners).
    let off = (strip.pieces.iter())
        .flat_map(|p| p.corners.iter())
        .fold(0.0f64, |m, v| m.max((v.y / cos * sin).abs()));
    let meets =
        |(l, h): (DVec3, DVec3), (bl, bh): (DVec3, DVec3)| l.cmple(bh).all() && bl.cmple(h).all();
    // The edges of `g` within ρ of the strip, in P (the search's reach):
    // among g's coedges within `off` + ρ of f's plane, which every point
    // within ρ of the strip's feet is, those whose boxes meet the strip's
    // box padded by ρ (each with a margin of ρ more, for rounding), found
    // by the face's tree over its edges' boxes.
    let reach = box_of(&on_g, 2.0 * rho);
    let mut coedges = Vec::new();
    gf.near_surface_in(
        &shape.surfaces[face.surface],
        off + 2.0 * rho,
        reach,
        &mut coedges,
    );
    budget.spend(coedges.len())?;
    let mut edges: Vec<usize> = (coedges.iter())
        .map(|&k| {
            let (li, i) = gf.coedge_at(k);
            gf.loops[li][i].0
        })
        .collect();
    edges.sort_unstable();
    edges.dedup();
    let mut near = Vec::new();
    for e in edges {
        if !meets(theirs.edge_boxes[e], reach) {
            continue;
        }
        if theirs.straight[e] {
            let (a, b) = theirs.edge_points(e);
            near.push((in_p(a), in_p(b)));
        } else if strip.reached(theirs, e, &in_p, rho + tau, budget)? {
            return Err(curved_beside(on_f(strip.pieces[0].corners[0])));
        }
    }
    // The edges of f that can enter the strip: among f's coedges within
    // `off` of P, those whose boxes meet the strip's box, found likewise.
    let within = box_of(&on_f, eps);
    face.near_surface_in(
        &theirs.surfaces[gf.surface],
        off + tau,
        within,
        &mut coedges,
    );
    budget.spend(coedges.len())?;
    let mut edges: Vec<usize> = (coedges.iter())
        .map(|&k| {
            let (li, i) = face.coedge_at(k);
            face.loops[li][i].0
        })
        .collect();
    edges.sort_unstable();
    edges.dedup();
    let mut cuts = Vec::new();
    for e in edges {
        if !meets(shape.edge_boxes[e], within) {
            continue;
        }
        if shape.straight[e] {
            let (a, b) = shape.edge_points(e);
            cuts.push((foot(a), foot(b)));
        } else if strip.reached(shape, e, &foot, -tau, budget)? {
            return Err(curved_beside(on_f(strip.pieces[0].corners[0])));
        }
    }
    // How a point is placed against each face: a flat face by its loops
    // in its plane (a step for each side, once, when its outline is first
    // built), a face with curved edges by its chart's outline, which
    // follows them (found by review: a flat face's outline takes a curved
    // edge by its chord, so a round cap held no point, and a round hole
    // every point its outer loop did).
    let mut against = |op: Operand, f: usize| -> Result<Against, KernelError> {
        match &outlines[slot(op)][f] {
            Some(o) => Ok(Against::Chart(o)),
            None => {
                let (o, fresh) = flat.get(pair, op, f);
                if fresh {
                    budget.spend(o.len())?;
                }
                Ok(Against::Flat(o))
            }
        }
    };
    let (against_f, against_g) = (against(op, f)?, against(them, g)?);
    // The points within the reach of a place where the clustering put g's
    // boundary on the chain are the contact's there: a cluster spreads at
    // most twice the tolerance about its place, so they lie within three
    // times it, measured in P as the search measures.
    let spots: Vec<DVec2> = spots.iter().map(|&p| in_p(p)).collect();
    let around = (&spots[..], 3.0 * eps);
    let found = held_point(&strip, &cuts, &near, around, rho, budget, |v, budget| {
        Ok(place(&against_g, gf, on_g(v), eps, budget)? == Some(true)
            && place(&against_f, face, on_f(v), eps, budget)? != Some(false))
    })?;
    match found {
        Some(v) => Err(KernelError::InvalidInput(format!(
            "{}, near {}",
            super::RUNS_WITHIN,
            show(on_f(v))
        ))),
        None => Ok(()),
    }
}

/// The outline [`place`] places a point against a face of a plane by.
enum Against<'s, 'a> {
    /// A flat face's loops in its plane ([`Face::to2`]).
    Flat(&'s FlatOutline),
    /// The chart's outline of a face with curved edges ([`outlines`]),
    /// which follows them.
    Chart(&'s Outline<'a>),
}

/// Where `p`, on the plane of face `f`, lies against the face by its
/// outline `against`: `Some(true)` inside, `Some(false)` outside, `None`
/// within `eps` of its boundary. A flat face's is decided in its plane
/// (within `eps` of a side, or else by the parity of a ray across its
/// sides), a face with curved edges as [`locate`] decides it. A step, and
/// one for each loop, side or chord looked at.
fn place(
    against: &Against,
    f: &Face,
    p: DVec3,
    eps: f64,
    budget: &mut Budget,
) -> Result<Option<bool>, KernelError> {
    let mut seen = 0;
    let at = match against {
        Against::Flat(outline) => {
            let q = f.to2(p);
            if outline.near(q, eps, &mut seen) {
                None
            } else {
                Some(outline.inside(q, &mut seen))
            }
        }
        Against::Chart(o) => o.locate_counting(o.chart.map(p), eps, eps / 16.0, &mut seen),
    };
    budget.spend(1 + seen)?;
    Ok(at)
}

/// The refusal of a strip ([`held_strip`]) that a curved edge may reach,
/// near `p`.
fn curved_beside(p: DVec3) -> KernelError {
    KernelError::Unsupported(format!(
        "booleans where a face of one body runs within the tolerance of the other's surface \
         beside a curved edge (design note 0012 §20) near {}",
        show(p)
    ))
}

/// The refusal of an edge of one body lying on a face of the other that
/// leaves the face where the tolerance cannot place the point, near `p`.
fn unplaced_leaving(p: DVec3) -> KernelError {
    KernelError::InvalidInput(format!("{} near {}", super::LEAVES_UNPLACED, show(p)))
}

/// A run of consecutive straight edges of a face's loop lying on a plane
/// of the other operand ([`lying_runs`]).
struct Run {
    /// Its edges, in loop order.
    edges: Vec<usize>,
    /// The places of its vertices, in loop order (one more than its edges).
    points: Vec<DVec3>,
    /// The clusters at its first vertex and at its last.
    ends: (usize, usize),
    /// True if the run is the whole interval of its face on the line, from
    /// the cluster at one end of it to the cluster at the other.
    exact: bool,
}

/// A face's runs lying on the other face's plane, for the checks of its
/// stretches ([`plane_meeting`]): each with its stretch along the planes'
/// common line, those stretches in order along the line, and the runs'
/// ends, so a stretch finds the runs beside it, and whether it is a run,
/// by searches.
struct Runs {
    runs: Vec<IndexedRun>,
    along: Spans,
    ends: Vec<(usize, usize)>,
}

/// A run with the clusters on its edges (sorted), found when a point
/// beside the run first needs them.
struct IndexedRun {
    run: Run,
    stretch: (f64, f64),
    on: OnceCell<Vec<usize>>,
}

impl IndexedRun {
    /// Whether cluster `c` lies on one of the run's edges (of `op`).
    fn holds(&self, c: usize, pts: &Points, op: Operand) -> bool {
        let on = self.on.get_or_init(|| {
            let mut on: Vec<usize> = (self.run.edges.iter())
                .flat_map(|&e| pts.on_edge(op, e).iter().copied())
                .collect();
            on.sort_unstable();
            on.dedup();
            on
        });
        on.binary_search(&c).is_ok()
    }
}

impl Runs {
    fn new(runs: Vec<Run>, dir: DVec3) -> Runs {
        let mut ends: Vec<(usize, usize)> = (runs.iter())
            .map(|r| (r.ends.0.min(r.ends.1), r.ends.0.max(r.ends.1)))
            .collect();
        ends.sort_unstable();
        let runs: Vec<IndexedRun> = (runs.into_iter())
            .map(|run| IndexedRun {
                stretch: run.stretch(dir),
                run,
                on: OnceCell::new(),
            })
            .collect();
        let along = Spans::new(runs.iter().map(|r| r.stretch));
        Runs { runs, along, ends }
    }

    /// Whether the interval `iv` is a run, from the cluster at one end of
    /// it to the cluster at the other: the whole interval of its face on
    /// the line.
    fn exact(&self, iv: &Interval) -> bool {
        let (c, d) = (iv.lo.1, iv.hi.1);
        self.ends.binary_search(&(c.min(d), c.max(d))).is_ok()
    }

    /// The runs whose stretches reach within 2ε of position `t` along the
    /// line (those within ε, and more for rounding), in order of their
    /// least positions.
    fn beside(&self, t: f64, eps: f64) -> impl Iterator<Item = &IndexedRun> {
        let (over, _) = self.along.meeting(t - 2.0 * eps, t + 2.0 * eps);
        over.map(|i| &self.runs[i])
    }
}

/// Spans along a line in order of their least positions, each with the
/// greatest end of those up to it, so the spans meeting a stretch of the
/// line are found by two searches and a pass over those starting before
/// its end and after the last span that ends before its start.
struct Spans {
    /// Each span's least and greatest position, and its place in the list
    /// given (a stable sort keeps that order among equal starts).
    along: Vec<(f64, f64, usize)>,
    /// The greatest end of the spans up to each in `along`.
    reach: Vec<f64>,
}

impl Spans {
    fn new(spans: impl Iterator<Item = (f64, f64)>) -> Spans {
        let mut along: Vec<(f64, f64, usize)> =
            spans.enumerate().map(|(i, (s, t))| (s, t, i)).collect();
        along.sort_by(|p, q| p.0.total_cmp(&q.0));
        let mut reach = Vec::with_capacity(along.len());
        for s in &along {
            let m = reach.last().copied().unwrap_or(f64::NEG_INFINITY);
            reach.push(s.1.max(m));
        }
        Spans { along, reach }
    }

    /// The places in the list given of the spans meeting `lo`..`hi`, in
    /// order of their least positions, and how many spans the pass over
    /// them takes.
    fn meeting(&self, lo: f64, hi: f64) -> (impl Iterator<Item = usize> + '_, usize) {
        let first = self.reach.partition_point(|&m| m < lo);
        let last = self.along.partition_point(|s| s.0 <= hi);
        let passed = self.along.get(first..last).unwrap_or(&[]);
        (
            (passed.iter()).filter(move |s| s.1 >= lo).map(|s| s.2),
            passed.len(),
        )
    }
}

impl Run {
    /// The least and greatest positions of its vertices along `dir`.
    fn stretch(&self, dir: DVec3) -> (f64, f64) {
        (self.points.iter()).fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), p| {
            let t = p.dot(dir);
            (lo.min(t), hi.max(t))
        })
    }

    /// The ends of its edge nearest `p`, in loop order (the first of
    /// those equally near).
    fn nearest(&self, p: DVec3) -> (DVec3, DVec3) {
        let ends = |k: usize| (self.points[k], self.points[k + 1]);
        let d = |k: usize| {
            let (a, b) = ends(k);
            points::segment_distance(p, a, b)
        };
        ends((0..self.edges.len()).fold(0, |best, k| if d(k) < d(best) { k } else { best }))
    }
}

/// The runs of face `f` of `op` lying on plane `surface` of the other
/// operand: consecutive straight edges of one of its loops whose vertices
/// are all within the tolerance of the plane, each marked exact where one
/// of the intervals `ivs` of the face on the planes' common line runs
/// along it from the cluster at one end to the cluster at the other. Only
/// the coedges near the plane are looked at.
fn lying_runs(
    pair: &Pair,
    pts: &Points,
    op: Operand,
    f: usize,
    surface: usize,
    ivs: &[Interval],
) -> Vec<Run> {
    let holds = |c: usize| {
        let (v0, v1) = pts.vertices_of(c);
        match op {
            Operand::First => v0.is_some(),
            Operand::Second => v1.is_some(),
        }
    };
    // The ends of the intervals whose two end clusters hold vertices of
    // the face's operand, each the lesser first, sorted, so a run finds
    // whether it is one by a search (found by review: a scan of them all
    // for each run took time growing with the square of the runs on a face
    // touching the line along many of them).
    let mut held: Vec<(usize, usize)> = (ivs.iter())
        .map(|iv| (iv.lo.1, iv.hi.1))
        .filter(|&(c, d)| holds(c) && holds(d))
        .map(|(c, d)| (c.min(d), c.max(d)))
        .collect();
    held.sort_unstable();
    let mut runs = Vec::new();
    let shape = pair.get(op);
    let face = &shape.faces[f];
    let lies = |li: usize, i: usize| {
        let (e, _) = face.loops[li][i];
        let vs = &face.verts[li];
        shape.straight[e]
            && pts.sign(op, vs[i], surface) == 0
            && pts.sign(op, vs[(i + 1) % vs.len()], surface) == 0
    };
    let mut near = Vec::new();
    face.near_surface(&pair.get(other(op)).surfaces[surface], pts.eps, &mut near);
    #[cfg(test)]
    super::work::add(&super::work::RUNS, near.len());
    for k in near {
        let (li, i) = face.coedge_at(k);
        let n = face.loops[li].len();
        // A run starts at a lying coedge whose predecessor does not lie (a
        // loop lying wholly on the plane has none: the faces would be
        // coplanar).
        if !lies(li, i) || lies(li, (i + n - 1) % n) {
            continue;
        }
        let vs = &face.verts[li];
        let mut edges = Vec::new();
        let mut points = vec![shape.points[vs[i]]];
        let mut j = i;
        while edges.len() < n && lies(li, j) {
            edges.push(face.loops[li][j].0);
            j = (j + 1) % n;
            points.push(shape.points[vs[j]]);
        }
        let (s, t) = (pts.vertex_cluster(op, vs[i]), pts.vertex_cluster(op, vs[j]));
        let exact = held.binary_search(&(s.min(t), s.max(t))).is_ok();
        runs.push(Run {
            edges,
            points,
            ends: (s, t),
            exact,
        });
    }
    runs
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
    // A line or a ray that misses the box where faces of both bodies can
    // be bounds no part of the result: its stretches in the face matter to
    // nothing. Beyond that box the tolerance need not decide where it
    // leaves the face (a plane within the tolerance of a cone along a
    // ruling there, the face's boundary grazing the cone: no rule refuses
    // it, as no result depends on it), so its cuts may not pair up (found
    // by the fuzz of second cuts: §17).
    if let Track::Line { origin, dir } | Track::Ray { origin, dir, .. } = *track {
        let (t0, t1) = track.ends();
        if !line_meets_box(origin, dir, t0, t1, region) {
            return Ok(Vec::new());
        }
    }
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

/// Whether the line `origin + t dir` for `t` in `[t0, t1]` meets the box
/// `(lo, hi)`.
fn line_meets_box(origin: DVec3, dir: DVec3, t0: f64, t1: f64, (lo, hi): (DVec3, DVec3)) -> bool {
    let (mut a, mut b) = (t0, t1);
    for k in 0..3 {
        if dir[k] == 0.0 {
            if origin[k] < lo[k] || origin[k] > hi[k] {
                return false;
            }
            continue;
        }
        let (x, y) = ((lo[k] - origin[k]) / dir[k], (hi[k] - origin[k]) / dir[k]);
        a = a.max(x.min(y));
        b = b.min(x.max(y));
    }
    a <= b
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

/// A convex piece of a strip ([`strip_pieces`]): the points `x` with
/// `n · x ≥ c` for each of its sides `(n, c)` (`n` a unit vector), and its
/// corners in order around it.
struct Piece {
    sides: Vec<(DVec2, f64)>,
    corners: Vec<DVec2>,
    /// The ends of the segment of the chain whose sweep it is a part of.
    along: (DVec2, DVec2),
}

impl Piece {
    /// The least and greatest x of its corners.
    fn span(&self) -> (f64, f64) {
        (self.corners.iter()).fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), v| {
            (lo.min(v.x), hi.max(v.x))
        })
    }

    /// Its part where `n · x ≥ c` (`n` a unit vector), if that has an area
    /// (Sutherland and Hodgman 1974). A side with the same normal is
    /// replaced by the tighter of the two, so a piece split along L again
    /// and again keeps its four sides and one for each meeting that cuts
    /// it (found by review: carried sides made each test of a piece beside
    /// n meetings take time in proportion to n).
    fn cut(&self, n: DVec2, c: f64) -> Option<Piece> {
        let k = self.corners.len();
        let mut corners: Vec<DVec2> = Vec::with_capacity(k + 1);
        for i in 0..k {
            let (u, v) = (self.corners[i], self.corners[(i + 1) % k]);
            let (du, dv) = (n.dot(u) - c, n.dot(v) - c);
            if du >= 0.0 {
                corners.push(u);
            }
            if (du >= 0.0) != (dv >= 0.0) {
                corners.push(u + (v - u) * (du / (du - dv)));
            }
        }
        corners.dedup();
        if corners.len() > 1 && corners[0] == corners[corners.len() - 1] {
            corners.pop();
        }
        let m = corners.len();
        let area: f64 = (0..m)
            .map(|i| corners[i].perp_dot(corners[(i + 1) % m]))
            .sum();
        if m < 3 || area == 0.0 {
            return None;
        }
        let mut sides = self.sides.clone();
        match sides.iter_mut().find(|side| side.0 == n) {
            Some(side) => side.1 = side.1.max(c),
            None => sides.push((n, c)),
        }
        Some(Piece {
            sides,
            corners,
            along: self.along,
        })
    }

    /// True if `v` is in the piece or within `tau` of it (outside each side
    /// by at most `tau`).
    fn contains(&self, v: DVec2, tau: f64) -> bool {
        self.sides.iter().all(|&(n, c)| n.dot(v) - c >= -tau)
    }

    /// The parameters of the segment `p`–`q` (0 at `p`, 1 at `q`) whose
    /// points are in the piece or within `tau` of it (Cyrus and Beck 1978),
    /// if they are not empty.
    fn clip(&self, p: DVec2, q: DVec2, tau: f64) -> Option<(f64, f64)> {
        let d = q - p;
        let (mut t0, mut t1) = (0.0f64, 1.0f64);
        for &(n, c) in &self.sides {
            // n · (p + t d) − c + tau ≥ 0.
            let (v, dv) = (n.dot(p) - c + tau, n.dot(d));
            if dv == 0.0 {
                if v < 0.0 {
                    return None;
                }
            } else if dv > 0.0 {
                t0 = t0.max(-v / dv);
            } else {
                t1 = t1.min(-v / dv);
            }
        }
        (t0 < t1).then_some((t0, t1))
    }
}

/// The pieces of the strip of a chain of `kind` whose points are at
/// `feet`, in the coordinates of P of [`held_strip`] (x along L, y across
/// it, each point of the face's plane at its foot on P), on the `sides` of
/// L listed (the signs of y): the sweep of each segment
/// ([`strip_pieces`]) or of a lone vertex ([`point_pieces`]) less its
/// points within `eps` of L, or a flap's piece ([`flap_piece`]), cut by
/// the meeting `lines`, the pair's stretches at their feet
/// ([`beyond_meeting`]; `cos` the cosine of the angle between the planes).
fn strip_of(
    feet: &[DVec2],
    kind: Kind,
    sides: &[f64],
    mut lines: Vec<(DVec2, DVec2)>,
    eps: f64,
    cos: f64,
    budget: &mut Budget,
) -> Result<Vec<Piece>, KernelError> {
    let (band, tau) = (eps * cos, eps * (1.0 + DEEP) * DEEP * cos);
    let mut pieces = Vec::new();
    match (kind, feet) {
        (Kind::Flap, &[x0, v, x1]) => {
            if sides.contains(&v.y.signum()) {
                pieces.extend(flap_piece(x0, v, x1, eps, cos));
            }
        }
        (_, &[v]) => point_pieces(v, band, tau, sides, &mut pieces),
        _ => {
            for w in feet.windows(2) {
                strip_pieces(w[0], w[1], band, tau, sides, &mut pieces);
            }
        }
    }
    // The stretches by their least position along L, each with the
    // greatest end of those up to it, to find those beside a piece.
    for l in &mut lines {
        if l.1.x < l.0.x {
            *l = (l.1, l.0);
        }
    }
    lines.sort_by(|p, q| p.0.x.total_cmp(&q.0.x));
    let reach: Vec<f64> = (lines.iter())
        .scan(f64::NEG_INFINITY, |m, l| {
            *m = m.max(l.1.x);
            Some(*m)
        })
        .collect();
    budget.spend(lines.len())?;
    let mut cut = Vec::new();
    for piece in pieces {
        let (x0, x1) = piece.span();
        let (lo, hi) = (
            reach.partition_point(|&m| m <= x0),
            lines.partition_point(|l| l.0.x < x1),
        );
        budget.spend(1 + hi.saturating_sub(lo))?;
        // The parts in order along L; those before `done` end at or before
        // the start of the line being cut, and so of every later one, which
        // keeps each whole (found by review: cutting every part by every
        // line took steps growing with the square of the meetings beside
        // one segment).
        let mut parts = vec![piece];
        let mut done = 0;
        for &line in lines.get(lo..hi).unwrap_or(&[]) {
            while done < parts.len() && parts[done].span().1 <= line.0.x {
                done += 1;
            }
            budget.spend(1 + parts.len() - done)?;
            let rest = parts.split_off(done);
            parts.extend((rest.into_iter()).flat_map(|p| beyond_meeting(p, line, eps, cos)));
        }
        cut.extend(parts);
    }
    Ok(cut)
}

/// The convex pieces of the sweep of the segment `a`–`b` to the x axis
/// (each point to its foot on it) that are at least `band` from the axis,
/// on each of the `sides` of it (the signs of y, +1 or −1) where the
/// segment reaches past `band` with its face on its left toward the axis,
/// into `out`. The segment runs as its face's loop
/// runs it, the face on its left (the x axis along L, the y axis across it
/// in the face's plane, the face's outward normal toward the viewer); on a
/// side where the face is on the segment's far side from the axis, the
/// sweep lies outside the face and bounds no strip of it. Each piece is a
/// trapezoid (or a triangle) between the segment's part beyond the band
/// and that part's foot on the band's edge. Its sides are the band's edge,
/// the two lines square to the axis through its ends, and the segment's
/// own line, so no side is taken from a short edge, whose direction
/// rounding could turn. A segment square to the axis to `tau` (its ends
/// no farther apart along it) sweeps no area to speak of and gives none.
fn strip_pieces(a: DVec2, b: DVec2, band: f64, tau: f64, sides: &[f64], out: &mut Vec<Piece>) {
    let e = b - a;
    if e.x.abs() <= tau {
        return;
    }
    for &side in sides {
        let (ya, yb) = (side * a.y, side * b.y);
        if (ya <= band && yb <= band) || side * e.x > 0.0 {
            continue;
        }
        let edge = side * band;
        let cut = |t: f64| DVec2::new(a.x + e.x * t, edge);
        let p = if ya < band {
            cut((band - ya) / (yb - ya))
        } else {
            a
        };
        let q = if yb < band {
            cut((band - ya) / (yb - ya))
        } else {
            b
        };
        let (x0, x1) = (p.x.min(q.x), p.x.max(q.x));
        if x0 >= x1 {
            continue;
        }
        // The side of the segment's line toward the axis.
        let n = e.perp() / e.length() * (-side * e.x.signum());
        let sides = vec![
            (DVec2::X, x0),
            (-DVec2::X, -x1),
            (DVec2::new(0.0, side), band),
            (n, n.dot(a)),
        ];
        let mut corners = vec![p, q, DVec2::new(q.x, edge), DVec2::new(p.x, edge)];
        corners.dedup();
        if corners.len() > 1 && corners[0] == corners[corners.len() - 1] {
            corners.pop();
        }
        out.push(Piece {
            sides,
            corners,
            along: (a, b),
        });
    }
}

/// The strip of a lone vertex at `v` ([`strip_runs`]): the segment from
/// it to its foot on the x axis, less its points within `band` of the
/// axis, as a piece `half` wide on each side of it, into `out` (none if
/// `v` is within `band` of the axis or on none of the `sides` of it, the
/// signs of y).
fn point_pieces(v: DVec2, band: f64, half: f64, sides: &[f64], out: &mut Vec<Piece>) {
    let side = v.y.signum();
    if v.y.abs() <= band || !sides.contains(&side) {
        return;
    }
    let (x0, x1, edge) = (v.x - half, v.x + half, side * band);
    out.push(Piece {
        sides: vec![
            (DVec2::X, x0),
            (-DVec2::X, -x1),
            (DVec2::new(0.0, side), band),
            (DVec2::new(0.0, -side), -side * v.y),
        ],
        corners: vec![
            DVec2::new(x0, v.y),
            DVec2::new(x1, v.y),
            DVec2::new(x1, edge),
            DVec2::new(x0, edge),
        ],
        along: (DVec2::new(x0, v.y), DVec2::new(x1, v.y)),
    });
}

/// The flap of a corner across the x axis at `v` whose edges cross the
/// axis at `x0` and `x1` ([`strip_runs`]), in the coordinates of
/// [`strip_of`]: the triangle between the edges and the axis, less its
/// points within `eps` of either edge in the face's plane (its coordinates
/// across the axis scaled by 1 / `cos`), a triangle with its base on the
/// axis or nothing. The axis carries no meeting there, so its points on
/// the axis are kept. A half-plane in the face's plane is one in P (the
/// map between them is linear), so each side is taken there and mapped.
fn flap_piece(x0: DVec2, v: DVec2, x1: DVec2, eps: f64, cos: f64) -> Option<Piece> {
    let side = v.y.signum();
    let (x0, x1) = (DVec2::new(x0.x, 0.0), DVec2::new(x1.x, 0.0));
    let in_face = |u: DVec2| DVec2::new(u.x, u.y / cos);
    let vf = in_face(v);
    // The half-plane `n · u ≥ c` of the face's plane (n a unit vector
    // there), in P's coordinates.
    let to_p = |n: DVec2, c: f64| {
        let m = DVec2::new(n.x, n.y / cos);
        let l = m.length();
        (m / l, c / l)
    };
    // The unit normal of the edge from `v` to `x` toward `other`, in the
    // face's plane.
    let inward = |x: DVec2, other: DVec2| {
        let e = in_face(x) - vf;
        let n = e.perp() / e.length();
        if n.dot(in_face(other) - vf) < 0.0 {
            -n
        } else {
            n
        }
    };
    let (n0, n1) = (inward(x0, x1), inward(x1, x0));
    let (lo, hi) = (x0.x.min(x1.x).min(v.x), x0.x.max(x1.x).max(v.x));
    if lo >= hi {
        return None;
    }
    let whole = Piece {
        sides: vec![
            (DVec2::new(0.0, side), 0.0),
            to_p(n0, n0.dot(vf)),
            to_p(n1, n1.dot(vf)),
        ],
        corners: vec![v, x0, x1],
        along: (DVec2::new(lo, v.y), DVec2::new(hi, v.y)),
    };
    let (m0, c0) = to_p(n0, n0.dot(vf) + eps);
    let (m1, c1) = to_p(n1, n1.dot(vf) + eps);
    whole.cut(m0, c0)?.cut(m1, c1)
}

/// The parts of a strip's piece ([`strip_pieces`]) left by the meeting
/// `a`–`b`, a stretch of the pair (each end at its foot on P, `a` the
/// nearer along L; `cos` the cosine of the angle between the planes),
/// over the stretch's span along L: those on the side of the stretch's
/// line where the piece's own chain is, farther than `eps` from the line
/// in the face's plane. Over a part of the span where the chain is beyond
/// the line from L, that is all of the piece but its points within `eps`
/// of the line; where the line is beyond L or beyond the chain, the
/// points of the piece within `eps` of it, near L or near the chain. The
/// piece is split where the chain crosses the line, each side kept on its
/// own; where the chain lies on the line nothing is kept. Outside the span
/// the piece is kept whole.
fn beyond_meeting(piece: Piece, (a, b): (DVec2, DVec2), eps: f64, cos: f64) -> Vec<Piece> {
    let (x0, x1) = piece.span();
    if b.x <= a.x || x1 <= a.x || b.x <= x0 {
        return vec![piece];
    }
    // The line in the face's plane (its coordinates across L scaled back by
    // 1 / cos), and the signed distance from it there, as a linear function
    // of the coordinates in P.
    let (af, bf) = (DVec2::new(a.x, a.y / cos), DVec2::new(b.x, b.y / cos));
    let nf = (bf - af).normalize().perp();
    let n2 = DVec2::new(nf.x, nf.y / cos);
    let g = |v: DVec2| n2.dot(v - a);
    let (p, q) = piece.along;
    let (gp, gq) = (g(p), g(q));
    let mut splits = vec![a.x, b.x];
    if (gp < 0.0) != (gq < 0.0) {
        splits.push(p.x + (q.x - p.x) * (gp / (gp - gq)));
    }
    splits.retain(|&x| x > x0 && x < x1);
    splits.sort_by(f64::total_cmp);
    splits.dedup();
    let mut parts = Vec::with_capacity(splits.len() + 1);
    let mut rest = Some(piece);
    for &x in &splits {
        let Some(r) = rest else { break };
        parts.extend(r.cut(-DVec2::X, -x));
        rest = r.cut(DVec2::X, x);
    }
    parts.extend(rest);
    let mut out = Vec::with_capacity(parts.len());
    for part in parts {
        let (u0, u1) = part.span();
        let x = (u0 + u1) * 0.5;
        if x <= a.x || x >= b.x {
            out.push(part);
            continue;
        }
        // The side of the line the chain is on over this part.
        let side = g(p + (q - p) * ((x - p.x) / (q.x - p.x)));
        if side == 0.0 {
            continue;
        }
        let m = n2 * side.signum();
        let l = m.length();
        out.extend(part.cut(m / l, (eps + m.dot(a)) / l));
    }
    out
}

/// A strip's pieces, with a tree over their boxes (padded by `tau`).
struct Strip {
    pieces: Vec<Piece>,
    tree: Bvh,
    tau: f64,
    /// The box of every piece.
    bounds: (DVec2, DVec2),
}

impl Strip {
    fn new(pieces: Vec<Piece>, tau: f64) -> Strip {
        let boxes: Vec<(DVec2, DVec2)> = (pieces.iter())
            .map(|p| {
                (p.corners.iter()).fold(
                    (DVec2::splat(f64::INFINITY), DVec2::splat(f64::NEG_INFINITY)),
                    |(lo, hi), &v| (lo.min(v), hi.max(v)),
                )
            })
            .collect();
        let bounds = (boxes.iter()).fold(
            (DVec2::splat(f64::INFINITY), DVec2::splat(f64::NEG_INFINITY)),
            |(lo, hi), &(l, h)| (lo.min(l), hi.max(h)),
        );
        let tree = Bvh::new(
            (boxes.iter())
                .map(|&(lo, hi)| flat(lo - DVec2::splat(tau), hi + DVec2::splat(tau)))
                .collect(),
        );
        Strip {
            pieces,
            tree,
            tau,
            bounds,
        }
    }

    /// The pieces whose (padded) boxes meet the box `lo`–`hi`, into `out`.
    fn meeting(&self, lo: DVec2, hi: DVec2, out: &mut Vec<usize>) {
        let (lo, hi) = flat(lo, hi);
        self.tree.overlapping(lo, hi, out);
    }

    /// True if `v` is in a piece, to `tau`.
    fn contains(&self, v: DVec2, budget: &mut Budget) -> Result<bool, KernelError> {
        let mut found = Vec::new();
        self.meeting(v, v, &mut found);
        budget.spend(1 + found.len())?;
        Ok(found.iter().any(|&k| self.pieces[k].contains(v, self.tau)))
    }

    /// True unless edge `e` of `shape`, curved, is farther than `gap` from
    /// every piece (a negative `gap`: enters none deeper than `-gap`), by
    /// the convex hulls of its rational Bézier pieces' control points,
    /// which hold it, taken into the strip's plane by `map`.
    fn reached(
        &self,
        shape: &Shape,
        e: usize,
        map: &dyn Fn(DVec3) -> DVec2,
        gap: f64,
        budget: &mut Budget,
    ) -> Result<bool, KernelError> {
        let (curve, t0, t1) = shape.curve(e);
        let Some(hulls) = quadric::hull_points(curve, t0, t1) else {
            return Ok(true);
        };
        let mut found = Vec::new();
        for hull in hulls {
            let hull: Vec<DVec2> = hull.into_iter().map(map).collect();
            let (lo, hi) = (hull.iter()).fold(
                (DVec2::splat(f64::INFINITY), DVec2::splat(f64::NEG_INFINITY)),
                |(lo, hi), &v| (lo.min(v), hi.max(v)),
            );
            let pad = DVec2::splat(gap.max(0.0));
            self.meeting(lo - pad, hi + pad, &mut found);
            budget.spend(1 + found.len())?;
            if (found.iter()).any(|&k| !apart(&hull, &self.pieces[k], gap)) {
                return Ok(true);
            }
        }
        Ok(false)
    }
}

/// A box of the plane as a box of space, at height zero.
fn flat(lo: DVec2, hi: DVec2) -> (DVec3, DVec3) {
    (lo.extend(0.0), hi.extend(0.0))
}

/// True if the convex hull of `points` and the piece are farther apart
/// than `gap` along some direction (a negative `gap`: overlap by less than
/// `-gap` along it): the separating axis test over the piece's sides, the
/// lines through two of the points, and the directions from each point to
/// each corner of the piece. Those hold the direction of least overlap of
/// two convex polygons and of their greatest separation, so the test is
/// exact up to rounding.
fn apart(points: &[DVec2], piece: &Piece, gap: f64) -> bool {
    fn span(vs: &[DVec2], n: DVec2) -> (f64, f64) {
        (vs.iter()).fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), v| {
            (lo.min(v.dot(n)), hi.max(v.dot(n)))
        })
    }
    let separates = |n: DVec2| {
        let ((a0, a1), (b0, b1)) = (span(points, n), span(&piece.corners, n));
        b0 - a1 > gap || a0 - b1 > gap
    };
    let pairs = (points.iter().enumerate())
        .flat_map(|(i, &p)| points[i + 1..].iter().map(move |&q| (q - p).perp()));
    let to_corners = (points.iter()).flat_map(|&p| piece.corners.iter().map(move |&c| c - p));
    (piece.sides.iter().map(|&(n, _)| n))
        .chain(pairs)
        .chain(to_corners)
        .filter_map(|n| n.try_normalize())
        .any(separates)
}

/// The parameters of the line `p + t (q − p)` whose points are within `r`
/// of the segment `a`–`b`, if any: the union of its meetings with the discs
/// about the segment's ends and with the band along it, which is one
/// interval since their union (a capsule) is convex.
fn capsule(p: DVec2, q: DVec2, (a, b): (DVec2, DVec2), r: f64) -> Option<(f64, f64)> {
    let d = q - p;
    let dd = d.length_squared();
    let (mut lo, mut hi) = (f64::INFINITY, f64::NEG_INFINITY);
    for c in [a, b] {
        let f = p - c;
        let half = f.dot(d);
        let disc = half * half - dd * (f.length_squared() - r * r);
        if disc >= 0.0 {
            let s = disc.sqrt();
            lo = lo.min((-half - s) / dd);
            hi = hi.max((-half + s) / dd);
        }
    }
    let e = b - a;
    let len = e.length();
    if len > 0.0 {
        let (u, n) = (e / len, e.perp() / len);
        let (mut t0, mut t1) = (f64::NEG_INFINITY, f64::INFINITY);
        let mut empty = false;
        // Each as v + t dv ≥ 0.
        for (v, dv) in [
            ((p - a).dot(u), d.dot(u)),
            ((b - p).dot(u), -d.dot(u)),
            (r - (p - a).dot(n), -d.dot(n)),
            (r + (p - a).dot(n), d.dot(n)),
        ] {
            if dv == 0.0 {
                empty |= v < 0.0;
            } else if dv > 0.0 {
                t0 = t0.max(-v / dv);
            } else {
                t1 = t1.min(-v / dv);
            }
        }
        if !empty && t0 <= t1 {
            lo = lo.min(t0);
            hi = hi.max(t1);
        }
    }
    (lo <= hi).then_some((lo, hi))
}

/// The parameters in (0, 1) of the segment `p`–`q` where the segment
/// `a`–`b` crosses or touches it, into `out`: a proper crossing, and each
/// end of `a`–`b` within `tau` of its line (so that a part of `p`–`q`
/// running along `a`–`b` is split from the rest at its ends).
fn crossings(p: DVec2, q: DVec2, (a, b): (DVec2, DVec2), tau: f64, out: &mut Vec<f64>) {
    let d = q - p;
    let dd = d.length_squared();
    let len = dd.sqrt();
    let mut put = |t: f64| {
        if t > 0.0 && t < 1.0 {
            out.push(t);
        }
    };
    for c in [a, b] {
        if d.perp_dot(c - p).abs() <= tau * len {
            put((c - p).dot(d) / dd);
        }
    }
    let e = b - a;
    let den = d.perp_dot(e);
    if den != 0.0 {
        let w = a - p;
        let (t, u) = (w.perp_dot(e) / den, w.perp_dot(d) / den);
        if (0.0..=1.0).contains(&u) {
            put(t);
        }
    }
}

/// Sorted intervals merged where they meet or overlap.
fn merge(mut list: Vec<(f64, f64)>) -> Vec<(f64, f64)> {
    list.sort_by(|x, y| x.0.total_cmp(&y.0).then(x.1.total_cmp(&y.1)));
    let mut out: Vec<(f64, f64)> = Vec::with_capacity(list.len());
    for (a, b) in list {
        match out.last_mut() {
            Some(last) if a <= last.1 => last.1 = last.1.max(b),
            _ => out.push((a, b)),
        }
    }
    out
}

/// A point of a strip ([`Strip`]) that `holds` accepts, found among the
/// points where a region of the strip held deeper than `rho` must reach
/// (design note 0012 §20, *Finding a held point*), in the plane of the
/// other face: `cuts` are the edges of the strip's face that can enter it
/// (taken there), `near` the edges of the other face within `rho` of it
/// (straight), and `around` places with a reach, whose discs are left out
/// of the strip (§21).
///
/// Let D be the points of the other face farther than `rho` from its
/// boundary and farther than the reach from every place, and K a part of
/// D in the strip's interior, apart from the rest. K's boundary runs along
/// the strip's boundary (the pieces' sides and the cuts) or along D's,
/// which is made of the lines at `rho` from the near edges, arcs of radius
/// `rho` about their ends and arcs of the reach about the places; so it
/// has a stretch on one of those lines, or it is bounded by arcs alone and
/// has a corner where two of the circles cross or touch. The search tests
/// (1) on each of those lines, the parts in the strip farther than `rho`
/// from every near edge (a capsule about each, [`capsule`]) and farther
/// than the reach from every place, split where a cut crosses or touches
/// them, each at its middle: no near edge comes within `rho` of a part and
/// no cut crosses it, so each face holds all of it or none of it; and (2)
/// for each two of the circles that meet (and a margin for rounding), the
/// two crossings (one point where they touch) that no place's disc holds
/// (to the margin): a corner of K is one, and lies in K's closure, at
/// least `rho` from every edge of the other face. So it finds a point of
/// every such K (to the strip's `tau`).
fn held_point(
    strip: &Strip,
    cuts: &[(DVec2, DVec2)],
    near: &[(DVec2, DVec2)],
    (spots, reach): (&[DVec2], f64),
    rho: f64,
    budget: &mut Budget,
    mut holds: impl FnMut(DVec2, &mut Budget) -> Result<bool, KernelError>,
) -> Result<Option<DVec2>, KernelError> {
    let tau = strip.tau;
    let segment_box = |(a, b): (DVec2, DVec2), pad: f64| {
        flat(a.min(b) - DVec2::splat(pad), a.max(b) + DVec2::splat(pad))
    };
    let near_tree = Bvh::new(near.iter().map(|&s| segment_box(s, rho)).collect());
    let cut_tree = Bvh::new(cuts.iter().map(|&s| segment_box(s, tau)).collect());
    let spot_tree = Bvh::new(spots.iter().map(|&v| segment_box((v, v), reach)).collect());
    budget.spend(strip.pieces.len() + near.len() + cuts.len() + spots.len())?;
    // (1) The lines: the pieces' sides, the cuts, and the lines at `rho`
    // on each side of each near edge, with the edge they come from.
    let mut lines: Vec<(DVec2, DVec2, Option<usize>)> = Vec::new();
    for piece in &strip.pieces {
        let c = &piece.corners;
        lines.extend((0..c.len()).map(|k| (c[k], c[(k + 1) % c.len()], None)));
    }
    lines.extend(cuts.iter().map(|&(a, b)| (a, b, None)));
    for (k, &(a, b)) in near.iter().enumerate() {
        let n = (b - a).perp().normalize_or_zero() * rho;
        if n != DVec2::ZERO {
            lines.push((a + n, b + n, Some(k)));
            lines.push((a - n, b - n, Some(k)));
        }
    }
    let mut found = Vec::new();
    let mut splits = Vec::new();
    for &(p, q, own) in &lines {
        let d = q - p;
        if d == DVec2::ZERO {
            continue;
        }
        let (lo, hi) = (p.min(q), p.max(q));
        strip.meeting(lo, hi, &mut found);
        budget.spend(1 + found.len())?;
        let inside = merge(
            (found.iter())
                .filter_map(|&k| strip.pieces[k].clip(p, q, tau))
                .collect(),
        );
        if inside.is_empty() {
            continue;
        }
        let (lo3, hi3) = flat(lo, hi);
        near_tree.overlapping(lo3, hi3, &mut found);
        budget.spend(found.len())?;
        let mut out: Vec<(f64, f64)> = (found.iter())
            .filter(|&&k| Some(k) != own)
            .filter_map(|&k| capsule(p, q, near[k], rho))
            .collect();
        spot_tree.overlapping(lo3, hi3, &mut found);
        budget.spend(found.len())?;
        out.extend((found.iter()).filter_map(|&k| capsule(p, q, (spots[k], spots[k]), reach)));
        let out = merge(out);
        cut_tree.overlapping(lo3, hi3, &mut found);
        budget.spend(found.len())?;
        splits.clear();
        for &k in &found {
            crossings(p, q, cuts[k], tau, &mut splits);
        }
        splits.sort_by(f64::total_cmp);
        // The parts of `inside` clear of `out`, split at `splits`.
        let mut x = 0;
        for &(a, b) in &inside {
            let mut from = a;
            while x < out.len() && out[x].1 <= a {
                x += 1;
            }
            let mut ends = Vec::new();
            for &(c, e) in out[x..].iter().take_while(|iv| iv.0 < b) {
                if c > from {
                    ends.push((from, c));
                }
                from = from.max(e);
            }
            ends.push((from, b));
            for (s, t) in ends {
                let first = splits.partition_point(|&u| u <= s);
                let mut from = s;
                for u in (splits[first..].iter().copied().take_while(|&u| u < t)).chain([t]) {
                    if u > from {
                        let m = p + d * ((from + u) * 0.5);
                        if holds(m, budget)? {
                            return Ok(Some(m));
                        }
                    }
                    from = u;
                }
            }
        }
    }
    // (2) The crossings of the circles about the ends of near edges, of
    // radius `rho`, and about the places, of the reach.
    let mut ends: Vec<DVec2> = near.iter().flat_map(|&(a, b)| [a, b]).collect();
    ends.sort_by(|u, v| u.x.total_cmp(&v.x).then(u.y.total_cmp(&v.y)));
    ends.dedup();
    let circles: Vec<(DVec2, f64)> = (ends.iter().map(|&v| (v, rho)))
        .chain(spots.iter().map(|&v| (v, reach)))
        .collect();
    let circle_tree = Bvh::new(
        (circles.iter())
            .map(|&(v, r)| segment_box((v, v), r))
            .collect(),
    );
    let mut inside = Vec::new();
    for (i, &(v, r)) in circles.iter().enumerate() {
        // Boxes padded by each radius meet where the centres are within
        // the two radii in each coordinate.
        let (lo, hi) = segment_box((v, v), r);
        circle_tree.overlapping(lo, hi, &mut found);
        budget.spend(1 + found.len())?;
        for &j in found.iter().filter(|&&j| j > i) {
            let (w, s) = circles[j];
            let l = (w - v).length();
            // Circles that touch bound a corner too (a cusp); a margin
            // keeps rounding from skipping them.
            if l == 0.0 || l > (r + s) * (1.0 + DEEP) || l < (r - s).abs() {
                continue;
            }
            let (e, u) = ((w - v) / l, (w - v).perp() / l);
            let a = (r * r - s * s + l * l) / (2.0 * l);
            let h = (r * r - a * a).max(0.0).sqrt();
            for side in [1.0, -1.0] {
                let m = v + e * a + u * (side * h);
                // A crossing in a place's disc is no corner of D.
                let (lo, hi) = flat(m, m);
                spot_tree.overlapping(lo, hi, &mut inside);
                budget.spend(inside.len())?;
                let clear =
                    (inside.iter()).all(|&k| (m - spots[k]).length() >= reach * (1.0 - DEEP));
                if clear && strip.contains(m, budget)? && holds(m, budget)? {
                    return Ok(Some(m));
                }
            }
        }
    }
    Ok(None)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    /// A deterministic pseudo-random sequence in [0, 1).
    fn lcg(state: &mut u64) -> f64 {
        *state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (*state >> 11) as f64 / (1u64 << 53) as f64
    }

    /// A polygon star-shaped about `c`, `n` corners at radii in
    /// `r0`..`r1`.
    fn star(s: &mut u64, c: DVec2, n: usize, r0: f64, r1: f64) -> Vec<DVec2> {
        (0..n)
            .map(|k| {
                let a = 2.0 * PI * (k as f64 + 0.8 * lcg(s)) / n as f64;
                c + DVec2::new(a.cos(), a.sin()) * (r0 + (r1 - r0) * lcg(s))
            })
            .collect()
    }

    fn edges(loops: &[Vec<DVec2>]) -> Vec<(DVec2, DVec2)> {
        (loops.iter())
            .flat_map(|l| (0..l.len()).map(move |k| (l[k], l[(k + 1) % l.len()])))
            .collect()
    }

    fn distance(v: DVec2, edges: &[(DVec2, DVec2)]) -> f64 {
        (edges.iter()).fold(f64::INFINITY, |m, &(a, b)| {
            let t = ((v - a).dot(b - a) / (b - a).length_squared()).clamp(0.0, 1.0);
            m.min((v - (a + (b - a) * t)).length())
        })
    }

    /// Inside by the parity of the edges crossed by a ray to +x.
    fn inside(v: DVec2, edges: &[(DVec2, DVec2)]) -> bool {
        (edges.iter()).fold(false, |odd, &(a, b)| {
            let cross =
                (a.y > v.y) != (b.y > v.y) && a.x + (v.y - a.y) / (b.y - a.y) * (b.x - a.x) > v.x;
            odd != cross
        })
    }

    /// The strip search finds a point wherever a grid finds one of the
    /// strip, in the strip's face, more than `rho` inside the other face
    /// (and some margin, for the grid's rounding), and every point it finds
    /// passes the test it was given: on random strips (bent runs, some
    /// crossing the line) against random faces with holes, most with discs
    /// about a few places left out.
    #[test]
    fn the_strip_search_finds_every_deeply_held_region() {
        let (eps, mut s) = (1.0, 11);
        let rho = eps * (1.0 + DEEP);
        // Per kind of face, the cases held deeper than `rho` and the cases
        // where the search finds a point.
        let (mut held, mut found) = ([0; 2], [0; 2]);
        for case in 0..200 {
            let mut run = vec![DVec2::new(0.0, 1.0 + 4.0 * lcg(&mut s))];
            for _ in 0..1 + case % 4 {
                let last = *run.last().unwrap();
                run.push(DVec2::new(
                    last.x + 1.0 + 4.0 * lcg(&mut s),
                    6.0 * lcg(&mut s) - (case % 5) as f64,
                ));
            }
            // Half the runs go the other way: their face lies above them,
            // so that they sweep the side above the line.
            if (case / 2) % 2 == 1 {
                run.reverse();
            }
            let mut pieces = Vec::new();
            for w in run.windows(2) {
                strip_pieces(w[0], w[1], eps, rho * DEEP, &[1.0, -1.0], &mut pieces);
            }
            if pieces.is_empty() {
                continue;
            }
            let strip = Strip::new(pieces, rho * DEEP);
            let (lo, hi) = strip.bounds;
            let mid = (lo + hi) * 0.5;
            let size = (hi - lo).max_element();
            // Half the faces are a star with a few holes; the others are
            // riddled with holes about `rho` apart, leaving small pockets
            // deeper than `rho` between them, some bounded by arcs alone.
            let mut g = vec![star(&mut s, mid, 5 + case % 7, 0.2 * size, 0.9 * size)];
            let holes = if case % 2 == 0 {
                case % 4
            } else {
                ((hi - lo).x * (hi - lo).y / 2.5) as usize
            };
            for _ in 0..holes {
                let c = lo + (hi - lo) * DVec2::new(lcg(&mut s), lcg(&mut s));
                g.push(star(&mut s, c, 3 + case % 3, 0.3, 1.5));
            }
            // The strip's face: a star around the strip, or (every third
            // case) a small one inside it, with holes.
            let mut f = if case % 3 == 0 {
                let c = lo + (hi - lo) * DVec2::new(lcg(&mut s), lcg(&mut s));
                vec![star(&mut s, c, 4, 0.3, 1.5)]
            } else {
                vec![star(&mut s, mid, 6, 0.4 * size, 1.2 * size)]
            };
            for _ in 0..case % 3 {
                let c = lo + (hi - lo) * DVec2::new(lcg(&mut s), lcg(&mut s));
                f.push(star(&mut s, c, 3, 0.3, 1.5));
            }
            let (ge, fe) = (edges(&g), edges(&f));
            // Two cases in three leave out the discs about a few places
            // (§21), on the strip or near it, of a reach from `rho` to four
            // times it.
            let spots: Vec<DVec2> = (0..(case % 3) * (1 + case % 4))
                .map(|_| lo + (hi - lo) * DVec2::new(lcg(&mut s), lcg(&mut s)))
                .collect();
            let reach = rho * (1.0 + 3.0 * lcg(&mut s));
            let clear = |v: DVec2, r: f64| spots.iter().all(|&c| (v - c).length() > r);
            let holds = |v: DVec2| {
                inside(v, &ge)
                    && distance(v, &ge) > eps
                    && (inside(v, &fe) || distance(v, &fe) <= eps)
            };
            let mut budget = Budget::new(STRIP_STEPS);
            let around = (&spots[..], reach);
            let got = held_point(&strip, &fe, &ge, around, rho, &mut budget, |v, _| {
                Ok(holds(v))
            })
            .unwrap();
            if let Some(v) = got {
                assert!(holds(v), "case {case}: {v} fails the test");
                assert!(
                    clear(v, reach * (1.0 - 2.0 * DEEP)),
                    "case {case}: {v} in a disc"
                );
                assert!(strip.contains(v, &mut Budget::new(STRIP_STEPS)).unwrap());
                found[case % 2] += 1;
            }
            let n = 160;
            let deep = (0..=n * n).any(|k| {
                let v = lo + (hi - lo) * DVec2::new((k % n) as f64, (k / n) as f64) / n as f64;
                (strip.pieces.iter()).any(|p| p.contains(v, -1e-9))
                    && inside(v, &fe)
                    && distance(v, &fe) > 1e-9
                    && inside(v, &ge)
                    && distance(v, &ge) > rho + 1e-9
                    && clear(v, reach + 1e-9)
            });
            if deep {
                held[case % 2] += 1;
                assert!(got.is_some(), "case {case}: a held point was missed");
            }
        }
        // The cases are varied: some strips are held, some are not.
        assert!(held.iter().all(|&h| (10..90).contains(&h)), "{held:?}");
    }

    /// A small face inside the strip, all of it deep inside the other
    /// face, is found on its own edges (no line of the strip or of the
    /// other face passes through it).
    #[test]
    fn the_strip_search_finds_a_small_face_inside_the_strip() {
        let eps = 1.0;
        let rho = eps * (1.0 + DEEP);
        let mut pieces = Vec::new();
        strip_pieces(
            DVec2::new(10.0, 10.0),
            DVec2::new(0.0, 10.0),
            eps,
            rho * DEEP,
            &[1.0, -1.0],
            &mut pieces,
        );
        let strip = Strip::new(pieces, rho * DEEP);
        let ge = edges(&[vec![
            DVec2::new(-20.0, -20.0),
            DVec2::new(30.0, -20.0),
            DVec2::new(30.0, 30.0),
            DVec2::new(-20.0, 30.0),
        ]]);
        let fe = edges(&[vec![
            DVec2::new(4.0, 4.0),
            DVec2::new(5.0, 4.0),
            DVec2::new(4.5, 5.0),
        ]]);
        let holds = |v: DVec2| {
            inside(v, &ge) && distance(v, &ge) > eps && (inside(v, &fe) || distance(v, &fe) <= eps)
        };
        let mut budget = Budget::new(STRIP_STEPS);
        let got = held_point(&strip, &fe, &ge, (&[], 0.0), rho, &mut budget, |v, _| {
            Ok(holds(v))
        })
        .unwrap();
        let v = got.expect("the small face is found");
        assert!(distance(v, &fe) < 1e-9, "{v}");
    }

    /// A region of the other face deeper than `rho` bounded by arcs alone
    /// (three holes whose tips are closer than `2 rho` in pairs, around a
    /// point farther than `rho` from each) is found at a crossing of two
    /// of the tips' circles, where no line the search clips passes through
    /// it.
    #[test]
    fn the_strip_search_finds_a_region_bounded_by_arcs() {
        let eps = 1.0;
        let rho = eps * (1.0 + DEEP);
        let r = 1.1 * rho;
        let centre = DVec2::new(5.0, 3.0);
        let mut g = vec![vec![
            DVec2::new(-20.0, -20.0),
            DVec2::new(30.0, -20.0),
            DVec2::new(30.0, 30.0),
            DVec2::new(-20.0, 30.0),
        ]];
        for k in 0..3 {
            let axis = PI / 2.0 + 2.0 * PI * k as f64 / 3.0;
            let dir = |a: f64| DVec2::new(a.cos(), a.sin());
            let tip = centre + dir(axis) * r;
            // A wedge opening away from the centre, its sides at 60° from
            // its axis, clockwise as a hole.
            g.push(vec![
                tip,
                tip + dir(axis + PI / 3.0) * 12.0,
                tip + dir(axis - PI / 3.0) * 12.0,
            ]);
        }
        let mut pieces = Vec::new();
        strip_pieces(
            DVec2::new(10.0, 10.0),
            DVec2::new(0.0, 10.0),
            eps,
            rho * DEEP,
            &[1.0, -1.0],
            &mut pieces,
        );
        let strip = Strip::new(pieces, rho * DEEP);
        let ge = edges(&g);
        let f = [vec![
            DVec2::new(3.0, 1.5),
            DVec2::new(7.0, 1.5),
            DVec2::new(7.0, 4.5),
            DVec2::new(3.0, 4.5),
        ]];
        let fe = edges(&f);
        let holds = |v: DVec2| inside(v, &ge) && distance(v, &ge) > eps && inside(v, &fe);
        assert!(inside(centre, &ge) && distance(centre, &ge) > rho);
        let mut budget = Budget::new(STRIP_STEPS);
        let got = held_point(&strip, &fe, &ge, (&[], 0.0), rho, &mut budget, |v, _| {
            Ok(holds(v))
        })
        .unwrap();
        let v = got.expect("the region bounded by arcs is found");
        assert!((v - centre).length() < 0.3 * rho, "{v}");
    }

    /// The region of [`the_strip_search_finds_a_region_bounded_by_arcs`]
    /// with a place's disc (§21) over each of its corners, so that what is
    /// left of it is bounded by arcs whose corners are crossings of a
    /// place's circle with a tip's: it is found at one of those, clear of
    /// the discs; and nothing is found when the discs grow to cover it.
    #[test]
    fn the_strip_search_finds_a_region_bounded_by_the_arcs_of_places() {
        let eps = 1.0;
        let rho = eps * (1.0 + DEEP);
        let r = 1.1 * rho;
        let centre = DVec2::new(5.0, 3.0);
        let dir = |a: f64| DVec2::new(a.cos(), a.sin());
        let mut g = vec![vec![
            DVec2::new(-20.0, -20.0),
            DVec2::new(30.0, -20.0),
            DVec2::new(30.0, 30.0),
            DVec2::new(-20.0, 30.0),
        ]];
        let tips = [0, 1, 2].map(|k| PI / 2.0 + 2.0 * PI * k as f64 / 3.0);
        for axis in tips {
            let tip = centre + dir(axis) * r;
            g.push(vec![
                tip,
                tip + dir(axis + PI / 3.0) * 12.0,
                tip + dir(axis - PI / 3.0) * 12.0,
            ]);
        }
        // The region's corners: the crossings of the tips' circles nearest
        // the centre.
        let corners = [0, 1, 2].map(|k| {
            let (u, w) = (
                centre + dir(tips[k]) * r,
                centre + dir(tips[(k + 1) % 3]) * r,
            );
            let m = (u + w) * 0.5;
            let half = (w - u).length() * 0.5;
            centre
                + (m - centre).normalize()
                    * ((m - centre).length() - (rho * rho - half * half).sqrt())
        });
        let mut pieces = Vec::new();
        strip_pieces(
            DVec2::new(10.0, 10.0),
            DVec2::new(0.0, 10.0),
            eps,
            rho * DEEP,
            &[1.0, -1.0],
            &mut pieces,
        );
        let strip = Strip::new(pieces, rho * DEEP);
        let ge = edges(&g);
        let fe = edges(&[vec![
            DVec2::new(3.0, 1.5),
            DVec2::new(7.0, 1.5),
            DVec2::new(7.0, 4.5),
            DVec2::new(3.0, 4.5),
        ]]);
        let holds = |v: DVec2| inside(v, &ge) && distance(v, &ge) > eps && inside(v, &fe);
        let mut budget = Budget::new(STRIP_STEPS);
        let reach = 0.1 * rho;
        let got = held_point(
            &strip,
            &fe,
            &ge,
            (&corners, reach),
            rho,
            &mut budget,
            |v, _| Ok(holds(v)),
        )
        .unwrap();
        let v = got.expect("the region bounded by the places' arcs is found");
        assert!((v - centre).length() < 0.3 * rho, "{v}");
        for c in corners {
            assert!(
                (v - c).length() >= reach * (1.0 - 2.0 * DEEP),
                "{v} in a disc"
            );
        }
        let got = held_point(
            &strip,
            &fe,
            &ge,
            (&corners, 0.3 * rho),
            rho,
            &mut budget,
            |v, _| Ok(holds(v)),
        )
        .unwrap();
        assert_eq!(got, None);
    }

    /// The strip of `feet` on `sides` with the meeting `lines`
    /// ([`strip_of`]), at ε = 1 and with the planes' cosine 1.
    fn strip(feet: &[DVec2], sides: &[f64], lines: &[(DVec2, DVec2)]) -> Vec<Piece> {
        let mut budget = Budget::new(STRIP_STEPS);
        let kind = if feet.len() == 1 {
            Kind::Inner
        } else {
            Kind::Run
        };
        strip_of(feet, kind, sides, lines.to_vec(), 1.0, 1.0, &mut budget).unwrap()
    }

    /// Whether a piece of `pieces` holds `v` (by more than rounding).
    fn holds(pieces: &[Piece], v: DVec2) -> bool {
        pieces.iter().any(|p| p.contains(v, -1e-9))
    }

    /// A strip leaves out its points within ε of the line, and covers only
    /// the sides of it asked for, a segment's on the side where its face
    /// lies toward the line, a vertex's on its own side.
    #[test]
    fn a_strip_leaves_out_the_band_and_keeps_to_its_sides() {
        let v = DVec2::new;
        let both = [1.0, -1.0];
        // The face below the segment, toward the line.
        let above = [v(10.0, 5.0), v(0.0, 5.0)];
        let s = strip(&above, &both, &[]);
        assert!(holds(&s, v(5.0, 1.5)) && holds(&s, v(5.0, 4.5)));
        assert!(!holds(&s, v(5.0, 0.5)) && !holds(&s, v(5.0, -2.0)));
        assert!(strip(&above, &[-1.0], &[]).is_empty());
        assert!(holds(&strip(&above, &[1.0], &[]), v(5.0, 1.5)));
        // The face above the segment, away from the line: no strip.
        assert!(strip(&[v(0.0, 5.0), v(10.0, 5.0)], &both, &[]).is_empty());
        let below = strip(&[v(0.0, -5.0), v(10.0, -5.0)], &both, &[]);
        assert!(holds(&below, v(5.0, -1.5)) && !holds(&below, v(5.0, -0.5)));
        // A vertex.
        let apex = [v(5.0, 5.0)];
        let s = strip(&apex, &both, &[]);
        assert!(holds(&s, v(5.0, 3.0)) && !holds(&s, v(5.0, 0.5)));
        assert!(strip(&apex, &[-1.0], &[]).is_empty());
        assert!(holds(&strip(&[v(5.0, -5.0)], &[-1.0], &[]), v(5.0, -3.0)));
        assert!(strip(&[v(5.0, 0.5)], &both, &[]).is_empty());
    }

    /// A flap's piece ([`flap_piece`]) is the triangle between the corner's
    /// edges and L on the corner's side, less its points within ε of
    /// either edge, its points on L kept, in P's coordinates (the face's
    /// across L scaled by the planes' cosine); and nothing where no point
    /// of the triangle is farther than ε from both edges.
    #[test]
    fn a_flap_is_its_triangle_less_the_band_of_its_edges() {
        let v = DVec2::new;
        // To rounding: L bounds the piece, and its points on L are in it.
        let holds = |p: &Option<Piece>, u: DVec2| p.as_ref().is_some_and(|p| p.contains(u, 1e-12));
        // A right angle 4 from L, its bisector square to L: points to 2.83
        // from both edges.
        let flap = flap_piece(v(-4.0, 0.0), v(0.0, 4.0), v(4.0, 0.0), 1.0, 1.0);
        assert!(holds(&flap, v(0.0, 1.5)) && holds(&flap, v(0.0, 0.0)));
        assert!(holds(&flap, v(1.0, 0.0)) && holds(&flap, v(-1.0, 0.5)));
        // Within ε of an edge, across L, or beyond the corner.
        for u in [
            v(0.0, 3.0),
            v(-2.5, 0.5),
            v(2.5, 0.5),
            v(0.0, -0.5),
            v(0.0, 4.5),
        ] {
            assert!(!holds(&flap, u), "{u}");
        }
        // The corner 1.3 from L: no point farther than 0.92 from both.
        assert!(flap_piece(v(-1.3, 0.0), v(0.0, 1.3), v(1.3, 0.0), 1.0, 1.0).is_none());
        // Below L, and with the bisector slanting (the edges to x = −1 and
        // 12 from the corner at (0, 4)).
        let below = flap_piece(v(-4.0, 0.0), v(0.0, -4.0), v(4.0, 0.0), 1.0, 1.0);
        assert!(holds(&below, v(0.0, -1.5)) && !holds(&below, v(0.0, 1.5)));
        let slant = flap_piece(v(-1.0, 0.0), v(0.0, 4.0), v(12.0, 0.0), 1.0, 1.0);
        assert!(holds(&slant, v(1.5, 1.0)) && !holds(&slant, v(0.0, 1.0)));
        // A face at 60° to P: the right angle 4 from L in the face is 2
        // from it in P, and its points 1.5 and 3 from L in the face are 0.75
        // and 1.5 from it in P.
        let steep = flap_piece(v(-4.0, 0.0), v(0.0, 2.0), v(4.0, 0.0), 1.0, 0.5);
        assert!(holds(&steep, v(0.0, 0.75)) && !holds(&steep, v(0.0, 1.5)));
    }

    /// The parts of a whole run's span that no stretch covers, each
    /// stretch grown by the tolerance: none where one stretch or several
    /// overlapping ones (in any order) cover it, the gaps between and
    /// beyond stretches otherwise, and the whole span with no stretch.
    #[test]
    fn a_whole_run_is_tested_where_no_stretch_covers_it() {
        // Binary fractions, so the sums are exact.
        let eps = 0.0625;
        let gaps = |spans: &[(f64, f64)]| uncovered((0.0, 1.0), spans.iter().copied(), eps);
        assert_eq!(gaps(&[]), vec![(0.0, 1.0)]);
        assert_eq!(gaps(&[(0.0, 1.0)]), vec![]);
        assert_eq!(gaps(&[(0.0625, 0.9375)]), vec![]);
        assert_eq!(gaps(&[(0.5, 1.0), (0.0, 0.5)]), vec![]);
        assert_eq!(
            gaps(&[(0.625, 0.96875), (0.0, 0.25)]),
            vec![(0.3125, 0.5625)]
        );
        assert_eq!(gaps(&[(0.5, 1.0)]), vec![(0.0, 0.4375)]);
        assert_eq!(
            gaps(&[(0.3125, 0.5), (0.25, 0.375)]),
            vec![(0.0, 0.1875), (0.5625, 1.0)]
        );
        assert_eq!(gaps(&[(-0.5, -0.25), (1.25, 1.5)]), vec![(0.0, 1.0)]);
    }

    /// The chains and notch tips ([`strip_runs`]) of the bottom face of a
    /// prism, its profile `pts` (w in tolerances where marked) on the
    /// plane through the x axis tilted 10° toward +y, against the top
    /// z = 0 of a box: each chain's kind and points, and each tip's place.
    fn corners(pts: &[(f64, f64, bool)]) -> (Vec<(Kind, Vec<DVec3>)>, Vec<DVec3>) {
        use crate::math::Frame;
        use crate::{BrepKernel, ExtrudeSpec, Kernel, Profile};
        let kernel = BrepKernel::default();
        let eps = kernel.tolerance().linear;
        let t = 10f64.to_radians();
        let tilted = Frame {
            origin: DVec3::ZERO,
            x: DVec3::X,
            y: DVec3::new(0.0, t.cos(), t.sin()),
            z: DVec3::new(0.0, -t.sin(), t.cos()),
        };
        let profile: Vec<DVec2> = (pts.iter())
            .map(|&(u, w, tol)| DVec2::new(u, if tol { w * eps } else { w }))
            .collect();
        let spec = ExtrudeSpec {
            distance: 0.5,
            symmetric: false,
        };
        let prism = (kernel.extrude(&Profile::polygon(tilted, &profile), &spec)).expect("a prism");
        let block = kernel
            .make_box(
                &Frame {
                    origin: DVec3::new(-0.25, -1.0, -1.0),
                    ..Frame::WORLD
                },
                DVec3::new(1.5, 2.0, 1.0),
            )
            .expect("a box");
        let pair = Pair {
            first: super::super::shape::shape(&block.brep, "first").expect("a shape"),
            second: super::super::shape::shape(&prism.brep, "second").expect("a shape"),
        };
        let pts = Points::new(&pair, eps);
        let facing = |s: &Shape, n: DVec3| {
            (s.faces.iter())
                .position(|f| f.normal().dot(n) > 1.0 - 1e-12)
                .expect("a face")
        };
        let (g, f) = (
            facing(&pair.first, DVec3::Z),
            facing(&pair.second, -tilted.z),
        );
        let surface = pair.first.faces[g].surface;
        let (chains, tips) = strip_runs(&pair, &pts, Operand::Second, f, surface, &[]);
        (
            chains.into_iter().map(|c| (c.kind, c.points)).collect(),
            tips.into_iter().map(|v| pair.second.points[v]).collect(),
        )
    }

    /// [`strip_runs`] sorts each corner within the tolerance of the plane
    /// once (design note 0012 §21): a reflex corner across L from its
    /// neighbours is a notch tip; a convex one, a flap from its edges'
    /// crossings of the plane through it; a reflex corner whose edges run
    /// away from the plane on its side, an apex; and nothing else is either
    /// (a convex corner pointing at L, a corner on L, a corner whose edges
    /// leave it to both sides of the plane).
    #[test]
    fn the_corners_within_the_tolerance_of_the_plane_are_sorted_once() {
        let eps = crate::Tolerance::default().linear;
        let t = 10f64.to_radians();
        let at = |u: f64, w: f64| DVec3::new(u, w * eps * t.cos(), w * eps * t.sin());
        let notch = |w: f64| {
            corners(&[
                (-0.5, -0.5, false),
                (2.0, -0.5, false),
                (2.0, 0.5, false),
                (0.6, 0.5, false),
                (0.5, w, true),
                (0.4, 0.5, false),
                (-0.5, 0.5, false),
            ])
        };
        let spike = |w: f64| {
            corners(&[
                (-0.5, -0.5, false),
                (2.0, -0.5, false),
                (2.0, -0.2, false),
                (0.6, -0.2, false),
                (0.5, w, true),
                (0.4, -0.2, false),
                (-0.5, -0.2, false),
            ])
        };
        let close = |p: DVec3, q: DVec3| (p - q).length() < 1e-15;
        // Across: a notch tip, and a flap whose ends are on the plane.
        let (chains, tips) = notch(-3.0);
        assert!(chains.is_empty() && tips.len() == 1 && close(tips[0], at(0.5, -3.0)));
        let (chains, tips) = spike(3.0);
        assert!(tips.is_empty() && chains.len() == 1, "{chains:?}");
        let (kind, points) = &chains[0];
        assert_eq!(*kind, Kind::Flap);
        assert!(points.len() == 3 && close(points[1], at(0.5, 3.0)));
        // The edges cross the plane 3 / 0.2 tolerances' worth of 0.1 along
        // x from the corner, one to either side.
        let reach = 0.1 * 3.0 * eps / (3.0 * eps + 0.2);
        let mut ends = [points[0], points[2]];
        ends.sort_by(|p, q| p.x.total_cmp(&q.x));
        assert!(ends.iter().all(|p| p.z.abs() < 1e-18), "{ends:?}");
        assert!((ends[0].x - (0.5 - reach)).abs() < 1e-15, "{ends:?}");
        assert!((ends[1].x - (0.5 + reach)).abs() < 1e-15, "{ends:?}");
        // An apex: the notch's tip short of L.
        let (chains, tips) = notch(3.0);
        assert!(tips.is_empty() && chains.len() == 1);
        assert_eq!(chains[0].0, Kind::Inner);
        assert!(chains[0].1.len() == 1 && close(chains[0].1[0], at(0.5, 3.0)));
        // Neither: a spike pointing at L short of it, a tip on L, and a
        // corner whose edges leave it to both sides.
        let pointing = corners(&[
            (-0.5, 0.2, false),
            (0.4, 0.2, false),
            (0.5, 3.0, true),
            (0.6, 0.2, false),
            (2.0, 0.2, false),
            (2.0, 0.5, false),
            (-0.5, 0.5, false),
        ]);
        let bent = corners(&[
            (-0.5, -0.5, false),
            (2.0, -0.5, false),
            (2.0, 0.5, false),
            (-0.5, 0.5, false),
            (-0.4, -3.0, true),
        ]);
        for (chains, tips) in [pointing, notch(0.0), spike(0.0), bent] {
            assert!(chains.is_empty() && tips.is_empty(), "{chains:?} {tips:?}");
        }
    }

    /// A strip ends ε short of a meeting between its chain and the line,
    /// over the meeting's span along the line; is split where its chain
    /// crosses the meeting's line, each side kept on its own side of it;
    /// and has nothing where its chain lies on that line.
    #[test]
    fn a_strip_ends_short_of_a_meeting() {
        let v = DVec2::new;
        let both = [1.0, -1.0];
        let chain = [v(10.0, 5.0), v(0.0, 5.0)];
        for line in [(v(-1.0, 3.0), v(11.0, 3.0)), (v(11.0, 3.0), v(-1.0, 3.0))] {
            let s = strip(&chain, &both, &[line]);
            assert!(holds(&s, v(5.0, 4.5)));
            assert!(!holds(&s, v(5.0, 3.5)) && !holds(&s, v(5.0, 2.0)));
        }
        // A meeting over part of the span.
        let s = strip(&chain, &both, &[(v(4.0, 3.0), v(6.0, 3.0))]);
        assert!(holds(&s, v(2.0, 2.0)) && holds(&s, v(8.0, 2.0)));
        assert!(!holds(&s, v(5.0, 2.0)) && holds(&s, v(5.0, 4.5)));
        // The chain crossing the meeting's line at x = 20 / 3: beyond it
        // to the left, between it and the line to the right.
        let line = (v(-1.0, 3.0), v(11.0, 3.0));
        let s = strip(&[v(10.0, 1.5), v(0.0, 6.0)], &both, &[line]);
        assert!(holds(&s, v(2.0, 5.0)) && holds(&s, v(8.0, 1.5)));
        assert!(!holds(&s, v(2.0, 3.5)) && !holds(&s, v(8.0, 2.2)));
        // The chain on the meeting's line, on either side of L (the face on
        // its left, toward L).
        for (y, x) in [(3.0, 10.0), (-3.0, 0.0)] {
            let line = (v(-1.0, y), v(11.0, y));
            let chain = [v(x, y), v(10.0 - x, y)];
            assert!(!strip(&chain, &both, &[]).is_empty());
            assert!(strip(&chain, &both, &[line]).is_empty());
        }
    }

    /// A segment beside `n` meetings, each over its own stretch of L
    /// between the segment and L, slanting a little: cutting its strip
    /// takes steps in proportion to `n`, and each piece keeps at most five
    /// sides (its four and its meeting's), at 256 and 1,024 meetings (found
    /// by review: every part was cut by every meeting, about n² steps, and
    /// each part carried a side for every split before it).
    #[test]
    fn a_strip_beside_many_meetings_takes_steps_in_proportion_to_them() {
        let v = DVec2::new;
        for n in [256, 1024] {
            let w = 10.0 / n as f64;
            let lines: Vec<(DVec2, DVec2)> = (0..n)
                .map(|i| (v(w * i as f64, 3.0), v(w * (i as f64 + 0.5), 3.0 + 0.1 * w)))
                .collect();
            let mut budget = Budget::new(STRIP_STEPS);
            let pieces = strip_of(
                &[v(10.0, 5.0), v(0.0, 5.0)],
                Kind::Run,
                &[1.0, -1.0],
                lines,
                1.0,
                1.0,
                &mut budget,
            )
            .unwrap();
            assert_eq!(pieces.len(), 2 * n);
            let steps = STRIP_STEPS - budget.left;
            assert!(steps <= 8 * n, "{n}: {steps} steps");
            let most = pieces.iter().map(|p| p.sides.len()).max();
            assert_eq!(most, Some(5), "{n}");
        }
    }

    /// [`apart`] against the separation of the points' hull and a piece
    /// along 7,200 directions, on random points and pieces: where some
    /// direction separates them by more than `gap` (a negative `gap`: they
    /// overlap along it by less than `-gap`), `apart` says they are apart,
    /// and where it says so, the best direction sampled comes within the
    /// sampling's error of `gap`.
    #[test]
    fn apart_is_the_separation_of_a_hull_and_a_piece() {
        let mut s = 11;
        let dirs: Vec<DVec2> = (0..7200)
            .map(|i| DVec2::from_angle(i as f64 * TAU / 7200.0))
            .collect();
        let span = |vs: &[DVec2], n: DVec2| {
            (vs.iter()).fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), v| {
                (lo.min(v.dot(n)), hi.max(v.dot(n)))
            })
        };
        let (mut said, mut not) = (0, 0);
        for case in 0..600 {
            let mut at = || DVec2::new(6.0 * lcg(&mut s) - 3.0, 2.0 * lcg(&mut s) + 1.0);
            let (a, b) = (at(), at());
            let (a, b) = if a.x > b.x { (a, b) } else { (b, a) };
            let mut out = Vec::new();
            strip_pieces(a, b, 1.0, 1e-9, &[1.0], &mut out);
            let Some(piece) = out.pop() else {
                continue;
            };
            let points: Vec<DVec2> = (0..2 + case % 4)
                .map(|_| DVec2::new(6.0 * lcg(&mut s) - 3.0, 6.0 * lcg(&mut s) - 2.0))
                .collect();
            let gap = 2.0 * lcg(&mut s) - 1.0;
            let best = (dirs.iter())
                .map(|&n| {
                    let ((a0, a1), (b0, b1)) = (span(&points, n), span(&piece.corners, n));
                    (b0 - a1).max(a0 - b1)
                })
                .fold(f64::NEG_INFINITY, f64::max);
            let got = apart(&points, &piece, gap);
            if best > gap + 1e-12 {
                assert!(got, "case {case}: {best} against {gap}");
            }
            if got {
                assert!(best > gap - 0.01, "case {case}: {best} against {gap}");
                said += 1;
            } else {
                not += 1;
            }
        }
        assert!(said > 100 && not > 100, "{said} apart, {not} not");
    }

    /// The strip check of a run of a prism's face lying on the plane of a
    /// block's top (the tolerance 10⁻³; the planes 0.25° apart, the run 0.2
    /// off the line where they cross, the strip between), the face having a
    /// hole inside the strip and the block's top lying in the strip: a
    /// diamond hole holding the middles of the sides of the top's region
    /// deeper than ρ, and a rectangle holding all of that region.
    /// Where the corners of that region are left in the face, held and
    /// bounded only by the hole's edges, the check finds them through the
    /// face's own edges that enter the strip (design note 0012 §20,
    /// *Finding a held point*): each line it tests elsewhere has its middle
    /// outside the top or in the hole.
    #[test]
    fn the_strip_check_takes_its_faces_edges_that_enter_the_strip() {
        use crate::math::Frame;
        use crate::{BrepKernel, ExtrudeSpec, Kernel, Profile, Tolerance};
        let tol = Tolerance {
            linear: 1e-3,
            ..Tolerance::default()
        };
        let kernel = BrepKernel::with_tolerance(tol);
        let t = 0.25f64.to_radians();
        let tilted = Frame {
            origin: DVec3::ZERO,
            x: DVec3::X,
            y: DVec3::new(0.0, t.cos(), t.sin()),
            z: DVec3::new(0.0, -t.sin(), t.cos()),
        };
        let v = |pts: &[(f64, f64)]| -> Vec<DVec2> {
            pts.iter().map(|&(x, y)| DVec2::new(x, y)).collect()
        };
        let outer = v(&[
            (-0.5, -0.5),
            (0.0, -0.2),
            (1.0, -0.2),
            (1.5, -0.5),
            (2.0, -0.5),
            (2.0, 0.5),
            (-0.5, 0.5),
        ]);
        let block = kernel
            .make_box(
                &Frame {
                    origin: DVec3::new(0.3, -0.17, -1.0),
                    ..Frame::WORLD
                },
                DVec3::new(0.4, 0.14, 1.0),
            )
            .expect("a block");
        let check = |hole: &[(f64, f64)]| {
            let prism = kernel
                .extrude(
                    &Profile::polygon(tilted, &outer).with(Profile::polygon(tilted, &v(hole))),
                    &ExtrudeSpec {
                        distance: 0.5,
                        symmetric: false,
                    },
                )
                .expect("a prism");
            let pair = Pair {
                first: super::super::shape::shape(&block.brep, "first").expect("a shape"),
                second: super::super::shape::shape(&prism.brep, "second").expect("a shape"),
            };
            let pts = Points::new(&pair, tol.linear);
            let outlines = [
                outlines(&pair.first).expect("outlines"),
                outlines(&pair.second).expect("outlines"),
            ];
            let facing = |s: &Shape, n: DVec3| {
                (s.faces.iter())
                    .position(|f| f.normal().dot(n) > 1.0 - 1e-12)
                    .expect("a face")
            };
            let (g, f) = (
                facing(&pair.first, DVec3::Z),
                facing(&pair.second, -tilted.z),
            );
            let dir = (pair.first.faces[g].normal())
                .cross(pair.second.faces[f].normal())
                .normalize();
            let chain = Chain {
                points: [1.0, 0.0].map(|x| tilted.x * x - tilted.y * 0.2).to_vec(),
                kind: Kind::Run,
                edges: Vec::new(),
                corner: None,
            };
            held_strip(
                &pair,
                &pts,
                &outlines,
                &mut StripSearch::new(STRIP_STEPS),
                dir,
                [(Operand::Second, f), (Operand::First, g)],
                &chain,
                &[],
                None,
                &[],
            )
        };
        let diamond = [(0.5, -0.175), (0.72, -0.1), (0.5, -0.025), (0.28, -0.1)];
        match check(&diamond) {
            Err(KernelError::InvalidInput(m)) => assert!(m.starts_with(super::super::RUNS_WITHIN)),
            other => panic!("{other:?}"),
        }
        let rectangle = [(0.29, -0.18), (0.71, -0.18), (0.71, -0.02), (0.29, -0.02)];
        assert!(check(&rectangle).is_ok());
    }

    /// Cutting each piece by only the meetings over its span gives the same
    /// pieces as cutting it by all of them, in the same order, on random
    /// chains and meetings; and the meetings beside no piece cost nothing
    /// beyond their sorting.
    #[test]
    fn a_strip_takes_only_the_meetings_over_each_piece() {
        let mut s = 7;
        let both = [1.0, -1.0];
        let tau = (1.0 + DEEP) * DEEP;
        let mut cut = 0;
        for case in 0..400 {
            let mut feet = vec![DVec2::new(
                12.0 * lcg(&mut s) - 1.0,
                12.0 * lcg(&mut s) - 6.0,
            )];
            for _ in 0..1 + case % 4 {
                feet.push(DVec2::new(
                    12.0 * lcg(&mut s) - 1.0,
                    12.0 * lcg(&mut s) - 6.0,
                ));
            }
            let lines: Vec<(DVec2, DVec2)> = (0..case % 30)
                .map(|_| {
                    let a = DVec2::new(14.0 * lcg(&mut s) - 2.0, 12.0 * lcg(&mut s) - 6.0);
                    let b = a + DVec2::new(0.1 + 5.0 * lcg(&mut s), 4.0 * lcg(&mut s) - 2.0);
                    if lcg(&mut s) < 0.5 {
                        (a, b)
                    } else {
                        (b, a)
                    }
                })
                .collect();
            let got = strip(&feet, &both, &lines);
            let mut all = Vec::new();
            for w in feet.windows(2) {
                strip_pieces(w[0], w[1], 1.0, tau, &both, &mut all);
            }
            let mut sorted: Vec<(DVec2, DVec2)> = (lines.iter())
                .map(|&(a, b)| if a.x <= b.x { (a, b) } else { (b, a) })
                .collect();
            sorted.sort_by(|p, q| p.0.x.total_cmp(&q.0.x));
            let before = all.len();
            for line in sorted {
                all = (all.into_iter())
                    .flat_map(|p| beyond_meeting(p, line, 1.0, 1.0))
                    .collect();
            }
            cut += usize::from(all.len() != before);
            assert_eq!(got.len(), all.len(), "case {case}");
            for (p, q) in got.iter().zip(&all) {
                assert_eq!(
                    (&p.corners, &p.sides),
                    (&q.corners, &q.sides),
                    "case {case}"
                );
            }
        }
        assert!(cut > 100, "{cut}");
        let v = DVec2::new;
        let lines: Vec<(DVec2, DVec2)> = (0..1000)
            .map(|i| {
                let x = if i % 2 == 0 { -30.0 } else { 20.0 } + 0.01 * i as f64;
                (v(x, 3.0), v(x + 1.0, 3.0))
            })
            .collect();
        let mut budget = Budget::new(STRIP_STEPS);
        let pieces = strip_of(
            &[v(10.0, 5.0), v(0.0, 5.0)],
            Kind::Run,
            &both,
            lines.clone(),
            1.0,
            1.0,
            &mut budget,
        )
        .unwrap();
        assert_eq!(pieces.len(), 1);
        assert_eq!(STRIP_STEPS - budget.left, 1000 + 1);
        // One step fewer left, and the boolean is refused as too large.
        let refused = strip_of(
            &[v(10.0, 5.0), v(0.0, 5.0)],
            Kind::Run,
            &both,
            lines,
            1.0,
            1.0,
            &mut Budget::new(1000),
        );
        assert!(
            matches!(&refused, Err(KernelError::TooLarge(m)) if m.contains(STRIP_WORK)),
            "{:?}",
            refused.map(|p| p.len())
        );
    }
}
