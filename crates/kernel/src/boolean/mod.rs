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
use crate::geom::Curve3;
use crate::origin::{Operand, Origins};
use crate::topo::Brep;
use crate::{BooleanOp, KernelError, Tolerance};

mod assemble;
mod bernstein;
mod chart;
mod classify;
#[cfg(test)]
mod dd;
mod imprint;
mod meet;
mod planar;
mod points;
mod quadric;
mod shape;
#[cfg(test)]
pub(crate) mod testing;

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
    /// Along curved track `id`, from parameter `ta` at `a` to `tb` at `b`
    /// (less than a turn apart on a closed track).
    Track { id: usize, ta: f64, tb: f64 },
}

/// The courses of curved segments: the operands' curved edges and the
/// tracks where their surfaces meet.
#[derive(Clone, Copy)]
struct Courses<'a> {
    pair: &'a Pair<'a>,
    tracks: &'a [Track],
}

impl Courses<'_> {
    /// The point of `course` at parameter `t`.
    fn point(&self, course: Course, t: f64) -> DVec3 {
        match course {
            Course::Straight => unreachable!("a straight segment has no parameter"),
            Course::Edge { op, e, .. } => self.pair.get(op).curve(e).0.point(t),
            Course::Track { id, .. } => self.tracks[id].point(t),
        }
    }

    /// The point of `course` at the fraction `f` of its parameters' way
    /// from `ta` to `tb`.
    fn along(&self, course: Course, f: f64) -> DVec3 {
        match course {
            Course::Straight => unreachable!("a straight segment has no parameter"),
            Course::Edge { ta, tb, .. } | Course::Track { ta, tb, .. } => {
                self.point(course, ta + (tb - ta) * f)
            }
        }
    }

    /// The point of `course` as far from one end as from the other: at
    /// least half its chord from each, wherever its parameter crowds (the
    /// point halfway along a hyperbola's NURBS parameter can be within
    /// 10⁻⁴ of the stretch of an end, review measured). By bisection on
    /// the difference of the squares of the distances, which is affine in
    /// the point, so changes sign once along a conic's stretch.
    fn middle(&self, course: Course) -> DVec3 {
        let (ta, tb) = match course {
            Course::Straight => unreachable!("a straight segment has no parameter"),
            Course::Edge { ta, tb, .. } | Course::Track { ta, tb, .. } => (ta, tb),
        };
        let (a, b) = (self.point(course, ta), self.point(course, tb));
        let side = |t: f64| {
            let x = self.point(course, t);
            (x - a).length_squared() - (x - b).length_squared()
        };
        let (mut lo, mut hi) = (ta, tb);
        for _ in 0..60 {
            let m = 0.5 * (lo + hi);
            if m == lo || m == hi {
                break;
            }
            if side(m) < 0.0 {
                lo = m;
            } else {
                hi = m;
            }
        }
        self.point(course, 0.5 * (lo + hi))
    }

    /// The point of `course` where its tangent has turned the fraction `f`
    /// of its whole turning from `ta` to `tb`, for an ellipse, a conic or
    /// a NURBS curve (`None` for a circle, whose parameter is its angle):
    /// on a conic track in closed form ([`meet::conic_at_turn`]), else by
    /// [`turned`], an ellipse's stretch broken into pieces of at most a
    /// quarter of its parameter's period and a NURBS curve's at its knots
    /// ([`ellipse_breaks`], [`knot_breaks`]).
    fn at_turn(&self, course: Course, f: f64) -> Option<DVec3> {
        let t = match course {
            Course::Straight => unreachable!("a straight segment has no parameter"),
            Course::Track { id, ta, tb } => match self.tracks[id] {
                Track::Conic { c, .. } => {
                    let (fa, fb) = (meet::conic_turn(c, ta), meet::conic_turn(c, tb));
                    meet::conic_at_turn(c, fa + (fb - fa) * f)
                }
                track @ Track::Ellipse { .. } => {
                    turned(|t| track.deriv(t), &ellipse_breaks(ta, tb), f)
                }
                _ => return None,
            },
            Course::Edge { op, e, ta, tb } => {
                let curve = self.pair.get(op).curve(e).0;
                let Curve3::Nurbs(n) = curve else {
                    return None;
                };
                turned(|t| curve.deriv(t), &knot_breaks(n, ta, tb), f)
            }
        };
        Some(self.point(course, t))
    }

    /// The distance from `p` (a point near the curve) to the stretch of
    /// `course` between its parameters: to the foot of `p` on the curve
    /// where that is in the stretch, else to the nearer end.
    fn distance(&self, course: Course, p: DVec3) -> f64 {
        match course {
            Course::Straight => unreachable!("a straight segment has no parameter"),
            Course::Edge { op, e, ta, tb } => {
                points::foot_on(self.pair.get(op).curve(e).0, ta, tb, p).1
            }
            Course::Track { id, ta, tb } => {
                let track = self.tracks[id];
                let (lo, hi) = if ta <= tb { (ta, tb) } else { (tb, ta) };
                let (t, d) = track.nearest(p);
                let t = if track.closed() {
                    lo + (t - lo).rem_euclid(TAU)
                } else {
                    t
                };
                if (lo..=hi).contains(&t) {
                    d
                } else {
                    (track.point(lo) - p)
                        .length()
                        .min((track.point(hi) - p).length())
                }
            }
        }
    }
}

/// The parameter between the first and the last of `breaks` (in order
/// along the stretch) where the direction of `deriv` has turned the
/// fraction `f` of its whole turning along them. The breaks are chosen
/// so that between consecutive ones a plane curve that turns one way (a
/// conic's stretch) turns less than a half turn ([`ellipse_breaks`],
/// [`knot_breaks`]), so its turning there is the angle between its
/// directions at the two, and grows steadily from the first: the turnings
/// are summed to find the piece, then the parameter is found in it by
/// bisection on the angle from its first direction. Along another curve,
/// some parameter between the ends. (The angle from the stretch's first
/// direction alone falls back after a half turn, which crowded the probes
/// of an ellipse's longer stretches into their first part; found by
/// review.)
fn turned(deriv: impl Fn(f64) -> DVec3, breaks: &[f64], f: f64) -> f64 {
    let angle = |d0: DVec3, d: DVec3| d0.cross(d).length().atan2(d0.dot(d));
    let dirs: Vec<DVec3> = breaks.iter().map(|&t| deriv(t)).collect();
    let turns: Vec<f64> = dirs.windows(2).map(|w| angle(w[0], w[1])).collect();
    let mut goal = f * turns.iter().sum::<f64>();
    let mut k = 0;
    while k + 1 < turns.len() && goal > turns[k] {
        goal -= turns[k];
        k += 1;
    }
    let (d0, mut lo, mut hi) = (dirs[k], breaks[k], breaks[k + 1]);
    // Within 10⁻¹² of the piece's parameters: far finer than a probe
    // needs.
    for _ in 0..40 {
        let m = 0.5 * (lo + hi);
        if angle(d0, deriv(m)) < goal {
            lo = m;
        } else {
            hi = m;
        }
    }
    0.5 * (lo + hi)
}

/// Breaks of an ellipse's stretch from `ta` to `tb` for [`turned`]: equal
/// pieces of at most a quarter of its parameter's period. The tangents at
/// parameters half a period apart are opposite (the ellipse is symmetric
/// about its centre) and turn steadily, so a piece shorter than half a
/// period turns less than a half turn.
fn ellipse_breaks(ta: f64, tb: f64) -> Vec<f64> {
    let k = ((tb - ta).abs() / std::f64::consts::FRAC_PI_2)
        .ceil()
        .max(1.0) as usize;
    (0..=k)
        .map(|i| {
            if i == k {
                tb
            } else {
                ta + (tb - ta) * (i as f64 / k as f64)
            }
        })
        .collect()
}

/// Breaks of the stretch of the NURBS curve `n` from `ta` to `tb` for
/// [`turned`]: its ends, the knots between them, and the halvings of any
/// piece between them not shown to turn less than a half turn
/// ([`within_half_turn`]). A rational quadratic piece of a conic, the
/// edges booleans write, is always shown so; a piece of higher degree
/// can turn further (review: an arc of 238° of an ellipse as one cubic
/// span, its weights positive), and is halved until its parts are.
fn knot_breaks(n: &opendatum_nurbs::NurbsCurve, ta: f64, tb: f64) -> Vec<f64> {
    /// Halvings of a piece before it is kept as it is: a smooth curve's
    /// pieces pass long before (its control points close in on its
    /// tangent), a cusp's never do.
    const DEPTH: usize = 30;
    let (lo, hi) = (ta.min(tb), ta.max(tb));
    let segments = n.bezier_segments();
    // The pieces' ends from `lo` up, each piece halved until it passes.
    let mut breaks = vec![lo];
    let mut stack: Vec<(f64, f64, usize)> = quadric::spans(&segments, lo, hi)
        .map(|(a, b, _)| (a, b, 0))
        .collect();
    stack.reverse();
    while let Some((a, b, depth)) = stack.pop() {
        let m = 0.5 * (a + b);
        let pass = depth >= DEPTH
            || m <= a
            || m >= b
            || quadric::spans(&segments, a, b).all(|(_, _, pw)| within_half_turn(&pw));
        if pass {
            breaks.push(b);
        } else {
            stack.push((m, b, depth + 1));
            stack.push((a, m, depth + 1));
        }
    }
    if breaks.len() == 1 {
        breaks.push(hi);
    }
    // The stretch's own ends, exactly.
    *breaks.last_mut().expect("two breaks") = hi;
    if ta > tb {
        breaks.reverse();
    }
    breaks
}

/// Whether every tangent of the rational Bézier piece with homogeneous
/// control points `pw` (positive weights) lies on one side of a plane
/// through the origin, so that the piece, turning one way, turns less
/// than a half turn. Its tangent is a combination with nonnegative
/// coefficients of the differences `Pⱼ − Pᵢ` (`i < j`) of its control
/// points: by the quotient rule (Piegl & Tiller, *The NURBS Book*, ch.
/// 4) `X' = Σᵢ<ⱼ wᵢ wⱼ (Bᵢ Bⱼ' − Bᵢ' Bⱼ) (Pⱼ − Pᵢ) / W²`, and `Bₖ' / Bₖ =
/// k / s − (n − k) / (1 − s)` grows with `k`. So
/// the test asks whether each difference has a positive component along
/// the sum of the directions of the first and the last leg: true of a
/// quadratic whose legs are not opposite, as a conic's arc of less than
/// a half turn has. A sufficient test, not a necessary one.
fn within_half_turn(pw: &[glam::DVec4]) -> bool {
    let p: Vec<DVec3> = pw.iter().map(|q| q.truncate() / q.w).collect();
    let n = p.len() - 1;
    let first = (1..=n).map(|j| p[j] - p[0]).find(|d| *d != DVec3::ZERO);
    let last = (0..n)
        .rev()
        .map(|i| p[n] - p[i])
        .find(|d| *d != DVec3::ZERO);
    let (Some(first), Some(last)) = (first, last) else {
        // A point.
        return true;
    };
    let u = first.normalize() + last.normalize();
    (0..n).all(|i| {
        (i + 1..=n).all(|j| {
            let d = p[j] - p[i];
            d == DVec3::ZERO || u.dot(d) > 0.0
        })
    })
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

    // Face pairs that can meet, and how; a crossing pair with whether each
    // face's boundary has a mark on the other's surface (`Relation`).
    let mut crossing: Vec<(usize, usize, [bool; 2])> = Vec::new();
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
                Relation::Crossing { zero } => {
                    pts.add_crossings(&pair, Operand::First, fa, b.surface)?;
                    pts.add_crossings(&pair, Operand::Second, fb, a.surface)?;
                    if pts.raw_count() > limits.points {
                        return Err(too_large("candidate points", limits.points));
                    }
                    crossing.push((fa, fb, zero));
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
    let courses = Courses {
        pair: &pair,
        tracks: &contacts.tracks,
    };
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
                            (0..along.len() - 1)
                                .map(|k| {
                                    let (ta, tb) = (ts[k], ts[k + 1]);
                                    segs.add_curved(
                                        along[k],
                                        along[k + 1],
                                        Course::Edge { op, e, ta, tb },
                                        Source::Edge(op, e),
                                        &courses,
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
    // Where a face's boundary touches a face of the other body.
    let mut touched: Vec<imprint::Touched> = Vec::new();
    for &(fa, fb, zero) in &crossing {
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
            for (lo, hi) in imprint::plane_meeting(&pair, &pts, fa, fb)? {
                for w in pts.split_segment(lo, hi)?.windows(2) {
                    imprint(segs.add(w[0], w[1], Course::Straight, Source::Meet(fa, fb)));
                }
            }
            imprint::touches_in(&pair, &pts, &outlines, fa, fb, zero, &mut touched)?;
            continue;
        }
        let (parts, full) = imprint::stretches(&pair, &pts, &outlines, &contacts, fa, fb)?;
        whole.extend(full.into_iter().map(|id| (fa, fb, id)));
        for st in parts {
            let track = contacts.tracks[st.track];
            if track.straight() {
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
                    Source::Meet(fa, fb),
                    &courses,
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
                    Source::Meet(fa, fb),
                    &courses,
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
                        let m = match seg.course {
                            Course::Straight => (pts.rep[seg.a] + pts.rep[seg.b]) * 0.5,
                            course => courses.middle(course),
                        };
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
    let isolated = settle_touches(&pair, &pts, &segs, &face_segs, touched)?;

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
    let mut contacts_in: Vec<assemble::Contact> = Vec::new();
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
            let lone = isolated[slot(op)].get(&fi);
            let points: Vec<usize> =
                lone.map_or_else(Vec::new, |l| l.iter().map(|x| x.0).collect());
            let split =
                planar::regions(&edges, &pos, &reach, eps, &points).map_err(|p| match p {
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
                    planar::Problem::Far(c) => KernelError::InvalidInput(format!(
                        "{} near {}",
                        planar::FAR,
                        show(pts.rep[c])
                    )),
                })?;
            if clear && !split.loose.is_empty() {
                return Err(KernelError::InvalidResult(format!(
                    "boolean: face {fi} of the {} body: a contact in a face clear of the \
                     other body",
                    name(op)
                )));
            }
            // The contacts each region holds (design note 0012 §18.2): as
            // places for the rule after selection, as obstacles its test
            // point keeps away from, and with the sides and the wedge the
            // cross-checks compare its class with.
            let base = pieces.len();
            let mut obstacles: Vec<Vec<classify::Obstacle>> =
                (0..split.regions.len()).map(|_| Vec::new()).collect();
            let mut held: Vec<Vec<Held>> = (0..split.regions.len()).map(|_| Vec::new()).collect();
            for l in &split.loose {
                let holder = base + l.region;
                match &l.kind {
                    planar::LooseKind::Bridge {
                        seg,
                        path,
                        band,
                        free,
                    } => {
                        let s = &segs.list[*seg];
                        let m = match s.course {
                            Course::Straight => (pts.rep[s.a] + pts.rep[s.b]) * 0.5,
                            course => courses.middle(course),
                        };
                        let e = s
                            .sources
                            .iter()
                            .find_map(|&src| match src {
                                Source::Edge(o, e) if o == other(op) => Some(e),
                                _ => None,
                            })
                            .ok_or_else(|| edgeless_contact(&segs, &pts, *seg, op, m, eps))?;
                        let n = outward(shape, f, m).normalize();
                        let (sigma, convex) = wedge(them, e, n, m, tol)?;
                        contacts_in.push(assemble::Contact {
                            place: assemble::Place::Seg(*seg),
                            holder,
                            sigma,
                            at: m,
                        });
                        for (k, c) in [s.a, s.b].into_iter().enumerate() {
                            if free[k] {
                                contacts_in.push(assemble::Contact {
                                    place: assemble::Place::Point(c),
                                    holder,
                                    sigma,
                                    at: pts.rep[c],
                                });
                            }
                        }
                        obstacles[l.region].push(classify::Obstacle {
                            path: path.clone(),
                            band: *band,
                        });
                        held[l.region].push(Held {
                            at: m,
                            sigma,
                            convex,
                        });
                    }
                    &planar::LooseKind::Point(c) => {
                        let sigma = lone
                            .and_then(|l| l.iter().find(|x| x.0 == c))
                            .map(|x| x.1)
                            .ok_or_else(|| {
                                KernelError::InvalidResult(format!(
                                    "boolean: face {fi} of the {} body: an isolated point that \
                                     was not touched near {}",
                                    name(op),
                                    show(pts.rep[c])
                                ))
                            })?;
                        contacts_in.push(assemble::Contact {
                            place: assemble::Place::Point(c),
                            holder,
                            sigma,
                            at: pts.rep[c],
                        });
                        obstacles[l.region].push(classify::Obstacle {
                            path: vec![pos(c)],
                            band: 0.0,
                        });
                        held[l.region].push(Held {
                            at: pts.rep[c],
                            sigma,
                            convex: None,
                        });
                    }
                }
            }
            for (ri, region) in split.regions.into_iter().enumerate() {
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
                let (inner, clearance, refined) = interior(&region, &edges, &chart, fi, eps)?;
                if clearance <= eps {
                    return Err(KernelError::InvalidInput(format!(
                        "the bodies leave a sliver thinner than the tolerance near {}",
                        show(chart.unmap(inner))
                    )));
                }
                // A region holding contacts is tested at points away from
                // them, best first.
                let tests = if held[ri].is_empty() {
                    vec![chart.unmap(inner)]
                } else {
                    avoiding(
                        refined.as_ref().unwrap_or(&region),
                        &obstacles[ri],
                        &chart,
                        fi,
                        eps,
                        held[ri][0].at,
                    )?
                };
                let at = tests[0];
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
                let boundary = &boundaries[slot(other(op))];
                let class = match class {
                    Some(c) => c,
                    None if held[ri].is_empty() => boundary.classify(at, eps)?,
                    None => holder_class(boundary, &tests, eps, held[ri][0].at)?,
                };
                if !held[ri].is_empty() {
                    cross_check(boundary, shape, f, class, &held[ri], eps)?;
                }
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
    assemble::assemble(
        &geometry,
        &pts,
        segs,
        pieces,
        &contacts_in,
        &coplanar,
        op,
        tol,
    )
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

/// A contact held by a region, for the cross-checks of its class: a point
/// of it, the side of the face the other body's touching boundary lies on
/// (+1 outside), and for a lying edge whether it is convex, if its faces
/// are not too nearly flat to tell.
struct Held {
    at: DVec3,
    sigma: i8,
    convex: Option<bool>,
}

/// Test points a region holding contacts tries, at most.
const HOLDER_TRIES: usize = 4;

/// The error for a contact edge (a bridge, segment `seg`, its middle `m`)
/// that carries no edge of the other body. Where a straight piece of an
/// edge of the other body leaves one end of the bridge and ends within
/// twice the tolerance of the bridge's line without reaching its other
/// end, the edge lies on the face and leaves it where the tolerance cannot
/// place the point: the edge's crossing of the face's boundary and the end
/// of the curve where the faces meet are different points, a few
/// tolerances apart where the face's boundary meets its neighbour at a
/// slant (found by the contact fuzz). That is refused as input; anything
/// else is a kernel error, a crossing curve left dangling.
fn edgeless_contact(
    segs: &Segments,
    pts: &Points,
    seg: usize,
    op: Operand,
    m: DVec3,
    eps: f64,
) -> KernelError {
    let s = &segs.list[seg];
    let (a, b) = (pts.rep[s.a], pts.rep[s.b]);
    let line = (b - a).normalize_or_zero();
    let off = |q: DVec3| (q - a).reject_from_normalized(line).length();
    let leaves = s.course == Course::Straight
        && segs.list.iter().any(|t| {
            let from = [t.a, t.b].iter().position(|&c| c == s.a || c == s.b);
            t.course == Course::Straight
                && from.is_some_and(|k| {
                    let end = [t.b, t.a][k];
                    end != s.a && end != s.b && off(pts.rep[end]) <= 2.0 * eps
                })
                && t.sources
                    .iter()
                    .any(|x| matches!(*x, Source::Edge(o, _) if o == other(op)))
        });
    if leaves {
        KernelError::InvalidInput(format!(
            "an edge of one body lies on a face of the other and leaves it at a point the \
             tolerance cannot place near {}",
            show(b)
        ))
    } else {
        KernelError::InvalidResult(format!(
            "boolean: a contact edge carries no edge of the other body near {}",
            show(m)
        ))
    }
}

/// For an edge `e` of `shape` lying on a face of the other operand whose
/// unit outward normal at `m` (a point of the edge) is `n`: the side of
/// that face the edge's two faces leave it to (+1 outward), and whether
/// the edge is convex, `None` where its faces are too nearly flat to tell.
/// The direction into each face across the edge is its outward normal
/// crossed with the edge's direction as the face's loop runs it (material
/// on the left); one within the angular tolerance of the touched face's
/// plane is refused.
fn wedge(
    shape: &Shape,
    e: usize,
    n: DVec3,
    m: DVec3,
    tol: &Tolerance,
) -> Result<(i8, Option<bool>), KernelError> {
    let along = if shape.straight[e] {
        let (a, b) = shape.edge_points(e);
        (b - a).normalize()
    } else {
        let (curve, t0, t1) = shape.curve(e);
        let (t, _) = points::foot_on(curve, t0, t1, m);
        curve.deriv(t).normalize() * (t1 - t0).signum()
    };
    let &[(f1, fw1), (f2, fw2)] = shape.edge_faces(e) else {
        return Err(KernelError::InvalidResult(format!(
            "boolean: an edge lying on a face is not on two faces near {}",
            show(m)
        )));
    };
    let mut sides = [0i8; 2];
    let mut into = [DVec3::ZERO; 2];
    let mut normals = [DVec3::ZERO; 2];
    for (k, (fi, fw)) in [(f1, fw1), (f2, fw2)].into_iter().enumerate() {
        let nk = outward(shape, &shape.faces[fi], m).normalize();
        let r = nk.cross(if fw { along } else { -along });
        let d = n.dot(r);
        if !(d.abs() >= tol.angular) {
            return Err(KernelError::InvalidInput(format!(
                "the bodies touch at an angle too shallow for the tolerance near {}",
                show(m)
            )));
        }
        sides[k] = if d > 0.0 { 1 } else { -1 };
        into[k] = r;
        normals[k] = nk;
    }
    if sides[0] != sides[1] {
        return Err(KernelError::InvalidResult(format!(
            "boolean: the faces at a lying edge leave the face it lies on to both sides near {}",
            show(m)
        )));
    }
    let c = normals[0].dot(into[1]);
    Ok((sides[0], (c.abs() >= tol.angular).then_some(c < 0.0)))
}

/// Test points in `region` of face `face` away from its contacts
/// (`obstacles`), best first, each farther than `eps` in space from the
/// region's boundary and from the contacts
/// ([`classify::interior_points_avoiding`]). With none, the contacts leave
/// no room, and the boolean is refused at the contact `near`.
fn avoiding(
    region: &planar::Region,
    obstacles: &[classify::Obstacle],
    chart: &Chart,
    face: usize,
    eps: f64,
    near: DVec3,
) -> Result<Vec<DVec3>, KernelError> {
    let band = region.band();
    let found = classify::interior_points_avoiding(&region.polygons(), obstacles, HOLDER_TRIES)
        .map_err(|m| KernelError::InvalidResult(format!("boolean: face {face}: {m}")))?;
    let tests: Vec<DVec3> = found
        .iter()
        .filter(|&&(q, _, score)| chart.space(q, score - band) > eps)
        .map(|&(q, _, _)| chart.unmap(q))
        .collect();
    if tests.is_empty() {
        return Err(KernelError::InvalidInput(format!(
            "a contact of the bodies lies within the tolerance of a narrow piece of a face \
             near {}",
            show(near)
        )));
    }
    Ok(tests)
}

/// The class of a region holding contacts, at the first of its test points
/// that the other body's surface is not within the tolerance of: those
/// can be near the touching body, at a distance from the contact about the
/// tolerance over the sine of the angle the body rises at. With none, the
/// angle is too shallow to classify, and the boolean is refused at the
/// contact `near`.
fn holder_class(
    boundary: &classify::Boundary,
    tests: &[DVec3],
    eps: f64,
    near: DVec3,
) -> Result<Class, KernelError> {
    for &p in tests {
        match boundary.classify(p, eps) {
            Err(KernelError::InvalidInput(_)) => continue,
            found => return found,
        }
    }
    Err(KernelError::InvalidInput(format!(
        "the bodies touch at an angle too shallow for the tolerance near {}",
        show(near)
    )))
}

/// Checks the class of a region of face `f` of `shape` holding the
/// contacts `held` against two derivations that do not use its test point
/// (design note 0012 §18.2): it is not on a face of the other body; a point
/// 4 tolerances across the face from each contact, on the side away from
/// the touching body, has the same class against the other body (skipped
/// where that point is within the tolerance of its surface); and at a
/// lying edge, the region is outside the other body exactly when the edge
/// is convex. Any disagreement is a kernel error.
fn cross_check(
    boundary: &classify::Boundary,
    shape: &Shape,
    f: &shape::Face,
    class: Class,
    held: &[Held],
    eps: f64,
) -> Result<(), KernelError> {
    if !matches!(class, Class::Inside | Class::Outside) {
        return Err(KernelError::InvalidResult(format!(
            "boolean: a contact lies inside a piece on a face of the other body near {}",
            show(held[0].at)
        )));
    }
    for h in held {
        let disagree = || {
            KernelError::InvalidResult(format!(
                "boolean: a contact's pieces disagree on their side near {}",
                show(h.at)
            ))
        };
        let n = outward(shape, f, h.at).normalize();
        let across = h.at - n * (4.0 * eps * f64::from(h.sigma));
        match boundary.classify(across, eps) {
            Ok(c) if c != class => return Err(disagree()),
            Ok(_) | Err(KernelError::InvalidInput(_)) => {}
            Err(e) => return Err(e),
        }
        if h.convex
            .is_some_and(|convex| convex != (class == Class::Outside))
        {
            return Err(disagree());
        }
    }
    Ok(())
}

/// The isolated points of each face that the other body's boundary
/// touches (design note 0012 §18.2), per operand by face, each with its
/// cluster and the side the touching boundary lies on (+1 outside), by
/// cluster. A point that ends a segment of the touched face is left to the
/// face's graph and dropped: a point of its boundary or of a contact
/// segment, or a vertex where a face crosses through a reflex corner
/// (whose corners the vertex fans decide); but a station of an edge where
/// the meeting of one of the edge's own faces with the touched face ends
/// is a crossing tangent to the edge, not built yet. A point that ends no
/// segment, touched from both sides, is an inconsistency.
fn settle_touches(
    pair: &Pair,
    pts: &Points,
    segs: &Segments,
    face_segs: &[Vec<Vec<FaceSeg>>; 2],
    mut touched: Vec<imprint::Touched>,
) -> Result<[FxHashMap<usize, Vec<(usize, i8)>>; 2], KernelError> {
    let key = |t: &imprint::Touched| (slot(t.op), t.face, t.cluster, t.mark, t.side);
    touched.sort_by_key(key);
    touched.dedup_by_key(|t| key(t));
    let mut out = [FxHashMap::default(), FxHashMap::default()];
    // The ends of the segments of the face last looked at, with the
    // segments.
    let mut ends: Vec<(usize, usize)> = Vec::new();
    let mut ends_of = None;
    let mut i = 0;
    while i < touched.len() {
        let t = touched[i];
        let same = |x: &imprint::Touched| (x.op, x.face, x.cluster) == (t.op, t.face, t.cluster);
        let group = &touched[i..i + touched[i..].iter().take_while(|x| same(x)).count()];
        i += group.len();
        if ends_of != Some((t.op, t.face)) {
            ends.clear();
            for fs in &face_segs[slot(t.op)][t.face] {
                let s = &segs.list[fs.seg];
                ends.extend([(s.a, fs.seg), (s.b, fs.seg)]);
            }
            ends.sort_unstable();
            ends.dedup();
            ends_of = Some((t.op, t.face));
        }
        let from = ends.partition_point(|&(c, _)| c < t.cluster);
        let at: Vec<usize> = ends[from..]
            .iter()
            .take_while(|&&(c, _)| c == t.cluster)
            .map(|&(_, s)| s)
            .collect();
        if !at.is_empty() {
            let them = pair.get(other(t.op));
            for x in group {
                let imprint::Mark::Station(e, _) = x.mark else {
                    continue;
                };
                let meets_own = at.iter().any(|&s| {
                    segs.list[s].sources.iter().any(|&src| match (src, t.op) {
                        (Source::Meet(g, h), Operand::First)
                        | (Source::Meet(h, g), Operand::Second) => {
                            g == t.face && them.edge_faces(e).iter().any(|&(k, _)| k == h)
                        }
                        _ => false,
                    })
                });
                if meets_own {
                    return Err(KernelError::Unsupported(format!(
                        "booleans where an edge of one body touches a curve where the faces \
                         meet (a crossing tangent to the edge; design note 0012 §18) near {}",
                        show(pts.rep[t.cluster])
                    )));
                }
            }
            continue;
        }
        if group.iter().any(|x| x.side != t.side) {
            return Err(KernelError::InvalidResult(format!(
                "boolean: one point touched from both sides near {}",
                show(pts.rep[t.cluster])
            )));
        }
        out[slot(t.op)]
            .entry(t.face)
            .or_insert_with(Vec::new)
            .push((t.cluster, t.side));
    }
    Ok(out)
}

/// Points on a piece's boundary above which [`interior`] probes the piece
/// ([`classify::probe_point`]) before triangulating it.
const PROBE_ABOVE: usize = 64;

/// A point inside a region of face `face` and its distance in space from
/// the region's boundary at least: inside the polygon of the region's
/// chords, as far from its sides as can be found, less the chords' bands,
/// as a distance in space ([`Chart::space`]).
/// While that is within the tolerance, the chords whose sags are widest
/// are halved, until every sag is within a sixteenth of it or down to its
/// own rounding ([`planar::Sag::wide`]). A sag that rounding holds over
/// the sixteenth is refused ([`planar::FAR`]): the bodies are too far from
/// the origin for the tolerance (found by review: without the floors the
/// halving went on forever, the chords doubling each round). Also returns
/// the region as refined, if it was.
fn interior(
    region: &planar::Region,
    edges: &[planar::Edge],
    chart: &Chart,
    face: usize,
    eps: f64,
) -> Result<(DVec2, f64, Option<planar::Region>), KernelError> {
    let fine = eps / 16.0;
    let mut refined: Option<planar::Region> = None;
    for _ in 0..planar::MAX_ROUNDS {
        let r = refined.as_ref().unwrap_or(region);
        let polygons = r.polygons();
        let band = r.band();
        // A large piece is probed first: a triangulation of it costs far
        // more, most of all when its points lie on circles.
        if polygons.iter().map(Vec::len).sum::<usize>() > PROBE_ABOVE {
            if let Some((inner, clearance)) = classify::probe_point(&polygons) {
                let space = chart.space(inner, clearance - band);
                if space > eps {
                    return Ok((inner, space, refined));
                }
            }
        }
        let (inner, clearance) = classify::interior_point(&polygons)
            .map_err(|m| KernelError::InvalidResult(format!("boolean: face {face}: {m}")))?;
        let space = chart.space(inner, clearance - band);
        if space > eps {
            return Ok((inner, space, refined));
        }
        let wide = r.widest(fine);
        if !wide.is_finite() {
            return Err(KernelError::InvalidResult(format!(
                "boolean: face {face}: a curved boundary's band could not be bounded"
            )));
        }
        if wide == 0.0 {
            if r.stuck(fine) {
                return Err(KernelError::InvalidInput(format!(
                    "{} near {}",
                    planar::FAR,
                    show(chart.unmap(inner))
                )));
            }
            return Ok((inner, space, refined));
        }
        refined = Some(r.refined(edges, wide / 4.0, fine));
    }
    Err(KernelError::InvalidResult(format!(
        "boolean: face {face}: a curved boundary could not be resolved into chords"
    )))
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
            });
            self.list.len() - 1
        });
        self.source(i, s);
        i
    }

    /// The curved segment from `from` to `to` along `course` (given from
    /// `from`), from source `s`: a new one, or one already between the two
    /// clusters along the same curve within the tolerance `eps` (two pieces
    /// of one course, an operand's edge and the curve where faces meet
    /// along it), which gains the source. The same curve: the points a
    /// quarter, half and three quarters of the way along each course, by
    /// its parameter and, for an ellipse, a conic or a NURBS curve, by the
    /// turning of its tangent ([`Courses::at_turn`]), are within `eps` of
    /// the other's stretch (design note 0012 §16.3). Its parameters and
    /// the other's need not agree (a conic's NURBS edge and its track).
    /// Two different conics between the same ends can pass within the
    /// tolerance of three points of each and stray far between them where
    /// the points crowd: by parameter, on a hyperbola's arms (missing its
    /// vertex); by turning, near its vertex. Review measured two coplanar
    /// conic tracks apart by 4.4 × 10⁵ times the largest distance at their
    /// parameter probes, and by at most 2.3 times it at both sets.
    fn add_curved(
        &mut self,
        from: usize,
        to: usize,
        course: Course,
        s: Source,
        courses: &Courses,
        eps: f64,
    ) -> usize {
        let k = key(from, to);
        // Every probe of `a` within the tolerance of `b`'s stretch (the
        // fractions are symmetric, so a course run either way gives the
        // same points).
        let near = |a: Course, b: Course| {
            [0.25, 0.5, 0.75].iter().all(|&f| {
                courses.distance(b, courses.along(a, f)) <= eps
                    && courses
                        .at_turn(a, f)
                        .is_none_or(|p| courses.distance(b, p) <= eps)
            })
        };
        let found = self.curved.get(&k).and_then(|l| {
            l.iter().copied().find(|&i| {
                let other = self.list[i].course;
                near(course, other) && near(other, course)
            })
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

    /// Points on the ellipse a tilted cut leaves on a rod's side, 1 m to
    /// 3,000 km from the origin, located against that face: on its
    /// boundary, found at once. Far out the sags of the ellipse's chords
    /// stop narrowing at their rounding; locating halved every chord
    /// within its band of the point to the parameter's last bit, the
    /// chords doubling (found by review: from about 145 km).
    #[test]
    fn points_on_a_far_conic_edge_are_located_at_once() {
        use crate::{ExtrudeSpec, Profile};
        let tol = Tolerance::default();
        let eps = tol.linear;
        let tilt = 20f64.to_radians();
        let n = DVec3::new(tilt.sin(), 0.0, tilt.cos());
        for far in [1.0, 3e4, 3e5, 3e6] {
            let at = DVec3::new(far, 0.3 * far, 0.0);
            let k = BrepKernel::default();
            let rod = k
                .extrude(
                    &Profile::circle(
                        Frame {
                            origin: at,
                            ..Frame::WORLD
                        },
                        DVec2::ZERO,
                        0.05,
                    ),
                    &ExtrudeSpec {
                        distance: 0.15,
                        symmetric: false,
                    },
                )
                .expect("a rod");
            let f = Frame::from_normal(at + DVec3::new(0.0, 0.0, 0.1), n, DVec3::Y);
            let slab = k
                .make_box(
                    &Frame {
                        origin: f.origin - (f.x + f.y) * 0.15,
                        ..f
                    },
                    DVec3::new(0.3, 0.3, 0.2),
                )
                .expect("a box");
            let (cut, _) =
                boolean(&rod.brep, &slab.brep, BooleanOp::Subtract, &tol).expect("a cut rod");
            let shape = shape::shape(&cut, "cut").expect("a shape");
            let mut located = 0;
            for face in &shape.faces {
                let outline = chart::Outline::new(&shape, face).expect("an outline");
                if !matches!(outline.chart, Chart::Cylinder { .. }) {
                    continue;
                }
                for &(e, _) in face.loops.iter().flatten() {
                    let (curve, t0, t1) = shape.curve(e);
                    if !matches!(curve, Curve3::Nurbs(_)) {
                        continue;
                    }
                    for f in [0.3, 0.5, 0.7] {
                        let q = outline.chart.map(curve.point(t0 + (t1 - t0) * f));
                        let start = std::time::Instant::now();
                        assert_eq!(outline.locate(q, eps, eps / 16.0), None, "{far}");
                        let took = start.elapsed().as_secs_f64();
                        assert!(took < 1.0, "{far}: {took} s");
                        located += 1;
                    }
                }
            }
            assert!(located >= 3, "{far}: {located}");
        }
    }

    /// A course whose sag does not narrow when its chords are halved.
    struct Stuck;

    impl planar::Bend for Stuck {
        fn at(&self, t: f64) -> DVec2 {
            DVec2::new(1.0 - t, 2e-8)
        }
        fn sag(&self, _: f64, _: f64) -> planar::Sag {
            planar::Sag {
                bound: 1e-8,
                floor: 1e-8,
                settles: false,
            }
        }
        fn slack(&self) -> f64 {
            0.0
        }
    }

    /// A sliver whose curved side's sag stays over a sixteenth of the
    /// tolerance, as the rounding of a course far from the origin keeps
    /// it: refused, its floor showing that halving cannot narrow it (it was
    /// halved forever, the chords doubling each round; found by review).
    #[test]
    fn a_sag_rounding_holds_up_is_refused() {
        let chord = |p: DVec2, edge: usize, curved: bool| planar::RChord {
            p,
            edge,
            t0: 0.0,
            t1: if curved { 1.0 } else { 0.0 },
            band: if curved { 1e-8 } else { 0.0 },
            sag: if curved {
                planar::Bend::sag(&Stuck, 0.0, 1.0)
            } else {
                planar::Sag::ZERO
            },
        };
        let region = planar::Region {
            loops: vec![vec![(0, true), (1, true), (2, true), (3, true)]],
            chords: vec![vec![
                chord(DVec2::ZERO, 0, false),
                chord(DVec2::new(1.0, 0.0), 1, false),
                chord(DVec2::new(1.0, 2e-8), 2, true),
                chord(DVec2::new(0.0, 2e-8), 3, false),
            ]],
        };
        let edges: Vec<planar::Edge> = (0..4)
            .map(|k| planar::Edge {
                seg: k,
                a: k,
                b: (k + 1) % 4,
                side: Some(true),
                bend: (k == 2).then(|| (Box::new(Stuck) as Box<dyn planar::Bend>, vec![0.0, 1.0])),
            })
            .collect();
        let chart = Chart::Plane(Frame::WORLD);
        let got = interior(&region, &edges, &chart, 0, 1e-7).map(|(q, d, _)| (q, d));
        let Err(KernelError::InvalidInput(m)) = got else {
            panic!("{got:?}")
        };
        assert!(m.contains("too far from the origin"), "{m}");
    }

    /// Two hyperbolas with the same asymptotes, their vertices ten
    /// tolerances apart, between ends a fiftieth of the tolerance apart:
    /// every probe of either by parameter falls on its arms, within the
    /// tolerance of the other, and only the probes by turning, near the
    /// vertices, tell them apart (design note 0012 §16.3; found by
    /// review). A hyperbola from either end is still one segment.
    #[test]
    fn courses_apart_only_at_a_vertex_are_two_segments() {
        let box_ = cube(DVec3::ZERO, 1.0);
        let pair = Pair {
            first: shape::shape(&box_, "first").expect("a shape"),
            second: shape::shape(&box_, "second").expect("a shape"),
        };
        let (p, q) = (1e-3, 1e-3 * (1.0 + 1e-3));
        let tracks = [
            Track::Conic {
                vertex: DVec3::ZERO,
                x: DVec3::X,
                y: DVec3::Y,
                p,
                c: -1.0,
            },
            Track::Conic {
                vertex: DVec3::X * (q - p),
                x: DVec3::X,
                y: DVec3::Y,
                p: q,
                c: -1.0,
            },
        ];
        let courses = Courses {
            pair: &pair,
            tracks: &tracks,
        };
        let eps = 1e-7;
        let first = Course::Track {
            id: 0,
            ta: -0.5 / p,
            tb: 0.3 / p,
        };
        let second = Course::Track {
            id: 1,
            ta: -0.5 / q,
            tb: 0.3 / q,
        };
        for (a, b) in [(first, second), (second, first)] {
            for f in [0.0, 0.25, 0.5, 0.75, 1.0] {
                assert!(courses.distance(b, courses.along(a, f)) <= eps / 4.0);
            }
            let vertex = courses.at_turn(a, 0.5).expect("a conic");
            assert!(courses.distance(b, vertex) > 9.0 * eps);
        }
        let mut segs = Segments::default();
        let one = segs.add_curved(0, 1, first, Source::Meet(0, 0), &courses, eps);
        let two = segs.add_curved(0, 1, second, Source::Meet(0, 1), &courses, eps);
        assert_ne!(one, two);
        let again = segs.add_curved(1, 0, first.reversed(), Source::Meet(1, 0), &courses, eps);
        assert_eq!(again, one);
        assert_eq!(segs.list.len(), 2);
    }

    /// Probes by turning along stretches that turn more than a half turn,
    /// either way: an ellipse track's, and as an operand's edge a NURBS
    /// ellipse in four quarters and an arc of 238° of another written as
    /// one cubic span with positive weights. They sit at a quarter, a half
    /// and three quarters of the turning, by the ellipse's closed form
    /// (they crowded into the stretch's first part; found by review, the
    /// cubic by the review after).
    #[test]
    fn probes_by_turning_follow_stretches_past_a_half_turn() {
        use std::f64::consts::FRAC_1_SQRT_2;
        // The tangent's direction at angle t of the ellipse with semi-axes
        // a and b along x and y, lifted: within a quarter turn of the
        // circle's, t + π/2.
        let turn = |(a, b): (f64, f64), t: f64| {
            let (s, c) = t.sin_cos();
            t + ((a - b) * s * c).atan2(a * s * s + b * c * c)
        };
        // The ellipse's angle of a point of it, lifted to the stretch
        // from angle `from` the way `up` says.
        let angle = |(a, b): (f64, f64), p: DVec3, from: f64, up: bool| {
            let t = (p.y / b).atan2(p.x / a);
            if up {
                from + (t - from).rem_euclid(TAU)
            } else {
                from - (from - t).rem_euclid(TAU)
            }
        };
        let check = |ab: (f64, f64), at: &dyn Fn(f64) -> DVec3, ta: f64, tb: f64, up: bool| {
            let from = angle(ab, at(0.0), ta, up);
            let to = angle(ab, at(1.0), from + if up { 1e-9 } else { -1e-9 }, up);
            let whole = turn(ab, to) - turn(ab, from);
            assert!(whole.abs() > PI, "{ta} {tb}: {whole}");
            for f in [0.25, 0.5, 0.75] {
                let t = angle(ab, at(f), from, up);
                let got = (turn(ab, t) - turn(ab, from)) / whole;
                assert!((got - f).abs() <= 1e-9, "{ta} {tb}: {got} for {f}");
            }
        };
        let (a, b) = (4.0, 1.0);
        let box_ = cube(DVec3::ZERO, 1.0);
        let pair = Pair {
            first: shape::shape(&box_, "first").expect("a shape"),
            second: shape::shape(&box_, "second").expect("a shape"),
        };
        let tracks = [Track::Ellipse {
            center: DVec3::ZERO,
            x: DVec3::X,
            y: DVec3::Y,
            a,
            b,
        }];
        let courses = Courses {
            pair: &pair,
            tracks: &tracks,
        };
        for (ta, tb) in [(0.0, 1.5 * PI), (1.5 * PI, 0.0), (0.3, 0.3 + 1.9 * PI)] {
            let course = Course::Track { id: 0, ta, tb };
            let at = |f: f64| match f {
                0.0 => tracks[0].point(ta),
                1.0 => tracks[0].point(tb),
                f => courses.at_turn(course, f).expect("an ellipse"),
            };
            check((a, b), &at, ta, tb, tb > ta);
        }
        // The ellipse in four quarters, each a rational quadratic.
        let w = FRAC_1_SQRT_2;
        let corner = |x: f64, y: f64| DVec3::new(a * x, b * y, 0.0);
        let quarters = opendatum_nurbs::NurbsCurve::new(
            2,
            vec![
                0.0, 0.0, 0.0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1.0, 1.0, 1.0,
            ],
            vec![
                corner(1.0, 0.0),
                corner(1.0, 1.0),
                corner(0.0, 1.0),
                corner(-1.0, 1.0),
                corner(-1.0, 0.0),
                corner(-1.0, -1.0),
                corner(0.0, -1.0),
                corner(1.0, -1.0),
                corner(1.0, 0.0),
            ],
            Some(vec![1.0, w, 1.0, w, 1.0, w, 1.0, w, 1.0]),
        )
        .expect("an ellipse");
        // The arc of the ellipse with semi-axes 1.5 and 1 from −119° to
        // 119° through 0°: a rational quadratic with the middle weight
        // cos 119° < 0 (the arc of a circle stretched), raised to a cubic
        // (Piegl & Tiller §5.5), whose weights 1, 0.0101, 0.0101, 1 are
        // positive.
        let (s, c) = 119f64.to_radians().sin_cos();
        let (a2, b2) = (1.5, 1.0);
        let q = [
            glam::DVec4::new(a2 * c, -b2 * s, 0.0, 1.0),
            glam::DVec4::new(a2, 0.0, 0.0, c),
            glam::DVec4::new(a2 * c, b2 * s, 0.0, 1.0),
        ];
        let cubic = opendatum_nurbs::NurbsCurve::from_homogeneous(
            3,
            vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0],
            vec![
                q[0],
                (q[0] + q[1] * 2.0) / 3.0,
                (q[1] * 2.0 + q[2]) / 3.0,
                q[2],
            ],
        )
        .expect("an arc");
        // The quadratics are not halved; the cubic's one span is.
        assert_eq!(
            knot_breaks(&quarters, 0.05, 0.85),
            [0.05, 0.25, 0.5, 0.75, 0.85]
        );
        assert!(knot_breaks(&cubic, 0.0, 1.0).len() > 2);
        for (n, ab, stretches) in [
            (quarters, (a, b), [(0.05, 0.85), (0.85, 0.05), (0.7, 0.05)]),
            (cubic, (a2, b2), [(0.0, 1.0), (1.0, 0.0), (0.1, 0.9)]),
        ] {
            // The curve as edge 0 of the first operand.
            let curve = Curve3::Nurbs(n);
            let mut brep = box_.clone();
            brep.edges[0].curve = curve.clone();
            let pair = Pair {
                first: shape::shape(&brep, "first").expect("a shape"),
                second: shape::shape(&box_, "second").expect("a shape"),
            };
            let courses = Courses {
                pair: &pair,
                tracks: &[],
            };
            for (ta, tb) in stretches {
                let course = Course::Edge {
                    op: Operand::First,
                    e: 0,
                    ta,
                    tb,
                };
                let at = |f: f64| match f {
                    0.0 => curve.point(ta),
                    1.0 => curve.point(tb),
                    f => courses.at_turn(course, f).expect("a NURBS ellipse"),
                };
                check(ab, &at, ta, tb, tb > ta);
            }
        }
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
