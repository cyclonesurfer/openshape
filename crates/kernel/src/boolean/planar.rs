//! The pieces of a face cut along segments (design note 0012 §4.2 and
//! §10.4).
//!
//! The segments join given points (clusters), so nothing is intersected
//! here. Two segments crossing, or a point within the tolerance of a
//! segment it does not end, means the bodies come closer than the
//! tolerance somewhere the earlier stages found no contact (seen in the
//! face's chart rather than in space), and the boolean is refused there.
//!
//! The split is made in a chart of the face's surface. A point's
//! tolerance there is its [`reach`](super::chart::Chart::reach), the chart
//! distance beyond which points are farther than the tolerance from it in
//! space (the tolerance itself in a plane's or a developed surface's
//! chart, more toward the rim of a sphere's); bands are measured in the
//! chart, which no chart makes shorter than in space. A straight segment is one chord; a curved
//! one is a chain of chords through points of its course, each with a
//! *band*, a bound on how far the course strays from the chord. Chords are
//! halved wherever a band could hide a contact (another chord or a point
//! within the band, plus the tolerance for a point) until every such band
//! is under a sixteenth of the tolerance (where the rounding of a bound
//! keeps it wider, the bodies lying too far from the origin, the face is
//! refused), so the chords form a straight-line
//! graph with the same structure as the segments: two courses that cross
//! have chords that cross, and a decision taken on the chords is off only
//! within that sixteenth (a refusal for points within 1.06 rather than 1
//! tolerances, never a wrong region).
//!
//! The chords form a planar graph; walking it with the face on the left of
//! every half-edge gives its cycles. Counter-clockwise cycles bound regions,
//! clockwise ones are the outsides of connected components and become holes
//! of the smallest region of another component around them. A region is
//! inside the face when it runs along the face's own edges in their
//! direction (material on the left); regions bounded only by segments
//! imprinted inside the face are inside it too.
//!
//! What separates nothing is reported apart from the regions (design note
//! 0012 §18.2): a *bridge*, a segment with the same cycle on both sides (an
//! edge of one body lying on the face, from a free end, between two loops
//! or on its own), and an isolated point (a corner of one body resting on
//! the face). Neither enters a region's loops; each is given with the
//! region that holds it, for the booleans' rule on contacts.

use std::cell::OnceCell;

use glam::{DVec2, DVec3};
use rustc_hash::FxHashMap;

use super::flat_outline::{FlatOutline, LONG_LOOP};
use crate::bvh::Bvh;

/// Why a face could not be split, at a point (a cluster).
#[derive(Debug)]
pub(super) enum Problem {
    /// A point within the tolerance of a segment it does not end.
    Close(usize),
    /// Two segments crossing.
    Cross(usize),
    /// The earlier stages disagree: a kernel error.
    Inconsistent(&'static str, usize),
    /// A chord whose band could hide a contact cannot be narrowed: its
    /// sag's own rounding holds it over the tolerance's sixteenth (it is
    /// over that and not [`halvable`]), the bodies being too far from the
    /// origin (design note 0012 §16.2).
    Far(usize),
}

impl Problem {
    /// What went wrong, and where (a cluster).
    pub fn describe(&self) -> (&'static str, usize) {
        match *self {
            Problem::Close(c) => (
                "a point within the tolerance of a segment it does not end",
                c,
            ),
            Problem::Cross(c) => ("two segments crossing", c),
            Problem::Inconsistent(m, c) => (m, c),
            Problem::Far(c) => (FAR, c),
        }
    }
}

/// Why a curved course that rounding keeps from being resolved is refused.
pub(super) const FAR: &str = "the bodies lie too far from the origin for the tolerance: \
                              rounding keeps a curved edge's course from being resolved \
                              within it";

/// A bound on how far a course strays from the chord of a stretch of it
/// ([`Bend::sag`]).
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub(super) struct Sag {
    /// The bound.
    pub bound: f64,
    /// The part of `bound` that its terms' allowances for their own
    /// rounding make, which halving the chord does not narrow: as the chord
    /// shrinks, the bound comes to at most about three times this (the
    /// allowances, and as much again of rounding in the terms they allow
    /// for), or to this itself where it [`Sag::settles`]. Zero where the
    /// bound narrows to nothing.
    pub floor: f64,
    /// Whether the bound comes down to its floor itself as the chord
    /// shrinks: its other part is computed from the chord's parameters
    /// with rounding only in proportion to it, and narrows to nothing (a
    /// conic track's, design note 0012 §16.2).
    pub settles: bool,
}

impl Sag {
    /// A straight chord's.
    pub const ZERO: Sag = Sag {
        bound: 0.0,
        floor: 0.0,
        settles: false,
    };

    /// Whether halving the chord can still narrow the bound toward `fine`:
    /// it is over `fine` and over four times its floor (an unbounded one
    /// is, to be halved toward a bound). A bound that settles on its floor
    /// is halved while it is over `fine` and over its floor by more than
    /// [`SETTLED`] of it: each halving quarters that part, so the chord
    /// comes within `fine` wherever the floor is under it by that much.
    pub fn wide(self, fine: f64) -> bool {
        if self.settles {
            self.bound > fine && self.bound > self.floor * (1.0 + SETTLED)
        } else {
            self.bound > fine.max(4.0 * self.floor)
        }
    }
}

/// How far over its floor a bound that settles on it ([`Sag::settles`])
/// is halved: a thousandth of it, about six halvings past four floors.
pub(super) const SETTLED: f64 = 1.0 / 1024.0;

/// Whether the chord of a course on `[t0, t1]` with `sag` is to be halved
/// toward `fine`: [`Sag::wide`], and its parameters have a middle between
/// them.
pub(super) fn halvable(sag: Sag, t0: f64, t1: f64, fine: f64) -> bool {
    let tm = 0.5 * (t0 + t1);
    sag.wide(fine) && tm != t0 && tm != t1
}

/// The course of a curved segment in a face's chart.
pub(super) trait Bend {
    /// The chart point at parameter `t`.
    fn at(&self, t: f64) -> DVec2;
    /// A bound on the distance from the course on `[t0, t1]` (either
    /// order) to the chord between its ends, less [`Bend::slack`], from the
    /// course's bend: halving the chord narrows it, by about four times
    /// where the course is smooth and the bound from its derivatives, by
    /// about two where it is from control points, down to a few times its
    /// floor ([`Sag`]).
    fn sag(&self, t0: f64, t1: f64) -> Sag;
    /// What every such bound adds to the sag and no halving narrows: an
    /// allowance for the rounding of the chart's points (design note 0012
    /// §12.1).
    fn slack(&self) -> f64;
}

/// A segment of a face's graph.
pub(super) struct Edge<'a> {
    /// The segment's index, which the regions are given in.
    pub seg: usize,
    /// Its ends (clusters), `a` the lesser.
    pub a: usize,
    pub b: usize,
    /// Material on the left of `a → b` (`Some(true)`), of `b → a`
    /// (`Some(false)`), or on both sides (`None`, an imprint inside the
    /// face).
    pub side: Option<bool>,
    /// For a curved segment, its course and its parameters from `a` to
    /// `b`, at least its two ends: a first division into chords.
    pub bend: Option<(Box<dyn Bend + 'a>, Vec<f64>)>,
}

/// A chord of a region's boundary: from `p` to the next chord's `p`, along
/// an edge (by its position in the list given), between its parameters
/// `t0` and `t1` (both 0 for a straight edge), straying at most `band`.
#[derive(Clone, Copy, Debug)]
pub(super) struct RChord {
    pub p: DVec2,
    pub edge: usize,
    pub t0: f64,
    pub t1: f64,
    pub band: f64,
    /// The part of `band` from the course's bend ([`Band`]).
    pub sag: Sag,
}

/// A region inside the face.
pub(super) struct Region {
    /// Its cycles of segments (outer first, counter-clockwise; holes
    /// clockwise), each segment with whether it is run from `a` to `b`.
    /// Each cycle starts at a cluster.
    pub loops: Vec<Vec<(usize, bool)>>,
    /// Per cycle, its chords in order.
    pub chords: Vec<Vec<RChord>>,
}

impl Region {
    /// The cycles as polygons of chart points.
    pub fn polygons(&self) -> Vec<Vec<DVec2>> {
        self.chords
            .iter()
            .map(|l| l.iter().map(|c| c.p).collect())
            .collect()
    }

    /// The widest band of its chords.
    pub fn band(&self) -> f64 {
        self.chords.iter().flatten().fold(0.0, |m, c| m.max(c.band))
    }

    /// The widest sag of its chords that halving can still narrow toward
    /// `fine` ([`halvable`]), zero if none can be.
    pub fn widest(&self, fine: f64) -> f64 {
        self.chords
            .iter()
            .flatten()
            .filter(|c| halvable(c.sag, c.t0, c.t1, fine))
            .fold(0.0, |m, c| m.max(c.sag.bound))
    }

    /// Whether a chord's sag is over `fine` though halving can no longer
    /// narrow it ([`halvable`]): rounding holds it up, the course lying too
    /// far from the origin for the tolerance.
    pub fn stuck(&self, fine: f64) -> bool {
        self.chords
            .iter()
            .flatten()
            .any(|c| c.sag.bound > fine && !halvable(c.sag, c.t0, c.t1, fine))
    }

    /// The region with every chord halved whose sag is over `above` and
    /// can still be narrowed toward `fine` ([`halvable`]).
    pub fn refined(&self, edges: &[Edge], above: f64, fine: f64) -> Region {
        let chords = self
            .chords
            .iter()
            .map(|l| {
                let mut out = Vec::with_capacity(2 * l.len());
                for (i, c) in l.iter().enumerate() {
                    let next = l[(i + 1) % l.len()].p;
                    match &edges[c.edge].bend {
                        Some((bend, _))
                            if c.sag.bound > above && halvable(c.sag, c.t0, c.t1, fine) =>
                        {
                            let tm = 0.5 * (c.t0 + c.t1);
                            let pm = bend.at(tm);
                            let first = Band::of(bend.as_ref(), c.t0, tm, c.p, pm);
                            let second = Band::of(bend.as_ref(), tm, c.t1, pm, next);
                            out.push(RChord {
                                p: c.p,
                                edge: c.edge,
                                t0: c.t0,
                                t1: tm,
                                band: first.band,
                                sag: first.sag,
                            });
                            out.push(RChord {
                                p: pm,
                                edge: c.edge,
                                t0: tm,
                                t1: c.t1,
                                band: second.band,
                                sag: second.sag,
                            });
                        }
                        _ => out.push(*c),
                    }
                }
                out
            })
            .collect();
        Region {
            loops: self.loops.clone(),
            chords,
        }
    }
}

/// A face's regions, and what lies in them without bounding them.
pub(super) struct Split {
    /// The regions inside the face.
    pub regions: Vec<Region>,
    /// Bridges, then isolated points, each with the region holding it.
    pub loose: Vec<Loose>,
}

/// A bridge or an isolated point of a face's graph, in a region.
#[derive(Debug)]
pub(super) struct Loose {
    /// The holding region, by its position in [`Split::regions`].
    pub region: usize,
    pub kind: LooseKind,
}

/// What a [`Loose`] item is.
#[derive(Debug)]
pub(super) enum LooseKind {
    /// A segment with the same region on both sides.
    Bridge {
        /// The segment's index, as given.
        seg: usize,
        /// Its chords' points in the chart, from its `a` to its `b`.
        path: Vec<DVec2>,
        /// The widest band of its chords.
        band: f64,
        /// Whether its `a` and its `b` are free: no segment of the
        /// region's loops ends there.
        free: [bool; 2],
    },
    /// A point (a cluster) on no segment.
    Point(usize),
}

/// The band of the chord from `p` to `q` along a course on `[t0, t1]`.
#[derive(Clone, Copy, Debug)]
struct Band {
    /// How far the course strays from the chord: its band about its own
    /// chord, plus how far the chord's ends are moved from the course's
    /// (to the clusters' places, at the ends of a segment).
    band: f64,
    /// The course's band about its own chord ([`Bend::sag`]). The rest of
    /// the band, the rounding allowance and an end moved to its cluster (by
    /// up to the tolerance), stays with the chord however short it gets, so
    /// halving is decided by the sag alone: halving by the whole band never
    /// ended at a cluster moved more than `fine` from the course, nor, far
    /// from the origin, where the allowance is wider than `fine` (found by
    /// review).
    sag: Sag,
}

impl Band {
    fn of(bend: &dyn Bend, t0: f64, t1: f64, p: DVec2, q: DVec2) -> Band {
        let moved = (bend.at(t0) - p).length().max((bend.at(t1) - q).length());
        let sag = bend.sag(t0, t1);
        Band {
            band: sag.bound + bend.slack() + moved,
            sag,
        }
    }
}

/// Rounds of halving after which a curved segment is taken as unresolved,
/// and the most halvings of a chord at a node ([`split_ends`]) or of the
/// chords about a point sought inside a region. A halving narrows a sag
/// about four times where it is bounded from the course's derivatives and
/// about two times where from control points, down to a few times its
/// floor ([`Bend::sag`]); sixty halvings shorten a chord 2⁶⁰ ≈ 10¹⁸ times,
/// past the last bit of its parameters, so a sag stops narrowing at its
/// floor or at the parameter's last bit ([`halvable`]) before this is
/// reached.
pub(super) const MAX_ROUNDS: usize = 60;

/// A chord of the graph.
#[derive(Clone, Copy, Debug)]
struct Chord {
    from: usize,
    to: usize,
    /// The unique edge it is part of.
    edge: usize,
    t0: f64,
    t1: f64,
    band: f64,
    /// The part of `band` from the course's bend ([`Band`]).
    sag: Sag,
}

/// The chords of the edges, given each curved edge's parameters: nodes
/// are the clusters first, then a node per inner parameter of each edge.
struct Graph {
    p: Vec<DVec2>,
    /// Per node, its cluster (for an inner node, its edge's start).
    cluster: Vec<usize>,
    /// The number of cluster nodes (the first ones).
    real: usize,
    chords: Vec<Chord>,
}

fn graph(
    edges: &[&Edge],
    params: &[Vec<f64>],
    ends: &[(usize, usize)],
    nodes: &[usize],
    pos: &[DVec2],
) -> Graph {
    let mut p = pos.to_vec();
    let mut cluster = nodes.to_vec();
    let mut chords = Vec::with_capacity(edges.len());
    for (i, e) in edges.iter().enumerate() {
        let (na, nb) = ends[i];
        match &e.bend {
            None => chords.push(Chord {
                from: na,
                to: nb,
                edge: i,
                t0: 0.0,
                t1: 0.0,
                band: 0.0,
                sag: Sag::ZERO,
            }),
            Some((bend, _)) => {
                let ts = &params[i];
                let k = ts.len() - 1;
                let mut prev = na;
                for j in 0..k {
                    let next = if j + 1 == k {
                        nb
                    } else {
                        p.push(bend.at(ts[j + 1]));
                        cluster.push(nodes[na]);
                        p.len() - 1
                    };
                    let b = Band::of(bend.as_ref(), ts[j], ts[j + 1], p[prev], p[next]);
                    chords.push(Chord {
                        from: prev,
                        to: next,
                        edge: i,
                        t0: ts[j],
                        t1: ts[j + 1],
                        band: b.band,
                        sag: b.sag,
                    });
                    prev = next;
                }
            }
        }
    }
    Graph {
        p,
        cluster,
        real: nodes.len(),
        chords,
    }
}

/// Halves the chords whose bands could hide a contact until none can:
/// a chord whose sag is over `eps / 16` is halved when another chord comes
/// within the two bands of it, or a cluster within its band and the
/// cluster's reach `near`. Where such a chord's sag is held over `eps / 16`
/// by its own rounding (it is not [`halvable`]), the face is refused
/// ([`Problem::Far`]) rather than decided on a band wider than allowed.
/// Consecutive chords of one edge share a point of its course and are not
/// tested against each other.
///
/// Two edges meeting at a node always come within their bands there, so an
/// edge's end chord at a node it shares with another edge is halved until
/// its band is narrow, one halving a round; those halvings are made first,
/// all at once ([`split_ends`]). After that, a round tests only the chords
/// new in it against the others (marking either of a close pair that is
/// wide): a pair of chords that both stood unchanged through the round
/// before was tested then, and they were not close.
fn refine(
    edges: &[&Edge],
    params: &mut [Vec<f64>],
    ends: &[(usize, usize)],
    nodes: &[usize],
    pos: &[DVec2],
    near: &[f64],
    eps: f64,
) -> Result<Graph, Problem> {
    let fine = eps / 16.0;
    // The widest reach of a cluster, for the box searches.
    let widest = near.iter().fold(0.0f64, |m, &r| m.max(r));
    split_ends(edges, params, ends, pos, fine);
    // Per parameter of each edge, the round it was made in.
    let mut born: Vec<Vec<usize>> = params.iter().map(|ts| vec![0; ts.len()]).collect();
    for round in 0..MAX_ROUNDS {
        let g = graph(edges, params, ends, nodes, pos);
        // The chords new in this round, in the order `graph` makes them.
        let mut fresh = Vec::with_capacity(g.chords.len());
        for (i, e) in edges.iter().enumerate() {
            match e.bend {
                None => fresh.push(round == 0),
                Some(_) => {
                    let b = &born[i];
                    fresh.extend(b.windows(2).map(|w| w[0] == round || w[1] == round));
                }
            }
        }
        debug_assert_eq!(fresh.len(), g.chords.len());
        let boxes: Vec<(DVec3, DVec3)> = g
            .chords
            .iter()
            .map(|c| {
                let (a, b) = (g.p[c.from], g.p[c.to]);
                let r = DVec2::splat(c.band);
                ((a.min(b) - r).extend(0.0), (a.max(b) + r).extend(0.0))
            })
            .collect();
        let tree = Bvh::new(boxes.clone());
        let points: Vec<DVec3> = g.p[..g.real].iter().map(|q| q.extend(0.0)).collect();
        let point_tree = Bvh::of_points(&points);
        let mut split: Vec<bool> = vec![false; g.chords.len()];
        // A chord whose band could hide a contact is halved; one whose sag
        // rounding holds over `fine` cannot be, and the face is refused.
        let halve = |k: usize, split: &mut [bool]| {
            let c = g.chords[k];
            if c.sag.bound > fine && !split[k] {
                if !halvable(c.sag, c.t0, c.t1, fine) {
                    return Err(Problem::Far(g.cluster[c.from]));
                }
                split[k] = true;
            }
            Ok(())
        };
        let mut close = Vec::new();
        for i in (0..g.chords.len()).filter(|&i| fresh[i]) {
            let c = g.chords[i];
            let (a, b) = (g.p[c.from], g.p[c.to]);
            tree.overlapping(boxes[i].0, boxes[i].1, &mut close);
            for &j in &close {
                let d = g.chords[j];
                let consecutive = d.edge == c.edge && (d.to == c.from || d.from == c.to);
                // Nothing to decide when neither could be halved anew.
                if j == i
                    || consecutive
                    || (split[i] || c.sag.bound <= fine) && (split[j] || d.sag.bound <= fine)
                {
                    continue;
                }
                if segment_distance(a, b, g.p[d.from], g.p[d.to]) <= c.band + d.band {
                    halve(i, &mut split)?;
                    halve(j, &mut split)?;
                }
            }
            if c.sag.bound > fine && !split[i] {
                let r = DVec3::splat(c.band + widest);
                point_tree.overlapping(boxes[i].0 - r, boxes[i].1 + r, &mut close);
                if close.iter().any(|&n| {
                    n != c.from && n != c.to && point_segment(g.p[n], a, b) <= c.band + near[n]
                }) {
                    halve(i, &mut split)?;
                }
            }
        }
        if !split.iter().any(|&s| s) {
            return Ok(g);
        }
        for (i, c) in g.chords.iter().enumerate() {
            if split[i] {
                let ts = &mut params[c.edge];
                let at = ts
                    .iter()
                    .position(|&t| t == c.t1)
                    .expect("a chord ends at a parameter of its edge");
                ts.insert(at, 0.5 * (c.t0 + c.t1));
                born[c.edge].insert(at, round + 1);
            }
        }
    }
    let first = edges.first().map_or(0, |e| e.a);
    Err(Problem::Inconsistent(
        "a curved segment could not be resolved into chords",
        first,
    ))
}

/// Halves the end chords of each curved edge at a node it shares with
/// another edge until their sags are at most `fine`, or rounding holds
/// them over it ([`halvable`]), or [`MAX_ROUNDS`] halvings are made: the
/// halvings [`refine`] would make there anyway,
/// one a round, because chords meeting at a node are within their bands
/// of each other.
fn split_ends(
    edges: &[&Edge],
    params: &mut [Vec<f64>],
    ends: &[(usize, usize)],
    pos: &[DVec2],
    fine: f64,
) {
    // Per node, how many edge ends meet there.
    let mut degree = vec![0usize; pos.len()];
    for &(a, b) in ends {
        degree[a] += 1;
        degree[b] += 1;
    }
    for (i, e) in edges.iter().enumerate() {
        let Some((bend, _)) = &e.bend else { continue };
        let (na, nb) = ends[i];
        // An edge from a node back to it meets another edge there only if
        // a third end is there.
        let own = if na == nb { 2 } else { 1 };
        for (node, at_start) in [(na, true), (nb, false)] {
            if degree[node] <= own {
                continue;
            }
            let ts = &mut params[i];
            let (t_end, t_next) = if at_start {
                (ts[0], ts[1])
            } else {
                (ts[ts.len() - 1], ts[ts.len() - 2])
            };
            let mut inner = Vec::new();
            let mut t = t_next;
            for _ in 0..MAX_ROUNDS {
                if !halvable(bend.sag(t_end, t), t_end, t, fine) {
                    break;
                }
                t = 0.5 * (t_end + t);
                inner.push(t);
            }
            // `inner` runs from the end's neighbour toward the end.
            if at_start {
                inner.reverse();
                ts.splice(1..1, inner);
            } else {
                let k = ts.len() - 1;
                ts.splice(k..k, inner);
            }
        }
    }
}

/// The regions inside the face whose segments are `edges`, the clusters
/// at chart points `pos`, and its bridges and the isolated `points`
/// (clusters on no segment) in them; `reach` gives the tolerance `eps` in
/// the chart at a chart point.
///
/// The half-edges are walked once; a chord whose two half-edges lie in one
/// cycle is a bridge (in a plane graph an edge has the same face on both
/// sides exactly when removing it disconnects its component), and with
/// bridges the walk is made again without them. A bridge lies in the
/// region whose loops hold the other half-edges of its first cycle, or,
/// in a cycle of bridges alone (a tree), in the region around one of its
/// points, as an isolated point does. A bridge or a point in no region, or
/// in a hole, is off the face and dropped.
/// With no bridges and no points the result is the one walk's.
pub(super) fn regions(
    edges: &[Edge],
    pos: &dyn Fn(usize) -> DVec2,
    reach: &dyn Fn(DVec2) -> f64,
    eps: f64,
    points: &[usize],
) -> Result<Split, Problem> {
    // Unique segments, in order of first appearance, with their sides.
    let mut unique: Vec<usize> = Vec::new();
    let mut slot: FxHashMap<usize, usize> = FxHashMap::default();
    let mut side: Vec<Option<bool>> = Vec::new();
    for (i, e) in edges.iter().enumerate() {
        match slot.get(&e.seg) {
            None => {
                slot.insert(e.seg, unique.len());
                unique.push(i);
                side.push(e.side);
            }
            Some(&k) => match (side[k], e.side) {
                (_, None) => {}
                (None, Some(x)) => side[k] = Some(x),
                (Some(x), Some(y)) if x == y => {}
                _ => {
                    return Err(Problem::Inconsistent(
                        "an edge with material on both sides",
                        e.a,
                    ))
                }
            },
        }
    }
    if unique.is_empty() {
        return Ok(Split {
            regions: Vec::new(),
            loose: Vec::new(),
        });
    }
    let list: Vec<&Edge> = unique.iter().map(|&i| &edges[i]).collect();
    // Nodes: the segments' ends, then the isolated points.
    let mut node_of: FxHashMap<usize, usize> = FxHashMap::default();
    let mut nodes: Vec<usize> = Vec::new();
    for e in &list {
        for c in [e.a, e.b] {
            node_of.entry(c).or_insert_with(|| {
                nodes.push(c);
                nodes.len() - 1
            });
        }
    }
    let ends_count = nodes.len();
    for &c in points {
        if node_of.contains_key(&c) {
            return Err(Problem::Inconsistent(
                "a point touched off the segments ends one",
                c,
            ));
        }
        node_of.insert(c, nodes.len());
        nodes.push(c);
    }
    let pos_nodes: Vec<DVec2> = nodes.iter().map(|&c| pos(c)).collect();
    let near: Vec<f64> = pos_nodes.iter().map(|&q| reach(q)).collect();
    let ends: Vec<(usize, usize)> = list
        .iter()
        .map(|e| (node_of[&e.a], node_of[&e.b]))
        .collect();
    let mut params: Vec<Vec<f64>> = list
        .iter()
        .map(|e| e.bend.as_ref().map(|(_, t)| t.clone()).unwrap_or_default())
        .collect();
    let g = if list.iter().any(|e| e.bend.is_some()) {
        refine(&list, &mut params, &ends, &nodes, &pos_nodes, &near, eps)?
    } else {
        graph(&list, &params, &ends, &nodes, &pos_nodes)
    };
    let p = &g.p;
    let chord_ends: Vec<(usize, usize)> = g.chords.iter().map(|c| (c.from, c.to)).collect();

    consistent(&chord_ends, p, g.real, &g.cluster, &near)?;

    // Half-edges: 2s runs from → to, 2s + 1 runs to → from.
    let nh = 2 * chord_ends.len();
    let from = |h: usize| {
        if h % 2 == 0 {
            chord_ends[h / 2].0
        } else {
            chord_ends[h / 2].1
        }
    };
    let to = |h: usize| {
        if h % 2 == 0 {
            chord_ends[h / 2].1
        } else {
            chord_ends[h / 2].0
        }
    };
    let mut out: Vec<Vec<usize>> = vec![Vec::new(); p.len()];
    for h in 0..nh {
        out[from(h)].push(h);
    }
    let angle = |h: usize| {
        let d = p[to(h)] - p[from(h)];
        d.y.atan2(d.x)
    };
    for (n, list) in out.iter_mut().enumerate() {
        list.sort_by(|&x, &y| angle(x).total_cmp(&angle(y)));
        for w in list.windows(2) {
            if angle(w[0]) == angle(w[1]) {
                return Err(Problem::Close(g.cluster[n]));
            }
        }
    }
    let unclosed =
        |h: usize| Problem::Inconsistent("a half-edge walk did not close", g.cluster[from(h)]);
    let (first_cycles, first_of) = walk(&out, nh, &from).map_err(unclosed)?;
    // Bridges: chords with one cycle on both sides.
    let bridge: Vec<bool> = (0..chord_ends.len())
        .map(|c| first_of[2 * c] == first_of[2 * c + 1])
        .collect();
    let bridged = bridge.iter().any(|&b| b);
    let (cycles, cycle_of, first) = if bridged {
        for (c, &b) in bridge.iter().enumerate() {
            if b && side[g.chords[c].edge].is_some() {
                return Err(Problem::Inconsistent(
                    "an edge of the face bounds no region",
                    g.cluster[chord_ends[c].0],
                ));
            }
        }
        for list in &mut out {
            list.retain(|&h| !bridge[h / 2]);
        }
        let (cycles, of) = walk(&out, nh, &from).map_err(unclosed)?;
        (cycles, of, Some((first_cycles, first_of)))
    } else {
        (first_cycles, first_of, None)
    };
    debug_assert!(
        cycles.iter().all(|c| prune(c.clone()).len() == c.len()),
        "a half-edge followed by its twin after the bridges are removed"
    );
    let area = |c: &[usize]| -> f64 {
        c.iter()
            .map(|&h| p[from(h)].perp_dot(p[to(h)]))
            .sum::<f64>()
            * 0.5
    };

    // Components, by union–find over nodes, through chords that are not
    // bridges.
    let mut comp: Vec<usize> = (0..p.len()).collect();
    fn find(c: &mut [usize], mut i: usize) -> usize {
        while c[i] != i {
            c[i] = c[c[i]];
            i = c[i];
        }
        i
    }
    for (k, &(a, b)) in chord_ends.iter().enumerate() {
        if bridge[k] {
            continue;
        }
        let (ra, rb) = (find(&mut comp, a), find(&mut comp, b));
        if ra != rb {
            comp[ra.max(rb)] = ra.min(rb);
        }
    }

    // Which side of the face each cycle is on, from its boundary half-edges.
    let inside_of = |c: &[usize]| -> Result<Option<bool>, Problem> {
        let mut verdict = None;
        for &h in c {
            if let Some(s) = side[g.chords[h / 2].edge] {
                let along = (h % 2 == 0) == s;
                match verdict {
                    None => verdict = Some(along),
                    Some(v) if v != along => {
                        return Err(Problem::Inconsistent(
                            "a region both inside and outside the face",
                            g.cluster[from(h)],
                        ))
                    }
                    _ => {}
                }
            }
        }
        Ok(verdict)
    };

    struct Cell {
        cycle: Vec<usize>,
        area: f64,
        comp: usize,
        inside: bool,
        holes: Vec<Vec<usize>>,
    }
    let mut cells: Vec<Cell> = Vec::new();
    // Per cycle, the cell it bounds (as its outer cycle or a hole).
    let mut cell_of: Vec<Option<usize>> = vec![None; cycles.len()];
    let mut outsides: Vec<(usize, Vec<usize>, f64)> = Vec::new();
    for (ci, c) in cycles.into_iter().enumerate() {
        if c.is_empty() {
            continue;
        }
        let a = area(&c);
        if a > 0.0 {
            let inside = inside_of(&c)?.unwrap_or(true);
            let k = find(&mut comp, from(c[0]));
            cell_of[ci] = Some(cells.len());
            cells.push(Cell {
                cycle: c,
                area: a,
                comp: k,
                inside,
                holes: Vec::new(),
            });
        } else if a < 0.0 {
            outsides.push((ci, c, a));
        } else {
            return Err(Problem::Inconsistent(
                "a cycle of zero area",
                g.cluster[from(c[0])],
            ));
        }
    }
    // Bridges, by segment: their chords in order and their first cycle.
    let mut bridges: Vec<(usize, Vec<usize>)> = Vec::new();
    for (k, &b) in bridge.iter().enumerate() {
        if !b {
            continue;
        }
        let e = g.chords[k].edge;
        match bridges.last_mut() {
            Some((last, chords)) if *last == e => chords.push(k),
            _ => bridges.push((e, vec![k])),
        }
    }
    // The bridges whose first cycle holds no other half-edge (trees),
    // placed by a point as isolated points are.
    let first_cycle = |chord: usize| -> &[usize] {
        let (cycles, of) = first.as_ref().expect("bridges were found in a first walk");
        &cycles[of[2 * chord]]
    };
    let tree = |chords: &[usize]| first_cycle(chords[0]).iter().all(|&h| bridge[h / 2]);
    let placed_by_point = bridges.iter().any(|(_, ch)| tree(ch)) || nodes.len() > ends_count;
    // The regions around a point are among those whose boxes hold it
    // (padded against the rounding of the containment test).
    let pad = 1e-9 * scale(p);
    let boxes = (!outsides.is_empty() || placed_by_point).then(|| {
        Bvh::new(
            cells
                .iter()
                .map(|r| {
                    let (lo, hi) = r.cycle.iter().fold(
                        (DVec2::splat(f64::INFINITY), DVec2::splat(f64::NEG_INFINITY)),
                        |(lo, hi), &h| (lo.min(p[from(h)]), hi.max(p[from(h)])),
                    );
                    (lo.extend(0.0), hi.extend(0.0))
                })
                .collect(),
        )
    });
    // Whether `q` is inside cell `i`'s outer cycle: a long cycle through
    // its outline, with a tree of its sides built when a point is first
    // placed against it, so a point meets the sides near it rather than
    // every side (found by profiling: a face crossing the other's plane
    // thousands of times, whose cells' cycles each run along all the
    // crossings).
    let outlines: Vec<OnceCell<FlatOutline>> = cells.iter().map(|_| OnceCell::new()).collect();
    let in_cycle = |i: usize, q: DVec2, cells: &[Cell]| -> bool {
        let c = &cells[i].cycle;
        if c.len() <= LONG_LOOP {
            return contains(c, q, p, &from);
        }
        let outline = outlines[i].get_or_init(|| {
            let corner = |k: usize| p[from(c[k % c.len()])];
            let sides = (0..c.len()).map(|k| (corner(k), corner(k + 1))).collect();
            FlatOutline::new(sides, vec![0, c.len()])
        });
        outline.inside(q, &mut 0)
    };
    // The smallest cell around `q` (farther than its reach from every
    // chord), other than of component `skip`.
    let around = |q: DVec2, skip: Option<usize>, cells: &[Cell]| -> Option<usize> {
        let mut near: Vec<usize> = Vec::new();
        if let Some(t) = &boxes {
            let q3 = q.extend(0.0);
            t.overlapping(q3 - DVec3::splat(pad), q3 + DVec3::splat(pad), &mut near);
        }
        let mut best: Option<usize> = None;
        for &i in &near {
            let r = &cells[i];
            if Some(r.comp) == skip || !in_cycle(i, q, cells) {
                continue;
            }
            if best.is_none_or(|b| r.area < cells[b].area) {
                best = Some(i);
            }
        }
        best
    };
    for (ci, c, _) in outsides {
        let k = find(&mut comp, from(c[0]));
        // A cluster of the cycle: farther than the tolerance from every
        // segment it does not end.
        let first = c
            .iter()
            .map(|&h| from(h))
            .find(|&n| n < g.real)
            .expect("a cycle passes through a cluster");
        let best = around(p[first], Some(k), &cells);
        let verdict = inside_of(&c)?;
        match best {
            Some(i) => {
                if !cells[i].inside || verdict == Some(false) {
                    return Err(Problem::Inconsistent(
                        "a hole outside the face",
                        g.cluster[from(c[0])],
                    ));
                }
                cells[i].holes.push(c);
                cell_of[ci] = Some(i);
            }
            None => {
                if verdict != Some(false) {
                    return Err(Problem::Inconsistent(
                        "an outer boundary inside the face",
                        g.cluster[from(c[0])],
                    ));
                }
            }
        }
    }
    // The regions are the cells inside the face, in order.
    let mut region_of: Vec<Option<usize>> = vec![None; cells.len()];
    let mut count = 0;
    for (i, r) in cells.iter().enumerate() {
        if r.inside {
            region_of[i] = Some(count);
            count += 1;
        }
    }
    // The region around a point of no cycle: the smallest cell around it,
    // if inside the face and the point is in none of its holes.
    let holding = |q: DVec2| -> Option<usize> {
        let i = around(q, None, &cells)?;
        let in_hole = cells[i].holes.iter().any(|h| contains(h, q, p, &from));
        if in_hole {
            None
        } else {
            region_of[i]
        }
    };
    let mut loose: Vec<Loose> = Vec::new();
    // Per node, whether a chord of the regions' loops reaches it.
    let mut looped = vec![false; p.len()];
    for (k, &(a, b)) in chord_ends.iter().enumerate() {
        if !bridge[k] {
            looped[a] = true;
            looped[b] = true;
        }
    }
    for (e, chords) in &bridges {
        let (na, nb) = ends[*e];
        let region = if tree(chords) {
            holding(p[na])
        } else {
            first_cycle(chords[0])
                .iter()
                .find(|&&h| !bridge[h / 2])
                .and_then(|&h| cell_of[cycle_of[h]])
                .and_then(|i| region_of[i])
        };
        // A bridge off the face is dropped: the imprint of a coincident
        // face takes the pieces of its edges within the tolerance of this
        // face, and one can reach just outside it (a copy of a body moved
        // by just over the tolerance, found by the motions test).
        let Some(region) = region else {
            continue;
        };
        let mut path = Vec::with_capacity(chords.len() + 1);
        path.push(p[chord_ends[chords[0]].0]);
        path.extend(chords.iter().map(|&k| p[chord_ends[k].1]));
        loose.push(Loose {
            region,
            kind: LooseKind::Bridge {
                seg: list[*e].seg,
                path,
                band: chords.iter().fold(0.0f64, |m, &k| m.max(g.chords[k].band)),
                free: [!looped[na], !looped[nb]],
            },
        });
    }
    loose.sort_by_key(|l| match l.kind {
        LooseKind::Bridge { seg, .. } => seg,
        LooseKind::Point(_) => unreachable!("bridges come first"),
    });
    let mut isolated: Vec<usize> = (ends_count..nodes.len()).collect();
    isolated.sort_by_key(|&n| nodes[n]);
    for n in isolated {
        if let Some(region) = holding(p[n]) {
            loose.push(Loose {
                region,
                kind: LooseKind::Point(nodes[n]),
            });
        }
    }
    // Each cycle as segments and chords, from its first cluster.
    let as_region = |cycles: Vec<&Vec<usize>>| -> Region {
        let mut loops = Vec::with_capacity(cycles.len());
        let mut chords = Vec::with_capacity(cycles.len());
        for c in cycles {
            let start = c
                .iter()
                .position(|&h| from(h) < g.real)
                .expect("a cycle passes through a cluster");
            let mut segs: Vec<(usize, bool)> = Vec::new();
            let mut cs: Vec<RChord> = Vec::with_capacity(c.len());
            for k in 0..c.len() {
                let h = c[(start + k) % c.len()];
                let ch = g.chords[h / 2];
                let forward = h % 2 == 0;
                let e = list[ch.edge];
                if segs.last() != Some(&(e.seg, forward)) || from(h) < g.real {
                    segs.push((e.seg, forward));
                }
                let (t0, t1) = if forward {
                    (ch.t0, ch.t1)
                } else {
                    (ch.t1, ch.t0)
                };
                cs.push(RChord {
                    p: p[from(h)],
                    edge: unique[ch.edge],
                    t0,
                    t1,
                    band: ch.band,
                    sag: ch.sag,
                });
            }
            loops.push(segs);
            chords.push(cs);
        }
        Region { loops, chords }
    };
    let regions = cells
        .iter()
        .filter(|r| r.inside)
        .map(|r| as_region(std::iter::once(&r.cycle).chain(&r.holes).collect()))
        .collect();
    Ok(Split { regions, loose })
}

/// The cycles of the half-edges in the angle-sorted lists `out` (per node,
/// the half-edges leaving it; `start` gives a half-edge's node), each
/// half-edge followed by the one before its twin about the twin's node:
/// the face on the left of every half-edge. Half-edges in no list are in no
/// cycle. Returns the cycles, in order of their least half-edge, and per
/// half-edge its cycle (`usize::MAX` for none), or a half-edge whose walk
/// did not close.
fn walk(
    out: &[Vec<usize>],
    nh: usize,
    start: &dyn Fn(usize) -> usize,
) -> Result<(Vec<Vec<usize>>, Vec<usize>), usize> {
    let mut rank = vec![usize::MAX; nh];
    for list in out {
        for (r, &h) in list.iter().enumerate() {
            rank[h] = r;
        }
    }
    let next = |h: usize| {
        let t = h ^ 1;
        let list = &out[start(t)];
        list[(rank[t] + list.len() - 1) % list.len()]
    };
    let mut of = vec![usize::MAX; nh];
    let mut cycles: Vec<Vec<usize>> = Vec::new();
    for h0 in 0..nh {
        if of[h0] != usize::MAX || rank[h0] == usize::MAX {
            continue;
        }
        let mut cyc = Vec::new();
        let mut h = h0;
        while of[h] == usize::MAX {
            of[h] = cycles.len();
            cyc.push(h);
            h = next(h);
        }
        if h != h0 {
            return Err(h0);
        }
        cycles.push(cyc);
    }
    Ok((cycles, of))
}

/// Segments below this count are checked pair by pair; above it, through
/// bounding-volume trees. Both find the same first problem.
const SMALL: usize = 24;

/// The largest coordinate magnitude, for padding boxes against rounding.
fn scale(p: &[DVec2]) -> f64 {
    p.iter()
        .fold(0.0f64, |m, q| m.max(q.x.abs()).max(q.y.abs()))
}

/// No two chords cross, and no cluster (the first `real` nodes) lies within
/// its reach `near` of a chord it does not end. Problems are reported in
/// the order of a full pairwise scan (by chord, then crossings before
/// close points, each by index), whether the candidates come from the scan
/// or from the trees.
fn consistent(
    ends: &[(usize, usize)],
    p: &[DVec2],
    real: usize,
    cluster: &[usize],
    near: &[f64],
) -> Result<(), Problem> {
    let boxes: Vec<(DVec3, DVec3)> = ends
        .iter()
        .map(|&(a, b)| (p[a].min(p[b]).extend(0.0), p[a].max(p[b]).extend(0.0)))
        .collect();
    // A computed distance within a reach puts the point within that reach
    // of the chord's box, up to rounding.
    let widest = near.iter().fold(0.0f64, |m, &r| m.max(r));
    let r = widest + 1e-9 * (widest + scale(p));
    let trees = (ends.len() > SMALL).then(|| {
        let points: Vec<DVec3> = p[..real].iter().map(|q| q.extend(0.0)).collect();
        (Bvh::new(boxes.clone()), Bvh::of_points(&points))
    });
    let mut segs: Vec<usize> = Vec::new();
    let mut points: Vec<usize> = Vec::new();
    for (i, &(a, b)) in ends.iter().enumerate() {
        let (lo, hi) = boxes[i];
        match &trees {
            Some((st, pt)) => {
                // Chords that cross share a point, so their boxes meet.
                st.overlapping(lo, hi, &mut segs);
                pt.overlapping(lo - DVec3::splat(r), hi + DVec3::splat(r), &mut points);
            }
            None => {
                segs.clear();
                segs.extend(0..ends.len());
                points.clear();
                points.extend(0..real);
            }
        }
        for &j in &segs {
            let (c, d) = ends[j];
            if j <= i || a == c || a == d || b == c || b == d {
                continue;
            }
            let o1 = orient(p[a], p[b], p[c]);
            let o2 = orient(p[a], p[b], p[d]);
            let o3 = orient(p[c], p[d], p[a]);
            let o4 = orient(p[c], p[d], p[b]);
            if o1 * o2 < 0.0 && o3 * o4 < 0.0 {
                return Err(Problem::Cross(cluster[a]));
            }
        }
        for &n in &points {
            if n != a && n != b && point_segment(p[n], p[a], p[b]) <= near[n] {
                return Err(Problem::Close(cluster[n]));
            }
        }
    }
    Ok(())
}

/// Removes antennas: a half-edge followed by its twin. Only checked: with
/// the bridges removed, no cycle has one.
fn prune(mut c: Vec<usize>) -> Vec<usize> {
    loop {
        let n = c.len();
        if n < 2 {
            return c;
        }
        let Some(i) = (0..n).find(|&i| c[(i + 1) % n] == c[i] ^ 1) else {
            return c;
        };
        let j = (i + 1) % n;
        if j == 0 {
            c.remove(i);
            c.remove(0);
        } else {
            c.remove(j);
            c.remove(i);
        }
    }
}

fn orient(a: DVec2, b: DVec2, c: DVec2) -> f64 {
    robust::orient2d(
        robust::Coord { x: a.x, y: a.y },
        robust::Coord { x: b.x, y: b.y },
        robust::Coord { x: c.x, y: c.y },
    )
}

pub(super) fn point_segment(q: DVec2, a: DVec2, b: DVec2) -> f64 {
    let d = b - a;
    let len2 = d.length_squared();
    if len2 == 0.0 {
        return (q - a).length();
    }
    let t = ((q - a).dot(d) / len2).clamp(0.0, 1.0);
    (q - (a + d * t)).length()
}

/// The distance between the segments `a`–`b` and `c`–`d`: zero if they
/// cross, else the least distance from an end of one to the other.
fn segment_distance(a: DVec2, b: DVec2, c: DVec2, d: DVec2) -> f64 {
    let o1 = orient(a, b, c);
    let o2 = orient(a, b, d);
    let o3 = orient(c, d, a);
    let o4 = orient(c, d, b);
    if o1 * o2 <= 0.0 && o3 * o4 <= 0.0 {
        return 0.0;
    }
    point_segment(a, c, d)
        .min(point_segment(b, c, d))
        .min(point_segment(c, a, b))
        .min(point_segment(d, a, b))
}

/// True if `q` is inside the polygon of half-edge cycle `c` (`q` is not on
/// it: the caller's points keep the tolerance from every segment).
fn contains(c: &[usize], q: DVec2, p: &[DVec2], from: &dyn Fn(usize) -> usize) -> bool {
    let mut inside = false;
    for (i, &h) in c.iter().enumerate() {
        let a = p[from(h)];
        let b = p[from(c[(i + 1) % c.len()])];
        if (a.y > q.y) != (b.y > q.y) {
            let x = a.x + (q.y - a.y) / (b.y - a.y) * (b.x - a.x);
            if x > q.x {
                inside = !inside;
            }
        }
    }
    inside
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-7;

    /// A face's graph on the points `at` (cluster `i` at `at[i]`), with
    /// straight segments numbered in order.
    struct Face {
        at: Vec<DVec2>,
        segs: Vec<(usize, usize, Option<bool>)>,
    }

    impl Face {
        /// The square (0, 0)–(4, 4), counter-clockwise, as clusters 0 to 3.
        fn square() -> Face {
            let mut f = Face {
                at: vec![
                    DVec2::new(0.0, 0.0),
                    DVec2::new(4.0, 0.0),
                    DVec2::new(4.0, 4.0),
                    DVec2::new(0.0, 4.0),
                ],
                segs: Vec::new(),
            };
            for i in 0..4 {
                f.side(i, (i + 1) % 4);
            }
            f
        }

        fn point(&mut self, x: f64, y: f64) -> usize {
            self.at.push(DVec2::new(x, y));
            self.at.len() - 1
        }

        /// A segment of the face's own boundary, material on the left of
        /// `from → to`.
        fn side(&mut self, from: usize, to: usize) -> usize {
            let s = (from.min(to), from.max(to), Some(from < to));
            self.segs.push(s);
            self.segs.len() - 1
        }

        /// A segment imprinted inside the face.
        fn imprint(&mut self, a: usize, b: usize) -> usize {
            self.segs.push((a.min(b), a.max(b), None));
            self.segs.len() - 1
        }

        /// The square hole (1, 1)–(2, 2), as clusters from 4 on.
        fn hole(&mut self) -> [usize; 4] {
            let h = [
                self.point(1.0, 1.0),
                self.point(2.0, 1.0),
                self.point(2.0, 2.0),
                self.point(1.0, 2.0),
            ];
            for i in 0..4 {
                // Clockwise: the face's material outside the hole.
                self.side(h[(i + 1) % 4], h[i]);
            }
            h
        }

        fn split(&self, points: &[usize]) -> Result<Split, Problem> {
            let edges: Vec<Edge> = self
                .segs
                .iter()
                .enumerate()
                .map(|(seg, &(a, b, side))| Edge {
                    seg,
                    a,
                    b,
                    side,
                    bend: None,
                })
                .collect();
            regions(&edges, &|c| self.at[c], &|_| EPS, EPS, points)
        }
    }

    fn loops(s: &Split) -> Vec<Vec<Vec<(usize, bool)>>> {
        s.regions.iter().map(|r| r.loops.clone()).collect()
    }

    fn bridge(l: &Loose) -> (usize, usize, Vec<DVec2>, [bool; 2]) {
        match &l.kind {
            LooseKind::Bridge {
                seg, path, free, ..
            } => (l.region, *seg, path.clone(), *free),
            LooseKind::Point(c) => panic!("a point, {c}, where a bridge was expected"),
        }
    }

    fn inconsistent(r: Result<Split, Problem>) -> &'static str {
        match r {
            Err(Problem::Inconsistent(m, _)) => m,
            Err(p) => panic!("{p:?}"),
            Ok(s) => panic!("split: {} regions, {:?}", s.regions.len(), s.loose),
        }
    }

    #[test]
    fn a_face_alone_is_one_region_with_nothing_loose() {
        let s = Face::square().split(&[]).expect("a split");
        assert_eq!(
            loops(&s),
            vec![vec![vec![(0, true), (1, true), (2, true), (3, false)]]]
        );
        assert!(s.loose.is_empty());
    }

    #[test]
    fn an_antenna_is_a_bridge_and_leaves_the_region_as_it_was() {
        let plain = Face::square().split(&[]).expect("a split");
        let mut f = Face::square();
        let inner = f.point(1.0, 1.0);
        let seg = f.imprint(0, inner);
        let s = f.split(&[]).expect("a split");
        assert_eq!(loops(&s), loops(&plain));
        assert_eq!(s.loose.len(), 1);
        assert_eq!(
            bridge(&s.loose[0]),
            (0, seg, vec![DVec2::ZERO, DVec2::ONE], [false, true])
        );
    }

    #[test]
    fn a_segment_on_its_own_is_a_bridge_free_at_both_ends() {
        let mut f = Face::square();
        let (a, b) = (f.point(1.0, 1.0), f.point(3.0, 2.0));
        let seg = f.imprint(a, b);
        let s = f.split(&[]).expect("a split");
        assert_eq!(s.regions.len(), 1);
        assert_eq!(
            bridge(&s.loose[0]),
            (0, seg, vec![DVec2::ONE, DVec2::new(3.0, 2.0)], [true, true])
        );
    }

    #[test]
    fn a_tree_is_bridges_in_the_region_around_it() {
        let mut f = Face::square();
        let (a, b, c, d) = (
            f.point(1.0, 1.0),
            f.point(2.0, 1.0),
            f.point(2.0, 2.0),
            f.point(3.0, 1.0),
        );
        let segs = [f.imprint(a, b), f.imprint(b, c), f.imprint(b, d)];
        let s = f.split(&[]).expect("a split");
        assert_eq!(s.regions.len(), 1);
        let got: Vec<_> = s.loose.iter().map(bridge).collect();
        assert_eq!(got.len(), 3);
        for (g, seg) in got.iter().zip(segs) {
            assert_eq!((g.0, g.1, g.3), (0, seg, [true, true]));
        }
    }

    #[test]
    fn a_connector_between_a_hole_and_the_outer_loop_is_a_bridge() {
        let mut plain = Face::square();
        plain.hole();
        let plain = plain.split(&[]).expect("a split");
        assert_eq!(plain.regions[0].loops.len(), 2);
        let mut f = Face::square();
        let h = f.hole();
        let seg = f.imprint(0, h[0]);
        let s = f.split(&[]).expect("a split");
        // Found by review: the walk's pruning kept such a connector as a
        // slit in the region's loops.
        assert_eq!(loops(&s), loops(&plain));
        assert_eq!(
            bridge(&s.loose[0]),
            (0, seg, vec![DVec2::ZERO, DVec2::ONE], [false, false])
        );
    }

    #[test]
    fn a_chord_across_the_face_is_no_bridge() {
        let mut f = Face::square();
        let (a, b) = (f.point(2.0, 0.0), f.point(2.0, 4.0));
        // The bottom and top sides end at the chord's ends.
        f.segs.clear();
        for (x, y) in [(0, a), (a, 1), (1, 2), (2, b), (b, 3), (3, 0)] {
            f.side(x, y);
        }
        f.imprint(a, b);
        let s = f.split(&[]).expect("a split");
        assert_eq!(s.regions.len(), 2);
        assert!(s.loose.is_empty());
    }

    #[test]
    fn points_inside_the_face_are_kept_and_others_dropped() {
        let mut f = Face::square();
        f.hole();
        let inside = f.point(0.5, 3.0);
        let in_hole = f.point(1.5, 1.5);
        let outside = f.point(5.0, 5.0);
        let s = f.split(&[outside, in_hole, inside]).expect("a split");
        assert_eq!(s.regions.len(), 1);
        assert_eq!(s.loose.len(), 1);
        assert!(matches!(s.loose[0].kind, LooseKind::Point(c) if c == inside));
        assert_eq!(s.loose[0].region, 0);
    }

    #[test]
    fn a_point_within_the_tolerance_of_a_side_is_refused() {
        let mut f = Face::square();
        let near = f.point(2.0, 0.5 * EPS);
        assert!(matches!(f.split(&[near]), Err(Problem::Close(_))));
    }

    #[test]
    fn a_point_on_a_segment_end_is_an_inconsistency() {
        assert_eq!(
            inconsistent(Face::square().split(&[2])),
            "a point touched off the segments ends one"
        );
    }

    #[test]
    fn a_side_of_the_face_bounding_no_region_is_an_inconsistency() {
        let mut f = Face::square();
        let inner = f.point(1.0, 1.0);
        f.side(0, inner);
        assert_eq!(
            inconsistent(f.split(&[])),
            "an edge of the face bounds no region"
        );
    }

    #[test]
    fn a_bridge_off_the_face_is_dropped() {
        // Outside the outer loop, on its own or from a corner; and inside
        // a hole, from its corner.
        let mut away = Face::square();
        let (a, b) = (away.point(5.0, 5.0), away.point(6.0, 5.0));
        away.imprint(a, b);
        let mut spur = Face::square();
        let out = spur.point(5.0, 5.0);
        spur.imprint(2, out);
        let mut into = Face::square();
        let h = into.hole();
        let mid = into.point(1.5, 1.5);
        into.imprint(h[0], mid);
        for f in [away, spur, into] {
            let s = f.split(&[]).expect("a split");
            assert_eq!(s.regions.len(), 1);
            assert!(s.loose.is_empty(), "{:?}", s.loose);
        }
    }

    #[test]
    fn loose_items_come_bridges_first_by_segment_then_points_by_cluster() {
        let mut f = Face::square();
        let (a, b) = (f.point(1.0, 1.0), f.point(1.0, 2.0));
        let (c, d) = (f.point(3.0, 1.0), f.point(3.0, 2.0));
        let (p, q) = (f.point(2.0, 3.0), f.point(2.0, 3.5));
        let late = f.imprint(c, d);
        let early = f.imprint(a, b);
        assert!(early > late);
        let s = f.split(&[q, p]).expect("a split");
        let order: Vec<String> = s
            .loose
            .iter()
            .map(|l| match l.kind {
                LooseKind::Bridge { seg, .. } => format!("seg {seg}"),
                LooseKind::Point(c) => format!("point {c}"),
            })
            .collect();
        assert_eq!(
            order,
            [
                format!("seg {late}"),
                format!("seg {early}"),
                format!("point {p}"),
                format!("point {q}")
            ]
        );
    }
}
