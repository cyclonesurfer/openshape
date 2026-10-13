//! The pieces of a face cut along segments (design note 0012 §4.2 and
//! §10.4).
//!
//! The segments join given points (clusters), so nothing is intersected
//! here. Two segments crossing, or a point within the tolerance of a
//! segment it does not end, means the bodies come closer than the
//! tolerance somewhere the earlier stages found no contact (seen in the
//! face's chart rather than in space), and the boolean is refused there.
//!
//! The split is made in a chart of the face's surface, in which lengths
//! are lengths along the surface. A straight segment is one chord; a curved
//! one is a chain of chords through points of its course, each with a
//! *band*, a bound on how far the course strays from the chord. Chords are
//! halved wherever a band could hide a contact (another chord or a point
//! within the band, plus the tolerance for a point) until every such band
//! is under a sixteenth of the tolerance, so the chords form a straight-line
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

use glam::{DVec2, DVec3};
use rustc_hash::FxHashMap;

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
        }
    }
}

/// The course of a curved segment in a face's chart.
pub(super) trait Bend {
    /// The chart point at parameter `t`.
    fn at(&self, t: f64) -> DVec2;
    /// A bound on the distance from the course on `[t0, t1]` (either
    /// order) to the chord between its ends.
    fn band(&self, t0: f64, t1: f64) -> f64;
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

    /// The region with every chord whose band is over `fine` halved.
    pub fn refined(&self, edges: &[Edge], fine: f64) -> Region {
        let chords = self
            .chords
            .iter()
            .map(|l| {
                let mut out = Vec::with_capacity(2 * l.len());
                for (i, c) in l.iter().enumerate() {
                    let next = l[(i + 1) % l.len()].p;
                    match &edges[c.edge].bend {
                        Some((bend, _)) if c.band > fine => {
                            let tm = 0.5 * (c.t0 + c.t1);
                            let pm = bend.at(tm);
                            out.push(RChord {
                                p: c.p,
                                edge: c.edge,
                                t0: c.t0,
                                t1: tm,
                                band: chord_band(bend.as_ref(), c.t0, tm, c.p, pm),
                            });
                            out.push(RChord {
                                p: pm,
                                edge: c.edge,
                                t0: tm,
                                t1: c.t1,
                                band: chord_band(bend.as_ref(), tm, c.t1, pm, next),
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

/// The band of the chord from `p` to `q` along a course on `[t0, t1]`:
/// the course's band about its own chord, plus how far the chord's ends
/// are moved (to the clusters' places, at the ends of a segment).
fn chord_band(bend: &dyn Bend, t0: f64, t1: f64, p: DVec2, q: DVec2) -> f64 {
    let moved = (bend.at(t0) - p).length().max((bend.at(t1) - q).length());
    bend.band(t0, t1) + moved
}

/// Rounds of halving after which a curved segment is taken as unresolved:
/// each round quarters the bands it halves, so this many rounds take a
/// chord across a kilometre down to well under a sixteenth of the default
/// tolerance.
const MAX_ROUNDS: usize = 60;

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
                    chords.push(Chord {
                        from: prev,
                        to: next,
                        edge: i,
                        t0: ts[j],
                        t1: ts[j + 1],
                        band: chord_band(bend.as_ref(), ts[j], ts[j + 1], p[prev], p[next]),
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
/// a chord whose band is over `eps / 16` is halved when another chord comes
/// within the two bands of it, or a cluster within its band and `eps`.
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
    eps: f64,
) -> Result<Graph, Problem> {
    let fine = eps / 16.0;
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
        let mut near = Vec::new();
        for i in (0..g.chords.len()).filter(|&i| fresh[i]) {
            let c = g.chords[i];
            let (a, b) = (g.p[c.from], g.p[c.to]);
            tree.overlapping(boxes[i].0, boxes[i].1, &mut near);
            for &j in &near {
                let d = g.chords[j];
                let consecutive = d.edge == c.edge && (d.to == c.from || d.from == c.to);
                // Nothing to decide when neither could be halved anew.
                if j == i
                    || consecutive
                    || (split[i] || c.band <= fine) && (split[j] || d.band <= fine)
                {
                    continue;
                }
                if segment_distance(a, b, g.p[d.from], g.p[d.to]) <= c.band + d.band {
                    split[i] |= c.band > fine;
                    split[j] |= d.band > fine;
                }
            }
            if c.band > fine && !split[i] {
                let r = DVec3::splat(c.band + eps);
                point_tree.overlapping(boxes[i].0 - r, boxes[i].1 + r, &mut near);
                split[i] = near.iter().any(|&n| {
                    n != c.from && n != c.to && point_segment(g.p[n], a, b) <= c.band + eps
                });
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
/// another edge until their bands are at most `fine` (or [`MAX_ROUNDS`]
/// halvings are made): the halvings [`refine`] would make there anyway,
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
                if chord_band(bend.as_ref(), t_end, t, pos[node], bend.at(t)) <= fine {
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

/// The regions inside the face.
pub(super) fn regions(
    edges: &[Edge],
    pos: &dyn Fn(usize) -> DVec2,
    eps: f64,
) -> Result<Vec<Region>, Problem> {
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
        return Ok(Vec::new());
    }
    let list: Vec<&Edge> = unique.iter().map(|&i| &edges[i]).collect();
    // Nodes.
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
    let pos_nodes: Vec<DVec2> = nodes.iter().map(|&c| pos(c)).collect();
    let ends: Vec<(usize, usize)> = list
        .iter()
        .map(|e| (node_of[&e.a], node_of[&e.b]))
        .collect();
    let mut params: Vec<Vec<f64>> = list
        .iter()
        .map(|e| e.bend.as_ref().map(|(_, t)| t.clone()).unwrap_or_default())
        .collect();
    let g = if list.iter().any(|e| e.bend.is_some()) {
        refine(&list, &mut params, &ends, &nodes, &pos_nodes, eps)?
    } else {
        graph(&list, &params, &ends, &nodes, &pos_nodes)
    };
    let p = &g.p;
    let chord_ends: Vec<(usize, usize)> = g.chords.iter().map(|c| (c.from, c.to)).collect();

    consistent(&chord_ends, p, g.real, &g.cluster, eps)?;

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
    let mut rank = vec![0usize; nh];
    for (n, list) in out.iter_mut().enumerate() {
        list.sort_by(|&x, &y| angle(x).total_cmp(&angle(y)));
        for w in list.windows(2) {
            if angle(w[0]) == angle(w[1]) {
                return Err(Problem::Close(g.cluster[n]));
            }
        }
        for (r, &h) in list.iter().enumerate() {
            rank[h] = r;
        }
    }
    let next = |h: usize| {
        let t = h ^ 1;
        let v = from(t);
        let list = &out[v];
        list[(rank[t] + list.len() - 1) % list.len()]
    };

    // Cycles, with antennas (segments with the same cycle on both sides)
    // removed.
    let mut seen = vec![false; nh];
    let mut cycles: Vec<Vec<usize>> = Vec::new();
    for h0 in 0..nh {
        if seen[h0] {
            continue;
        }
        let mut cyc = Vec::new();
        let mut h = h0;
        while !seen[h] {
            seen[h] = true;
            cyc.push(h);
            h = next(h);
        }
        if h != h0 {
            return Err(Problem::Inconsistent(
                "a half-edge walk did not close",
                g.cluster[from(h0)],
            ));
        }
        cycles.push(prune(cyc));
    }
    let area = |c: &[usize]| -> f64 {
        c.iter()
            .map(|&h| p[from(h)].perp_dot(p[to(h)]))
            .sum::<f64>()
            * 0.5
    };

    // Components, by union–find over nodes.
    let mut comp: Vec<usize> = (0..p.len()).collect();
    fn find(c: &mut [usize], mut i: usize) -> usize {
        while c[i] != i {
            c[i] = c[c[i]];
            i = c[i];
        }
        i
    }
    for &(a, b) in &chord_ends {
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
    let mut outsides: Vec<(Vec<usize>, f64)> = Vec::new();
    for c in cycles {
        if c.is_empty() {
            continue;
        }
        let a = area(&c);
        if a > 0.0 {
            let inside = inside_of(&c)?.unwrap_or(true);
            let k = find(&mut comp, from(c[0]));
            cells.push(Cell {
                cycle: c,
                area: a,
                comp: k,
                inside,
                holes: Vec::new(),
            });
        } else if a < 0.0 {
            outsides.push((c, a));
        } else {
            return Err(Problem::Inconsistent(
                "a cycle of zero area",
                g.cluster[from(c[0])],
            ));
        }
    }
    // The regions around a point are among those whose boxes hold it
    // (padded against the rounding of the containment test).
    let pad = 1e-9 * scale(p);
    let tree = (!outsides.is_empty()).then(|| {
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
    let mut near: Vec<usize> = Vec::new();
    for (c, _) in outsides {
        let k = find(&mut comp, from(c[0]));
        // A cluster of the cycle: farther than the tolerance from every
        // segment it does not end.
        let first = c
            .iter()
            .map(|&h| from(h))
            .find(|&n| n < g.real)
            .expect("a cycle passes through a cluster");
        let q = p[first];
        let mut best: Option<usize> = None;
        if let Some(t) = &tree {
            let q3 = q.extend(0.0);
            t.overlapping(q3 - DVec3::splat(pad), q3 + DVec3::splat(pad), &mut near);
        }
        for &i in &near {
            let r = &cells[i];
            if r.comp == k || !contains(&r.cycle, q, p, &from) {
                continue;
            }
            if best.is_none_or(|b| r.area < cells[b].area) {
                best = Some(i);
            }
        }
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
                });
            }
            loops.push(segs);
            chords.push(cs);
        }
        Region { loops, chords }
    };
    Ok(cells
        .iter()
        .filter(|r| r.inside)
        .map(|r| as_region(std::iter::once(&r.cycle).chain(&r.holes).collect()))
        .collect())
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
/// `eps` of a chord it does not end. Problems are reported in the order of
/// a full pairwise scan (by chord, then crossings before close points, each
/// by index), whether the candidates come from the scan or from the trees.
fn consistent(
    ends: &[(usize, usize)],
    p: &[DVec2],
    real: usize,
    cluster: &[usize],
    eps: f64,
) -> Result<(), Problem> {
    let boxes: Vec<(DVec3, DVec3)> = ends
        .iter()
        .map(|&(a, b)| (p[a].min(p[b]).extend(0.0), p[a].max(p[b]).extend(0.0)))
        .collect();
    // A computed distance within `eps` puts the point within `eps` of the
    // chord's box, up to rounding.
    let r = eps + 1e-9 * (eps + scale(p));
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
            if n != a && n != b && point_segment(p[n], p[a], p[b]) <= eps {
                return Err(Problem::Close(cluster[n]));
            }
        }
    }
    Ok(())
}

/// Removes antennas: a half-edge followed by its twin.
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
