//! The candidate vertices of a boolean's result and their incidences with
//! the operands' edges (design note 0012 §4.1: one decision per fact).
//!
//! Raw points are the operands' vertices, the points where an edge of one
//! operand crosses a face plane of the other (found from the signs of its
//! ends, decided once per vertex and plane), and the points where two
//! edges touch. Raw points within the tolerance of each other are one
//! *cluster*, a candidate vertex. Whether a cluster lies on an edge is then
//! decided once, by its distance from the edge, and every later stage asks
//! this module rather than testing again. Two clusters on one edge within
//! the tolerance of each other along it are one point too, and are joined
//! before the edges are split (see [`Points::incidences`]).

use glam::DVec3;
use rustc_hash::FxHashMap;

use super::meet::Track;
use super::quadric::Quadric;
use super::{show, slot, Pair};
use crate::bvh::Bvh;
use crate::geom::Curve3;
use crate::origin::Operand;
use crate::KernelError;

/// What a raw point is.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Raw {
    /// A vertex of an operand.
    Vertex(Operand, usize),
    /// An edge of an operand crossing a surface of the other.
    Crossing,
    /// A station of an edge within the tolerance of a surface of the
    /// other operand, where the edge reaches the surface from one side
    /// ([`Station`]).
    Station,
    /// Two edges touching.
    Contact,
    /// A point put on a closed track that no edge reaches, to split it
    /// into edges with ends.
    Seam,
}

/// A station of an edge against a surface of the other operand (design
/// note 0012 §10.2): one of the edge's ends, or a split of it (a turning
/// point of the surface's implicit function along the edge, or a joint of
/// its rational Bézier pieces). Between consecutive stations the edge
/// crosses the surface at most once, exactly when their signs are strictly
/// opposite.
#[derive(Clone, Copy, Debug, PartialEq)]
pub(super) struct Station {
    /// The edge's parameter.
    pub t: f64,
    /// The side of the surface: −1, 0 (within the tolerance) or +1, for an
    /// end the sign of its vertex.
    pub sign: i8,
    /// True at a turning point, where the edge may touch the surface.
    pub turn: bool,
}

/// An edge of an operand and a surface of the other.
type EdgeSurface = (Operand, usize, usize);

pub(super) struct Points {
    pub eps: f64,
    /// Per operand, its vertices and the other operand's surfaces, for the
    /// sides of the surfaces the vertices are on.
    verts: [Vec<DVec3>; 2],
    surfaces: [Vec<Quadric>; 2],
    pos: Vec<DVec3>,
    raw: Vec<Raw>,
    /// Raw index of each operand's first vertex.
    offset: [usize; 2],
    /// Per edge and surface of the other operand, its stations, from the
    /// edge's start to its end, made when first asked for.
    stations: FxHashMap<EdgeSurface, Vec<Station>>,
    /// (edge and surface, stretch between stations `k` and `k + 1`) → raw
    /// crossing point.
    crossings: FxHashMap<(EdgeSurface, usize), usize>,
    /// (edge and surface, station) → raw point of an interior station on
    /// the surface.
    station_points: FxHashMap<(EdgeSurface, usize), usize>,
    /// Raw point → cluster.
    cluster: Vec<usize>,
    /// The union–find forest over the raw points that the clusters are
    /// the trees of.
    parent: Vec<usize>,
    /// Per cluster, its position: an operand vertex if it has one (the
    /// first operand's first), else the mean of its raw points.
    pub rep: Vec<DVec3>,
    /// Per cluster, its raw points.
    members: Vec<Vec<usize>>,
    /// Per cluster, whether it holds a vertex of the first and the second
    /// operand.
    has_vertex: Vec<[bool; 2]>,
    /// Per operand, per edge: the clusters on it from its start to its end.
    on_edges: [Vec<Vec<usize>>; 2],
    /// Per operand, per curved edge: the parameters of those clusters (the
    /// edge's own at its ends, of their nearest points between); empty for
    /// a straight edge.
    on_params: [Vec<Vec<f64>>; 2],
    /// The clusters' positions, indexed.
    index: Bvh,
}

impl Points {
    pub fn new(pair: &Pair, eps: f64) -> Points {
        let verts = [pair.first.points.clone(), pair.second.points.clone()];
        let surfaces = [
            pair.get(Operand::Second).surfaces.clone(),
            pair.get(Operand::First).surfaces.clone(),
        ];
        let mut pos = Vec::new();
        let mut raw = Vec::new();
        let mut offset = [0, 0];
        for op in [Operand::First, Operand::Second] {
            offset[slot(op)] = pos.len();
            for (v, p) in pair.get(op).points.iter().enumerate() {
                pos.push(*p);
                raw.push(Raw::Vertex(op, v));
            }
        }
        Points {
            eps,
            verts,
            surfaces,
            pos,
            raw,
            offset,
            stations: FxHashMap::default(),
            crossings: FxHashMap::default(),
            station_points: FxHashMap::default(),
            cluster: Vec::new(),
            parent: Vec::new(),
            rep: Vec::new(),
            members: Vec::new(),
            has_vertex: Vec::new(),
            on_edges: [Vec::new(), Vec::new()],
            on_params: [Vec::new(), Vec::new()],
            index: Bvh::new(Vec::new()),
        }
    }

    /// The number of raw points so far.
    pub fn raw_count(&self) -> usize {
        self.pos.len()
    }

    /// The side of `surface` (of the other operand) that vertex `v` of
    /// `op` is on: −1, 0 (within the tolerance) or +1, by its signed
    /// distance. A function of the vertex and the surface alone, computed
    /// the same way every time it is asked, so every stage sees one answer
    /// (the decision is made once in effect, without a table of every
    /// vertex against every surface).
    pub fn sign(&self, op: Operand, v: usize, surface: usize) -> i8 {
        let d = self.surfaces[slot(op)][surface].distance(self.verts[slot(op)][v]);
        if d.abs() <= self.eps {
            0
        } else {
            d.signum() as i8
        }
    }

    /// The side of `surface` (of the other operand) at distance `d`.
    fn side(&self, d: f64) -> i8 {
        if d.abs() <= self.eps {
            0
        } else {
            d.signum() as i8
        }
    }

    /// The stations of edge `e` of `op` against `surface` of the other
    /// operand, from the edge's start to its end, each decided once. A
    /// straight edge has none but its ends against a plane (the distance is
    /// linear along it).
    pub fn stations(
        &mut self,
        pair: &Pair,
        op: Operand,
        e: usize,
        surface: usize,
    ) -> Result<&[Station], KernelError> {
        let key = (op, e, surface);
        if !self.stations.contains_key(&key) {
            let list = self.station_list(pair, op, e, surface)?;
            self.stations.insert(key, list);
        }
        Ok(&self.stations[&key])
    }

    /// The stations of edge `e` of `op` against `surface`, as
    /// [`Points::stations`] finds them.
    fn station_list(
        &self,
        pair: &Pair,
        op: Operand,
        e: usize,
        surface: usize,
    ) -> Result<Vec<Station>, KernelError> {
        let shape = pair.get(op);
        let q = self.surfaces[slot(op)][surface];
        let (va, vb) = shape.edges[e];
        let (curve, t0, t1) = shape.curve(e);
        let mut list = vec![Station {
            t: t0,
            sign: self.sign(op, va, surface),
            turn: false,
        }];
        if !(shape.straight[e] && matches!(q, Quadric::Plane { .. })) {
            for split in q.splits(curve, t0, t1)? {
                list.push(Station {
                    t: split.t,
                    sign: self.side(q.distance(curve.point(split.t))),
                    turn: split.turn,
                });
            }
        }
        list.push(Station {
            t: t1,
            sign: self.sign(op, vb, surface),
            turn: false,
        });
        Ok(list)
    }

    /// The side of `surface` (of the other operand) that edge `e` of `op`
    /// leaves its start (or, if `!from_start`, its end) to: the sign of its
    /// next station (the implicit function is monotone between stations),
    /// 0 if that is on the surface.
    pub fn leaving(
        &self,
        pair: &Pair,
        op: Operand,
        e: usize,
        surface: usize,
        from_start: bool,
    ) -> Result<i8, KernelError> {
        let found;
        let st = match self.cached_stations(op, e, surface) {
            Some(st) => st,
            None => {
                found = self.station_list(pair, op, e, surface)?;
                &found
            }
        };
        Ok(if from_start {
            st[1].sign
        } else {
            st[st.len() - 2].sign
        })
    }

    /// The stations of edge `e` of `op` against `surface`, if made.
    pub fn cached_stations(&self, op: Operand, e: usize, surface: usize) -> Option<&[Station]> {
        self.stations.get(&(op, e, surface)).map(Vec::as_slice)
    }

    /// Records where the edges of face `f` of `op` cross `surface` of the
    /// other operand (between stations of strictly opposite signs), and
    /// the stations inside them on the surface that the edges reach it at
    /// from one side.
    pub fn add_crossings(
        &mut self,
        pair: &Pair,
        op: Operand,
        f: usize,
        surface: usize,
    ) -> Result<(), KernelError> {
        let shape = pair.get(op);
        let face = &shape.faces[f];
        let q = pair.get(super::other(op)).surfaces[surface];
        let mut near = Vec::new();
        face.near_surface(&q, self.eps, &mut near);
        for k in near {
            let (li, i) = face.coedge_at(k);
            let e = face.loops[li][i].0;
            if shape.straight[e] && matches!(q, Quadric::Plane { .. }) {
                // Step 1's crossing of a line and a plane, where the
                // distance is linear.
                let (va, vb) = shape.edges[e];
                let (sa, sb) = (self.sign(op, va, surface), self.sign(op, vb, surface));
                if sa * sb >= 0 || self.crossings.contains_key(&((op, e, surface), 0)) {
                    continue;
                }
                let (a, b) = (shape.points[va], shape.points[vb]);
                let (da, db) = (q.distance(a), q.distance(b));
                let p = a + (b - a) * (da / (da - db));
                self.crossings.insert(((op, e, surface), 0), self.pos.len());
                self.pos.push(p);
                self.raw.push(Raw::Crossing);
                continue;
            }
            let st = self.stations(pair, op, e, surface)?.to_vec();
            let (curve, _, _) = shape.curve(e);
            for k in 0..st.len() - 1 {
                let key = ((op, e, surface), k);
                if st[k].sign * st[k + 1].sign >= 0 || self.crossings.contains_key(&key) {
                    continue;
                }
                let p = curve.point(q.crossing(curve, st[k].t, st[k + 1].t));
                self.crossings.insert(key, self.pos.len());
                self.pos.push(p);
                self.raw.push(Raw::Crossing);
            }
            for i in 1..st.len() - 1 {
                let key = ((op, e, surface), i);
                let reached = st[i - 1].sign != 0 || st[i + 1].sign != 0;
                if st[i].sign != 0 || !reached || self.station_points.contains_key(&key) {
                    continue;
                }
                self.station_points.insert(key, self.pos.len());
                self.pos.push(curve.point(st[i].t));
                self.raw.push(Raw::Station);
            }
        }
        Ok(())
    }

    /// Records where an edge of one operand touches an edge of the other
    /// (closer than the tolerance, not parallel): crossings within a common
    /// plane that no edge/plane crossing finds, such as two seams. Fails if
    /// more than `limit` pairs of edges come within the tolerance's reach of
    /// each other.
    pub fn add_edge_contacts(&mut self, pair: &Pair, limit: usize) -> Result<(), usize> {
        let eps = self.eps;
        let boxes = |op: Operand| -> Vec<(DVec3, DVec3)> {
            let shape = pair.get(op);
            (0..shape.edges.len())
                .map(|e| {
                    let (a, b) = shape.edge_points(e);
                    (a.min(b), a.max(b))
                })
                .collect()
        };
        let ba = boxes(Operand::First);
        let bb = boxes(Operand::Second);
        let tree = Bvh::new(bb.clone());
        let pad = super::pad(eps);
        let mut near = Vec::new();
        let mut examined = 0;
        for (ea, xa) in ba.iter().enumerate() {
            tree.overlapping(xa.0 - pad, xa.1 + pad, &mut near);
            examined += near.len();
            if examined > limit {
                return Err(limit);
            }
            if !pair.first.straight[ea] {
                continue;
            }
            for &eb in &near {
                let xb = &bb[eb];
                if !pair.second.straight[eb] || !super::boxes_overlap(*xa, *xb, eps) {
                    continue;
                }
                let (p0, p1) = pair.first.edge_points(ea);
                let (q0, q1) = pair.second.edge_points(eb);
                let (u, v) = (p1 - p0, q1 - q0);
                if u.cross(v).length() <= 1e-9 * u.length() * v.length() {
                    continue;
                }
                let (s, t) = closest_parameters(p0, u, q0, v);
                let (x, y) = (p0 + u * s, q0 + v * t);
                if (x - y).length() <= eps {
                    self.pos.push((x + y) * 0.5);
                    self.raw.push(Raw::Contact);
                }
            }
        }
        Ok(())
    }

    /// Clusters the raw points (union–find within the tolerance) and checks
    /// that no operand has two vertices in one cluster and that no cluster
    /// spreads wider than twice the tolerance.
    pub fn cluster(&mut self) -> Result<(), KernelError> {
        let eps = self.eps;
        let n = self.pos.len();
        let mut parent: Vec<usize> = (0..n).collect();
        let h = 2.0 * eps;
        let cell = |p: DVec3| {
            (
                (p.x / h).floor() as i64,
                (p.y / h).floor() as i64,
                (p.z / h).floor() as i64,
            )
        };
        let mut grid: FxHashMap<(i64, i64, i64), Vec<usize>> = FxHashMap::default();
        for i in 0..n {
            let c = cell(self.pos[i]);
            for dx in -1..=1 {
                for dy in -1..=1 {
                    for dz in -1..=1 {
                        if let Some(list) = grid.get(&(c.0 + dx, c.1 + dy, c.2 + dz)) {
                            for &j in list {
                                if (self.pos[i] - self.pos[j]).length() <= eps {
                                    union(&mut parent, i, j);
                                }
                            }
                        }
                    }
                }
            }
            grid.entry(c).or_default().push(i);
        }
        self.parent = parent;
        self.settle(false)
    }

    /// Makes the clusters from the union–find forest: their members,
    /// places and vertices, and the index of their places. `joined` says
    /// that clusters were joined along edges since the first clustering:
    /// two vertices of one body then met through points between them, not
    /// by being within the tolerance of each other.
    fn settle(&mut self, joined: bool) -> Result<(), KernelError> {
        let eps = self.eps;
        let n = self.pos.len();
        let mut id = vec![usize::MAX; n];
        self.cluster = vec![0; n];
        self.members.clear();
        self.rep.clear();
        self.has_vertex.clear();
        for i in 0..n {
            let r = find(&mut self.parent, i);
            if id[r] == usize::MAX {
                id[r] = self.members.len();
                self.members.push(Vec::new());
            }
            self.cluster[i] = id[r];
            self.members[id[r]].push(i);
        }
        for m in &self.members {
            let mut count = [0usize; 2];
            let mut at = None;
            for &i in m {
                if let Raw::Vertex(op, _) = self.raw[i] {
                    count[slot(op)] += 1;
                    if at.is_none() {
                        at = Some(self.pos[i]);
                    }
                }
            }
            let rep = at
                .unwrap_or_else(|| m.iter().map(|&i| self.pos[i]).sum::<DVec3>() / m.len() as f64);
            for op in [Operand::First, Operand::Second] {
                if count[slot(op)] > 1 {
                    return Err(KernelError::InvalidInput(if joined {
                        format!(
                            "the bodies meet in points closer than the tolerance near {}",
                            show(rep)
                        )
                    } else {
                        format!(
                            "the {} body has vertices closer than the tolerance near {}",
                            super::name(op),
                            show(rep)
                        )
                    }));
                }
            }
            if m.iter().any(|&i| (self.pos[i] - rep).length() > 2.0 * eps) {
                return Err(KernelError::InvalidInput(format!(
                    "the bodies meet in points closer than the tolerance near {}",
                    show(rep)
                )));
            }
            self.rep.push(rep);
            self.has_vertex.push([count[0] > 0, count[1] > 0]);
        }
        self.index = Bvh::of_points(&self.rep);
        Ok(())
    }

    /// New clusters at `ps`, after the clustering, each of one point: the
    /// seams of closed tracks that no edge reaches, which the caller keeps
    /// clear of every other point. Returns the first new cluster; the rest
    /// follow it in order.
    pub fn add_seams(&mut self, ps: &[DVec3]) -> usize {
        let first = self.rep.len();
        for &p in ps {
            let i = self.pos.len();
            self.pos.push(p);
            self.raw.push(Raw::Seam);
            self.parent.push(i);
            self.cluster.push(self.rep.len());
            self.members.push(vec![i]);
            self.rep.push(p);
            self.has_vertex.push([false, false]);
        }
        if !ps.is_empty() {
            self.index = Bvh::of_points(&self.rep);
        }
        first
    }

    /// The cluster of vertex `v` of `op`.
    pub fn vertex_cluster(&self, op: Operand, v: usize) -> usize {
        self.cluster[self.offset[slot(op)] + v]
    }

    /// The vertices of the first and the second operand in cluster `c`.
    pub fn vertices_of(&self, c: usize) -> (Option<usize>, Option<usize>) {
        let mut out = (None, None);
        for &i in &self.members[c] {
            match self.raw[i] {
                Raw::Vertex(Operand::First, v) => out.0 = Some(v),
                Raw::Vertex(Operand::Second, v) => out.1 = Some(v),
                _ => {}
            }
        }
        out
    }

    /// The cluster of the point where edge `e` of `op` crosses `surface` of
    /// the other operand between its stations `k` and `k + 1`, if recorded.
    pub fn crossing_cluster(
        &self,
        op: Operand,
        e: usize,
        surface: usize,
        k: usize,
    ) -> Option<usize> {
        self.crossings
            .get(&((op, e, surface), k))
            .map(|&i| self.cluster[i])
    }

    /// The cluster of station `i` of edge `e` of `op` against `surface`:
    /// its vertex's at an end, else its own if recorded.
    pub fn station_cluster(
        &self,
        pair: &Pair,
        op: Operand,
        e: usize,
        surface: usize,
        i: usize,
    ) -> Option<usize> {
        let n = self.stations.get(&(op, e, surface))?.len();
        let (va, vb) = pair.get(op).edges[e];
        if i == 0 {
            Some(self.vertex_cluster(op, va))
        } else if i + 1 == n {
            Some(self.vertex_cluster(op, vb))
        } else {
            self.station_points
                .get(&((op, e, surface), i))
                .map(|&r| self.cluster[r])
        }
    }

    /// Decides, once per cluster and edge, which clusters lie on each edge
    /// of both operands (within the tolerance, between its ends), and
    /// orders them along the edge.
    ///
    /// First, two clusters of derived points (crossings and contacts, no
    /// vertex of either body) on one edge within the tolerance of each
    /// other along it are joined into one: they are one point of the
    /// result. Where an edge meets a plane at a shallow angle, the
    /// tolerance across the plane is a long stretch along the edge, so a
    /// crossing found from each side (the edge of one body crossing a
    /// plane of the other, and the other's edge crossing the first's
    /// plane) can land just over the tolerance apart in space while within
    /// it along the edge; refusing them refused copies of a body moved by
    /// a fraction of the tolerance. The bodies' vertices are their own
    /// data, placed exactly: two of them are one point only within the
    /// tolerance, and a cluster holding one is not joined (joining copies
    /// of a body's vertices moved by just over the tolerance made an
    /// inconsistent result instead of a refusal). Joining moves a
    /// cluster's place, so it is repeated until nothing is joined (a few
    /// rounds); whatever is still that close after them is refused below,
    /// as is a joined cluster spreading wider than twice the tolerance.
    pub fn incidences(&mut self, pair: &Pair) -> Result<(), KernelError> {
        // Each round finds the clusters on every edge once; the last
        // round's lists, made after the last join, are the edges' lists.
        let mut rounds = 0;
        let lists = loop {
            let lists = [Operand::First, Operand::Second].map(|op| {
                let shape = pair.get(op);
                shape
                    .edges
                    .iter()
                    .enumerate()
                    .map(|(e, &(va, vb))| {
                        let (c0, c1) = (self.vertex_cluster(op, va), self.vertex_cluster(op, vb));
                        if shape.straight[e] {
                            let (a, b) = shape.edge_points(e);
                            self.on_segment(a, b, c0, c1, Some(op))
                        } else {
                            self.on_curve(pair, op, e, c0, c1)
                        }
                    })
                    .collect::<Vec<_>>()
            });
            let derived = |c: usize| self.has_vertex[c] == [false, false];
            let joins: Vec<(usize, usize)> = lists
                .iter()
                .flatten()
                .flat_map(|on| on.windows(2))
                .filter(|w| w[1].0 - w[0].0 <= self.eps && derived(w[0].1) && derived(w[1].1))
                .map(|w| (w[0].1, w[1].1))
                .collect();
            if joins.is_empty() || rounds == JOIN_ROUNDS {
                break lists;
            }
            for (c, d) in joins {
                let (i, j) = (self.members[c][0], self.members[d][0]);
                union(&mut self.parent, i, j);
            }
            self.settle(true)?;
            rounds += 1;
        };
        for (op, ons) in [Operand::First, Operand::Second].into_iter().zip(lists) {
            let shape = pair.get(op);
            // A point on two edges of one body away from its vertices (near
            // a corner, within the tolerance of both edges) belongs to
            // neither more than the other: refused rather than guessed.
            let mut edge_of: FxHashMap<usize, usize> = FxHashMap::default();
            let mut out = Vec::with_capacity(ons.len());
            let mut params = Vec::with_capacity(ons.len());
            for (e, on) in ons.into_iter().enumerate() {
                let (_, b) = shape.edge_points(e);
                let list = self.spaced(on, b)?;
                for &c in &list[1..list.len() - 1] {
                    if edge_of.insert(c, e).is_some() {
                        return Err(KernelError::InvalidInput(format!(
                            "a point where the bodies meet lies within the tolerance of two \
                             edges of the {} body near {}",
                            super::name(op),
                            show(self.rep[c])
                        )));
                    }
                }
                params.push(if shape.straight[e] {
                    Vec::new()
                } else {
                    let (curve, t0, t1) = shape.curve(e);
                    let n = list.len();
                    (0..n)
                        .map(|k| match k {
                            0 => t0,
                            _ if k + 1 == n => t1,
                            _ => foot_on(curve, t0, t1, self.rep[list[k]]).0,
                        })
                        .collect()
                });
                out.push(list);
            }
            self.on_edges[slot(op)] = out;
            self.on_params[slot(op)] = params;
        }
        Ok(())
    }

    /// The clusters on curved edge `e` of `op` from cluster `c0` (at its
    /// start) to `c1` (at its end), skipping clusters holding a vertex of
    /// `op`: those within the tolerance of a point of the edge strictly
    /// between its ends, nearest in turn, each with its distance from the
    /// start along the chords between them.
    fn on_curve(
        &self,
        pair: &Pair,
        op: Operand,
        e: usize,
        c0: usize,
        c1: usize,
    ) -> Vec<(f64, usize)> {
        let eps = self.eps;
        let shape = pair.get(op);
        let (curve, t0, t1) = shape.curve(e);
        let (lo, hi) = shape.edge_boxes[e];
        let mut near = Vec::new();
        self.index
            .overlapping(lo - DVec3::splat(eps), hi + DVec3::splat(eps), &mut near);
        let mut found: Vec<(f64, usize)> = Vec::new();
        for c in near {
            if c == c0 || c == c1 || self.has_vertex[c][slot(op)] {
                continue;
            }
            let (t, d) = foot_on(curve, t0, t1, self.rep[c]);
            if d <= eps && t > t0 && t < t1 {
                found.push((t, c));
            }
        }
        found.sort_by(|x, y| x.0.total_cmp(&y.0).then(x.1.cmp(&y.1)));
        let mut out = Vec::with_capacity(found.len() + 2);
        out.push((0.0, c0));
        let (mut at, mut prev) = (0.0, curve.point(t0));
        for (t, c) in found.into_iter().chain([(t1, c1)]) {
            let p = curve.point(t);
            at += (p - prev).length();
            prev = p;
            out.push((at, c));
        }
        out
    }

    /// The parameters of the clusters on curved edge `e` of `op`
    /// ([`Points::on_edge`]).
    pub fn on_edge_params(&self, op: Operand, e: usize) -> &[f64] {
        &self.on_params[slot(op)][e]
    }

    /// The clusters on the curved track `track` strictly between parameters
    /// `ta` (cluster `lo`) and `tb > ta` (cluster `hi`), less than a turn
    /// apart on a closed track: those within the tolerance of it there, in
    /// order, each with its parameter (unrolled from `ta`), both ends
    /// included. Two within the tolerance of each other along it are
    /// refused.
    pub fn on_track(
        &self,
        track: &Track,
        (lo, ta): (usize, f64),
        (hi, tb): (usize, f64),
    ) -> Result<Vec<(usize, f64)>, KernelError> {
        let eps = self.eps;
        // A box holding the track's stretch.
        let (box_lo, box_hi) = track.bounds(ta, tb);
        let mut near = Vec::new();
        self.index.overlapping(
            box_lo - DVec3::splat(eps),
            box_hi + DVec3::splat(eps),
            &mut near,
        );
        let mut found: Vec<(f64, usize)> = Vec::new();
        for c in near {
            if c == lo || c == hi {
                continue;
            }
            let (t, d) = track.nearest(self.rep[c]);
            if d > eps {
                continue;
            }
            let u = if track.closed() {
                ta + (t - ta).rem_euclid(std::f64::consts::TAU)
            } else {
                t
            };
            if u > ta && u < tb {
                found.push((u, c));
            }
        }
        found.sort_by(|x, y| x.0.total_cmp(&y.0).then(x.1.cmp(&y.1)));
        let mut out = Vec::with_capacity(found.len() + 2);
        out.push((lo, ta));
        out.extend(found.into_iter().map(|(u, c)| (c, u)));
        out.push((hi, tb));
        for w in out.windows(2) {
            if (track.point(w[1].1) - track.point(w[0].1)).length() <= eps {
                return Err(KernelError::InvalidInput(format!(
                    "the bodies meet in points closer than the tolerance near {}",
                    show(self.rep[w[1].0])
                )));
            }
        }
        Ok(out)
    }

    /// The clusters on segment `a`–`b` from cluster `c0` (at `a`) to `c1`
    /// (at `b`), skipping clusters holding a vertex of `skip`, each with its
    /// distance from `a` along the segment, in order.
    fn on_segment(
        &self,
        a: DVec3,
        b: DVec3,
        c0: usize,
        c1: usize,
        skip: Option<Operand>,
    ) -> Vec<(f64, usize)> {
        let eps = self.eps;
        let d = b - a;
        let len2 = d.length_squared();
        let (lo, hi) = (a.min(b) - DVec3::splat(eps), a.max(b) + DVec3::splat(eps));
        let mut found: Vec<(f64, usize)> = Vec::new();
        let mut near = Vec::new();
        self.index.overlapping(lo, hi, &mut near);
        for c in near {
            let p = &self.rep[c];
            if c == c0 || c == c1 {
                continue;
            }
            if let Some(op) = skip {
                if self.has_vertex[c][slot(op)] {
                    continue;
                }
            }
            let t = (*p - a).dot(d) / len2;
            if t <= 0.0 || t >= 1.0 {
                continue;
            }
            if (*p - (a + d * t)).length() <= eps {
                found.push((t, c));
            }
        }
        found.sort_by(|x, y| x.0.total_cmp(&y.0).then(x.1.cmp(&y.1)));
        let len = len2.sqrt();
        let mut out = Vec::with_capacity(found.len() + 2);
        out.push((0.0, c0));
        out.extend(found.into_iter().map(|(t, c)| (t * len, c)));
        out.push((len, c1));
        out
    }

    /// The clusters on segment `a`–`b` from cluster `c0` (at `a`) to `c1`
    /// (at `b`), skipping clusters holding a vertex of `skip`; two of them
    /// within the tolerance of each other along it are refused.
    fn along(
        &self,
        a: DVec3,
        b: DVec3,
        c0: usize,
        c1: usize,
        skip: Option<Operand>,
    ) -> Result<Vec<usize>, KernelError> {
        self.spaced(self.on_segment(a, b, c0, c1, skip), b)
    }

    /// The clusters of a list from [`Points::on_segment`] of a segment
    /// ending at `b`; two within the tolerance of each other along it are
    /// refused.
    fn spaced(&self, on: Vec<(f64, usize)>, b: DVec3) -> Result<Vec<usize>, KernelError> {
        for w in on.windows(2) {
            if w[1].0 - w[0].0 <= self.eps {
                // At the far end, name the end.
                let last = on.last().expect("a segment's list holds its ends").1;
                let at = if w[1].1 == last { b } else { self.rep[w[1].1] };
                return Err(KernelError::InvalidInput(format!(
                    "the bodies meet in points closer than the tolerance near {}",
                    show(at)
                )));
            }
        }
        Ok(on.into_iter().map(|(_, c)| c).collect())
    }

    /// The clusters on edge `e` of `op`, from its start vertex to its end.
    pub fn on_edge(&self, op: Operand, e: usize) -> &[usize] {
        &self.on_edges[slot(op)][e]
    }

    /// The clusters on the segment between clusters `lo` and `hi`, in
    /// order, both ends included.
    pub fn split_segment(&self, lo: usize, hi: usize) -> Result<Vec<usize>, KernelError> {
        self.along(self.rep[lo], self.rep[hi], lo, hi, None)
    }
}

/// Rounds of joining clusters close along an edge
/// ([`Points::incidences`]): each joins every such pair at once, and a
/// second round is needed only where a join moved a cluster close to
/// another.
const JOIN_ROUNDS: usize = 4;

/// The root of `i` in a union–find forest, halving the path.
fn find(parent: &mut [usize], mut i: usize) -> usize {
    while parent[i] != i {
        parent[i] = parent[parent[i]];
        i = parent[i];
    }
    i
}

/// Joins the trees of `i` and `j`, under the lesser root (so the forest
/// does not depend on the order of the joins).
fn union(parent: &mut [usize], i: usize, j: usize) {
    let (ri, rj) = (find(parent, i), find(parent, j));
    if ri != rj {
        parent[ri.max(rj)] = ri.min(rj);
    }
}

/// The parameters `(s, t)` of the closest points of the lines `p + s u`
/// and `q + t v` (not parallel), clamped to the segments `[0, 1]`.
fn closest_parameters(p: DVec3, u: DVec3, q: DVec3, v: DVec3) -> (f64, f64) {
    let w = p - q;
    let (a, b, c) = (u.dot(u), u.dot(v), v.dot(v));
    let (d, e) = (u.dot(w), v.dot(w));
    let den = a * c - b * b;
    let mut s = ((b * e - c * d) / den).clamp(0.0, 1.0);
    let mut t = ((b * s + e) / c).clamp(0.0, 1.0);
    s = ((b * t - d) / a).clamp(0.0, 1.0);
    t = ((b * s + e) / c).clamp(0.0, 1.0);
    (s, t)
}

/// The parameter of the point of the stretch of `curve` between `t0` and
/// `t1` (either order) nearest `p`, and its distance from `p`: the foot of
/// `p` on the curve where that is in the stretch (for a circle, its angle
/// unrolled from the stretch's start), else the nearer end; a NURBS curve
/// is searched within the stretch only (design note 0012 §16.4), so a
/// nearer point of the curve outside it is never taken.
pub(super) fn foot_on(curve: &Curve3, t0: f64, t1: f64, p: DVec3) -> (f64, f64) {
    let (lo, hi) = if t0 <= t1 { (t0, t1) } else { (t1, t0) };
    let t = match curve {
        Curve3::Line { origin, dir } => {
            ((p - *origin).dot(*dir) / dir.length_squared()).clamp(lo, hi)
        }
        Curve3::Circle { center, x, y, .. } => {
            let d = p - *center;
            let a = d.dot(*y).atan2(d.dot(*x));
            let t = lo + (a - lo).rem_euclid(std::f64::consts::TAU);
            if t <= hi {
                t
            } else if (curve.point(lo) - p).length() <= (curve.point(hi) - p).length() {
                lo
            } else {
                hi
            }
        }
        Curve3::Nurbs(n) => n.closest_point_in(p, lo, hi).0,
    };
    (t, (curve.point(t) - p).length())
}

/// Distance from `p` to the segment `a`–`b`.
pub(super) fn segment_distance(p: DVec3, a: DVec3, b: DVec3) -> f64 {
    (p - segment_closest(p, a, b)).length()
}

/// The point of the segment `a`–`b` nearest `p`.
pub(super) fn segment_closest(p: DVec3, a: DVec3, b: DVec3) -> DVec3 {
    let d = b - a;
    let len2 = d.length_squared();
    if len2 == 0.0 {
        return a;
    }
    let t = ((p - a).dot(d) / len2).clamp(0.0, 1.0);
    a + d * t
}
