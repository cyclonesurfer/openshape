//! The candidate vertices of a boolean's result and their incidences with
//! the operands' edges (design note 0012 §4.1: one decision per fact).
//!
//! Raw points are the operands' vertices, the points where an edge of one
//! operand crosses a face plane of the other (found from the signs of its
//! ends, decided once per vertex and plane), and the points where two
//! edges touch. Raw points within the tolerance of each other are one
//! *cluster*, a candidate vertex. Whether a cluster lies on an edge is then
//! decided once, by its distance from the edge, and every later stage asks
//! this module rather than testing again.

use std::collections::HashMap;

use glam::DVec3;

use super::{show, slot, Pair};
use crate::origin::Operand;
use crate::KernelError;

/// What a raw point is.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Raw {
    /// A vertex of an operand.
    Vertex(Operand, usize),
    /// An edge of an operand crossing a plane of the other.
    Crossing,
    /// Two edges touching.
    Contact,
}

pub(super) struct Points {
    pub eps: f64,
    /// Per operand, per vertex, per plane of the other operand: the side
    /// of the plane the vertex is on (−1, 0 within the tolerance, +1).
    signs: [Vec<i8>; 2],
    /// The number of planes of the other operand, per operand.
    stride: [usize; 2],
    pos: Vec<DVec3>,
    raw: Vec<Raw>,
    /// Raw index of each operand's first vertex.
    offset: [usize; 2],
    /// (operand, edge, plane of the other operand) → raw crossing point.
    crossings: HashMap<(Operand, usize, usize), usize>,
    /// Raw point → cluster.
    cluster: Vec<usize>,
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
}

impl Points {
    pub fn new(pair: &Pair, eps: f64) -> Points {
        let mut signs = [Vec::new(), Vec::new()];
        let mut stride = [0, 0];
        for op in [Operand::First, Operand::Second] {
            let me = pair.get(op);
            let them = pair.get(super::other(op));
            let n = them.planes.len();
            stride[slot(op)] = n;
            let mut s = Vec::with_capacity(me.points.len() * n);
            for p in &me.points {
                for q in &them.planes {
                    let d = (*p - q.point).dot(q.normal);
                    s.push(if d.abs() <= eps { 0 } else { d.signum() as i8 });
                }
            }
            signs[slot(op)] = s;
        }
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
            signs,
            stride,
            pos,
            raw,
            offset,
            crossings: HashMap::new(),
            cluster: Vec::new(),
            rep: Vec::new(),
            members: Vec::new(),
            has_vertex: Vec::new(),
            on_edges: [Vec::new(), Vec::new()],
        }
    }

    /// The side of `plane` (of the other operand) that vertex `v` of `op`
    /// is on: −1, 0 (within the tolerance) or +1.
    pub fn sign(&self, op: Operand, v: usize, plane: usize) -> i8 {
        self.signs[slot(op)][v * self.stride[slot(op)] + plane]
    }

    /// Records where the edges of face `f` of `op` cross `plane` of the
    /// other operand (their ends strictly on opposite sides).
    pub fn add_crossings(&mut self, pair: &Pair, op: Operand, f: usize, plane: usize) {
        let poly = pair.get(op);
        let q = &pair.get(super::other(op)).planes[plane];
        for l in &poly.faces[f].loops {
            for &(e, _) in l {
                if self.crossings.contains_key(&(op, e, plane)) {
                    continue;
                }
                let (va, vb) = poly.edges[e];
                let (sa, sb) = (self.sign(op, va, plane), self.sign(op, vb, plane));
                if sa * sb >= 0 {
                    continue;
                }
                let (a, b) = (poly.points[va], poly.points[vb]);
                let (da, db) = ((a - q.point).dot(q.normal), (b - q.point).dot(q.normal));
                let p = a + (b - a) * (da / (da - db));
                self.crossings.insert((op, e, plane), self.pos.len());
                self.pos.push(p);
                self.raw.push(Raw::Crossing);
            }
        }
    }

    /// Records where an edge of one operand touches an edge of the other
    /// (closer than the tolerance, not parallel): crossings within a common
    /// plane that no edge/plane crossing finds, such as two seams.
    pub fn add_edge_contacts(&mut self, pair: &Pair) {
        let eps = self.eps;
        let boxes = |op: Operand| -> Vec<(DVec3, DVec3)> {
            let poly = pair.get(op);
            (0..poly.edges.len())
                .map(|e| {
                    let (a, b) = poly.edge_points(e);
                    (a.min(b), a.max(b))
                })
                .collect()
        };
        let (ba, bb) = (boxes(Operand::First), boxes(Operand::Second));
        for (ea, xa) in ba.iter().enumerate() {
            for (eb, xb) in bb.iter().enumerate() {
                if !super::boxes_overlap(*xa, *xb, eps) {
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
    }

    /// Clusters the raw points (union–find within the tolerance) and checks
    /// that no operand has two vertices in one cluster and that no cluster
    /// spreads wider than twice the tolerance.
    pub fn cluster(&mut self) -> Result<(), KernelError> {
        let eps = self.eps;
        let n = self.pos.len();
        let mut parent: Vec<usize> = (0..n).collect();
        fn find(parent: &mut [usize], mut i: usize) -> usize {
            while parent[i] != i {
                parent[i] = parent[parent[i]];
                i = parent[i];
            }
            i
        }
        let h = 2.0 * eps;
        let cell = |p: DVec3| {
            (
                (p.x / h).floor() as i64,
                (p.y / h).floor() as i64,
                (p.z / h).floor() as i64,
            )
        };
        let mut grid: HashMap<(i64, i64, i64), Vec<usize>> = HashMap::new();
        for i in 0..n {
            let c = cell(self.pos[i]);
            for dx in -1..=1 {
                for dy in -1..=1 {
                    for dz in -1..=1 {
                        if let Some(list) = grid.get(&(c.0 + dx, c.1 + dy, c.2 + dz)) {
                            for &j in list {
                                if (self.pos[i] - self.pos[j]).length() <= eps {
                                    let (ri, rj) = (find(&mut parent, i), find(&mut parent, j));
                                    if ri != rj {
                                        parent[ri.max(rj)] = ri.min(rj);
                                    }
                                }
                            }
                        }
                    }
                }
            }
            grid.entry(c).or_default().push(i);
        }
        let mut id = vec![usize::MAX; n];
        self.cluster = vec![0; n];
        for i in 0..n {
            let r = find(&mut parent, i);
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
            for op in [Operand::First, Operand::Second] {
                if count[slot(op)] > 1 {
                    return Err(KernelError::InvalidInput(format!(
                        "the {} body has vertices closer than the tolerance near {}",
                        super::name(op),
                        show(at.unwrap_or(self.pos[m[0]]))
                    )));
                }
            }
            let rep = at
                .unwrap_or_else(|| m.iter().map(|&i| self.pos[i]).sum::<DVec3>() / m.len() as f64);
            if m.iter().any(|&i| (self.pos[i] - rep).length() > 2.0 * eps) {
                return Err(KernelError::InvalidInput(format!(
                    "the bodies meet in points closer than the tolerance near {}",
                    show(rep)
                )));
            }
            self.rep.push(rep);
            self.has_vertex.push([count[0] > 0, count[1] > 0]);
        }
        Ok(())
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

    /// The cluster of the point where edge `e` of `op` crosses `plane` of
    /// the other operand, if recorded.
    pub fn crossing_cluster(&self, op: Operand, e: usize, plane: usize) -> Option<usize> {
        self.crossings
            .get(&(op, e, plane))
            .map(|&i| self.cluster[i])
    }

    /// Decides, once per cluster and edge, which clusters lie on each edge
    /// of both operands (within the tolerance, between its ends), and
    /// orders them along the edge.
    pub fn incidences(&mut self, pair: &Pair) -> Result<(), KernelError> {
        for op in [Operand::First, Operand::Second] {
            let poly = pair.get(op);
            let mut lists = Vec::with_capacity(poly.edges.len());
            // A point on two edges of one body away from its vertices (near
            // a corner, within the tolerance of both edges) belongs to
            // neither more than the other: refused rather than guessed.
            let mut edge_of: HashMap<usize, usize> = HashMap::new();
            for (e, &(va, vb)) in poly.edges.iter().enumerate() {
                let (c0, c1) = (self.vertex_cluster(op, va), self.vertex_cluster(op, vb));
                let (a, b) = poly.edge_points(e);
                let list = self.along(a, b, c0, c1, Some(op))?;
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
                lists.push(list);
            }
            self.on_edges[slot(op)] = lists;
        }
        Ok(())
    }

    /// The clusters on segment `a`–`b` from cluster `c0` (at `a`) to `c1`
    /// (at `b`), skipping clusters holding a vertex of `skip`.
    fn along(
        &self,
        a: DVec3,
        b: DVec3,
        c0: usize,
        c1: usize,
        skip: Option<Operand>,
    ) -> Result<Vec<usize>, KernelError> {
        let eps = self.eps;
        let d = b - a;
        let len2 = d.length_squared();
        let (lo, hi) = (a.min(b) - DVec3::splat(eps), a.max(b) + DVec3::splat(eps));
        let mut found: Vec<(f64, usize)> = Vec::new();
        for (c, p) in self.rep.iter().enumerate() {
            if c == c0 || c == c1 || p.cmplt(lo).any() || p.cmpgt(hi).any() {
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
        out.push(c0);
        let mut last = 0.0;
        for (t, c) in found {
            if (t - last) * len <= eps {
                return Err(KernelError::InvalidInput(format!(
                    "the bodies meet in points closer than the tolerance near {}",
                    show(self.rep[c])
                )));
            }
            out.push(c);
            last = t;
        }
        if (1.0 - last) * len <= eps {
            return Err(KernelError::InvalidInput(format!(
                "the bodies meet in points closer than the tolerance near {}",
                show(b)
            )));
        }
        out.push(c1);
        Ok(out)
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

/// Distance from `p` to the segment `a`–`b`.
pub(super) fn segment_distance(p: DVec3, a: DVec3, b: DVec3) -> f64 {
    let d = b - a;
    let len2 = d.length_squared();
    if len2 == 0.0 {
        return (p - a).length();
    }
    let t = ((p - a).dot(d) / len2).clamp(0.0, 1.0);
    (p - (a + d * t)).length()
}
