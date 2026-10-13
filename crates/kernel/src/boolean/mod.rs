//! Booleans of solids bounded by planes: design note 0012 §8 step 1 and
//! ADR 0020.
//!
//! The stages of design note 0012 §4, for polyhedra:
//!
//! 1. **Points** (`points`): every vertex of either body, every point where
//!    an edge of one crosses a face plane of the other, and every point
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
//!    or outside it by the winding number at a point well inside the piece.
//! 5. **Select, sew, merge** (`assemble`): pieces are kept by the table of
//!    §4.4, sewn along their segments, coplanar neighbours and collinear
//!    edges joined, shells and solids rebuilt, and the origins of ADR 0020
//!    recorded. The kernel's validator then checks the result.
//!
//! Anything the stages cannot decide within the tolerance is refused with
//! a message and a location rather than guessed.

use std::collections::HashMap;

use glam::{DVec2, DVec3};

use crate::geom::Curve3;
use crate::math::Frame;
use crate::origin::{Operand, Origins};
use crate::topo::Brep;
use crate::{BooleanOp, KernelError, Tolerance};

mod assemble;
mod classify;
mod planar;
mod points;

use classify::Class;
use points::Points;

/// A face of an operand: a polygon, possibly with holes, in its plane.
struct PolyFace {
    /// The plane's frame: origin on the plane, `x × y` the outward normal.
    frame: Frame,
    /// The face's plane group in its operand (faces on one plane share it).
    plane: usize,
    /// Per loop (outer first), the coedges as (edge, forward).
    loops: Vec<Vec<(usize, bool)>>,
    /// Per loop, the vertex each coedge starts at.
    verts: Vec<Vec<usize>>,
    /// Bounding box of the face's vertices.
    lo: DVec3,
    hi: DVec3,
}

impl PolyFace {
    fn normal(&self) -> DVec3 {
        self.frame.z
    }

    fn to2(&self, p: DVec3) -> DVec2 {
        let d = p - self.frame.origin;
        DVec2::new(d.dot(self.frame.x), d.dot(self.frame.y))
    }

    fn to3(&self, q: DVec2) -> DVec3 {
        self.frame.point2(q)
    }
}

/// A plane shared by one or more faces of an operand, with the outward
/// normal of the first of them.
struct Plane {
    point: DVec3,
    normal: DVec3,
}

/// An operand bounded by planes, with straight edges.
struct Poly<'a> {
    brep: &'a Brep,
    faces: Vec<PolyFace>,
    planes: Vec<Plane>,
    /// Per edge, its start and end vertex.
    edges: Vec<(usize, usize)>,
    /// Vertex positions.
    points: Vec<DVec3>,
}

impl Poly<'_> {
    fn edge_points(&self, e: usize) -> (DVec3, DVec3) {
        let (a, b) = self.edges[e];
        (self.points[a], self.points[b])
    }
}

/// Relative tolerance for recognising faces built on one plane (rounding
/// only), as the tessellator's seam test does.
const SAME_PLANE: f64 = 1e-10;

fn poly<'a>(b: &'a Brep, which: &str) -> Result<Poly<'a>, KernelError> {
    if b.edges
        .iter()
        .any(|e| !matches!(e.curve, Curve3::Line { .. }))
    {
        return Err(KernelError::Unsupported(format!(
            "booleans of bodies with curved edges (the {which} body has some; \
             design note 0012 §8 step 2)"
        )));
    }
    let points: Vec<DVec3> = b.vertices.iter().map(|v| v.point).collect();
    let edges: Vec<(usize, usize)> = b.edges.iter().map(|e| (e.v0.idx(), e.v1.idx())).collect();
    let mut planes: Vec<Plane> = Vec::new();
    let mut faces = Vec::with_capacity(b.faces.len());
    for f in &b.faces {
        let (p, n) = f.surface.plane().ok_or_else(|| {
            KernelError::Unsupported(format!(
                "booleans of bodies with curved faces (the {which} body has some; \
                 design note 0012 §8 step 2)"
            ))
        })?;
        let n = if f.forward { n } else { -n };
        let frame = Frame::from_normal(p, n, n.any_orthonormal_vector());
        let scale = 1.0 + p.length();
        let plane = match planes.iter().position(|q| {
            q.normal.dot(n).abs() >= 1.0 - SAME_PLANE
                && q.normal.dot(p - q.point).abs() <= SAME_PLANE * scale
        }) {
            Some(i) => i,
            None => {
                planes.push(Plane {
                    point: p,
                    normal: n,
                });
                planes.len() - 1
            }
        };
        let mut loops = Vec::with_capacity(f.loops.len());
        let mut verts = Vec::with_capacity(f.loops.len());
        let (mut lo, mut hi) = (DVec3::splat(f64::MAX), DVec3::splat(f64::MIN));
        for l in &f.loops {
            let mut cs = Vec::with_capacity(l.coedges.len());
            let mut vs = Vec::with_capacity(l.coedges.len());
            for c in &l.coedges {
                let v = b.coedge_start(*c).idx();
                cs.push((c.edge.idx(), c.forward));
                vs.push(v);
                lo = lo.min(points[v]);
                hi = hi.max(points[v]);
            }
            loops.push(cs);
            verts.push(vs);
        }
        faces.push(PolyFace {
            frame,
            plane,
            loops,
            verts,
            lo,
            hi,
        });
    }
    Ok(Poly {
        brep: b,
        faces,
        planes,
        edges,
        points,
    })
}

/// The two operands, indexed by [`Operand`].
struct Pair<'a> {
    first: Poly<'a>,
    second: Poly<'a>,
}

impl<'a> Pair<'a> {
    fn get(&self, op: Operand) -> &Poly<'a> {
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

/// How two faces of different operands lie.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Relation {
    /// One face lies wholly on one side of the other's plane.
    Apart,
    /// On one plane (either orientation).
    Coplanar,
    /// The faces' planes cross, and each face reaches the other's plane.
    Crossing,
}

/// A segment of the result's candidate edges: two clusters, the lesser
/// first.
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

/// A segment in one face's graph: material on the left of `a → b`
/// (`Some(true)`), of `b → a` (`Some(false)`), or on both sides (`None`,
/// an imprint inside the face).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct FaceSeg {
    key: Key,
    side: Option<bool>,
}

/// A piece of an operand face, as cycles of clusters (outer first),
/// counter-clockwise about the face's outward normal for the outer loop.
struct Piece {
    operand: Operand,
    face: usize,
    loops: Vec<Vec<usize>>,
    class: Class,
    /// For a piece on a coincident face of the other operand, that face.
    partner: Option<usize>,
}

pub(crate) fn boolean(
    a: &Brep,
    b: &Brep,
    op: BooleanOp,
    tol: &Tolerance,
) -> Result<(Brep, Origins), KernelError> {
    let pair = Pair {
        first: poly(a, "first")?,
        second: poly(b, "second")?,
    };
    let eps = tol.linear;
    let mut pts = Points::new(&pair, eps);

    // Face pairs that can meet, and how.
    let mut crossing: Vec<(usize, usize)> = Vec::new();
    let mut coplanar: Vec<(usize, usize)> = Vec::new();
    for (fa, a) in pair.first.faces.iter().enumerate() {
        for (fb, b) in pair.second.faces.iter().enumerate() {
            if !boxes_overlap((a.lo, a.hi), (b.lo, b.hi), eps) {
                continue;
            }
            match relation(&pair, &pts, fa, fb)? {
                Relation::Apart => {}
                Relation::Coplanar => coplanar.push((fa, fb)),
                Relation::Crossing => {
                    pts.add_crossings(&pair, Operand::First, fa, b.plane);
                    pts.add_crossings(&pair, Operand::Second, fb, a.plane);
                    crossing.push((fa, fb));
                }
            }
        }
    }
    pts.add_edge_contacts(&pair);
    pts.cluster()?;
    pts.incidences(&pair)?;

    // Each face's segments: its own edges' pieces, then imprints.
    let mut segs = Segments::default();
    let mut face_segs: [Vec<Vec<FaceSeg>>; 2] = [
        vec![Vec::new(); pair.first.faces.len()],
        vec![Vec::new(); pair.second.faces.len()],
    ];
    for op in [Operand::First, Operand::Second] {
        let poly = pair.get(op);
        for (fi, f) in poly.faces.iter().enumerate() {
            for l in &f.loops {
                for &(e, forward) in l {
                    let along = pts.on_edge(op, e);
                    for w in along.windows(2) {
                        let (from, to) = if forward { (w[0], w[1]) } else { (w[1], w[0]) };
                        let k = segs.add(from, to, Source::Edge(op, e));
                        face_segs[slot(op)][fi].push(FaceSeg {
                            key: k,
                            side: Some(from == k.0),
                        });
                    }
                }
            }
        }
    }
    for &(fa, fb) in &crossing {
        for (lo, hi) in meeting(&pair, &pts, fa, fb)? {
            for w in pts.split_segment(lo, hi)?.windows(2) {
                let k = segs.add(w[0], w[1], Source::Meet(fa, fb));
                face_segs[0][fa].push(FaceSeg { key: k, side: None });
                face_segs[1][fb].push(FaceSeg { key: k, side: None });
            }
        }
    }
    for &(fa, fb) in &coplanar {
        for (op, f, g) in [(Operand::First, fa, fb), (Operand::Second, fb, fa)] {
            // The pieces of g's edges inside f (closed).
            let face = &pair.get(op).faces[f];
            let gp = pair.get(other(op));
            for l in &gp.faces[g].loops {
                for &(e, _) in l {
                    let along = pts.on_edge(other(op), e);
                    for w in along.windows(2) {
                        let m = (pts.rep[w[0]] + pts.rep[w[1]]) * 0.5;
                        if in_closed_face(pair.get(op), face, m, eps) {
                            let k = key(w[0], w[1]);
                            face_segs[slot(op)][f].push(FaceSeg { key: k, side: None });
                        }
                    }
                }
            }
        }
    }

    // Split every face into pieces and classify them.
    let mut partners: [HashMap<usize, Vec<usize>>; 2] = [HashMap::new(), HashMap::new()];
    for &(fa, fb) in &coplanar {
        partners[0].entry(fa).or_default().push(fb);
        partners[1].entry(fb).or_default().push(fa);
    }
    let windings = [
        classify::Winding::new(pair.first.brep, tol)?,
        classify::Winding::new(pair.second.brep, tol)?,
    ];
    let mut pieces: Vec<Piece> = Vec::new();
    for op in [Operand::First, Operand::Second] {
        let poly = pair.get(op);
        let them = pair.get(other(op));
        for (fi, f) in poly.faces.iter().enumerate() {
            let list = &face_segs[slot(op)][fi];
            let pos = |c: usize| f.to2(pts.rep[c]);
            let regions = planar::regions(list, &pos, eps).map_err(|p| match p {
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
            for loops in regions {
                let poly2: Vec<Vec<DVec2>> = loops
                    .iter()
                    .map(|l| l.iter().map(|&c| pos(c)).collect())
                    .collect();
                let (inner, clearance) = classify::interior_point(&poly2)
                    .map_err(|m| KernelError::InvalidResult(format!("boolean: face {fi}: {m}")))?;
                let at = f.to3(inner);
                if clearance <= eps {
                    return Err(KernelError::InvalidInput(format!(
                        "the bodies leave a sliver thinner than the tolerance near {}",
                        show(at)
                    )));
                }
                let mut class = None;
                let mut partner = None;
                for &g in partners[slot(op)].get(&fi).into_iter().flatten() {
                    let gf = &them.faces[g];
                    if in_closed_face(them, gf, at, eps) {
                        let same = gf.normal().dot(f.normal()) > 0.0;
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
                    None => windings[slot(other(op))].classify(at, eps)?,
                };
                pieces.push(Piece {
                    operand: op,
                    face: fi,
                    loops,
                    class,
                    partner,
                });
            }
        }
    }
    assemble::assemble(&pair, &pts, &segs, pieces, op, tol)
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
    sources: HashMap<Key, Vec<Source>>,
}

impl Segments {
    fn add(&mut self, a: usize, b: usize, s: Source) -> Key {
        let k = key(a, b);
        let list = self.sources.entry(k).or_default();
        if !list.contains(&s) {
            list.push(s);
        }
        k
    }
}

fn relation(pair: &Pair, pts: &Points, fa: usize, fb: usize) -> Result<Relation, KernelError> {
    let a = &pair.first.faces[fa];
    let b = &pair.second.faces[fb];
    let signs = |op: Operand, f: &PolyFace, plane: usize| -> (bool, bool, bool) {
        let (mut pos, mut neg, mut zero) = (false, false, false);
        for &v in f.verts.iter().flatten() {
            match pts.sign(op, v, plane) {
                1 => pos = true,
                -1 => neg = true,
                _ => zero = true,
            }
        }
        (pos, neg, zero)
    };
    let (ap, an, az) = signs(Operand::First, a, b.plane);
    let (bp, bn, bz) = signs(Operand::Second, b, a.plane);
    let a_on = az && !ap && !an;
    let b_on = bz && !bp && !bn;
    if a_on && b_on {
        return Ok(Relation::Coplanar);
    }
    if a_on || b_on {
        let p = pair.first.points[a.verts[0][0]];
        return Err(KernelError::InvalidInput(format!(
            "two faces of the bodies nearly coincide but are not on one plane, \
             within the tolerance, near {}",
            show(p)
        )));
    }
    let apart = |p: bool, n: bool, z: bool| !z && (p != n);
    if apart(ap, an, az) || apart(bp, bn, bz) {
        return Ok(Relation::Apart);
    }
    Ok(Relation::Crossing)
}

/// A closed interval of the common line of two face planes, between two
/// clusters, with their positions along the line.
#[derive(Clone, Copy, Debug)]
struct Span {
    lo: (f64, usize),
    hi: (f64, usize),
}

/// The stretches of the common line of two crossing faces that lie in both
/// (closed) faces, as pairs of clusters.
fn meeting(
    pair: &Pair,
    pts: &Points,
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
    let sa = face_line(pair, pts, Operand::First, fa, b.plane, dir)?;
    let sb = face_line(pair, pts, Operand::Second, fb, a.plane, dir)?;
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
    Ok(out)
}

fn compare(
    x: (f64, usize),
    y: (f64, usize),
    eps: f64,
    pts: &Points,
) -> Result<std::cmp::Ordering, KernelError> {
    if x.1 == y.1 {
        return Ok(std::cmp::Ordering::Equal);
    }
    if (x.0 - y.0).abs() <= eps {
        return Err(KernelError::InvalidInput(format!(
            "the bodies meet in points closer than the tolerance near {}",
            show(pts.rep[x.1])
        )));
    }
    Ok(x.0.total_cmp(&y.0))
}

fn later(
    x: (f64, usize),
    y: (f64, usize),
    eps: f64,
    pts: &Points,
) -> Result<(f64, usize), KernelError> {
    Ok(match compare(x, y, eps, pts)? {
        std::cmp::Ordering::Less => y,
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
        std::cmp::Ordering::Greater => y,
        _ => x,
    })
}

/// The closed face `f` of `op` cut by `plane` of the other operand: the
/// intervals of the planes' common line (direction `dir`) inside the face.
///
/// The crossings of the face's boundary with the plane come from the signs
/// of its vertices against the plane (decided once, `Points::sign`). Two
/// half-open rules, counting a vertex on the plane as on its positive side
/// and then as on its negative side, give the interior together with the
/// boundary edges lying in the plane from either side; their union is the
/// closed face's intersection with the line.
fn face_line(
    pair: &Pair,
    pts: &Points,
    op: Operand,
    f: usize,
    plane: usize,
    dir: DVec3,
) -> Result<Vec<Span>, KernelError> {
    let poly = pair.get(op);
    let face = &poly.faces[f];
    let eps = pts.eps;
    let mut spans: Vec<Span> = Vec::new();
    for zero_positive in [true, false] {
        let side = |s: i8| if zero_positive { s >= 0 } else { s > 0 };
        let mut cuts: Vec<(f64, usize)> = Vec::new();
        for (l, vs) in face.loops.iter().zip(&face.verts) {
            for (k, &(e, _)) in l.iter().enumerate() {
                let (va, vb) = (vs[k], vs[(k + 1) % vs.len()]);
                let (sa, sb) = (pts.sign(op, va, plane), pts.sign(op, vb, plane));
                if side(sa) == side(sb) {
                    continue;
                }
                let c = if sa == 0 {
                    pts.vertex_cluster(op, va)
                } else if sb == 0 {
                    pts.vertex_cluster(op, vb)
                } else {
                    pts.crossing_cluster(op, e, plane).ok_or_else(|| {
                        KernelError::InvalidResult(format!(
                            "boolean: no crossing recorded for edge {e} of the {} body",
                            name(op)
                        ))
                    })?
                };
                cuts.push((pts.rep[c].dot(dir), c));
            }
        }
        cuts.sort_by(|x, y| x.0.total_cmp(&y.0).then(x.1.cmp(&y.1)));
        for w in cuts.windows(2) {
            if w[0].1 != w[1].1 && (w[1].0 - w[0].0) <= eps {
                return Err(KernelError::InvalidInput(format!(
                    "the bodies meet in points closer than the tolerance near {}",
                    show(pts.rep[w[0].1])
                )));
            }
        }
        if cuts.len() % 2 != 0 {
            return Err(KernelError::InvalidResult(format!(
                "boolean: face {f} of the {} body crosses a plane an odd number of times",
                name(op)
            )));
        }
        for c in cuts.chunks(2) {
            spans.push(Span { lo: c[0], hi: c[1] });
        }
    }
    // Union of the two rules' intervals.
    spans.sort_by(|x, y| x.lo.0.total_cmp(&y.lo.0));
    let mut out: Vec<Span> = Vec::new();
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

/// True if `p` (on the face's plane) lies in the closed face: within the
/// tolerance of one of its edges, or inside it.
fn in_closed_face(poly: &Poly, face: &PolyFace, p: DVec3, eps: f64) -> bool {
    for l in &face.loops {
        for &(e, _) in l {
            let (a, b) = poly.edge_points(e);
            if points::segment_distance(p, a, b) <= eps {
                return true;
            }
        }
    }
    let q = face.to2(p);
    let mut inside = false;
    for vs in &face.verts {
        for k in 0..vs.len() {
            let a = face.to2(poly.points[vs[k]]);
            let b = face.to2(poly.points[vs[(k + 1) % vs.len()]]);
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
