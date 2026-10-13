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

use glam::DVec3;
use rustc_hash::{FxHashMap, FxHashSet};

use crate::bvh::Bvh;
use crate::origin::{Operand, Origins};
use crate::topo::Brep;
use crate::{BooleanOp, KernelError, Tolerance};

mod assemble;
// The geometry layer of step 2 (design note 0012 §10.7, item 1), tested on
// its own; the pipeline uses it as it gains curved faces.
#[allow(dead_code)]
mod bernstein;
mod classify;
#[allow(dead_code)]
mod meet;
mod planar;
mod points;
#[allow(dead_code)]
mod quadric;
#[allow(dead_code)]
mod shape;
#[cfg(test)]
mod testing;

use classify::Class;
use points::Points;
use shape::{Face, Shape};

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
    for (s, which) in [(&pair.first, "first"), (&pair.second, "second")] {
        if !s.flat {
            return Err(KernelError::Unsupported(format!(
                "booleans of bodies with curved faces or edges (the {which} body has some; \
                 design note 0012 §10)"
            )));
        }
    }
    let eps = tol.linear;
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
            match relation(&pair, &pts, fa, fb)? {
                Relation::Apart => {}
                Relation::Coplanar => coplanar.push((fa, fb)),
                Relation::Crossing => {
                    pts.add_crossings(&pair, Operand::First, fa, b.surface);
                    pts.add_crossings(&pair, Operand::Second, fb, a.surface);
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
                        edge_segs[slot(op)][e] = along
                            .windows(2)
                            .map(|w| segs.add(w[0], w[1], Course::Straight, Source::Edge(op, e)))
                            .collect();
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
    for &(fa, fb) in &crossing {
        for (lo, hi) in meeting(&pair, &pts, fa, fb)? {
            for w in pts.split_segment(lo, hi)?.windows(2) {
                let s = segs.add(w[0], w[1], Course::Straight, Source::Meet(fa, fb));
                face_segs[0][fa].push(FaceSeg { seg: s, side: None });
                face_segs[1][fb].push(FaceSeg { seg: s, side: None });
            }
        }
    }
    for &(fa, fb) in &coplanar {
        for (op, f, g) in [(Operand::First, fa, fb), (Operand::Second, fb, fa)] {
            // The pieces of g's edges inside f (closed).
            let face = &pair.get(op).faces[f];
            let gs = pair.get(other(op));
            for l in &gs.faces[g].loops {
                for &(e, _) in l {
                    for &s in &edge_segs[slot(other(op))][e] {
                        let seg = &segs.list[s];
                        let m = (pts.rep[seg.a] + pts.rep[seg.b]) * 0.5;
                        if in_closed_face(pair.get(op), face, m, eps) {
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
        classify::Boundary::new(&pair.first)?,
        classify::Boundary::new(&pair.second)?,
    ];
    let boxes = [pair.first.extent(), pair.second.extent()];
    let mut pieces: Vec<Piece> = Vec::new();
    for op in [Operand::First, Operand::Second] {
        let shape = pair.get(op);
        let them = pair.get(other(op));
        for (fi, f) in shape.faces.iter().enumerate() {
            // A face farther than the tolerance from the other body's box
            // is untouched by it (no point where the bodies meet is on its
            // edges, and no face of the other body is on its plane) and
            // wholly outside it: its one piece, its own loops, needs no
            // point inside it and no test.
            let clear = !boxes_overlap((f.lo, f.hi), boxes[slot(other(op))], eps);
            let edges: Vec<planar::Edge> = face_segs[slot(op)][fi]
                .iter()
                .map(|fs| {
                    let s = &segs.list[fs.seg];
                    planar::Edge {
                        seg: fs.seg,
                        a: s.a,
                        b: s.b,
                        side: fs.side,
                        bend: None,
                    }
                })
                .collect();
            let pos = |c: usize| f.to2(pts.rep[c]);
            let regions = planar::regions(&edges, &pos, eps).map_err(|p| match p {
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
                let (inner, clearance) = classify::interior_point(&region.polygons())
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
                    if in_face(them, gf, at) {
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
    assemble::assemble(&pair, &pts, segs, pieces, &coplanar, op, tol)
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
}

impl Segments {
    /// The segment from `from` to `to` along `course` from `source`: a new
    /// one, or the one already there (which gains the source).
    fn add(&mut self, from: usize, to: usize, course: Course, s: Source) -> usize {
        let k = key(from, to);
        let i = match course {
            Course::Straight => *self.straight.entry(k).or_insert_with(|| {
                self.list.push(Seg {
                    a: k.0,
                    b: k.1,
                    course,
                    sources: Vec::new(),
                });
                self.list.len() - 1
            }),
        };
        let list = &mut self.list[i].sources;
        if !list.contains(&s) {
            list.push(s);
        }
        i
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

fn relation(pair: &Pair, pts: &Points, fa: usize, fb: usize) -> Result<Relation, KernelError> {
    let a = &pair.first.faces[fa];
    let b = &pair.second.faces[fb];
    // Which sides of the plane the face's vertices are on: those of the
    // coedges near the plane one by one, the others by their boxes.
    let signs = |op: Operand, f: &Face, surface: usize| -> (bool, bool, bool) {
        let mut near = Vec::new();
        let q = &pair.get(other(op)).surfaces[surface];
        let (mut pos, mut neg) = f.near_surface(q, pts.eps, &mut near);
        let mut zero = false;
        for k in near {
            let (li, i) = f.coedge_at(k);
            match pts.sign(op, f.verts[li][i], surface) {
                1 => pos = true,
                -1 => neg = true,
                _ => zero = true,
            }
        }
        (pos, neg, zero)
    };
    let (ap, an, az) = signs(Operand::First, a, b.surface);
    let (bp, bn, bz) = signs(Operand::Second, b, a.surface);
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

/// The closed face `f` of `op` cut by the plane `surface` of the other
/// operand: the intervals of the planes' common line (direction `dir`)
/// inside the face.
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
    surface: usize,
    dir: DVec3,
) -> Result<Vec<Span>, KernelError> {
    let shape = pair.get(op);
    let face = &shape.faces[f];
    let eps = pts.eps;
    let mut spans: Vec<Span> = Vec::new();
    // Only coedges near the plane can cross it under either rule.
    let mut near = Vec::new();
    face.near_surface(&pair.get(other(op)).surfaces[surface], eps, &mut near);
    for zero_positive in [true, false] {
        let side = |s: i8| if zero_positive { s >= 0 } else { s > 0 };
        let mut cuts: Vec<(f64, usize)> = Vec::new();
        for &c in &near {
            let (li, k) = face.coedge_at(c);
            let (vs, e) = (&face.verts[li], face.loops[li][k].0);
            let (va, vb) = (vs[k], vs[(k + 1) % vs.len()]);
            let (sa, sb) = (pts.sign(op, va, surface), pts.sign(op, vb, surface));
            if side(sa) == side(sb) {
                continue;
            }
            let c = if sa == 0 {
                pts.vertex_cluster(op, va)
            } else if sb == 0 {
                pts.vertex_cluster(op, vb)
            } else {
                pts.crossing_cluster(op, e, surface).ok_or_else(|| {
                    KernelError::InvalidResult(format!(
                        "boolean: no crossing recorded for edge {e} of the {} body",
                        name(op)
                    ))
                })?
            };
            cuts.push((pts.rep[c].dot(dir), c));
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
fn in_closed_face(shape: &Shape, face: &Face, p: DVec3, eps: f64) -> bool {
    for l in &face.loops {
        for &(e, _) in l {
            let (a, b) = shape.edge_points(e);
            if points::segment_distance(p, a, b) <= eps {
                return true;
            }
        }
    }
    in_face(shape, face, p)
}

/// True if `p` (on the face's plane) lies inside the face, by the parity
/// of the face's edges crossed by a ray from it; for a point off its edges.
fn in_face(shape: &Shape, face: &Face, p: DVec3) -> bool {
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
