//! Selecting, sewing and simplifying the pieces into the result (design
//! note 0012 §4.4 and §4.5), and its origins (ADR 0020).

use std::collections::{BTreeMap, BTreeSet};

use glam::{DVec2, DVec3};
use rustc_hash::{FxHashMap, FxHashSet};
use spade::Triangulation;

use super::chart::Chart;
use super::classify::{self, Class};
use super::meet::Track;
use super::points::{segment_distance, Points};
use super::quadric::Quadric;
use super::{
    corners, imprint, key, planar, planar_edges, shape, show, slot, Course, Key, Pair, Piece,
    Segments, Source,
};
use crate::geom::{Curve3, Surface};
use crate::math::Frame;
use crate::origin::{
    EdgeIx, EdgeOrigin, FaceIx, FaceOrigin, Merged, Operand, Origins, VertexIx, VertexOrigin,
};
use crate::tess;
use crate::topo::{Brep, Coedge, EdgeId, Loop, ShellId, VertexId};
use crate::{BooleanOp, KernelError, Tolerance};

/// The geometry the assembly reads: the operands, the curves where their
/// surfaces meet, and their faces' charts.
pub(super) struct Geometry<'a> {
    pub pair: &'a Pair<'a>,
    pub tracks: &'a [Track],
    pub charts: &'a [Vec<Chart>; 2],
}

impl Geometry<'_> {
    /// True if the operand face `f` of `op` is on a plane.
    fn planar(&self, op: Operand, f: usize) -> bool {
        let shape = self.pair.get(op);
        matches!(
            shape.surfaces[shape.faces[f].surface],
            Quadric::Plane { .. }
        )
    }

    /// The chart of a face of the result, in which its loops run
    /// counter-clockwise.
    fn chart(&self, f: &Out) -> Chart {
        let c = self.charts[slot(f.operand)][f.face];
        if f.reversed {
            c.reversed()
        } else {
            c
        }
    }

    /// The point of `course` at parameter `t`.
    fn point(&self, course: Course, t: f64) -> DVec3 {
        match course {
            Course::Straight => unreachable!("a straight segment has no parameter"),
            Course::Edge { op, e, .. } => self.pair.get(op).curve(e).0.point(t),
            Course::Track { id, .. } => self.tracks[id].point(t),
        }
    }
}

/// A face of the result being assembled.
struct Out {
    /// The operand face whose surface it lies on (the least of its
    /// constituents).
    operand: Operand,
    face: usize,
    /// True if its outward normal is the opposite of that face's.
    reversed: bool,
    /// Cycles of segments, outer first, counter-clockwise about the
    /// outward normal, each segment with whether it is run from its `a`.
    loops: Vec<Vec<(usize, bool)>>,
    /// Every operand face it is made of or coincides with.
    sources: BTreeSet<(Operand, usize)>,
    /// True if it joins several pieces (of one face or of faces on one
    /// plane).
    joined: bool,
}

/// Whether a piece is kept, and if so reversed (design note 0012 §4.4).
fn kept(op: BooleanOp, p: &Piece) -> Option<bool> {
    use BooleanOp::*;
    use Class::*;
    match (op, p.operand, p.class) {
        (Union, Operand::First, Outside | OnSame) => Some(false),
        (Union, Operand::Second, Outside) => Some(false),
        (Intersect, Operand::First, Inside | OnSame) => Some(false),
        (Intersect, Operand::Second, Inside) => Some(false),
        (Subtract, Operand::First, Outside | OnOpposite) => Some(false),
        (Subtract, Operand::Second, Inside) => Some(true),
        _ => None,
    }
}

pub(super) fn assemble(
    geometry: &Geometry,
    pts: &Points,
    mut segs: Segments,
    pieces: Vec<Piece>,
    coplanar: &FxHashSet<(usize, usize)>,
    op: BooleanOp,
    tol: &Tolerance,
) -> Result<(Brep, Origins), KernelError> {
    let eps = tol.linear;
    let mut faces: Vec<Out> = Vec::new();
    for p in pieces {
        let Some(reversed) = kept(op, &p) else {
            continue;
        };
        let mut sources = BTreeSet::from([(p.operand, p.face)]);
        if p.class == Class::OnSame {
            if let Some(g) = p.partner {
                sources.insert((super::other(p.operand), g));
            }
        }
        let mut loops = p.loops;
        if reversed {
            // Run the other way, from the same cluster as the cycle's last
            // segment started at.
            for l in &mut loops {
                l.reverse();
                for e in l.iter_mut() {
                    e.1 = !e.1;
                }
                l.rotate_left(1);
            }
        }
        faces.push(Out {
            operand: p.operand,
            face: p.face,
            reversed,
            loops,
            sources,
            joined: false,
        });
    }
    if faces.is_empty() {
        return Err(KernelError::InvalidInput(
            "the result would be empty".into(),
        ));
    }
    check_closed(&faces, pts, &segs)?;
    let faces = merge_faces(geometry, pts, &segs, faces, coplanar, eps)?;
    let faces = merge_edges(geometry, faces, pts, &mut segs, eps);
    check_closed(&faces, pts, &segs)?;
    check_vertices(&faces, pts, &segs)?;
    let shells = shells(&faces);
    let pair = geometry.pair;
    // Faces on planes bounded by straight segments: the shells' volumes
    // and the voids' places are exact on the clusters' places.
    let flat = faces.iter().all(|f| {
        geometry.planar(f.operand, f.face)
            && directed(f).all(|(s, _)| segs.list[s].course == Course::Straight)
    });
    let solids = if flat {
        Some(solids(pair, pts, &segs, &faces, &shells)?)
    } else {
        None
    };
    let (mut brep, org, shell_ids) = build(geometry, pts, &segs, &faces, &shells)?;
    let solids = match solids {
        Some(s) => s,
        None => curved_solids(&brep, &shell_ids, tol)?,
    };
    for (outer, voids) in solids {
        let mut list = vec![shell_ids[outer]];
        list.extend(voids.iter().map(|&v| shell_ids[v]));
        brep.add_solid(list);
    }
    Ok((brep, org))
}

/// The segments of a face's loops, each with the way it is run.
fn directed(f: &Out) -> impl Iterator<Item = (usize, bool)> + '_ {
    f.loops.iter().flatten().copied()
}

/// Every segment must be used once in each direction: more uses are a
/// non-manifold result, fewer an inconsistency. Problems are reported in
/// order of the segments' ends.
fn check_closed(faces: &[Out], pts: &Points, segs: &Segments) -> Result<(), KernelError> {
    let mut uses: BTreeMap<(Key, usize), (usize, usize)> = BTreeMap::new();
    for f in faces {
        for (s, fw) in directed(f) {
            let seg = &segs.list[s];
            let e = uses.entry((key(seg.a, seg.b), s)).or_default();
            if fw {
                e.0 += 1;
            } else {
                e.1 += 1;
            }
        }
    }
    for ((k, _), (fw, bw)) in uses {
        let at = (pts.rep[k.0] + pts.rep[k.1]) * 0.5;
        if fw + bw > 2 {
            return Err(KernelError::InvalidInput(format!(
                "the result would be non-manifold: the bodies touch along an edge near {}",
                show(at)
            )));
        }
        if fw != 1 || bw != 1 {
            return Err(KernelError::InvalidResult(format!(
                "boolean: the kept pieces do not close along an edge near {}",
                show(at)
            )));
        }
    }
    Ok(())
}

/// The frame a face of the result lies in, with `z` its outward normal, so
/// its loops run counter-clockwise in `x`, `y`.
fn plane_of(pair: &Pair, f: &Out) -> Frame {
    let fr = pair.get(f.operand).faces[f.face].frame;
    if f.reversed {
        Frame {
            origin: fr.origin,
            x: fr.x,
            y: -fr.y,
            z: -fr.z,
        }
    } else {
        fr
    }
}

fn to2(fr: &Frame, p: DVec3) -> DVec2 {
    let d = p - fr.origin;
    DVec2::new(d.dot(fr.x), d.dot(fr.y))
}

fn find(p: &mut [usize], mut i: usize) -> usize {
    while p[i] != i {
        p[i] = p[p[i]];
        i = p[i];
    }
    i
}

fn unite(p: &mut [usize], i: usize, j: usize) {
    let (ri, rj) = (find(p, i), find(p, j));
    if ri != rj {
        p[ri.max(rj)] = ri.min(rj);
    }
}

/// Joins faces that meet along a segment, lie on one plane and face the
/// same way (design note 0012 §4.5): pieces of one face cut and kept on
/// both sides, and coplanar faces of the two operands.
///
/// Whether two faces lie on one plane was decided before, and is not
/// tested again here: faces of one operand do if they are on one of its
/// planes (built on one plane, up to rounding), and faces of the two
/// operands if they were found to coincide (`coplanar`, pairs of the first
/// operand's face and the second's, every vertex of each within the
/// tolerance of the other's plane). A test on the faces' normals instead
/// would either join two faces of one body meeting at a slight crease,
/// moving corners of the joined face off its plane, or refuse to join
/// small coinciding faces tilted within the tolerance.
fn merge_faces(
    geometry: &Geometry,
    pts: &Points,
    segs: &Segments,
    faces: Vec<Out>,
    coplanar: &FxHashSet<(usize, usize)>,
    eps: f64,
) -> Result<Vec<Out>, KernelError> {
    let pair = geometry.pair;
    let n = faces.len();
    let mut parent: Vec<usize> = (0..n).collect();
    let planar: Vec<bool> = faces
        .iter()
        .map(|f| geometry.planar(f.operand, f.face))
        .collect();
    let frames: Vec<Frame> = faces
        .iter()
        .zip(&planar)
        .map(|(f, &p)| if p { plane_of(pair, f) } else { Frame::WORLD })
        .collect();
    let on_one_plane = |x: (Operand, usize), y: (Operand, usize)| match (x.0, y.0) {
        (Operand::First, Operand::Second) => coplanar.contains(&(x.1, y.1)),
        (Operand::Second, Operand::First) => coplanar.contains(&(y.1, x.1)),
        _ => {
            let shape = pair.get(x.0);
            shape.faces[x.1].surface == shape.faces[y.1].surface
        }
    };
    // Faces on a curved surface are joined only with pieces of the same
    // operand face facing the same way: those span no more than the face,
    // at most half a turn (ADR 0008).
    let same = |i: usize, j: usize| {
        if !planar[i] || !planar[j] {
            let (x, y) = (&faces[i], &faces[j]);
            return !planar[i]
                && !planar[j]
                && (x.operand, x.face, x.reversed) == (y.operand, y.face, y.reversed);
        }
        frames[i].z.dot(frames[j].z) > 0.0
            && faces[i]
                .sources
                .iter()
                .any(|&x| faces[j].sources.iter().any(|&y| on_one_plane(x, y)))
    };
    let mut by_seg: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
    for (i, f) in faces.iter().enumerate() {
        for (s, _) in directed(f) {
            by_seg.entry(s).or_default().push(i);
        }
    }
    for fs in by_seg.values() {
        if let [i, j] = fs[..] {
            if i != j && same(i, j) {
                unite(&mut parent, i, j);
            }
        }
    }
    let mut groups: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
    for i in 0..n {
        let r = find(&mut parent, i);
        groups.entry(r).or_default().push(i);
    }
    let mut faces: Vec<Option<Out>> = faces.into_iter().map(Some).collect();
    let mut out = Vec::new();
    for members in groups.into_values() {
        let mut parts: Vec<Out> = members
            .iter()
            .map(|&i| faces[i].take().expect("each face is in one group"))
            .collect();
        if parts.len() == 1 {
            out.append(&mut parts);
            continue;
        }
        // The least operand face carries the merged face.
        let primary = parts
            .iter()
            .min_by_key(|f| (f.operand, f.face))
            .expect("a group has members");
        let chart = geometry.chart(primary);
        let (operand, face, reversed) = (primary.operand, primary.face, primary.reversed);
        let sources: BTreeSet<(Operand, usize)> = parts
            .iter()
            .flat_map(|f| f.sources.iter().copied())
            .collect();
        let mut all: BTreeSet<(usize, bool)> = BTreeSet::new();
        for f in &parts {
            all.extend(directed(f));
        }
        // The uses not cancelled by a use the other way, in order of their
        // ends as run.
        let mut boundary: Vec<(Key, usize, bool)> = all
            .iter()
            .filter(|(s, fw)| !all.contains(&(*s, !*fw)))
            .map(|&(s, fw)| ((segs.start(s, fw), segs.end(s, fw)), s, fw))
            .collect();
        boundary.sort_unstable_by_key(|x| (x.0, x.1));
        let edges = planar_edges(
            pair,
            geometry.tracks,
            segs,
            &chart,
            boundary.iter().map(|&(_, s, fw)| (s, Some(fw))),
        )?;
        let reach = |q: DVec2| chart.reach(q, eps);
        let regions =
            planar::regions(&edges, &|c| chart.map(pts.rep[c]), &reach, eps).map_err(|p| {
                let (m, c) = p.describe();
                KernelError::InvalidResult(format!(
                    "boolean: merging faces: {m} near {}",
                    super::show(pts.rep[c])
                ))
            })?;
        for region in regions {
            out.push(Out {
                operand,
                face,
                reversed,
                loops: region.loops,
                sources: sources.clone(),
                joined: true,
            });
        }
    }
    Ok(out)
}

/// Joins two segments end to end where they meet at a point no other
/// segment reaches, lie on one line and bound the same two faces.
///
/// Joins are made in rounds: each round takes, in order of their points,
/// the joins whose three points no earlier join of the round touches. The
/// segments at each point and the faces along each segment are kept up to
/// date as joins are made, and a round looks again only at the points whose
/// segments changed or that an earlier join of the last round touched, so
/// a long run of collinear segments costs time in proportion to its length
/// rather than to the whole result for every round.
fn merge_edges(
    geometry: &Geometry,
    mut faces: Vec<Out>,
    pts: &Points,
    segs: &mut Segments,
    eps: f64,
) -> Vec<Out> {
    let mut at: FxHashMap<usize, BTreeSet<usize>> = FxHashMap::default();
    let mut seg_faces: FxHashMap<usize, BTreeSet<usize>> = FxHashMap::default();
    for (i, f) in faces.iter().enumerate() {
        for (s, _) in directed(f) {
            let seg = &segs.list[s];
            at.entry(seg.a).or_default().insert(s);
            at.entry(seg.b).or_default().insert(s);
            seg_faces.entry(s).or_default().insert(i);
        }
    }
    let mut look: BTreeSet<usize> = at.keys().copied().collect();
    let mut gone: FxHashSet<usize> = FxHashSet::default();
    // Segments joined into another, and the one they joined.
    let mut joined: FxHashMap<usize, usize> = FxHashMap::default();
    loop {
        let mut touched: FxHashSet<usize> = FxHashSet::default();
        let mut joins: Vec<(usize, usize, usize, usize, usize, Course)> = Vec::new();
        let mut again: BTreeSet<usize> = BTreeSet::new();
        for &v in &look {
            let Some(ss) = at.get(&v) else {
                continue;
            };
            let [s0, s1] = ss.iter().copied().collect::<Vec<usize>>()[..] else {
                continue;
            };
            let (g0, g1) = (&segs.list[s0], &segs.list[s1]);
            let u = if g0.a == v { g0.b } else { g0.a };
            let w = if g1.a == v { g1.b } else { g1.a };
            if u == w || seg_faces[&s0] != seg_faces[&s1] {
                continue;
            }
            // A point on the axis of a face's surface of revolution stays a
            // vertex: the surface is singular there, and an edge running
            // through it would have no parameters there on the face (its
            // angle jumps by half a turn).
            if seg_faces[&s0]
                .iter()
                .any(|&i| on_axis(geometry.pair, &faces[i], pts.rep[v], eps))
            {
                continue;
            }
            // The courses from u to v and from v to w.
            let c0 = if g0.a == u {
                g0.course
            } else {
                g0.course.reversed()
            };
            let c1 = if g1.a == v {
                g1.course
            } else {
                g1.course.reversed()
            };
            let joined = match (c0, c1) {
                (Course::Straight, Course::Straight) => {
                    if segment_distance(pts.rep[v], pts.rep[u], pts.rep[w]) > eps {
                        continue;
                    }
                    Course::Straight
                }
                // Pieces of one curve: its parameters at v agree (an
                // operand's edge), or differ by whole turns (a track).
                (
                    Course::Edge { op, e, ta, tb },
                    Course::Edge {
                        op: op1,
                        e: e1,
                        ta: tv,
                        tb: tw,
                    },
                ) if (op, e) == (op1, e1) => Course::Edge {
                    op,
                    e,
                    ta,
                    tb: tw + (tb - tv),
                },
                (
                    Course::Track { id, ta, tb },
                    Course::Track {
                        id: id1,
                        ta: tv,
                        tb: tw,
                    },
                ) if id == id1 => Course::Track {
                    id,
                    ta,
                    tb: tw + (tb - tv),
                },
                _ => continue,
            };
            if [u, v, w].iter().any(|c| touched.contains(c)) {
                again.insert(v);
                continue;
            }
            touched.extend([u, v, w]);
            joins.push((u, v, w, s0, s1, joined));
        }
        if joins.is_empty() {
            break;
        }
        for (u, v, w, s0, s1, course) in joins {
            let mut sources: BTreeSet<Source> = segs.list[s0].sources.iter().copied().collect();
            sources.extend(segs.list[s1].sources.iter().copied());
            let mid = match course {
                Course::Straight => None,
                Course::Edge { ta, tb, .. } | Course::Track { ta, tb, .. } => {
                    Some(geometry.point(course, 0.5 * (ta + tb)))
                }
            };
            let mut add = |x: Source| match mid {
                None => segs.add(u, w, Course::Straight, x),
                Some(m) => segs.add_curved(u, w, course, m, x, eps),
            };
            let mut sources = sources.into_iter();
            let first = sources.next().expect("a segment has a source");
            let s = add(first);
            for x in sources {
                add(x);
            }
            let fs = seg_faces.remove(&s0).unwrap_or_default();
            seg_faces.remove(&s1);
            seg_faces.entry(s).or_default().extend(fs);
            at.remove(&v);
            for (x, old) in [(u, s0), (w, s1)] {
                if let Some(ss) = at.get_mut(&x) {
                    ss.remove(&old);
                    ss.insert(s);
                }
            }
            joined.insert(s0, s);
            joined.insert(s1, s);
            gone.insert(v);
            again.extend([u, w]);
        }
        look = again;
    }
    if !gone.is_empty() {
        let resolve = |mut s: usize| {
            while let Some(&t) = joined.get(&s) {
                s = t;
            }
            s
        };
        for f in &mut faces {
            for l in &mut f.loops {
                let old = std::mem::take(l);
                for (s, fw) in old {
                    let start = segs.start(s, fw);
                    if gone.contains(&start) {
                        continue;
                    }
                    let t = resolve(s);
                    l.push((t, segs.list[t].a == start));
                }
            }
        }
    }
    faces
}

/// True if `p` lies within `eps` of the axis of the surface of `f`'s
/// operand face, if that is a surface of revolution: its singular points
/// (a pole, an apex, the centre of a flat cap) are there.
fn on_axis(pair: &Pair, f: &Out, p: DVec3, eps: f64) -> bool {
    match &pair.get(f.operand).brep.faces[f.face].surface {
        Surface::Revolution { frame, .. } => (p - frame.origin).cross(frame.z).length() <= eps,
        _ => false,
    }
}

/// Around every vertex the faces must form one fan: two fans are bodies
/// touching at a point, which is non-manifold.
fn check_vertices(faces: &[Out], pts: &Points, segs: &Segments) -> Result<(), KernelError> {
    // Corners: (vertex, incoming segment, outgoing segment), each segment
    // with the way it is run.
    type Corner = (usize, (usize, bool), (usize, bool));
    let mut corners: Vec<Corner> = Vec::new();
    for f in faces {
        for l in &f.loops {
            let n = l.len();
            for i in 0..n {
                let (p, q) = (l[(i + n - 1) % n], l[i]);
                corners.push((segs.start(q.0, q.1), p, q));
            }
        }
    }
    let by_in: FxHashMap<(usize, bool), usize> =
        corners.iter().enumerate().map(|(i, c)| (c.1, i)).collect();
    let mut seen = vec![false; corners.len()];
    let mut fans: BTreeMap<usize, usize> = BTreeMap::new();
    for i in 0..corners.len() {
        if seen[i] {
            continue;
        }
        *fans.entry(corners[i].0).or_default() += 1;
        let mut c = i;
        while !seen[c] {
            seen[c] = true;
            let (v, _, (s, fw)) = corners[c];
            match by_in.get(&(s, !fw)) {
                Some(&d) => c = d,
                None => {
                    return Err(KernelError::InvalidResult(format!(
                        "boolean: a corner with no neighbour near {}",
                        show(pts.rep[v])
                    )))
                }
            }
        }
    }
    for (v, n) in fans {
        if n > 1 {
            return Err(KernelError::InvalidInput(format!(
                "the result would be non-manifold: the bodies touch at a point near {}",
                show(pts.rep[v])
            )));
        }
    }
    Ok(())
}

/// Connected sets of faces, in order of their least face.
fn shells(faces: &[Out]) -> Vec<Vec<usize>> {
    let mut parent: Vec<usize> = (0..faces.len()).collect();
    let mut first: FxHashMap<usize, usize> = FxHashMap::default();
    for (i, f) in faces.iter().enumerate() {
        for (s, _) in directed(f) {
            match first.get(&s) {
                Some(&j) => unite(&mut parent, i, j),
                None => {
                    first.insert(s, i);
                }
            }
        }
    }
    let mut groups: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
    for i in 0..faces.len() {
        let r = find(&mut parent, i);
        groups.entry(r).or_default().push(i);
    }
    groups.into_values().collect()
}

/// The enclosed volume of a shell of planar faces: a third of the sum of
/// each face's area times its plane's distance from the origin along its
/// outward normal.
fn shell_volume(pair: &Pair, pts: &Points, segs: &Segments, faces: &[Out], shell: &[usize]) -> f64 {
    let mut v = 0.0;
    for &i in shell {
        let f = &faces[i];
        let fr = plane_of(pair, f);
        let area: f64 = f
            .loops
            .iter()
            .map(|l| {
                let q: Vec<DVec2> = corners(segs, l)
                    .into_iter()
                    .map(|c| to2(&fr, pts.rep[c]))
                    .collect();
                (0..q.len())
                    .map(|k| q[k].perp_dot(q[(k + 1) % q.len()]))
                    .sum::<f64>()
                    * 0.5
            })
            .sum();
        v += fr.origin.dot(fr.z) * area / 3.0;
    }
    v
}

/// Outer shells (enclosing positive volume), each with the void shells
/// (negative volume) inside it: the smallest outer shell around each void.
fn solids(
    pair: &Pair,
    pts: &Points,
    segs: &Segments,
    faces: &[Out],
    shells: &[Vec<usize>],
) -> Result<Vec<(usize, Vec<usize>)>, KernelError> {
    let vols: Vec<f64> = shells
        .iter()
        .map(|s| shell_volume(pair, pts, segs, faces, s))
        .collect();
    let mut solids: Vec<(usize, Vec<usize>)> = Vec::new();
    for (i, v) in vols.iter().enumerate() {
        if *v > 0.0 {
            solids.push((i, Vec::new()));
        } else if *v == 0.0 {
            return Err(KernelError::InvalidResult(
                "boolean: a shell encloses no volume".into(),
            ));
        }
    }
    let voids: Vec<usize> = (0..shells.len()).filter(|&i| vols[i] < 0.0).collect();
    if voids.is_empty() {
        return Ok(solids);
    }
    let meshes: Vec<Vec<[DVec3; 3]>> = solids
        .iter()
        .map(|(s, _)| shell_triangles(pair, pts, segs, faces, &shells[*s]))
        .collect::<Result<_, _>>()?;
    for v in voids {
        let (s, fw) = faces[shells[v][0]].loops[0][0];
        let q = pts.rep[segs.start(s, fw)];
        let mut best: Option<usize> = None;
        for (k, (s, _)) in solids.iter().enumerate() {
            if winding(&meshes[k], q) > 0.5 && best.is_none_or(|b| vols[*s] < vols[solids[b].0]) {
                best = Some(k);
            }
        }
        match best {
            Some(k) => solids[k].1.push(v),
            None => {
                return Err(KernelError::InvalidResult(format!(
                    "boolean: a void near {} is inside no solid",
                    show(q)
                )))
            }
        }
    }
    Ok(solids)
}

fn shell_triangles(
    pair: &Pair,
    pts: &Points,
    segs: &Segments,
    faces: &[Out],
    shell: &[usize],
) -> Result<Vec<[DVec3; 3]>, KernelError> {
    let mut out = Vec::new();
    for &i in shell {
        let f = &faces[i];
        let fr = plane_of(pair, f);
        let cycles: Vec<Vec<usize>> = f.loops.iter().map(|l| corners(segs, l)).collect();
        let loops: Vec<Vec<DVec2>> = cycles
            .iter()
            .map(|l| l.iter().map(|&c| to2(&fr, pts.rep[c])).collect())
            .collect();
        let (cdt, index) = tess::polygon_cdt(&loops)
            .map_err(|e| KernelError::InvalidResult(format!("boolean: triangulation: {e}")))?;
        // Each vertex's first cluster.
        let mut at = vec![usize::MAX; cdt.num_vertices()];
        for (&i, &c) in index.iter().zip(cycles.iter().flatten()) {
            if at[i] == usize::MAX {
                at[i] = c;
            }
        }
        for t in tess::inside_triangles(&cdt) {
            let mut tri = t.map(|k| pts.rep[at[k]]);
            // Triangles of a counter-clockwise boundary are
            // counter-clockwise; keep them facing outward.
            let n = (tri[1] - tri[0]).cross(tri[2] - tri[0]);
            if n.dot(fr.z) < 0.0 {
                tri.swap(1, 2);
            }
            out.push(tri);
        }
    }
    Ok(out)
}

fn winding(tris: &[[DVec3; 3]], p: DVec3) -> f64 {
    let mut total = 0.0;
    for t in tris {
        let (a, b, c) = (t[0] - p, t[1] - p, t[2] - p);
        let (la, lb, lc) = (a.length(), b.length(), c.length());
        let num = a.dot(b.cross(c));
        let den = la * lb * lc + a.dot(b) * lc + a.dot(c) * lb + b.dot(c) * la;
        total += 2.0 * num.atan2(den);
    }
    total / (4.0 * std::f64::consts::PI)
}

fn edge_origin(s: Source) -> EdgeOrigin {
    match s {
        Source::Edge(operand, e) => EdgeOrigin::Operand {
            operand,
            edge: EdgeIx(e as u32),
        },
        Source::Meet(a, b) => EdgeOrigin::Meet {
            first: FaceIx(a as u32),
            second: FaceIx(b as u32),
        },
    }
}

/// The B-rep and its origin table.
/// The edge geometry of segment `seg`: its curve, its parameter range,
/// and whether the curve runs from the segment's `b` to its `a`.
fn edge_curve(
    geometry: &Geometry,
    pts: &Points,
    seg: &super::Seg,
) -> Result<(Curve3, f64, f64, bool), KernelError> {
    let (a, b) = (pts.rep[seg.a], pts.rep[seg.b]);
    Ok(match seg.course {
        Course::Straight => (
            Curve3::Line {
                origin: a,
                dir: b - a,
            },
            0.0,
            1.0,
            false,
        ),
        Course::Edge { op, e, ta, tb } => {
            let curve = geometry.pair.get(op).curve(e).0.clone();
            (curve, ta.min(tb), ta.max(tb), ta > tb)
        }
        Course::Track { id, ta, tb } => {
            let (lo, hi) = (ta.min(tb), ta.max(tb));
            let (curve, t0, t1) = match geometry.tracks[id] {
                // A circle's frame turned so that the arc's middle is at a
                // quarter turn: its parameters are then within the kernel's
                // range for circles, [−π/2, 3π/2).
                Track::Circle {
                    center,
                    x,
                    y,
                    radius,
                } => {
                    let shift = 0.5 * (lo + hi) - std::f64::consts::FRAC_PI_2;
                    let (sn, cs) = shift.sin_cos();
                    (
                        Curve3::Circle {
                            center,
                            x: x * cs + y * sn,
                            y: y * cs - x * sn,
                            radius,
                        },
                        lo - shift,
                        hi - shift,
                    )
                }
                track => track.edge(lo, hi)?,
            };
            (curve, t0, t1, ta > tb)
        }
    })
}

/// The B-rep's vertices, edges, faces and shells, its origin table, and
/// its shells' ids (in the order of `shells`).
fn build(
    geometry: &Geometry,
    pts: &Points,
    segs: &Segments,
    faces: &[Out],
    shells: &[Vec<usize>],
) -> Result<(Brep, Origins, Vec<ShellId>), KernelError> {
    let pair = geometry.pair;
    let mut brep = Brep::default();
    let mut org = Origins::default();
    let mut merged: Vec<Merged> = Vec::new();
    let mut vertex_of: FxHashMap<usize, VertexId> = FxHashMap::default();
    // Per segment, its edge and whether the edge runs from its `b`.
    let mut edge_of: FxHashMap<usize, (EdgeId, bool)> = FxHashMap::default();
    let mut face_ids = Vec::with_capacity(faces.len());
    for f in faces {
        let mut loops = Vec::with_capacity(f.loops.len());
        for l in &f.loops {
            let mut coedges = Vec::with_capacity(l.len());
            for &(s, fw) in l {
                let seg = &segs.list[s];
                for c in [seg.a, seg.b] {
                    if let std::collections::hash_map::Entry::Vacant(e) = vertex_of.entry(c) {
                        let id = brep.add_vertex(pts.rep[c]);
                        e.insert(id);
                        let (fv, sv) = pts.vertices_of(c);
                        let at = VertexIx(id.0);
                        let first = |op, v: usize| VertexOrigin::Operand {
                            operand: op,
                            vertex: VertexIx(v as u32),
                        };
                        match (fv, sv) {
                            (Some(a), b) => {
                                org.vertices.push(first(Operand::First, a));
                                if let Some(b) = b {
                                    merged.push(Merged::Vertex(at, first(Operand::Second, b)));
                                }
                            }
                            (None, Some(b)) => org.vertices.push(first(Operand::Second, b)),
                            (None, None) => org.vertices.push(VertexOrigin::Meet),
                        }
                    }
                }
                let (edge, flipped) = match edge_of.get(&s) {
                    Some(&e) => e,
                    None => {
                        let (a, b) = (pts.rep[seg.a], pts.rep[seg.b]);
                        let (curve, t0, t1, flipped) = edge_curve(geometry, pts, seg)?;
                        let (v0, v1) = if flipped {
                            (vertex_of[&seg.b], vertex_of[&seg.a])
                        } else {
                            (vertex_of[&seg.a], vertex_of[&seg.b])
                        };
                        let id = brep.add_edge(curve, t0, t1, v0, v1);
                        edge_of.insert(s, (id, flipped));
                        let mut src = seg
                            .sources
                            .iter()
                            .copied()
                            .collect::<BTreeSet<_>>()
                            .into_iter();
                        let primary = src.next().ok_or_else(|| {
                            KernelError::InvalidResult(format!(
                                "boolean: an edge with no origin near {}",
                                show((a + b) * 0.5)
                            ))
                        })?;
                        org.edges.push(edge_origin(primary));
                        for s in src {
                            merged.push(Merged::Edge(EdgeIx(id.0), edge_origin(s)));
                        }
                        (id, flipped)
                    }
                };
                coedges.push(Coedge {
                    edge,
                    forward: fw != flipped,
                });
            }
            loops.push(Loop { coedges });
        }
        // A face joined from pieces on a plane recorded as a surface of
        // revolution (a revolve's cap, in faces of at most half a turn, ADR
        // 0008) can go around the axis: it is recorded as a plane.
        let src = &pair.get(f.operand).brep.faces[f.face];
        let id = if f.joined
            && geometry.planar(f.operand, f.face)
            && !matches!(src.surface, Surface::Plane(_))
        {
            brep.add_face(Surface::Plane(plane_of(pair, f)), true, loops)
        } else {
            brep.add_face(src.surface.clone(), src.forward != f.reversed, loops)
        };
        org.faces.push(FaceOrigin::Operand {
            operand: f.operand,
            face: FaceIx(f.face as u32),
        });
        for &(o, g) in &f.sources {
            if (o, g) != (f.operand, f.face) {
                merged.push(Merged::Face(
                    FaceIx(id.0),
                    FaceOrigin::Operand {
                        operand: o,
                        face: FaceIx(g as u32),
                    },
                ));
            }
        }
        face_ids.push(id);
    }
    let shell_ids: Vec<ShellId> = shells
        .iter()
        .map(|s| brep.add_shell(s.iter().map(|&i| face_ids[i]).collect()))
        .collect();
    merged.sort();
    org.merged = merged;
    Ok((brep, org, shell_ids))
}

/// Outer shells (enclosing positive volume) of a result with curved faces
/// or edges, each with the void shells inside it, as [`solids`] finds them
/// for flat ones. A lone shell is the outer shell of the one solid (were it
/// inside out, validation would refuse the body). With several, a shell's
/// volume is the one validation measures orientation by
/// ([`crate::check::shell_signed_volumes`]: from a coarse mesh, or from the
/// mass properties where the mesh's error could reach its sign); each void
/// goes in the smallest outer shell holding one of its vertices (by ray
/// parity on the shell's own faces).
fn curved_solids(
    b: &Brep,
    shells: &[ShellId],
    tol: &Tolerance,
) -> Result<Vec<(usize, Vec<usize>)>, KernelError> {
    if shells.len() == 1 {
        return Ok(vec![(0, Vec::new())]);
    }
    let eps = tol.linear;
    let signed = crate::check::shell_signed_volumes(b, tol)
        .map_err(|e| KernelError::InvalidResult(format!("boolean: the shells' volumes: {e}")))?;
    let vols: Vec<f64> = shells.iter().map(|s| signed[s.idx()]).collect();
    let bodies: Vec<Brep> = shells
        .iter()
        .map(|&s| crate::check::shell_brep(b, s))
        .collect();
    let body = |i: usize| &bodies[i];
    let mut solids: Vec<(usize, Vec<usize>)> = Vec::new();
    for (i, v) in vols.iter().enumerate() {
        if *v > 0.0 {
            solids.push((i, Vec::new()));
        } else if *v == 0.0 {
            return Err(KernelError::InvalidResult(
                "boolean: a shell encloses no volume".into(),
            ));
        }
    }
    for v in (0..shells.len()).filter(|&i| vols[i] < 0.0) {
        // A corner of the void shell itself (its body shares the whole
        // result's vertices and edges).
        let shell = &b.shells[shells[v].idx()];
        let first = b.faces[shell.faces[0].idx()].loops[0].coedges[0].edge;
        let q = b.vertices[b.edges[first.idx()].v0.idx()].point;
        let mut best: Option<usize> = None;
        for (k, (s, _)) in solids.iter().enumerate() {
            let outer = body(*s);
            let shape = shape::shape(outer, "result")?;
            let outlines = imprint::outlines(&shape)?;
            let boundary = classify::Boundary::new(&shape, &outlines)?;
            if boundary.classify(q, eps)? == Class::Inside
                && best.is_none_or(|b| vols[*s] < vols[solids[b].0])
            {
                best = Some(k);
            }
        }
        match best {
            Some(k) => solids[k].1.push(v),
            None => {
                return Err(KernelError::InvalidResult(format!(
                    "boolean: a void near {} is inside no solid",
                    show(q)
                )))
            }
        }
    }
    Ok(solids)
}
