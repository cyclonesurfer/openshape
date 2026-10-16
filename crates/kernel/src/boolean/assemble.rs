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
    corners, imprint, key, planar, planar_edges, shape, show, slot, Course, Courses, Key, Pair,
    Piece, Segments, Source,
};
use crate::geom::{Curve3, Surface};
use crate::math::Frame;
use crate::origin::{
    EdgeIx, EdgeOrigin, FaceIx, FaceOrigin, Merged, Operand, Origins, VertexIx, VertexOrigin,
};
use crate::topo::{Brep, Coedge, EdgeId, Loop, ShellId, VertexId};
use crate::{check, tess};
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

    /// The courses of curved segments.
    fn courses(&self) -> Courses<'_> {
        Courses {
            pair: self.pair,
            tracks: self.tracks,
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

/// Where the bodies touch without crossing inside a face (design note
/// 0012 §18.2): along a segment, or at a cluster.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub(super) enum Place {
    Seg(usize),
    Point(usize),
}

/// A contact inside a piece: its place, the piece holding it (by index in
/// the pieces), the side of the piece's face the other body's touching
/// boundary lies on (+1 outside the piece's body), and a point of it (a
/// segment's middle).
pub(super) struct Contact {
    pub place: Place,
    pub holder: usize,
    pub sigma: i8,
    pub at: DVec3,
}

/// Decides every contact inside a face (design note 0012 §18.2): a place
/// pinches, making the result non-manifold, if and only if kept pieces of
/// both operands are at it: a piece holding it, and a piece of the other
/// operand *incident* to it (with its segment in its loops, or a corner at
/// its cluster). The least pinching place is refused, a segment before a
/// point. Cross-checks, each a kernel error: every place has an incident
/// piece, and its incident pieces are all outside the holder's body where
/// the touching boundary lies outside it (σ = +1) and all inside where it
/// lies inside.
fn check_contacts(
    op: BooleanOp,
    pieces: &[Piece],
    contacts: &[Contact],
    segs: &Segments,
) -> Result<(), KernelError> {
    if contacts.is_empty() {
        return Ok(());
    }
    let keep = |i: usize| kept(op, &pieces[i]).is_some();
    // Per place, in order: its holders, side and point.
    let mut places: BTreeMap<Place, (Vec<usize>, i8, DVec3)> = BTreeMap::new();
    for c in contacts {
        let e = places
            .entry(c.place)
            .or_insert_with(|| (Vec::new(), c.sigma, c.at));
        if e.1 != c.sigma {
            return Err(KernelError::InvalidResult(format!(
                "boolean: a contact's pieces disagree on their side near {}",
                show(e.2)
            )));
        }
        e.0.push(c.holder);
    }
    let places: Vec<(Place, (Vec<usize>, i8, DVec3))> = places.into_iter().collect();
    let mut of_seg: FxHashMap<usize, usize> = FxHashMap::default();
    let mut of_point: FxHashMap<usize, usize> = FxHashMap::default();
    for (k, (place, _)) in places.iter().enumerate() {
        match *place {
            Place::Seg(s) => of_seg.insert(s, k),
            Place::Point(c) => of_point.insert(c, k),
        };
    }
    // The pieces of the other operand at each place, in one pass over the
    // loops.
    let mut incident: Vec<Vec<usize>> = vec![Vec::new(); places.len()];
    for (i, p) in pieces.iter().enumerate() {
        for &(s, fw) in p.loops.iter().flatten() {
            for k in [of_seg.get(&s), of_point.get(&segs.start(s, fw))]
                .into_iter()
                .flatten()
            {
                if pieces[places[*k].1 .0[0]].operand != p.operand {
                    incident[*k].push(i);
                }
            }
        }
    }
    let mut pinch: Option<usize> = None;
    for (k, (_, (holders, sigma, at))) in places.iter().enumerate() {
        let inc = &mut incident[k];
        inc.sort_unstable();
        inc.dedup();
        let Some(&first) = inc.first() else {
            return Err(KernelError::InvalidResult(format!(
                "boolean: a contact with no piece of the other body at it near {}",
                show(*at)
            )));
        };
        let class = pieces[first].class;
        let side = match class {
            Class::Outside => 1,
            Class::Inside => -1,
            _ => 0,
        };
        if side != *sigma || inc.iter().any(|&i| pieces[i].class != class) {
            return Err(KernelError::InvalidResult(format!(
                "boolean: a contact's pieces disagree on their side near {}",
                show(*at)
            )));
        }
        if pinch.is_none() && holders.iter().any(|&i| keep(i)) && inc.iter().any(|&i| keep(i)) {
            pinch = Some(k);
        }
    }
    match pinch {
        None => Ok(()),
        Some(k) => {
            let (place, (_, _, at)) = &places[k];
            Err(KernelError::InvalidInput(match place {
                Place::Seg(_) => format!(
                    "the result would be non-manifold: the bodies touch along an edge near {}",
                    show(*at)
                ),
                Place::Point(_) => format!(
                    "the result would be non-manifold: the bodies touch at a point near {}",
                    show(*at)
                ),
            }))
        }
    }
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
    contacts: &[Contact],
    coplanar: &FxHashSet<(usize, usize)>,
    op: BooleanOp,
    tol: &Tolerance,
) -> Result<(Brep, Origins), KernelError> {
    let eps = tol.linear;
    check_contacts(op, &pieces, contacts, &segs)?;
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
    if flat {
        check_apart(pair, &faces, pts, &segs, eps)?;
    }
    let solids = if flat {
        Some(solids(pair, pts, &segs, &faces, &shells)?)
    } else {
        None
    };
    let (mut brep, org, shell_ids) = build(geometry, pts, &segs, &faces, &shells)?;
    check_axes(&brep, tol)?;
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

/// Every place of a flat result keeps clear of the faces it is not on
/// (design note 0012 §19): each end of a segment by more than the
/// tolerance from every face without a corner there, and its middle from
/// every face it does not bound, else the result would touch itself and is
/// refused as non-manifold there. Two faces are let off for the middle of
/// a short segment: one with a corner at one end of a segment shorter
/// than four tolerances (a sharp corner there: an edge leaving a corner of
/// a face at a slight angle stays within the tolerance of the face for the
/// tolerance over the angle's sine), and one with a corner at either end
/// of a segment no longer than two tolerances, whose middle is within the
/// tolerance of both ends. A face with corners at both ends of a longer
/// segment and its middle on it holds the segment lying on it. The rule
/// after selection (`check_contacts`) decides the contacts the imprint
/// records;
/// this decides, from the result alone, those that decisions at the edge
/// of the tolerance let through (found by the contact fuzz, on main too):
/// an edge of one body a tolerance from a face of the other, its vertices
/// called off the face's plane while the points where it leaves the face
/// are within the tolerance of the face's edges and placed on them, so the
/// result has the edge lying on the face; and a piece of an edge of one
/// body a few tolerances long, from a corner of the other within the
/// tolerance of the edge to where the edge crosses a face at that corner,
/// lying on that face. One query of a tree of the
/// faces' boxes per point; a face within the tolerance of a point's plane
/// is tested by parity of its loops in its plane.
fn check_apart(
    pair: &Pair,
    faces: &[Out],
    pts: &Points,
    segs: &Segments,
    eps: f64,
) -> Result<(), KernelError> {
    let planes: Vec<Frame> = faces.iter().map(|f| plane_of(pair, f)).collect();
    let polys: Vec<Vec<Vec<DVec2>>> = (faces.iter().zip(&planes))
        .map(|(f, fr)| {
            (f.loops.iter())
                .map(|l| {
                    corners(segs, l)
                        .into_iter()
                        .map(|c| to2(fr, pts.rep[c]))
                        .collect()
                })
                .collect()
        })
        .collect();
    let boxes: Vec<(DVec3, DVec3)> = (faces.iter())
        .map(|f| {
            directed(f).fold(
                (DVec3::splat(f64::INFINITY), DVec3::splat(f64::NEG_INFINITY)),
                |(lo, hi), (s, fw)| {
                    let p = pts.rep[segs.start(s, fw)];
                    (lo.min(p), hi.max(p))
                },
            )
        })
        .map(|(lo, hi)| (lo - DVec3::splat(eps), hi + DVec3::splat(eps)))
        .collect();
    let tree = crate::bvh::Bvh::new(boxes);
    let mut along: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
    let mut at: FxHashMap<usize, Vec<usize>> = FxHashMap::default();
    for (i, f) in faces.iter().enumerate() {
        for (s, fw) in directed(f) {
            along.entry(s).or_default().push(i);
            at.entry(segs.start(s, fw)).or_default().push(i);
        }
    }
    // Whether `p` is on face `g`: within the tolerance of its plane, and
    // inside its loops or within the tolerance of them.
    let on = |g: usize, p: DVec3| {
        let fr = &planes[g];
        if (p - fr.origin).dot(fr.z).abs() > eps {
            return false;
        }
        let q = to2(fr, p);
        let mut inside = false;
        for l in &polys[g] {
            for k in 0..l.len() {
                let (a, b) = (l[k], l[(k + 1) % l.len()]);
                if (a.y > q.y) != (b.y > q.y) && q.x < a.x + (q.y - a.y) / (b.y - a.y) * (b.x - a.x)
                {
                    inside = !inside;
                }
                let d = b - a;
                let t =
                    ((q - a).dot(d) / d.length_squared().max(f64::MIN_POSITIVE)).clamp(0.0, 1.0);
                if (q - (a + d * t)).length() <= eps {
                    return true;
                }
            }
        }
        inside
    };
    let mut near = Vec::new();
    for (&s, by) in &along {
        let seg = &segs.list[s];
        let (a, b) = (pts.rep[seg.a], pts.rep[seg.b]);
        let ends = |c: usize| at.get(&c).map_or(&[][..], Vec::as_slice);
        let len = (b - a).length();
        // Whether the middle may be near face `g` at a corner at the ends.
        let let_off = |g: usize| {
            let (x, y) = (ends(seg.a).contains(&g), ends(seg.b).contains(&g));
            (len <= 2.0 * eps && (x || y)) || (len < 4.0 * eps && x != y)
        };
        let samples = [
            (a, ends(seg.a), true),
            (b, ends(seg.b), true),
            ((a + b) * 0.5, by.as_slice(), false),
        ];
        for (p, allowed, point) in samples {
            tree.overlapping(p, p, &mut near);
            near.sort_unstable();
            for &g in &near {
                let skip = allowed.contains(&g) || (!point && let_off(g));
                if skip || !on(g, p) {
                    continue;
                }
                return Err(KernelError::InvalidInput(if point {
                    format!(
                        "the result would be non-manifold: the bodies touch at a point near {}",
                        show(p)
                    )
                } else {
                    format!(
                        "the result would be non-manifold: the bodies touch along an edge near {}",
                        show(p)
                    )
                }));
            }
        }
    }
    Ok(())
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
        // Every edge passed has a side, so the split has no bridge (one
        // would be an inconsistency) and no isolated point.
        let split =
            planar::regions(&edges, &|c| chart.map(pts.rep[c]), &reach, eps, &[]).map_err(|p| {
                let (m, c) = p.describe();
                match p {
                    planar::Problem::Far(_) => {
                        KernelError::InvalidInput(format!("{m} near {}", super::show(pts.rep[c])))
                    }
                    _ => KernelError::InvalidResult(format!(
                        "boolean: merging faces: {m} near {}",
                        super::show(pts.rep[c])
                    )),
                }
            })?;
        debug_assert!(split.loose.is_empty(), "a merged face with loose items");
        for region in split.regions {
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
        let mut joins: Vec<Join> = Vec::new();
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
            let along: BTreeSet<(Operand, usize)> = (seg_faces[&s0].iter())
                .flat_map(|&i| faces[i].sources.iter().copied())
                .collect();
            let Some(sources) = join_sources(geometry, pts, segs, &along, [s0, s1], [u, v, w], eps)
            else {
                continue;
            };
            if [u, v, w].iter().any(|c| touched.contains(c)) {
                again.insert(v);
                continue;
            }
            touched.extend([u, v, w]);
            joins.push(Join {
                ends: [u, v, w],
                segs: [s0, s1],
                course: joined,
                sources,
            });
        }
        if joins.is_empty() {
            break;
        }
        for Join {
            ends: [u, v, w],
            segs: [s0, s1],
            course,
            sources,
        } in joins
        {
            let courses = geometry.courses();
            let mut add = |x: Source| match course {
                Course::Straight => segs.add(u, w, Course::Straight, x),
                _ => segs.add_curved(u, w, course, x, &courses, eps),
            };
            let mut sources = sources.into_iter();
            let first = sources.next().expect("a join has a source along all of it");
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

/// A join of segments `segs[0]`, from `ends[0]` to `ends[1]`, and
/// `segs[1]`, from `ends[1]` to `ends[2]`, end to end along `course`, with
/// its sources.
struct Join {
    ends: [usize; 3],
    segs: [usize; 2],
    course: Course,
    sources: Vec<Source>,
}

/// The sources of the segment joining segment `s0`, from `u` to `v`, and
/// segment `s1`, from `v` to `w`, end to end: those that run along all of
/// it, as every origin of an entity does (ADR 0020). `None` if none of
/// them bounds the faces `along` it (ADR 0024): the point then stays a
/// vertex.
///
/// A source of both pieces runs along all of the join. A source of one
/// piece does if the other piece's far end and middle are within the
/// tolerance of it: of a straight operand edge's line, or of the surfaces
/// of a curved edge's faces or of the faces that meet. A straight join
/// allows a bend of the tolerance where the pieces meet, so without this an
/// operand's edge along a piece a few tolerances long, meeting a long piece
/// that leaves its line at a slight angle there, was named the long edge's
/// origin.
fn join_sources(
    geometry: &Geometry,
    pts: &Points,
    segs: &Segments,
    along: &BTreeSet<(Operand, usize)>,
    [s0, s1]: [usize; 2],
    [u, v, w]: [usize; 3],
    eps: f64,
) -> Option<Vec<Source>> {
    let pair = geometry.pair;
    let courses = geometry.courses();
    let (g0, g1) = (&segs.list[s0], &segs.list[s1]);
    // The far end and the middle of the piece `g`, which runs from `v`.
    let probes = |g: &super::Seg, far: usize| {
        let mid = match g.course {
            Course::Straight => (pts.rep[v] + pts.rep[far]) * 0.5,
            c => courses.middle(c),
        };
        [pts.rep[far], mid]
    };
    let mut sources: BTreeSet<Source> = BTreeSet::new();
    for (g, h, far) in [(g0, g1, w), (g1, g0, u)] {
        let mut at: Option<[DVec3; 2]> = None;
        for &x in &g.sources {
            if sources.contains(&x) {
                continue;
            }
            let on = h.sources.contains(&x) || {
                let p = *at.get_or_insert_with(|| probes(h, far));
                p.iter().all(|&q| source_distance(pair, x, q) <= eps)
            };
            if on {
                sources.insert(x);
            }
        }
    }
    (sources.iter().any(|&x| bounds(pair, along, x))).then(|| sources.into_iter().collect())
}

/// How far `p` is from what source `x` runs along: a straight operand
/// edge's line, the surfaces of a curved edge's faces, or those of the
/// faces that meet.
fn source_distance(pair: &Pair, x: Source, p: DVec3) -> f64 {
    let surface = |op: Operand, f: usize| {
        let shape = pair.get(op);
        shape.surfaces[shape.faces[f].surface].distance(p).abs()
    };
    match x {
        Source::Edge(op, e) => {
            let shape = pair.get(op);
            if let (Curve3::Line { .. }, _, _) = shape.curve(e) {
                let (a, b) = shape.edge_points(e);
                (p - a).reject_from(b - a).length()
            } else {
                (shape.edge_faces(e).iter())
                    .map(|&(f, _)| surface(op, f))
                    .fold(0.0, f64::max)
            }
        }
        Source::Meet(fa, fb) => surface(Operand::First, fa).max(surface(Operand::Second, fb)),
    }
}

/// Whether one of the operand faces in `along` (the sources of the
/// result's faces along a segment) bounds what source `x` runs along: the
/// operand's edge, or both faces that meet (ADR 0024).
fn bounds(pair: &Pair, along: &BTreeSet<(Operand, usize)>, x: Source) -> bool {
    match x {
        Source::Edge(op, e) => {
            (pair.get(op).edge_faces(e).iter()).any(|&(g, _)| along.contains(&(op, g)))
        }
        Source::Meet(fa, fb) => {
            along.contains(&(Operand::First, fa)) && along.contains(&(Operand::Second, fb))
        }
    }
}

/// Refuses a result with an edge passing near the axis of its face's
/// surface of revolution between its ends, within the tolerance the
/// validator holds edges to: a flat cap's centre, a cone's apex or a
/// sphere's pole, where the surface is singular and a face keeps a vertex
/// (design note 0012 §15). A plane a few tolerances off a turned rod's
/// axis cuts its ends so; a plane through the axis meets it there, where
/// the edges keep a vertex.
fn check_axes(brep: &Brep, tol: &Tolerance) -> Result<(), KernelError> {
    let near = check::on_tolerance(tol);
    for f in &brep.faces {
        let Surface::Revolution { frame, .. } = &f.surface else {
            continue;
        };
        for c in f.loops.iter().flat_map(|l| &l.coedges) {
            let e = brep.edge(c.edge);
            if let Some(p) = check::through_axis(&e.curve, e.t0, e.t1, frame, near) {
                return Err(KernelError::Unsupported(format!(
                    "booleans whose result has an edge passing within the tolerance of a turned \
                     face's centre, a cone's apex or a sphere's pole (design note 0012 §15), \
                     near {}",
                    show(p)
                )));
            }
        }
    }
    Ok(())
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
/// each face's area times its plane's distance along its outward normal
/// from a corner of the shell. The distance is taken from that corner, not
/// from the world's origin: the faces' areas are of the clusters' places
/// projected onto each plane, which lie within the tolerance of it, so the
/// projected faces fail to close by up to the tolerance times their
/// perimeter, and the sum is off by that times the distance. From the
/// world's origin, that gave the wrong sign to a shell a few tolerances
/// across a metre from it (found by the contact fuzz: a corner dipping four
/// tolerances into a face). The planes, unlike the places, hold the volume
/// of a sliver a tolerance or so thick between two faces with one loop of
/// places: a fan over the places, the other formula tried, gave it none
/// (found by the fuzz of motions by fractions of the tolerance, run a
/// hundred times longer).
fn shell_volume(pair: &Pair, pts: &Points, segs: &Segments, faces: &[Out], shell: &[usize]) -> f64 {
    let (s, fw) = faces[shell[0]].loops[0][0];
    let c = pts.rep[segs.start(s, fw)];
    let mut v = 0.0;
    for &i in shell {
        let f = &faces[i];
        let fr = plane_of(pair, f);
        let area: f64 = (f.loops.iter())
            .map(|l| {
                let q: Vec<DVec2> = (corners(segs, l).into_iter())
                    .map(|k| to2(&fr, pts.rep[k]))
                    .collect();
                (0..q.len())
                    .map(|k| q[k].perp_dot(q[(k + 1) % q.len()]))
                    .sum::<f64>()
                    * 0.5
            })
            .sum();
        v += (fr.origin - c).dot(fr.z) * area / 3.0;
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
    // An origin names only operand entities that bound the result entity
    // (ADR 0024, amending ADR 0020): the sources of the result's faces
    // along each segment and at each cluster, merged faces and coplanar
    // partners included.
    let mut along: FxHashMap<usize, BTreeSet<(Operand, usize)>> = FxHashMap::default();
    let mut around: FxHashMap<usize, BTreeSet<(Operand, usize)>> = FxHashMap::default();
    for f in faces {
        for (s, fw) in directed(f) {
            along.entry(s).or_default().extend(&f.sources);
            around
                .entry(segs.start(s, fw))
                .or_default()
                .extend(&f.sources);
        }
    }
    // Whether a face in `u` holds vertex `v` of `op`.
    let holds = |u: &BTreeSet<(Operand, usize)>, op: Operand, v: usize| {
        (pair.get(op).edges_at(v).iter()).any(|&e| bounds(pair, u, Source::Edge(op, e)))
    };
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
                        let u = &around[&c];
                        let fv = fv.filter(|&v| holds(u, Operand::First, v));
                        let sv = sv.filter(|&v| holds(u, Operand::Second, v));
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
                        let u = &along[&s];
                        let mut src = (seg.sources.iter().copied())
                            .filter(|&x| bounds(pair, u, x))
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
    if vols.iter().all(|v| *v > 0.0) {
        return Ok(solids);
    }
    // Each outer shell's boundary, built once for all the voids.
    let shapes = solids
        .iter()
        .map(|(s, _)| shape::shape(body(*s), "result"))
        .collect::<Result<Vec<_>, _>>()?;
    let outlines = shapes
        .iter()
        .map(imprint::outlines)
        .collect::<Result<Vec<_>, _>>()?;
    let boundaries = shapes
        .iter()
        .zip(&outlines)
        .map(|(shape, outlines)| classify::Boundary::new(shape, outlines))
        .collect::<Result<Vec<_>, _>>()?;
    for v in (0..shells.len()).filter(|&i| vols[i] < 0.0) {
        // A corner of the void shell itself (its body shares the whole
        // result's vertices and edges).
        let shell = &b.shells[shells[v].idx()];
        let first = b.faces[shell.faces[0].idx()].loops[0].coedges[0].edge;
        let q = b.vertices[b.edges[first.idx()].v0.idx()].point;
        let mut best: Option<usize> = None;
        for (k, (s, _)) in solids.iter().enumerate() {
            if boundaries[k].classify(q, eps)? == Class::Inside
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

#[cfg(test)]
mod tests {
    use super::*;

    /// Two triangles' segments: the loop of a piece of the other operand
    /// through clusters 10, 11 and 12 (segments 0, 1 and 2), and the loop
    /// of a holder through clusters 20, 21 and 22.
    fn triangles() -> (Segments, Vec<(usize, bool)>, Vec<(usize, bool)>) {
        let mut segs = Segments::default();
        let mut cycle = |c: [usize; 3]| -> Vec<(usize, bool)> {
            (0..3)
                .map(|i| {
                    let (a, b) = (c[i], c[(i + 1) % 3]);
                    let s = segs.add(a, b, Course::Straight, Source::Edge(Operand::Second, a));
                    (s, segs.start(s, true) == a)
                })
                .collect()
        };
        let (l, h) = (cycle([10, 11, 12]), cycle([20, 21, 22]));
        (segs, l, h)
    }

    fn piece(operand: Operand, class: Class, loops: Vec<Vec<(usize, bool)>>) -> Piece {
        Piece {
            operand,
            face: 0,
            loops,
            class,
            partner: None,
        }
    }

    fn other(op: Operand) -> Operand {
        super::super::other(op)
    }

    /// The four ways of design note 0012 §18.1: a piece of X holding a
    /// contact, of class τ against Y, Y's touching boundary on side σ of
    /// it, so Y's pieces at the contact are outside X for σ = +1 and inside
    /// for σ = −1. Each pinches in exactly one of X ∪ Y, X ∩ Y, X − Y and
    /// Y − X, whichever operand X is, at a segment or at a point.
    #[test]
    fn each_way_of_touching_pinches_in_one_operation() {
        use BooleanOp::*;
        let (segs, l, h) = triangles();
        for (sigma, tau, pinches) in [
            (1, Class::Outside, "X ∪ Y"),
            (1, Class::Inside, "Y − X"),
            (-1, Class::Outside, "X − Y"),
            (-1, Class::Inside, "X ∩ Y"),
        ] {
            let y_class = if sigma > 0 {
                Class::Outside
            } else {
                Class::Inside
            };
            for x in [Operand::First, Operand::Second] {
                let pieces = [
                    piece(x, tau, vec![h.clone()]),
                    piece(other(x), y_class, vec![l.clone()]),
                ];
                let first = x == Operand::First;
                for (place, text) in [
                    (
                        Place::Seg(l[0].0),
                        "along an edge near (1000.000000, 0.000000, 0.000000) mm",
                    ),
                    (
                        Place::Point(11),
                        "at a point near (1000.000000, 0.000000, 0.000000) mm",
                    ),
                ] {
                    let contacts = [Contact {
                        place,
                        holder: 0,
                        sigma,
                        at: DVec3::X,
                    }];
                    for (op, name) in [
                        (Union, "X ∪ Y"),
                        (Intersect, "X ∩ Y"),
                        (Subtract, if first { "X − Y" } else { "Y − X" }),
                    ] {
                        let got = check_contacts(op, &pieces, &contacts, &segs);
                        if name == pinches {
                            let Err(KernelError::InvalidInput(m)) = got else {
                                panic!("{name}, σ {sigma}, τ {tau:?}: {got:?}");
                            };
                            assert_eq!(
                                m,
                                format!(
                                    "the result would be non-manifold: the bodies touch {text}"
                                )
                            );
                        } else {
                            assert!(got.is_ok(), "{name}, σ {sigma}, τ {tau:?}: {got:?}");
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn contacts_are_cross_checked() {
        let (segs, l, h) = triangles();
        let at = |place, holder, sigma| Contact {
            place,
            holder,
            sigma,
            at: DVec3::ZERO,
        };
        let internal = |r: Result<(), KernelError>, what: &str| match r {
            Err(KernelError::InvalidResult(m)) => assert!(m.contains(what), "{m}"),
            other => panic!("expected {what}: {other:?}"),
        };
        let holder = || piece(Operand::First, Class::Outside, vec![h.clone()]);
        let union = BooleanOp::Union;
        // No piece of the other operand at the place: a piece of the
        // holder's own operand does not count.
        let alone = [
            holder(),
            piece(Operand::First, Class::Outside, vec![l.clone()]),
        ];
        internal(
            check_contacts(union, &alone, &[at(Place::Point(11), 0, 1)], &segs),
            "no piece of the other body",
        );
        // An incident piece inside the holder's body where the touching
        // boundary lies outside it, or on a coincident face.
        for class in [Class::Inside, Class::OnSame, Class::OnOpposite] {
            let wrong = [holder(), piece(Operand::Second, class, vec![l.clone()])];
            internal(
                check_contacts(union, &wrong, &[at(Place::Seg(l[0].0), 0, 1)], &segs),
                "disagree on their side",
            );
        }
        // Incident pieces of different classes.
        let mixed = [
            holder(),
            piece(Operand::Second, Class::Outside, vec![l.clone()]),
            piece(Operand::Second, Class::Inside, vec![l.clone()]),
        ];
        internal(
            check_contacts(union, &mixed, &[at(Place::Point(10), 0, 1)], &segs),
            "disagree on their side",
        );
        // Two holders of one place on different sides.
        let two = [
            holder(),
            holder(),
            piece(Operand::Second, Class::Outside, vec![l.clone()]),
        ];
        internal(
            check_contacts(
                union,
                &two,
                &[at(Place::Point(10), 0, 1), at(Place::Point(10), 1, -1)],
                &segs,
            ),
            "disagree on their side",
        );
    }

    /// The least pinching place is reported: segments before points, then
    /// by index; a place that does not pinch is passed over.
    #[test]
    fn the_least_pinching_place_is_reported() {
        let (segs, l, h) = triangles();
        let pieces = [
            piece(Operand::First, Class::Outside, vec![h.clone()]),
            piece(Operand::Second, Class::Outside, vec![l.clone()]),
            piece(Operand::First, Class::Inside, vec![h.clone()]),
        ];
        let at = |place, holder, x: f64| Contact {
            place,
            holder,
            sigma: 1,
            at: DVec3::new(x, 0.0, 0.0),
        };
        let report =
            |contacts: &[Contact]| match check_contacts(BooleanOp::Union, &pieces, contacts, &segs)
            {
                Err(KernelError::InvalidInput(m)) => m,
                other => panic!("{other:?}"),
            };
        let (s0, s2) = (l[0].0, l[2].0);
        // Point 10 is held by an inside piece, which a union drops.
        let m = report(&[
            at(Place::Point(10), 2, 0.001),
            at(Place::Point(12), 0, 0.002),
            at(Place::Seg(s2), 0, 0.003),
            at(Place::Seg(s0), 2, 0.004),
        ]);
        assert!(m.contains("along an edge near (3.000000, "), "{m}");
        let m = report(&[
            at(Place::Point(12), 0, 0.002),
            at(Place::Point(11), 0, 0.001),
        ]);
        assert!(m.contains("at a point near (1.000000, "), "{m}");
    }
}
