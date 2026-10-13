//! Selecting, sewing and simplifying the pieces into the result (design
//! note 0012 §4.4 and §4.5), and its origins (ADR 0020).

use std::collections::{BTreeMap, BTreeSet, HashMap};

use glam::{DVec2, DVec3};
use spade::{ConstrainedDelaunayTriangulation, Point2, Triangulation};

use super::classify::Class;
use super::points::{segment_distance, Points};
use super::{key, planar, show, FaceSeg, Key, Pair, Piece, Segments, Source, SAME_PLANE};
use crate::geom::Curve3;
use crate::math::Frame;
use crate::origin::{
    EdgeIx, EdgeOrigin, FaceIx, FaceOrigin, Merged, Operand, Origins, VertexIx, VertexOrigin,
};
use crate::tess;
use crate::topo::{Brep, Coedge, EdgeId, Loop, ShellId, VertexId};
use crate::{BooleanOp, KernelError, Tolerance};

/// A face of the result being assembled.
struct Out {
    /// The operand face whose surface it lies on (the least of its
    /// constituents).
    operand: Operand,
    face: usize,
    /// True if its outward normal is the opposite of that face's.
    reversed: bool,
    /// Cycles of clusters, outer first, counter-clockwise about the
    /// outward normal.
    loops: Vec<Vec<usize>>,
    /// Every operand face it is made of or coincides with.
    sources: BTreeSet<(Operand, usize)>,
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
    pair: &Pair,
    pts: &Points,
    segs: &Segments,
    pieces: Vec<Piece>,
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
            for l in &mut loops {
                l.reverse();
            }
        }
        faces.push(Out {
            operand: p.operand,
            face: p.face,
            reversed,
            loops,
            sources,
        });
    }
    if faces.is_empty() {
        return Err(KernelError::InvalidInput(
            "the result would be empty".into(),
        ));
    }
    check_closed(&faces, pts)?;
    let faces = merge_faces(pair, pts, faces, eps)?;
    let mut edge_sources: HashMap<Key, BTreeSet<Source>> = segs
        .sources
        .iter()
        .map(|(k, v)| (*k, v.iter().copied().collect()))
        .collect();
    let faces = merge_edges(faces, pts, &mut edge_sources, eps);
    check_closed(&faces, pts)?;
    check_vertices(&faces, pts)?;
    let shells = shells(&faces);
    let solids = solids(pair, pts, &faces, &shells)?;
    build(pair, pts, &faces, &shells, &solids, &edge_sources)
}

/// The directed segments of a face's loops.
fn directed(f: &Out) -> impl Iterator<Item = (usize, usize)> + '_ {
    f.loops
        .iter()
        .flat_map(|l| (0..l.len()).map(move |i| (l[i], l[(i + 1) % l.len()])))
}

/// Every segment must be used once in each direction: more uses are a
/// non-manifold result, fewer an inconsistency.
fn check_closed(faces: &[Out], pts: &Points) -> Result<(), KernelError> {
    let mut uses: BTreeMap<Key, (usize, usize)> = BTreeMap::new();
    for f in faces {
        for (u, v) in directed(f) {
            let e = uses.entry(key(u, v)).or_default();
            if u < v {
                e.0 += 1;
            } else {
                e.1 += 1;
            }
        }
    }
    for (k, (fw, bw)) in uses {
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
fn merge_faces(
    pair: &Pair,
    pts: &Points,
    faces: Vec<Out>,
    eps: f64,
) -> Result<Vec<Out>, KernelError> {
    let n = faces.len();
    let mut parent: Vec<usize> = (0..n).collect();
    let frames: Vec<Frame> = faces.iter().map(|f| plane_of(pair, f)).collect();
    let same = |i: usize, j: usize| {
        let (a, b) = (&frames[i], &frames[j]);
        a.z.dot(b.z) >= 1.0 - SAME_PLANE && a.z.dot(b.origin - a.origin).abs() <= eps
    };
    let mut by_seg: BTreeMap<Key, Vec<usize>> = BTreeMap::new();
    for (i, f) in faces.iter().enumerate() {
        for (u, v) in directed(f) {
            by_seg.entry(key(u, v)).or_default().push(i);
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
        let fr = plane_of(pair, primary);
        let (operand, face, reversed) = (primary.operand, primary.face, primary.reversed);
        let sources: BTreeSet<(Operand, usize)> = parts
            .iter()
            .flat_map(|f| f.sources.iter().copied())
            .collect();
        let mut all: BTreeSet<(usize, usize)> = BTreeSet::new();
        for f in &parts {
            all.extend(directed(f));
        }
        let boundary: Vec<FaceSeg> = all
            .iter()
            .filter(|(u, v)| !all.contains(&(*v, *u)))
            .map(|&(u, v)| FaceSeg {
                key: key(u, v),
                side: Some(u < v),
            })
            .collect();
        let regions = planar::regions(&boundary, &|c| to2(&fr, pts.rep[c]), eps).map_err(|p| {
            let (m, c) = p.describe();
            KernelError::InvalidResult(format!(
                "boolean: merging faces: {m} near {}",
                super::show(pts.rep[c])
            ))
        })?;
        for loops in regions {
            out.push(Out {
                operand,
                face,
                reversed,
                loops,
                sources: sources.clone(),
            });
        }
    }
    Ok(out)
}

/// Joins two segments end to end where they meet at a point no other
/// segment reaches, lie on one line and bound the same two faces.
fn merge_edges(
    mut faces: Vec<Out>,
    pts: &Points,
    sources: &mut HashMap<Key, BTreeSet<Source>>,
    eps: f64,
) -> Vec<Out> {
    loop {
        let mut at: BTreeMap<usize, BTreeSet<Key>> = BTreeMap::new();
        let mut seg_faces: HashMap<Key, BTreeSet<usize>> = HashMap::new();
        for (i, f) in faces.iter().enumerate() {
            for (u, v) in directed(f) {
                let k = key(u, v);
                at.entry(u).or_default().insert(k);
                at.entry(v).or_default().insert(k);
                seg_faces.entry(k).or_default().insert(i);
            }
        }
        let mut touched: BTreeSet<usize> = BTreeSet::new();
        let mut joins: Vec<(usize, usize, usize)> = Vec::new();
        for (&v, ks) in &at {
            let [k0, k1] = ks.iter().copied().collect::<Vec<Key>>()[..] else {
                continue;
            };
            let u = if k0.0 == v { k0.1 } else { k0.0 };
            let w = if k1.0 == v { k1.1 } else { k1.0 };
            if u == w
                || seg_faces[&k0] != seg_faces[&k1]
                || segment_distance(pts.rep[v], pts.rep[u], pts.rep[w]) > eps
                || [u, v, w].iter().any(|c| touched.contains(c))
            {
                continue;
            }
            touched.extend([u, v, w]);
            joins.push((u, v, w));
        }
        if joins.is_empty() {
            return faces;
        }
        for (u, v, w) in joins {
            for f in &mut faces {
                for l in &mut f.loops {
                    l.retain(|&c| c != v);
                }
            }
            let mut s = sources.remove(&key(u, v)).unwrap_or_default();
            s.extend(sources.remove(&key(v, w)).unwrap_or_default());
            sources.entry(key(u, w)).or_default().extend(s);
        }
    }
}

/// Around every vertex the faces must form one fan: two fans are bodies
/// touching at a point, which is non-manifold.
fn check_vertices(faces: &[Out], pts: &Points) -> Result<(), KernelError> {
    // Corners: (vertex, incoming segment, outgoing segment).
    type Corner = (usize, (usize, usize), (usize, usize));
    let mut corners: Vec<Corner> = Vec::new();
    for f in faces {
        for l in &f.loops {
            let n = l.len();
            for i in 0..n {
                let (p, v, q) = (l[(i + n - 1) % n], l[i], l[(i + 1) % n]);
                corners.push((v, (p, v), (v, q)));
            }
        }
    }
    let by_in: HashMap<(usize, usize), usize> =
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
            let (v, _, (_, q)) = corners[c];
            match by_in.get(&(q, v)) {
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
    let mut first: HashMap<Key, usize> = HashMap::new();
    for (i, f) in faces.iter().enumerate() {
        for (u, v) in directed(f) {
            match first.get(&key(u, v)) {
                Some(&j) => unite(&mut parent, i, j),
                None => {
                    first.insert(key(u, v), i);
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
fn shell_volume(pair: &Pair, pts: &Points, faces: &[Out], shell: &[usize]) -> f64 {
    let mut v = 0.0;
    for &i in shell {
        let f = &faces[i];
        let fr = plane_of(pair, f);
        let area: f64 = f
            .loops
            .iter()
            .map(|l| {
                let q: Vec<DVec2> = l.iter().map(|&c| to2(&fr, pts.rep[c])).collect();
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
    faces: &[Out],
    shells: &[Vec<usize>],
) -> Result<Vec<(usize, Vec<usize>)>, KernelError> {
    let vols: Vec<f64> = shells
        .iter()
        .map(|s| shell_volume(pair, pts, faces, s))
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
        .map(|(s, _)| shell_triangles(pair, pts, faces, &shells[*s]))
        .collect::<Result<_, _>>()?;
    for v in voids {
        let q = pts.rep[faces[shells[v][0]].loops[0][0]];
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
    faces: &[Out],
    shell: &[usize],
) -> Result<Vec<[DVec3; 3]>, KernelError> {
    let mut out = Vec::new();
    for &i in shell {
        let f = &faces[i];
        let fr = plane_of(pair, f);
        let mut cdt: ConstrainedDelaunayTriangulation<Point2<f64>> =
            ConstrainedDelaunayTriangulation::new();
        let mut at: Vec<usize> = Vec::new();
        for l in &f.loops {
            let mut hs = Vec::with_capacity(l.len());
            for &c in l {
                let q = to2(&fr, pts.rep[c]);
                let h = cdt.insert(Point2::new(q.x, q.y)).map_err(|e| {
                    KernelError::InvalidResult(format!("boolean: triangulation: {e:?}"))
                })?;
                if h.index() == at.len() {
                    at.push(c);
                }
                hs.push(h);
            }
            for k in 0..hs.len() {
                let (a, b) = (hs[k], hs[(k + 1) % hs.len()]);
                if a != b {
                    cdt.add_constraint(a, b);
                }
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
fn build(
    pair: &Pair,
    pts: &Points,
    faces: &[Out],
    shells: &[Vec<usize>],
    solids: &[(usize, Vec<usize>)],
    edge_sources: &HashMap<Key, BTreeSet<Source>>,
) -> Result<(Brep, Origins), KernelError> {
    let mut brep = Brep::default();
    let mut org = Origins::default();
    let mut merged: Vec<Merged> = Vec::new();
    let mut vertex_of: HashMap<usize, VertexId> = HashMap::new();
    let mut edge_of: HashMap<Key, EdgeId> = HashMap::new();
    let mut face_ids = Vec::with_capacity(faces.len());
    for f in faces {
        let mut loops = Vec::with_capacity(f.loops.len());
        for l in &f.loops {
            let mut coedges = Vec::with_capacity(l.len());
            for i in 0..l.len() {
                let (u, v) = (l[i], l[(i + 1) % l.len()]);
                let k = key(u, v);
                for c in [k.0, k.1] {
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
                let edge = match edge_of.get(&k) {
                    Some(&e) => e,
                    None => {
                        let (a, b) = (pts.rep[k.0], pts.rep[k.1]);
                        let id = brep.add_edge(
                            Curve3::Line {
                                origin: a,
                                dir: b - a,
                            },
                            0.0,
                            1.0,
                            vertex_of[&k.0],
                            vertex_of[&k.1],
                        );
                        edge_of.insert(k, id);
                        let mut src = edge_sources.get(&k).into_iter().flatten().copied();
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
                        id
                    }
                };
                coedges.push(Coedge {
                    edge,
                    forward: u == k.0,
                });
            }
            loops.push(Loop { coedges });
        }
        let src = &pair.get(f.operand).brep.faces[f.face];
        let id = brep.add_face(src.surface.clone(), src.forward != f.reversed, loops);
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
    for (outer, voids) in solids {
        let mut list = vec![shell_ids[*outer]];
        list.extend(voids.iter().map(|&v| shell_ids[v]));
        brep.add_solid(list);
    }
    merged.sort();
    org.merged = merged;
    Ok((brep, org))
}
