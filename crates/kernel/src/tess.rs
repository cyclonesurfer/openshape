//! Tolerance-controlled tessellation of B-rep faces.
//!
//! Approach (see `docs/design/0002-tessellation.md`):
//! 1. Every edge is discretized once, so the two faces that share an edge
//!    share exactly the same boundary vertices and the mesh is watertight.
//! 2. Each face's boundary samples are mapped into the face's parameter
//!    space. Singular boundary points (cone apex, sphere pole) become a run
//!    of parameter-space points along the free parameter, so the
//!    parameter-space boundary stays a simple polygon.
//! 3. The parameter-space polygon (outer loop plus holes) is triangulated with
//!    a constrained Delaunay triangulation, in coordinates scaled so the
//!    mapping is roughly isometric (on a cylinder or a cone, squeezed along
//!    the rulings).
//! 4. Triangles are refined, by inserting their parameter-space centroid or
//!    bisecting their longest edge, until every triangle is within the
//!    chord tolerance of the surface (measured at thirteen points and
//!    estimated between them by the quartic through them), faces the way
//!    the surface does, and its normals turn by less than the angle
//!    tolerance. A triangle facing against the surface beside a boundary
//!    segment has the segment bisected and the face meshed again; a sliver
//!    too small to split that faces against it is mended by an edge flip.
//!    A face that cannot be refined, or would pass the vertex budget, is an
//!    error.

use std::f64::consts::FRAC_PI_4;

use glam::{DVec2, DVec3};
use rustc_hash::{FxHashMap, FxHashSet};
use spade::{ConstrainedDelaunayTriangulation, Point2, Triangulation};

use crate::geom::{Curve3, Surface};
use crate::topo::{Brep, Face};
use crate::{KernelError, Tolerance};

/// Tessellation controls.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TessellationOptions {
    /// Maximum distance between the mesh and the exact surface (metres).
    pub chord_tolerance: f64,
    /// Maximum angle between adjacent mesh normals or segment directions (radians).
    pub max_angle: f64,
}

impl Default for TessellationOptions {
    fn default() -> Self {
        Self {
            chord_tolerance: 1e-4,
            max_angle: 15f64.to_radians(),
        }
    }
}

/// Mesh of one face. Positions are shared with neighbouring faces along edges
/// bit-for-bit, so welding by position gives a watertight mesh.
#[derive(Clone, Debug, Default)]
pub(crate) struct FaceTess {
    pub positions: Vec<DVec3>,
    pub normals: Vec<DVec3>,
    pub triangles: Vec<[u32; 3]>,
    /// Slivers too small to split that faced against the surface, mended
    /// by an edge flip (design note 0002 §4; counted for the tests).
    #[cfg_attr(not(test), allow(dead_code))]
    pub unfolded: usize,
}

/// Mesh of a whole B-rep.
#[derive(Clone, Debug, Default)]
pub(crate) struct BrepTess {
    pub faces: Vec<FaceTess>,
    pub edges: Vec<Vec<DVec3>>,
    /// Per edge: true for a seam (see [`seams`]).
    pub seams: Vec<bool>,
}

/// Per edge of `brep`: true where two pieces of one surface meet with the
/// same orientation, as where a periodic face is split (ADR 0008), where
/// one face uses the edge twice, or where a sweep's faces from two profile
/// curves meet along the track of a point at which the curves join with
/// equal first derivatives (pieces of one C¹ curve, as a closed spline
/// cut in two; [`swept_join`]). Such an edge bounds faces but is no crease
/// in the shape, so displays hide it.
///
/// The tests err towards showing an edge: surfaces they cannot prove to
/// continue one another count as different, so a real crease is never
/// hidden.
pub(crate) fn seams(brep: &Brep) -> Vec<bool> {
    let mut uses: Vec<Vec<usize>> = vec![Vec::new(); brep.edges.len()];
    for (fi, f) in brep.faces.iter().enumerate() {
        for lp in &f.loops {
            for c in &lp.coedges {
                uses[c.edge.idx()].push(fi);
            }
        }
    }
    uses.iter()
        .enumerate()
        .map(|(e, u)| match u.as_slice() {
            [a, b] if a == b => true,
            [a, b] => {
                let (fa, fb) = (&brep.faces[*a], &brep.faces[*b]);
                swept_join(brep, e, fa, fb).unwrap_or_else(|| same_surface(fa, fb))
            }
            _ => false,
        })
        .collect()
}

/// For an edge shared by two faces of one sweep (extrusions along one
/// direction, or revolutions about one axis) traced by a point where one
/// face's profile curve ends and the other's starts: `Some(true)` if the
/// curves have equal first derivatives there, so the surfaces continue one
/// another with no crease; `Some(false)` if not. `None` when the faces and
/// the edge are not such a pair (the caller then compares the surfaces).
///
/// Equal parametric derivatives (not just parallel tangents) are what the
/// pieces of one C¹ curve have where it was cut, so a tangent joint between
/// two different curves (a line and an arc, say) stays a visible edge, as
/// tangent edges are elsewhere.
fn swept_join(brep: &Brep, edge: usize, a: &Face, b: &Face) -> Option<bool> {
    if a.forward != b.forward {
        return None;
    }
    let e = &brep.edges[edge];
    let p = brep.vertices[e.v0.idx()].point;
    let ends = |c: &Curve3| match c {
        Curve3::Nurbs(n) => {
            let (t0, t1) = n.domain();
            let d0 = n.derivatives(t0, 1);
            let d1 = n.derivatives(t1, 1);
            Some([(d0[0], d0[1]), (d1[0], d1[1])])
        }
        _ => None,
    };
    match (&a.surface, &b.surface) {
        (
            Surface::Extrusion {
                curve: ca, dir: da, ..
            },
            Surface::Extrusion {
                curve: cb, dir: db, ..
            },
        ) => {
            if !parallel(*da, *db) {
                return None;
            }
            let (ea, eb) = (ends(ca)?, ends(cb)?);
            let axis = da.normalize_or_zero();
            let scale = 1.0 + p.length() + ea[0].0.length().max(eb[0].0.length());
            // On the edge's ruling through `p`.
            let on_edge = |q: DVec3| {
                let r = q - p;
                (r - axis * r.dot(axis)).length() <= SAME * scale
            };
            let joined = |(qa, _): (DVec3, DVec3), (qb, _): (DVec3, DVec3)| {
                on_edge(qa) && on_edge(qb) && near(qa, qb, scale)
            };
            let same = |x: DVec3, y: DVec3| (x - y).length() <= SAME * x.length().max(y.length());
            if joined(ea[1], eb[0]) {
                Some(same(ea[1].1, eb[0].1))
            } else if joined(eb[1], ea[0]) {
                Some(same(eb[1].1, ea[0].1))
            } else {
                None
            }
        }
        (
            Surface::Revolution {
                frame: fa,
                profile: pa,
            },
            Surface::Revolution {
                frame: fb,
                profile: pb,
            },
        ) => {
            let scale = 1.0 + fa.origin.length().max(fb.origin.length()) + p.length();
            let off = fb.origin - fa.origin;
            if !parallel(fa.z, fb.z) || (off - fa.z * off.dot(fa.z)).length() > SAME * scale {
                return None;
            }
            let (ea, eb) = (ends(pa)?, ends(pb)?);
            // A point and derivative in the frame that turns with the
            // meridian plane: (distance from the axis, height), and the
            // derivative's radial, axial and turning parts.
            let meridian = |(q, d): (DVec3, DVec3)| {
                let r = q - fa.origin;
                let h = r.dot(fa.z);
                let radial = r - fa.z * h;
                let rho = radial.length();
                if rho <= SAME * scale {
                    return None;
                }
                let x = radial / rho;
                Some((
                    DVec3::new(rho, h, 0.0),
                    DVec3::new(d.dot(x), d.dot(fa.z), d.dot(fa.z.cross(x))),
                ))
            };
            let (mp, _) = meridian((p, DVec3::ZERO))?;
            let joined = |x: (DVec3, DVec3), y: (DVec3, DVec3)| {
                let (Some((qa, da)), Some((qb, db))) = (meridian(x), meridian(y)) else {
                    return None;
                };
                (near(qa, mp, scale) && near(qb, mp, scale))
                    .then(|| (da - db).length() <= SAME * da.length().max(db.length()))
            };
            joined(ea[1], eb[0]).or_else(|| joined(eb[1], ea[0]))
        }
        _ => None,
    }
}

/// Relative tolerance for comparing surface data built from the same
/// definition by rigid motions (rounding only).
const SAME: f64 = 1e-10;

fn near(a: DVec3, b: DVec3, scale: f64) -> bool {
    (a - b).length() <= SAME * scale
}

/// True if `a` and `b` point the same way up to rounding: the sine of the
/// angle between them at most [`SAME`] (a test on the cosine,
/// `a · b ≥ 1 − SAME`, would pass angles up to 1.4e-5 rad and hide a slight
/// crease).
fn parallel(a: DVec3, b: DVec3) -> bool {
    let (a, b) = (a.normalize_or_zero(), b.normalize_or_zero());
    a.dot(b) > 0.0 && a.cross(b).length() <= SAME
}

/// True if two faces are pieces of one surface with one outward side:
///
/// * planes with the same outward normal and offset;
/// * extrusions of circles (cylinders) along the same direction with equal
///   radii, parallel circle planes and centres on one axis line;
/// * surfaces of revolution about the same axis whose profiles have the same
///   meridian (distance from and height along the axis): lines and splines
///   point for point, as the copies a revolve rotates about its axis, and
///   circles in planes through the axis by centre and radius;
/// * otherwise, identical surface data.
///
/// Except for planes, the faces must also have the same orientation flag;
/// the pieces a sweep makes are rigid copies, so theirs agree.
pub(crate) fn same_surface(a: &Face, b: &Face) -> bool {
    match (&a.surface, &b.surface) {
        (Surface::Plane(fa), Surface::Plane(fb)) => {
            let sign = |f: &Face| if f.forward { 1.0 } else { -1.0 };
            let (na, nb) = (fa.z * sign(a), fb.z * sign(b));
            let scale = 1.0 + fa.origin.length().max(fb.origin.length());
            parallel(na, nb) && na.dot(fb.origin - fa.origin).abs() <= SAME * scale
        }
        _ if a.forward != b.forward => false,
        (
            Surface::Extrusion {
                curve:
                    Curve3::Circle {
                        center: ca,
                        x: xa,
                        y: ya,
                        radius: ra,
                    },
                dir: da,
                ..
            },
            Surface::Extrusion {
                curve:
                    Curve3::Circle {
                        center: cb,
                        x: xb,
                        y: yb,
                        radius: rb,
                    },
                dir: db,
                ..
            },
        ) => {
            let scale = 1.0 + ca.length().max(cb.length()) + ra.max(*rb);
            let axis = da.normalize_or_zero();
            let off = *cb - *ca;
            (ra - rb).abs() <= SAME * scale
                && parallel(*da, *db)
                && xa.cross(*ya).cross(xb.cross(*yb)).length() <= SAME
                && (off - axis * off.dot(axis)).length() <= SAME * scale
        }
        (
            Surface::Revolution {
                frame: fa,
                profile: pa,
            },
            Surface::Revolution {
                frame: fb,
                profile: pb,
            },
        ) => {
            let scale = 1.0 + fa.origin.length().max(fb.origin.length());
            let off = fb.origin - fa.origin;
            if !parallel(fa.z, fb.z) || (off - fa.z * off.dot(fa.z)).length() > SAME * scale {
                return false;
            }
            // Distance from the axis and height along it.
            let meridian = |p: DVec3| {
                let r = p - fa.origin;
                let h = r.dot(fa.z);
                DVec3::new((r - fa.z * h).length(), h, 0.0)
            };
            let same_at = |ts: &[f64]| {
                ts.iter()
                    .all(|&t| near(meridian(pa.point(t)), meridian(pb.point(t)), scale))
            };
            match (pa, pb) {
                (Curve3::Line { .. }, Curve3::Line { .. }) => same_at(&[0.0, 1.0]),
                // Circles in planes through the axis sweep the same torus
                // (or sphere) when their centres have the same meridian and
                // their radii agree, whatever their parameterization.
                (
                    Curve3::Circle {
                        center: ca,
                        x: xa,
                        y: ya,
                        radius: ra,
                    },
                    Curve3::Circle {
                        center: cb,
                        x: xb,
                        y: yb,
                        radius: rb,
                    },
                ) => {
                    let in_meridian_plane = |c: DVec3, n: DVec3| {
                        let r = c - fa.origin;
                        let radial = r - fa.z * r.dot(fa.z);
                        n.dot(fa.z).abs() <= SAME && n.dot(radial).abs() <= SAME * scale
                    };
                    (ra - rb).abs() <= SAME * scale
                        && near(meridian(*ca), meridian(*cb), scale)
                        && in_meridian_plane(*ca, xa.cross(*ya))
                        && in_meridian_plane(*cb, xb.cross(*yb))
                }
                (Curve3::Nurbs(na), Curve3::Nurbs(nb)) => {
                    na.degree() == nb.degree()
                        && na.knots() == nb.knots()
                        && na.weights() == nb.weights()
                        && na
                            .control_points()
                            .iter()
                            .zip(nb.control_points())
                            .all(|(p, q)| near(meridian(*p), meridian(q), scale))
                }
                _ => false,
            }
        }
        (sa, sb) => sa == sb,
    }
}

/// Rounds of boundary refinement before giving up (design note 0002 §5).
const REFINE_ROUNDS: usize = 60;

/// The most vertices a face's mesh, or an edge's samples, may take: work
/// that grows past it is refused with a message, not left to exhaust the
/// machine.
const VERTEX_BUDGET: usize = 2_000_000;

/// Tessellates every face and edge of `brep`. Edges are discretized to the
/// chord and angle tolerances, then refined where a face's boundary needs
/// it (design note 0002 §5); no segment shorter than `tol` is split.
pub(crate) fn tessellate(
    brep: &Brep,
    opts: &TessellationOptions,
    tol: &Tolerance,
) -> Result<BrepTess, KernelError> {
    if !(opts.chord_tolerance > 0.0 && opts.max_angle > 0.0) {
        return Err(KernelError::InvalidInput(
            "tessellation tolerances must be positive".into(),
        ));
    }
    let mut params: Vec<Vec<f64>> = brep
        .edges
        .iter()
        .enumerate()
        .map(|(i, e)| {
            e.curve
                .discretize(
                    e.t0,
                    e.t1,
                    opts.chord_tolerance,
                    opts.max_angle,
                    VERTEX_BUDGET,
                )
                .map_err(|m| KernelError::Tessellation(format!("edge {i}: {m}")))
        })
        .collect::<Result<_, _>>()?;
    let points = |params: &[Vec<f64>]| -> Vec<Vec<DVec3>> {
        brep.edges
            .iter()
            .zip(params)
            .map(|(e, ts)| {
                let mut pts: Vec<DVec3> = ts.iter().map(|&t| e.curve.point(t)).collect();
                // Pin the ends to the vertices so loops close exactly.
                let n = pts.len();
                pts[0] = brep.vertex(e.v0).point;
                pts[n - 1] = brep.vertex(e.v1).point;
                pts
            })
            .collect()
    };
    let cusps = find_cusps(brep, opts.max_angle);
    // Matches the cusps' samples until none is added, returning which
    // edges changed. Within a group every sample keeps its distance along
    // the group's direction, so a distance has crossed the longest chain of
    // cusps after as many passes as there are cusps.
    let rematch = |params: &mut Vec<Vec<f64>>, edges: &mut Vec<Vec<DVec3>>| {
        let mut changed = vec![false; brep.edges.len()];
        for _ in 0..=cusps.len() {
            let pass = match_cusps(brep, &cusps, params, edges, tol.linear);
            if !pass.contains(&true) {
                break;
            }
            *edges = points(params);
            for (c, p) in changed.iter_mut().zip(pass) {
                *c |= p;
            }
        }
        changed
    };
    let mut edges = points(&params);
    rematch(&mut params, &mut edges);
    let face_loops = |i: usize, edges: &[Vec<DVec3>]| {
        loops_uv(&brep.faces[i], edges, opts.max_angle)
            .map_err(|m| KernelError::Tessellation(format!("face {i}: {m}")))
    };
    let mut loops: Vec<Vec<Vec<Sample>>> = (0..brep.faces.len())
        .map(|i| face_loops(i, &edges))
        .collect::<Result<_, _>>()?;
    // Refinement is local; a boundary growing past this many samples is
    // chasing something bisection cannot resolve (the triangulation then
    // reports it).
    let budget = 8 * params.iter().map(Vec::len).sum::<usize>() + 100_000;
    // The cusps of each face, and the cusps beside each edge, by index.
    let mut cusps_of_face: Vec<Vec<usize>> = vec![Vec::new(); brep.faces.len()];
    let mut cusps_of_edge: Vec<Vec<usize>> = vec![Vec::new(); brep.edges.len()];
    for (k, c) in cusps.iter().enumerate() {
        cusps_of_face[c.face].push(k);
        for (e, _) in c.sides {
            if cusps_of_edge[e].last() != Some(&k) {
                cusps_of_edge[e].push(k);
            }
        }
    }
    // Bisects the edge segments `split` names (per edge, the indices of
    // their first samples) longer than the model tolerance, keeping the
    // cusps matched; returns which edges changed.
    let bisect = |split: &mut Vec<Vec<usize>>,
                  params: &mut Vec<Vec<f64>>,
                  edges: &mut Vec<Vec<DVec3>>| {
        let runs: Vec<CuspRuns> = cusps
            .iter()
            .map(|c| cusp_runs(brep, c, params, edges, tol.linear))
            .collect();
        let mut changed = vec![false; brep.edges.len()];
        for (e, segs) in split.iter_mut().enumerate() {
            if segs.is_empty() {
                continue;
            }
            segs.sort_unstable();
            segs.dedup();
            let ts = &params[e];
            let pts = &edges[e];
            let curve = &brep.edges[e].curve;
            let mut out: Vec<f64> = Vec::with_capacity(ts.len() + segs.len());
            let mut next = segs.iter().peekable();
            for i in 0..ts.len() {
                out.push(ts[i]);
                if next.peek() == Some(&&i) {
                    next.next();
                    // In a cusp's matched run, halve the distance along its
                    // direction instead, so that both sides split at the
                    // same abscissa and stay matched.
                    let cusp = cusps_of_edge[e]
                        .iter()
                        .find_map(|&k| runs[k].has((e, i)).then_some(&cusps[k]));
                    let mid = cusp
                        .and_then(|c| {
                            let d = 0.5 * (pts[i] - c.v + pts[i + 1] - c.v).dot(c.dir);
                            solve_monotone(
                                |t| (curve.point(t) - c.v).dot(c.dir) - d,
                                ts[i],
                                ts[i + 1],
                            )
                        })
                        .unwrap_or(0.5 * (ts[i] + ts[i + 1]));
                    if (pts[i + 1] - pts[i]).length() > tol.linear && mid > ts[i] && mid < ts[i + 1]
                    {
                        out.push(mid);
                        changed[e] = true;
                    }
                }
            }
            params[e] = out;
        }
        if changed.contains(&true) {
            *edges = points(params);
            // Bisection near a cusp breaks its matching; restore it.
            for (c, r) in changed.iter_mut().zip(rematch(params, edges)) {
                *c |= r;
            }
        }
        changed
    };
    // The faces whose loops use an edge that changed.
    let touched = |changed: &[bool]| -> Vec<usize> {
        (0..brep.faces.len())
            .filter(|&i| {
                brep.faces[i]
                    .loops
                    .iter()
                    .flat_map(|lp| &lp.coedges)
                    .any(|c| changed[c.edge.idx()])
            })
            .collect()
    };
    for _ in 0..REFINE_ROUNDS {
        if params.iter().map(Vec::len).sum::<usize>() > budget {
            break;
        }
        let runs: Vec<CuspRuns> = cusps
            .iter()
            .map(|c| cusp_runs(brep, c, &params, &edges, tol.linear))
            .collect();
        let mut split: Vec<Vec<usize>> = vec![Vec::new(); brep.edges.len()];
        for (fi, (face, l)) in brep.faces.iter().zip(&loops).enumerate() {
            let matched: Vec<&CuspRuns> = cusps_of_face[fi].iter().map(|&k| &runs[k]).collect();
            for (e, i) in boundary_conflicts(brep, face, l, &params, &matched) {
                split[e].push(i);
            }
        }
        let changed = bisect(&mut split, &mut params, &mut edges);
        if !changed.contains(&true) {
            break;
        }
        for i in touched(&changed) {
            loops[i] = face_loops(i, &edges)?;
        }
    }
    // The faces, each meshed again after an edge of it is bisected where
    // its triangles cannot face the way the surface does (design note 0002
    // §5).
    let mut faces: Vec<FaceTess> = vec![FaceTess::default(); brep.faces.len()];
    let mut todo: Vec<usize> = (0..brep.faces.len()).collect();
    for round in 0.. {
        let mut split: Vec<Vec<usize>> = vec![Vec::new(); brep.edges.len()];
        // The faces that asked, with the edges they asked to bisect.
        let mut asked: Vec<(usize, Vec<usize>)> = Vec::new();
        for &i in &todo {
            match tessellate_face(&brep.faces[i], loops[i].clone(), opts, tol.linear) {
                Ok(t) => faces[i] = t,
                Err(Unmeshed::Split(segs)) => {
                    asked.push((i, segs.iter().map(|&(e, _)| e).collect()));
                    for (e, k) in segs {
                        split[e].push(k);
                    }
                }
                Err(Unmeshed::Failed(m)) => {
                    return Err(KernelError::Tessellation(format!("face {i}: {m}")))
                }
            }
        }
        if asked.is_empty() {
            break;
        }
        let refuse = |i: usize, why: &str| {
            KernelError::Tessellation(format!(
                "face {i}: triangles beside its boundary face against the surface, and {why}"
            ))
        };
        if round == REFINE_ROUNDS || params.iter().map(Vec::len).sum::<usize>() > budget {
            return Err(refuse(
                asked[0].0,
                "bisecting the boundary did not settle them",
            ));
        }
        let changed = bisect(&mut split, &mut params, &mut edges);
        for (i, wanted) in &asked {
            if !wanted.iter().any(|&e| changed[e]) {
                return Err(refuse(*i, "its segments there are too short to bisect"));
            }
        }
        todo = touched(&changed);
        for &i in &todo {
            loops[i] = face_loops(i, &edges)?;
        }
    }
    Ok(BrepTess {
        faces,
        edges,
        seams: seams(brep),
    })
}

/// Why a face was not meshed.
enum Unmeshed {
    /// These edge segments (edge, index of the segment's first sample)
    /// are to be bisected first: triangles beside them face against the
    /// surface, and cannot be made to face along it inside the face.
    Split(Vec<(usize, usize)>),
    /// Refused, for this reason.
    Failed(String),
}

impl From<String> for Unmeshed {
    fn from(m: String) -> Self {
        Unmeshed::Failed(m)
    }
}

/// How a triangle of a face's refinement fares.
enum Judged {
    /// Within the tolerances.
    Fine,
    /// To be split at this point (in the face's parameters).
    Insert(DVec2),
    /// Beside this edge segment (edge, index of its first sample), to be
    /// bisected first.
    Bisect((usize, usize)),
    /// Too small to split, and facing against the surface: to be unfolded
    /// once refinement is done.
    Folded,
}

/// A vertex of a face's loop at which the two boundary curves leave in
/// directions less than the angle tolerance apart: a cusp (design note
/// 0002 §5). Near it the curves run side by side, closer than any
/// deviation bound, so their samples are matched instead.
struct Cusp {
    /// The face.
    face: usize,
    /// The vertex.
    v: DVec3,
    /// Unit direction along which both sides are matched: halfway between
    /// the two curves' directions there, or, for cusps that share a side
    /// (a vertex where the edges of several faces leave together), the
    /// mean direction of all their sides, one for the whole group.
    dir: DVec3,
    /// The two edges leaving the vertex, each with `true` if the vertex is
    /// at the edge's start.
    sides: [(usize, bool); 2],
    /// Per side, `true` if its edge also has a cusp at its other end: its
    /// run then stops at the edge's parameter midpoint.
    half: [bool; 2],
}

/// The cusps of every face of `brep`.
fn find_cusps(brep: &Brep, max_angle: f64) -> Vec<Cusp> {
    let mut out = Vec::new();
    // Each cusp's sides with the directions they leave in.
    let mut leaves: Vec<[((usize, bool), DVec3); 2]> = Vec::new();
    for (fi, f) in brep.faces.iter().enumerate() {
        for lp in &f.loops {
            let n = lp.coedges.len();
            for k in 0..n {
                let prev = &lp.coedges[(k + n - 1) % n];
                let next = &lp.coedges[k];
                let (ep, en) = (brep.edge(prev.edge), brep.edge(next.edge));
                // Each side as (edge, vertex at its start), and the
                // direction in which its curve leaves the vertex.
                let leave = |e: &crate::topo::Edge, at_start: bool| {
                    let d = if at_start {
                        e.curve.deriv(e.t0)
                    } else {
                        -e.curve.deriv(e.t1)
                    };
                    d.normalize_or_zero()
                };
                let sp = (prev.edge.idx(), !prev.forward);
                let sn = (next.edge.idx(), next.forward);
                let (a, b) = (leave(ep, sp.1), leave(en, sn.1));
                if a == DVec3::ZERO || b == DVec3::ZERO || a.angle_between(b) > max_angle {
                    continue;
                }
                let v = if next.forward { en.v0 } else { en.v1 };
                out.push(Cusp {
                    face: fi,
                    v: brep.vertex(v).point,
                    dir: (a + b).normalize(),
                    sides: [sp, sn],
                    half: [false; 2],
                });
                leaves.push([(sp, a), (sn, b)]);
            }
        }
    }
    // Cusps sharing a side are matched along one direction: along their
    // own, the samples one adds to the side they share are at distances
    // the other has not matched, and matching around a vertex where three
    // faces' edges leave together never settles.
    let mut group: Vec<usize> = (0..out.len()).collect();
    fn root(group: &mut [usize], mut k: usize) -> usize {
        while group[k] != k {
            group[k] = group[group[k]];
            k = group[k];
        }
        k
    }
    let mut sides: Vec<((usize, bool), usize, DVec3)> = leaves
        .iter()
        .enumerate()
        .flat_map(|(k, l)| l.map(|(s, d)| (s, k, d)))
        .collect();
    sides.sort_by_key(|&(s, k, _)| (s, k));
    for w in sides.windows(2) {
        if w[0].0 == w[1].0 {
            let (r0, r1) = (root(&mut group, w[0].1), root(&mut group, w[1].1));
            group[r0.max(r1)] = r0.min(r1);
        }
    }
    let mut sum = vec![DVec3::ZERO; out.len()];
    for (i, &(s, k, d)) in sides.iter().enumerate() {
        if i == 0 || sides[i - 1].0 != s {
            sum[root(&mut group, k)] += d;
        }
    }
    // An edge with cusps at both ends is matched from each end to its
    // middle only, so the two ends' matchings never meet.
    let mut ends = vec![0u8; brep.edges.len()];
    for &(s, _, _) in &sides {
        ends[s.0] |= if s.1 { 1 } else { 2 };
    }
    for (k, c) in out.iter_mut().enumerate() {
        let dir = sum[root(&mut group, k)].normalize_or_zero();
        if dir != DVec3::ZERO {
            c.dir = dir;
        }
        c.half = c.sides.map(|(e, _)| ends[e] == 3);
    }
    out
}

/// The samples of one side of a cusp over which its curve runs away from
/// the vertex as a graph over the cusp's direction: (sample index, distance
/// along the direction), outward, the vertex itself excluded; with `half`,
/// only those before the edge's parameter midpoint (from its end: those at
/// or after it).
fn cusp_run(
    brep: &Brep,
    c: &Cusp,
    ((e, at_start), half): ((usize, bool), bool),
    params: &[f64],
    pts: &[DVec3],
) -> Vec<(usize, f64)> {
    let curve = &brep.edges[e].curve;
    let n = pts.len();
    let sign = if at_start { 1.0 } else { -1.0 };
    let mid = 0.5 * (params[0] + params[n - 1]);
    let mut out = Vec::new();
    let mut last = 0.0;
    for k in 1..n {
        let j = if at_start { k } else { n - 1 - k };
        if half && (params[j] >= mid) == at_start {
            break;
        }
        let d = (pts[j] - c.v).dot(c.dir);
        // Each segment turns by less than the angle tolerance, so tangents
        // within 60° of the direction at both ends keep it monotonic.
        let tangent = (curve.deriv(params[j]) * sign).normalize_or_zero();
        if d <= last || tangent.dot(c.dir) < 0.5 {
            break;
        }
        out.push((j, d));
        last = d;
    }
    out
}

/// The segments of each side of a cusp within its matched run: (edge,
/// first segment, end segment), segment `i` running from sample `i` to
/// `i + 1`.
struct CuspRuns {
    sides: [(usize, usize, usize); 2],
}

fn on_side((e, i): (usize, usize), (edge, lo, hi): (usize, usize, usize)) -> bool {
    e == edge && i >= lo && i < hi
}

impl CuspRuns {
    /// True if edge segments `x` and `y` lie on opposite sides of the cusp,
    /// both within its matched run.
    fn pairs(&self, x: (usize, usize), y: (usize, usize)) -> bool {
        let [a, b] = self.sides;
        (on_side(x, a) && on_side(y, b)) || (on_side(x, b) && on_side(y, a))
    }

    /// True if edge segment `x` lies within the matched run of either side.
    fn has(&self, x: (usize, usize)) -> bool {
        self.sides.iter().any(|&s| on_side(x, s))
    }
}

/// The matched runs of a cusp: each side's samples out to where the
/// shorter of the two runs ends.
fn cusp_runs(brep: &Brep, c: &Cusp, params: &[Vec<f64>], pts: &[Vec<DVec3>], tol: f64) -> CuspRuns {
    let run = |k: usize| {
        let e = c.sides[k].0;
        cusp_run(brep, c, (c.sides[k], c.half[k]), &params[e], &pts[e])
    };
    let (ra, rb) = (run(0), run(1));
    let dmax = match (ra.last(), rb.last()) {
        (Some(a), Some(b)) => a.1.min(b.1) + tol,
        _ => 0.0,
    };
    let side = |(e, at_start): (usize, bool), r: &[(usize, f64)]| {
        let k = r.iter().take_while(|(_, d)| *d <= dmax).count();
        let n = pts[e].len();
        if at_start {
            (e, 0, k)
        } else {
            (e, n - 1 - k, n - 1)
        }
    };
    CuspRuns {
        sides: [side(c.sides[0], &ra), side(c.sides[1], &rb)],
    }
}

/// Gives the two sides of every cusp samples at the same distances along
/// its direction, out to where the shorter run ends (design note 0002 §5):
/// two curves that are graphs over one line, sampled at the same abscissae,
/// have polygons that never cross, however close the curves run. A sample
/// closer than `tol` to an existing one is not added. Returns, per edge,
/// whether its samples changed.
fn match_cusps(
    brep: &Brep,
    cusps: &[Cusp],
    params: &mut [Vec<f64>],
    pts: &[Vec<DVec3>],
    tol: f64,
) -> Vec<bool> {
    let mut add: Vec<Vec<f64>> = vec![Vec::new(); params.len()];
    for c in cusps {
        let runs = [0, 1].map(|k| {
            let e = c.sides[k].0;
            cusp_run(brep, c, (c.sides[k], c.half[k]), &params[e], &pts[e])
        });
        let (Some(la), Some(lb)) = (runs[0].last(), runs[1].last()) else {
            continue;
        };
        let dmax = la.1.min(lb.1) + tol;
        for (from, to) in [(0, 1), (1, 0)] {
            let (e, at_start) = c.sides[to];
            let curve = &brep.edges[e].curve;
            let ts = &params[e];
            let vertex = if at_start { 0 } else { ts.len() - 1 };
            // The target run with the vertex in front, at distance zero.
            let mut target = vec![(vertex, 0.0)];
            target.extend(runs[to].iter().copied());
            for &(_, d) in runs[from].iter().take_while(|(_, d)| *d <= dmax) {
                let k = target.partition_point(|(_, x)| *x < d);
                if k == 0 || k == target.len() {
                    continue;
                }
                let ((j0, d0), (j1, d1)) = (target[k - 1], target[k]);
                if d - d0 <= tol || d1 - d <= tol {
                    continue;
                }
                let f = |t: f64| (curve.point(t) - c.v).dot(c.dir) - d;
                if let Some(t) = solve_monotone(f, ts[j0], ts[j1]) {
                    add[e].push(t);
                }
            }
        }
    }
    let mut changed = vec![false; params.len()];
    for (e, mut extra) in add.into_iter().enumerate() {
        if extra.is_empty() {
            continue;
        }
        extra.sort_by(f64::total_cmp);
        let curve = &brep.edges[e].curve;
        let mut out: Vec<f64> = Vec::with_capacity(params[e].len() + extra.len());
        let old = &params[e];
        let mut next = extra.into_iter().peekable();
        for (i, &t) in old.iter().enumerate() {
            if i > 0 {
                while let Some(&x) = next.peek() {
                    if x >= t {
                        break;
                    }
                    next.next();
                    let prev = *out.last().unwrap_or(&old[0]);
                    let p = curve.point(x);
                    if x > prev
                        && (p - curve.point(prev)).length() > tol
                        && (curve.point(t) - p).length() > tol
                    {
                        out.push(x);
                        changed[e] = true;
                    }
                }
            }
            out.push(t);
        }
        params[e] = out;
    }
    changed
}

/// The root of `f` between `a` and `b` (in either order), if `f` changes
/// sign there: safeguarded regula falsi (Illinois), to rounding.
fn solve_monotone(f: impl Fn(f64) -> f64, a: f64, b: f64) -> Option<f64> {
    let (mut a, mut b) = (a.min(b), a.max(b));
    let (mut fa, mut fb) = (f(a), f(b));
    if fa == 0.0 {
        return Some(a);
    }
    if fb == 0.0 {
        return Some(b);
    }
    if fa.signum() == fb.signum() {
        return None;
    }
    let mut side = 0;
    for _ in 0..200 {
        let mut x = (a * fb - b * fa) / (fb - fa);
        if !(x > a && x < b) {
            x = 0.5 * (a + b);
        }
        let fx = f(x);
        if fx == 0.0 || b - a <= 4.0 * f64::EPSILON * a.abs().max(b.abs()) {
            return Some(x);
        }
        if fx.signum() == fa.signum() {
            a = x;
            fa = fx;
            if side == -1 {
                fb *= 0.5;
            }
            side = -1;
        } else {
            b = x;
            fb = fx;
            if side == 1 {
                fa *= 0.5;
            }
            side = 1;
        }
    }
    Some(0.5 * (a + b))
}

/// Edge segments (edge, index of the segment's first point) to bisect so
/// that `face`'s boundary polygon has the topology of its exact boundary
/// (design note 0002 §5), in the face's scaled parameter space:
/// * every segment of a loop with fewer than three samples;
/// * of two segments that cross or touch other than at a shared end, or
///   whose curves may come closer than their chords show (the chords are
///   within the sum of the curves' deviations from them), the one that
///   deviates more (both when within 10%), or failing that the longer;
/// * of two segments sharing an end, one whose curve may pass the other's
///   far end (the far end is within its deviation of its chord).
///
/// A segment's deviation is 1.25 times the largest distance from its chord
/// of the curve's points at a quarter, half and three quarters of the
/// segment, mapped into the parameter space. Two segments on the two sides
/// of a cusp, within its matched runs (`cusps`), conflict only if they
/// cross or touch: the curves there come closer than any deviation bound,
/// and matching keeps their polygons apart. Nor do deviations alone make
/// two segments of one loop conflict across a monotone stretch
/// ([`monotone_between`]).
fn boundary_conflicts(
    brep: &Brep,
    face: &Face,
    loops: &[Vec<Sample>],
    params: &[Vec<f64>],
    cusps: &[&CuspRuns],
) -> Vec<(usize, usize)> {
    let mut out = Vec::new();
    for l in loops {
        if l.len() < 3 {
            out.extend(l.iter().filter_map(|smp| smp.seg));
        }
    }
    if !out.is_empty() {
        return out;
    }
    let s = &face.surface;
    let Ok(scale) = param_scale(s, loops) else {
        return out;
    };
    let is_planar = matches!(s, Surface::Plane(_));
    struct BSeg {
        a: DVec2,
        b: DVec2,
        lo: DVec2,
        hi: DVec2,
        seg: Option<(usize, usize)>,
        chord: f64,
        dev: f64,
        lp: usize,
        idx: usize,
        len: usize,
    }
    let mut segs: Vec<BSeg> = Vec::new();
    for (li, l) in loops.iter().enumerate() {
        for k in 0..l.len() {
            let (p, q) = (l[k], l[(k + 1) % l.len()]);
            let (a, b) = (p.uv * scale, q.uv * scale);
            let dev = match p.seg {
                Some((e, i)) => {
                    let (t0, t1) = (params[e][i], params[e][i + 1]);
                    let c = &brep.edges[e].curve;
                    let mut d = 0.0f64;
                    for f in [0.25, 0.5, 0.75] {
                        let x = c.point(t0 + (t1 - t0) * f);
                        let w = if is_planar {
                            let (u, v) = s.project(x);
                            DVec2::new(u, v)
                        } else {
                            // From the chord's point, so the parameter
                            // stays on the segment's side of any closure
                            // point or seam.
                            project_near(s, x, p.uv.lerp(q.uv, f))
                        };
                        d = d.max(point_segment(w * scale, a, b));
                    }
                    // A smooth segment turning less than the angle
                    // tolerance strays from its chord by a fraction of its
                    // length; more is a projection artefact.
                    (1.25 * d).min(0.5 * (b - a).length())
                }
                None => 0.0,
            };
            segs.push(BSeg {
                a,
                b,
                lo: a.min(b) - DVec2::splat(dev),
                hi: a.max(b) + DVec2::splat(dev),
                seg: p.seg,
                chord: (q.p - p.p).length(),
                dev,
                lp: li,
                idx: k,
                len: l.len(),
            });
        }
    }
    // Of a conflicting pair, the segment(s) to bisect.
    let pick = |x: &BSeg, y: &BSeg, out: &mut Vec<(usize, usize)>| {
        let (kx, ky) = if x.dev + y.dev > 0.0 {
            (x.dev, y.dev)
        } else {
            (x.chord, y.chord)
        };
        if let Some(sx) = x.seg.filter(|_| kx >= 0.9 * ky || y.seg.is_none()) {
            out.push(sx);
        }
        if let Some(sy) = y.seg.filter(|_| ky >= 0.9 * kx || x.seg.is_none()) {
            out.push(sy);
        }
    };
    // The matched runs beside each edge.
    let mut runs_of_edge: FxHashMap<usize, Vec<&CuspRuns>> = FxHashMap::default();
    for &c in cusps {
        for (e, _, _) in c.sides {
            let list = runs_of_edge.entry(e).or_default();
            if !list.iter().any(|r| std::ptr::eq(*r, c)) {
                list.push(c);
            }
        }
    }
    segs.sort_by(|x, y| x.lo.x.total_cmp(&y.lo.x));
    for i in 0..segs.len() {
        let x = &segs[i];
        for y in &segs[i + 1..] {
            if y.lo.x > x.hi.x {
                break;
            }
            if y.lo.y > x.hi.y || x.lo.y > y.hi.y {
                continue;
            }
            let consecutive =
                x.lp == y.lp && ((x.idx + 1) % x.len == y.idx || (y.idx + 1) % y.len == x.idx);
            let shared = [(x.a, y.b), (x.a, y.a), (x.b, y.a), (x.b, y.b)]
                .into_iter()
                .find(|(p, q)| p == q);
            let matched = match (x.seg, y.seg) {
                (Some(sx), Some(sy)) => runs_of_edge
                    .get(&sx.0)
                    .is_some_and(|rs| rs.iter().any(|c| c.pairs(sx, sy))),
                _ => false,
            };
            if matched {
                if !(consecutive || shared.is_some()) && segments_touch(x.a, x.b, y.a, y.b) {
                    pick(x, y, &mut out);
                }
                continue;
            }
            if consecutive || shared.is_some() {
                // The far ends of each, against the other's deviation.
                let far = |z: &BSeg, w: &BSeg| if z.a == w.a || z.a == w.b { z.b } else { z.a };
                if let Some(sx) = x.seg {
                    if x.dev > 0.0 && point_segment(far(y, x), x.a, x.b) <= x.dev {
                        out.push(sx);
                    }
                }
                if let Some(sy) = y.seg {
                    if y.dev > 0.0 && point_segment(far(x, y), y.a, y.b) <= y.dev {
                        out.push(sy);
                    }
                }
                continue;
            }
            if segments_touch(x.a, x.b, y.a, y.b)
                || (segment_distance(x.a, x.b, y.a, y.b) <= x.dev + y.dev
                    && !(x.lp == y.lp && monotone_between(&loops[x.lp], scale, x.idx, y.idx)))
            {
                pick(x, y, &mut out);
            }
        }
    }
    out
}

/// Segments of boundary between two segments of one loop beyond which
/// [`monotone_between`] gives up.
const MONOTONE_SPAN: usize = 32;

/// True if the boundary polygon `l` (scaled by `scale`) runs from segment
/// `i` to segment `j` (the shorter way round, both included, at most
/// [`MONOTONE_SPAN`] segments) within 45° of one direction. Each segment's
/// curve turns by less than the angle tolerance, so the curve along that
/// stretch is a graph over the direction too: two segments of it cannot
/// meet, however close their deviation bands come where short segments
/// lie between them (design note 0002 §5).
fn monotone_between(l: &[Sample], scale: DVec2, i: usize, j: usize) -> bool {
    let n = l.len();
    let (fwd, back) = ((j + n - i) % n, (i + n - j) % n);
    let (start, span) = if fwd <= back { (i, fwd) } else { (j, back) };
    if span > MONOTONE_SPAN {
        return false;
    }
    let at = |k: usize| l[k % n].uv * scale;
    let dir = at(start + span + 1) - at(start);
    (start..=start + span).all(|k| {
        let d = at(k + 1) - at(k);
        d.dot(dir) >= std::f64::consts::FRAC_1_SQRT_2 * d.length() * dir.length()
    })
}

/// The parameters of the point of `s` nearest to `p`, by Gauss–Newton
/// steps from `guess` (enough for points near the guess's surface point).
fn project_near(s: &Surface, p: DVec3, guess: DVec2) -> DVec2 {
    let mut uv = guess;
    for _ in 0..6 {
        let (su, sv) = s.derivs(uv.x, uv.y);
        let r = p - s.point(uv.x, uv.y);
        let (a, b, c) = (su.dot(su), su.dot(sv), sv.dot(sv));
        let det = a * c - b * b;
        let (r0, r1) = (su.dot(r), sv.dot(r));
        let step = if det.abs() > 1e-14 * a.max(c).max(1e-300).powi(2) {
            DVec2::new((c * r0 - b * r1) / det, (a * r1 - b * r0) / det)
        } else if a < c {
            DVec2::new(0.0, r1 / c.max(1e-300))
        } else {
            DVec2::new(r0 / a.max(1e-300), 0.0)
        };
        if !step.is_finite() {
            break;
        }
        uv += step;
        if step.length() <= 1e-12 * (1.0 + uv.length()) {
            break;
        }
    }
    uv
}

/// Distance from `p` to the segment `ab`.
fn point_segment(p: DVec2, a: DVec2, b: DVec2) -> f64 {
    let ab = b - a;
    let l2 = ab.length_squared();
    let s = if l2 > 0.0 {
        ((p - a).dot(ab) / l2).clamp(0.0, 1.0)
    } else {
        0.0
    };
    (a + ab * s - p).length()
}

/// Distance between two segments that do not cross.
fn segment_distance(a: DVec2, b: DVec2, c: DVec2, d: DVec2) -> f64 {
    point_segment(a, c, d)
        .min(point_segment(b, c, d))
        .min(point_segment(c, a, b))
        .min(point_segment(d, a, b))
}

fn orient(a: DVec2, b: DVec2, c: DVec2) -> f64 {
    robust::orient2d(
        robust::Coord { x: a.x, y: a.y },
        robust::Coord { x: b.x, y: b.y },
        robust::Coord { x: c.x, y: c.y },
    )
}

/// True if the closed segments `ab` and `cd` cross or touch (exact
/// orientation tests, Shewchuk).
fn segments_touch(a: DVec2, b: DVec2, c: DVec2, d: DVec2) -> bool {
    let o1 = orient(a, b, c);
    let o2 = orient(a, b, d);
    let o3 = orient(c, d, a);
    let o4 = orient(c, d, b);
    if ((o1 > 0.0 && o2 < 0.0) || (o1 < 0.0 && o2 > 0.0))
        && ((o3 > 0.0 && o4 < 0.0) || (o3 < 0.0 && o4 > 0.0))
    {
        return true;
    }
    let on = |p: DVec2, q: DVec2, r: DVec2| {
        r.x >= p.x.min(q.x) && r.x <= p.x.max(q.x) && r.y >= p.y.min(q.y) && r.y <= p.y.max(q.y)
    };
    (o1 == 0.0 && on(a, b, c))
        || (o2 == 0.0 && on(a, b, d))
        || (o3 == 0.0 && on(c, d, a))
        || (o4 == 0.0 && on(c, d, b))
}

#[derive(Clone, Copy)]
struct Sample {
    uv: DVec2,
    p: DVec3,
    /// The edge segment (edge, index of its first point) from this sample
    /// to the next along the loop; `None` between the two parameter-space
    /// points of a split singular point.
    seg: Option<(usize, usize)>,
}

fn singular_dir(s: &Surface, uv: DVec2) -> Option<bool> {
    let (su, sv) = s.derivs(uv.x, uv.y);
    let (a, b) = (su.length(), sv.length());
    if a < 1e-9 * b.max(1e-300) {
        Some(true) // u is ambiguous
    } else if b < 1e-9 * a.max(1e-300) {
        Some(false) // v is ambiguous
    } else {
        None
    }
}

/// The parameter direction (`true` for u) and domain of a surface swept
/// from a closed NURBS curve (its two ends one point): an extrusion of one,
/// or a revolution of one.
fn closed_direction(s: &Surface) -> Option<(bool, f64, f64)> {
    let (along_u, curve) = match s {
        Surface::Extrusion {
            curve: Curve3::Nurbs(n),
            ..
        } => (true, n),
        Surface::Revolution {
            profile: Curve3::Nurbs(n),
            ..
        } => (false, n),
        _ => return None,
    };
    let (a, b) = curve.domain();
    // The test profile normalization uses to split closed splines (at the
    // default tolerance), so every surface built from a split curve is seen.
    let closed = (curve.point(a) - curve.point(b)).length() <= Tolerance::default().linear;
    closed.then_some((along_u, a, b))
}

/// On a surface swept from a closed curve, boundary points over the
/// curve's closure point project to either end of its domain. A face never
/// runs across the closure point (it is a vertex of the profile), so each
/// such point belongs at the end of the domain on the face's side: the end
/// nearer the mean parameter of the loop's other points. Without this the
/// parameter-space boundary can jump across the whole domain and cross
/// itself.
fn unwrap_closure(s: &Surface, samples: &mut [Sample]) {
    let Some((along_u, a, b)) = closed_direction(s) else {
        return;
    };
    let get = |uv: DVec2| if along_u { uv.x } else { uv.y };
    let eps = 1e-6 * (b - a);
    let at_end = |t: f64| (t - a).abs() <= eps || (b - t).abs() <= eps;
    let (mut sum, mut n) = (0.0, 0usize);
    for smp in samples.iter() {
        let t = get(smp.uv);
        if !at_end(t) {
            sum += t;
            n += 1;
        }
    }
    if n == 0 {
        return;
    }
    let end = if sum / n as f64 > 0.5 * (a + b) { b } else { a };
    for smp in samples.iter_mut() {
        if at_end(get(smp.uv)) {
            if along_u {
                smp.uv.x = end;
            } else {
                smp.uv.y = end;
            }
        }
    }
}

/// The boundary loops of `face` in its parameter space, each sample tagged
/// with the edge segment that starts at it. A singular point (a cone's
/// apex, a sphere's pole) becomes a run along the ambiguous parameter, on
/// a surface of revolution split into steps of at most `max_angle` (and
/// π/4) of its angle: the triangles meeting there
/// then span no more of it than their neighbours do (two points spanning
/// half a turn made fans the refinement took tens of passes to resolve).
/// Fails past [`VERTEX_BUDGET`] samples (a tiny `max_angle` asks for
/// hundreds of millions at a pole: found by review, as an allocation
/// that aborted the process).
fn loops_uv(face: &Face, edges: &[Vec<DVec3>], max_angle: f64) -> Result<Vec<Vec<Sample>>, String> {
    let s = &face.surface;
    let mut out = Vec::new();
    let mut count = 0usize;
    for lp in &face.loops {
        let mut pts: Vec<(DVec3, (usize, usize))> = Vec::new();
        for c in &lp.coedges {
            let e = c.edge.idx();
            let ep = &edges[e];
            let n = ep.len();
            if c.forward {
                pts.extend((0..n - 1).map(|i| (ep[i], (e, i))));
            } else {
                // Walking back from point i, the segment is [i − 1, i].
                pts.extend((1..n).rev().map(|i| (ep[i], (e, i - 1))));
            }
        }
        let mut raw: Vec<Sample> = pts
            .iter()
            .map(|&(p, seg)| {
                let (u, v) = s.project(p);
                Sample {
                    uv: DVec2::new(u, v),
                    p,
                    seg: Some(seg),
                }
            })
            .collect();
        unwrap_closure(s, &mut raw);
        // Split singular points into two parameter-space points.
        let n = raw.len();
        let mut fixed: Vec<Sample> = Vec::with_capacity(n + 4);
        for i in 0..n {
            let cur = raw[i];
            match singular_dir(s, cur.uv) {
                None => fixed.push(cur),
                Some(u_ambiguous) => {
                    let prev = raw[(i + n - 1) % n].uv;
                    let next = raw[(i + 1) % n].uv;
                    let (a, b) = if u_ambiguous {
                        (DVec2::new(prev.x, cur.uv.y), DVec2::new(next.x, cur.uv.y))
                    } else {
                        (DVec2::new(cur.uv.x, prev.y), DVec2::new(cur.uv.x, next.y))
                    };
                    if (a - b).length() > 1e-12 {
                        // On a surface of revolution u is the angle.
                        let steps = if u_ambiguous && matches!(s, Surface::Revolution { .. }) {
                            ((b.x - a.x).abs() / max_angle.min(FRAC_PI_4))
                                .ceil()
                                .max(1.0)
                        } else {
                            1.0
                        };
                        if !(steps + (count + fixed.len()) as f64 <= VERTEX_BUDGET as f64) {
                            return Err("tessellation exceeded the vertex budget".into());
                        }
                        for k in 0..steps as usize {
                            fixed.push(Sample {
                                uv: a + (b - a) * (k as f64 / steps),
                                p: cur.p,
                                seg: None,
                            });
                        }
                        fixed.push(Sample { uv: b, ..cur });
                    } else {
                        fixed.push(Sample { uv: a, ..cur });
                    }
                }
            }
        }
        count += fixed.len();
        out.push(fixed);
    }
    if count > VERTEX_BUDGET {
        return Err("tessellation exceeded the vertex budget".into());
    }
    Ok(out)
}

/// The scale of a face's parameters for its triangulation on a cylinder
/// or a cone of revolution, straight along `v`, given `param_scale`'s and
/// the face's extent `dv` in `v`: there a triangle's chord error depends
/// only on its width across `u` (an arc of radius at most `su`, the
/// scale of `u`), which the tolerance holds to `h = √(8 tol su)`, so `v`
/// is shrunk until a triangle `h` wide spans the face along `v` (at most
/// 256 times). The triangulation then follows the rulings, as on a face
/// bounded by its parameter lines, instead of filling a face cut across
/// them (by an ellipse, a conic or a line off the rulings) with
/// triangles `h` wide each way. The refinement measures the error in
/// space, so the mesh meets the tolerance either way.
fn ruled_scale(s: &Surface, scale: DVec2, dv: f64, tol: f64) -> DVec2 {
    let ruled = matches!(
        s,
        Surface::Revolution {
            profile: Curve3::Line { .. },
            ..
        } | Surface::Extrusion {
            curve: Curve3::Circle { .. },
            ..
        }
    );
    if !ruled || dv.is_nan() || dv <= 0.0 {
        return scale;
    }
    let h = (8.0 * tol * scale.x).sqrt();
    let shrink = (h / (scale.y * dv)).clamp(1.0 / 256.0, 1.0);
    DVec2::new(scale.x, scale.y * shrink)
}

/// Parameter-space scaling for a roughly isometric triangulation: the
/// largest `|S_u|` and `|S_v|` over a grid on the loops' bounding box.
fn param_scale(s: &Surface, loops: &[Vec<Sample>]) -> Result<DVec2, String> {
    let (mut lo, mut hi) = (DVec2::splat(f64::MAX), DVec2::splat(f64::MIN));
    for l in loops {
        for smp in l {
            lo = lo.min(smp.uv);
            hi = hi.max(smp.uv);
        }
    }
    let mut su = 0.0f64;
    let mut sv = 0.0f64;
    for fu in [0.25, 0.5, 0.75] {
        for fv in [0.25, 0.5, 0.75] {
            let uv = lo + (hi - lo) * DVec2::new(fu, fv);
            let (a, b) = s.derivs(uv.x, uv.y);
            su = su.max(a.length());
            sv = sv.max(b.length());
        }
    }
    if !(su > 0.0 && sv > 0.0 && su.is_finite() && sv.is_finite()) {
        return Err("degenerate surface parameterization".into());
    }
    Ok(DVec2::new(su, sv))
}

/// The constrained Delaunay triangulation the kernel uses.
pub(crate) type Cdt = ConstrainedDelaunayTriangulation<Point2<f64>>;

/// The constrained Delaunay triangulation of a polygon: constraint edges
/// around each of `loops` (closed cycles of points), and, per point in loop
/// order (the first loop's points, then the second's, …), its vertex index
/// in the triangulation. Exactly equal points are one vertex, numbered in
/// order of first appearance; an error if two boundary edges cross.
///
/// The triangulation is bulk loaded, which costs O(n log n) Delaunay tests
/// whatever the input. Inserting a polygon's vertices one at a time around
/// its boundary puts each one outside the hull of the earlier ones, and for
/// vertices on a circle (any regular polygon) that costs a number of tests
/// growing with the square of their count.
///
/// The result is the same on every run: the bulk load of the published
/// triangulation library inserts the vertices it defers in the iteration
/// order of a hash set seeded at random, which changes the triangulation
/// of points on a common circle, and the workspace uses a copy patched to
/// defer them in order (vendor/README.md).
pub(crate) fn polygon_cdt(loops: &[Vec<DVec2>]) -> Result<(Cdt, Vec<usize>), String> {
    let n: usize = loops.iter().map(Vec::len).sum();
    let mut first: FxHashMap<[u64; 2], usize> =
        FxHashMap::with_capacity_and_hasher(n, Default::default());
    let mut verts: Vec<Point2<f64>> = Vec::with_capacity(n);
    let mut index = Vec::with_capacity(n);
    // Coordinates this far below the polygon's size are rounding of a
    // zero (a cap's centre mapped through a frame whose axes carry
    // rounding came out at 8e-49), and the triangulation refuses nonzero
    // coordinates under about 1.8e-43: they are taken as zero.
    let size = loops
        .iter()
        .flatten()
        .fold(0.0f64, |m, q| m.max(q.x.abs()).max(q.y.abs()));
    let tiny = (size * 1e-30).max(1e-42);
    let flush = |x: f64| if x.abs() < tiny { 0.0 } else { x };
    for q in loops.iter().flatten() {
        let q = DVec2::new(flush(q.x), flush(q.y));
        // `+ 0.0` makes a negative zero positive, so equal points share a key.
        let key = [(q.x + 0.0).to_bits(), (q.y + 0.0).to_bits()];
        let i = *first.entry(key).or_insert_with(|| {
            verts.push(Point2::new(q.x, q.y));
            verts.len() - 1
        });
        index.push(i);
    }
    let mut seen: FxHashSet<[usize; 2]> = FxHashSet::default();
    let mut edges: Vec<[usize; 2]> = Vec::new();
    let mut start = 0;
    for l in loops {
        for k in 0..l.len() {
            let (a, b) = (index[start + k], index[start + (k + 1) % l.len()]);
            if a != b && seen.insert([a.min(b), a.max(b)]) {
                edges.push([a, b]);
            }
        }
        start += l.len();
    }
    let mut crossing = false;
    let cdt = Cdt::try_bulk_load_cdt(verts, edges, |_| crossing = true)
        .map_err(|e| format!("triangulation failed: {e:?}"))?;
    if crossing {
        return Err("a boundary intersects itself".into());
    }
    Ok((cdt, index))
}

/// Whether a closed polygon is convex and counter-clockwise: `Some(true)`
/// if every turn is strictly to the left, `Some(false)` if some vertices
/// are straight on (on one line with their neighbours, between them), and
/// `None` for any other polygon.
///
/// Decided exactly: every turn is to the left or straight on (never back),
/// at least one is to the left, and the edges' directions sweep one full
/// turn, which they do if their x components change sign at most twice;
/// more would wind the boundary around more than once (Schorn and Fisher,
/// "Testing the convexity of a polygon", *Graphics Gems IV*, 1994).
pub(crate) fn convexity(l: &[DVec2]) -> Option<bool> {
    let n = l.len();
    if n < 3 {
        return None;
    }
    let sign = |a: f64, b: f64| (b > a) as i8 - (b < a) as i8;
    let dir = |k: usize| {
        let (a, b) = (l[k], l[(k + 1) % n]);
        (sign(a.x, b.x), sign(a.y, b.y))
    };
    // Changes of sign of the edges' x components around the loop: from
    // the last edge's that is not zero.
    let mut x_was = (0..n).rev().map(|k| dir(k).0).find(|&x| x != 0)?;
    let (mut left, mut strict, mut flips) = (false, true, 0);
    for k in 0..n {
        let (d0, d1) = (dir(k), dir((k + 1) % n));
        if d0 == (0, 0) {
            return None;
        }
        let o = orient(l[k], l[(k + 1) % n], l[(k + 2) % n]);
        if o < 0.0 {
            return None;
        }
        if o > 0.0 {
            left = true;
        } else if (d0.0 != 0 && d0.0 == -d1.0) || (d0.1 != 0 && d0.1 == -d1.1) {
            // Straight back along the same line.
            return None;
        } else {
            strict = false;
        }
        if d0.0 != 0 && d0.0 != x_was {
            flips += 1;
            x_was = d0.0;
        }
    }
    (left && flips <= 2).then_some(strict)
}

/// The triangles, by position in the loop, that split a convex polygon of
/// `n` vertices by halving the boundary's runs level by level: the first
/// spans the polygon, and each further one a run of the boundary and the
/// chord across it. Each triangle's box is about the size of its run (the
/// slivers of a fan all reach across the polygon).
pub(crate) fn balanced(n: usize) -> Vec<[usize; 3]> {
    let mut out = vec![[0, n / 3, 2 * n / 3]];
    let mut runs = vec![(0, n / 3), (n / 3, 2 * n / 3), (2 * n / 3, n)];
    while let Some((i, j)) = runs.pop() {
        if j - i >= 2 {
            let m = (i + j) / 2;
            out.push([i, m, j % n]);
            runs.push((m, j));
            runs.push((i, m));
        }
    }
    out
}

/// The points of the quartic Lagrange triangle other than its corners, in
/// quarters of their barycentric coordinates: a quarter, half and three
/// quarters of the way along each edge `01`, `12` and `20`, then the three
/// inside, halfway from the centroid to each corner.
const QUARTERS: [[usize; 3]; 12] = [
    [3, 1, 0],
    [2, 2, 0],
    [1, 3, 0],
    [0, 3, 1],
    [0, 2, 2],
    [0, 1, 3],
    [1, 0, 3],
    [2, 0, 2],
    [3, 0, 1],
    [2, 1, 1],
    [1, 2, 1],
    [1, 1, 2],
];

/// No quartic that is zero at a triangle's corners is larger anywhere in it
/// than this many times its largest value at the [`QUARTERS`] (the sum of
/// the magnitudes of their Lagrange functions peaks at 3.3536, near the
/// corners; checked by a test).
const LEBESGUE: f64 = 3.36;

/// The quartic Lagrange functions of the [`QUARTERS`] at barycentric
/// coordinates `l`: products of `ℓₘ(x) = ∏_{a<m} (4x − a)/(a + 1)`, one
/// for each coordinate.
fn quartic_basis(l: [f64; 3]) -> [f64; 12] {
    let f = l.map(|x| {
        let y = 4.0 * x;
        let mut v = [1.0; 5];
        for m in 1..5 {
            v[m] = v[m - 1] * (y - (m - 1) as f64) / m as f64;
        }
        v
    });
    QUARTERS.map(|[i, j, k]| f[0][i] * f[1][j] * f[2][k])
}

/// The steps along each side of [`GRID`], and its points.
const STEPS: usize = 12;
const POINTS: usize = (STEPS + 1) * (STEPS + 2) / 2;

/// A point of [`GRID`]: its barycentric coordinates, the quartic Lagrange
/// functions there, and its neighbours' indices.
type GridPoint = ([f64; 3], [f64; 12], Vec<usize>);

/// A grid of twelfths over a triangle (twelfths hold the thirds, where a
/// cubic error along an edge peaks, and the quarters).
static GRID: std::sync::LazyLock<Vec<GridPoint>> = std::sync::LazyLock::new(|| {
    let n = STEPS as i32;
    let index = |i: i32, j: i32| -> Option<usize> {
        (i >= 0 && j >= 0 && i + j <= n).then(|| {
            // Rows of `i` hold `n + 1 - i` points.
            (0..i).map(|r| (n + 1 - r) as usize).sum::<usize>() + j as usize
        })
    };
    let mut out = Vec::new();
    for i in 0..=n {
        for j in 0..=n - i {
            let l = [i, j, n - i - j].map(|x| x as f64 / n as f64);
            let near = [(1, 0), (-1, 0), (0, 1), (0, -1), (1, -1), (-1, 1)]
                .into_iter()
                .filter_map(|(a, b)| index(i + a, j + b))
                .collect();
            out.push((l, quartic_basis(l), near));
        }
    }
    out
});

/// The value of largest magnitude, over a triangle, of the quartic that is
/// zero at its corners and `at` the [`QUARTERS`] (the quartic Lagrange
/// interpolant of the distance to the surface, the estimate between the
/// points tested; after R. E. Bank and R. K. Smith, "Mesh smoothing using a
/// posteriori error estimates", SIAM J. Numer. Anal. 34, 1997, who estimate
/// a linear element's error from the quadratic one), and where it is, in
/// barycentric coordinates. The edges `01`, `12` and `20` that `skip`
/// marks are left out. Found on a grid of twelfths, then by a local
/// search, down to 1/384, from each of the grid's local peaks within a
/// tenth of its largest (design note 0002 §4). Given `stop`, the search
/// ends at the first point found past it, and is not refined when the
/// grid's largest is under it by more than a tenth.
fn peak(at: &[f64; 12], skip: [bool; 3], stop: Option<f64>) -> (f64, [f64; 3]) {
    let refine = |top: f64| stop.is_none_or(|s| top * 1.1 >= s);
    let stop = stop.unwrap_or(f64::INFINITY);
    let value = |b: &[f64; 12]| b.iter().zip(at).map(|(x, y)| x * y).sum::<f64>();
    let left_out = |l: &[f64; 3]| (0..3).any(|e| skip[e] && l[(e + 2) % 3] < 1e-9);
    let grid = &*GRID;
    let mut v = [0.0f64; POINTS];
    let mut top = 0.0f64;
    for (k, (l, b, _)) in grid.iter().enumerate() {
        if !left_out(l) {
            v[k] = value(b);
            if v[k].abs() > stop {
                return (v[k], *l);
            }
            top = top.max(v[k].abs());
        }
    }
    let mut best = (0.0f64, [1.0 / 3.0; 3]);
    if top == 0.0 {
        return best;
    }
    for (k, (l, _, near)) in grid.iter().enumerate() {
        if v[k].abs() * 1.1 < top || near.iter().any(|&m| v[m].abs() > v[k].abs()) {
            continue;
        }
        let mut here = (v[k], *l);
        if refine(top) {
            let mut step = 1.0 / 24.0;
            for _ in 0..5 {
                let (_, l0) = here;
                for (i, j) in [(0, 1), (1, 0), (1, 2), (2, 1), (2, 0), (0, 2)] {
                    let mut l = l0;
                    l[i] += step;
                    l[j] -= step;
                    if l[j] > -1e-12 {
                        l[j] = l[j].max(0.0);
                        if !left_out(&l) {
                            let x = value(&quartic_basis(l));
                            if x.abs() > here.0.abs() {
                                here = (x, l);
                                if x.abs() > stop {
                                    return here;
                                }
                            }
                        }
                    }
                }
                step *= 0.5;
            }
        }
        if here.0.abs() > best.0.abs() {
            best = here;
        }
    }
    best
}

/// Triangles inside the constrained boundary, found by flood fill from the
/// outer face, toggling inside/outside each time a constraint edge is
/// crossed. Purely topological, so robust to nearly collinear boundaries.
pub(crate) fn inside_triangles(
    cdt: &ConstrainedDelaunayTriangulation<Point2<f64>>,
) -> Vec<[usize; 3]> {
    let mut parity: Vec<Option<bool>> = vec![None; cdt.num_all_faces()];
    let mut queue = std::collections::VecDeque::new();
    for f in cdt.inner_faces() {
        for e in f.adjacent_edges() {
            if e.rev().face().is_outer() {
                let inside = e.as_undirected().data().is_constraint_edge();
                let i = f.fix().index();
                if parity[i].is_none() {
                    parity[i] = Some(inside);
                    queue.push_back(f.fix());
                }
            }
        }
    }
    while let Some(fx) = queue.pop_front() {
        let f = cdt.face(fx);
        let here = parity[fx.index()].unwrap_or(false);
        for e in f.adjacent_edges() {
            if let Some(g) = e.rev().face().as_inner() {
                let gi = g.fix().index();
                if parity[gi].is_none() {
                    parity[gi] = Some(here ^ e.as_undirected().data().is_constraint_edge());
                    queue.push_back(g.fix());
                }
            }
        }
    }
    cdt.inner_faces()
        .filter(|f| parity[f.fix().index()] == Some(true))
        .map(|f| {
            let v = f.vertices();
            [v[0].fix().index(), v[1].fix().index(), v[2].fix().index()]
        })
        .collect()
}

fn tessellate_face(
    face: &Face,
    loops: Vec<Vec<Sample>>,
    opts: &TessellationOptions,
    linear: f64,
) -> Result<FaceTess, Unmeshed> {
    let s = &face.surface;
    if loops.iter().any(|l| l.len() < 3) {
        return Err(Unmeshed::Failed(
            "loop has fewer than three boundary samples".into(),
        ));
    }
    let (mut lo, mut hi) = (DVec2::splat(f64::MAX), DVec2::splat(f64::MIN));
    for l in &loops {
        for smp in l {
            lo = lo.min(smp.uv);
            hi = hi.max(smp.uv);
        }
    }
    let scale = ruled_scale(
        s,
        param_scale(s, &loops)?,
        hi.y - lo.y,
        opts.chord_tolerance,
    );

    let is_planar = matches!(s, Surface::Plane(_));
    // Spade rejects non-zero coordinates below ~1e-42; such values are
    // round-off (e.g. an angle of 1e-60 rad) and are flushed to zero.
    let flush = |x: f64| if x.abs() < 1e-30 { 0.0 } else { x };
    let place = |uv: DVec2| {
        let q = uv * scale;
        DVec2::new(flush(q.x), flush(q.y))
    };
    // Per CDT vertex: parameter value and, for boundary vertices, the exact 3D point.
    let mut verts: Vec<(DVec2, Option<DVec3>)> = Vec::new();
    // Per CDT vertex of a curved face: its surface point and normal,
    // evaluated once.
    let mut at: Vec<(DVec3, Option<DVec3>)> = Vec::new();
    // At a singular point (an apex, a pole, a disc's centre) `S_u × S_v`
    // is the rounding of a vanishing derivative and can point either way:
    // such a point has no normal of its own.
    let evaluate = |uv: DVec2, p: Option<DVec3>| {
        let sp = p.unwrap_or_else(|| s.point(uv.x, uv.y));
        let n = singular_dir(s, uv).is_none().then(|| s.normal(uv.x, uv.y));
        (sp, n.flatten())
    };
    let placed: Vec<Vec<DVec2>> = loops
        .iter()
        .map(|l| l.iter().map(|smp| place(smp.uv)).collect())
        .collect();
    // A flat face bounded by one strictly convex loop is split by halving
    // the boundary's runs: no point is added on a flat face, every
    // triangulation of a convex polygon on its own vertices is valid (no
    // three of them are on a line), and this one needs no search.
    let mut unfolded = 0;
    let tris = if is_planar && placed.len() == 1 && convexity(&placed[0]) == Some(true) {
        verts = loops[0].iter().map(|smp| (smp.uv, Some(smp.p))).collect();
        balanced(placed[0].len())
    } else {
        let (mut cdt, index) =
            polygon_cdt(&placed).map_err(|e| format!("face boundary in parameter space: {e}"))?;
        // Vertices are numbered by the triangulation; a repeated boundary point
        // keeps its first sample.
        let mut filled = vec![false; cdt.num_vertices()];
        verts.resize(cdt.num_vertices(), (DVec2::ZERO, None));
        if !is_planar {
            at.resize(cdt.num_vertices(), (DVec3::ZERO, None));
        }
        for (&i, smp) in index.iter().zip(loops.iter().flatten()) {
            if !filled[i] {
                filled[i] = true;
                verts[i] = (smp.uv, Some(smp.p));
                if !is_planar {
                    at[i] = evaluate(smp.uv, Some(smp.p));
                }
            }
        }
        let insert = |cdt: &mut Cdt,
                      verts: &mut Vec<(DVec2, Option<DVec3>)>,
                      at: &mut Vec<(DVec3, Option<DVec3>)>,
                      uv: DVec2,
                      hint: usize|
         -> Result<(), String> {
            let q = place(uv);
            let h = cdt
                .insert_with_hint(
                    Point2::new(q.x, q.y),
                    spade::handles::FixedVertexHandle::from_index(hint),
                )
                .map_err(|e| format!("triangulation insert failed: {e:?}"))?;
            if h.index() == verts.len() {
                verts.push((uv, None));
                if !is_planar {
                    at.push(evaluate(uv, None));
                }
            }
            Ok(())
        };

        let boundary_constraints = cdt.num_constraints();
        if is_planar {
            inside_triangles(&cdt)
        } else {
            let min_area = (opts.chord_tolerance * 0.01).powi(2);
            // The boundary's segments, each with the edge segment it is
            // (`None` between the points a singular point is split into).
            // The edges' samples hold their chords to the tolerance, so
            // their midpoints are not tested (no point added inside the
            // face would change them). Refinement adds no constraint, so
            // they are found once.
            let mut boundary: FxHashMap<[usize; 2], Option<(usize, usize)>> = FxHashMap::default();
            let mut start = 0;
            for l in &loops {
                for (j, smp) in l.iter().enumerate() {
                    let (a, b) = (index[start + j], index[start + (j + 1) % l.len()]);
                    if a != b {
                        boundary.entry([a.min(b), a.max(b)]).or_insert(smp.seg);
                    }
                }
                start += l.len();
            }
            // How triangle `t` fares (design note 0002 §4): fine, to be
            // split at a point, or beside a boundary segment to bisect
            // first.
            let judge = |t: [usize; 3],
                         verts: &[(DVec2, Option<DVec3>)],
                         at: &[(DVec3, Option<DVec3>)],
                         boundary: &FxHashMap<[usize; 2], Option<(usize, usize)>>,
                         min_area: f64|
             -> Result<Judged, Unmeshed> {
                let uv = t.map(|i| verts[i].0);
                let q = uv.map(|x| x * scale);
                let area = 0.5 * (q[1] - q[0]).perp_dot(q[2] - q[0]).abs();
                let lens = [0, 1, 2].map(|e| (q[(e + 1) % 3] - q[e]).length_squared());
                let longest = lens[0].max(lens[1]).max(lens[2]);
                let side = |e: usize| {
                    let (a, b) = (t[e], t[(e + 1) % 3]);
                    boundary.get(&[a.min(b), a.max(b)])
                };
                let p = t.map(|i| at[i].0);
                let tri_n = (p[1] - p[0]).cross(p[2] - p[0]);
                if tri_n.length() <= 0.0 {
                    return Ok(Judged::Fine);
                }
                // Its corners in a line within the model tolerance, a
                // triangle has no side to face (by an apex, three points
                // on nearly one ruling).
                let side_3d = [0, 1, 2].map(|e| (p[(e + 1) % 3] - p[e]).length());
                let thin = tri_n.length() < linear * side_3d[0].max(side_3d[1]).max(side_3d[2]);
                let tri_n = tri_n.normalize();
                let c = (uv[0] + uv[1] + uv[2]) / 3.0;
                // A triangle facing away from the surface (a long sliver
                // folded over), or a normal turning too far.
                let nc = s.normal(c.x, c.y);
                let folded = !thin && nc.is_some_and(|nc| tri_n.dot(nc) <= 0.0);
                let turned = folded
                    || nc.is_some_and(|nc| {
                        t.iter().any(|&i| {
                            at[i]
                                .1
                                .is_some_and(|nw| nw.angle_between(nc) > opts.max_angle)
                        })
                    });
                if turned {
                    // Beside a boundary segment that is its longest side,
                    // points added inside close in on the segment's middle
                    // in the parameters, which can lie off its chord: the
                    // segment is bisected (a segment Delaunay refinement
                    // splits rather than crowd: J. Ruppert, "A Delaunay
                    // refinement algorithm for quality 2-dimensional mesh
                    // generation", J. Algorithms 18, 1995; found by review,
                    // as triangles folded over by a cut near a ball's pole
                    // and by an apex, settled once too small to split).
                    let e = (0..3)
                        .max_by(|&a, &b| lens[a].total_cmp(&lens[b]))
                        .unwrap_or(0);
                    if let Some(Some(seg)) = side(e) {
                        return Ok(Judged::Bisect(*seg));
                    }
                }
                // Tiny triangles and slivers are not split (their centroid
                // could land on a boundary constraint within round-off).
                // One folded over is unfolded once refinement is done; one
                // whose normals turn too far is refused.
                if area < min_area || 2.0 * area < 1e-8 * longest {
                    if folded {
                        return Ok(Judged::Folded);
                    }
                    if turned {
                        return Err(
                            "a triangle's normals turn past the angle tolerance where it is \
                             too small to split"
                                .to_string()
                                .into(),
                        );
                    }
                    return Ok(Judged::Fine);
                }
                // Chord error: the distance to the surface from the
                // triangle's centroid, the points halfway from it to the
                // corners, and the points at a quarter, half and three
                // quarters of its edges, and the largest value inside it of
                // the quartic through those points (zero at the corners).
                // The worst failing edge is split at its midpoint if it is
                // the triangle's longest (Rivara's bisection); otherwise
                // the centroid is added. Estimated, not bounded: where the
                // surface bends faster than a quartic between the points,
                // the mesh can stray a little further (design note 0002
                // §4).
                let tol = opts.chord_tolerance;
                // The distance from `q`, a point of the triangle, to the
                // surface, from the surface point at `w` under it: along
                // the surface's normal there (an offset along the surface
                // is the parameters' slack, not the mesh's error), signed
                // by the side of the triangle it is on.
                let gap = |w: DVec2, q: DVec3| {
                    let off = s.point(w.x, w.y) - q;
                    match s.normal(w.x, w.y) {
                        Some(n) => off.dot(n) * n.dot(tri_n).signum(),
                        None => off.dot(tri_n),
                    }
                };
                // The parameter of a corner that is free: at an apex or a
                // pole, always a boundary sample.
                let free = |i: usize| verts[i].1.and_then(|_| singular_dir(s, verts[i].0));
                // The parameters under the point `f` of the way along the
                // chord from `a` to `b`: those as far between theirs,
                // except that from an apex or a pole the chord runs along
                // the other end's ruling or meridian, so the free
                // parameter is that end's (the parameters as far between
                // are a point beside the chord, where the distance along
                // the normal overstates the chord's: review measured a
                // polar cap at 10,424 triangles that way, against 7,262).
                let under = |a: usize, b: usize, f: f64| {
                    let (x, y) = (verts[a].0, verts[b].0);
                    let m = x.lerp(y, f);
                    match (free(a), free(b)) {
                        (Some(true), None) => DVec2::new(y.x, m.y),
                        (Some(false), None) => DVec2::new(m.x, y.y),
                        (None, Some(true)) => DVec2::new(x.x, m.y),
                        (None, Some(false)) => DVec2::new(m.x, x.y),
                        _ => m,
                    }
                };
                // The centroid and the points halfway from it to the
                // corners.
                let inner = [
                    [1.0 / 3.0; 3],
                    [0.5, 0.25, 0.25],
                    [0.25, 0.5, 0.25],
                    [0.25, 0.25, 0.5],
                ]
                .map(|w| {
                    let x = uv[0] * w[0] + uv[1] * w[1] + uv[2] * w[2];
                    gap(x, p[0] * w[0] + p[1] * w[1] + p[2] * w[2])
                });
                let mut split = inner.iter().any(|x| x.abs() > tol).then_some(c);
                // The distances at the quartic triangle's points
                // ([`QUARTERS`]).
                let mut d = [0.0; 12];
                d[9..].copy_from_slice(&inner[1..]);
                let mut worst = tol;
                for e in 0..3 {
                    let (a, b) = (t[e], t[(e + 1) % 3]);
                    let along =
                        [0.25, 0.5, 0.75].map(|f| gap(under(a, b, f), at[a].0.lerp(at[b].0, f)));
                    if side(e).is_some() {
                        // The edges' samples hold a boundary chord to the
                        // tolerance; the parameters between its ends only
                        // approximate the curve's, so the distances are
                        // capped there and not tested.
                        for (x, y) in d[3 * e..3 * e + 3].iter_mut().zip(along) {
                            *x = y.clamp(-tol, tol);
                        }
                        continue;
                    }
                    d[3 * e..3 * e + 3].copy_from_slice(&along);
                    let far = along.into_iter().fold(0.0f64, |m, x| m.max(x.abs()));
                    if far > worst {
                        worst = far;
                        let m = (verts[a].0 + verts[b].0) * 0.5;
                        split = Some(if lens[e] >= longest { m } else { c });
                    }
                }
                // Between the points: the quartic through them, searched
                // only when it could pass the tolerance, and not on the
                // boundary's segments (no point added inside moves them).
                // Past it on an edge that is the triangle's longest, that
                // edge is bisected; elsewhere the centroid is added.
                let large = d.iter().fold(0.0f64, |m, x| m.max(x.abs()));
                if std::env::var_os("NOPEAK").is_none() && split.is_none() && large * LEBESGUE > tol
                {
                    let on = |e: usize, l: [f64; 3]| l[(e + 2) % 3] < 1e-9;
                    let (v, l) = peak(&d, [0, 1, 2].map(|e| side(e).is_some()), Some(tol));
                    if v.abs() > tol {
                        split = Some(match (0..3).find(|&e| on(e, l)) {
                            Some(e) if lens[e] >= longest => {
                                (verts[t[e]].0 + verts[t[(e + 1) % 3]].0) * 0.5
                            }
                            _ => c,
                        });
                    }
                }
                if split.is_none() && turned {
                    split = Some(c);
                }
                Ok(split.map_or(Judged::Fine, Judged::Insert))
            };
            // Whether triangle `t` (counter-clockwise) is still one of the
            // triangulation's.
            let present = |cdt: &Cdt, t: [usize; 3]| {
                let h = spade::handles::FixedVertexHandle::from_index;
                cdt.get_edge_from_neighbors(h(t[0]), h(t[1]))
                    .and_then(|e| e.face().as_inner())
                    .is_some_and(|f| f.vertices().iter().any(|v| v.fix().index() == t[2]))
            };
            // Each pass tests the triangles made by the last (around the
            // points it added: inserting a point re-triangulates only the
            // triangles about it) and splits those that fail, until none
            // does: a triangle's test depends only on its corners. A pass
            // that adds no point while triangles fail, or a mesh past the
            // vertex budget, is refused.
            let mut work = inside_triangles(&cdt);
            let mut splits: Vec<(usize, usize)> = Vec::new();
            let mut folded: Vec<[usize; 3]> = Vec::new();
            while !work.is_empty() {
                let mut new_pts: Vec<(DVec2, usize)> = Vec::new();
                let mut failing: Vec<[usize; 3]> = Vec::new();
                for &t in &work {
                    if !present(&cdt, t) {
                        continue;
                    }
                    match judge(t, &verts, &at, &boundary, min_area)? {
                        Judged::Fine => {}
                        Judged::Insert(w) => {
                            new_pts.push((w, t[0]));
                            failing.push(t);
                        }
                        Judged::Bisect(seg) => splits.push(seg),
                        Judged::Folded => folded.push(t),
                    }
                }
                if !splits.is_empty() {
                    splits.sort_unstable();
                    splits.dedup();
                    return Err(Unmeshed::Split(splits));
                }
                if new_pts.is_empty() {
                    break;
                }
                if verts.len() + new_pts.len() > VERTEX_BUDGET {
                    return Err("tessellation exceeded the vertex budget".to_string().into());
                }
                let before = verts.len();
                for (c, hint) in new_pts {
                    insert(&mut cdt, &mut verts, &mut at, c, hint)?;
                }
                if verts.len() == before {
                    return Err("refining the mesh added no point where triangles fail"
                        .to_string()
                        .into());
                }
                // The triangles about the new points, and any that failed
                // and are still there; each once, from its least corner (in
                // order, so the result does not depend on a hash table's).
                let mut next: Vec<[usize; 3]> = Vec::new();
                for v in before..verts.len() {
                    let h = spade::handles::FixedVertexHandle::from_index(v);
                    for e in cdt.vertex(h).out_edges() {
                        if let Some(f) = e.face().as_inner() {
                            next.push(f.vertices().map(|x| x.fix().index()));
                        }
                    }
                }
                next.extend(failing.into_iter().filter(|&t| present(&cdt, t)));
                for t in &mut next {
                    let m = (0..3).min_by_key(|&i| t[i]).unwrap_or(0);
                    t.rotate_left(m);
                }
                next.sort_unstable();
                next.dedup();
                work = next;
            }
            if cdt.num_constraints() != boundary_constraints {
                return Err("refinement split a boundary edge".to_string().into());
            }
            // Unfolding: a sliver too small to split can face against the
            // surface (three points nearly in a line in the parameters,
            // the surface's bend carrying the middle one across the chord
            // of the others in space; found by the near-tangent apex fuzz,
            // design note 0002 §4). The triangle and the one across its
            // longest edge that allows it are re-split along their other
            // diagonal (an edge flip, as in C. L. Lawson, "Software for C¹
            // surface interpolation", 1977), kept only if both new
            // triangles turn counter-clockwise in the parameters and pass
            // the test above; one that no flip mends is refused.
            let mut tris = inside_triangles(&cdt);
            if !folded.is_empty() {
                let canon = |t: [usize; 3]| {
                    let m = (0..3).min_by_key(|&i| t[i]).unwrap_or(0);
                    let mut t = t;
                    t.rotate_left(m);
                    t
                };
                let mut slot: FxHashMap<[usize; 3], usize> = FxHashMap::default();
                let mut across: FxHashMap<(usize, usize), usize> = FxHashMap::default();
                for (i, &t) in tris.iter().enumerate() {
                    slot.insert(canon(t), i);
                    for e in 0..3 {
                        across.insert((t[e], t[(e + 1) % 3]), i);
                    }
                }
                folded.sort_unstable();
                folded.dedup();
                let ccw = |t: [usize; 3]| {
                    let q = t.map(|i| verts[i].0 * scale);
                    (q[1] - q[0]).perp_dot(q[2] - q[0]) > 0.0
                };
                for f in folded {
                    let Some(&i) = slot.get(&canon(f)) else {
                        continue;
                    };
                    let t = tris[i];
                    if !matches!(judge(t, &verts, &at, &boundary, min_area)?, Judged::Folded) {
                        continue;
                    }
                    let length = |e: usize| {
                        (verts[t[(e + 1) % 3]].0 * scale - verts[t[e]].0 * scale).length_squared()
                    };
                    let mut order = [0, 1, 2];
                    order.sort_by(|&a, &b| length(b).total_cmp(&length(a)));
                    let mut mended = false;
                    for e in order {
                        let (x, y, z) = (t[e], t[(e + 1) % 3], t[(e + 2) % 3]);
                        if boundary.contains_key(&[x.min(y), x.max(y)]) {
                            continue;
                        }
                        let Some(&j) = across.get(&(y, x)) else {
                            continue;
                        };
                        let u = tris[j];
                        let Some(&d) = u.iter().find(|&&w| w != x && w != y) else {
                            continue;
                        };
                        let (a, b) = ([x, d, z], [d, y, z]);
                        let fine = |n: [usize; 3]| -> Result<bool, Unmeshed> {
                            Ok(ccw(n)
                                && matches!(
                                    judge(n, &verts, &at, &boundary, min_area)?,
                                    Judged::Fine
                                ))
                        };
                        if fine(a)? && fine(b)? {
                            for (k, old) in [(i, t), (j, u)] {
                                slot.remove(&canon(old));
                                for e in 0..3 {
                                    across.remove(&(old[e], old[(e + 1) % 3]));
                                }
                                let new = if k == i { a } else { b };
                                tris[k] = new;
                                slot.insert(canon(new), k);
                            }
                            for k in [i, j] {
                                let n = tris[k];
                                for e in 0..3 {
                                    across.insert((n[e], n[(e + 1) % 3]), k);
                                }
                            }
                            mended = true;
                            unfolded += 1;
                            break;
                        }
                    }
                    if !mended {
                        return Err(
                            "a triangle faces against the surface where it is too small \
                                    to split or flip"
                                .to_string()
                                .into(),
                        );
                    }
                }
            }
            tris
        }
    };

    // Emit vertices. Boundary vertices with identical 3D points (split
    // singular points) are welded.
    let mut out = FaceTess {
        unfolded,
        ..FaceTess::default()
    };
    let mut index_of: Vec<Option<u32>> = vec![None; verts.len()];
    let mut by_pos: FxHashMap<[u64; 3], u32> = FxHashMap::default();
    let center = (lo + hi) * 0.5;
    let sign = if face.forward { 1.0 } else { -1.0 };
    // A singular point has no normal of its own: each of its points in the
    // parameters takes the normal a step from it toward the middle of the
    // face, along its meridian or ruling, and the welded vertex their mean
    // (an apex's normal was one meridian's, facing away from triangles
    // across a half-turn face of a slender cone: found by review).
    let mut singular: Vec<u32> = Vec::new();
    for t in &tris {
        for &i in t {
            if index_of[i].is_some() {
                continue;
            }
            let (uv, exact) = verts[i];
            let p = exact.unwrap_or_else(|| {
                if is_planar {
                    s.point(uv.x, uv.y)
                } else {
                    at[i].0
                }
            });
            let mut n = if is_planar {
                s.normal(uv.x, uv.y)
            } else {
                at[i].1
            };
            let stepped = n.is_none();
            let mut k = 1e-6;
            while n.is_none() && k < 0.5 {
                let w = uv + (center - uv) * k;
                n = s.normal(w.x, w.y);
                k *= 4.0;
            }
            let n = n.unwrap_or(DVec3::Z) * sign;
            let key = [p.x.to_bits(), p.y.to_bits(), p.z.to_bits()];
            if exact.is_some() {
                if let Some(&j) = by_pos.get(&key) {
                    if stepped {
                        out.normals[j as usize] += n;
                    }
                    index_of[i] = Some(j);
                    continue;
                }
            }
            let j = out.positions.len() as u32;
            out.positions.push(p);
            out.normals.push(n);
            if stepped {
                singular.push(j);
            }
            if exact.is_some() {
                by_pos.insert(key, j);
            }
            index_of[i] = Some(j);
        }
        let a = index_of[t[0]].unwrap_or(0);
        let b = index_of[t[1]].unwrap_or(0);
        let c = index_of[t[2]].unwrap_or(0);
        if a == b || b == c || a == c {
            continue;
        }
        if face.forward {
            out.triangles.push([a, b, c]);
        } else {
            out.triangles.push([a, c, b]);
        }
    }
    for j in singular {
        let n = &mut out.normals[j as usize];
        *n = n.try_normalize().unwrap_or(DVec3::Z * sign);
    }
    if out.triangles.is_empty() {
        return Err(Unmeshed::Failed("face produced no triangles".into()));
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Body, BooleanOp, BrepKernel, ExtrudeSpec, Kernel, Profile, ProfileCurve};
    use crate::{Frame, RevolveSpec};
    use std::f64::consts::{FRAC_PI_2, PI, TAU};

    fn k() -> BrepKernel {
        BrepKernel::default()
    }

    /// A solid turned a whole turn about z from `(ρ, z)` points.
    fn turned(curves: Vec<ProfileCurve>) -> Body {
        turned_by(curves, TAU)
    }

    /// A solid turned by `angle` about z from `(ρ, z)` points.
    fn turned_by(curves: Vec<ProfileCurve>, angle: f64) -> Body {
        let plane = Frame {
            origin: DVec3::ZERO,
            x: DVec3::X,
            y: DVec3::Z,
            z: -DVec3::Y,
        };
        k().revolve(
            &Profile { plane, curves },
            &RevolveSpec {
                axis_origin: DVec2::ZERO,
                axis_dir: DVec2::Y,
                angle,
            },
        )
        .expect("a revolve")
    }

    fn polygon(pts: &[(f64, f64)]) -> Vec<ProfileCurve> {
        (0..pts.len())
            .map(|i| {
                let (a, b) = (pts[i], pts[(i + 1) % pts.len()]);
                ProfileCurve::Line {
                    start: DVec2::new(a.0, a.1),
                    end: DVec2::new(b.0, b.1),
                }
            })
            .collect()
    }

    /// A box `across` square and very deep on the side of the plane
    /// through `point` that `normal` points to.
    fn slab(point: DVec3, normal: DVec3) -> Body {
        let f = Frame::from_normal(point, normal, DVec3::Y);
        k().make_box(
            &Frame {
                origin: point - (f.x + f.y) * 500.0,
                ..f
            },
            DVec3::splat(1000.0),
        )
        .expect("a box")
    }

    /// The farthest a curved face's mesh strays from its surface, in chord
    /// tolerances, over a grid of fifteen points of each triangle (its
    /// corners, the quarter points of its sides and three inside), against
    /// the surface's closest points; and the triangles facing against their
    /// vertices' normals.
    fn strays(body: &Body, opts: &TessellationOptions) -> (f64, usize) {
        let mesh = tessellate(&body.brep, opts, &Tolerance::default()).expect("a mesh");
        let (mut worst, mut folded) = (0.0f64, 0);
        for (f, m) in body.brep.faces.iter().zip(&mesh.faces) {
            if matches!(f.surface, Surface::Plane(_)) {
                continue;
            }
            for t in &m.triangles {
                let p = t.map(|i| m.positions[i as usize]);
                for i in 0..=4 {
                    for j in 0..=4 - i {
                        let (a, b) = (0.25 * i as f64, 0.25 * j as f64);
                        let q = p[0] * (1.0 - a - b) + p[1] * a + p[2] * b;
                        let (u, v) = f.surface.project(q);
                        worst =
                            worst.max((f.surface.point(u, v) - q).length() / opts.chord_tolerance);
                    }
                }
                let n = (p[1] - p[0]).cross(p[2] - p[0]);
                folded += t
                    .iter()
                    .filter(|&&i| n.dot(m.normals[i as usize]) < 0.0)
                    .count()
                    .min(1);
            }
        }
        (worst, folded)
    }

    /// Meshes of cones, a ball and faces cut across their rulings stay
    /// within the chord tolerance and face out, sampled at fifteen points
    /// of every triangle (its corners, the quarter points of its sides and
    /// three inside) against the surface's closest points (found by
    /// review, design note 0002 §4: refinement stopped after twelve passes,
    /// leaving a pointed cone's mesh 16 tolerances off by its apex at the
    /// defaults and 159 at a chord of 1e-5; it never tested a diagonal
    /// between two boundary samples, which faces meshed along their rulings
    /// are made of, so cut faces strayed by up to a quarter more; long
    /// slivers on a frustum folded over; and a ball strayed 5.2%, 8.2% at
    /// a chord of 1e-5, mostly 10° to 30° from its poles, and a polar cap
    /// 45% at 1e-5).
    #[test]
    fn meshes_stay_within_the_chord_tolerance() {
        let tilt = 79.3f64.to_radians();
        let cut = |body: Body, point: DVec3, normal: DVec3| {
            k().boolean(
                &body,
                &slab(point, normal.normalize()),
                BooleanOp::Intersect,
            )
            .expect("a cut")
        };
        let ball = turned(vec![
            ProfileCurve::Arc {
                center: DVec2::ZERO,
                radius: 1.0,
                start_angle: -FRAC_PI_2,
                sweep_angle: PI,
            },
            ProfileCurve::Line {
                start: DVec2::new(0.0, 1.0),
                end: DVec2::new(0.0, -1.0),
            },
        ]);
        // The finest tolerance meshes a polar cap of the ball, not all of
        // it (over a million triangles).
        let cap = cut(ball.clone(), DVec3::new(0.0, 0.0, 0.9), DVec3::Z);
        let mut bodies = vec![
            turned(polygon(&[(0.0, 0.0), (0.5, 0.0), (0.0, 1.0)])),
            turned(polygon(&[(0.0, -1.0), (1.0, 0.0), (0.0, 1.0)])),
            cap,
            cut(
                turned(polygon(&[(0.0, 0.0), (0.5, 0.0), (0.025, 1.0), (0.0, 1.0)])),
                DVec3::new(0.065, 0.0, 0.37),
                DVec3::new(tilt.sin(), 0.31 * tilt.sin(), tilt.cos()),
            ),
            cut(
                turned(polygon(&[(0.0, 0.0), (2.0, 0.0), (0.1, 1.0), (0.0, 1.0)])),
                DVec3::new(0.26, 0.0, 0.37),
                DVec3::new(0.8984, 0.2785, 0.3395),
            ),
            cut(
                k().extrude(
                    &Profile::circle(Frame::WORLD, DVec2::ZERO, 2.0),
                    &ExtrudeSpec {
                        distance: 1.0,
                        symmetric: false,
                    },
                )
                .expect("a cylinder"),
                DVec3::new(0.26, 0.0, 0.37),
                DVec3::new(0.954, 0.296, 0.045),
            ),
        ];
        // Frustums cut by planes at any angle.
        let mut rng = crate::boolean::testing::Rng(0x9E37_79B9_7F4A_7C15);
        for _ in 0..40 {
            let rad = rng.range(0.3, 2.0);
            let (top, h) = (rad * rng.range(0.02, 0.8), rng.range(0.3, 1.5));
            let (a, b) = (rng.range(0.0, TAU), rng.range(0.05, 1.5));
            let n = DVec3::new(b.sin() * a.cos(), b.sin() * a.sin(), b.cos());
            let point = DVec3::new(
                rng.range(-0.3, 0.3) * rad,
                rng.range(-0.3, 0.3) * rad,
                rng.range(0.2, 0.8) * h,
            );
            let frustum = turned(polygon(&[(0.0, 0.0), (rad, 0.0), (top, h), (0.0, h)]));
            if let Ok(x) = k().boolean(&frustum, &slab(point, n), BooleanOp::Intersect) {
                bodies.push(x);
            }
        }
        assert!(bodies.len() > 40, "{}", bodies.len());
        for opts in [
            TessellationOptions::default(),
            TessellationOptions {
                chord_tolerance: 1e-3,
                max_angle: 0.5,
            },
            TessellationOptions {
                chord_tolerance: 1e-5,
                max_angle: 0.5,
            },
        ] {
            let fine = opts.chord_tolerance < 1e-4;
            for (i, b) in bodies.iter().chain((!fine).then_some(&ball)).enumerate() {
                let (worst, folded) = strays(b, &opts);
                assert!(
                    worst <= 1.0 + 1e-6,
                    "body {i} {opts:?}: {worst} tolerances off"
                );
                assert_eq!(folded, 0, "body {i} {opts:?}");
            }
        }
    }

    /// A ball of radius `r` about the origin.
    fn ball(r: f64) -> Body {
        turned(vec![
            ProfileCurve::Arc {
                center: DVec2::ZERO,
                radius: r,
                start_angle: -FRAC_PI_2,
                sweep_angle: PI,
            },
            ProfileCurve::Line {
                start: DVec2::new(0.0, r),
                end: DVec2::new(0.0, -r),
            },
        ])
    }

    /// Faces whose boundary runs where the parameters' straight segments
    /// lie far off the chords (a cut passing a ball's pole, a cone cut
    /// steeply near its apex) mesh without a triangle facing against the
    /// surface: the segment is bisected (found by review: refinement drew a
    /// triangle's free corner to the segment's middle in the parameters,
    /// beside the chord, and settled it folded over once too small to
    /// split; 2 triangles of a polar cap at the defaults, 3.96e-3 m² on a
    /// ball of radius 0.928 at a chord of 1e-2).
    #[test]
    fn triangles_by_a_pole_or_an_apex_face_out() {
        let o = |chord_tolerance: f64, max_angle: f64| TessellationOptions {
            chord_tolerance,
            max_angle,
        };
        let tilted = DVec3::new(1.4f64.sin(), 0.2 * 1.4f64.sin(), 1.4f64.cos()).normalize();
        let cap = k()
            .boolean(
                &ball(1.0),
                &slab(DVec3::Z - tilted * 0.01, tilted),
                BooleanOp::Intersect,
            )
            .expect("a cap");
        let piece = k()
            .boolean(
                &ball(0.927_848_640_067_481_9),
                &slab(
                    DVec3::new(
                        -0.187_570_959_161_979_9,
                        0.243_362_344_506_771_45,
                        -0.820_777_160_938_975,
                    ),
                    DVec3::new(
                        -0.586_601_633_680_040_9,
                        0.780_907_379_632_937_2,
                        0.214_667_621_682_280_57,
                    )
                    .normalize(),
                ),
                BooleanOp::Intersect,
            )
            .expect("a piece of a ball");
        let frustum = k()
            .boolean(
                &turned(polygon(&[
                    (0.0, 0.0),
                    (0.475_446_144_542_218_06, 0.0),
                    (0.095_854_760_067_169_53, 0.789_033_539_859_123_3),
                    (0.0, 0.789_033_539_859_123_3),
                ])),
                &slab(
                    DVec3::new(
                        0.118_846_628_415_245_27,
                        -0.052_971_486_938_529_475,
                        0.306_540_880_015_415_13,
                    ),
                    DVec3::new(
                        0.527_975_390_656_859_4,
                        -0.674_893_649_018_185,
                        0.515_519_688_640_167_8,
                    )
                    .normalize(),
                ),
                BooleanOp::Intersect,
            )
            .expect("a cut frustum");
        for (name, body, opts) in [
            ("cap", &cap, TessellationOptions::default()),
            ("cap", &cap, o(1e-4, 0.3)),
            ("cap", &cap, o(1e-4, 1.0)),
            ("piece", &piece, o(1e-2, 1.0)),
            ("frustum", &frustum, o(1e-3, 0.05)),
        ] {
            let (worst, folded) = strays(body, &opts);
            assert_eq!(folded, 0, "{name} {opts:?}");
            assert!(worst <= 1.0, "{name} {opts:?}: {worst}");
        }
    }

    /// Thin wedges of a ball and a torus mesh within the tolerance: their
    /// refinement moves in from the boundary a row at a time, and needed
    /// up to 186 passes where a limit of 64 refused them (found by review).
    #[test]
    fn thin_wedges_of_a_ball_and_a_torus_mesh() {
        let wedge = |r: f64, degrees: f64| {
            turned_by(
                vec![
                    ProfileCurve::Arc {
                        center: DVec2::ZERO,
                        radius: r,
                        start_angle: -FRAC_PI_2,
                        sweep_angle: PI,
                    },
                    ProfileCurve::Line {
                        start: DVec2::new(0.0, r),
                        end: DVec2::new(0.0, -r),
                    },
                ],
                degrees.to_radians(),
            )
        };
        let fine = TessellationOptions {
            chord_tolerance: 1e-5,
            max_angle: 0.5,
        };
        let ring = turned_by(
            vec![ProfileCurve::Circle {
                center: DVec2::new(1.0, 0.0),
                radius: 0.3,
            }],
            0.5f64.to_radians(),
        );
        for (name, body, opts) in [
            ("unit, 1°", wedge(1.0, 1.0), fine),
            (
                "radius 10, 1°",
                wedge(10.0, 1.0),
                TessellationOptions::default(),
            ),
            (
                "radius 3, 1.7°",
                wedge(3.0, 1.7),
                TessellationOptions::default(),
            ),
            ("torus, 0.5°", ring, fine),
        ] {
            let (worst, folded) = strays(&body, &opts);
            assert!(worst <= 1.0, "{name}: {worst}");
            assert_eq!(folded, 0, "{name}");
        }
    }

    /// Angle tolerances asking for hundreds of millions of samples are
    /// refused at once (found by review: a pole's run of samples was
    /// allocated whole and aborted the process).
    #[test]
    fn tiny_angle_tolerances_are_refused() {
        for body in [
            turned(polygon(&[(0.0, 0.0), (1.0, 0.0), (0.0, 1.0)])),
            turned(polygon(&[(0.0, 0.0), (1.0, 0.0), (1.0, 0.1), (0.0, 0.1)])),
        ] {
            for max_angle in [1e-6, 1e-8] {
                let opts = TessellationOptions {
                    chord_tolerance: 1e-4,
                    max_angle,
                };
                match tessellate(&body.brep, &opts, &Tolerance::default()) {
                    Err(KernelError::Tessellation(m)) => assert!(m.contains("more than"), "{m}"),
                    r => panic!("{max_angle}: {:?}", r.map(|m| m.faces.len())),
                }
            }
        }
    }

    /// A slender pin's apex takes the mean of its face's normals there,
    /// which faces along every triangle at the apex (found by review: it
    /// took one meridian's normal, 128° from a triangle across the face).
    #[test]
    fn a_slender_pins_apex_normal_faces_its_triangles() {
        let pin = |at: DVec3| {
            let plane = Frame {
                origin: at,
                x: DVec3::X,
                y: DVec3::Z,
                z: -DVec3::Y,
            };
            k().revolve(
                &Profile {
                    plane,
                    curves: polygon(&[(0.0, 0.0), (0.05, 0.0), (0.0, 1.0)]),
                },
                &RevolveSpec {
                    axis_origin: DVec2::ZERO,
                    axis_dir: DVec2::Y,
                    angle: TAU,
                },
            )
            .expect("a pin")
        };
        for (at, opts) in [
            (
                DVec3::new(10.0, 3.0, 0.0),
                TessellationOptions {
                    chord_tolerance: 1e-5,
                    max_angle: 0.5,
                },
            ),
            (DVec3::ZERO, TessellationOptions::default()),
        ] {
            let body = pin(at);
            let mesh = tessellate(&body.brep, &opts, &Tolerance::default()).expect("a mesh");
            let apex = at + DVec3::Z;
            let mut seen = 0;
            for m in &mesh.faces {
                for t in &m.triangles {
                    let p = t.map(|i| m.positions[i as usize]);
                    let n = (p[1] - p[0]).cross(p[2] - p[0]).normalize();
                    for &i in t {
                        if (m.positions[i as usize] - apex).length() < 1e-12 {
                            seen += 1;
                            let d = n.dot(m.normals[i as usize]);
                            assert!(d > 0.0, "{at}: {d}");
                        }
                    }
                }
            }
            assert!(seen > 0);
        }
    }

    /// A sliver too small to split that faces against the surface is
    /// unfolded by an edge flip, not refused: a frustum a plane cuts near
    /// its apex, nearly tangent to it (the near-tangent apex fuzz's case 19,
    /// whose cut failed validation: "a triangle faces against the surface
    /// where it is too small to split"). Its coarse validation mesh flips
    /// slivers and has no triangle facing against its corners' normals.
    #[test]
    fn folded_slivers_are_unfolded() {
        let axis = Frame::from_normal(DVec3::ZERO, DVec3::ONE.normalize(), DVec3::X);
        let (za, zt) = (1.0, 0.999);
        let plane = Frame {
            origin: axis.origin,
            x: axis.x,
            y: axis.z,
            z: -axis.y,
        };
        let pts =
            [(0.0, 0.0), (1.0, 0.0), (1.0 - zt, zt), (0.0, zt)].map(|(x, y)| DVec2::new(x, y));
        let frustum = k()
            .revolve(
                &Profile::polygon(plane, &pts),
                &RevolveSpec {
                    axis_origin: DVec2::ZERO,
                    axis_dir: DVec2::Y,
                    angle: TAU,
                },
            )
            .expect("a frustum");
        let psi = FRAC_PI_4 * (1.0 + 1e-5);
        let n = axis.x * psi.cos() + axis.z * psi.sin();
        let q = axis.origin + axis.z * za - n * 1e-6;
        let f = Frame::from_normal(q, n, DVec3::Y);
        let block = k()
            .make_box(
                &Frame {
                    origin: q - (f.x + f.y) * 6.0,
                    ..f
                },
                DVec3::splat(12.0),
            )
            .expect("a slab");
        let cut = k()
            .boolean(&frustum, &block, BooleanOp::Subtract)
            .expect("a cut");
        let opts = TessellationOptions {
            chord_tolerance: crate::check::bounding_size(&cut.brep) * 1e-2,
            max_angle: 0.6,
        };
        let mesh = tessellate(&cut.brep, &opts, &Tolerance::default()).expect("a mesh");
        let unfolded: usize = mesh.faces.iter().map(|f| f.unfolded).sum();
        assert!(unfolded > 0);
        assert_eq!(strays(&cut, &opts).1, 0);
    }

    /// Meshes stay within the chord tolerance between the points refinement
    /// tests, against exact distances over a grid of 24ths of every
    /// triangle: a ball (`|1 − |q||`), polar caps of it, a nearly closed
    /// torus, and solids turned from S-shaped splines (the distance to the
    /// profile in the plane through the axis). Found by review: the error
    /// between the points peaks off them (a cubic along an edge at a third
    /// of it); tested at the points and by the quadratic through the edges'
    /// midpoints, an S-curve's mesh strayed 3.6% and a torus's 1.3%, and
    /// without that quadratic a ball's 2.4% (sampled on a grid of
    /// twelfths).
    #[test]
    fn meshes_stay_within_the_tolerance_between_the_points_tested() {
        let o = |c: f64, a: f64| TessellationOptions {
            chord_tolerance: c,
            max_angle: a,
        };
        let worst = |body: &Body, opts: TessellationOptions, dist: &dyn Fn(DVec3) -> f64| {
            let mesh = tessellate(&body.brep, &opts, &Tolerance::default()).expect("a mesh");
            let mut worst = 0.0f64;
            for (f, m) in body.brep.faces.iter().zip(&mesh.faces) {
                // Flat faces (a turned profile's ends are surfaces of
                // revolution) are exact and left out.
                let z = m.positions[0].z;
                if matches!(f.surface, Surface::Plane(_)) || m.positions.iter().all(|p| p.z == z) {
                    continue;
                }
                for t in &m.triangles {
                    let p = t.map(|i| m.positions[i as usize]);
                    let n = 24;
                    for i in 0..=n {
                        for j in 0..=n - i {
                            let (a, b) = (i as f64 / n as f64, j as f64 / n as f64);
                            let q = p[0] * (1.0 - a - b) + p[1] * a + p[2] * b;
                            worst = worst.max(dist(q) / opts.chord_tolerance);
                        }
                    }
                }
            }
            worst
        };
        let sphere = |q: DVec3| (1.0 - q.length()).abs();
        let tilted = DVec3::new(1.4f64.sin(), 0.2 * 1.4f64.sin(), 1.4f64.cos()).normalize();
        let caps = [
            k().boolean(
                &ball(1.0),
                &slab(DVec3::new(0.0, 0.0, 0.9), DVec3::Z),
                BooleanOp::Intersect,
            ),
            k().boolean(
                &ball(1.0),
                &slab(DVec3::Z - tilted * 0.01, tilted),
                BooleanOp::Intersect,
            ),
        ]
        .map(|c| c.expect("a cap"));
        let mut cases: Vec<(String, f64)> = vec![
            ("ball".into(), worst(&ball(1.0), o(1e-3, 0.5), &sphere)),
            ("cap".into(), worst(&caps[0], o(1e-3, 0.5), &sphere)),
            ("tilted cap".into(), worst(&caps[1], o(1e-3, 0.5), &sphere)),
        ];
        let torus = turned(vec![ProfileCurve::Circle {
            center: DVec2::new(1.0, 0.0),
            radius: 0.999,
        }]);
        let ring = |q: DVec3| (DVec2::new(q.truncate().length() - 1.0, q.z).length() - 0.999).abs();
        cases.push(("torus".into(), worst(&torus, o(1e-2, 0.3), &ring)));
        // S-shaped profiles, their height rising along them: the distance
        // to a polyline of 20,000 points of the curve (off it by under
        // 1e-6), searched out from the point's height.
        let k5 = [0.0, 0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0, 1.0];
        for pts in [
            [(0.5, 0.0), (1.2, 0.2), (0.1, 0.5), (1.0, 0.8), (0.5, 1.0)],
            [(0.5, 0.0), (1.5, 0.2), (0.1, 0.5), (1.5, 0.8), (0.5, 1.0)],
        ] {
            let points: Vec<DVec2> = pts.iter().map(|&(x, y)| DVec2::new(x, y)).collect();
            let curve = opendatum_nurbs::NurbsCurve::new(
                3,
                k5.to_vec(),
                points.iter().map(|p| p.extend(0.0)).collect(),
                None,
            )
            .expect("a spline");
            let line: Vec<DVec2> = (0..=20_000)
                .map(|i| curve.point(i as f64 / 20_000.0).truncate())
                .collect();
            assert!(line.windows(2).all(|w| w[1].y > w[0].y));
            let along = |r: DVec2, k: usize| {
                let (a, b) = (line[k], line[k + 1]);
                let t = ((r - a).dot(b - a) / (b - a).length_squared()).clamp(0.0, 1.0);
                (a + (b - a) * t - r).length()
            };
            let profile = |q: DVec3| {
                let r = DVec2::new(q.truncate().length(), q.z);
                let i = line.partition_point(|p| p.y < r.y).clamp(1, line.len() - 1) - 1;
                let mut best = f64::INFINITY;
                for k in (0..=i).rev() {
                    best = best.min(along(r, k));
                    if r.y - line[k].y > best {
                        break;
                    }
                }
                for k in i..line.len() - 1 {
                    best = best.min(along(r, k));
                    if line[k].y - r.y > best {
                        break;
                    }
                }
                best
            };
            let mut curves = vec![ProfileCurve::Line {
                start: DVec2::new(0.0, 0.0),
                end: DVec2::new(0.5, 0.0),
            }];
            curves.push(ProfileCurve::Spline {
                degree: 3,
                knots: k5.to_vec(),
                points,
                weights: None,
            });
            curves.extend(
                polygon(&[(0.5, 1.0), (0.0, 1.0), (0.0, 0.0)])
                    .into_iter()
                    .take(2),
            );
            cases.push((
                format!("S {pts:?}"),
                worst(&turned(curves), o(1e-2, 1.0), &profile),
            ));
        }
        for (name, w) in &cases {
            eprintln!("{name}: {w:.4}");
            assert!(*w <= 1.0, "{name}: {w} tolerances off");
        }
    }

    /// The quartic's peak, the estimate between the points tested (found by
    /// review to have no test of its own), against known peaks and a fine
    /// grid of its values; and the bound that decides when it is searched.
    #[test]
    fn the_quartics_peak_is_its_largest_value() {
        let at = |g: &dyn Fn([f64; 3]) -> f64| QUARTERS.map(|q| g(q.map(|x| x as f64 / 4.0)));
        // A cubic error along an edge peaks at a third of it, 4/27, above
        // all the points tested (9/64 at three quarters).
        let edge = at(&|l| l[1] * l[1] * l[0]);
        let (v, l) = peak(&edge, [false; 3], None);
        assert!(
            (v - 4.0 / 27.0).abs() < 1e-15 && l[2] == 0.0,
            "{v} at {l:?}"
        );
        assert!(edge.iter().all(|&x| x <= 9.0 / 64.0));
        // Left out, the edge's points give way to the largest off it,
        // found by the search within 1/384 of it.
        let (v, l) = peak(&edge, [true, false, false], None);
        assert!(l[2] > 0.0 && v < 4.0 / 27.0 && v > 0.146, "{v} at {l:?}");
        // The bubble peaks at the centroid; a quadratic's equal midpoints
        // at four thirds of them; nothing anywhere is zero.
        let (v, _) = peak(&at(&|l| 27.0 * l[0] * l[1] * l[2]), [false; 3], None);
        assert!((v - 1.0).abs() < 1e-15);
        let (v, _) = peak(
            &at(&|l| 4.0 * (l[0] * l[1] + l[1] * l[2] + l[2] * l[0])),
            [false; 3],
            None,
        );
        assert!((v - 4.0 / 3.0).abs() < 1e-12);
        assert_eq!(peak(&[0.0; 12], [false; 3], None).0, 0.0);
        assert_eq!(GRID.len(), POINTS);
        // Asked whether it passes a value, it stops at a point past it.
        let (v, _) = peak(&edge, [false; 3], Some(0.147));
        assert!(v > 0.147);
        // Random quartics: against the largest on a grid of 150ths (itself
        // a little under the peak), above it by no more than that grid
        // misses, and below it by no more than 0.1% for cubic errors and 5%
        // for random data (several peaks of nearly one height).
        let fine = |at: &[f64; 12]| {
            let n = 150;
            let mut best = 0.0f64;
            for i in 0..=n {
                for j in 0..=n - i {
                    let l = [i, j, n - i - j].map(|x| x as f64 / n as f64);
                    let x: f64 = quartic_basis(l).iter().zip(at).map(|(b, y)| b * y).sum();
                    best = best.max(x.abs());
                }
            }
            best
        };
        let mut rng = crate::boolean::testing::Rng(0x0BEA_C0DE);
        for k in 0..60 {
            let at = if k % 2 == 0 {
                let c: Vec<f64> = (0..10).map(|_| rng.range(-1.0, 1.0)).collect();
                at(&|[a, b, d]: [f64; 3]| {
                    c[0] * a * b
                        + c[1] * b * d
                        + c[2] * d * a
                        + c[3] * a * a * b
                        + c[4] * a * b * b
                        + c[5] * b * b * d
                        + c[6] * b * d * d
                        + c[7] * d * d * a
                        + c[8] * d * a * a
                        + c[9] * a * b * d
                })
            } else {
                std::array::from_fn(|_| rng.range(-1.0, 1.0))
            };
            let (v, _) = peak(&at, [false; 3], None);
            let best = fine(&at);
            let floor = if k % 2 == 0 { 0.999 } else { 0.95 };
            assert!(v.abs() <= best * 1.001, "{k}: {v} > {best}");
            assert!(v.abs() >= best * floor, "{k}: {v} < {best}");
        }
        // The bound: the magnitudes of the Lagrange functions sum to no
        // more than LEBESGUE anywhere (their peak, 3.3536, is near the
        // corners).
        let n = 240;
        let mut most = 0.0f64;
        for i in 0..=n {
            for j in 0..=n - i {
                let l = [i, j, n - i - j].map(|x| x as f64 / n as f64);
                most = most.max(quartic_basis(l).iter().map(|x| x.abs()).sum());
            }
        }
        assert!(most < LEBESGUE && most > 3.35, "{most}");
    }

    fn poly(pts: &[(f64, f64)]) -> Vec<DVec2> {
        pts.iter().map(|&(x, y)| DVec2::new(x, y)).collect()
    }

    /// Coordinates far below a polygon's size are rounding of a zero and
    /// taken as zero: the triangulation refuses nonzero ones under about
    /// 1.8e-43 (a cap's centre came out of its chart at 8e-49).
    #[test]
    fn negligible_coordinates_are_zero() {
        let square = poly(&[(8e-49, 6.8e-33), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]);
        let (cdt, index) = polygon_cdt(&[square]).expect("a triangulation");
        assert_eq!(cdt.num_vertices(), 4);
        assert_eq!(index, vec![0, 1, 2, 3]);
        let tiny = poly(&[(0.0, 0.0), (1e-20, 0.0), (1e-20, 1e-20), (0.0, 1e-20)]);
        assert_eq!(polygon_cdt(&[tiny]).expect("kept").0.num_vertices(), 4);
    }

    #[test]
    fn convexity_is_decided_exactly() {
        let square = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)];
        assert_eq!(convexity(&poly(&square)), Some(true));
        assert_eq!(convexity(&poly(&square[..3])), Some(true));
        let on_side = [(0.0, 0.0), (0.5, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)];
        assert_eq!(convexity(&poly(&on_side)), Some(false));
        let pentagram: Vec<DVec2> = (0..5)
            .map(|k| {
                let a = std::f64::consts::TAU * (2 * k) as f64 / 5.0;
                DVec2::new(a.cos(), a.sin())
            })
            .collect();
        let others = [
            poly(&[(0.0, 0.0), (0.0, 1.0), (1.0, 1.0), (1.0, 0.0)]),
            poly(&[
                (0.0, 0.0),
                (2.0, 0.0),
                (2.0, 1.0),
                (1.0, 1.0),
                (1.0, 2.0),
                (0.0, 2.0),
            ]),
            poly(&[(0.0, 0.0), (2.0, 0.0), (3.0, 0.0), (2.0, 0.0), (2.0, 2.0)]),
            poly(&[(0.0, 0.0), (1.0, 0.0), (1.0, 0.0), (1.0, 1.0)]),
            poly(&[(0.0, 0.0), (1.0, 0.0), (2.0, 0.0)]),
            poly(&[(0.0, 0.0), (1.0, 0.0)]),
            pentagram,
        ];
        for l in &others {
            assert_eq!(convexity(l), None, "{l:?}");
        }
    }

    /// The balanced triangles of a convex polygon cover it once: `n − 2`
    /// counter-clockwise triangles whose areas sum to the polygon's, each
    /// side of the polygon in one triangle and each chord in two.
    #[test]
    fn balanced_triangles_cover_the_polygon() {
        for n in 3..300 {
            let l: Vec<DVec2> = (0..n)
                .map(|k| {
                    let a = std::f64::consts::TAU * k as f64 / n as f64;
                    DVec2::new(a.cos(), 0.5 * a.sin())
                })
                .collect();
            let tris = balanced(n);
            assert_eq!(tris.len(), n - 2);
            let mut area = 0.0;
            let mut uses: FxHashMap<(usize, usize), (u32, u32)> = FxHashMap::default();
            for t in &tris {
                let [a, b, c] = t.map(|k| l[k]);
                let twice = (b - a).perp_dot(c - a);
                assert!(twice > 0.0, "n {n}: {t:?}");
                area += 0.5 * twice;
                for (u, v) in [(t[0], t[1]), (t[1], t[2]), (t[2], t[0])] {
                    let e = uses.entry((u.min(v), u.max(v))).or_default();
                    if u < v {
                        e.0 += 1;
                    } else {
                        e.1 += 1;
                    }
                }
            }
            let shoelace: f64 = (0..n).map(|k| 0.5 * l[k].perp_dot(l[(k + 1) % n])).sum();
            assert!((area - shoelace).abs() < 1e-12, "n {n}");
            for ((u, v), (fw, bw)) in uses {
                let side = v - u == 1 || (u == 0 && v == n - 1);
                if side {
                    assert_eq!(fw + bw, 1, "n {n}: side {u}-{v}");
                } else {
                    assert_eq!((fw, bw), (1, 1), "n {n}: chord {u}-{v}");
                }
            }
        }
    }
}
