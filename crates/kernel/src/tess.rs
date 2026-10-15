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
//!    estimated between them by the quartic through them, bounded by its
//!    Bernstein-Bézier coefficients), faces the way the surface does at
//!    its corners and centroid (unless its corners lie in a line within
//!    the model tolerance, when it has no side to face), and, if a chord
//!    long or longer, its normals turn by less than the angle tolerance. A
//!    triangle facing against the surface beside a boundary segment that
//!    is its longest side has the segment bisected and the face meshed
//!    again; one too small to split that faces against it is mended by an
//!    edge flip, or else split at its centroid and refined again. A sliver
//!    too thin for either, or a face that cannot be refined or would pass
//!    the vertex budget, is an error.

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
    /// Triangles too small to split that faced against the surface,
    /// mended by an edge flip (design note 0002 §4; counted for the tests).
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
    // §4).
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

/// A point measured against a surface: its offset to the surface point
/// under it, and the surface's unit normal there (`None` at a singular
/// point).
type Probe = (DVec3, Option<DVec3>);

/// A mesh vertex on a curved face: its point on the surface, the
/// surface's unit normal there (`None` at a singular point), and which
/// parameter is free there (`Some(true)` for `u`; [`singular_dir`]).
type Corner = (DVec3, Option<DVec3>, Option<bool>);

/// The points measured along mesh edges, each edge's three ([`Probe`]s at
/// a quarter, half and three quarters of it from its lesser corner), kept
/// for the other triangle beside it and for a triangle judged again: a
/// direct-mapped cache, one slot per edge by a hash of its corners,
/// overwritten on a collision, and sized to the face's vertices up to
/// [`EdgeProbes::MOST`] slots. What a slot holds depends only on its edge,
/// so a miss changes no result. A hash table of every edge measured
/// (until 2026-10-15) grew to over a hundred thousand entries a face on a
/// unit ball, tens of megabytes; this cache, of 754 kB, held 42% of the
/// lookups of the ball's default mesh, of the 49% that can be held
/// (docs/PERFORMANCE.md, fifth pass).
struct EdgeProbes {
    slots: Vec<([usize; 2], [Probe; 3])>,
}

impl EdgeProbes {
    /// The most slots: 754 kB.
    const MOST: usize = 1 << 12;

    fn new() -> Self {
        Self { slots: Vec::new() }
    }

    /// Room for a face of `vertices` vertices: about two slots a vertex
    /// (a mesh has about three edges a vertex), cleared when it grows.
    fn fit(&mut self, vertices: usize) {
        let want = (2 * vertices).next_power_of_two().clamp(64, Self::MOST);
        if self.slots.len() < want {
            self.slots = vec![([usize::MAX; 2], [(DVec3::ZERO, None); 3]); want];
        }
    }

    /// The probes of edge `key` (its corners, lesser first), measured by
    /// `measure` unless held.
    fn get(&mut self, key: [usize; 2], measure: impl FnOnce() -> [Probe; 3]) -> [Probe; 3] {
        let h = (key[0] as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15)
            ^ (key[1] as u64).wrapping_mul(0xC2B2_AE3D_27D4_EB4F);
        let n = self.slots.len();
        let slot = &mut self.slots[(h >> 32) as usize & (n - 1)];
        if slot.0 != key {
            *slot = (key, measure());
        }
        slot.1
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
    singular_of(su, sv)
}

/// Which parameter, if either, is free where a surface's partial
/// derivatives are `su` and `sv` ([`singular_dir`]).
fn singular_of(su: DVec3, sv: DVec3) -> Option<bool> {
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
                        let total = steps + (count + fixed.len()) as f64;
                        if total.is_nan() || total > VERTEX_BUDGET as f64 {
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
/// space, so the mesh is held to the tolerance, as estimated there, either
/// way.
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

/// The angle between unit vectors `a` and `b`, from above and without an
/// inverse cosine: twice `asin(s)` for `s` half their distance, which is
/// at most `s / √(1 − s²)` (its tangent; about `s³/3` above it for small
/// `s`), capped at π.
fn turn(a: DVec3, b: DVec3) -> f64 {
    let s2 = 0.25 * (a - b).length_squared();
    if s2 >= 1.0 {
        return std::f64::consts::PI;
    }
    2.0 * (s2 / (1.0 - s2)).sqrt().min(std::f64::consts::FRAC_PI_2)
}

/// The points of the quartic Lagrange triangle other than its corners, in
/// quarters of their barycentric coordinates: a quarter, half and three
/// quarters of the way along each edge `01`, `12` and `20`, then the three
/// inside, at (½, ¼, ¼) and its permutations: a quarter of the way from the
/// centroid to each corner.
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
#[cfg(test)]
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

/// The exponents of the quartic Bernstein polynomials of a triangle,
/// `4!/(i! j! k!) λ₀ⁱ λ₁ʲ λ₂ᵏ`: its corners', then the [`QUARTERS`]'. A
/// quarter of each is where its coefficient sits.
const EXPONENTS: [[usize; 3]; 15] = [
    [4, 0, 0],
    [0, 4, 0],
    [0, 0, 4],
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

/// The quartic Bernstein polynomials, as [`EXPONENTS`], at barycentric
/// coordinates `l`.
fn bernstein(l: [f64; 3]) -> [f64; 15] {
    const FACTORIAL: [f64; 5] = [1.0, 1.0, 2.0, 6.0, 24.0];
    EXPONENTS.map(|e| (0..3).fold(24.0, |x, a| x * l[a].powi(e[a] as i32) / FACTORIAL[e[a]]))
}

/// The corners of the four triangles that splitting one at its edges'
/// middles makes, in its barycentric coordinates.
const HALVES: [[[f64; 3]; 3]; 4] = [
    [[1.0, 0.0, 0.0], [0.5, 0.5, 0.0], [0.5, 0.0, 0.5]],
    [[0.5, 0.5, 0.0], [0.0, 1.0, 0.0], [0.0, 0.5, 0.5]],
    [[0.5, 0.0, 0.5], [0.0, 0.5, 0.5], [0.0, 0.0, 1.0]],
    [[0.0, 0.5, 0.5], [0.5, 0.0, 0.5], [0.5, 0.5, 0.0]],
];

/// The Bernstein-Bézier form of the quartic estimate: its coefficients
/// from its values at the [`QUARTERS`] (zero at the corners), and each of
/// the [`HALVES`]' coefficients from its parent's. Stored by column, each
/// input's contribution to every output, so the products vectorize.
struct Bezier {
    from_values: [[f64; 15]; 12],
    halves: [[[f64; 15]; 15]; 4],
}

impl Bezier {
    /// The coefficients of the quartic that is zero at the corners and
    /// `at` the [`QUARTERS`].
    fn coefficients(&self, at: &[f64; 12]) -> [f64; 15] {
        product(&self.from_values, at)
    }

    /// The coefficients of half `h` ([`HALVES`]) of the quartic whose
    /// coefficients are `c`.
    fn half(&self, h: usize, c: &[f64; 15]) -> [f64; 15] {
        product(&self.halves[h], c)
    }
}

/// The matrix whose columns are `columns` times `x`.
fn product<const N: usize>(columns: &[[f64; 15]; N], x: &[f64; N]) -> [f64; 15] {
    let mut out = [0.0; 15];
    for (column, &x) in columns.iter().zip(x) {
        for (o, m) in out.iter_mut().zip(column) {
            *o += m * x;
        }
    }
    out
}

static BEZIER: std::sync::LazyLock<Bezier> = std::sync::LazyLock::new(|| {
    // The Bernstein polynomials at the points their coefficients sit at,
    // inverted by Gauss-Jordan elimination with partial pivoting (a fixed
    // matrix, well conditioned; a test checks the result).
    let node = |e: [usize; 3]| e.map(|x| x as f64 / 4.0);
    let mut a: Vec<[f64; 30]> = EXPONENTS
        .iter()
        .enumerate()
        .map(|(r, &e)| {
            let mut row = [0.0; 30];
            row[..15].copy_from_slice(&bernstein(node(e)));
            row[15 + r] = 1.0;
            row
        })
        .collect();
    for c in 0..15 {
        let p = (c..15)
            .max_by(|&x, &y| a[x][c].abs().total_cmp(&a[y][c].abs()))
            .unwrap_or(c);
        a.swap(c, p);
        let pivot = a[c][c];
        for x in a[c].iter_mut() {
            *x /= pivot;
        }
        let row = a[c];
        for (r, other) in a.iter_mut().enumerate() {
            if r != c {
                let f = other[c];
                for (x, y) in other.iter_mut().zip(row) {
                    *x -= f * y;
                }
            }
        }
    }
    let inverse: [[f64; 15]; 15] = std::array::from_fn(|r| std::array::from_fn(|c| a[r][15 + c]));
    let column = QUARTERS.map(|q| EXPONENTS.iter().position(|&e| e == q).unwrap_or(0));
    let from_values = column.map(|c| std::array::from_fn(|r| inverse[r][c]));
    let halves = HALVES.map(|corners| {
        // The parent's Bernstein polynomials at the points of the half's
        // coefficients.
        let at: Vec<[f64; 15]> = EXPONENTS
            .iter()
            .map(|&e| {
                let l = node(e);
                bernstein(std::array::from_fn(|a| {
                    (0..3).map(|v| l[v] * corners[v][a]).sum()
                }))
            })
            .collect();
        std::array::from_fn(|c| {
            std::array::from_fn(|r| (0..15).map(|n| inverse[r][n] * at[n][c]).sum())
        })
    });
    Bezier {
        from_values,
        halves,
    }
});

/// How many parts [`past`] splits a triangle into, at most, before it
/// takes a quartic it can neither bound under the tolerance nor find past
/// it as past it.
const PARTS: usize = 256;

/// Where, in barycentric coordinates, the quartic that is zero at a
/// triangle's corners and `at` the [`QUARTERS`] passes `tol` in magnitude,
/// or `None` where it is proven not to: the quartic Lagrange interpolant
/// of the distance to the surface, the estimate between the points tested
/// (after R. E. Bank and R. K. Smith, "Mesh smoothing using a posteriori
/// error estimates", SIAM J. Numer. Anal. 34, 1997, who estimate a linear
/// element's error from the quadratic one).
///
/// A quartic lies within the hull of its Bernstein-Bézier coefficients
/// (G. Farin, "Triangular Bernstein-Bézier patches", Computer Aided
/// Geometric Design 3, 1986), so none past `tol` proves it under it. Where
/// one is, the triangle is split at its edges' middles, the part with the
/// largest coefficients searched first, until a corner of a part is past
/// `tol` (the coefficients at the corners are the values there) or every
/// part is proven under it. Past [`PARTS`] parts the last one's centroid
/// is taken as past, so the answer errs only toward refining (design note
/// 0002 §4). On the edges `01`, `12` and `20` that `skip` marks (segments
/// of the boundary, which no point added inside moves) the values and
/// coefficients are left out: inside, the quartic then passes `tol` by no
/// more than it does on those edges, less away from them.
fn past(at: &[f64; 12], skip: [bool; 3], tol: f64) -> Option<[f64; 3]> {
    let bz = &*BEZIER;
    let skipping = skip.contains(&true);
    // Whether coefficient `m` of a part with these corners sits on an edge
    // left out (the corners' coordinates are sums of halves, so exact).
    let left_out = |m: usize, corners: &[[f64; 3]; 3]| {
        skipping
            && (0..3).any(|e| {
                skip[e] && (0..3).all(|v| EXPONENTS[m][v] == 0 || corners[v][(e + 2) % 3] == 0.0)
            })
    };
    let bound = |c: &[f64; 15], corners: &[[f64; 3]; 3]| {
        if skipping {
            (0..15)
                .filter(|&m| !left_out(m, corners))
                .fold(0.0f64, |b, m| b.max(c[m].abs()))
        } else {
            c.iter().fold(0.0f64, |b, x| b.max(x.abs()))
        }
    };
    let mut stack: Vec<([f64; 15], [[f64; 3]; 3])> = Vec::with_capacity(16);
    stack.push((
        bz.coefficients(at),
        [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
    ));
    let mut parts = 0;
    while let Some((c, corners)) = stack.pop() {
        if let Some(v) = (0..3).find(|&v| c[v].abs() > tol && !left_out(v, &corners)) {
            return Some(corners[v]);
        }
        if bound(&c, &corners) <= tol {
            continue;
        }
        parts += 1;
        if parts > PARTS {
            return Some(std::array::from_fn(|a| {
                (corners[0][a] + corners[1][a] + corners[2][a]) / 3.0
            }));
        }
        let mut halves: [(f64, [f64; 15], [[f64; 3]; 3]); 4] = std::array::from_fn(|h| {
            let ch = bz.half(h, &c);
            let corners: [[f64; 3]; 3] = std::array::from_fn(|v| {
                std::array::from_fn(|a| (0..3).map(|w| HALVES[h][v][w] * corners[w][a]).sum())
            });
            (bound(&ch, &corners), ch, corners)
        });
        // The largest last, to be searched first.
        halves.sort_by(|x, y| x.0.total_cmp(&y.0));
        stack.extend(halves.map(|(_, c, corners)| (c, corners)));
    }
    None
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
    // Per CDT vertex of a curved face: its surface point, its normal, and
    // which parameter is free there ([`singular_dir`]), from one
    // evaluation.
    let mut at: Vec<Corner> = Vec::new();
    // At a singular point (an apex, a pole, a disc's centre) `S_u × S_v`
    // is the rounding of a vanishing derivative and can point either way:
    // such a point has no normal of its own.
    let evaluate = |uv: DVec2, p: Option<DVec3>| -> Corner {
        let (q, su, sv) = s.point_derivs(uv.x, uv.y);
        let free = singular_of(su, sv);
        let n = free.is_none().then(|| crate::geom::unit_normal(su, sv));
        (p.unwrap_or(q), n.flatten(), free)
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
            at.resize(cdt.num_vertices(), (DVec3::ZERO, None, None));
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
                      at: &mut Vec<Corner>,
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
            // A corner's normal is past the angle tolerance from another
            // where their dot product is under this.
            let cos_max = opts.max_angle.min(std::f64::consts::PI).cos();
            // How triangle `t` fares (design note 0002 §4): fine, to be
            // split at a point, or beside a boundary segment to bisect
            // first.
            let judge = |t: [usize; 3],
                         verts: &[(DVec2, Option<DVec3>)],
                         at: &[Corner],
                         boundary: &FxHashMap<[usize; 2], Option<(usize, usize)>>,
                         min_area: f64,
                         probes: &mut EdgeProbes|
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
                // A triangle facing away from the surface's normal at its
                // centroid or at a corner (a long sliver folded over, or
                // one standing on edge across the surface), or a normal
                // turning too far. Facing is judged at the corners too: a
                // sliver standing nearly square to a cone faced along the
                // normal at its centroid in the parameters and against two
                // of its corners' and the one at its closest point (found
                // by review, 0.45 mm high in a validation mesh). The angle
                // tolerance is held on triangles as long as the chord
                // tolerance or longer, as on the edges' curves (§1): a
                // surface turning through more than it within the chord
                // tolerance (nearly a crease) is followed by the chord
                // alone there (found by review: a revolved rational cubic
                // with a weight of 30 had its seam bisected for 35 rounds,
                // 96,000 triangles and 18 s, and with a weight of 60 at 5°
                // for minutes and then refused).
                let (pc, su, sv) = s.point_derivs(c.x, c.y);
                let nc = crate::geom::unit_normal(su, sv);
                let folded = !thin
                    && (nc.is_some_and(|nc| tri_n.dot(nc) <= 0.0)
                        || t.iter()
                            .any(|&i| at[i].1.is_some_and(|n| tri_n.dot(n) <= 0.0)));
                let long = side_3d[0].max(side_3d[1]).max(side_3d[2]) >= opts.chord_tolerance;
                let turned = folded
                    || long
                        && nc.is_some_and(|nc| {
                            t.iter()
                                .any(|&i| at[i].1.is_some_and(|nw| nw.dot(nc) < cos_max))
                        });
                if folded {
                    // Beside a boundary segment that is its longest side,
                    // points added inside close in on the segment's middle
                    // in the parameters, which can lie off its chord: the
                    // segment is bisected (a segment Delaunay refinement
                    // splits rather than crowd: J. Ruppert, "A Delaunay
                    // refinement algorithm for quality 2-dimensional mesh
                    // generation", J. Algorithms 18, 1995; found by review,
                    // as triangles folded over by a cut near a ball's pole
                    // and by an apex, settled once too small to split).
                    // Only a folded triangle asks: each bisection meshes the
                    // face again, and triangles whose normals merely turn
                    // past the angle tolerance asked round after round
                    // beside a seam (review: 35 rounds, 96,000 triangles and
                    // 18 s for a revolved rational cubic with a weight of
                    // 30, which now meshes to 13,830 triangles in 82 ms).
                    let e = (0..3)
                        .max_by(|&a, &b| lens[a].total_cmp(&lens[b]))
                        .unwrap_or(0);
                    if let Some(Some(seg)) = side(e) {
                        return Ok(Judged::Bisect(*seg));
                    }
                }
                // Tiny triangles and slivers are not split (their centroid
                // could land on a boundary constraint within round-off).
                // One folded over is unfolded once refinement is done.
                if area < min_area || 2.0 * area < 1e-8 * longest {
                    return Ok(if folded { Judged::Folded } else { Judged::Fine });
                }
                // Chord error: the distance to the surface from the points
                // at a quarter, half and three quarters of the triangle's
                // edges, its centroid and the points a quarter of the way
                // from it to the corners, and the largest value inside it
                // of the quartic through those points (zero at the
                // corners). The longest edge, failing, is split at its
                // midpoint (after Rivara's longest-edge bisection, though
                // the point goes into the Delaunay triangulation, so his
                // bound on angles does not strictly hold); a failure
                // elsewhere adds the centroid. Estimated, not bounded:
                // where the surface bends faster than a quartic between the
                // points, the mesh can stray further, by an amount nothing
                // here bounds (design note 0002 §4).
                let tol = opts.chord_tolerance;
                // How fast the surface bends across the triangle: the
                // largest turn of its normal per unit length between its
                // corners, and from them to its centroid.
                let bend = {
                    let mid = (p[0] + p[1] + p[2]) / 3.0;
                    let mut bend = 0.0f64;
                    for i in 0..3 {
                        let Some(a) = at[t[i]].1 else {
                            continue;
                        };
                        for j in i + 1..3 {
                            if let Some(b) = at[t[j]].1 {
                                bend = bend.max(turn(a, b) / (p[i] - p[j]).length());
                            }
                        }
                        if let Some(nc) = nc {
                            bend = bend.max(turn(a, nc) / (p[i] - mid).length());
                        }
                    }
                    bend
                };
                // The distance to the surface from a point of the triangle,
                // from its offset `off` to the surface point under it and
                // the surface's normal there, signed by the side of the
                // triangle the surface is on: the offset along the normal
                // and, for the offset along the surface (the parameters'
                // slack, not the mesh's error), how far the surface can
                // bend toward or away from the point across it, `bend`
                // times half its square. Along the normal alone (until
                // 2026-10-15), the distance came out short on a saddle,
                // which bends toward the side the slack lies on: review
                // found a torus nearly closed at its middle meshed 1.7%
                // past the tolerance.
                let gap = |&(off, n): &Probe| match n {
                    Some(n) => {
                        let h = off.dot(n);
                        let slack = (off.length_squared() - h * h).max(0.0);
                        (h.abs() + 0.5 * bend * slack).copysign(h * n.dot(tri_n))
                    }
                    None => off.dot(tri_n),
                };
                // The parameter of a corner that is free: at an apex or a
                // pole, always a boundary sample.
                let free = |i: usize| verts[i].1.and(at[i].2);
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
                let probe = |w: DVec2, q: DVec3| -> Probe {
                    let (x, su, sv) = s.point_derivs(w.x, w.y);
                    (x - q, crate::geom::unit_normal(su, sv))
                };
                // The distances at the quartic triangle's points
                // ([`QUARTERS`]), the edges' first, the longest first. An
                // edge's are measured from its lesser corner and kept for
                // the triangle beside it ([`EdgeProbes`]); the first failing
                // edge decides the split without the rest.
                let mut d = [0.0; 12];
                let mut order = [0, 1, 2];
                order.sort_by(|&a, &b| lens[b].total_cmp(&lens[a]));
                for e in order {
                    let (a, b) = (t[e], t[(e + 1) % 3]);
                    let (x, y) = (a.min(b), a.max(b));
                    let mut along = probes
                        .get([x, y], || {
                            [0.25, 0.5, 0.75]
                                .map(|f| probe(under(x, y, f), at[x].0.lerp(at[y].0, f)))
                        })
                        .map(|pr| gap(&pr));
                    if a > b {
                        along.reverse();
                    }
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
                    if along.iter().any(|x| x.abs() > tol) {
                        let m = (verts[a].0 + verts[b].0) * 0.5;
                        return Ok(Judged::Insert(if lens[e] >= longest { m } else { c }));
                    }
                }
                // The centroid (measured already), then the points a
                // quarter of the way from it to the corners.
                if gap(&(pc - (p[0] + p[1] + p[2]) / 3.0, nc)).abs() > tol {
                    return Ok(Judged::Insert(c));
                }
                for (k, w) in [[0.5, 0.25, 0.25], [0.25, 0.5, 0.25], [0.25, 0.25, 0.5]]
                    .into_iter()
                    .enumerate()
                {
                    let x = uv[0] * w[0] + uv[1] * w[1] + uv[2] * w[2];
                    d[9 + k] = gap(&probe(x, p[0] * w[0] + p[1] * w[1] + p[2] * w[2]));
                }
                if d[9..].iter().any(|x| x.abs() > tol) {
                    return Ok(Judged::Insert(c));
                }
                let mut split = None;
                // Between the points: the quartic through them, bounded
                // where it could pass the tolerance, leaving out the
                // boundary's segments (no point added inside moves them).
                // Past it on an edge that is the triangle's longest, that
                // edge is bisected; elsewhere the centroid is added.
                let large = d.iter().fold(0.0f64, |m, x| m.max(x.abs()));
                if large * LEBESGUE > tol {
                    let on = |e: usize, l: [f64; 3]| l[(e + 2) % 3] < 1e-9;
                    if let Some(l) = past(&d, [0, 1, 2].map(|e| side(e).is_some()), tol) {
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
            let mut probes = EdgeProbes::new();
            let mut splits: Vec<(usize, usize)> = Vec::new();
            let mut folded: Vec<[usize; 3]> = Vec::new();
            // The triangles about the points from `from` to `to` (inserting
            // a point re-triangulates only the triangles about it), and
            // those of `failing` still there; each once, from its least
            // corner, in the order found (the points' order, then each
            // one's edges' about it: a hash set only says which are taken,
            // so the result does not depend on its order). Triangles about
            // one point come together and share their edges about it,
            // which [`EdgeProbes`] then holds (sorted by their corners,
            // until 2026-10-15, the two triangles beside an edge came far
            // apart: on a ball's default mesh, repeated measurements came
            // to a third of its edges, an eighth now).
            let about = |cdt: &Cdt, from: usize, to: usize, failing: Vec<[usize; 3]>| {
                let mut next: Vec<[usize; 3]> = Vec::new();
                let mut taken: FxHashSet<[usize; 3]> = FxHashSet::default();
                let mut take = |t: [usize; 3]| {
                    let m = (0..3).min_by_key(|&i| t[i]).unwrap_or(0);
                    let t = [t[m], t[(m + 1) % 3], t[(m + 2) % 3]];
                    if taken.insert(t) {
                        next.push(t);
                    }
                };
                for v in from..to {
                    let h = spade::handles::FixedVertexHandle::from_index(v);
                    for e in cdt.vertex(h).out_edges() {
                        if let Some(f) = e.face().as_inner() {
                            take(f.vertices().map(|x| x.fix().index()));
                        }
                    }
                }
                for t in failing {
                    if present(cdt, t) {
                        take(t);
                    }
                }
                next
            };
            // Refinement, then unfolding; a folded triangle that no flip
            // mends is split and refinement resumes.
            let tris = loop {
                while !work.is_empty() {
                    probes.fit(verts.len());
                    let mut new_pts: Vec<(DVec2, usize)> = Vec::new();
                    let mut failing: Vec<[usize; 3]> = Vec::new();
                    for &t in &work {
                        if !present(&cdt, t) {
                            continue;
                        }
                        match judge(t, &verts, &at, &boundary, min_area, &mut probes)? {
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
                    work = about(&cdt, before, verts.len(), failing);
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
                // the test above. One that no flip mends is split at its
                // centroid, unless it is a sliver thinner than 1e-8 of its
                // length in the parameters, which is refused: tiny by area
                // alone, an ear between two boundary segments by a cone's
                // apex, 12 mm long at the coarse chord of a result's check,
                // could be neither flipped nor split, and the valid union was
                // refused (found by review: 23 of 300 unions of a frustum with
                // a slab nearly tangent to it by its apex).
                let mut tris = inside_triangles(&cdt);
                probes.fit(verts.len());
                let mut flips = 0;
                let mut stuck: Vec<(DVec2, usize)> = Vec::new();
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
                    for &f in &folded {
                        let Some(&i) = slot.get(&canon(f)) else {
                            continue;
                        };
                        // Still there, it is still folded: a triangle's
                        // judgement depends only on its corners, which
                        // never move. (It was judged again here until
                        // 2026-10-15, in the triangulation's rotation of
                        // its corners rather than the one judged, and
                        // skipped unless found folded again: review found
                        // a tie between side lengths or an area within
                        // rounding of a limit could then have let it
                        // through unmended.)
                        let t = tris[i];
                        let length = |e: usize| {
                            (verts[t[(e + 1) % 3]].0 * scale - verts[t[e]].0 * scale)
                                .length_squared()
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
                            let mut fine = |n: [usize; 3]| -> Result<bool, Unmeshed> {
                                Ok(ccw(n)
                                    && matches!(
                                        judge(n, &verts, &at, &boundary, min_area, &mut probes)?,
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
                                flips += 1;
                                break;
                            }
                        }
                        if !mended {
                            let q = t.map(|i| verts[i].0 * scale);
                            let twice = (q[1] - q[0]).perp_dot(q[2] - q[0]).abs();
                            if twice < 1e-8 * (0..3).map(length).fold(0.0, f64::max) {
                                return Err(
                                    "a triangle faces against the surface where it is too thin \
                                            to split or flip"
                                        .to_string()
                                        .into(),
                                );
                            }
                            stuck.push((
                                (verts[t[0]].0 + verts[t[1]].0 + verts[t[2]].0) / 3.0,
                                t[0],
                            ));
                        }
                    }
                }
                if stuck.is_empty() {
                    unfolded += flips;
                    break tris;
                }
                if verts.len() + stuck.len() > VERTEX_BUDGET {
                    return Err("tessellation exceeded the vertex budget".to_string().into());
                }
                let before = verts.len();
                for (c, hint) in stuck {
                    insert(&mut cdt, &mut verts, &mut at, c, hint)?;
                }
                if verts.len() == before {
                    return Err("splitting a folded triangle added no point"
                        .to_string()
                        .into());
                }
                work = about(&cdt, before, verts.len(), Vec::new());
            };
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
        // Every cut is done (design note 0002 §4 counts them).
        assert_eq!(bodies.len(), 46);
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

    /// A solid turned from one rational cubic whose middle weights are far
    /// apart bends nearly to a crease. It meshes within the chord
    /// tolerance, facing out, the angle tolerance held on triangles a chord
    /// long or longer (found by review: the angle tolerance was held at any
    /// size, and triangles beside the seam whose normals turned past it
    /// had the seam bisected, the face meshed again each time, round after
    /// round: 35 rounds and 96,000 triangles in 18 s at a weight of 30, and
    /// at 60 and 5° a refusal after minutes; now 13,830 and 20,600
    /// triangles, in 82 and 100 ms. The time guards keep 60 at 5°, 92,264
    /// triangles).
    #[test]
    fn a_nearly_creased_rational_revolve_meshes() {
        let pts = [(0.2, 0.0), (1.2, 0.2), (1.2, 1.8), (0.2, 2.0)].map(|(x, y)| DVec2::new(x, y));
        for w in [30.0, 60.0] {
            let body = turned(vec![
                ProfileCurve::Spline {
                    degree: 3,
                    knots: vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0],
                    points: pts.to_vec(),
                    weights: Some(vec![1.0, w, 0.2, 1.0]),
                },
                ProfileCurve::Line {
                    start: pts[3],
                    end: pts[0],
                },
            ]);
            let opts = TessellationOptions {
                chord_tolerance: 1e-3,
                max_angle: 15f64.to_radians(),
            };
            let mesh = tessellate(&body.brep, &opts, &Tolerance::default()).expect("a mesh");
            let triangles: usize = mesh.faces.iter().map(|f| f.triangles.len()).sum();
            assert!(triangles < 25_000, "{w}: {triangles} triangles");
            let (worst, folded) = strays(&body, &opts);
            assert!(worst <= 1.0, "{w}: {worst}");
            assert_eq!(folded, 0, "{w}");
        }
    }

    /// An extruded zig-zag spline, whose bends are tighter than the chord
    /// tolerance, meshes with about as many triangles as its samples need,
    /// within the tolerance and facing out (found by review: its edges'
    /// pieces shorter than the chord were let turn past the angle
    /// tolerance, and the triangles along the rulings beside them were
    /// split until tiny: 19,784 triangles where 216 do at 1e-2 and 0.5
    /// rad, 379,510 in 11 s at 1e-3 and 0.1, and a refusal at 1e-2 and
    /// 0.1).
    #[test]
    fn an_extruded_zig_zag_meshes_to_its_samples() {
        let n = 10;
        let mut points = vec![DVec2::new(0.0, 0.5)];
        for i in 1..n {
            let y = if i % 2 == 1 { 1.5 } else { 0.05 };
            points.push(DVec2::new(i as f64 / n as f64, y));
        }
        points.push(DVec2::new(1.0, 0.5));
        let mut knots = vec![0.0; 4];
        knots.extend((1..n - 2).map(|i| i as f64 / (n - 2) as f64));
        knots.extend([1.0; 4]);
        let corners =
            [(1.0, 0.5), (1.0, -0.5), (0.0, -0.5), (0.0, 0.5)].map(|(x, y)| DVec2::new(x, y));
        let mut curves = vec![ProfileCurve::Spline {
            degree: 3,
            knots,
            points,
            weights: None,
        }];
        curves.extend((0..3).map(|i| ProfileCurve::Line {
            start: corners[i],
            end: corners[i + 1],
        }));
        let body = k()
            .extrude(
                &Profile {
                    plane: Frame::WORLD,
                    curves,
                },
                &ExtrudeSpec {
                    distance: 0.7,
                    symmetric: false,
                },
            )
            .expect("an extrusion");
        for (chord_tolerance, max_angle, most) in [
            (1e-2, 0.5, 400),
            (1e-2, 0.2, 800),
            (1e-2, 0.1, 1_500),
            (1e-3, 0.1, 2_000),
        ] {
            let opts = TessellationOptions {
                chord_tolerance,
                max_angle,
            };
            let mesh = tessellate(&body.brep, &opts, &Tolerance::default()).expect("a mesh");
            let curved = body
                .brep
                .faces
                .iter()
                .zip(&mesh.faces)
                .filter(|(f, _)| !matches!(f.surface, Surface::Plane(_)))
                .map(|(_, m)| m.triangles.len());
            for triangles in curved {
                assert!(
                    triangles <= most,
                    "{chord_tolerance} {max_angle}: {triangles}"
                );
            }
            let (worst, folded) = strays(&body, &opts);
            assert!(worst <= 1.0, "{chord_tolerance} {max_angle}: {worst}");
            assert_eq!(folded, 0, "{chord_tolerance} {max_angle}");
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

    /// Slabs nearly tangent to a frustum along a ruling, by its apex. A
    /// union's own check meshes it at a hundredth of its size, where an ear
    /// between two boundary segments, 12 mm long but tiny in area, folded
    /// over and could be neither split nor flipped, so the valid union was
    /// refused as invalid; and cuts' meshes held slivers standing on edge,
    /// facing along the normal at their centroids in the parameters but
    /// against their corners' (found by review: 23 of 300 such unions
    /// refused, 10 slivers in 1,800 meshes). Each operation is done or
    /// refused as unsupported, and in its meshes at a hundredth of its
    /// size and at the defaults every triangle thicker than the model
    /// tolerance faces along its corners' normals (thinner ones have no
    /// side to face, §4).
    #[test]
    fn slivers_by_a_frustums_apex_face_out() {
        // A frustum about `axis`, `s` across at its base, its top a
        // hundredth of its apex's height below it, and a slab whose face is
        // tilted `f` of the cone's angle from tangent along the ruling at
        // azimuth `phi`, `d` off the apex.
        let pair = |axis: Frame, slope: f64, phi: f64, f: f64, d: f64, s: f64| {
            let top = 0.99 / slope;
            let plane = Frame {
                origin: axis.origin,
                x: axis.x,
                y: axis.z,
                z: -axis.y,
            };
            let pts = [(0.0, 0.0), (1.0, 0.0), (1.0 - slope * top, top), (0.0, top)]
                .map(|(r, h)| DVec2::new(r, h) * s);
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
            let apex = axis.origin + axis.z * (s / slope);
            let psi = slope.atan() * (1.0 + f);
            let n = (axis.x * phi.cos() + axis.y * phi.sin()) * psi.cos() + axis.z * psi.sin();
            let q = apex - n * d;
            let g = Frame::from_normal(q, n, DVec3::Y);
            let slab = k()
                .make_box(
                    &Frame {
                        origin: q - (g.x + g.y) * (6.0 * s),
                        ..g
                    },
                    DVec3::splat(12.0 * s),
                )
                .expect("a slab");
            (frustum, slab)
        };
        let facing_in = |body: &Body, opts: &TessellationOptions| {
            let mesh = tessellate(&body.brep, opts, &Tolerance::default()).expect("a mesh");
            let mut count = 0;
            for m in &mesh.faces {
                for t in &m.triangles {
                    let p = t.map(|i| m.positions[i as usize]);
                    let n = (p[1] - p[0]).cross(p[2] - p[0]);
                    let longest = (0..3)
                        .map(|e| (p[(e + 1) % 3] - p[e]).length())
                        .fold(0.0, f64::max);
                    if n.length() >= 1e-7 * longest {
                        count += t.iter().any(|&i| n.dot(m.normals[i as usize]) <= 0.0) as usize;
                    }
                }
            }
            count
        };
        let check = |a: &Body, b: &Body, what: &str| {
            for op in [BooleanOp::Union, BooleanOp::Subtract, BooleanOp::Intersect] {
                match k().boolean(a, b, op) {
                    Ok(r) => {
                        let size = crate::check::bounding_size(&r.brep);
                        for opts in [
                            TessellationOptions {
                                chord_tolerance: size * 1e-2,
                                max_angle: 0.6,
                            },
                            TessellationOptions::default(),
                        ] {
                            let n = facing_in(&r, &opts);
                            assert_eq!(n, 0, "{what} {op:?} at {}", opts.chord_tolerance);
                        }
                    }
                    Err(KernelError::Unsupported(_)) => {}
                    // A slab clear of the frustum leaves nothing in common,
                    // and the cut takes nothing from it.
                    Err(KernelError::InvalidInput(m))
                        if op == BooleanOp::Intersect && m.contains("empty") =>
                    {
                        let volume =
                            |x: &Body| k().mass_properties(x, 1e-9).expect("a volume").volume;
                        if let Ok(cut) = k().boolean(a, b, BooleanOp::Subtract) {
                            let whole = volume(a);
                            assert!((volume(&cut) - whole).abs() < 1e-6 * whole, "{what}");
                        }
                    }
                    Err(e) => panic!("{what} {op:?}: {e}"),
                }
            }
        };
        // Review's two: the union refused, and the cut with a sliver.
        let axis = Frame::from_normal(DVec3::ZERO, DVec3::ONE.normalize(), DVec3::X);
        let (a, b) = pair(axis, 0.5, 0.1936177480131406, -1e-4, 1e-6, 1.0);
        let union = k().boolean(&a, &b, BooleanOp::Union).expect("a union");
        let volume = |x: &Body| k().mass_properties(x, 1e-9).expect("a volume").volume;
        let both = k()
            .boolean(&a, &b, BooleanOp::Intersect)
            .expect("an intersection");
        let expected = volume(&a) + volume(&b) - volume(&both);
        assert!((volume(&union) - expected).abs() < 1e-9 * expected);
        check(&a, &b, "union");
        let tilted = Frame {
            origin: DVec3::ZERO,
            x: DVec3::new(0.1886813858826233, -0.238360305074152, -0.9526718740397293),
            y: DVec3::new(
                -0.8699696966508363,
                0.40946832331301936,
                -0.27475155888998387,
            ),
            z: DVec3::new(0.4555788203271353, 0.8806361661709003, -0.13010719158122622),
        };
        let (a, b) = pair(tilted, 1.0, 4.679463002356427, -1e-6, 1e-5, 1.0);
        check(&a, &b, "cut");
        // And more of their kind, at unit and ten times the size.
        let mut rng = crate::boolean::testing::Rng(0x0F01_D5ED);
        for i in 0..48 {
            let dir = DVec3::new(
                rng.range(-1.0, 1.0),
                rng.range(-1.0, 1.0),
                rng.range(-1.0, 1.0),
            );
            let axis = Frame::from_normal(DVec3::ZERO, dir.normalize(), DVec3::X);
            let slope = [0.5, 1.0, 2.0][i % 3];
            let phi = rng.range(0.0, TAU);
            let f = [-1e-6, -1e-4, -1e-3, 1e-4][i % 4];
            let s = [1.0, 10.0][i % 2];
            let (a, b) = pair(axis, slope, phi, f, [1e-6, 1e-5, 1e-4][i % 3] * s, s);
            check(&a, &b, &format!("{i}"));
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
    /// midpoints, an S-curve's mesh strayed 3.6% at a chord of 1e-2 and
    /// 5.2% at 1e-3, and a torus's 1.3% at 1e-2, and without that
    /// quadratic a ball's 2.4% at 1e-3 (sampled on a grid of twelfths).
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
                for (k, p) in line.iter().enumerate().take(line.len() - 1).skip(i) {
                    best = best.min(along(r, k));
                    if p.y - r.y > best {
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

    /// The turn between two unit normals is never less than their angle,
    /// and close to it while small (under `0.72 s³` above it for an angle
    /// under half a radian, `s` the sine of half of it): random pairs at
    /// every angle, nearly equal ones, opposite ones.
    #[test]
    fn a_turn_bounds_the_angle_from_above() {
        let mut rng = crate::boolean::testing::Rng(0x7E57_7C2A);
        let angle = |a: DVec3, b: DVec3| a.cross(b).length().atan2(a.dot(b));
        for _ in 0..20_000 {
            let a = rng.unit().normalize();
            // b at an angle up to π from a, from tiny to wide.
            let t = PI * rng.next().powi(1 + rng.below(4) as i32);
            let axis = a.any_orthonormal_vector();
            let b = (a * t.cos() + axis * t.sin()).normalize();
            let (up, exact) = (turn(a, b), angle(a, b));
            // Within rounding: both come from coordinates of unit vectors.
            assert!(up >= exact - 1e-15, "{t}: {up} < {exact}");
            assert!(up <= PI);
            if exact < 0.5 {
                let s = (0.5 * exact).sin();
                assert!(
                    up - exact <= 2.0 * s.powi(3) * 0.36 + 1e-15,
                    "{t}: {up} {exact}"
                );
            }
        }
        assert_eq!(turn(DVec3::Z, -DVec3::Z), PI);
        assert_eq!(turn(DVec3::X, DVec3::X), 0.0);
    }

    /// The quartic estimate between the points tested, bounded by its
    /// Bernstein-Bézier coefficients: their conversion and halving against
    /// its values; known peaks; random quartics against a fine grid of
    /// their values; the two that review found a grid search let pass (one
    /// 1.8% past the tolerance and 12% above the grid's largest, for which
    /// a gate skipped the search, and one past it on an edge between the
    /// grid's points); and the bound that decides when it is made.
    #[test]
    fn the_quartic_estimate_is_bounded() {
        let at = |g: &dyn Fn([f64; 3]) -> f64| QUARTERS.map(|q| g(q.map(|x| x as f64 / 4.0)));
        let value = |at: &[f64; 12], l: [f64; 3]| -> f64 {
            quartic_basis(l).iter().zip(at).map(|(b, y)| b * y).sum()
        };
        let dot =
            |a: &[f64; 15], b: &[f64; 15]| -> f64 { a.iter().zip(b).map(|(x, y)| x * y).sum() };
        let mut rng = crate::boolean::testing::Rng(0x0BEA_C0DE);
        let bz = &*BEZIER;
        for _ in 0..20 {
            let d: [f64; 12] = std::array::from_fn(|_| rng.range(-1.0, 1.0));
            let c = bz.coefficients(&d);
            for _ in 0..10 {
                let (x, y) = (rng.range(0.0, 1.0), rng.range(0.0, 1.0));
                let l = if x + y <= 1.0 {
                    [x, y, 1.0 - x - y]
                } else {
                    [1.0 - x, 1.0 - y, x + y - 1.0]
                };
                let b = dot(&bernstein(l), &c);
                assert!((b - value(&d, l)).abs() < 1e-12, "{b} {}", value(&d, l));
                for (h, corners) in HALVES.iter().enumerate() {
                    let ch = bz.half(h, &c);
                    let p = std::array::from_fn(|a| (0..3).map(|v| l[v] * corners[v][a]).sum());
                    assert!((dot(&bernstein(l), &ch) - dot(&bernstein(p), &c)).abs() < 1e-12);
                }
            }
        }
        let none = [false; 3];
        // A cubic error along an edge peaks at a third of it, 4/27, above
        // all the points tested (9/64 at three quarters): found there.
        let edge = at(&|l| l[1] * l[1] * l[0]);
        assert!(edge.iter().all(|&x| x <= 9.0 / 64.0));
        let l = past(&edge, none, 0.148).expect("past");
        assert!(l[2] == 0.0 && value(&edge, l) > 0.148, "{l:?}");
        assert_eq!(past(&edge, none, 4.0 / 27.0 * 1.001), None);
        // Left out, the edge counts by its values alone (its coefficient
        // off it is 1/12): past a value under them, a point off it is.
        assert_eq!(past(&edge, [true, false, false], 0.084), None);
        let l = past(&edge, [true, false, false], 0.05).expect("past");
        assert!(l[2] > 0.0 && value(&edge, l) > 0.05, "{l:?}");
        // The bubble peaks at the centroid; a quadratic's equal midpoints
        // at four thirds of them; nothing anywhere is zero.
        let bubble = at(&|l| 27.0 * l[0] * l[1] * l[2]);
        assert!(past(&bubble, none, 0.999).is_some());
        assert_eq!(past(&bubble, none, 1.001), None);
        let quadratic = at(&|l| 4.0 * (l[0] * l[1] + l[1] * l[2] + l[2] * l[0]));
        assert!(past(&quadratic, none, 4.0 / 3.0 - 1e-3).is_some());
        assert_eq!(past(&quadratic, none, 4.0 / 3.0 + 1e-3), None);
        assert_eq!(past(&[0.0; 12], none, 0.0), None);
        // Review's two, peaking at 1.01795 and at 1.02397 on edge 12.
        let gate = [
            -0.3611585035471566,
            0.9094083989247623,
            0.37904350143913534,
            -1.202217860587422,
            2.0429570521603937,
            -1.1339729795303184,
            -0.716947599259962,
            -0.6800351169042179,
            -0.1892340169872107,
            1.5752034001628834,
            -0.23594576132925216,
            0.22848214460088379,
        ]
        .map(|x| x * 0.37);
        let search = [
            -0.831191981437047,
            0.6169252596250212,
            0.9885251333825471,
            0.8316466780302958,
            -1.2713504551309718,
            0.7339366431706363,
            -1.128768640394832,
            -1.2368455436308736,
            -0.5254951362597207,
            -0.5814529199266456,
            0.5680130903624987,
            -0.397719588040232,
        ]
        .map(|x| x * 0.58);
        for at in [gate, search] {
            assert!(at.iter().all(|x: &f64| x.abs() <= 1.0));
            let l = past(&at, none, 1.0).expect("past");
            assert!(value(&at, l).abs() > 1.0, "{l:?}");
        }
        // Random quartics, against the largest on a grid of 300ths: past
        // a thousandth under it, proven under a hundredth over it.
        let fine = |at: &[f64; 12]| {
            let n = 300;
            let mut best = 0.0f64;
            for i in 0..=n {
                for j in 0..=n - i {
                    let l = [i, j, n - i - j].map(|x| x as f64 / n as f64);
                    best = best.max(value(at, l).abs());
                }
            }
            best
        };
        for k in 0..100 {
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
            let best = fine(&at);
            let l = past(&at, none, best * 0.999).expect("past");
            assert!(value(&at, l).abs() > best * 0.999, "{k}: {l:?}");
            assert_eq!(past(&at, none, best * 1.01), None, "{k}");
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
