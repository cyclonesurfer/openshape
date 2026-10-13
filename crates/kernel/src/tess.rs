//! Tolerance-controlled tessellation of B-rep faces.
//!
//! Approach (see `docs/design/0002-tessellation.md`):
//! 1. Every edge is discretized once, so the two faces that share an edge
//!    share exactly the same boundary vertices and the mesh is watertight.
//! 2. Each face's boundary samples are mapped into the face's parameter
//!    space. Singular boundary points (cone apex, sphere pole) are split into
//!    two parameter-space points so the parameter-space boundary stays a
//!    simple polygon.
//! 3. The parameter-space polygon (outer loop plus holes) is triangulated with
//!    a constrained Delaunay triangulation, in coordinates scaled so the
//!    mapping is roughly isometric.
//! 4. Triangles are refined by inserting their parameter-space centroid until
//!    the surface deviates from every triangle by less than the chord
//!    tolerance and the normal varies by less than the angle tolerance.

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
        .map(|e| {
            e.curve
                .discretize(e.t0, e.t1, opts.chord_tolerance, opts.max_angle)
        })
        .collect();
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
    let mut edges = points(&params);
    if match_cusps(brep, &cusps, &mut params, &edges, tol.linear).contains(&true) {
        edges = points(&params);
    }
    let mut loops: Vec<Vec<Vec<Sample>>> = brep.faces.iter().map(|f| loops_uv(f, &edges)).collect();
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
            edges = points(&params);
            // Bisection near a cusp breaks its matching; restore it.
            let rematched = match_cusps(brep, &cusps, &mut params, &edges, tol.linear);
            if rematched.contains(&true) {
                edges = points(&params);
                for (c, r) in changed.iter_mut().zip(rematched) {
                    *c |= r;
                }
            }
        } else {
            break;
        }
        for (face, l) in brep.faces.iter().zip(&mut loops) {
            let touched = face
                .loops
                .iter()
                .flat_map(|lp| &lp.coedges)
                .any(|c| changed[c.edge.idx()]);
            if touched {
                *l = loops_uv(face, &edges);
            }
        }
    }
    let mut faces = Vec::with_capacity(brep.faces.len());
    for (i, (face, l)) in brep.faces.iter().zip(loops).enumerate() {
        faces.push(
            tessellate_face(face, l, opts)
                .map_err(|e| KernelError::Tessellation(format!("face {i}: {e}")))?,
        );
    }
    Ok(BrepTess {
        faces,
        edges,
        seams: seams(brep),
    })
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
    /// Unit direction halfway between the two curves' directions there.
    dir: DVec3,
    /// The two edges leaving the vertex, each with `true` if the vertex is
    /// at the edge's start.
    sides: [(usize, bool); 2],
}

/// The cusps of every face of `brep`.
fn find_cusps(brep: &Brep, max_angle: f64) -> Vec<Cusp> {
    let mut out = Vec::new();
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
                });
            }
        }
    }
    out
}

/// The samples of one side of a cusp over which its curve runs away from
/// the vertex as a graph over the cusp's direction: (sample index, distance
/// along the direction), outward, the vertex itself excluded.
fn cusp_run(
    brep: &Brep,
    c: &Cusp,
    (e, at_start): (usize, bool),
    params: &[f64],
    pts: &[DVec3],
) -> Vec<(usize, f64)> {
    let curve = &brep.edges[e].curve;
    let n = pts.len();
    let sign = if at_start { 1.0 } else { -1.0 };
    let mut out = Vec::new();
    let mut last = 0.0;
    for k in 1..n {
        let j = if at_start { k } else { n - 1 - k };
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
    let run = |s: (usize, bool)| cusp_run(brep, c, s, &params[s.0], &pts[s.0]);
    let (ra, rb) = (run(c.sides[0]), run(c.sides[1]));
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
        let runs = c
            .sides
            .map(|s| cusp_run(brep, c, s, &params[s.0], &pts[s.0]));
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
/// with the edge segment that starts at it.
fn loops_uv(face: &Face, edges: &[Vec<DVec3>]) -> Vec<Vec<Sample>> {
    let s = &face.surface;
    let mut out = Vec::new();
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
                        fixed.push(Sample {
                            uv: a,
                            p: cur.p,
                            seg: None,
                        });
                        fixed.push(Sample { uv: b, ..cur });
                    } else {
                        fixed.push(Sample { uv: a, ..cur });
                    }
                }
            }
        }
        out.push(fixed);
    }
    out
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
) -> Result<FaceTess, String> {
    let s = &face.surface;
    if loops.iter().any(|l| l.len() < 3) {
        return Err("loop has fewer than three boundary samples".into());
    }
    let (mut lo, mut hi) = (DVec2::splat(f64::MAX), DVec2::splat(f64::MIN));
    for l in &loops {
        for smp in l {
            lo = lo.min(smp.uv);
            hi = hi.max(smp.uv);
        }
    }
    let scale = param_scale(s, &loops)?;

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
    let evaluate = |uv: DVec2, p: Option<DVec3>| {
        let sp = p.unwrap_or_else(|| s.point(uv.x, uv.y));
        (sp, s.normal(uv.x, uv.y))
    };
    let placed: Vec<Vec<DVec2>> = loops
        .iter()
        .map(|l| l.iter().map(|smp| place(smp.uv)).collect())
        .collect();
    // A flat face bounded by one strictly convex loop is split by halving
    // the boundary's runs: no point is added on a flat face, every
    // triangulation of a convex polygon on its own vertices is valid (no
    // three of them are on a line), and this one needs no search.
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
        if !is_planar {
            let min_area = (opts.chord_tolerance * 0.01).powi(2);
            // Triangles found fine (or too small to refine) in an earlier pass:
            // the test depends only on their vertices, so it is not repeated.
            let mut settled: FxHashSet<[usize; 3]> = FxHashSet::default();
            for _pass in 0..12 {
                let tris = inside_triangles(&cdt);
                let mut new_pts: Vec<(DVec2, usize)> = Vec::new();
                for t in &tris {
                    let mut key = *t;
                    key.sort_unstable();
                    if settled.contains(&key) {
                        continue;
                    }
                    let uv = t.map(|i| verts[i].0);
                    let q = uv.map(|x| x * scale);
                    let area = 0.5 * (q[1] - q[0]).perp_dot(q[2] - q[0]).abs();
                    let longest = (q[1] - q[0])
                        .length_squared()
                        .max((q[2] - q[1]).length_squared())
                        .max((q[0] - q[2]).length_squared());
                    // Skip tiny triangles and slivers (whose centroid could
                    // land on a boundary constraint within round-off).
                    if area < min_area || 2.0 * area < 1e-8 * longest {
                        settled.insert(key);
                        continue;
                    }
                    let p = t.map(|i| at[i].0);
                    let tri_n = (p[1] - p[0]).cross(p[2] - p[0]);
                    if tri_n.length() <= 0.0 {
                        settled.insert(key);
                        continue;
                    }
                    let tri_n = tri_n.normalize();
                    // Chord error: distance of the surface from the triangle's
                    // plane at the centroid and at interior-edge midpoints.
                    let c = (uv[0] + uv[1] + uv[2]) / 3.0;
                    let off = |w: DVec2| (s.point(w.x, w.y) - p[0]).dot(tri_n).abs();
                    let mut refine = off(c) > opts.chord_tolerance
                        || (0..3).any(|e| {
                            let (a, b) = (t[e], t[(e + 1) % 3]);
                            (verts[a].1.is_none() || verts[b].1.is_none())
                                && off((verts[a].0 + verts[b].0) * 0.5) > opts.chord_tolerance
                        });
                    if !refine {
                        if let Some(nc) = s.normal(c.x, c.y) {
                            refine = t.iter().any(|&i| {
                                at[i]
                                    .1
                                    .map(|nw| nw.angle_between(nc) > opts.max_angle)
                                    .unwrap_or(false)
                            });
                        }
                    }
                    if refine {
                        new_pts.push((c, t[0]));
                    } else {
                        settled.insert(key);
                    }
                }
                if new_pts.is_empty() {
                    break;
                }
                if verts.len() + new_pts.len() > 2_000_000 {
                    return Err("tessellation exceeded the vertex budget".into());
                }
                for (c, hint) in new_pts {
                    insert(&mut cdt, &mut verts, &mut at, c, hint)?;
                }
            }
            if cdt.num_constraints() != boundary_constraints {
                return Err("refinement split a boundary edge".into());
            }
        }

        inside_triangles(&cdt)
    };

    // Emit vertices. Boundary vertices with identical 3D points (split
    // singular points) are welded.
    let mut out = FaceTess::default();
    let mut index_of: Vec<Option<u32>> = vec![None; verts.len()];
    let mut by_pos: FxHashMap<[u64; 3], u32> = FxHashMap::default();
    let center = (lo + hi) * 0.5;
    let sign = if face.forward { 1.0 } else { -1.0 };
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
            let key = [p.x.to_bits(), p.y.to_bits(), p.z.to_bits()];
            if exact.is_some() {
                if let Some(&j) = by_pos.get(&key) {
                    index_of[i] = Some(j);
                    continue;
                }
            }
            let mut n = if is_planar {
                s.normal(uv.x, uv.y)
            } else {
                at[i].1
            };
            let mut k = 1e-6;
            while n.is_none() && k < 0.5 {
                let w = uv + (center - uv) * k;
                n = s.normal(w.x, w.y);
                k *= 4.0;
            }
            let j = out.positions.len() as u32;
            out.positions.push(p);
            out.normals.push(n.unwrap_or(DVec3::Z) * sign);
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
    if out.triangles.is_empty() {
        return Err("face produced no triangles".into());
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

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
