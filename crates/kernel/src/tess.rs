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

use std::collections::HashMap;

use glam::{DVec2, DVec3};
use spade::{ConstrainedDelaunayTriangulation, Point2, Triangulation};

use crate::geom::{Curve3, Surface};
use crate::topo::{Brep, Face};
use crate::KernelError;

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
/// same orientation, as where a periodic face is split (ADR 0008), or where
/// one face uses the edge twice. Such an edge bounds faces but is no crease
/// in the shape, so displays hide it.
///
/// The test ([`same_surface`]) errs towards showing an edge: surfaces it
/// cannot prove identical count as different, so a real crease is never
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
        .map(|u| match u.as_slice() {
            [a, b] => a == b || same_surface(&brep.faces[*a], &brep.faces[*b]),
            _ => false,
        })
        .collect()
}

/// Relative tolerance for comparing surface data built from the same
/// definition by rigid motions (rounding only).
const SAME: f64 = 1e-10;

fn near(a: DVec3, b: DVec3, scale: f64) -> bool {
    (a - b).length() <= SAME * scale
}

fn parallel(a: DVec3, b: DVec3) -> bool {
    let (a, b) = (a.normalize_or_zero(), b.normalize_or_zero());
    a.dot(b) >= 1.0 - SAME
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
                && xa.cross(*ya).dot(xb.cross(*yb)).abs() >= 1.0 - SAME
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

/// Tessellates every face and edge of `brep`.
pub(crate) fn tessellate(brep: &Brep, opts: &TessellationOptions) -> Result<BrepTess, KernelError> {
    if !(opts.chord_tolerance > 0.0 && opts.max_angle > 0.0) {
        return Err(KernelError::InvalidInput(
            "tessellation tolerances must be positive".into(),
        ));
    }
    let edges: Vec<Vec<DVec3>> = brep
        .edges
        .iter()
        .map(|e| {
            let ts = e
                .curve
                .discretize(e.t0, e.t1, opts.chord_tolerance, opts.max_angle);
            let mut pts: Vec<DVec3> = ts.iter().map(|&t| e.curve.point(t)).collect();
            // Pin the ends to the vertices so loops close exactly.
            let n = pts.len();
            pts[0] = brep.vertex(e.v0).point;
            pts[n - 1] = brep.vertex(e.v1).point;
            pts
        })
        .collect();
    let mut faces = Vec::with_capacity(brep.faces.len());
    for (i, face) in brep.faces.iter().enumerate() {
        faces.push(
            tessellate_face(brep, face, &edges, opts)
                .map_err(|e| KernelError::Tessellation(format!("face {i}: {e}")))?,
        );
    }
    Ok(BrepTess {
        faces,
        edges,
        seams: seams(brep),
    })
}

#[derive(Clone, Copy)]
struct Sample {
    uv: DVec2,
    p: DVec3,
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

fn face_loops_uv(
    brep: &Brep,
    face: &Face,
    edges: &[Vec<DVec3>],
) -> Result<Vec<Vec<Sample>>, String> {
    let s = &face.surface;
    let mut out = Vec::new();
    for lp in &face.loops {
        let mut pts: Vec<DVec3> = Vec::new();
        for c in &lp.coedges {
            let ep = &edges[c.edge.idx()];
            if c.forward {
                pts.extend_from_slice(&ep[..ep.len() - 1]);
            } else {
                pts.extend(ep.iter().rev().take(ep.len() - 1));
            }
        }
        let _ = brep;
        if pts.len() < 3 {
            return Err("loop has fewer than three boundary samples".into());
        }
        let raw: Vec<Sample> = pts
            .iter()
            .map(|&p| {
                let (u, v) = s.project(p);
                Sample {
                    uv: DVec2::new(u, v),
                    p,
                }
            })
            .collect();
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
                    fixed.push(Sample { uv: a, p: cur.p });
                    if (a - b).length() > 1e-12 {
                        fixed.push(Sample { uv: b, p: cur.p });
                    }
                }
            }
        }
        out.push(fixed);
    }
    Ok(out)
}

/// Triangles inside the constrained boundary, found by flood fill from the
/// outer face, toggling inside/outside each time a constraint edge is
/// crossed. Purely topological, so robust to nearly collinear boundaries.
fn inside_triangles(cdt: &ConstrainedDelaunayTriangulation<Point2<f64>>) -> Vec<[usize; 3]> {
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
    brep: &Brep,
    face: &Face,
    edges: &[Vec<DVec3>],
    opts: &TessellationOptions,
) -> Result<FaceTess, String> {
    let s = &face.surface;
    let loops = face_loops_uv(brep, face, edges)?;

    // Parameter-space scaling for a roughly isometric triangulation.
    let (mut lo, mut hi) = (DVec2::splat(f64::MAX), DVec2::splat(f64::MIN));
    for l in &loops {
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
    let scale = DVec2::new(su, sv);

    let mut cdt: ConstrainedDelaunayTriangulation<Point2<f64>> =
        ConstrainedDelaunayTriangulation::new();
    // Per CDT vertex: parameter value and, for boundary vertices, the exact 3D point.
    let mut verts: Vec<(DVec2, Option<DVec3>)> = Vec::new();
    let insert = |cdt: &mut ConstrainedDelaunayTriangulation<Point2<f64>>,
                  verts: &mut Vec<(DVec2, Option<DVec3>)>,
                  uv: DVec2,
                  p: Option<DVec3>|
     -> Result<spade::handles::FixedVertexHandle, String> {
        // Spade rejects non-zero coordinates below ~1e-42; such values are
        // round-off (e.g. an angle of 1e-60 rad) and are flushed to zero.
        let flush = |x: f64| if x.abs() < 1e-30 { 0.0 } else { x };
        let q = uv * scale;
        let q = DVec2::new(flush(q.x), flush(q.y));
        let h = cdt
            .insert(Point2::new(q.x, q.y))
            .map_err(|e| format!("triangulation insert failed: {e:?}"))?;
        if h.index() == verts.len() {
            verts.push((uv, p));
        }
        Ok(h)
    };
    for l in &loops {
        let mut handles = Vec::with_capacity(l.len());
        for smp in l {
            handles.push(insert(&mut cdt, &mut verts, smp.uv, Some(smp.p))?);
        }
        for i in 0..handles.len() {
            let a = handles[i];
            let b = handles[(i + 1) % handles.len()];
            if a != b && cdt.try_add_constraint(a, b).is_empty() && !cdt.exists_constraint(a, b) {
                return Err("face boundary intersects itself in parameter space".into());
            }
        }
    }

    let boundary_constraints = cdt.num_constraints();
    let is_planar = matches!(s, Surface::Plane(_));
    if !is_planar {
        let min_area = (opts.chord_tolerance * 0.01).powi(2);
        for _pass in 0..12 {
            let tris = inside_triangles(&cdt);
            let mut new_pts: Vec<DVec2> = Vec::new();
            for t in &tris {
                let uv: Vec<DVec2> = t.iter().map(|&i| verts[i].0).collect();
                let q: Vec<DVec2> = uv.iter().map(|x| *x * scale).collect();
                let area = 0.5 * (q[1] - q[0]).perp_dot(q[2] - q[0]).abs();
                let longest = (q[1] - q[0])
                    .length_squared()
                    .max((q[2] - q[1]).length_squared())
                    .max((q[0] - q[2]).length_squared());
                // Skip tiny triangles and slivers (whose centroid could land
                // on a boundary constraint within round-off).
                if area < min_area || 2.0 * area < 1e-8 * longest {
                    continue;
                }
                let p: Vec<DVec3> = t
                    .iter()
                    .map(|&i| {
                        verts[i]
                            .1
                            .unwrap_or_else(|| s.point(verts[i].0.x, verts[i].0.y))
                    })
                    .collect();
                let tri_n = (p[1] - p[0]).cross(p[2] - p[0]);
                if tri_n.length() <= 0.0 {
                    continue;
                }
                let tri_n = tri_n.normalize();
                // Chord error: distance of the surface from the triangle's
                // plane at the centroid and at interior-edge midpoints.
                let c = (uv[0] + uv[1] + uv[2]) / 3.0;
                let mut probes = vec![c];
                for e in 0..3 {
                    let (a, b) = (t[e], t[(e + 1) % 3]);
                    if verts[a].1.is_none() || verts[b].1.is_none() {
                        probes.push((verts[a].0 + verts[b].0) * 0.5);
                    }
                }
                let mut refine = probes.iter().any(|w| {
                    let sp = s.point(w.x, w.y);
                    (sp - p[0]).dot(tri_n).abs() > opts.chord_tolerance
                });
                if !refine {
                    if let Some(nc) = s.normal(c.x, c.y) {
                        refine = uv.iter().any(|w| {
                            s.normal(w.x, w.y)
                                .map(|nw| nw.angle_between(nc) > opts.max_angle)
                                .unwrap_or(false)
                        });
                    }
                }
                if refine {
                    new_pts.push(c);
                }
            }
            if new_pts.is_empty() {
                break;
            }
            if verts.len() + new_pts.len() > 2_000_000 {
                return Err("tessellation exceeded the vertex budget".into());
            }
            for c in new_pts {
                insert(&mut cdt, &mut verts, c, None)?;
            }
        }
        if cdt.num_constraints() != boundary_constraints {
            return Err("refinement split a boundary edge".into());
        }
    }

    // Emit vertices. Boundary vertices with identical 3D points (split
    // singular points) are welded.
    let tris = inside_triangles(&cdt);
    let mut out = FaceTess::default();
    let mut index_of: Vec<Option<u32>> = vec![None; verts.len()];
    let mut by_pos: HashMap<[u64; 3], u32> = HashMap::new();
    let center = (lo + hi) * 0.5;
    let sign = if face.forward { 1.0 } else { -1.0 };
    for t in &tris {
        for &i in t {
            if index_of[i].is_some() {
                continue;
            }
            let (uv, exact) = verts[i];
            let p = exact.unwrap_or_else(|| s.point(uv.x, uv.y));
            let key = [p.x.to_bits(), p.y.to_bits(), p.z.to_bits()];
            if exact.is_some() {
                if let Some(&j) = by_pos.get(&key) {
                    index_of[i] = Some(j);
                    continue;
                }
            }
            let mut n = s.normal(uv.x, uv.y);
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
