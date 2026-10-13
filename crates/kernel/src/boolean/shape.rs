//! The operands of a boolean as its stages see them (design note 0012
//! §10.1): each face on a surface recognised exactly (a plane, cylinder,
//! cone or sphere), faces on one surface sharing it, and boxes holding
//! every edge and face.

use std::f64::consts::{PI, TAU};

use glam::DVec3;

use super::quadric::{self, Quadric};
use crate::bvh::Bvh;
use crate::geom::{Curve3, Surface};
use crate::math::Frame;
use crate::topo::Brep;
use crate::KernelError;

/// A face of an operand.
pub(super) struct Face {
    /// For a face on a plane, the plane's frame: origin on the plane,
    /// `x × y` the face's outward normal. Unused for a curved face.
    pub frame: Frame,
    /// The face's surface in its operand (faces on one surface share it).
    pub surface: usize,
    /// True if the face's outward normal points the way the surface's
    /// signed distance grows (for a plane, along the normal it was recorded
    /// with: that of the first face on it).
    pub outward: bool,
    /// True for a face on a plane whose edges are all straight: a face of
    /// step 1.
    pub flat: bool,
    /// Per loop (outer first), the coedges as (edge, forward).
    pub loops: Vec<Vec<(usize, bool)>>,
    /// Per loop, the vertex each coedge starts at.
    pub verts: Vec<Vec<usize>>,
    /// For a face on a sphere, the unit direction from the centre to the
    /// middle of the half of the sphere that holds the face (longitude π/2
    /// of its surface's frame, on the equator; ADR 0008); zero otherwise.
    pub middle: DVec3,
    /// A box holding the face: of its edges' boxes (a face on a plane, a
    /// cylinder or a cone lies in the convex hull of its boundary), and for
    /// a face on a sphere also the sphere's extreme points along the axes
    /// that can be in it.
    pub lo: DVec3,
    pub hi: DVec3,
    /// For a face with many edges, its coedges (in loop order, all loops
    /// in turn) indexed by their boxes.
    index: Option<Bvh>,
    /// The position in loop order of each loop's first coedge.
    starts: Vec<usize>,
}

/// Faces with more coedges than this get an index of them.
const INDEXED_FACE: usize = 32;

/// Relative tolerance for recognising faces built on one plane (rounding
/// only), as the tessellator's seam test does: the sine of the angle
/// between their normals, and their offset relative to their distance from
/// the origin.
pub(super) const SAME_PLANE: f64 = 1e-10;

/// True if unit normals `n` and `m` are parallel up to rounding, facing
/// either way. The test is on the sine of the angle between them: one on
/// the cosine (`|n · m| ≥ 1 − SAME_PLANE`) would pass angles up to
/// `√(2 · SAME_PLANE)`, 1.4e-5 rad, and so take two faces of a body that
/// meet at a slight crease for one plane, putting the far corners of the
/// second microns off the first's plane.
pub(super) fn parallel(n: DVec3, m: DVec3) -> bool {
    n.cross(m).length() <= SAME_PLANE
}

impl Face {
    /// The loop and position in it of the `k`-th coedge in loop order.
    pub fn coedge_at(&self, k: usize) -> (usize, usize) {
        let li = self.starts.partition_point(|&s| s <= k) - 1;
        (li, k - self.starts[li])
    }

    /// The number of coedges in all loops.
    pub fn coedges(&self) -> usize {
        self.loops.iter().map(Vec::len).sum()
    }

    /// The coedges (by position in loop order, increasing) whose edges are
    /// not wholly farther than `eps` on one side of `surface`: the only ones
    /// that can meet it or change side across it. Also whether some coedge
    /// left out lies above the surface (on the side its distance grows to),
    /// and whether some lies below it.
    pub fn near_surface(&self, surface: &Quadric, eps: f64, out: &mut Vec<usize>) -> (bool, bool) {
        out.clear();
        let Some(index) = &self.index else {
            out.extend(0..self.coedges());
            return (false, false);
        };
        let (mut above, mut below) = (false, false);
        // A margin far above the rounding of the vertices' signed
        // distances, so no vertex within `eps` is left out.
        let scale = surface.magnitude() + self.lo.abs().max(self.hi.abs()).max_element();
        let reach = eps + 1e-9 * (1.0 + scale);
        let mut side = |d: f64, r: f64| {
            if d - r > reach {
                above = true;
                false
            } else if d + r < -reach {
                below = true;
                false
            } else {
                true
            }
        };
        match *surface {
            Quadric::Plane { point, normal } => index.visit(
                |lo, hi| {
                    let centre = (lo + hi) * 0.5;
                    let half = (hi - lo) * 0.5;
                    side((centre - point).dot(normal), half.dot(normal.abs()))
                },
                |k| {
                    out.push(k);
                    true
                },
            ),
            // The signed distance changes by at most the distance moved,
            // so a box is on one side if the distance at its centre is more
            // than half its diagonal.
            _ => index.visit(
                |lo, hi| side(surface.distance((lo + hi) * 0.5), (hi - lo).length() * 0.5),
                |k| {
                    out.push(k);
                    true
                },
            ),
        };
        out.sort_unstable();
        (above, below)
    }

    /// The face's outward normal, for a face on a plane.
    pub fn normal(&self) -> DVec3 {
        self.frame.z
    }

    pub fn to2(&self, p: DVec3) -> glam::DVec2 {
        let d = p - self.frame.origin;
        glam::DVec2::new(d.dot(self.frame.x), d.dot(self.frame.y))
    }
}

/// An operand.
pub(super) struct Shape<'a> {
    pub brep: &'a Brep,
    pub faces: Vec<Face>,
    pub surfaces: Vec<Quadric>,
    /// Per edge, its start and end vertex.
    pub edges: Vec<(usize, usize)>,
    /// Per edge, true if it is straight.
    pub straight: Vec<bool>,
    /// Per edge, a box holding it.
    pub edge_boxes: Vec<(DVec3, DVec3)>,
    /// Vertex positions.
    pub points: Vec<DVec3>,
}

impl Shape<'_> {
    /// The ends of edge `e` (straight or not).
    pub fn edge_points(&self, e: usize) -> (DVec3, DVec3) {
        let (a, b) = self.edges[e];
        (self.points[a], self.points[b])
    }

    /// The curve of edge `e` and its parameters at its start and end.
    pub fn curve(&self, e: usize) -> (&Curve3, f64, f64) {
        let edge = &self.brep.edges[e];
        (&edge.curve, edge.t0, edge.t1)
    }

    /// The box of the whole operand: of its faces' boxes.
    pub fn extent(&self) -> (DVec3, DVec3) {
        self.faces.iter().fold(
            (DVec3::splat(f64::MAX), DVec3::splat(f64::MIN)),
            |(lo, hi), f| (lo.min(f.lo), hi.max(f.hi)),
        )
    }
}

/// How far past the half of its sphere that holds it (in the cosine of
/// the angle from that half's middle) a face's boundary may reach before
/// the face is refused: a face lies in its half by construction (ADR
/// 0008), and its chart serves it (at a scale of at most 4) anywhere
/// within 120° of the middle.
const HALF: f64 = 0.5;

/// The middle of the half of the sphere (centre `center`, radius `radius`)
/// that holds the face with loops `loops` on the surface `surface`: the
/// direction of longitude π/2 on the equator of the surface's frame. A
/// face's points are within longitudes `[0, π]` of that frame (ADR 0008),
/// so within a quarter turn of it; this is checked on the boundary, where
/// the least cosine of the angle from it is found in closed form along each
/// edge (the face does not hold the antipode, so its least is on its
/// boundary).
fn sphere_middle(
    b: &Brep,
    loops: &[Vec<(usize, bool)>],
    surface: &Surface,
    center: DVec3,
    radius: f64,
    which: &str,
) -> Result<DVec3, KernelError> {
    let refuse = || {
        unsupported(
            which,
            "faces on spheres outside the half of the sphere their surface's frame gives",
        )
    };
    let Surface::Revolution { frame, .. } = surface else {
        return Err(refuse());
    };
    let middle = frame.y.normalize();
    let least = loops
        .iter()
        .flatten()
        .map(|&(e, _)| {
            let edge = &b.edges[e];
            let along = |p: DVec3| (p - center).dot(middle);
            match &edge.curve {
                Curve3::Line { .. } => {
                    along(edge.curve.point(edge.t0)).min(along(edge.curve.point(edge.t1)))
                }
                Curve3::Circle {
                    center: o,
                    x,
                    y,
                    radius: r,
                } => along(*o) + r * least_of_wave(x.dot(middle), y.dot(middle), edge.t0, edge.t1),
                // Within the hull of the edge's control points.
                Curve3::Nurbs(_) => quadric::hull_points(&edge.curve, edge.t0, edge.t1).map_or(
                    f64::NEG_INFINITY,
                    |pieces| {
                        pieces
                            .iter()
                            .flatten()
                            .map(|&p| along(p))
                            .fold(f64::INFINITY, f64::min)
                    },
                ),
            }
        })
        .fold(f64::INFINITY, f64::min);
    if least >= -HALF * radius {
        Ok(middle)
    } else {
        Err(refuse())
    }
}

/// The least of `a cos t + b sin t` for `t` between `t0` and `t1`: at an
/// end, or at the wave's trough `atan2(b, a) + π` (mod 2π) if that lies
/// between them.
pub(super) fn least_of_wave(a: f64, b: f64, t0: f64, t1: f64) -> f64 {
    let (lo, hi) = if t0 <= t1 { (t0, t1) } else { (t1, t0) };
    let wave = |t: f64| a * t.cos() + b * t.sin();
    let low = b.atan2(a) + PI;
    let k = ((lo - low) / TAU).ceil();
    if low + k * TAU <= hi {
        -a.hypot(b)
    } else {
        wave(lo).min(wave(hi))
    }
}

fn unsupported(which: &str, what: &str) -> KernelError {
    KernelError::Unsupported(format!(
        "booleans of bodies with {what} (the {which} body has some; design note 0012 §10)"
    ))
}

/// The operand `b` (`which` names it in messages).
pub(super) fn shape<'a>(b: &'a Brep, which: &str) -> Result<Shape<'a>, KernelError> {
    let points: Vec<DVec3> = b.vertices.iter().map(|v| v.point).collect();
    let mut edges = Vec::with_capacity(b.edges.len());
    let mut straight = Vec::with_capacity(b.edges.len());
    let mut edge_boxes = Vec::with_capacity(b.edges.len());
    for e in &b.edges {
        let (v0, v1) = (e.v0.idx(), e.v1.idx());
        if v0 == v1 {
            return Err(unsupported(which, "closed edges"));
        }
        let (p, q) = (points[v0], points[v1]);
        let line = matches!(e.curve, Curve3::Line { .. });
        let (lo, hi) = if line {
            (p.min(q), p.max(q))
        } else {
            let (lo, hi) = quadric::curve_box(&e.curve, e.t0, e.t1)
                .ok_or_else(|| unsupported(which, "edges of degree above five"))?;
            (lo.min(p).min(q), hi.max(p).max(q))
        };
        edges.push((v0, v1));
        straight.push(line);
        edge_boxes.push((lo, hi));
    }
    let mut surfaces: Vec<Quadric> = Vec::new();
    let mut faces = Vec::with_capacity(b.faces.len());
    for f in &b.faces {
        let mut loops = Vec::with_capacity(f.loops.len());
        let mut verts = Vec::with_capacity(f.loops.len());
        let (mut lo, mut hi) = (DVec3::splat(f64::MAX), DVec3::splat(f64::MIN));
        let mut middle = DVec3::ZERO;
        let mut all_straight = true;
        for l in &f.loops {
            let mut cs = Vec::with_capacity(l.coedges.len());
            let mut vs = Vec::with_capacity(l.coedges.len());
            for c in &l.coedges {
                let e = c.edge.idx();
                cs.push((e, c.forward));
                vs.push(b.coedge_start(*c).idx());
                lo = lo.min(edge_boxes[e].0);
                hi = hi.max(edge_boxes[e].1);
                all_straight &= straight[e];
            }
            loops.push(cs);
            verts.push(vs);
        }
        let (surface, outward, frame, planar) = if let Some((p, n)) = f.surface.plane() {
            let n = if f.forward { n } else { -n };
            let frame = Frame::from_normal(p, n, n.any_orthonormal_vector());
            let scale = 1.0 + p.length();
            let found = surfaces.iter().position(|q| match *q {
                Quadric::Plane {
                    point: qp,
                    normal: qn,
                } => parallel(qn, n) && qn.dot(p - qp).abs() <= SAME_PLANE * scale,
                _ => false,
            });
            let i = found.unwrap_or_else(|| {
                surfaces.push(Quadric::Plane {
                    point: p,
                    normal: n,
                });
                surfaces.len() - 1
            });
            let outward = match surfaces[i] {
                Quadric::Plane { normal, .. } => normal.dot(n) > 0.0,
                _ => unreachable!("a plane group holds a plane"),
            };
            (i, outward, frame, true)
        } else {
            let q = Quadric::of(&f.surface).ok_or_else(|| {
                unsupported(
                    which,
                    "faces other than planes, cylinders, cones and spheres",
                )
            })?;
            let i = surfaces.iter().position(|s| s.same(&q)).unwrap_or_else(|| {
                surfaces.push(q);
                surfaces.len() - 1
            });
            if loops.iter().all(Vec::is_empty) {
                return Err(KernelError::InvalidInput(format!(
                    "a face of the {which} body has no edges"
                )));
            }
            // Which way the face faces, at the most regular of its corners
            // and the middles of its edges: where `S_u × S_v` is largest. A
            // cone's apex and a sphere's poles, often corners, have no
            // normal, and the parameterisation's derivatives there are
            // rounding, whose cross product points any way.
            let corners = verts.iter().flatten().map(|&v| points[v]);
            let middles = loops.iter().flatten().map(|&(e, _)| {
                let edge = &b.edges[e];
                edge.curve.point(0.5 * (edge.t0 + edge.t1))
            });
            let (p, n) = corners
                .chain(middles)
                .map(|p| {
                    let (u, v) = f.surface.project(p);
                    let (su, sv) = f.surface.derivs(u, v);
                    (p, su.cross(sv))
                })
                .fold(None, |best: Option<(DVec3, DVec3)>, (p, n)| match best {
                    Some((_, m)) if m.length() >= n.length() => best,
                    _ => Some((p, n)),
                })
                .and_then(|(p, n)| Some((p, n.try_normalize()?)))
                .ok_or_else(|| {
                    KernelError::InvalidInput(format!(
                        "a face of the {which} body has no regular point on its boundary"
                    ))
                })?;
            let n = if f.forward { n } else { -n };
            let g = surfaces[i].gradient(p);
            let along = n.dot(g);
            if along.abs() < 0.5 {
                return Err(KernelError::InvalidInput(format!(
                    "a face of the {which} body is not on the surface it records"
                )));
            }
            if let Quadric::Sphere { center, radius } = q {
                middle = sphere_middle(b, &loops, &f.surface, center, radius, which)?;
                for e in [
                    DVec3::X,
                    DVec3::Y,
                    DVec3::Z,
                    -DVec3::X,
                    -DVec3::Y,
                    -DVec3::Z,
                ] {
                    if e.dot(middle) >= -HALF {
                        lo = lo.min(center + e * radius);
                        hi = hi.max(center + e * radius);
                    }
                }
            }
            (i, along > 0.0, Frame::WORLD, false)
        };
        let starts: Vec<usize> = loops
            .iter()
            .scan(0, |at, l| {
                let s = *at;
                *at += l.len();
                Some(s)
            })
            .collect();
        let count: usize = loops.iter().map(Vec::len).sum();
        let index = (count > INDEXED_FACE).then(|| {
            Bvh::new(
                loops
                    .iter()
                    .flatten()
                    .map(|&(e, _)| edge_boxes[e])
                    .collect(),
            )
        });
        faces.push(Face {
            frame,
            middle,
            surface,
            outward,
            flat: planar && all_straight,
            loops,
            verts,
            lo,
            hi,
            index,
            starts,
        });
    }
    Ok(Shape {
        brep: b,
        faces,
        surfaces,
        edges,
        straight,
        edge_boxes,
        points,
    })
}
