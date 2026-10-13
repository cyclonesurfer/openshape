//! Classifying the pieces of the faces (design note 0012 §4.3).
//!
//! A piece that lies on a face of the other body (the two faces on one
//! plane) is *on* it, with the same or the opposite outward normal. Any
//! other piece is wholly inside or wholly outside the other body, since
//! every place where that body's surface crosses the face is a cut between
//! pieces; it is classified at one point well inside it (for a convex
//! piece, the better of its centroid and the incentres of a few triangles
//! of it; otherwise the incentre of the triangle with the largest
//! inscribed circle in a triangulation of it).
//!
//! A point within the tolerance of the other body's surface is on it, so a
//! piece whose point is there without lying on a coincident face has no
//! side the tolerance can decide, and the boolean is refused. Any other
//! point is decided exactly, by ray parity: the other body's faces are
//! triangulated on its own vertices (so the triangles close up exactly,
//! edge to edge), and a segment from the point to beyond the body counts
//! the triangles it passes through, +1 leaving the body and −1 entering
//! it. Every test is an exact orientation predicate (Shewchuk, "Adaptive
//! precision floating-point arithmetic and fast robust geometric
//! predicates", 1997), so the count is the winding number of the
//! triangulated surface about the point: 1 inside, 0 outside, and anything
//! else an error. A segment that passes exactly through an edge or a
//! corner of a triangle is not counted; the test is repeated along another
//! fixed direction. A bounding-volume tree over the triangles keeps each
//! test to the triangles near the segment.
//!
//! Faces that are not flat (on a cylinder, or on a plane with curved
//! edges) are met in closed form instead (design note 0012 §10.5): the ray
//! meets a face's surface at the roots of a polynomial of degree one or
//! two, and a root counts if it lies in the face, decided in the face's
//! chart against its outline. A ray that grazes such a surface, or meets
//! it within the tolerance of the face's boundary, is abandoned for the
//! next direction, as one through a triangle's edge is.
//!
//! This replaces the propagation of design note 0012 §4.3 (classify faces
//! at the intersection edges, then flood fill): each piece is decided by
//! one exact test, and no decision depends on another.

use glam::{DVec2, DVec3};
use spade::Triangulation;

use super::chart::Outline;
use super::quadric::Quadric;
use super::shape::{Face, Shape};
use super::show;
use crate::bvh::{self, Bvh};
use crate::tess;
use crate::KernelError;

/// Where a piece lies with respect to the other body.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum Class {
    Inside,
    Outside,
    /// On a face of the other body with the same outward normal.
    OnSame,
    /// On a face of the other body with the opposite outward normal.
    OnOpposite,
}

/// Directions of the parity test, tried in turn until one passes through
/// no triangle's edge or corner: fixed, so the result is reproducible, and
/// generic (none along an axis, a face diagonal or a body diagonal, so not
/// along the edges and faces models are usually built from), so that the
/// first one nearly always serves.
const DIRECTIONS: [[f64; 3]; 8] = [
    [0.538, 0.632, 0.557],
    [-0.707, 0.309, 0.636],
    [0.224, -0.866, 0.447],
    [-0.412, -0.520, -0.748],
    [0.831, 0.173, -0.529],
    [-0.141, 0.917, -0.374],
    [0.616, -0.436, -0.656],
    [-0.911, -0.265, 0.316],
];

/// The cosine of the angle between a ray and a curved surface's normal
/// below which the ray is taken to graze the surface and is abandoned.
const GRAZE: f64 = 1e-3;

/// A face of a body that is not flat.
struct Curved<'a> {
    surface: Quadric,
    outward: bool,
    outline: &'a Outline<'a>,
}

impl Curved<'_> {
    /// The signed number of times the ray from `p` along the unit `dir`
    /// passes through the face before `len` (+1 leaving the body, −1
    /// entering it), or `None` if it grazes the surface or meets it within
    /// `eps` of the face's boundary.
    fn hits(&self, p: DVec3, dir: DVec3, len: f64, eps: f64) -> Option<i64> {
        let roots: Vec<f64> = match self.surface {
            Quadric::Plane { point, normal } => {
                let c = dir.dot(normal);
                if c.abs() < GRAZE {
                    return None;
                }
                vec![-(p - point).dot(normal) / c]
            }
            Quadric::Cylinder {
                point,
                axis,
                radius,
            } => {
                // |w + s u|² = r² across the axis.
                let w = p - point;
                let w = w - axis * w.dot(axis);
                let u = dir - axis * dir.dot(axis);
                let (a, b, c) = (u.dot(u), w.dot(u), w.dot(w) - radius * radius);
                if a < GRAZE * GRAZE {
                    return None;
                }
                let disc = b * b - a * c;
                if disc < 0.0 {
                    return Some(0);
                }
                // The cosine of the ray with the surface's normal at either
                // root is √disc / r.
                let sq = disc.sqrt();
                if sq < GRAZE * radius {
                    return None;
                }
                let q = -(b + sq.copysign(b));
                vec![q / a, c / q]
            }
            _ => return None,
        };
        let mut total = 0;
        for s in roots {
            if !(s > 0.0 && s < len) {
                continue;
            }
            let h = p + dir * s;
            match self
                .outline
                .locate(self.outline.chart.map(h), eps, eps / 16.0)
            {
                None => return None,
                Some(false) => {}
                Some(true) => {
                    let g = self.surface.gradient(h);
                    let n = if self.outward { g } else { -g };
                    total += if n.dot(dir) > 0.0 { 1 } else { -1 };
                }
            }
        }
        Some(total)
    }

    /// True if `p` is within `eps` of the face: of a point of its surface
    /// within twice `eps` of the face (`p`'s foot on the surface is within
    /// twice `eps` of the face if `p` is within `eps` of it).
    fn near(&self, p: DVec3, eps: f64) -> bool {
        if self.surface.distance(p).abs() > eps {
            return false;
        }
        match self.surface.foot(p) {
            None => true,
            Some(f) => {
                let q = self.outline.chart.map(f);
                self.outline.locate(q, 2.0 * eps, eps / 16.0) != Some(false)
            }
        }
    }
}

/// A body's boundary: its flat faces as triangles on its own vertices,
/// outward facing, and its other faces on their surfaces.
pub(super) struct Boundary<'a> {
    points: &'a [DVec3],
    tris: Vec<[u32; 3]>,
    tree: Bvh,
    curved: Vec<Curved<'a>>,
    curved_tree: Bvh,
    /// The centre and half the diagonal of the body's box.
    centre: DVec3,
    reach: f64,
    /// Padding of the boxes in the segment test, against rounding.
    pad: f64,
}

impl<'a> Boundary<'a> {
    /// The boundary of `shape`, whose faces that are not flat have the
    /// outlines `outlines`.
    pub fn new(
        shape: &'a Shape,
        outlines: &'a [Option<Outline<'a>>],
    ) -> Result<Boundary<'a>, KernelError> {
        let points = &shape.points[..];
        let mut tris: Vec<[u32; 3]> = Vec::new();
        let mut curved = Vec::new();
        let mut curved_boxes = Vec::new();
        for (fi, f) in shape.faces.iter().enumerate() {
            if let Some(outline) = &outlines[fi] {
                curved.push(Curved {
                    surface: shape.surfaces[f.surface],
                    outward: f.outward,
                    outline,
                });
                curved_boxes.push((f.lo, f.hi));
                continue;
            }
            triangulate(shape, f, &mut tris).map_err(|m| {
                KernelError::InvalidResult(format!("boolean: face {fi} of an operand: {m}"))
            })?;
        }
        // A triangle on three collinear points (a face with collinear
        // vertices, triangulated in its rounded plane) covers no area: a
        // segment through it passes through its neighbours' shared edges,
        // which the test sees, so it is dropped rather than tested.
        tris.retain(|t| !collinear(t.map(|i| points[i as usize])));
        let boxes: Vec<(DVec3, DVec3)> = tris
            .iter()
            .map(|t| {
                let [a, b, c] = t.map(|i| points[i as usize]);
                (a.min(b).min(c), a.max(b).max(c))
            })
            .collect();
        let tree = Bvh::new(boxes);
        let curved_tree = Bvh::new(curved_boxes);
        let (lo, hi) = match (tree.bounds(), curved_tree.bounds()) {
            (Some(a), Some(b)) => (a.0.min(b.0), a.1.max(b.1)),
            (Some(a), None) | (None, Some(a)) => a,
            (None, None) => (DVec3::ZERO, DVec3::ZERO),
        };
        let centre = (lo + hi) * 0.5;
        let reach = (hi - lo).length() * 0.5;
        let pad = 1e-9 * (reach + centre.abs().max_element()) + f64::MIN_POSITIVE;
        Ok(Boundary {
            points,
            tris,
            tree,
            curved,
            curved_tree,
            centre,
            reach,
            pad,
        })
    }

    fn triangle(&self, t: usize) -> [DVec3; 3] {
        self.tris[t].map(|i| self.points[i as usize])
    }

    /// True if `p` is within `eps` of the boundary.
    fn near(&self, p: DVec3, eps: f64) -> bool {
        let r = eps * (1.0 + 1e-9) + self.pad;
        !self.tree.visit(
            |lo, hi| bvh::near_box(p, lo, hi, r),
            |t| (closest_on_triangle(p, &self.triangle(t)) - p).length() > eps,
        ) || !self.curved_tree.visit(
            |lo, hi| bvh::near_box(p, lo, hi, r),
            |i| !self.curved[i].near(p, eps),
        )
    }

    /// The signed number of times the segment `p`–`q` passes through the
    /// boundary (+1 leaving the body, −1 entering it), or `None` if it
    /// passes through an edge or a corner of a triangle, lies in a
    /// triangle's plane across it, grazes a curved face's surface or meets
    /// it within `eps` of the face's boundary. `p` is off the boundary; `q`
    /// is outside the body's box.
    fn crossings(&self, p: DVec3, q: DVec3, eps: f64) -> Option<i64> {
        let mut total = 0i64;
        let mut clean = true;
        let pad = self.pad + 1e-9 * (p - self.centre).length();
        self.tree.visit(
            |lo, hi| bvh::segment_meets_box(p, q, lo, hi, pad),
            |t| {
                let [a, b, c] = self.triangle(t);
                let (sp, sq) = (orient(a, b, c, p), orient(a, b, c, q));
                if sp == 0.0 && sq == 0.0 {
                    // In the triangle's plane: whether it runs across the
                    // triangle is not worth deciding; try another
                    // direction.
                    clean = false;
                    return false;
                }
                // An end on the plane is off the triangle (p is off the
                // boundary, q outside the box), so only a strict change of
                // side can cross it.
                if sp == 0.0 || sq == 0.0 || (sp > 0.0) == (sq > 0.0) {
                    return true;
                }
                let e = [orient(p, q, a, b), orient(p, q, b, c), orient(p, q, c, a)];
                if e.iter().all(|&x| x > 0.0) || e.iter().all(|&x| x < 0.0) {
                    // `p` behind the outward-facing triangle: leaving.
                    total += if sp > 0.0 { 1 } else { -1 };
                    true
                } else if e.iter().all(|&x| x >= 0.0) || e.iter().all(|&x| x <= 0.0) {
                    // Through an edge or a corner.
                    clean = false;
                    false
                } else {
                    true
                }
            },
        );
        if !clean {
            return None;
        }
        let (dir, len) = ((q - p).normalize(), (q - p).length());
        self.curved_tree.visit(
            |lo, hi| bvh::segment_meets_box(p, q, lo, hi, pad),
            |i| match self.curved[i].hits(p, dir, len, eps) {
                None => {
                    clean = false;
                    false
                }
                Some(n) => {
                    total += n;
                    true
                }
            },
        );
        clean.then_some(total)
    }

    /// Inside or outside, at a point farther than `eps` from the boundary;
    /// a point within it is on the boundary, and a piece of a face that is
    /// on it without coinciding with a face is refused.
    pub fn classify(&self, p: DVec3, eps: f64) -> Result<Class, KernelError> {
        if self.near(p, eps) {
            return Err(KernelError::InvalidInput(format!(
                "a face of one body runs within the tolerance of the other's surface \
                 without lying on one of its faces, near {}",
                show(p)
            )));
        }
        let length = 2.0 * ((p - self.centre).length() + self.reach) + 1.0;
        for d in DIRECTIONS {
            let q = p + DVec3::from_array(d).normalize() * length;
            match self.crossings(p, q, eps) {
                None => continue,
                Some(1) => return Ok(Class::Inside),
                Some(0) => return Ok(Class::Outside),
                Some(w) => {
                    return Err(KernelError::InvalidResult(format!(
                        "boolean: the other body's surface winds {w} times around a point \
                         near {}",
                        show(p)
                    )))
                }
            }
        }
        Err(KernelError::InvalidResult(format!(
            "boolean: every test ray from a point near {} passes through an edge or a \
             corner of the other body",
            show(p)
        )))
    }
}

/// Appends the triangles of face `f` on its operand's vertices, facing
/// along its outward normal.
fn triangulate(shape: &Shape, f: &Face, out: &mut Vec<[u32; 3]>) -> Result<(), String> {
    let q2 = |v: usize| f.to2(shape.points[v]);
    // A strictly convex face without holes, counter-clockwise: triangles
    // that halve the boundary's runs level by level. Each triangle spans a
    // run of the boundary and its box is about the size of that run, so a
    // ray meets a few boxes per level; the slivers of a fan all reach
    // across the face, and a ray near it meets most of them.
    if let [vs] = &f.verts[..] {
        let ring: Vec<DVec2> = vs.iter().map(|&v| q2(v)).collect();
        if tess::convexity(&ring) == Some(true) {
            out.extend(
                tess::balanced(vs.len())
                    .into_iter()
                    .map(|t| t.map(|k| vs[k] as u32)),
            );
            return Ok(());
        }
    }
    let loops: Vec<Vec<DVec2>> = f
        .verts
        .iter()
        .map(|l| l.iter().map(|&v| q2(v)).collect())
        .collect();
    let (cdt, index) = tess::polygon_cdt(&loops)?;
    if cdt.num_vertices() != index.len() {
        return Err("two vertices at one point of the face's plane".into());
    }
    let mut at = vec![0; index.len()];
    for (&i, &v) in index.iter().zip(f.verts.iter().flatten()) {
        at[i] = v;
    }
    for t in tess::inside_triangles(&cdt) {
        out.push(t.map(|k| at[k] as u32));
    }
    Ok(())
}

fn orient(a: DVec3, b: DVec3, c: DVec3, d: DVec3) -> f64 {
    let k = |p: DVec3| robust::Coord3D {
        x: p.x,
        y: p.y,
        z: p.z,
    };
    robust::orient3d(k(a), k(b), k(c), k(d))
}

/// True if the three points are exactly on one line: collinear in every
/// coordinate plane.
fn collinear([a, b, c]: [DVec3; 3]) -> bool {
    let xy = |p: DVec3| DVec2::new(p.x, p.y);
    let yz = |p: DVec3| DVec2::new(p.y, p.z);
    let zx = |p: DVec3| DVec2::new(p.z, p.x);
    orient2(xy(a), xy(b), xy(c)) == 0.0
        && orient2(yz(a), yz(b), yz(c)) == 0.0
        && orient2(zx(a), zx(b), zx(c)) == 0.0
}

fn orient2(a: DVec2, b: DVec2, c: DVec2) -> f64 {
    let k = |p: DVec2| robust::Coord { x: p.x, y: p.y };
    robust::orient2d(k(a), k(b), k(c))
}

/// A point well inside the polygon `loops` (outer first) and its distance
/// from the boundary at least. For a convex polygon without holes, the
/// better of its centroid and the incentres of a balanced triangulation of
/// it ([`convex_interior`]); otherwise the incentre of the triangle with
/// the largest inscribed circle in a constrained triangulation of the
/// polygon.
pub(super) fn interior_point(loops: &[Vec<DVec2>]) -> Result<(DVec2, f64), String> {
    if let [l] = loops {
        if let Some(found) = convex_interior(l) {
            return Ok(found);
        }
    }
    let (cdt, index) = tess::polygon_cdt(loops).map_err(|e| format!("a piece's boundary: {e}"))?;
    let mut at = vec![DVec2::ZERO; cdt.num_vertices()];
    for (&i, q) in index.iter().zip(loops.iter().flatten()) {
        at[i] = *q;
    }
    let mut best: Option<(DVec2, f64)> = None;
    for t in tess::inside_triangles(&cdt) {
        let found = incircle(t.map(|i| at[i]));
        if found.is_some_and(|(_, r)| best.is_none_or(|(_, br)| r > br)) {
            best = found;
        }
    }
    best.ok_or_else(|| "a piece with no area".into())
}

/// Sides probed inward by [`probe_point`]: this many of the longest, and
/// this many more spread evenly along the boundary.
const PROBES: usize = 8;

/// A point inside the polygon `loops` (outer first, counter-clockwise;
/// holes clockwise: the inside on the left of every side) and its least
/// distance from the boundary, found without a triangulation, or `None`.
///
/// From the middle of each of a few sides ([`PROBES`] of the longest and
/// as many spread along the boundary), a ray goes inward, square to the
/// side, to the first side it meets; the candidate is halfway along it.
/// Each candidate is checked: its winding number, counted with exact
/// orientation tests, must be one (inside), and its distance from every
/// side is measured; the farthest from the boundary is returned. The work
/// is a few passes over the sides, where a triangulation of a polygon
/// whose points lie on circles (round holes) spends most of its time in
/// exact incircle tests. The point is good, not the best:
/// [`interior_point`] serves when it is not good enough.
pub(super) fn probe_point(loops: &[Vec<DVec2>]) -> Option<(DVec2, f64)> {
    let sides: Vec<(DVec2, DVec2)> = loops
        .iter()
        .flat_map(|l| (0..l.len()).map(move |i| (l[i], l[(i + 1) % l.len()])))
        .collect();
    let m = sides.len();
    if m < 3 {
        return None;
    }
    let mut order: Vec<usize> = (0..m).collect();
    order.sort_by(|&i, &j| {
        let len = |k: usize| (sides[k].1 - sides[k].0).length_squared();
        len(j).total_cmp(&len(i)).then(i.cmp(&j))
    });
    let mut chosen: Vec<usize> = order[..PROBES.min(m)].to_vec();
    chosen.extend((0..PROBES).map(|k| k * m / PROBES));
    chosen.sort_unstable();
    chosen.dedup();
    let mut best: Option<(DVec2, f64)> = None;
    for i in chosen {
        let (a, b) = sides[i];
        let d = b - a;
        if d.length_squared() == 0.0 {
            continue;
        }
        let (start, n) = (0.5 * (a + b), d.perp().normalize());
        // The nearest side the inward ray meets.
        let mut reach = f64::INFINITY;
        for (j, &(c, e)) in sides.iter().enumerate() {
            if j == i {
                continue;
            }
            let f = e - c;
            let den = n.perp_dot(f);
            if den == 0.0 {
                continue;
            }
            let w = c - start;
            let t = w.perp_dot(f) / den;
            let u = w.perp_dot(n) / den;
            if t > 0.0 && (0.0..=1.0).contains(&u) {
                reach = reach.min(t);
            }
        }
        if !reach.is_finite() {
            continue;
        }
        let q = start + n * (0.5 * reach);
        let clearance = sides
            .iter()
            .map(|&(c, e)| super::planar::point_segment(q, c, e))
            .fold(f64::INFINITY, f64::min);
        if clearance > 0.0 && best.is_none_or(|(_, r)| clearance > r) && winding(&sides, q) == 1 {
            best = Some((q, clearance));
        }
    }
    best
}

/// The winding number of the closed sides about `q` (off them), by the
/// sides crossing the ray from `q` toward +x: +1 upward with `q` on their
/// left, −1 downward with `q` on their right, each decided exactly.
fn winding(sides: &[(DVec2, DVec2)], q: DVec2) -> i64 {
    let mut w = 0;
    for &(a, b) in sides {
        if a.y <= q.y {
            if b.y > q.y && orient2(a, b, q) > 0.0 {
                w += 1;
            }
        } else if b.y <= q.y && orient2(a, b, q) < 0.0 {
            w -= 1;
        }
    }
    w
}

/// The centre and radius of the circle inscribed in a triangle, or `None`
/// for a triangle with no perimeter.
fn incircle([a, b, c]: [DVec2; 3]) -> Option<(DVec2, f64)> {
    let (la, lb, lc) = ((b - c).length(), (c - a).length(), (a - b).length());
    let per = la + lb + lc;
    let area = 0.5 * (b - a).perp_dot(c - a).abs();
    (per > 0.0).then(|| ((a * la + b * lb + c * lc) / per, 2.0 * area / per))
}

/// For a convex counter-clockwise polygon ([`tess::convexity`]), a point
/// inside it and its distance from the boundary at least, without a
/// triangulation; `None` for any other polygon.
///
/// The point is the better of two candidates. The centroid of the area is
/// at least a third of the polygon's width from every side (Minkowski and
/// Radon), so at least two thirds as far from the boundary as the centre
/// of the largest inscribed circle; and the incentres of the triangles of
/// [`tess::balanced`], of which a triangle's own incentre is the best
/// possible. A point's distance from a convex polygon's boundary is its
/// least distance from the sides' lines.
fn convex_interior(l: &[DVec2]) -> Option<(DVec2, f64)> {
    tess::convexity(l)?;
    let n = l.len();
    let o = l[0];
    let (mut twice_area, mut sum) = (0.0, DVec2::ZERO);
    for k in 1..n - 1 {
        let (p, q) = (l[k] - o, l[k + 1] - o);
        let w = p.perp_dot(q);
        twice_area += w;
        sum += (p + q) * w;
    }
    if !twice_area.is_finite() || twice_area <= 0.0 {
        return None;
    }
    // Least distance from the sides' lines, positive inside.
    let clearance = |p: DVec2| {
        (0..n)
            .map(|k| {
                let (a, b) = (l[k], l[(k + 1) % n]);
                (b - a).perp_dot(p - a) / (b - a).length()
            })
            .fold(f64::INFINITY, f64::min)
    };
    let centroid = o + sum / (3.0 * twice_area);
    let mut best = (centroid, clearance(centroid));
    for t in tess::balanced(n) {
        if let Some((c, r)) = incircle(t.map(|k| l[k])) {
            if r > best.1 {
                best = (c, r);
            }
        }
    }
    (best.1 > 0.0).then_some(best)
}

/// The point of triangle `t` closest to `p` (Ericson, *Real-Time Collision
/// Detection*, 2005, §5.1.5: by the Voronoi region of `p`).
fn closest_on_triangle(p: DVec3, t: &[DVec3; 3]) -> DVec3 {
    let [a, b, c] = *t;
    let (ab, ac, ap) = (b - a, c - a, p - a);
    let (d1, d2) = (ab.dot(ap), ac.dot(ap));
    if d1 <= 0.0 && d2 <= 0.0 {
        return a;
    }
    let bp = p - b;
    let (d3, d4) = (ab.dot(bp), ac.dot(bp));
    if d3 >= 0.0 && d4 <= d3 {
        return b;
    }
    let vc = d1 * d4 - d3 * d2;
    if vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0 {
        return a + ab * (d1 / (d1 - d3));
    }
    let cp = p - c;
    let (d5, d6) = (ab.dot(cp), ac.dot(cp));
    if d6 >= 0.0 && d5 <= d6 {
        return c;
    }
    let vb = d5 * d2 - d1 * d6;
    if vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0 {
        return a + ac * (d2 / (d2 - d6));
    }
    let va = d3 * d6 - d5 * d4;
    if va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0 {
        return b + (c - b) * ((d4 - d3) / ((d4 - d3) + (d5 - d6)));
    }
    let den = va + vb + vc;
    if den == 0.0 {
        // A degenerate triangle: its closest edge point.
        return [(a, b), (b, c), (c, a)]
            .iter()
            .map(|&(u, v)| {
                let d = v - u;
                let s =
                    ((p - u).dot(d) / d.length_squared().max(f64::MIN_POSITIVE)).clamp(0.0, 1.0);
                u + d * s
            })
            .min_by(|x, y| (*x - p).length().total_cmp(&(*y - p).length()))
            .unwrap_or(a);
    }
    a + ab * (vb / den) + ac * (vc / den)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn closest_points_on_a_triangle() {
        let t = [DVec3::ZERO, DVec3::X, DVec3::Y];
        let cases = [
            (DVec3::new(0.2, 0.2, 1.0), DVec3::new(0.2, 0.2, 0.0)),
            (DVec3::new(-1.0, -1.0, 0.5), DVec3::ZERO),
            (DVec3::new(2.0, -0.5, 0.0), DVec3::X),
            (DVec3::new(0.5, -1.0, 0.0), DVec3::new(0.5, 0.0, 0.0)),
            (DVec3::new(-1.0, 0.5, 0.0), DVec3::new(0.0, 0.5, 0.0)),
            (DVec3::new(1.0, 1.0, 3.0), DVec3::new(0.5, 0.5, 0.0)),
            (DVec3::new(-0.5, 3.0, 0.0), DVec3::Y),
        ];
        for (p, want) in cases {
            assert!((closest_on_triangle(p, &t) - want).length() < 1e-15, "{p}");
        }
        let flat = [DVec3::ZERO, DVec3::X, DVec3::X * 2.0];
        let q = closest_on_triangle(DVec3::new(1.5, 1.0, 0.0), &flat);
        assert!((q - DVec3::new(1.5, 0.0, 0.0)).length() < 1e-15);
    }

    fn poly(pts: &[(f64, f64)]) -> Vec<DVec2> {
        pts.iter().map(|&(x, y)| DVec2::new(x, y)).collect()
    }

    /// The distance from `p` to the polygon's boundary.
    fn boundary_distance(l: &[DVec2], p: DVec2) -> f64 {
        (0..l.len())
            .map(|k| {
                let (a, b) = (l[k], l[(k + 1) % l.len()]);
                let t = ((p - a).dot(b - a) / (b - a).length_squared()).clamp(0.0, 1.0);
                (a + (b - a) * t - p).length()
            })
            .fold(f64::INFINITY, f64::min)
    }

    /// True if `p` is strictly inside the counter-clockwise polygon `l`
    /// (by winding number).
    fn strictly_inside(l: &[DVec2], p: DVec2) -> bool {
        let mut w = 0;
        for k in 0..l.len() {
            let (a, b) = (l[k], l[(k + 1) % l.len()]);
            let o = orient2(a, b, p);
            if o == 0.0 && (p - a).dot(p - b) <= 0.0 {
                return false;
            }
            if a.y <= p.y && b.y > p.y && o > 0.0 {
                w += 1;
            } else if a.y > p.y && b.y <= p.y && o < 0.0 {
                w -= 1;
            }
        }
        w == 1
    }

    #[test]
    fn convex_pieces_get_a_point_without_a_triangulation() {
        // A square: its centre, half its side from the boundary.
        let square = poly(&[(0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 2.0)]);
        let (p, r) = convex_interior(&square).unwrap();
        assert!((p - DVec2::ONE).length() < 1e-15 && (r - 1.0).abs() < 1e-15);
        // With points on its sides (where other faces split them).
        let split = poly(&[
            (0.0, 0.0),
            (1.0, 0.0),
            (2.0, 0.0),
            (2.0, 2.0),
            (1.0, 2.0),
            (0.0, 2.0),
            (0.0, 1.0),
        ]);
        let (p, r) = convex_interior(&split).unwrap();
        assert!(strictly_inside(&split, p) && r > 0.9 && r <= boundary_distance(&split, p));
        // A thin triangle: its incentre, half as far again from the long
        // sides as its centroid.
        let thin = poly(&[(0.0, 0.0), (10.0, 0.0), (5.0, 0.1)]);
        let (p, r) = convex_interior(&thin).unwrap();
        let (c, ri) = incircle([thin[0], thin[1], thin[2]]).unwrap();
        assert_eq!((p, r), (c, ri));
        assert!((r - boundary_distance(&thin, p)).abs() < 1e-12);
        // Anything else is left to the triangulation.
        let refused = [
            // Clockwise.
            poly(&[(0.0, 0.0), (0.0, 1.0), (1.0, 1.0), (1.0, 0.0)]),
            // Not convex.
            poly(&[(0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (1.0, 1.0), (0.0, 2.0)]),
            // A spike straight back along a side.
            poly(&[(0.0, 0.0), (2.0, 0.0), (3.0, 0.0), (2.0, 0.0), (2.0, 2.0)]),
            // A repeated point.
            poly(&[(0.0, 0.0), (1.0, 0.0), (1.0, 0.0), (1.0, 1.0)]),
            // All on one line.
            poly(&[(0.0, 0.0), (1.0, 0.0), (2.0, 0.0)]),
            // A pentagram: every turn to the left, but around twice.
            (0..5)
                .map(|k| {
                    let a = std::f64::consts::TAU * (2 * k) as f64 / 5.0;
                    DVec2::new(a.cos(), a.sin())
                })
                .collect(),
        ];
        for l in &refused {
            assert_eq!(convex_interior(l), None, "{l:?}");
        }
    }

    /// A circle of radius `r` about `c` as `n` points, clockwise (a hole).
    fn hole(c: DVec2, r: f64, n: usize) -> Vec<DVec2> {
        (0..n)
            .map(|k| c + DVec2::from_angle(-std::f64::consts::TAU * k as f64 / n as f64) * r)
            .collect()
    }

    #[test]
    fn probes_find_a_point_among_round_holes() {
        // A plate with a 10 by 10 grid of round holes, the face a boolean
        // drilling it splits: the probes find a point between four holes,
        // as far from them as any point of the plate.
        let mut loops = vec![poly(&[(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0)])];
        for i in 0..10 {
            for j in 0..10 {
                loops.push(hole(DVec2::new(i as f64 + 0.5, j as f64 + 0.5), 0.3, 32));
            }
        }
        let (p, r) = probe_point(&loops).unwrap();
        let best = 0.5f64.sqrt() - 0.3;
        assert!(r > 0.99 * best && r <= best, "{p} {r}");
        assert!(strictly_inside(&loops[0], p));
        for h in &loops[1..] {
            let mut ccw = h.clone();
            ccw.reverse();
            assert!(!strictly_inside(&ccw, p));
        }
        // Too few sides to bound anything.
        assert_eq!(probe_point(&[poly(&[(0.0, 0.0), (1.0, 0.0)])]), None);
    }

    proptest::proptest! {
        #![proptest_config(proptest::prelude::ProptestConfig::with_cases(400))]

        /// On random polygons: a convex polygon always gets a point, and
        /// every point given is inside, its clearance a true lower bound
        /// and at least half the triangulation's; a polygon that is not
        /// convex never gets one.
        #[test]
        fn probe_points_are_inside(
            n in 3usize..40,
            holes in 0usize..12,
            seed in 0u64..u64::MAX,
        ) {
            let mut h = seed;
            let mut rnd = || {
                h = h.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
                (h >> 11) as f64 / (1u64 << 53) as f64
            };
            // A star-shaped outer loop about the origin, radii 1 to 2.
            let outer: Vec<DVec2> = (0..n)
                .map(|k| {
                    let a = std::f64::consts::TAU * k as f64 / n as f64;
                    DVec2::from_angle(a) * (1.0 + rnd())
                })
                .collect();
            let mut loops = vec![outer];
            let mut discs: Vec<(DVec2, f64)> = Vec::new();
            for _ in 0..holes {
                let c = DVec2::new(rnd() - 0.5, rnd() - 0.5) * 1.2;
                let r = 0.02 + 0.15 * rnd();
                if discs.iter().all(|&(d, s)| (c - d).length() > r + s + 1e-3)
                    && boundary_distance(&loops[0], c) > r + 1e-3
                {
                    discs.push((c, r));
                    loops.push(hole(c, r, 3 + (rnd() * 40.0) as usize));
                }
            }
            let (p, r) = probe_point(&loops).unwrap();
            proptest::prop_assert!(strictly_inside(&loops[0], p));
            let mut d = boundary_distance(&loops[0], p);
            for l in &loops[1..] {
                let mut ccw = l.clone();
                ccw.reverse();
                proptest::prop_assert!(!strictly_inside(&ccw, p));
                d = d.min(boundary_distance(l, p));
            }
            proptest::prop_assert!(r > 0.0 && (r - d).abs() <= 1e-12, "{r} against {d}");
        }

        #[test]
        fn convex_interior_points_are_inside(
            n in 3usize..40,
            seed in 0u64..u64::MAX,
            star in proptest::bool::ANY,
        ) {
            let mut h = seed;
            let mut rnd = || {
                h = h.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
                (h >> 11) as f64 / (1u64 << 53) as f64
            };
            let mut angles: Vec<f64> = (0..n).map(|_| rnd() * std::f64::consts::TAU).collect();
            angles.sort_by(f64::total_cmp);
            let (cx, cy, s) = (rnd() * 20.0 - 10.0, rnd() * 20.0 - 10.0, 0.01 + rnd() * 5.0);
            let l: Vec<DVec2> = angles
                .iter()
                .map(|&a| {
                    let r = if star { s * (0.3 + rnd()) } else { s };
                    DVec2::new(cx + r * a.cos(), cy + r * a.sin())
                })
                .collect();
            let convex = (0..n).all(|k| {
                (0..n).all(|j| orient2(l[k], l[(k + 1) % n], l[j]) >= 0.0)
            }) && (0..n).all(|k| l[k] != l[(k + 1) % n]);
            match convex_interior(&l) {
                Some((p, r)) => {
                    proptest::prop_assert!(convex, "not convex: {l:?}");
                    proptest::prop_assert!(strictly_inside(&l, p));
                    let d = boundary_distance(&l, p);
                    proptest::prop_assert!(r > 0.0 && r <= d * (1.0 + 1e-9) + 1e-15);
                    let (_, rc) = interior_point(std::slice::from_ref(&l)).unwrap();
                    let cdt = {
                        let (cdt, index) = tess::polygon_cdt(std::slice::from_ref(&l)).unwrap();
                        let mut at = vec![DVec2::ZERO; cdt.num_vertices()];
                        for (&i, q) in index.iter().zip(&l) {
                            at[i] = *q;
                        }
                        tess::inside_triangles(&cdt)
                            .into_iter()
                            .filter_map(|t| incircle(t.map(|i| at[i])))
                            .map(|(_, r)| r)
                            .fold(0.0, f64::max)
                    };
                    proptest::prop_assert_eq!(rc, r);
                    proptest::prop_assert!(r >= 0.5 * cdt, "{r} against {cdt}");
                }
                None => proptest::prop_assert!(!convex || angles.windows(2).any(|w| w[0] == w[1])),
            }
        }
    }
}
