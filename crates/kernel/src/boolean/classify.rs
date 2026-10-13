//! Classifying the pieces of the faces (design note 0012 §4.3).
//!
//! A piece that lies on a face of the other body (the two faces on one
//! plane) is *on* it, with the same or the opposite outward normal. Any
//! other piece is wholly inside or wholly outside the other body, since
//! every place where that body's surface crosses the face is a cut between
//! pieces; it is classified at one point well inside it (the incentre of
//! the piece's triangle with the largest inscribed circle).
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
//! This replaces the propagation of design note 0012 §4.3 (classify faces
//! at the intersection edges, then flood fill): each piece is decided by
//! one exact test, and no decision depends on another.

use glam::{DVec2, DVec3};
use spade::Triangulation;

use super::{show, Poly};
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

/// A body's boundary as triangles on its own vertices, outward facing.
pub(super) struct Boundary<'a> {
    points: &'a [DVec3],
    tris: Vec<[u32; 3]>,
    tree: Bvh,
    /// The centre and half the diagonal of the body's box.
    centre: DVec3,
    reach: f64,
    /// Padding of the boxes in the segment test, against rounding.
    pad: f64,
}

impl<'a> Boundary<'a> {
    pub fn new(poly: &'a Poly) -> Result<Boundary<'a>, KernelError> {
        let points = &poly.points[..];
        let mut tris: Vec<[u32; 3]> = Vec::new();
        for (fi, f) in poly.faces.iter().enumerate() {
            triangulate(poly, f, &mut tris).map_err(|m| {
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
        let (lo, hi) = tree.bounds().unwrap_or((DVec3::ZERO, DVec3::ZERO));
        let centre = (lo + hi) * 0.5;
        let reach = (hi - lo).length() * 0.5;
        let pad = 1e-9 * (reach + centre.abs().max_element()) + f64::MIN_POSITIVE;
        Ok(Boundary {
            points,
            tris,
            tree,
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
        )
    }

    /// The signed number of triangles the segment `p`–`q` passes through
    /// (+1 leaving the body, −1 entering it), or `None` if it passes
    /// through an edge or a corner of one, or lies in a triangle's plane
    /// across it. `p` is off the boundary; `q` is outside the body's box.
    fn crossings(&self, p: DVec3, q: DVec3) -> Option<i64> {
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
            match self.crossings(p, q) {
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
fn triangulate(poly: &Poly, f: &super::PolyFace, out: &mut Vec<[u32; 3]>) -> Result<(), String> {
    let q2 = |v: usize| f.to2(poly.points[v]);
    // A strictly convex face without holes, counter-clockwise: triangles
    // that halve the boundary's runs level by level. Each triangle spans a
    // run of the boundary and its box is about the size of that run, so a
    // ray meets a few boxes per level; the slivers of a fan all reach
    // across the face, and a ray near it meets most of them.
    if let [vs] = &f.verts[..] {
        let n = vs.len();
        let convex = n >= 3
            && (0..n).all(|k| orient2(q2(vs[k]), q2(vs[(k + 1) % n]), q2(vs[(k + 2) % n])) > 0.0);
        if convex {
            let tri = |a: usize, b: usize, c: usize| [vs[a], vs[b], vs[c % n]].map(|v| v as u32);
            out.push(tri(0, n / 3, 2 * n / 3));
            let mut runs = vec![(0, n / 3), (n / 3, 2 * n / 3), (2 * n / 3, n)];
            while let Some((i, j)) = runs.pop() {
                if j - i >= 2 {
                    let m = (i + j) / 2;
                    out.push(tri(i, m, j));
                    runs.push((m, j));
                    runs.push((i, m));
                }
            }
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
/// from the boundary at least: the incentre of the triangle with the
/// largest inscribed circle in a constrained triangulation of the polygon.
pub(super) fn interior_point(loops: &[Vec<DVec2>]) -> Result<(DVec2, f64), String> {
    let (cdt, index) = tess::polygon_cdt(loops).map_err(|e| format!("a piece's boundary: {e}"))?;
    let mut at = vec![DVec2::ZERO; cdt.num_vertices()];
    for (&i, q) in index.iter().zip(loops.iter().flatten()) {
        at[i] = *q;
    }
    let mut best: Option<(DVec2, f64)> = None;
    for t in tess::inside_triangles(&cdt) {
        let [a, b, c] = t.map(|i| at[i]);
        let (la, lb, lc) = ((b - c).length(), (c - a).length(), (a - b).length());
        let per = la + lb + lc;
        let area = 0.5 * (b - a).perp_dot(c - a).abs();
        if per <= 0.0 {
            continue;
        }
        let r = 2.0 * area / per;
        if best.is_none_or(|(_, br)| r > br) {
            best = Some(((a * la + b * lb + c * lc) / per, r));
        }
    }
    best.ok_or_else(|| "a piece with no area".into())
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
}
