//! Classifying the pieces of the faces (design note 0012 §4.3).
//!
//! A piece that lies on a face of the other body (the two faces on one
//! plane) is *on* it, with the same or the opposite outward normal. Any
//! other piece is wholly inside or wholly outside the other body, since
//! every place where that body's surface crosses the face is a cut between
//! pieces; it is classified by the generalised winding number of the other
//! body's (exact, planar) triangulation at a point well inside the piece:
//! the incentre of the piece's triangle with the largest inscribed circle.
//! The winding number is an integer (0 or 1) up to rounding at any point
//! not on the surface; a value far from both is an error, never a guess.
//! A point within the tolerance of the other body's surface is on it, so a
//! piece whose point is there without lying on a coincident face has no
//! side the tolerance can decide, and the boolean is refused.
//!
//! This replaces the propagation of design note 0012 §4.3 (classify faces
//! at the intersection edges, then flood fill): each piece is decided by
//! one test whose robustness is measured by its distance from the surface,
//! and no decision depends on another.

use std::f64::consts::PI;

use glam::{DVec2, DVec3};
use spade::{ConstrainedDelaunayTriangulation, Point2, Triangulation};

use super::show;
use crate::tess::{self, TessellationOptions};
use crate::topo::Brep;
use crate::{KernelError, Tolerance};

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

/// A body's boundary as triangles, for winding numbers.
pub(super) struct Winding {
    tris: Vec<[DVec3; 3]>,
}

impl Winding {
    pub fn new(b: &Brep, tol: &Tolerance) -> Result<Winding, KernelError> {
        // Planar faces with straight edges: the triangulation is exact and
        // needs no refinement, whatever the tolerances.
        let size = crate::check::bounding_size(b).max(tol.linear);
        let opts = TessellationOptions {
            chord_tolerance: size,
            max_angle: PI,
        };
        let t = tess::tessellate(b, &opts, tol)?;
        let mut tris = Vec::new();
        for f in &t.faces {
            for tri in &f.triangles {
                tris.push(tri.map(|i| f.positions[i as usize]));
            }
        }
        Ok(Winding { tris })
    }

    /// The winding number of the boundary about `p`.
    fn number(&self, p: DVec3) -> f64 {
        let mut total = 0.0;
        for t in &self.tris {
            let (a, b, c) = (t[0] - p, t[1] - p, t[2] - p);
            let (la, lb, lc) = (a.length(), b.length(), c.length());
            let num = a.dot(b.cross(c));
            let den = la * lb * lc + a.dot(b) * lc + a.dot(c) * lb + b.dot(c) * la;
            total += 2.0 * num.atan2(den);
        }
        total / (4.0 * PI)
    }

    /// The distance from `p` to the boundary.
    fn distance(&self, p: DVec3) -> f64 {
        self.tris
            .iter()
            .map(|t| (closest_on_triangle(p, t) - p).length())
            .fold(f64::INFINITY, f64::min)
    }

    /// Inside or outside, at a point farther than `eps` from the boundary;
    /// a point within it is on the boundary, and a piece of a face that is
    /// on it without coinciding with a face is refused.
    pub fn classify(&self, p: DVec3, eps: f64) -> Result<Class, KernelError> {
        if self.distance(p) <= eps {
            return Err(KernelError::InvalidInput(format!(
                "a face of one body runs within the tolerance of the other's surface \
                 without lying on one of its faces, near {}",
                show(p)
            )));
        }
        let w = self.number(p);
        if w.abs() < 0.25 {
            Ok(Class::Outside)
        } else if (w - 1.0).abs() < 0.25 {
            Ok(Class::Inside)
        } else {
            Err(KernelError::InvalidInput(format!(
                "cannot tell whether a piece near {} is inside the other body \
                 (winding number {w:.3})",
                show(p)
            )))
        }
    }
}

/// A point well inside the polygon `loops` (outer first) and its distance
/// from the boundary at least: the incentre of the triangle with the
/// largest inscribed circle in a constrained triangulation of the polygon.
pub(super) fn interior_point(loops: &[Vec<DVec2>]) -> Result<(DVec2, f64), String> {
    let mut cdt: ConstrainedDelaunayTriangulation<Point2<f64>> =
        ConstrainedDelaunayTriangulation::new();
    let mut at: Vec<DVec2> = Vec::new();
    for l in loops {
        let mut hs = Vec::with_capacity(l.len());
        for q in l {
            let h = cdt
                .insert(Point2::new(q.x, q.y))
                .map_err(|e| format!("triangulation insert failed: {e:?}"))?;
            if h.index() == at.len() {
                at.push(*q);
            }
            hs.push(h);
        }
        for i in 0..hs.len() {
            let (a, b) = (hs[i], hs[(i + 1) % hs.len()]);
            if a != b && cdt.try_add_constraint(a, b).is_empty() && !cdt.exists_constraint(a, b) {
                return Err("a piece's boundary crosses itself".into());
            }
        }
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
