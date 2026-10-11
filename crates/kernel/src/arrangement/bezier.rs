//! Rational Bézier pieces in homogeneous form, the working representation
//! for intersection, winding and area computations.
//!
//! References: G. Farin, *Curves and Surfaces for CAGD*, 5th ed., ch. 4–5
//! (de Casteljau algorithm, derivatives, subdivision) and ch. 13 (rational
//! Bézier curves); L. Piegl & W. Tiller, *The NURBS Book*, 2nd ed., §1.4
//! and §5.6 (decomposition of B-splines into Bézier segments via blossoms).

use glam::{DVec2, DVec3};

/// A rational Bézier curve on the local parameter `u ∈ [0, 1]`, with
/// homogeneous control points `(w·x, w·y, w)` and positive weights. By the
/// convex hull property the curve lies in the convex hull of the projected
/// control points.
#[derive(Clone, Debug)]
pub(super) struct Bez {
    /// Homogeneous control points.
    pub pts: Vec<DVec3>,
}

/// Projects a homogeneous point.
#[inline]
pub(super) fn project(h: DVec3) -> DVec2 {
    DVec2::new(h.x / h.z, h.y / h.z)
}

/// De Casteljau evaluation of a (vector-valued) Bernstein polynomial.
fn casteljau(pts: &[DVec3], u: f64) -> DVec3 {
    let mut b = pts.to_vec();
    let n = b.len();
    for k in 1..n {
        for i in 0..n - k {
            b[i] = b[i] * (1.0 - u) + b[i + 1] * u;
        }
    }
    b[0]
}

impl Bez {
    /// The straight segment from `a` to `b`.
    pub fn line(a: DVec2, b: DVec2) -> Bez {
        Bez {
            pts: vec![a.extend(1.0), b.extend(1.0)],
        }
    }

    /// The exact rational quadratic arc of the circle `(c, r)` from angle `f0`
    /// to angle `f1` (`|f1 − f0| < π`), Farin §13.6.
    pub fn arc(c: DVec2, r: f64, f0: f64, f1: f64) -> Bez {
        let h = 0.5 * (f1 - f0);
        let w = h.cos();
        let p0 = c + r * DVec2::from_angle(f0);
        let p2 = c + r * DVec2::from_angle(f1);
        let m = c + (r / w) * DVec2::from_angle(f0 + h);
        Bez {
            pts: vec![p0.extend(1.0), (m * w).extend(w), p2.extend(1.0)],
        }
    }

    /// Polynomial degree.
    pub fn degree(&self) -> usize {
        self.pts.len() - 1
    }

    /// Point at `u = 0`.
    pub fn start(&self) -> DVec2 {
        project(self.pts[0])
    }

    /// Point at `u = 1`.
    pub fn end(&self) -> DVec2 {
        project(self.pts[self.pts.len() - 1])
    }

    /// Projected control points (the convex hull vertices, unsorted).
    pub fn hull(&self) -> Vec<DVec2> {
        self.pts.iter().map(|&h| project(h)).collect()
    }

    /// True if all weights are equal (a polynomial curve).
    pub fn is_polynomial(&self) -> bool {
        self.pts.iter().all(|h| h.z == self.pts[0].z)
    }

    /// Point at `u`.
    pub fn point(&self, u: f64) -> DVec2 {
        project(casteljau(&self.pts, u))
    }

    /// Point and first and second derivatives with respect to `u`
    /// (quotient rule on the homogeneous curve, Farin §13.4).
    pub fn derivs(&self, u: f64) -> [DVec2; 3] {
        let n = self.degree();
        let h0 = casteljau(&self.pts, u);
        let d1: Vec<DVec3> = (0..n)
            .map(|i| (self.pts[i + 1] - self.pts[i]) * n as f64)
            .collect();
        let h1 = casteljau(&d1, u);
        let h2 = if n >= 2 {
            let d2: Vec<DVec3> = (0..n - 1)
                .map(|i| (d1[i + 1] - d1[i]) * (n - 1) as f64)
                .collect();
            casteljau(&d2, u)
        } else {
            DVec3::ZERO
        };
        let w = h0.z;
        let c = h0.truncate() / w;
        let c1 = (h1.truncate() - c * h1.z) / w;
        let c2 = (h2.truncate() - c1 * (2.0 * h1.z) - c * h2.z) / w;
        [c, c1, c2]
    }

    /// Splits at `u` into the pieces on `[0, u]` and `[u, 1]`, each
    /// reparametrized to `[0, 1]`.
    pub fn split(&self, u: f64) -> (Bez, Bez) {
        let n = self.pts.len();
        let mut b = self.pts.clone();
        let mut left = Vec::with_capacity(n);
        let mut right = vec![DVec3::ZERO; n];
        left.push(b[0]);
        right[n - 1] = b[n - 1];
        for k in 1..n {
            for i in 0..n - k {
                b[i] = b[i] * (1.0 - u) + b[i + 1] * u;
            }
            left.push(b[0]);
            right[n - 1 - k] = b[n - 1 - k];
        }
        (Bez { pts: left }, Bez { pts: right })
    }

    /// The piece on `[u0, u1]` (`0 ≤ u0 < u1 ≤ 1`), reparametrized to `[0, 1]`.
    pub fn sub(&self, u0: f64, u1: f64) -> Bez {
        let left = if u1 < 1.0 {
            self.split(u1).0
        } else {
            self.clone()
        };
        if u0 > 0.0 {
            left.split(u0 / u1).1
        } else {
            left
        }
    }

    /// The same curve traversed backwards.
    pub fn reversed(&self) -> Bez {
        Bez {
            pts: self.pts.iter().rev().copied().collect(),
        }
    }
}

/// Axis-aligned bounds of a point set.
pub(super) fn bounds(pts: &[DVec2]) -> (DVec2, DVec2) {
    let mut lo = DVec2::splat(f64::INFINITY);
    let mut hi = DVec2::splat(f64::NEG_INFINITY);
    for p in pts {
        lo = lo.min(*p);
        hi = hi.max(*p);
    }
    (lo, hi)
}

/// Distance from `p` to the segment `ab`.
pub(super) fn point_segment_distance(p: DVec2, a: DVec2, b: DVec2) -> f64 {
    let d = b - a;
    let l2 = d.length_squared();
    if l2 == 0.0 {
        return (p - a).length();
    }
    let s = ((p - a).dot(d) / l2).clamp(0.0, 1.0);
    (p - (a + d * s)).length()
}

/// Largest distance of the hull points from the chord between the first and
/// last of them: a bound on how far the curve strays from its chord.
pub(super) fn flatness(hull: &[DVec2]) -> f64 {
    let a = hull[0];
    let b = hull[hull.len() - 1];
    hull[1..hull.len() - 1]
        .iter()
        .map(|&p| point_segment_distance(p, a, b))
        .fold(0.0, f64::max)
}

/// Gauss–Legendre nodes and weights on `[-1, 1]` (Newton iteration on the
/// Legendre polynomial; Abramowitz & Stegun 25.4.29).
pub(super) fn gauss_legendre(n: usize) -> Vec<(f64, f64)> {
    let mut out = Vec::with_capacity(n);
    for i in 0..n {
        let mut x = (std::f64::consts::PI * (i as f64 + 0.75) / (n as f64 + 0.5)).cos();
        let mut dp = 1.0;
        for _ in 0..100 {
            let (mut p0, mut p1) = (1.0, x);
            for k in 2..=n {
                let p2 = ((2 * k - 1) as f64 * x * p1 - (k - 1) as f64 * p0) / k as f64;
                p0 = p1;
                p1 = p2;
            }
            dp = n as f64 * (x * p1 - p0) / (x * x - 1.0);
            let dx = p1 / dp;
            x -= dx;
            if dx.abs() < 1e-16 {
                break;
            }
        }
        out.push((x, 2.0 / ((1.0 - x * x) * dp * dp)));
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn arc_is_exact_and_split_preserves_geometry() {
        let b = Bez::arc(DVec2::new(1.0, 2.0), 3.0, 0.2, 1.4);
        for i in 0..=10 {
            let u = i as f64 / 10.0;
            let p = b.point(u);
            assert!(((p - DVec2::new(1.0, 2.0)).length() - 3.0).abs() < 1e-14);
        }
        let (l, r) = b.split(0.3);
        assert!((l.point(0.5) - b.point(0.15)).length() < 1e-14);
        assert!((r.point(0.5) - b.point(0.65)).length() < 1e-14);
        let s = b.sub(0.2, 0.7);
        assert!((s.point(0.5) - b.point(0.45)).length() < 1e-14);
        let [p, d1, d2] = b.derivs(0.4);
        let h = 1e-6;
        let fd1 = (b.point(0.4 + h) - b.point(0.4 - h)) / (2.0 * h);
        let fd2 = (b.point(0.4 + h) - 2.0 * p + b.point(0.4 - h)) / (h * h);
        assert!((d1 - fd1).length() < 1e-6);
        assert!((d2 - fd2).length() < 1e-3);
    }

    #[test]
    fn gauss_integrates_polynomials() {
        let g = gauss_legendre(5);
        let s: f64 = g.iter().map(|&(x, w)| w * x.powi(8)).sum();
        assert!((s - 2.0 / 9.0).abs() < 1e-14);
    }
}
