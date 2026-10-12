//! Curve and surface geometry carried by B-rep edges and faces.
//!
//! Analytic types are kept exact wherever possible (lines, circles, planes,
//! surfaces of revolution and linear extrusion), with NURBS for free-form
//! geometry. Every type supports evaluation, first derivatives, projection
//! (closest point) and affine transformation.

use glam::{DMat4, DVec3};
use openshape_nurbs::{NurbsCurve, NurbsSurface};
use serde::{Deserialize, Serialize};

use crate::math::{wrap_angle, Frame};

/// A 3D curve. Parameter ranges are carried by the edges that use it.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub(crate) enum Curve3 {
    /// `origin + t * dir` (dir need not be unit length).
    Line { origin: DVec3, dir: DVec3 },
    /// `center + r (cos t · x + sin t · y)` with orthonormal `x`, `y`.
    /// Parameters are angles in `[-π/2, 3π/2)`.
    Circle {
        center: DVec3,
        x: DVec3,
        y: DVec3,
        radius: f64,
    },
    /// Non-uniform rational B-spline.
    Nurbs(NurbsCurve),
}

/// Knot spans shorter than this fraction of a parameter range are merged
/// into their neighbours by [`Curve3::smooth_breaks`].
pub const KNOT_MERGE: f64 = 1e-12;

impl Curve3 {
    /// Point at parameter `t`.
    pub fn point(&self, t: f64) -> DVec3 {
        match self {
            Curve3::Line { origin, dir } => *origin + *dir * t,
            Curve3::Circle {
                center,
                x,
                y,
                radius,
            } => {
                let (s, c) = t.sin_cos();
                *center + (*x * c + *y * s) * *radius
            }
            Curve3::Nurbs(n) => n.point(t),
        }
    }

    /// First derivative at `t`.
    pub fn deriv(&self, t: f64) -> DVec3 {
        match self {
            Curve3::Line { dir, .. } => *dir,
            Curve3::Circle { x, y, radius, .. } => {
                let (s, c) = t.sin_cos();
                (*y * c - *x * s) * *radius
            }
            Curve3::Nurbs(n) => n.derivatives(t, 1)[1],
        }
    }

    /// Parameter of the closest point to `p`.
    pub fn project(&self, p: DVec3) -> f64 {
        match self {
            Curve3::Line { origin, dir } => (p - *origin).dot(*dir) / dir.length_squared(),
            Curve3::Circle { center, x, y, .. } => {
                let d = p - *center;
                wrap_angle(d.dot(*y).atan2(d.dot(*x)))
            }
            Curve3::Nurbs(n) => n.closest_point(p).0,
        }
    }

    /// The curve under an affine transform.
    pub fn transformed(&self, m: &DMat4) -> Curve3 {
        match self {
            Curve3::Line { origin, dir } => Curve3::Line {
                origin: m.transform_point3(*origin),
                dir: m.transform_vector3(*dir),
            },
            Curve3::Circle {
                center,
                x,
                y,
                radius,
            } => {
                // Rigid motions keep circles circular; the kernel only applies
                // rigid transforms to circles.
                let nx = m.transform_vector3(*x);
                let scale = nx.length();
                Curve3::Circle {
                    center: m.transform_point3(*center),
                    x: nx / scale,
                    y: m.transform_vector3(*y) / scale,
                    radius: radius * scale,
                }
            }
            Curve3::Nurbs(n) => Curve3::Nurbs(n.transformed(m)),
        }
    }

    /// Parameters splitting `[t0, t1]` into pieces on which the curve is
    /// one polynomial or rational piece (for quadrature): a line whole, a
    /// circle in pieces of at most π/4, a NURBS curve at its distinct knots
    /// inside the range. Knots within [`KNOT_MERGE`] of the range's length
    /// of a neighbouring break are skipped: such a knot span is a rounding
    /// artefact (a split one rounding step from an existing knot), its
    /// polynomial piece is numerically meaningless, and its contribution
    /// to any integral is negligible.
    pub fn smooth_breaks(&self, t0: f64, t1: f64) -> Vec<f64> {
        match self {
            Curve3::Line { .. } => vec![t0, t1],
            Curve3::Circle { .. } => {
                let n = (((t1 - t0).abs() / std::f64::consts::FRAC_PI_4).ceil() as usize).max(1);
                (0..=n)
                    .map(|i| t0 + (t1 - t0) * i as f64 / n as f64)
                    .collect()
            }
            Curve3::Nurbs(c) => {
                let (lo, hi) = (t0.min(t1), t0.max(t1));
                let gap = KNOT_MERGE * (hi - lo);
                let mut out = vec![lo];
                for &k in c.knots().as_slice() {
                    if k > lo + gap && k < hi - gap && out.last().is_some_and(|&l| k > l + gap) {
                        out.push(k);
                    }
                }
                out.push(hi);
                if t0 > t1 {
                    out.reverse();
                }
                out
            }
        }
    }

    /// Parameters subdividing `[t0, t1]` so that the polyline through them
    /// deviates from the curve by at most `chord_tol` and turns by at most
    /// `max_angle` per segment.
    pub fn discretize(&self, t0: f64, t1: f64, chord_tol: f64, max_angle: f64) -> Vec<f64> {
        match self {
            Curve3::Line { .. } => vec![t0, t1],
            Curve3::Circle { radius, .. } => {
                let sweep = t1 - t0;
                let by_chord = if chord_tol < *radius {
                    2.0 * (1.0 - chord_tol / radius).acos()
                } else {
                    std::f64::consts::PI
                };
                let step = by_chord.min(max_angle).max(1e-4);
                let n = ((sweep / step).ceil() as usize).max(1);
                (0..=n).map(|i| t0 + sweep * i as f64 / n as f64).collect()
            }
            Curve3::Nurbs(_) => {
                // Start from a uniform split, then bisect adaptively.
                let n0 = 8;
                let mut out = vec![t0];
                for i in 0..n0 {
                    let a = t0 + (t1 - t0) * i as f64 / n0 as f64;
                    let b = t0 + (t1 - t0) * (i + 1) as f64 / n0 as f64;
                    self.refine(a, b, chord_tol, max_angle, 0, &mut out);
                }
                out
            }
        }
    }

    fn refine(&self, a: f64, b: f64, tol: f64, max_angle: f64, depth: u32, out: &mut Vec<f64>) {
        let pa = self.point(a);
        let pb = self.point(b);
        let mut split = false;
        if depth < 16 {
            for f in [0.25, 0.5, 0.75] {
                let tm = a + (b - a) * f;
                let pm = self.point(tm);
                let chord = pa.lerp(pb, f);
                if (pm - chord).length() > tol {
                    split = true;
                    break;
                }
            }
            if !split {
                let da = self.deriv(a);
                let db = self.deriv(b);
                if da.length() > 0.0 && db.length() > 0.0 && da.angle_between(db) > max_angle {
                    split = true;
                }
            }
        }
        if split {
            let m = 0.5 * (a + b);
            self.refine(a, m, tol, max_angle, depth + 1, out);
            self.refine(m, b, tol, max_angle, depth + 1, out);
        } else {
            out.push(b);
        }
    }
}

/// Classification of a surface into analytic kinds, used by downstream
/// algorithms and exporters.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum SurfaceKind {
    /// Plane.
    Plane,
    /// Circular cylinder.
    Cylinder,
    /// Circular cone.
    Cone,
    /// Sphere.
    Sphere,
    /// Torus.
    Torus,
    /// General surface of revolution.
    Revolution,
    /// General linear extrusion.
    Extrusion,
    /// Free-form NURBS.
    Nurbs,
}

/// A surface. Parameterizations:
///
/// * `Plane`: `(u, v)` are frame x/y coordinates.
/// * `Extrusion`: `curve(u) + v · dir`.
/// * `Revolution`: `profile(v)` rotated by angle `u` about the frame's z axis;
///   the profile lies in the frame's x/z half-plane (x ≥ 0).
/// * `Nurbs`: native parameters.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub(crate) enum Surface {
    Plane(Frame),
    Extrusion {
        curve: Curve3,
        dir: DVec3,
        /// A point and normal of the plane the base curve lies in, used to
        /// invert the `v` parameter.
        base_origin: DVec3,
        base_normal: DVec3,
    },
    Revolution {
        frame: Frame,
        profile: Curve3,
    },
    Nurbs(NurbsSurface),
}

impl Surface {
    /// Point at `(u, v)`.
    pub fn point(&self, u: f64, v: f64) -> DVec3 {
        match self {
            Surface::Plane(f) => f.origin + f.x * u + f.y * v,
            Surface::Extrusion { curve, dir, .. } => curve.point(u) + *dir * v,
            Surface::Revolution { frame, profile } => {
                let p = profile.point(v) - frame.origin;
                let rho = p.dot(frame.x);
                let h = p.dot(frame.z);
                let (s, c) = u.sin_cos();
                frame.origin + (frame.x * c + frame.y * s) * rho + frame.z * h
            }
            Surface::Nurbs(n) => n.point(u, v),
        }
    }

    /// Partial derivatives `(S_u, S_v)` at `(u, v)`.
    pub fn derivs(&self, u: f64, v: f64) -> (DVec3, DVec3) {
        match self {
            Surface::Plane(f) => (f.x, f.y),
            Surface::Extrusion { curve, dir, .. } => (curve.deriv(u), *dir),
            Surface::Revolution { frame, profile } => {
                let p = profile.point(v) - frame.origin;
                let dp = profile.deriv(v);
                let rho = p.dot(frame.x);
                let drho = dp.dot(frame.x);
                let dh = dp.dot(frame.z);
                let (s, c) = u.sin_cos();
                let radial = frame.x * c + frame.y * s;
                let tangential = frame.y * c - frame.x * s;
                (tangential * rho, radial * drho + frame.z * dh)
            }
            Surface::Nurbs(n) => {
                let d = n.derivatives(u, v, 1);
                (d[1][0], d[0][1])
            }
        }
    }

    /// Unit normal `S_u × S_v` at `(u, v)`, or `None` at a singular point.
    pub fn normal(&self, u: f64, v: f64) -> Option<DVec3> {
        let (su, sv) = self.derivs(u, v);
        let n = su.cross(sv);
        let scale = su.length() * sv.length();
        if scale <= 0.0 || n.length() <= 1e-10 * scale.max(1e-300) || n.length() < 1e-300 {
            None
        } else {
            Some(n.normalize())
        }
    }

    /// Parameters of the closest point to `p`. Exact for points on the
    /// surface; a good approximation for nearby points.
    pub fn project(&self, p: DVec3) -> (f64, f64) {
        match self {
            Surface::Plane(f) => {
                let d = p - f.origin;
                (d.dot(f.x), d.dot(f.y))
            }
            Surface::Extrusion {
                curve,
                dir,
                base_origin,
                base_normal,
            } => {
                let dn = dir.dot(*base_normal);
                let v = if dn.abs() > 1e-15 {
                    (p - *base_origin).dot(*base_normal) / dn
                } else {
                    0.0
                };
                let u = curve.project(p - *dir * v);
                (u, v)
            }
            Surface::Revolution { frame, profile } => {
                let d = p - frame.origin;
                let lx = d.dot(frame.x);
                let ly = d.dot(frame.y);
                let rho = (lx * lx + ly * ly).sqrt();
                let u = if rho > 1e-14 {
                    wrap_angle(ly.atan2(lx))
                } else {
                    0.0
                };
                let back = frame.origin + frame.x * rho + frame.z * d.dot(frame.z);
                (u, profile.project(back))
            }
            Surface::Nurbs(n) => n.closest_point(p).0,
        }
    }

    /// The analytic kind of this surface.
    pub fn kind(&self) -> SurfaceKind {
        match self {
            Surface::Plane(_) => SurfaceKind::Plane,
            Surface::Extrusion { curve, .. } => match curve {
                Curve3::Line { .. } => SurfaceKind::Plane,
                Curve3::Circle { .. } => SurfaceKind::Cylinder,
                Curve3::Nurbs(_) => SurfaceKind::Extrusion,
            },
            Surface::Revolution { frame, profile } => match profile {
                Curve3::Line { dir, .. } => {
                    let along = dir.normalize().dot(frame.z).abs();
                    if along > 1.0 - 1e-12 {
                        SurfaceKind::Cylinder
                    } else if along < 1e-12 {
                        SurfaceKind::Plane
                    } else {
                        SurfaceKind::Cone
                    }
                }
                Curve3::Circle { center, .. } => {
                    let d = *center - frame.origin;
                    if (d - frame.z * d.dot(frame.z)).length() < 1e-12 {
                        SurfaceKind::Sphere
                    } else {
                        SurfaceKind::Torus
                    }
                }
                Curve3::Nurbs(_) => SurfaceKind::Revolution,
            },
            Surface::Nurbs(_) => SurfaceKind::Nurbs,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    fn check_derivs(s: &Surface, u: f64, v: f64) {
        let h = 1e-6;
        let (su, sv) = s.derivs(u, v);
        let fu = (s.point(u + h, v) - s.point(u - h, v)) / (2.0 * h);
        let fv = (s.point(u, v + h) - s.point(u, v - h)) / (2.0 * h);
        assert!((su - fu).length() < 1e-6, "{su} vs {fu}");
        assert!((sv - fv).length() < 1e-6, "{sv} vs {fv}");
    }

    fn check_project(s: &Surface, u: f64, v: f64) {
        let p = s.point(u, v);
        let (pu, pv) = s.project(p);
        assert!((s.point(pu, pv) - p).length() < 1e-9);
    }

    #[test]
    fn circle_projection_and_derivative() {
        let c = Curve3::Circle {
            center: DVec3::new(1.0, 2.0, 3.0),
            x: DVec3::X,
            y: DVec3::Z,
            radius: 2.0,
        };
        for i in 0..20 {
            let t = -PI / 2.0 + i as f64 * 0.31;
            let p = c.point(t);
            assert!((c.project(p) - t).abs() < 1e-12);
            let h = 1e-6;
            let fd = (c.point(t + h) - c.point(t - h)) / (2.0 * h);
            assert!((fd - c.deriv(t)).length() < 1e-6);
        }
    }

    #[test]
    fn revolution_torus_and_cone() {
        let frame = Frame::WORLD;
        let torus = Surface::Revolution {
            frame,
            profile: Curve3::Circle {
                center: DVec3::new(3.0, 0.0, 0.0),
                x: DVec3::X,
                y: DVec3::Z,
                radius: 1.0,
            },
        };
        assert_eq!(torus.kind(), SurfaceKind::Torus);
        let cone = Surface::Revolution {
            frame,
            profile: Curve3::Line {
                origin: DVec3::new(0.0, 0.0, 2.0),
                dir: DVec3::new(1.0, 0.0, -2.0),
            },
        };
        assert_eq!(cone.kind(), SurfaceKind::Cone);
        for (u, v) in [(0.1, 0.2), (1.0, 2.5), (2.9, -0.4), (4.0, 1.0)] {
            check_derivs(&torus, u, v);
            check_project(&torus, u, v);
        }
        for (u, v) in [(0.1, 0.2), (1.0, 0.5), (2.9, 0.9)] {
            check_derivs(&cone, u, v);
            check_project(&cone, u, v);
        }
        assert!(cone.normal(0.3, 0.0).is_none());
    }

    #[test]
    fn extrusion_of_circle_is_cylinder() {
        let s = Surface::Extrusion {
            curve: Curve3::Circle {
                center: DVec3::ZERO,
                x: DVec3::X,
                y: DVec3::Y,
                radius: 0.5,
            },
            dir: DVec3::new(0.0, 0.0, 2.0),
            base_origin: DVec3::ZERO,
            base_normal: DVec3::Z,
        };
        assert_eq!(s.kind(), SurfaceKind::Cylinder);
        for (u, v) in [(0.1, 0.2), (1.0, 0.7), (2.9, 0.4)] {
            check_derivs(&s, u, v);
            check_project(&s, u, v);
            assert!(((s.point(u, v).truncate().length()) - 0.5).abs() < 1e-12);
        }
    }

    #[test]
    fn circle_discretization_meets_chord_tolerance() {
        let c = Curve3::Circle {
            center: DVec3::ZERO,
            x: DVec3::X,
            y: DVec3::Y,
            radius: 1.0,
        };
        let ts = c.discretize(0.0, PI, 1e-3, 1.0);
        for w in ts.windows(2) {
            let mid = c.point(0.5 * (w[0] + w[1]));
            let chord = 0.5 * (c.point(w[0]) + c.point(w[1]));
            assert!((mid - chord).length() <= 1e-3 + 1e-12);
        }
    }
}
