//! Curve and surface geometry carried by B-rep edges and faces.
//!
//! Analytic types are kept exact wherever possible (lines, circles, planes,
//! surfaces of revolution and linear extrusion), with NURBS for free-form
//! geometry. Every type supports evaluation, first derivatives, projection
//! (closest point) and affine transformation.

use glam::{DMat4, DVec3};
use opendatum_nurbs::{NurbsCurve, NurbsSurface};
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

/// A piece of a curve this fraction of the chord tolerance long that turns
/// further than the angle tolerance is a kink ([`Curve3::discretize`]).
const KINK: f64 = 1e-3;

impl Curve3 {
    /// Size of the coordinates the curve's points are computed from (its
    /// origin, its centre plus radius, or its control points' bounds):
    /// together with the size of the points themselves, this sets their
    /// rounding errors, about this times the machine epsilon.
    pub fn magnitude(&self) -> f64 {
        match self {
            Curve3::Line { origin, .. } => origin.abs().max_element(),
            Curve3::Circle { center, radius, .. } => center.abs().max_element() + radius,
            Curve3::Nurbs(n) => n
                .control_points()
                .iter()
                .map(|p| p.abs().max_element())
                .fold(0.0, f64::max),
        }
    }

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

    /// The point and first derivative at `t`, as [`Curve3::point`] and
    /// [`Curve3::deriv`] give them, from one evaluation of a NURBS.
    pub fn point_deriv(&self, t: f64) -> (DVec3, DVec3) {
        match self {
            Curve3::Line { origin, dir } => (*origin + *dir * t, *dir),
            Curve3::Circle {
                center,
                x,
                y,
                radius,
            } => {
                let (s, c) = t.sin_cos();
                (
                    *center + (*x * c + *y * s) * *radius,
                    (*y * c - *x * s) * *radius,
                )
            }
            Curve3::Nurbs(n) => {
                let d = n.derivatives(t, 1);
                (d[0], d[1])
            }
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
    /// `max_angle` per segment, in at most `limit` segments.
    ///
    /// A circle's steps meet `r (1 − cos(Δ/2)) ≤ chord_tol` exactly. A
    /// NURBS curve's segments are bisected until the curve's points at a
    /// quarter, half and three quarters of each are within `chord_tol` of
    /// its chord (an estimate, not a bound), and until its ends' tangents
    /// differ by at most `max_angle`, also on pieces shorter than
    /// `chord_tol`: a face's triangles beside a short piece can be long,
    /// along its rulings, and take the turn across them. A piece a
    /// thousandth of `chord_tol` long that still turns further is taken as
    /// a kink, where the tangent jumps and no bisection meets the angle.
    /// Fails if that takes more than `limit` segments, instead of returning
    /// a coarser polyline, and if a tolerance is not positive (or is NaN).
    pub fn discretize(
        &self,
        t0: f64,
        t1: f64,
        chord_tol: f64,
        max_angle: f64,
        limit: usize,
    ) -> Result<Vec<f64>, String> {
        let too_many = || format!("an edge needs more than {limit} segments at this tolerance");
        // A NaN or negative tolerance would pass every comparison below
        // the wrong way and come back as a coarse polyline (found by
        // review: a half circle as one chord for an angle of -1).
        if !(chord_tol > 0.0 && max_angle > 0.0) {
            return Err(format!(
                "tolerances must be positive (chord {chord_tol}, angle {max_angle})"
            ));
        }
        match self {
            Curve3::Line { .. } => Ok(vec![t0, t1]),
            Curve3::Circle { radius, .. } => {
                let sweep = t1 - t0;
                // 1 − cos(Δ/2) = 2 sin²(Δ/4), without the cancellation of
                // the cosine's form when the tolerance is a tiny part of
                // the radius.
                let by_chord = if chord_tol < *radius {
                    4.0 * (0.5 * chord_tol / radius).sqrt().asin()
                } else {
                    std::f64::consts::PI
                };
                let step = by_chord.min(max_angle);
                let n = (sweep.abs() / step).ceil();
                // NaN too (a NaN parameter) is refused.
                if n.is_nan() || n > limit as f64 {
                    return Err(too_many());
                }
                let n = (n as usize).max(1);
                Ok((0..=n).map(|i| t0 + sweep * i as f64 / n as f64).collect())
            }
            Curve3::Nurbs(_) => {
                // Start from a uniform split, then bisect adaptively.
                let n0 = 8;
                let mut out = vec![t0];
                for i in 0..n0 {
                    let a = t0 + (t1 - t0) * i as f64 / n0 as f64;
                    let b = t0 + (t1 - t0) * (i + 1) as f64 / n0 as f64;
                    self.refine(a, b, chord_tol, max_angle, limit, &mut out)
                        .ok_or_else(too_many)?;
                }
                Ok(out)
            }
        }
    }

    /// Bisects `[a, b]` as [`Curve3::discretize`] says, appending the ends
    /// of its pieces after `a` to `out`; `None` past `limit` of them (or
    /// sixty levels of halving, where the parameters reach their
    /// rounding).
    fn refine(
        &self,
        a: f64,
        b: f64,
        tol: f64,
        max_angle: f64,
        limit: usize,
        out: &mut Vec<f64>,
    ) -> Option<()> {
        // Depth first, by an explicit stack of the pieces still to test
        // (the later half on top of the earlier).
        let mut stack = vec![(b, 0u32)];
        let mut at = a;
        while let Some((end, depth)) = stack.pop() {
            let (pa, pb) = (self.point(at), self.point(end));
            let off = [0.25, 0.5, 0.75].into_iter().any(|f| {
                let pm = self.point(at + (end - at) * f);
                (pm - pa.lerp(pb, f)).length() > tol
            });
            let turns = || {
                let (da, db) = (self.deriv(at), self.deriv(end));
                (pb - pa).length() > KINK * tol
                    && da.length() > 0.0
                    && db.length() > 0.0
                    && da.angle_between(db) > max_angle
            };
            if off || turns() {
                if depth == 60 {
                    return None;
                }
                stack.push((end, depth + 1));
                stack.push((0.5 * (at + end), depth + 1));
            } else {
                out.push(end);
                if out.len() > limit {
                    return None;
                }
                at = end;
            }
        }
        Some(())
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

/// The unit normal `S_u × S_v` from a surface's partial derivatives, or
/// `None` where they are parallel or vanish (a singular point).
pub(crate) fn unit_normal(su: DVec3, sv: DVec3) -> Option<DVec3> {
    let n = su.cross(sv);
    let scale = su.length() * sv.length();
    if scale <= 0.0 || n.length() <= 1e-10 * scale.max(1e-300) || n.length() < 1e-300 {
        None
    } else {
        Some(n.normalize())
    }
}

impl Surface {
    /// Size of the coordinates the surface's points are computed from, as
    /// [`Curve3::magnitude`].
    pub fn magnitude(&self) -> f64 {
        match self {
            Surface::Plane(f) => f.origin.abs().max_element(),
            Surface::Extrusion { curve, .. } => curve.magnitude(),
            Surface::Revolution { frame, profile } => {
                frame.origin.abs().max_element() + profile.magnitude()
            }
            Surface::Nurbs(n) => n
                .control_points()
                .iter()
                .flatten()
                .map(|p| p.abs().max_element())
                .fold(0.0, f64::max),
        }
    }

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

    /// The point at `(u, v)` and the partial derivatives `(S_u, S_v)`
    /// there, as [`Surface::point`] and [`Surface::derivs`] give them,
    /// sharing their work.
    pub fn point_derivs(&self, u: f64, v: f64) -> (DVec3, DVec3, DVec3) {
        match self {
            Surface::Plane(f) => (f.origin + f.x * u + f.y * v, f.x, f.y),
            Surface::Extrusion { curve, dir, .. } => {
                let (p, d) = curve.point_deriv(u);
                (p + *dir * v, d, *dir)
            }
            Surface::Revolution { frame, profile } => {
                let (p, dp) = profile.point_deriv(v);
                let p = p - frame.origin;
                let (rho, h) = (p.dot(frame.x), p.dot(frame.z));
                let (drho, dh) = (dp.dot(frame.x), dp.dot(frame.z));
                let (s, c) = u.sin_cos();
                let radial = frame.x * c + frame.y * s;
                let tangential = frame.y * c - frame.x * s;
                (
                    frame.origin + radial * rho + frame.z * h,
                    tangential * rho,
                    radial * drho + frame.z * dh,
                )
            }
            Surface::Nurbs(n) => {
                let d = n.derivatives(u, v, 1);
                (d[0][0], d[1][0], d[0][1])
            }
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
        unit_normal(su, sv)
    }

    /// The unit normal at `(u, v)` as [`Surface::normal`], or at a point
    /// of a surface of revolution's axis where the surface is smooth (a
    /// sphere's pole), where `S_u` vanishes: there the normal is along the
    /// axis, pointing as `S_u × S_v` does beside the point. `None` where
    /// the profile meets the axis at a slant (a cone's apex, where the
    /// surface has no tangent plane).
    pub fn normal_or_limit(&self, u: f64, v: f64) -> Option<DVec3> {
        if let Some(n) = self.normal(u, v) {
            return Some(n);
        }
        let Surface::Revolution { frame, profile } = self else {
            return None;
        };
        let dp = profile.deriv(v);
        let len = dp.length();
        if len.is_nan() || len == 0.0 || dp.dot(frame.z).abs() > 1e-9 * len {
            return None;
        }
        // A step along the profile far above the rounding of its points and
        // far below its size, to the side of the axis it lies on (x ≥ 0 in
        // the frame); beside a smooth point the normal keeps its side.
        let step = 1e-8 * profile.magnitude().max(len) / len;
        [step, -step].into_iter().find_map(|dv| {
            let w = v + dv;
            let rho = (profile.point(w) - frame.origin).dot(frame.x);
            let n = self.normal(u, w).filter(|_| rho > 0.0)?;
            Some(frame.z * n.dot(frame.z).signum())
        })
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

    /// For a surface of kind [`SurfaceKind::Plane`], a point of the plane
    /// and the unit direction of `S_u × S_v` (the same at every regular
    /// point); `None` for every other kind.
    pub fn plane(&self) -> Option<(DVec3, DVec3)> {
        if self.kind() != SurfaceKind::Plane {
            return None;
        }
        match self {
            Surface::Plane(f) => Some((f.origin, f.x.cross(f.y).normalize())),
            Surface::Extrusion {
                curve: Curve3::Line { origin, dir: along },
                dir,
                ..
            } => {
                let n = along.cross(*dir);
                (n.length() > 0.0).then(|| (*origin, n.normalize()))
            }
            // The profile is a line perpendicular to the axis in the x/z
            // half-plane (ρ ≥ 0): S_u × S_v = −ρ ρ' z, with ρ' = dir · x.
            Surface::Revolution {
                frame,
                profile: Curve3::Line { origin, dir },
            } => {
                let drho = dir.dot(frame.x);
                (drho != 0.0).then(|| (*origin, -frame.z * drho.signum()))
            }
            _ => None,
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
                    // Sine and cosine of the line's angle to the axis.
                    let dir = dir.normalize();
                    let (across, along) = (dir.cross(frame.z).length(), dir.dot(frame.z).abs());
                    if across < 1e-12 {
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
        let (p, du, dv) = s.point_derivs(u, v);
        let near = |a: DVec3, b: DVec3| (a - b).length() <= 1e-12 * (1.0 + b.length());
        assert!(near(p, s.point(u, v)) && near(du, su) && near(dv, sv));
    }

    fn check_project(s: &Surface, u: f64, v: f64) {
        let p = s.point(u, v);
        let (pu, pv) = s.project(p);
        assert!((s.point(pu, pv) - p).length() < 1e-9);
    }

    #[test]
    fn planes_have_the_normal_of_their_derivatives() {
        let check = |s: &Surface, samples: &[(f64, f64)]| {
            let (o, n) = s.plane().expect("planar");
            for &(u, v) in samples {
                let p = s.point(u, v);
                assert!((p - o).dot(n).abs() < 1e-12, "{p} off the plane");
                let m = s.normal(u, v).unwrap();
                assert!((m - n).length() < 1e-12, "{m} vs {n}");
            }
        };
        let frame = Frame::from_normal(
            DVec3::new(1.0, -2.0, 0.5),
            DVec3::new(1.0, 2.0, 3.0),
            DVec3::X,
        );
        check(&Surface::Plane(frame), &[(0.0, 0.0), (1.5, -2.0)]);
        let flipped = Frame {
            y: -frame.y,
            ..frame
        };
        check(&Surface::Plane(flipped), &[(0.3, 0.7)]);
        for along in [DVec3::new(1.0, 1.0, 0.0), DVec3::new(-2.0, 0.5, 0.0)] {
            let s = Surface::Extrusion {
                curve: Curve3::Line {
                    origin: DVec3::new(0.5, 0.0, 0.0),
                    dir: along,
                },
                dir: DVec3::new(0.0, 0.0, 2.0),
                base_origin: DVec3::ZERO,
                base_normal: DVec3::Z,
            };
            check(&s, &[(0.0, 0.0), (0.4, 1.0), (1.0, -0.5)]);
        }
        // Discs and annuli of a revolve, profile running out or in.
        for (start, dir) in [(0.0, 2.0), (3.0, -1.5), (1.0, 1.0)] {
            let s = Surface::Revolution {
                frame,
                profile: Curve3::Line {
                    origin: frame.origin + frame.x * start + frame.z * 0.75,
                    dir: frame.x * dir,
                },
            };
            check(&s, &[(0.2, 0.5), (2.5, 0.9), (5.0, 0.1)]);
        }
        let cylinder = Surface::Revolution {
            frame,
            profile: Curve3::Line {
                origin: frame.origin + frame.x,
                dir: frame.z,
            },
        };
        assert_eq!(cylinder.plane(), None);
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
        let ts = c.discretize(0.0, PI, 1e-3, 1.0, 1000).expect("samples");
        for w in ts.windows(2) {
            let mid = c.point(0.5 * (w[0] + w[1]));
            let chord = 0.5 * (c.point(w[0]) + c.point(w[1]));
            assert!((mid - chord).length() <= 1e-3 + 1e-12);
        }
    }

    /// A kilometre circle at a micrometre: steps of 8.9e-5 rad, where a
    /// floor of 1e-4 rad left chords 1.25 tolerances off (found by
    /// review). Far finer than the limit allows is refused, not coarsened.
    #[test]
    fn large_circles_meet_fine_tolerances_or_are_refused() {
        let c = Curve3::Circle {
            center: DVec3::ZERO,
            x: DVec3::X,
            y: DVec3::Y,
            radius: 1e3,
        };
        let ts = c.discretize(0.0, PI, 1e-6, 0.3, 100_000).expect("samples");
        assert_eq!(ts.len(), 35_126);
        for w in ts.windows(2) {
            let mid = c.point(0.5 * (w[0] + w[1]));
            let chord = 0.5 * (c.point(w[0]) + c.point(w[1]));
            assert!((mid - chord).length() <= 1e-6 * (1.0 + 1e-6));
        }
        assert!(c.discretize(0.0, PI, 1e-12, 0.3, 100_000).is_err());
        assert!(c.discretize(0.0, PI, 1e-6, 1e-9, 100_000).is_err());
    }

    /// Tolerances that are not positive, or NaN, are refused on every
    /// kind of curve; they had come back as coarse polylines (found by
    /// review: a half circle as one chord for an angle of -1), and a NaN
    /// parameter is refused too.
    #[test]
    fn discretizing_refuses_tolerances_that_are_not_positive() {
        let circle = Curve3::Circle {
            center: DVec3::ZERO,
            x: DVec3::X,
            y: DVec3::Y,
            radius: 1.0,
        };
        let spline = Curve3::Nurbs(
            NurbsCurve::new(
                2,
                vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
                vec![DVec3::ZERO, DVec3::new(1.0, 1.0, 0.0), DVec3::X * 2.0],
                None,
            )
            .expect("a parabola"),
        );
        let line = Curve3::Line {
            origin: DVec3::ZERO,
            dir: DVec3::X,
        };
        for c in [&circle, &spline, &line] {
            for (chord, angle) in [
                (f64::NAN, 0.3),
                (1e-3, f64::NAN),
                (1e-3, -1.0),
                (-1e-3, 0.3),
                (0.0, 0.3),
                (1e-3, 0.0),
            ] {
                assert!(
                    c.discretize(0.0, 1.0, chord, angle, 1000).is_err(),
                    "{chord} {angle}"
                );
            }
            assert!(c.discretize(0.0, 1.0, 1e-3, 0.3, 1000).is_ok());
        }
        assert!(circle.discretize(0.0, f64::NAN, 1e-3, 0.3, 1000).is_err());
    }

    /// A NURBS curve's samples turn by at most the angle tolerance also
    /// where its bends are tighter than the chord tolerance (found by
    /// review: pieces under the chord were let turn further, 0.835 rad at
    /// 0.5 on this zig-zag, and its extrusion meshed to 19,784 triangles
    /// where 216 do at 1e-2 and 0.5 rad, 92 times as many). A kink, where
    /// the tangent jumps, is closed in on to a thousandth of the tolerance
    /// and no further.
    #[test]
    fn nurbs_samples_turn_within_the_angle_at_tight_bends() {
        let n = 10;
        let points: Vec<DVec3> = (0..=n)
            .map(|i| {
                let y = match i {
                    0 => 0.5,
                    _ if i == n => 0.5,
                    _ if i % 2 == 1 => 1.5,
                    _ => 0.05,
                };
                DVec3::new(i as f64 / n as f64, y, 0.0)
            })
            .collect();
        let mut knots = vec![0.0; 4];
        knots.extend((1..n - 2).map(|i| i as f64 / (n - 2) as f64));
        knots.extend([1.0; 4]);
        let zig = Curve3::Nurbs(NurbsCurve::new(3, knots, points, None).expect("a spline"));
        for (tol, max_angle) in [(1e-2, 0.5), (1e-2, 0.1), (1e-4, 0.3)] {
            let ts = zig
                .discretize(0.0, 1.0, tol, max_angle, 100_000)
                .expect("samples");
            let turn = ts
                .windows(2)
                .map(|w| zig.deriv(w[0]).angle_between(zig.deriv(w[1])))
                .fold(0.0, f64::max);
            assert!(turn <= max_angle, "{tol} {max_angle}: {turn}");
        }
        let kink = Curve3::Nurbs(
            NurbsCurve::new(
                1,
                vec![0.0, 0.0, 0.3, 1.0, 1.0],
                vec![DVec3::ZERO, DVec3::X, DVec3::new(1.0, 1.0, 0.0)],
                None,
            )
            .expect("a polyline"),
        );
        // Surfaces swept from the spline evaluate their points and
        // derivatives together as apart.
        for s in [
            Surface::Extrusion {
                curve: zig.clone(),
                dir: DVec3::Z,
                base_origin: DVec3::ZERO,
                base_normal: DVec3::Z,
            },
            Surface::Revolution {
                frame: Frame::WORLD,
                profile: zig.clone(),
            },
        ] {
            for (u, v) in [(0.13, 0.4), (0.71, 0.9), (0.5, 0.05)] {
                check_derivs(&s, u, v);
            }
        }
        let ts = kink.discretize(0.0, 1.0, 1e-2, 0.5, 1000).expect("samples");
        let near = ts
            .iter()
            .map(|&t| (kink.point(t) - DVec3::X).length())
            .fold(f64::INFINITY, f64::min);
        assert!(ts.len() < 50, "{}", ts.len());
        assert!(near < 2.0 * KINK * 1e-2, "{near}");
    }
}
