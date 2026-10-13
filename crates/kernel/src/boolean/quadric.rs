//! Planes, circular cylinders, cones and spheres: the surfaces step 2 of
//! the booleans handles (design note 0012 §10.1), and how an edge's
//! distance to one varies along it (§10.2).
//!
//! Each surface has a signed distance, positive outside (away from the
//! axis or centre; the side a plane's normal points to), which decides
//! every sign, and an implicit polynomial `G` of degree one or two with
//! the same sign wherever the distance is small. Along a rational Bézier
//! piece of an edge `G` is a ratio of polynomials, and the roots of its
//! derivative, found in Bernstein form, split the edge into stretches on
//! which `G` is monotone: on each, the edge crosses the surface at most
//! once, exactly when the distances at the stretch's ends have strictly
//! opposite signs.

use glam::{DVec3, DVec4};

use super::bernstein::Bernstein;
use crate::geom::{Curve3, Surface, SurfaceKind};
use crate::KernelError;

/// Relative tolerance for recognising one surface built twice (rounding
/// only), as step 1's test for planes: sines of angles, and offsets and
/// radii relative to the size of the coordinates.
const SAME: f64 = 1e-10;

/// A surface of step 2.
#[derive(Clone, Copy, Debug, PartialEq)]
pub(super) enum Quadric {
    /// The points `p` with `(p − point) · normal = 0`, `normal` a unit
    /// vector.
    Plane { point: DVec3, normal: DVec3 },
    /// The points at `radius` from the line through `point` along the unit
    /// vector `axis`.
    Cylinder {
        point: DVec3,
        axis: DVec3,
        radius: f64,
    },
    /// One nappe of a circular cone: the points whose distance from the
    /// line through `point` along the unit vector `axis` is `radius + slope
    /// · h` where that is not negative, `h` their position along the axis
    /// from `point`. `slope = tan α > 0` for the half-angle `α`; the apex
    /// is at `h = −radius / slope`. The point is kept near the faces on the
    /// cone, so that far apexes (slight tapers) cost no precision.
    Cone {
        point: DVec3,
        axis: DVec3,
        radius: f64,
        slope: f64,
    },
    /// The points at `radius` from `center`.
    Sphere { center: DVec3, radius: f64 },
}

/// Where an edge is split so that a surface's implicit function is
/// monotone between consecutive splits (and the edge's ends).
#[derive(Clone, Copy, Debug, PartialEq)]
pub(super) struct Split {
    /// The curve parameter.
    pub t: f64,
    /// True at a turning point of the implicit function (where the edge
    /// may touch the surface); false at a joint between two Bézier pieces
    /// of the edge.
    pub turn: bool,
}

/// One rational Bézier piece of an edge, with the map from its parameter
/// `s ∈ [0, 1]` to the curve's.
struct Piece {
    /// Homogeneous control points `(w x, w y, w z, w)`.
    pw: Vec<DVec4>,
    map: Map,
}

#[derive(Clone, Copy)]
enum Map {
    /// `t = a + (b − a) s`.
    Linear(f64, f64),
    /// A circular arc of the rational quadratic form with weights `1,
    /// cos(δ/2), 1` about its middle angle `mid`: `tan((t − mid)/2) =
    /// tan(δ/4) (2s − 1)`, with `half = tan(δ/4)`.
    Arc { mid: f64, half: f64 },
}

impl Map {
    fn t(&self, s: f64) -> f64 {
        match *self {
            Map::Linear(a, b) => a + (b - a) * s,
            Map::Arc { mid, half } => mid + 2.0 * (half * (2.0 * s - 1.0)).atan(),
        }
    }
}

/// The implicit function along a piece: `G(C(s)) = n(s) / w(s)^k`.
struct Along {
    n: Bernstein,
    w: Bernstein,
    k: f64,
    /// Size of the terms `n` is computed from, for its rounding.
    scale: f64,
}

/// Bisection steps after which a crossing is taken as found (a parameter
/// interval of `2⁻²⁰⁰` of the stretch; floats run out long before).
const BISECTIONS: usize = 200;

impl Quadric {
    /// The surface of a face, exactly from the kernel's surface record, or
    /// `None` for surfaces step 2 does not handle (tori, general surfaces
    /// of revolution and extrusion, NURBS, and an extrusion of a circle
    /// along a direction other than its axis, an elliptic cylinder).
    pub fn of(surface: &Surface) -> Option<Quadric> {
        if let Some((point, normal)) = surface.plane() {
            return Some(Quadric::Plane { point, normal });
        }
        match surface {
            Surface::Extrusion {
                curve:
                    Curve3::Circle {
                        center,
                        x,
                        y,
                        radius,
                    },
                dir,
                ..
            } => {
                let axis = dir.normalize();
                (axis.cross(x.cross(*y)).length() <= 1e-12).then_some(Quadric::Cylinder {
                    point: *center,
                    axis,
                    radius: *radius,
                })
            }
            Surface::Revolution {
                frame,
                profile: Curve3::Line { origin, dir },
            } => {
                // The profile lies in the frame's x/z half-plane (x ≥ 0).
                let w = *origin - frame.origin;
                let h = w.dot(frame.z);
                let rho = w.dot(frame.x);
                match surface.kind() {
                    SurfaceKind::Cylinder => Some(Quadric::Cylinder {
                        point: frame.origin,
                        axis: frame.z,
                        radius: (w - frame.z * h).length(),
                    }),
                    SurfaceKind::Cone => {
                        let slope = dir.dot(frame.x) / dir.dot(frame.z);
                        let (axis, slope) = if slope > 0.0 {
                            (frame.z, slope)
                        } else {
                            (-frame.z, -slope)
                        };
                        // A profile starting on the axis (a pointed cone's
                        // apex) has ρ of rounding size, either sign.
                        let on_axis = -1e-12 * (1.0 + w.length());
                        (slope.is_finite() && rho >= on_axis).then_some(Quadric::Cone {
                            point: frame.origin + frame.z * h,
                            axis,
                            radius: rho.max(0.0),
                            slope,
                        })
                    }
                    _ => None,
                }
            }
            Surface::Revolution {
                profile: Curve3::Circle { center, radius, .. },
                ..
            } if surface.kind() == SurfaceKind::Sphere => Some(Quadric::Sphere {
                center: *center,
                radius: *radius,
            }),
            _ => None,
        }
    }

    /// The signed distance from `p`: positive outside (away from the axis
    /// or centre; the side the normal points to). Exact up to rounding;
    /// behind a cone's apex, the distance to the apex.
    pub fn distance(&self, p: DVec3) -> f64 {
        match *self {
            Quadric::Plane { point, normal } => (p - point).dot(normal),
            Quadric::Cylinder {
                point,
                axis,
                radius,
            } => {
                let w = p - point;
                (w - axis * w.dot(axis)).length() - radius
            }
            Quadric::Sphere { center, radius } => (p - center).length() - radius,
            Quadric::Cone {
                point,
                axis,
                radius,
                slope,
            } => {
                let w = p - point;
                let h = w.dot(axis);
                let rho = (w - axis * h).length();
                let cos = 1.0 / (1.0 + slope * slope).sqrt();
                let sin = slope * cos;
                // Position along the axis from the apex, and of the foot
                // of `p` on the generator line in its half-plane, measured
                // from the apex: the foot is on the nappe when this is not
                // negative.
                let z = h + radius / slope;
                if rho * sin + z * cos >= 0.0 {
                    (rho - radius - slope * h) * cos
                } else {
                    (rho * rho + z * z).sqrt()
                }
            }
        }
    }

    /// The unit direction in which the signed distance grows at `p`: the
    /// normal of a plane, away from the axis or the centre, and behind a
    /// cone's apex away from the apex. Zero on an axis or at a centre.
    pub fn gradient(&self, p: DVec3) -> DVec3 {
        match *self {
            Quadric::Plane { normal, .. } => normal,
            Quadric::Cylinder { point, axis, .. } => {
                let w = p - point;
                (w - axis * w.dot(axis)).normalize_or_zero()
            }
            Quadric::Sphere { center, .. } => (p - center).normalize_or_zero(),
            Quadric::Cone {
                point,
                axis,
                radius,
                slope,
            } => {
                let w = p - point;
                let h = w.dot(axis);
                let across = w - axis * h;
                let cos = 1.0 / (1.0 + slope * slope).sqrt();
                let sin = slope * cos;
                let z = h + radius / slope;
                if across.length() * sin + z * cos >= 0.0 {
                    across.normalize_or_zero() * cos - axis * sin
                } else {
                    (w + axis * (radius / slope)).normalize_or_zero()
                }
            }
        }
    }

    /// The point of the surface nearest `p` (`None` where there is no one
    /// nearest: on a cylinder's or a cone's axis, at a sphere's centre).
    /// On a cone, the foot on the ruling in `p`'s half-plane, or the apex
    /// for a point behind it.
    pub fn foot(&self, p: DVec3) -> Option<DVec3> {
        match *self {
            Quadric::Plane { point, normal } => Some(p - normal * (p - point).dot(normal)),
            Quadric::Cylinder {
                point,
                axis,
                radius,
            } => {
                let w = p - point;
                let h = w.dot(axis);
                let out = (w - axis * h).try_normalize()?;
                Some(point + axis * h + out * radius)
            }
            Quadric::Sphere { center, radius } => {
                Some(center + (p - center).try_normalize()? * radius)
            }
            Quadric::Cone {
                point,
                axis,
                radius,
                slope,
            } => {
                let w = p - point;
                let h = w.dot(axis);
                let across = w - axis * h;
                let out = across.try_normalize()?;
                let cos = 1.0 / (1.0 + slope * slope).sqrt();
                let sin = slope * cos;
                // Slant position along the ruling from the reference
                // circle, no further back than the apex.
                let sigma = ((across.length() - radius) * sin + h * cos).max(-radius / sin);
                Some(point + out * (radius + sigma * sin) + axis * (sigma * cos))
            }
        }
    }

    /// Size of the coordinates the surface is given by, for its rounding.
    pub fn magnitude(&self) -> f64 {
        match *self {
            Quadric::Plane { point, .. } => point.abs().max_element(),
            Quadric::Cylinder { point, radius, .. } | Quadric::Cone { point, radius, .. } => {
                point.abs().max_element() + radius
            }
            Quadric::Sphere { center, radius } => center.abs().max_element() + radius,
        }
    }

    /// True if `self` and `other` are one surface up to rounding: two
    /// records of a surface the kernel built once (the two halves of a
    /// cylinder, ADR 0008), as step 1 groups faces on one plane. Planes
    /// match facing either way.
    pub fn same(&self, other: &Quadric) -> bool {
        let scale = 1.0 + self.magnitude().max(other.magnitude());
        let near = |d: f64| d.abs() <= SAME * scale;
        let parallel = |a: DVec3, b: DVec3| a.cross(b).length() <= SAME;
        // The distance of `q` from the line through `p` along `axis`.
        let off_axis = |p: DVec3, axis: DVec3, q: DVec3| {
            let w = q - p;
            (w - axis * w.dot(axis)).length()
        };
        match (*self, *other) {
            (
                Quadric::Plane { point, normal },
                Quadric::Plane {
                    point: q,
                    normal: m,
                },
            ) => parallel(normal, m) && near((q - point).dot(normal)),
            (
                Quadric::Cylinder {
                    point,
                    axis,
                    radius,
                },
                Quadric::Cylinder {
                    point: q,
                    axis: v,
                    radius: s,
                },
            ) => parallel(axis, v) && near(off_axis(point, axis, q)) && near(radius - s),
            (
                Quadric::Cone {
                    point,
                    axis,
                    radius,
                    slope,
                },
                Quadric::Cone {
                    point: q,
                    axis: v,
                    radius: s,
                    slope: k,
                },
            ) => {
                // The other's radius at this one's point.
                let at = s + k * (point - q).dot(v);
                axis.dot(v) > 0.0
                    && parallel(axis, v)
                    && near(off_axis(point, axis, q))
                    && near(radius - at)
                    && (slope - k).abs() <= SAME * (1.0 + slope)
            }
            (
                Quadric::Sphere { center, radius },
                Quadric::Sphere {
                    center: c,
                    radius: s,
                },
            ) => near((center - c).length()) && near(radius - s),
            _ => false,
        }
    }

    /// The same surface with its reference point moved near `p` (a plane's
    /// foot of `p`, an axis's point nearest `p`), so that the implicit
    /// function near `p` is computed from small numbers.
    fn anchored(&self, p: DVec3) -> Quadric {
        match *self {
            Quadric::Plane { point, normal } => Quadric::Plane {
                point: p - normal * (p - point).dot(normal),
                normal,
            },
            Quadric::Cylinder {
                point,
                axis,
                radius,
            } => Quadric::Cylinder {
                point: point + axis * (p - point).dot(axis),
                axis,
                radius,
            },
            Quadric::Cone {
                point,
                axis,
                radius,
                slope,
            } => {
                let h = (p - point).dot(axis);
                Quadric::Cone {
                    point: point + axis * h,
                    axis,
                    radius: radius + slope * h,
                    slope,
                }
            }
            sphere @ Quadric::Sphere { .. } => sphere,
        }
    }

    /// The implicit function along a rational Bézier piece, `None` if its
    /// degree is beyond the polynomials of `bernstein`.
    fn along(&self, pw: &[DVec4]) -> Option<Along> {
        let first = pw[0].truncate() / pw[0].w;
        let q = self.anchored(first);
        let reference = match q {
            Quadric::Plane { point, .. }
            | Quadric::Cylinder { point, .. }
            | Quadric::Cone { point, .. } => point,
            Quadric::Sphere { center, .. } => center,
        };
        let weights: Vec<f64> = pw.iter().map(|p| p.w).collect();
        let w = Bernstein::new(&weights)?;
        let wmax = weights.iter().fold(0.0f64, |m, x| m.max(x.abs()));
        let rel: Vec<DVec3> = pw.iter().map(|p| p.truncate() - reference * p.w).collect();
        let reach = rel.iter().fold(0.0f64, |m, x| m.max(x.length()));
        let coords = |f: &dyn Fn(DVec3) -> f64| {
            let c: Vec<f64> = rel.iter().map(|x| f(*x)).collect();
            Bernstein::new(&c)
        };
        // |X|², the squared length of the homogeneous position.
        let squared = || -> Option<Bernstein> {
            let x = coords(&|v| v.x)?;
            let y = coords(&|v| v.y)?;
            let z = coords(&|v| v.z)?;
            Some(
                x.mul(&x)?
                    .add_scaled(1.0, &y.mul(&y)?)
                    .add_scaled(1.0, &z.mul(&z)?),
            )
        };
        Some(match q {
            Quadric::Plane { normal, .. } => Along {
                n: coords(&|v| v.dot(normal))?,
                w,
                k: 1.0,
                scale: reach,
            },
            Quadric::Cylinder { axis, radius, .. } => {
                let h = coords(&|v| v.dot(axis))?;
                let n = squared()?
                    .add_scaled(-1.0, &h.mul(&h)?)
                    .add_scaled(-radius * radius, &w.mul(&w)?);
                Along {
                    n,
                    w,
                    k: 2.0,
                    scale: (reach + radius * wmax).powi(2),
                }
            }
            Quadric::Cone {
                axis,
                radius,
                slope,
                ..
            } => {
                let h = coords(&|v| v.dot(axis))?;
                let r = w.scaled(radius).add_scaled(slope, &h);
                let n = squared()?
                    .add_scaled(-1.0, &h.mul(&h)?)
                    .add_scaled(-1.0, &r.mul(&r)?);
                Along {
                    n,
                    w,
                    k: 2.0,
                    scale: ((reach + radius.abs() * wmax) * (1.0 + slope)).powi(2),
                }
            }
            Quadric::Sphere { radius, .. } => Along {
                n: squared()?.add_scaled(-radius * radius, &w.mul(&w)?),
                w,
                k: 2.0,
                scale: (reach + radius * wmax).powi(2),
            },
        })
    }

    /// The splits of the edge `curve` on `[t0, t1]` strictly between its
    /// ends, ordered from `t0` to `t1`: the turning points of the implicit
    /// function along it, and the joints of its rational Bézier pieces
    /// (quarter turns of a circle, the knot spans of a NURBS curve).
    /// Between consecutive splits the implicit function is monotone, so
    /// the edge crosses the surface there at most once (behind a cone's
    /// apex `G` also vanishes on the other nappe, where the distance stays
    /// positive; a monotone `G` still crosses zero once at most, so the
    /// distance changes sign at most once).
    pub fn splits(&self, curve: &Curve3, t0: f64, t1: f64) -> Result<Vec<Split>, KernelError> {
        let (lo, hi) = if t0 <= t1 { (t0, t1) } else { (t1, t0) };
        let unsupported = || {
            KernelError::Unsupported(
                "booleans with an edge of degree above five (design note 0012 §10.2)".into(),
            )
        };
        let pieces = pieces(curve, lo, hi).ok_or_else(unsupported)?;
        let mut out: Vec<Split> = Vec::new();
        let mut push = |s: Split| {
            if s.t > lo && s.t < hi && out.last().is_none_or(|l| s.t > l.t) {
                out.push(s);
            }
        };
        for (i, piece) in pieces.iter().enumerate() {
            if i > 0 {
                push(Split {
                    t: piece.map.t(0.0),
                    turn: false,
                });
            }
            let along = self.along(&piece.pw).ok_or_else(unsupported)?;
            // d/ds (n / wᵏ) = (n' w − k n w') / wᵏ⁺¹.
            let d = along
                .n
                .deriv()
                .mul(&along.w)
                .and_then(|a| Some(a.add_scaled(-along.k, &along.n.mul(&along.w.deriv())?)))
                .ok_or_else(unsupported)?;
            let wmax = along.w.bound();
            let floor = 1e-12 * (d.degree() + 1) as f64 * (along.k + 1.0) * along.scale * wmax;
            for s in d.cleaned(floor).roots() {
                push(Split {
                    t: piece.map.t(s),
                    turn: true,
                });
            }
        }
        if t0 > t1 {
            out.reverse();
        }
        Ok(out)
    }

    /// The parameter where the edge `curve` crosses the surface between
    /// `ta` and `tb`, given that the implicit function is monotone between
    /// them and the signed distances there have strictly opposite signs:
    /// bisection on the sign of the distance, to the last bit of the
    /// parameter.
    pub fn crossing(&self, curve: &Curve3, mut ta: f64, mut tb: f64) -> f64 {
        let mut da = self.distance(curve.point(ta));
        let mut db = self.distance(curve.point(tb));
        let below = da < 0.0;
        for _ in 0..BISECTIONS {
            let mid = 0.5 * (ta + tb);
            if mid == ta || mid == tb {
                break;
            }
            let d = self.distance(curve.point(mid));
            if d == 0.0 {
                return mid;
            }
            if (d < 0.0) == below {
                (ta, da) = (mid, d);
            } else {
                (tb, db) = (mid, d);
            }
        }
        if da.abs() <= db.abs() {
            ta
        } else {
            tb
        }
    }
}

/// A box holding the edge `curve` on `[t0, t1]` (either order): the box of
/// its rational Bézier pieces' control points, whose convex hull holds each
/// piece (the weights are positive). `None` for a NURBS curve of degree
/// above five.
pub(super) fn curve_box(curve: &Curve3, t0: f64, t1: f64) -> Option<(DVec3, DVec3)> {
    let (lo, hi) = if t0 <= t1 { (t0, t1) } else { (t1, t0) };
    let mut b = (DVec3::splat(f64::INFINITY), DVec3::splat(f64::NEG_INFINITY));
    for piece in pieces(curve, lo, hi)? {
        for q in &piece.pw {
            let x = q.truncate() / q.w;
            b = (b.0.min(x), b.1.max(x));
        }
    }
    Some(b)
}

/// The control points of the rational Bézier pieces of `curve` on `[t0,
/// t1]` (either order), dehomogenised: the convex hull of each piece's
/// points holds the piece. `None` for a NURBS curve of degree above five.
pub(super) fn hull_points(curve: &Curve3, t0: f64, t1: f64) -> Option<Vec<Vec<DVec3>>> {
    let (lo, hi) = if t0 <= t1 { (t0, t1) } else { (t1, t0) };
    Some(
        pieces(curve, lo, hi)?
            .iter()
            .map(|piece| piece.pw.iter().map(|q| q.truncate() / q.w).collect())
            .collect(),
    )
}

/// The rational Bézier pieces of `curve` on `[t0, t1]` (`t0 ≤ t1`), in
/// order: a line whole, a circle in equal pieces of at most a quarter
/// turn (Piegl & Tiller §7.5), a NURBS curve by its knot spans, clipped.
/// `None` for a NURBS curve of degree above five.
fn pieces(curve: &Curve3, t0: f64, t1: f64) -> Option<Vec<Piece>> {
    let one = |p: DVec3| p.extend(1.0);
    Some(match curve {
        Curve3::Line { .. } => vec![Piece {
            pw: vec![one(curve.point(t0)), one(curve.point(t1))],
            map: Map::Linear(t0, t1),
        }],
        Curve3::Circle {
            center,
            x,
            y,
            radius,
        } => {
            let span = t1 - t0;
            let k = ((span / std::f64::consts::FRAC_PI_2).ceil() as usize).max(1);
            let delta = span / k as f64;
            let w = (0.5 * delta).cos();
            (0..k)
                .map(|i| {
                    let a = t0 + delta * i as f64;
                    let b = if i + 1 == k { t1 } else { a + delta };
                    let mid = 0.5 * (a + b);
                    let (s, c) = mid.sin_cos();
                    let corner = *center + (*x * c + *y * s) * (*radius / w);
                    Piece {
                        pw: vec![
                            one(curve.point(a)),
                            (corner * w).extend(w),
                            one(curve.point(b)),
                        ],
                        map: Map::Arc {
                            mid,
                            half: (0.25 * (b - a)).tan(),
                        },
                    }
                })
                .collect()
        }
        Curve3::Nurbs(n) => {
            if n.degree() > 5 {
                return None;
            }
            n.bezier_segments()
                .into_iter()
                .filter(|seg| seg.b > t0 && seg.a < t1)
                .map(|seg| {
                    let (a, b) = (seg.a.max(t0), seg.b.min(t1));
                    let len = seg.b - seg.a;
                    let pw = clip(&seg.points, (a - seg.a) / len, (b - seg.a) / len);
                    Piece {
                        pw,
                        map: Map::Linear(a, b),
                    }
                })
                .collect()
        }
    })
}

/// The part of the Bézier curve with control points `pw` on `[sa, sb]`,
/// reparameterised to `[0, 1]`, by de Casteljau subdivision.
fn clip(pw: &[DVec4], sa: f64, sb: f64) -> Vec<DVec4> {
    // Keep [0, sb], then the part of that on [sa / sb, 1].
    let left = if sb < 1.0 {
        subdivide(pw, sb).0
    } else {
        pw.to_vec()
    };
    if sa > 0.0 {
        subdivide(&left, sa / sb).1
    } else {
        left
    }
}

/// The two parts of a Bézier curve at `s`.
fn subdivide(pw: &[DVec4], s: f64) -> (Vec<DVec4>, Vec<DVec4>) {
    let n = pw.len();
    let mut b = pw.to_vec();
    let mut left = Vec::with_capacity(n);
    let mut right = vec![DVec4::ZERO; n];
    left.push(b[0]);
    right[n - 1] = b[n - 1];
    for r in 1..n {
        for i in 0..n - r {
            b[i] = b[i] * (1.0 - s) + b[i + 1] * s;
        }
        left.push(b[0]);
        right[n - 1 - r] = b[n - 1 - r];
    }
    (left, right)
}

#[cfg(test)]
mod tests {
    use std::f64::consts::{FRAC_PI_2, PI, TAU};

    use glam::DVec2;
    use opendatum_nurbs::NurbsCurve;

    use super::super::testing::Rng;
    use super::*;
    use crate::math::Frame;
    use crate::profile::{Profile, ProfileCurve};
    use crate::sweep::{self, ExtrudeSpec, RevolveSpec};
    use crate::topo::Brep;
    use crate::Tolerance;

    /// A curve of a 2D profile curve in `frame`, as the kernel stores it.
    fn nurbs_in(frame: &Frame, c: &ProfileCurve) -> Curve3 {
        let ProfileCurve::Spline {
            degree,
            knots,
            points,
            weights,
        } = c
        else {
            panic!("a spline");
        };
        let p3 = points.iter().map(|p| frame.point2(*p)).collect();
        Curve3::Nurbs(NurbsCurve::new(*degree, knots.clone(), p3, weights.clone()).unwrap())
    }

    /// Points on every edge of every face, with the face's surface.
    fn face_points(b: &Brep) -> Vec<(&Surface, Vec<DVec3>)> {
        b.faces
            .iter()
            .map(|f| {
                let mut pts = Vec::new();
                for l in &f.loops {
                    for c in &l.coedges {
                        let e = b.edge(c.edge);
                        for i in 0..=8 {
                            pts.push(e.curve.point(e.t0 + (e.t1 - e.t0) * i as f64 / 8.0));
                        }
                    }
                }
                (&f.surface, pts)
            })
            .collect()
    }

    fn kind(q: &Quadric) -> &'static str {
        match q {
            Quadric::Plane { .. } => "plane",
            Quadric::Cylinder { .. } => "cylinder",
            Quadric::Cone { .. } => "cone",
            Quadric::Sphere { .. } => "sphere",
        }
    }

    #[test]
    fn the_surfaces_of_extrudes_and_revolves_are_recognised() {
        let tol = Tolerance::default();
        let plane = Frame::from_normal(
            DVec3::new(0.1, -0.2, 0.3),
            DVec3::new(1.0, 2.0, 2.5),
            DVec3::X,
        );
        // A rounded slot with a round hole, extruded.
        let r = 0.25;
        let slot = Profile {
            plane,
            curves: vec![
                ProfileCurve::Line {
                    start: DVec2::new(0.0, -r),
                    end: DVec2::new(1.0, -r),
                },
                ProfileCurve::Arc {
                    center: DVec2::new(1.0, 0.0),
                    radius: r,
                    start_angle: -FRAC_PI_2,
                    sweep_angle: PI,
                },
                ProfileCurve::Line {
                    start: DVec2::new(1.0, r),
                    end: DVec2::new(0.0, r),
                },
                ProfileCurve::Arc {
                    center: DVec2::ZERO,
                    radius: r,
                    start_angle: FRAC_PI_2,
                    sweep_angle: PI,
                },
                ProfileCurve::Circle {
                    center: DVec2::new(0.5, 0.0),
                    radius: 0.1,
                },
            ],
        };
        let spec = ExtrudeSpec {
            distance: 0.4,
            symmetric: false,
        };
        let (slot, _) = sweep::extrude(&slot, &spec, &tol).unwrap();
        // A revolve of a profile with a slanted line (a cone), lines along
        // and across the axis (a cylinder, discs) and a half circle centred
        // on the axis (a sphere).
        let turned = Profile {
            plane,
            curves: vec![
                ProfileCurve::Line {
                    start: DVec2::new(0.0, 0.0),
                    end: DVec2::new(0.5, 0.0),
                },
                ProfileCurve::Line {
                    start: DVec2::new(0.5, 0.0),
                    end: DVec2::new(0.8, 0.6),
                },
                ProfileCurve::Line {
                    start: DVec2::new(0.8, 0.6),
                    end: DVec2::new(0.8, 1.2),
                },
                ProfileCurve::Arc {
                    center: DVec2::new(0.0, 1.2),
                    radius: 0.8,
                    start_angle: 0.0,
                    sweep_angle: FRAC_PI_2,
                },
                ProfileCurve::Line {
                    start: DVec2::new(0.0, 2.0),
                    end: DVec2::new(0.0, 0.0),
                },
            ],
        };
        let axis = RevolveSpec {
            axis_origin: DVec2::ZERO,
            axis_dir: DVec2::Y,
            angle: TAU,
        };
        let (turned, _) = sweep::revolve(&turned, &axis, &tol).unwrap();
        // A torus is not a quadric.
        let ring = Profile::circle(plane, DVec2::new(1.0, 0.0), 0.3);
        let (ring, _) = sweep::revolve(&ring, &axis, &tol).unwrap();
        let mut seen = std::collections::BTreeMap::new();
        for b in [&slot, &turned] {
            for (surface, pts) in face_points(b) {
                let q = Quadric::of(surface).expect("a quadric");
                *seen.entry(kind(&q)).or_insert(0) += 1;
                let mag = 1.0 + q.magnitude();
                for p in pts {
                    let d = q.distance(p);
                    assert!(d.abs() <= 1e-14 * mag, "{} off by {d}", kind(&q));
                }
            }
        }
        // Slot: 2 planar sides, a half cylinder at each end, 2 for the
        // hole, 2 caps; revolve: 2 cones, 2 cylinders, 2 spheres, 2 discs.
        let want: Vec<(&str, i32)> =
            vec![("cone", 2), ("cylinder", 6), ("plane", 6), ("sphere", 2)];
        assert_eq!(seen.into_iter().collect::<Vec<_>>(), want);
        for f in &ring.faces {
            assert_eq!(Quadric::of(&f.surface), None);
        }
        // The two halves of a cylinder are one surface; different
        // cylinders are not.
        for (b, want) in [(&slot, 1), (&turned, 3)] {
            let cylinders: Vec<Quadric> = b
                .faces
                .iter()
                .filter_map(|f| Quadric::of(&f.surface))
                .filter(|q| !matches!(q, Quadric::Plane { .. }))
                .collect();
            let pairs = (0..cylinders.len())
                .flat_map(|i| (i + 1..cylinders.len()).map(move |j| (i, j)))
                .filter(|&(i, j)| cylinders[i].same(&cylinders[j]))
                .count();
            // The slot's hole; the revolve's cone, cylinder and sphere.
            assert_eq!(pairs, want);
        }
    }

    #[test]
    fn distances_are_signed_and_exact() {
        let mut rng = Rng(0xD157_0001);
        for _ in 0..2000 {
            let axis = rng.unit();
            let radial = axis.any_orthonormal_vector();
            let radial = (radial * rng.range(-1.0, 1.0)
                + axis.cross(radial) * rng.range(-1.0, 1.0))
            .normalize();
            let point = rng.point(2.0);
            let h = rng.range(-1.0, 1.0);
            let delta = rng.range(-0.05, 0.05);
            let radius = rng.range(0.1, 1.0);
            let near = |a: f64, b: f64| (a - b).abs() <= 1e-14 * 8.0;
            let plane = Quadric::Plane {
                point,
                normal: axis,
            };
            assert!(near(plane.distance(point + radial + axis * delta), delta));
            let cylinder = Quadric::Cylinder {
                point,
                axis,
                radius,
            };
            let p = point + axis * h + radial * (radius + delta);
            assert!(near(cylinder.distance(p), delta));
            let sphere = Quadric::Sphere {
                center: point,
                radius,
            };
            let dir = (axis * h + radial).normalize();
            assert!(near(sphere.distance(point + dir * (radius + delta)), delta));
            // A cone, measured along its normal from a point of its nappe.
            let slope = 10f64.powf(rng.range(-2.0, 1.0));
            let cone = Quadric::Cone {
                point,
                axis,
                radius,
                slope,
            };
            let at = (radius + slope * h).max(0.05);
            let h = (at - radius) / slope;
            let cos = 1.0 / (1.0 + slope * slope).sqrt();
            let sin = slope * cos;
            let p = point + axis * (h - delta * sin) + radial * (at + delta * cos);
            let d = cone.distance(p);
            assert!(
                (d - delta).abs() <= 1e-13 * (1.0 + h.abs()),
                "{d} vs {delta}"
            );
            // Behind the apex, the distance to the apex.
            let apex = point - axis * (radius / slope);
            let behind = apex - axis * 0.3 + radial * (0.15 / slope);
            assert!(
                (cone.distance(behind) - (behind - apex).length()).abs()
                    < 1e-12 * (1.0 + 1.0 / slope)
            );
            // Moving the reference point keeps the surface.
            let moved = cone.anchored(point + axis * rng.range(-1.0, 1.0) + radial);
            assert!(cone.same(&moved) && moved.same(&cone));
            assert!(near(moved.distance(p), d));
            let shifted = Quadric::Cone {
                point,
                axis,
                radius: radius + 1e-6,
                slope,
            };
            assert!(!cone.same(&shifted));
            let flipped = Quadric::Plane {
                point: point + radial,
                normal: -axis,
            };
            assert!(plane.same(&flipped));
            assert!(!plane.same(&cylinder) && !cylinder.same(&sphere));
        }
    }

    /// A random line, circular arc, elliptical arc, conic or cubic spline
    /// in the box `[−1.5, 1.5]³`, with a parameter range.
    fn random_edge(rng: &mut Rng) -> (Curve3, f64, f64) {
        let frame = Frame::from_normal(rng.point(1.0), rng.unit(), rng.unit());
        let sub = |rng: &mut Rng, (a, b): (f64, f64)| {
            let (u, v) = (rng.next(), rng.next());
            let (u, v) = (u.min(v) * 0.5, 0.5 + u.max(v) * 0.5);
            (a + (b - a) * u, a + (b - a) * v)
        };
        match rng.below(5) {
            0 => {
                let a = rng.point(1.5);
                (
                    Curve3::Line {
                        origin: a,
                        dir: rng.point(1.5) - a,
                    },
                    0.0,
                    1.0,
                )
            }
            1 => {
                let t0 = rng.range(-FRAC_PI_2, PI);
                let c = Curve3::Circle {
                    center: frame.origin,
                    x: frame.x,
                    y: frame.y,
                    radius: rng.range(0.2, 1.5),
                };
                (c, t0, t0 + rng.range(0.1, 1.9 * PI))
            }
            2 => {
                let major = DVec2::new(rng.range(0.2, 1.5), 0.0);
                let minor = DVec2::new(rng.range(-0.3, 0.3), rng.range(0.2, 1.5));
                let e = ProfileCurve::elliptical_arc(
                    DVec2::ZERO,
                    major,
                    minor,
                    rng.range(-PI, PI),
                    rng.range(0.2, TAU),
                )
                .unwrap();
                let c = nurbs_in(&frame, &e);
                let Curve3::Nurbs(n) = &c else { unreachable!() };
                let (a, b) = sub(rng, n.domain());
                (c, a, b)
            }
            3 => {
                let e = ProfileCurve::conic(
                    DVec2::new(rng.range(-1.0, 1.0), rng.range(-1.0, 1.0)),
                    DVec2::new(rng.range(-1.0, 1.0), rng.range(-1.0, 1.0)),
                    DVec2::new(rng.range(-1.0, 1.0), rng.range(-1.0, 1.0)),
                    10f64.powf(rng.range(-0.7, 0.7)),
                )
                .unwrap();
                (nurbs_in(&frame, &e), 0.0, 1.0)
            }
            _ => {
                let pts: Vec<DVec3> = (0..6).map(|_| rng.point(1.2)).collect();
                let n = NurbsCurve::interpolate(&pts, 3).unwrap();
                let (a, b) = sub(rng, n.domain());
                (Curve3::Nurbs(n), a, b)
            }
        }
    }

    fn random_quadric(rng: &mut Rng) -> Quadric {
        let (point, axis) = (rng.point(0.5), rng.unit());
        match rng.below(4) {
            0 => Quadric::Plane {
                point,
                normal: axis,
            },
            1 => Quadric::Cylinder {
                point,
                axis,
                radius: rng.range(0.2, 1.0),
            },
            2 => Quadric::Cone {
                point,
                axis,
                radius: rng.range(0.0, 1.0),
                slope: 10f64.powf(rng.range(-1.5, 1.0)),
            },
            _ => Quadric::Sphere {
                center: point,
                radius: rng.range(0.2, 1.2),
            },
        }
    }

    /// Checks the splits and crossings of an edge against dense sampling of
    /// the distance; returns the number of crossings.
    fn check(q: &Quadric, curve: &Curve3, t0: f64, t1: f64, case: usize) -> usize {
        let splits = q.splits(curve, t0, t1).unwrap();
        let mut ends = vec![t0];
        ends.extend(splits.iter().map(|s| s.t));
        ends.push(t1);
        assert!(
            ends.windows(2).all(|w| w[0] < w[1]),
            "case {case}: {ends:?}"
        );
        let mag = 1.0 + q.magnitude() + curve.magnitude();
        let noise = 1e-11 * mag;
        let d = |t: f64| q.distance(curve.point(t));
        let mut found = Vec::new();
        for w in ends.windows(2) {
            let (da, db) = (d(w[0]), d(w[1]));
            if (da < -noise && db > noise) || (da > noise && db < -noise) {
                let t = q.crossing(curve, w[0], w[1]);
                assert!(t >= w[0] && t <= w[1]);
                assert!(d(t).abs() <= 1e-12 * mag, "case {case}: {} off", d(t));
                found.push(t);
            }
        }
        // Sampled: no crossing is missed, and within a stretch the
        // distance changes sign at most once.
        let m = 3000;
        let ts: Vec<f64> = (0..=m)
            .map(|i| t0 + (t1 - t0) * i as f64 / m as f64)
            .collect();
        let ds: Vec<f64> = ts.iter().map(|&t| d(t)).collect();
        let strict = |x: f64| {
            if x > noise {
                1
            } else if x < -noise {
                -1
            } else {
                0
            }
        };
        for i in 0..m {
            let (a, b) = (strict(ds[i]), strict(ds[i + 1]));
            if a * b < 0 {
                assert!(
                    found.iter().any(|&t| t >= ts[i] && t <= ts[i + 1]),
                    "case {case}: a crossing in [{}, {}] missed; found {found:?}",
                    ts[i],
                    ts[i + 1]
                );
            }
        }
        for w in ends.windows(2) {
            let signs: Vec<i32> = ts
                .iter()
                .zip(&ds)
                .filter(|(t, _)| **t > w[0] && **t < w[1])
                .map(|(_, x)| strict(*x))
                .filter(|s| *s != 0)
                .collect();
            let changes = signs.windows(2).filter(|p| p[0] != p[1]).count();
            assert!(
                changes <= 1,
                "case {case}: {changes} sign changes in a stretch"
            );
        }
        // The implicit function, computed here independently, is monotone
        // between consecutive splits and has an extremum at each turn.
        let g = |t: f64| implicit(q, curve.point(t));
        let gnoise = 1e-10 * mag * mag;
        for w in ends.windows(2) {
            let gs: Vec<f64> = ts
                .iter()
                .filter(|t| **t > w[0] && **t < w[1])
                .map(|&t| g(t))
                .collect();
            let up = gs.windows(2).any(|p| p[1] > p[0] + gnoise);
            let down = gs.windows(2).any(|p| p[1] < p[0] - gnoise);
            assert!(
                !(up && down),
                "case {case}: not monotone on [{}, {}]",
                w[0],
                w[1]
            );
        }
        let h = (t1 - t0) * 1e-4;
        for s in splits.iter().filter(|s| s.turn) {
            let (a, m, b) = (g(s.t - h), g(s.t), g(s.t + h));
            let min = a >= m - gnoise && b >= m - gnoise;
            let max = a <= m + gnoise && b <= m + gnoise;
            assert!(
                min || max,
                "case {case}: no extremum at {}: {a} {m} {b}",
                s.t
            );
        }
        found.len()
    }

    /// The implicit function, from its formula.
    fn implicit(q: &Quadric, p: DVec3) -> f64 {
        match *q {
            Quadric::Plane { point, normal } => (p - point).dot(normal),
            Quadric::Cylinder {
                point,
                axis,
                radius,
            } => {
                let w = p - point;
                (w - axis * w.dot(axis)).length_squared() - radius * radius
            }
            Quadric::Cone {
                point,
                axis,
                radius,
                slope,
            } => {
                let w = p - point;
                let h = w.dot(axis);
                (w - axis * h).length_squared() - (radius + slope * h).powi(2)
            }
            Quadric::Sphere { center, radius } => (p - center).length_squared() - radius * radius,
        }
    }

    #[test]
    fn splits_and_crossings_match_sampling() {
        let mut rng = Rng(0xC055_0001);
        let mut crossings = 0;
        let mut turns = 0;
        for case in 0..3000 {
            let q = random_quadric(&mut rng);
            let (curve, t0, t1) = random_edge(&mut rng);
            crossings += check(&q, &curve, t0, t1, case);
            turns += q
                .splits(&curve, t0, t1)
                .unwrap()
                .iter()
                .filter(|s| s.turn)
                .count();
            // The same edge the other way round.
            let back = q.splits(&curve, t1, t0).unwrap();
            let fwd = q.splits(&curve, t0, t1).unwrap();
            assert!(back.iter().rev().eq(fwd.iter()));
        }
        // The test exercises what it should.
        assert!(
            crossings > 2000 && turns > 1000,
            "{crossings} crossings, {turns} turns"
        );
    }

    #[test]
    fn edges_lying_on_a_surface_have_no_turns() {
        let frame = Frame::from_normal(
            DVec3::new(0.3, -0.2, 0.1),
            DVec3::new(0.2, 1.0, -0.4),
            DVec3::X,
        );
        let (o, z) = (frame.origin, frame.z);
        let circle = |c: DVec3, r: f64| Curve3::Circle {
            center: c,
            x: frame.x,
            y: frame.y,
            radius: r,
        };
        let cylinder = Quadric::Cylinder {
            point: o,
            axis: z,
            radius: 0.4,
        };
        let cone = Quadric::Cone {
            point: o,
            axis: z,
            radius: 0.4,
            slope: 0.3,
        };
        let sphere = Quadric::Sphere {
            center: o,
            radius: 0.5,
        };
        let plane = Quadric::Plane {
            point: o,
            normal: z,
        };
        // An ellipse where a plane at 40° cuts the cylinder.
        let tilt = 40f64.to_radians();
        let major = (frame.x * tilt.cos() + z * tilt.sin()) * (0.4 / tilt.cos());
        let e = ProfileCurve::elliptical_arc(DVec2::ZERO, DVec2::X, DVec2::Y, 0.3, 5.0).unwrap();
        let Curve3::Nurbs(unit) = nurbs_in(&Frame::WORLD, &e) else {
            unreachable!()
        };
        let m = glam::DMat4::from_cols(
            major.extend(0.0),
            (frame.y * 0.4).extend(0.0),
            z.extend(0.0),
            o.extend(1.0),
        );
        let ellipse = Curve3::Nurbs(unit.transformed(&m));
        let cases: Vec<(Quadric, Curve3, f64, f64)> = vec![
            (cylinder, circle(o + z * 0.7, 0.4), 0.2, 3.0),
            (
                cylinder,
                Curve3::Line {
                    origin: o + frame.y * 0.4,
                    dir: z,
                },
                -1.0,
                2.0,
            ),
            (cylinder, ellipse.clone(), 0.0, 1.0),
            (cone, circle(o + z * 0.5, 0.55), -1.0, 4.0),
            (
                cone,
                Curve3::Line {
                    origin: o + frame.x * 0.4,
                    dir: z + frame.x * 0.3,
                },
                -0.5,
                2.0,
            ),
            (sphere, circle(o, 0.5), 0.0, 6.0),
            (sphere, circle(o + z * 0.3, 0.4), 1.0, 2.5),
            (plane, circle(o, 0.7), 0.0, 6.0),
            (
                plane,
                Curve3::Nurbs(unit.transformed(&frame.matrix())),
                0.0,
                1.0,
            ),
        ];
        for (i, (q, c, t0, t1)) in cases.iter().enumerate() {
            let splits = q.splits(c, *t0, *t1).unwrap();
            assert!(splits.iter().all(|s| !s.turn), "case {i}: {splits:?}");
            for k in 0..=50 {
                let t = t0 + (t1 - t0) * k as f64 / 50.0;
                assert!(q.distance(c.point(t)).abs() < 1e-14 * 8.0, "case {i}");
            }
        }
    }

    #[test]
    fn turning_points_at_joints_and_tangent_contacts() {
        // A half circle whose highest point is the joint of its two
        // quarter-turn pieces, cut below that point: two crossings.
        let half = Curve3::Circle {
            center: DVec3::ZERO,
            x: DVec3::X,
            y: DVec3::Y,
            radius: 1.0,
        };
        let plane = Quadric::Plane {
            point: DVec3::new(0.0, 0.5, 0.0),
            normal: DVec3::Y,
        };
        assert_eq!(check(&plane, &half, 0.0, PI, 0), 2);
        let splits = plane.splits(&half, 0.0, PI).unwrap();
        let ends: Vec<f64> = splits.iter().map(|s| s.t).collect();
        let mut found: Vec<f64> = Vec::new();
        let mut all = vec![0.0];
        all.extend(ends);
        all.push(PI);
        for w in all.windows(2) {
            let (a, b) = (
                plane.distance(half.point(w[0])),
                plane.distance(half.point(w[1])),
            );
            if a * b < 0.0 {
                found.push(plane.crossing(&half, w[0], w[1]));
            }
        }
        assert!((found[0] - PI / 6.0).abs() < 1e-15 && (found[1] - 5.0 * PI / 6.0).abs() < 1e-15);
        // A line tangent to a cylinder, and a circle touching a plane: one
        // turning point each, on the surface.
        let cylinder = Quadric::Cylinder {
            point: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 0.5,
        };
        let line = Curve3::Line {
            origin: DVec3::new(-1.0, 0.5, 0.3),
            dir: DVec3::new(2.0, 0.0, 0.1),
        };
        let s = cylinder.splits(&line, 0.0, 1.0).unwrap();
        assert_eq!(s.len(), 1);
        assert!(s[0].turn && cylinder.distance(line.point(s[0].t)).abs() < 1e-15);
        let floor = Quadric::Plane {
            point: DVec3::new(0.0, -1.0, 0.0),
            normal: DVec3::Y,
        };
        let s = floor.splits(&half, 3.0, 6.0).unwrap();
        let turns: Vec<&Split> = s.iter().filter(|s| s.turn).collect();
        assert_eq!(turns.len(), 1);
        assert!((turns[0].t - 1.5 * PI).abs() < 1e-7);
        assert!(floor.distance(half.point(turns[0].t)).abs() < 1e-15);
    }
}
