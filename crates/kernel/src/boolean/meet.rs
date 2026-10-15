//! Where two surfaces of step 2 meet, in closed form (design note 0012
//! §10.3 and §14): lines, circles, ellipses and conics.
//!
//! Two surfaces of revolution about one axis (a plane across the axis, a
//! cylinder, a cone, a sphere with its centre on the axis) meet in circles
//! about that axis, where their meridians (their sections by a half-plane
//! bounded by the axis: lines and circles) meet. Otherwise two planes meet
//! in a line, a plane and a cylinder in two lines along the axis or an
//! ellipse, two cylinders with parallel axes in two lines, and a plane
//! across a cone in a conic (§14): an ellipse, a parabola or one branch of
//! a hyperbola, or two rulings where the plane passes through the apex.
//! Other pairs meet in curves of degree four; they are refused as
//! unsupported.
//!
//! Every decision is made within the tolerance `eps`: surfaces that stay
//! within it of each other over the region where the faces are are one
//! surface; surfaces whose least gap is within it touch, which is refused
//! as a tangent contact until tangency is handled (design note 0011 §4.4).
//! A configuration within a sixteenth of the tolerance of a special one
//! over the region (axes parallel, a plane square to an axis, a centre on
//! an axis) is treated as that one, so that no curve is computed from a
//! near-zero angle: the curves returned lie within `eps / 16` of both
//! surfaces in the region, and exactly on them otherwise.

use std::f64::consts::{FRAC_PI_2, PI};

use glam::{DVec2, DVec3};
use opendatum_nurbs::NurbsCurve;

use super::quadric::Quadric;
use crate::geom::Curve3;
use crate::KernelError;

/// The part of space where two faces are: their common bounding box
/// (`lo` to `hi`), the ball holding it, and for a pair with a cone the
/// least height above its apex, along its axis, of the faces on it
/// (`lift`: beside an oblique axis the box can hold the apex although no
/// face comes near it).
#[derive(Clone, Copy, Debug, PartialEq)]
pub(super) struct Region {
    pub center: DVec3,
    pub radius: f64,
    pub lo: DVec3,
    pub hi: DVec3,
    pub lift: f64,
}

impl Region {
    /// The region of the box from `lo` to `hi`, faces on a cone reaching
    /// its apex.
    pub fn of_box(lo: DVec3, hi: DVec3) -> Region {
        Region {
            center: (lo + hi) * 0.5,
            radius: (hi - lo).length() * 0.5,
            lo,
            hi,
            lift: f64::NEG_INFINITY,
        }
    }

    /// True if `p` is within `margin` of the box.
    fn holds(&self, p: DVec3, margin: f64) -> bool {
        p.cmpge(self.lo - margin).all() && p.cmple(self.hi + margin).all()
    }

    /// True if faces on the pair's cone can come within `margin` of its
    /// apex `apex` (and within `height` of it along its axis), where faces
    /// on the other surface can be.
    fn near_apex(&self, apex: DVec3, margin: f64, height: f64) -> bool {
        self.lift <= height && self.holds(apex, margin)
    }

    /// The parameters from `start` (at least zero) where the ray `p + s ·
    /// dir` is in the box, if it reaches it.
    fn ray(&self, p: DVec3, dir: DVec3, start: f64) -> Option<(f64, f64)> {
        let (mut s0, mut s1) = (start, f64::INFINITY);
        for i in 0..3 {
            if dir[i] == 0.0 {
                if p[i] < self.lo[i] || p[i] > self.hi[i] {
                    return None;
                }
                continue;
            }
            let a = (self.lo[i] - p[i]) / dir[i];
            let b = (self.hi[i] - p[i]) / dir[i];
            s0 = s0.max(a.min(b));
            s1 = s1.min(a.max(b));
        }
        (s0 <= s1).then_some((s0, s1))
    }
}

/// A curve where two surfaces meet.
#[derive(Clone, Copy, Debug, PartialEq)]
pub(super) enum Track {
    /// `origin + t · dir`, `dir` a unit vector.
    Line { origin: DVec3, dir: DVec3 },
    /// A ruling of a cone from its apex, where a plane through the apex
    /// crosses it (§14): `origin + t · dir` for `t ≥ start`, `dir` a unit
    /// vector away from the apex. The apex is at `start` when `from_apex`;
    /// else the ray starts where it comes within [`CLIP`] tolerances of the
    /// region, the apex lying farther away (at infinity where it stays
    /// farther): before there no point is in faces on both surfaces.
    Ray {
        origin: DVec3,
        dir: DVec3,
        start: f64,
        from_apex: bool,
    },
    /// `center + radius (cos t · x + sin t · y)`, `x` and `y` orthonormal.
    Circle {
        center: DVec3,
        x: DVec3,
        y: DVec3,
        radius: f64,
    },
    /// `center + a cos t · x + b sin t · y`, `x` and `y` orthonormal.
    Ellipse {
        center: DVec3,
        x: DVec3,
        y: DVec3,
        a: f64,
        b: f64,
    },
    /// One branch of a hyperbola or a parabola, or the half of an ellipse
    /// at one vertex, in the vertex form (design note 0012 §14):
    /// `vertex + x · p t² / (1 + w) + y · p t` with `w = √(1 − c t²)`, the
    /// points `(ξ, η)` from the vertex along `x` and `y` with `η² = 2pξ −
    /// cξ²`. `p > 0` is the semi-latus rectum and `c = 1 − e²` (`e` the
    /// eccentricity): negative on a hyperbola, zero on a parabola, positive
    /// on an ellipse, whose half here is `|t| < 1/√c`. `x` (the way the
    /// curve opens) and `y` are orthonormal.
    Conic {
        vertex: DVec3,
        x: DVec3,
        y: DVec3,
        p: f64,
        c: f64,
    },
}

/// How two surfaces meet.
#[derive(Clone, Debug, PartialEq)]
pub(super) enum Meet {
    /// They do not meet.
    Apart,
    /// They are one surface within the tolerance over the region.
    Same,
    /// They cross along these curves.
    Curves(Vec<Track>),
    /// They touch, along a curve or at a point.
    Touch,
    /// Their curves are not handled in closed form yet, or not within the
    /// tolerance: why, ending with the design note's section.
    Unsupported(&'static str),
}

impl Track {
    /// The point at `t`.
    pub fn point(&self, t: f64) -> DVec3 {
        match *self {
            Track::Line { origin, dir } | Track::Ray { origin, dir, .. } => origin + dir * t,
            Track::Circle {
                center,
                x,
                y,
                radius,
            } => {
                let (s, c) = t.sin_cos();
                center + (x * c + y * s) * radius
            }
            Track::Ellipse { center, x, y, a, b } => {
                let (s, c) = t.sin_cos();
                center + x * (a * c) + y * (b * s)
            }
            Track::Conic {
                vertex, x, y, p, c, ..
            } => vertex + x * (p * t * t / (1.0 + conic_w(c, t))) + y * (p * t),
        }
    }

    /// The first derivative at `t`.
    pub fn deriv(&self, t: f64) -> DVec3 {
        match *self {
            Track::Line { dir, .. } | Track::Ray { dir, .. } => dir,
            Track::Circle { x, y, radius, .. } => {
                let (s, c) = t.sin_cos();
                (y * c - x * s) * radius
            }
            Track::Ellipse { x, y, a, b, .. } => {
                let (s, c) = t.sin_cos();
                y * (b * c) - x * (a * s)
            }
            Track::Conic { x, y, p, c, .. } => x * (p * t / conic_w(c, t)) + y * p,
        }
    }

    /// The second derivative at `t`.
    pub fn second(&self, t: f64) -> DVec3 {
        match *self {
            Track::Line { .. } | Track::Ray { .. } => DVec3::ZERO,
            Track::Circle { center, .. } | Track::Ellipse { center, .. } => center - self.point(t),
            Track::Conic { x, p, c, .. } => x * (p / conic_w(c, t).powi(3)),
        }
    }

    /// The parameter of a point on the curve (within rounding of it); for
    /// a closed curve, in `[−π/2, 3π/2)` like the kernel's circles.
    pub fn param(&self, p: DVec3) -> f64 {
        match *self {
            Track::Line { origin, dir } | Track::Ray { origin, dir, .. } => (p - origin).dot(dir),
            Track::Circle { center, x, y, .. } => {
                let d = p - center;
                wrap(d.dot(y).atan2(d.dot(x)))
            }
            Track::Ellipse { center, x, y, a, b } => {
                let d = p - center;
                wrap((d.dot(y) / b).atan2(d.dot(x) / a))
            }
            Track::Conic {
                vertex, y, p: semi, ..
            } => (p - vertex).dot(y) / semi,
        }
    }

    /// True for a closed curve (a circle or an ellipse), of period 2π.
    pub fn closed(&self) -> bool {
        matches!(self, Track::Circle { .. } | Track::Ellipse { .. })
    }

    /// True for a line or a ray.
    pub fn straight(&self) -> bool {
        matches!(self, Track::Line { .. } | Track::Ray { .. })
    }

    /// The parameters of an open curve's ends, infinite where it runs
    /// off: a line, a parabola and a branch of a hyperbola run off both
    /// ways, a ray from the apex, and the half of an ellipse ends on its
    /// minor axis (a closed curve gives both infinite).
    pub fn ends(&self) -> (f64, f64) {
        match *self {
            Track::Ray { start, .. } => (start, f64::INFINITY),
            Track::Conic { c, .. } if c > 0.0 => {
                let end = 1.0 / c.sqrt();
                (-end, end)
            }
            _ => (f64::NEG_INFINITY, f64::INFINITY),
        }
    }

    /// The parameter of the point of the curve nearest `p` and its
    /// distance from `p`: in closed form for a line, a ray and a circle; for an
    /// ellipse or a conic by Newton's method on `(C(t) − p) · C'(t) = 0`
    /// from the parameter of `p` (its eccentric angle, or its position
    /// across the conic's axis), which for a point near the curve is near
    /// the foot (the questions asked of it are about points within a few
    /// tolerances of the curve). On half an ellipse also by its eccentric
    /// angle, the nearer foot taken: near the ends of the half, on its
    /// minor axis, a conic's parameter runs ever faster (`ξ' = p t / w`,
    /// `w → 0`) and Newton's method in it can stop short of the foot
    /// (found by review: 360 times the distance, 10⁻⁸ from an end).
    pub fn nearest(&self, p: DVec3) -> (f64, f64) {
        let t = match *self {
            Track::Line { origin, dir } => (p - origin).dot(dir),
            Track::Ray {
                origin, dir, start, ..
            } => (p - origin).dot(dir).max(start),
            Track::Circle { .. } => self.param(p),
            Track::Ellipse { center, x, y, a, b } => {
                let d = p - center;
                wrap(self.foot(p, (d.dot(y) / b).atan2(d.dot(x) / a), 0.25, f64::INFINITY))
            }
            Track::Conic {
                vertex,
                x,
                y,
                p: semi,
                c,
            } => {
                // On half an ellipse, within its half.
                let limit = if c > 0.0 {
                    (1.0 - 1e-12) / c.sqrt()
                } else {
                    f64::INFINITY
                };
                let t = self.param(p).clamp(-limit, limit);
                let t = self.foot(p, t, 0.25 * (1.0 + t.abs()), limit);
                if c <= 0.0 {
                    t
                } else {
                    // The whole ellipse, `ξ = a (1 − cos θ)`, `η = b sin θ`
                    // with `a = p / c` and `b = p / √c`, so `t = sin θ / √c`;
                    // a foot on the other half (`cos θ < 0`) puts the
                    // nearest point of this half at its end on that side.
                    let (a, b) = (semi / c, semi / c.sqrt());
                    let whole = Track::Ellipse {
                        center: vertex + x * a,
                        x: -x,
                        y,
                        a,
                        b,
                    };
                    let (theta, _) = whole.nearest(p);
                    let (sin, cos) = theta.sin_cos();
                    let end = 1.0 / c.sqrt();
                    let s = if cos >= 0.0 {
                        (sin * end).clamp(-end, end)
                    } else {
                        end.copysign(sin)
                    };
                    if (self.point(s) - p).length() < (self.point(t) - p).length() {
                        s
                    } else {
                        t
                    }
                }
            }
        };
        (t, (self.point(t) - p).length())
    }

    /// Newton's method for the foot of `p` from `t`, each step at most
    /// `most`, the parameter kept within `±limit`.
    fn foot(&self, p: DVec3, mut t: f64, most: f64, limit: f64) -> f64 {
        for _ in 0..NEWTON {
            let (c, c1) = (self.point(t) - p, self.deriv(t));
            let f = c.dot(c1);
            let df = c1.dot(c1) + c.dot(self.second(t));
            if df <= 0.0 || !df.is_finite() {
                break;
            }
            let step = (f / df).clamp(-most, most);
            t = (t - step).clamp(-limit, limit);
            if step.abs() <= 1e-15 * (1.0 + t.abs()) {
                break;
            }
        }
        t
    }

    /// A box holding the curve between parameters `ta` and `tb` (either
    /// order), exact: the ends and every turning point of a coordinate
    /// between them (where a coordinate of a circle or an ellipse is a
    /// sinusoid, and one of a conic is convex or concave in `t`).
    pub fn bounds(&self, ta: f64, tb: f64) -> (DVec3, DVec3) {
        let (lo, hi) = (ta.min(tb), ta.max(tb));
        let mut ts = vec![lo, hi];
        let mut turning = |t: f64| {
            if t > lo && t < hi {
                ts.push(t);
            }
        };
        match *self {
            Track::Line { .. } | Track::Ray { .. } => {}
            Track::Circle {
                x, y, radius: a, ..
            }
            | Track::Ellipse { x, y, a, .. } => {
                let b = match *self {
                    Track::Ellipse { b, .. } => b,
                    _ => a,
                };
                for i in 0..3 {
                    // a cos t xᵢ + b sin t yᵢ turns where tan t = b yᵢ / (a xᵢ).
                    let t0 = (b * y[i]).atan2(a * x[i]);
                    let first = ((lo - t0) / PI).ceil() as i64;
                    let last = ((hi - t0) / PI).floor() as i64;
                    for j in first..=last.min(first + 4) {
                        turning(t0 + PI * j as f64);
                    }
                }
            }
            Track::Conic { x, y, c, .. } => {
                for i in 0..3 {
                    // p (t/w xᵢ + yᵢ) is zero where t/w = m = −yᵢ/xᵢ, at
                    // t = m / √(1 + c m²).
                    if x[i] == 0.0 {
                        continue;
                    }
                    let m = -y[i] / x[i];
                    let q = 1.0 + c * m * m;
                    if q > 0.0 {
                        turning(m / q.sqrt());
                    }
                }
            }
        }
        ts.iter().fold(
            (DVec3::splat(f64::INFINITY), DVec3::splat(f64::NEG_INFINITY)),
            |(l, h), &t| {
                let q = self.point(t);
                (l.min(q), h.max(q))
            },
        )
    }

    /// The edge geometry of the curve from `ta` to `tb` (`ta < tb`, less
    /// than a turn apart on a closed curve), with its parameter range:
    /// a line (a ray's piece too) or circle as itself, an elliptical arc as an exact rational
    /// quadratic B-spline in pieces of at most a quarter turn (Piegl &
    /// Tiller A7.1 under the affine map from the unit circle, as the
    /// sketch's ellipses).
    pub fn edge(&self, ta: f64, tb: f64) -> Result<(Curve3, f64, f64), KernelError> {
        Ok(match *self {
            Track::Line { origin, dir } | Track::Ray { origin, dir, .. } => {
                (Curve3::Line { origin, dir }, ta, tb)
            }
            Track::Circle {
                center,
                x,
                y,
                radius,
            } => (
                Curve3::Circle {
                    center,
                    x,
                    y,
                    radius,
                },
                ta,
                tb,
            ),
            Track::Ellipse { center, x, y, a, b } => {
                let sweep = tb - ta;
                let k = ((sweep / FRAC_PI_2).ceil() as usize).max(1);
                let delta = sweep / k as f64;
                let w = (0.5 * delta).cos();
                let at = |t: f64, scale: f64| {
                    let (s, c) = t.sin_cos();
                    center + (x * (a * c) + y * (b * s)) * scale
                };
                let mut points = Vec::with_capacity(2 * k + 1);
                let mut weights = Vec::with_capacity(2 * k + 1);
                for i in 0..k {
                    let t0 = ta + delta * i as f64;
                    points.push(at(t0, 1.0));
                    weights.push(1.0);
                    points.push(at(t0 + 0.5 * delta, 1.0 / w));
                    weights.push(w);
                }
                points.push(at(tb, 1.0));
                weights.push(1.0);
                let mut knots = vec![0.0; 3];
                for i in 1..k {
                    let u = i as f64 / k as f64;
                    knots.extend([u, u]);
                }
                knots.extend([1.0; 3]);
                let n = NurbsCurve::new(2, knots, points, Some(weights))
                    .map_err(|e| KernelError::InvalidResult(format!("an elliptical edge: {e}")))?;
                (Curve3::Nurbs(n), 0.0, 1.0)
            }
            Track::Conic {
                vertex, x, y, p, c, ..
            } => {
                // In s = t / (1 + w) the conic is vertex + 2p (s² x + s y) /
                // (1 + c s²), a rational quadratic (design note 0012 §14),
                // cut into pieces turning at most a quarter turn each at
                // equal turning of the tangent. Its weights are taken in t,
                // where they do not cancel far along a branch: 1 + c s² =
                // 2 / (1 + w), and between pieces' ends t₀, t₁ the middle
                // weight is N / ((1 + w₀)(1 + w₁)) with N = 1 + w₀ + w₁ +
                // w₀w₁ + c t₀t₁, the last two summed as (1 − c t₀² − c t₁²) /
                // (w₀w₁ − c t₀t₁) where they would cancel.
                let (fa, fb) = (conic_turn(c, ta), conic_turn(c, tb));
                let k = (((fb - fa).abs() / FRAC_PI_2).ceil() as usize).max(1);
                let ts: Vec<f64> = (0..=k)
                    .map(|i| match i {
                        0 => ta,
                        _ if i == k => tb,
                        _ => conic_at_turn(c, fa + (fb - fa) * i as f64 / k as f64),
                    })
                    .collect();
                let mut points = Vec::with_capacity(2 * k + 1);
                let mut weights = Vec::with_capacity(2 * k + 1);
                for i in 0..k {
                    let (t0, t1) = (ts[i], ts[i + 1]);
                    let (w0, w1) = (conic_w(c, t0), conic_w(c, t1));
                    if i == 0 {
                        points.push(self.point(t0));
                        weights.push(2.0 / (1.0 + w0));
                    }
                    let m = if c < 0.0 && c * t0 * t1 < 0.0 {
                        (1.0 - c * t0 * t0 - c * t1 * t1) / (w0 * w1 - c * t0 * t1)
                    } else {
                        w0 * w1 + c * t0 * t1
                    };
                    let n = 1.0 + w0 + w1 + m;
                    // The middle control point: vertex + 2p (s₀s₁ x + (s₀ +
                    // s₁)/2 y) over the middle weight.
                    points.push(
                        vertex
                            + (x * (2.0 * p * t0 * t1)
                                + y * (p * (t0 * (1.0 + w1) + t1 * (1.0 + w0))))
                                / n,
                    );
                    weights.push(n / ((1.0 + w0) * (1.0 + w1)));
                    points.push(self.point(t1));
                    weights.push(2.0 / (1.0 + w1));
                }
                let mut knots = vec![0.0; 3];
                for i in 1..k {
                    let u = i as f64 / k as f64;
                    knots.extend([u, u]);
                }
                knots.extend([1.0; 3]);
                let n = NurbsCurve::new(2, knots, points, Some(weights))
                    .map_err(|e| KernelError::InvalidResult(format!("a conic edge: {e}")))?;
                (Curve3::Nurbs(n), 0.0, 1.0)
            }
        })
    }
}

/// `w = √(1 − c t²)` of a conic track (zero beyond half an ellipse).
pub(super) fn conic_w(c: f64, t: f64) -> f64 {
    (1.0 - c * t * t).max(0.0).sqrt()
}

/// The angle of a conic track's tangent at `t` from its `y`, toward its
/// `x`: `atan2(t, w)`, increasing with `t`, within a half turn.
pub(super) fn conic_turn(c: f64, t: f64) -> f64 {
    t.atan2(conic_w(c, t))
}

/// The parameter of a conic track where its tangent is at angle `phi`
/// ([`conic_turn`]): `t / w = tan φ`, so `t = sin φ / √(cos² φ + c sin² φ)`.
pub(super) fn conic_at_turn(c: f64, phi: f64) -> f64 {
    let (s, co) = phi.sin_cos();
    s / (co * co + c * s * s).sqrt()
}

/// Newton steps for the foot of a point on an ellipse.
const NEWTON: usize = 30;

/// An angle wrapped into `[−π/2, 3π/2)`.
fn wrap(a: f64) -> f64 {
    crate::math::wrap_angle(a)
}

/// A meridian: a surface of revolution's section by a half-plane bounded
/// by its axis, in coordinates `(ρ, h)`: distance from the axis and
/// position along it. Only its points with `ρ ≥ 0` belong to the surface
/// (for a cone, its nappe).
#[derive(Clone, Copy, Debug)]
enum Meridian {
    /// The line through `point` along the unit vector `dir`.
    Line { point: DVec2, dir: DVec2 },
    /// The circle about `(0, h)` of radius `r`.
    Circle { h: f64, r: f64 },
}

/// A common axis of two surfaces of revolution.
#[derive(Clone, Copy, Debug)]
struct Axis {
    point: DVec3,
    dir: DVec3,
}

impl Axis {
    /// The meridian of `q` about this axis (on which `q` must be a surface
    /// of revolution, up to the snapping above).
    fn meridian(&self, q: &Quadric) -> Meridian {
        let h = |p: DVec3| (p - self.point).dot(self.dir);
        match *q {
            Quadric::Plane { point, normal } => Meridian::Line {
                point: DVec2::new(0.0, (point - self.point).dot(normal) / normal.dot(self.dir)),
                dir: DVec2::X,
            },
            Quadric::Cylinder { point, radius, .. } => Meridian::Line {
                point: DVec2::new(radius, h(point)),
                dir: DVec2::Y,
            },
            Quadric::Cone {
                point,
                axis,
                radius,
                slope,
            } => {
                let s = if axis.dot(self.dir) > 0.0 {
                    slope
                } else {
                    -slope
                };
                Meridian::Line {
                    point: DVec2::new(radius, h(point)),
                    dir: DVec2::new(s, 1.0).normalize(),
                }
            }
            Quadric::Sphere { center, radius } => Meridian::Circle {
                h: h(center),
                r: radius,
            },
        }
    }

    /// The circle about the axis through the meridian point `(ρ, h)`.
    fn circle(&self, at: DVec2) -> Track {
        let x = self.dir.any_orthonormal_vector();
        Track::Circle {
            center: self.point + self.dir * at.y,
            x,
            y: self.dir.cross(x),
            radius: at.x,
        }
    }
}

/// The axis and point of a cylinder or cone.
fn axis_of(q: &Quadric) -> Option<Axis> {
    match *q {
        Quadric::Cylinder { point, axis, .. } | Quadric::Cone { point, axis, .. } => {
            Some(Axis { point, dir: axis })
        }
        _ => None,
    }
}

/// The point of the line through `p` along the unit `dir` nearest `c`.
fn nearest_on(p: DVec3, dir: DVec3, c: DVec3) -> DVec3 {
    p + dir * (c - p).dot(dir)
}

/// How near a cone's axis, in tolerances, a curve where a plane crosses
/// it may pass where faces can be: beyond the hundred the validator holds
/// a face's edges to off its surface's singular points.
const ON_AXIS: f64 = 128.0;

/// How far beyond the region, in tolerances, a ray from a cone's apex
/// outside it starts (§14): the crossings of a face's boundary with it
/// before there are no business of the boolean, as no point there is in
/// faces on both surfaces, and near the apex they would lie on both rays.
pub(super) const CLIP: f64 = 256.0;

/// How `a` and `b` meet in `region`, deciding within `eps`.
pub(super) fn meet(a: &Quadric, b: &Quadric, region: &Region, eps: f64) -> Meet {
    // Order the pair: plane, cylinder, cone, sphere.
    let rank = |q: &Quadric| match q {
        Quadric::Plane { .. } => 0,
        Quadric::Cylinder { .. } => 1,
        Quadric::Cone { .. } => 2,
        Quadric::Sphere { .. } => 3,
    };
    let (a, b) = if rank(a) <= rank(b) { (a, b) } else { (b, a) };
    let snap = eps / 16.0;
    if let Some(axis) = common_axis(a, b, region, snap) {
        return coaxial(&axis, a, b, region, eps);
    }
    match (*a, *b) {
        (
            Quadric::Plane { point, normal },
            Quadric::Plane {
                point: q,
                normal: m,
            },
        ) => planes(point, normal, q, m, region, eps),
        (
            Quadric::Plane { point, normal },
            Quadric::Cylinder {
                point: a,
                axis,
                radius,
            },
        ) => plane_cylinder(point, normal, a, axis, radius, region, eps),
        (
            Quadric::Cylinder {
                point: a1,
                axis: u1,
                radius: r1,
            },
            Quadric::Cylinder {
                point: a2,
                axis: u2,
                radius: r2,
            },
        ) => {
            if u1.cross(u2).length() * region.radius > snap {
                return Meet::Unsupported(
                    "two cylinders at an angle meet in a quartic curve (design note 0012 §10.3)",
                );
            }
            parallel_cylinders(a1, u1, r1, a2, u2, r2, region, eps)
        }
        (
            Quadric::Plane { point, normal },
            Quadric::Cone {
                point: c,
                axis,
                radius,
                slope,
            },
        ) => plane_cone(point, normal, c, axis, radius, slope, region, eps),
        _ => Meet::Unsupported("these surfaces meet in a quartic curve (design note 0012 §10.3)"),
    }
}

/// The axis about which both surfaces are surfaces of revolution, if
/// there is one (within `snap` over the region).
fn common_axis(a: &Quadric, b: &Quadric, region: &Region, snap: f64) -> Option<Axis> {
    match (*a, *b) {
        (Quadric::Plane { .. }, Quadric::Plane { .. }) => None,
        (Quadric::Plane { normal, .. }, Quadric::Sphere { center, .. }) => Some(Axis {
            point: center,
            dir: normal,
        }),
        (Quadric::Plane { point, normal }, other) => {
            // A plane square to the axis: its tilt moves the circle's
            // points off it by at most the sine times the circle's radius.
            let axis = axis_of(&other)?;
            let sine = normal.cross(axis.dir).length();
            let off = match other {
                Quadric::Cone {
                    point: c,
                    axis: v,
                    radius,
                    slope,
                } => {
                    // The cone's radius ρ₀ where the plane crosses the axis.
                    // Tilted by θ, the plane meets the cone at radius ρ₀ / (1
                    // − k tan θ cos φ) (φ about the axis from the tilt's
                    // direction): within 2 ρ₀ tan θ (1 + k) of the circle
                    // while k tan θ ≤ 1/2, the section an ellipse; nearer a
                    // ruling's direction it is not near the circle, however
                    // small ρ₀ (a plane through the apex at a slant crosses
                    // the cone in two rulings).
                    let tan = sine / normal.dot(v).abs();
                    let steep = slope * tan;
                    if steep.is_nan() || steep > 0.5 {
                        return None;
                    }
                    let h = (point - c).dot(normal) / normal.dot(v);
                    2.0 * (1.0 + slope) * (radius + slope * h).abs() * tan
                }
                Quadric::Cylinder { radius, .. } => sine * radius,
                _ => return None,
            };
            (off <= snap).then_some(axis)
        }
        (Quadric::Sphere { center, .. }, Quadric::Sphere { center: c2, .. }) => {
            let d = c2 - center;
            let dir = if d.length() > snap {
                d.normalize()
            } else {
                DVec3::Z
            };
            Some(Axis { point: center, dir })
        }
        (other, Quadric::Sphere { center, .. }) => {
            let axis = axis_of(&other)?;
            let off = center - nearest_on(axis.point, axis.dir, center);
            (off.length() <= snap).then_some(axis)
        }
        (one, other) => {
            let (p, q) = (axis_of(&one)?, axis_of(&other)?);
            let sine = p.dir.cross(q.dir).length();
            let c = region.center;
            let off = nearest_on(q.point, q.dir, c) - nearest_on(p.point, p.dir, c);
            (sine * region.radius <= snap && off.length() <= snap).then_some(p)
        }
    }
}

/// Two surfaces of revolution about `axis`: circles where their meridians
/// meet.
fn coaxial(axis: &Axis, a: &Quadric, b: &Quadric, region: &Region, eps: f64) -> Meet {
    let (ma, mb) = (axis.meridian(a), axis.meridian(b));
    let parallel = |u: DVec2, v: DVec2| u.perp_dot(v).abs() * region.radius <= eps / 16.0;
    // Points where the meridians meet, or how they otherwise relate.
    let points: Vec<DVec2> = match (ma, mb) {
        (Meridian::Line { point: p, dir: u }, Meridian::Line { point: q, dir: v }) => {
            if parallel(u, v) {
                // One line within the tolerance, or two apart, measured
                // where the region is: from the point of the first nearest
                // the region's meridian point to the second.
                let w = region.center - axis.point;
                let h = w.dot(axis.dir);
                let at = DVec2::new((w - axis.dir * h).length(), h);
                let near = p + u * (at - p).dot(u);
                let gap = (near - q).perp_dot(v).abs();
                return if gap <= eps { Meet::Same } else { Meet::Apart };
            }
            let t = (q - p).perp_dot(v) / u.perp_dot(v);
            vec![p + u * t]
        }
        (Meridian::Line { point, dir }, Meridian::Circle { h, r })
        | (Meridian::Circle { h, r }, Meridian::Line { point, dir }) => {
            let c = DVec2::new(0.0, h);
            let foot = point + dir * (c - point).dot(dir);
            let d = (foot - c).length();
            if d > r + eps {
                return Meet::Apart;
            }
            if d >= r - eps {
                return if foot.x >= -eps {
                    Meet::Touch
                } else {
                    Meet::Apart
                };
            }
            let half = (r * r - d * d).sqrt();
            vec![foot - dir * half, foot + dir * half]
        }
        (Meridian::Circle { h: h1, r: r1 }, Meridian::Circle { h: h2, r: r2 }) => {
            let d = (h2 - h1).abs();
            if d <= eps / 16.0 {
                return if (r1 - r2).abs() <= eps {
                    Meet::Same
                } else {
                    Meet::Apart
                };
            }
            if d > r1 + r2 + eps || d < (r1 - r2).abs() - eps {
                return Meet::Apart;
            }
            if d >= r1 + r2 - eps || d <= (r1 - r2).abs() + eps {
                return Meet::Touch;
            }
            let x = (d * d + r1 * r1 - r2 * r2) / (2.0 * d);
            let h = h1 + x * (h2 - h1).signum();
            vec![DVec2::new((r1 * r1 - x * x).sqrt(), h)]
        }
    };
    let mut curves = Vec::new();
    for p in points {
        if p.x < -eps {
            continue;
        }
        if p.x <= eps {
            // They meet on the axis: a point contact (a cone's apex, the
            // pole of a sphere).
            return Meet::Touch;
        }
        if p.x <= ON_AXIS * eps && region.holds(axis.point + axis.dir * p.y, p.x + eps) {
            // A circle the result's faces would hold within the tolerance
            // of a pointed cone's apex or a sphere's pole, singular points
            // of their surfaces, where faces on both can be (§15).
            return Meet::Unsupported(
                "surfaces about one axis meet in a circle within the tolerance of the axis \
                 (design note 0012 §15)",
            );
        }
        curves.push(axis.circle(p));
    }
    if curves.is_empty() {
        Meet::Apart
    } else {
        Meet::Curves(curves)
    }
}

/// Two planes: a line, or one plane, or two apart.
fn planes(p: DVec3, n: DVec3, q: DVec3, m: DVec3, region: &Region, eps: f64) -> Meet {
    let cross = n.cross(m);
    let s = cross.length();
    let c = region.center;
    if s * region.radius <= eps / 16.0 {
        let foot = c - n * (c - p).dot(n);
        return if (foot - q).dot(m).abs() <= eps {
            Meet::Same
        } else {
            Meet::Apart
        };
    }
    // The point of both planes nearest the region's centre.
    let g = n.dot(m);
    let (r1, r2) = ((p - c).dot(n), (q - c).dot(m));
    let det = s * s;
    let alpha = (r1 - g * r2) / det;
    let beta = (r2 - g * r1) / det;
    Meet::Curves(vec![Track::Line {
        origin: c + n * alpha + m * beta,
        dir: cross / s,
    }])
}

/// A plane and a cylinder off square: two lines along the axis, or an
/// ellipse.
fn plane_cylinder(
    q: DVec3,
    n: DVec3,
    a: DVec3,
    u: DVec3,
    r: f64,
    region: &Region,
    eps: f64,
) -> Meet {
    let gamma = n.dot(u);
    if gamma.abs() * region.radius <= eps / 16.0 {
        // The axis along the plane: rulings of the cylinder at the
        // plane's distance from the axis, near the region.
        let ac = nearest_on(a, u, region.center);
        let delta = (ac - q).dot(n);
        if delta.abs() > r + eps {
            return Meet::Apart;
        }
        if delta.abs() >= r - eps {
            return Meet::Touch;
        }
        let across = (n - u * gamma).normalize();
        let side = u.cross(across);
        let half = (r * r - delta * delta).sqrt();
        let foot = ac - across * delta;
        return Meet::Curves(vec![
            Track::Line {
                origin: foot - side * half,
                dir: u,
            },
            Track::Line {
                origin: foot + side * half,
                dir: u,
            },
        ]);
    }
    // The ellipse: centred where the axis crosses the plane, its minor
    // semi-axis the radius across the axis, its major one the radius over
    // the cosine of the tilt, along the axis's direction in the plane.
    // y across the axis, then x = y × n = u − n γ (normalized) along the
    // axis's direction in the plane: both in the plane to rounding (y made
    // square to n once more, as the rounding of a short cross product is
    // not), even when the plane is nearly square to the axis.
    let center = a + u * ((q - a).dot(n) / gamma);
    let y = n.cross(u).normalize();
    let y = (y - n * y.dot(n)).normalize();
    let x = y.cross(n);
    Meet::Curves(vec![Track::Ellipse {
        center,
        x,
        y,
        a: r / gamma.abs(),
        b: r,
    }])
}

/// A plane (through `q`, unit normal `n`) across a cone (reference point
/// `a`, unit axis `u`, radius `r` there, slope `k`), not square to its
/// axis: a conic, or two rulings through the apex (design note 0012 §14).
#[allow(clippy::too_many_arguments)]
fn plane_cone(
    q: DVec3,
    n: DVec3,
    a: DVec3,
    u: DVec3,
    r: f64,
    k: f64,
    region: &Region,
    eps: f64,
) -> Meet {
    let snap = eps / 16.0;
    // The plane's frame: x along the axis's direction in the plane (made
    // square to n once more, as the rounding of a short difference is
    // not), y across the axis.
    let gamma = n.dot(u);
    let along = u - n * gamma;
    let g = along.length();
    if g == 0.0 {
        // Square to the axis: circles, found as coaxial.
        return Meet::Unsupported(
            "a plane square to a cone's axis was not taken as coaxial (design note 0012 §15)",
        );
    }
    let x = along / g;
    let x = (x - n * x.dot(n)).normalize();
    let y = n.cross(x);
    let d = (a - q).dot(n);
    let o = a - n * d;
    // The section A ξ² + η² − 2Bξ + C = 0 about o, and σ = k δ, δ the
    // apex's distance from the plane: B² − AC = σ².
    let big_a = (gamma.abs() - g * k) * (gamma.abs() + g * k);
    let big_b = g * (r * k - (1.0 + k * k) * gamma * d);
    let big_c = (d * g - r + k * gamma * d) * (d * g + r - k * gamma * d);
    let sigma = d * k - gamma * r;
    let apex = a - u * (r / k);
    let delta = sigma / k;
    // Whether the plane is within the tolerance of a plane tangent to the
    // cone along a ruling, over the region: their angle is |g − cos α| /
    // |γ|, with |g − cos α| = |A| cos² α / (cos α + g).
    let cos = 1.0 / (1.0 + k * k).sqrt();
    let tilt = big_a.abs() * cos * cos / ((cos + g) * gamma.abs());
    let tangent = tilt * ((region.center - apex).length() + region.radius) <= eps;
    if delta.abs() <= eps {
        // Through the apex within the tolerance. Along the ruling nearest
        // the plane, toward x, the cone leaves the plane by μ = |cos α |γ| −
        // sin α g| = cos α |A| / (|γ| + k g) per unit of length from the
        // apex (α its half-angle): it stays within the tolerance of the plane
        // (A > 0), or the plane cuts a sliver from it between two rulings
        // thinner than the tolerance (A < 0), out to `reach` from the apex,
        // and further the more nearly the plane is tangent to it. Where faces
        // on both can come that near the apex, the tolerance cannot decide
        // how they meet (a frustum's top within `reach` of the apex: found
        // by review, as two solids a nanometre apart and as an inconsistent
        // cut).
        let mu = cos * big_a.abs() / (gamma.abs() + k * g);
        let reach = ((eps + delta.abs()) / mu).max(eps);
        if region.near_apex(apex, reach, reach * cos) {
            return if big_a < 0.0 && !tangent {
                Meet::Unsupported(
                    "a plane through a cone's apex meets it in two rulings (design note 0012 §14)",
                )
            } else {
                Meet::Touch
            };
        }
    }
    if delta.abs() <= snap {
        return through_apex(
            n, x, y, a, u, r, k, gamma, g, big_a, tangent, apex, region, eps,
        );
    }
    // Within `tip` of the apex the cone is within 128 tolerances of its
    // axis. The curve comes nearest the axis by the apex, k |δ| / (|γ| +
    // k g) from it (k |δ| for a circle, |δ| for a plane along the axis); a
    // curve nearer than 128 tolerances, crossing the cone where faces can
    // reach the apex, turns there: a face of a pointed cone would get an
    // edge passing its singular point within the tolerance (§14).
    let tip = ON_AXIS * eps * (1.0 + k * k).sqrt() / k;
    if k * delta.abs() <= ON_AXIS * eps * (gamma.abs() + k * g)
        && region.near_apex(apex, tip, tip)
        && (big_a < 0.0 || sigma * gamma < 0.0)
    {
        return Meet::Unsupported(
            "a plane passes so near a cone's apex that their curve comes within the tolerance \
             of its axis (design note 0012 §15)",
        );
    }
    // The ruling nearest the plane: from the apex along r̂ = cos α u −
    // sign(γ) sin α n⊥ / g (n⊥ = n − γ u), the cone leaves the plane at
    // m = sign(γ) cos α (|γ| − k g) per unit of length, f(s) = δ + s m,
    // every other ruling at a rate further from zero on the side of γ.
    // With γ δ > 0 the cone lies on the apex's side of the plane wherever
    // this ruling does (the others leave toward that side, A ≥ 0, or come
    // toward the plane more slowly than it, A < 0), nearest the plane on
    // it: within the tolerance of the plane where faces on both can be,
    // without crossing, it touches (a frustum's rim within the tolerance
    // of a plane nearly tangent along this ruling, found by review as two
    // solids a few hundredths of a micrometre apart).
    if gamma * delta > 0.0 {
        let sign = gamma.signum();
        let ruling = u * cos - (n - u * gamma) * (sign * k * cos / g);
        let m = sign * cos * (gamma.abs() - k * g);
        if let Some((s0, s1)) = region.ray(apex, ruling, (region.lift / cos).max(0.0)) {
            let (f0, f1) = (delta + s0 * m, delta + s1 * m);
            if f0 * delta > 0.0 && f1 * delta > 0.0 && f0.abs().min(f1.abs()) <= eps {
                return Meet::Touch;
            }
        }
    }
    // The disk where the plane crosses the region.
    let w = region.center - o;
    let off = w.dot(n);
    if off.abs() > region.radius + eps {
        return Meet::Apart;
    }
    let r0 = (region.radius * region.radius - off * off).max(0.0).sqrt() + eps;
    let xi0 = w.dot(x);
    let p = sigma.abs();
    if big_a >= 0.0 && sigma * gamma >= 0.0 {
        // An ellipse or a parabola on the other nappe.
        return Meet::Apart;
    }
    if big_a > 0.0 && (big_b - big_a * xi0).abs() <= 16.0 * big_a * r0 {
        // An ellipse whose centre is near the region: closed.
        let center = o + x * (big_b / big_a);
        let (ea, eb) = (p / big_a, p / big_a.sqrt());
        let c = center - region.center;
        if c.dot(x).abs() > ea + r0 + region.radius || c.dot(y).abs() > eb + r0 + region.radius {
            return Meet::Apart;
        }
        return Meet::Curves(vec![Track::Ellipse {
            center,
            x,
            y,
            a: ea,
            b: eb,
        }]);
    }
    // The roots of A ξ² − 2Bξ + C, each in the form without cancellation:
    // the vertices of a conic's branch or of an ellipse.
    let greater = || {
        if big_b >= 0.0 {
            (big_b + p) / big_a
        } else {
            big_c / (big_b - p)
        }
    };
    let lesser = || {
        if big_b >= 0.0 {
            big_c / (big_b + p)
        } else {
            (big_b - p) / big_a
        }
    };
    // The vertex (by its ξ) and the way the curve opens from it.
    let (xi, opens) = if big_a > 0.0 && big_b - big_a * xi0 < 0.0 {
        // The half of an ellipse at its vertex of greater ξ, the centre
        // beyond the region the other way: opening along −x.
        (greater(), -1.0)
    } else {
        // A hyperbola's branch on the nappe, a parabola, or the half of an
        // ellipse at its vertex of lesser ξ: opening along +x.
        (lesser(), 1.0)
    };
    if !xi.is_finite() || (xi0 - xi) * opens < -r0 {
        // The curve lies beyond the region.
        return Meet::Apart;
    }
    let vertex = o + x * xi;
    // Points are found from the vertex: refuse where its distance could
    // make their rounding a sixty-fourth of the tolerance (§14).
    let reach = vertex.abs().max_element() + (vertex - region.center).length() + region.radius;
    if 8.0 * f64::EPSILON * reach > eps / 64.0 {
        return Meet::Unsupported(
            "a plane across a cone meets it in a conic whose vertex is too far away to compute \
             it within the tolerance (design note 0012 §14)",
        );
    }
    let mut curves = vec![Track::Conic {
        vertex,
        x: x * opens,
        y,
        p,
        c: big_a,
    }];
    if big_a > 0.0 {
        // The ellipse's other half, from its other vertex. It lies beyond
        // the region, but a face's boundary can cross it where the faces'
        // boxes reach further (a small tool on a large cone), and every such
        // crossing must lie on a track. Where its vertex is too far to
        // compute it within the tolerance, hundreds of kilometres away, no
        // face reaches it.
        let far = o + x * if opens > 0.0 { greater() } else { lesser() };
        let reach = far.abs().max_element() + (far - region.center).length() + region.radius;
        if far.is_finite() && 8.0 * f64::EPSILON * reach <= eps / 64.0 {
            curves.push(Track::Conic {
                vertex: far,
                x: x * -opens,
                y,
                p,
                c: big_a,
            });
        }
    }
    Meet::Curves(curves)
}

/// A plane within `eps / 16` of a cone's apex, which lies outside the
/// region's box (§14), if not always its ball: two rulings where it crosses the cone (`A < 0`), rays
/// away from the apex (from where they come within [`CLIP`] tolerances of
/// the region, or from the apex where it is that near), exact on the cone
/// at the region's height and within `eps / 16` of the plane; a tangent contact where it is within the tolerance of the plane
/// tangent along the ruling toward `x` over the region (`tangent`); else
/// apart.
#[allow(clippy::too_many_arguments)]
fn through_apex(
    n: DVec3,
    x: DVec3,
    y: DVec3,
    a: DVec3,
    u: DVec3,
    r: f64,
    k: f64,
    gamma: f64,
    g: f64,
    big_a: f64,
    tangent: bool,
    apex: DVec3,
    region: &Region,
    eps: f64,
) -> Meet {
    if tangent {
        // The ruling toward x, if it reaches the region.
        let w = region.center - apex;
        let s = w.dot(x).max(0.0);
        return if (w - x * s).length() <= region.radius + eps {
            Meet::Touch
        } else {
            Meet::Apart
        };
    }
    if big_a > 0.0 {
        // Only the apex.
        return Meet::Apart;
    }
    // The rulings square to n: unit (cos/g) x ± (√(−A) cos/g) y, and the
    // unit e from the axis toward each, ((γ²/g) x ± (√(−A)/g) y − γ n) / k,
    // at the region's height h; at as far beyond the apex where the
    // region's centre lies below it (beside an oblique axis), as the cone
    // is the one nappe.
    let h = (region.center - a).dot(u);
    let h = if r + k * h < 0.0 { -2.0 * r / k - h } else { h };
    let rho = r + k * h;
    let spread = (-big_a).sqrt();
    let mut curves = Vec::with_capacity(2);
    let mut near = false;
    for side in [-1.0, 1.0] {
        let dir = (x + y * (side * spread)).normalize();
        let e = (x * (gamma * gamma / g) + y * (side * spread / g) - n * gamma) / k;
        let origin = a + u * h + e * rho;
        // Away from the apex, which lies behind the region.
        let start = (apex - origin).dot(dir);
        let (dir, start) = if start > 0.0 {
            (-dir, -start)
        } else {
            (dir, start)
        };
        let w = region.center - origin;
        let along = w.dot(dir);
        near |= (w - dir * along.max(start)).length() <= region.radius + eps;
        // From where it comes within CLIP tolerances of the region; where it
        // stays farther, from infinity (no stretch of it matters).
        let reach = region.radius + CLIP * eps;
        let disc = reach * reach - w.length_squared() + along * along;
        let entry = if disc >= 0.0 {
            along - disc.sqrt()
        } else {
            f64::INFINITY
        };
        curves.push(Track::Ray {
            origin,
            dir,
            start: entry.max(start),
            from_apex: entry <= start,
        });
    }
    // Both rays where either comes near the region: a face's boundary can
    // cross the other beyond it, where the faces' boxes reach further, and
    // every such crossing must lie on a track.
    if near {
        Meet::Curves(curves)
    } else {
        Meet::Apart
    }
}

/// Two cylinders with parallel axes (within the snap), not coaxial: two
/// lines, or touching, or apart. The lines are rulings of the first.
#[allow(clippy::too_many_arguments)]
fn parallel_cylinders(
    a1: DVec3,
    u: DVec3,
    r1: f64,
    a2: DVec3,
    u2: DVec3,
    r2: f64,
    region: &Region,
    eps: f64,
) -> Meet {
    let c = region.center;
    let p1 = nearest_on(a1, u, c);
    let off = nearest_on(a2, u2, c) - p1;
    let off = off - u * off.dot(u);
    let d = off.length();
    if d > r1 + r2 + eps || d < (r1 - r2).abs() - eps {
        return Meet::Apart;
    }
    if d >= r1 + r2 - eps || d <= (r1 - r2).abs() + eps {
        return Meet::Touch;
    }
    let e = off / d;
    let f = u.cross(e);
    let x = (d * d + r1 * r1 - r2 * r2) / (2.0 * d);
    let y = (r1 * r1 - x * x).sqrt();
    Meet::Curves(vec![
        Track::Line {
            origin: p1 + e * x - f * y,
            dir: u,
        },
        Track::Line {
            origin: p1 + e * x + f * y,
            dir: u,
        },
    ])
}

#[cfg(test)]
mod tests {
    use std::f64::consts::{PI, TAU};

    use super::super::testing::Rng;
    use super::*;

    const EPS: f64 = 1e-7;

    /// The region of the box holding the ball about `center`.
    fn ball(center: DVec3, radius: f64) -> Region {
        let r = DVec3::splat(radius / 3f64.sqrt());
        Region::of_box(center - r, center + r)
    }

    /// An `(n + 1)²` grid of points of `q` around the region (`None` off
    /// a cone's nappe): plane coordinates about the region's centre, or an
    /// angle about the axis and a height, or longitude and latitude.
    fn grid(q: &Quadric, region: &Region, n: usize) -> Vec<Vec<Option<DVec3>>> {
        let frame = |n: DVec3| {
            let x = n.any_orthonormal_vector();
            (x, n.cross(x))
        };
        let f = |i: usize| i as f64 / n as f64;
        (0..=n)
            .map(|i| {
                (0..=n)
                    .map(|j| {
                        let (s, t) = (f(i), f(j));
                        match *q {
                            Quadric::Plane { point, normal } => {
                                let (x, y) = frame(normal);
                                let c =
                                    region.center - normal * (region.center - point).dot(normal);
                                let r = region.radius;
                                Some(c + x * (r * (2.0 * s - 1.0)) + y * (r * (2.0 * t - 1.0)))
                            }
                            Quadric::Cylinder {
                                point,
                                axis,
                                radius,
                            } => {
                                let (x, y) = frame(axis);
                                let u = TAU * s;
                                Some(
                                    point
                                        + axis * (4.0 * t - 2.0)
                                        + (x * u.cos() + y * u.sin()) * radius,
                                )
                            }
                            Quadric::Cone {
                                point,
                                axis,
                                radius,
                                slope,
                            } => {
                                let (x, y) = frame(axis);
                                let (u, v) = (TAU * s, 4.0 * t - 2.0);
                                let rho = radius + slope * v;
                                (rho >= 0.0)
                                    .then(|| point + axis * v + (x * u.cos() + y * u.sin()) * rho)
                            }
                            Quadric::Sphere { center, radius } => {
                                let (u, v) = (TAU * s, PI * (t - 0.5));
                                Some(
                                    center
                                        + (DVec3::X * u.cos() + DVec3::Y * u.sin())
                                            * (radius * v.cos())
                                        + DVec3::Z * (radius * v.sin()),
                                )
                            }
                        }
                    })
                    .collect()
            })
            .collect()
    }

    /// Neighbouring grid points, along both directions.
    fn chords(g: &[Vec<Option<DVec3>>]) -> Vec<(DVec3, DVec3)> {
        let n = g.len();
        let mut out = Vec::new();
        for i in 0..n {
            for j in 0..n {
                let p = g[i][j];
                for q in [g.get(i + 1).map(|r| r[j]), g[i].get(j + 1).copied()] {
                    if let (Some(p), Some(Some(q))) = (p, q) {
                        out.push((p, q));
                    }
                }
            }
        }
        out
    }

    /// Parameters of a track covering its part near the region: a whole
    /// closed curve; a line's stretch through the region (a ray's from
    /// its apex at most); a conic's
    /// stretch across it (by its position across the axis), within half
    /// an ellipse.
    fn span(t: &Track, region: &Region) -> (f64, f64) {
        match *t {
            Track::Circle { .. } | Track::Ellipse { .. } => (0.0, TAU),
            Track::Line { .. } => {
                let c = t.param(region.center);
                (c - region.radius, c + region.radius)
            }
            Track::Ray { start, .. } => {
                let c = t.param(region.center);
                ((c - region.radius).max(start), c + region.radius)
            }
            Track::Conic { p, c, .. } => {
                let m = t.param(region.center);
                let r = 2.0 * region.radius / p;
                let limit = if c > 0.0 {
                    0.999 / c.sqrt()
                } else {
                    f64::INFINITY
                };
                ((m - r).max(-limit), (m + r).min(limit))
            }
        }
    }

    /// The distance from `p` to the curve, for points near it.
    fn off_track(t: &Track, p: DVec3) -> f64 {
        match *t {
            Track::Line { origin, dir } => {
                let w = p - origin;
                (w - dir * w.dot(dir)).length()
            }
            Track::Ray {
                origin, dir, start, ..
            } => {
                let w = p - origin;
                (w - dir * w.dot(dir).max(start)).length()
            }
            _ => {
                // A few Newton steps on the squared distance, from the
                // closed-form guess.
                let mut s = t.param(p);
                for _ in 0..8 {
                    let d = t.point(s) - p;
                    let (d1, d2) = (t.deriv(s), t.second(s));
                    let g = d1.dot(d1) + d.dot(d2);
                    if g > 0.0 {
                        s -= d.dot(d1) / g;
                    }
                }
                (t.point(s) - p).length()
            }
        }
    }

    /// Checks a meeting: every curve lies on both surfaces, and every
    /// crossing of the second surface found by sampling the first around
    /// the region is on a curve. Returns the number of crossings sampled.
    fn check(a: &Quadric, b: &Quadric, region: &Region, case: &str) -> usize {
        let m = meet(a, b, region, EPS);
        let Meet::Curves(tracks) = &m else {
            panic!("{case}: {m:?}");
        };
        let mag = 1.0 + a.magnitude() + b.magnitude() + region.center.abs().max_element();
        for t in tracks.iter().filter(|t| t.ends().0 < f64::INFINITY) {
            let (lo, hi) = span(t, region);
            for i in 0..=64 {
                let p = t.point(lo + (hi - lo) * i as f64 / 64.0);
                if (p - region.center).length() > 2.0 * region.radius {
                    continue;
                }
                let (da, db) = (a.distance(p), b.distance(p));
                assert!(
                    da.abs() <= 1e-12 * mag && db.abs() <= 1e-12 * mag,
                    "{case}: {da} {db}"
                );
            }
        }
        let mut sampled = 0;
        for (p, q) in chords(&grid(a, region, 80)) {
            let (dp, dq) = (b.distance(p), b.distance(q));
            if dp * dq >= 0.0 || dp.abs() < 1e-9 || dq.abs() < 1e-9 {
                continue;
            }
            // The crossing along the chord, which is within the chord's
            // sag of the first surface.
            let (mut lo, mut hi) = (p, q);
            for _ in 0..60 {
                let mid = (lo + hi) * 0.5;
                if (b.distance(mid) < 0.0) == (dp < 0.0) {
                    lo = mid;
                } else {
                    hi = mid;
                }
            }
            let x = (lo + hi) * 0.5;
            // Off the first surface by the chord's sag, so off the curve by
            // about that over the sine of the angle the surfaces meet at.
            let sine = normal(a, x).cross(normal(b, x)).length();
            if sine < 1e-3 {
                continue;
            }
            let sag = (p - q).length_squared() + a.distance(x).abs();
            let near = tracks
                .iter()
                .map(|t| off_track(t, x))
                .fold(f64::MAX, f64::min);
            assert!(
                near <= 4.0 * sag / sine + 1e-9,
                "{case}: a crossing {near} from every curve"
            );
            sampled += 1;
        }
        sampled
    }

    /// The direction in which the distance from `q` grows fastest at `p`.
    fn normal(q: &Quadric, p: DVec3) -> DVec3 {
        let h = 1e-6;
        let d = |v: DVec3| q.distance(p + v * h) - q.distance(p - v * h);
        DVec3::new(d(DVec3::X), d(DVec3::Y), d(DVec3::Z)).normalize()
    }

    /// Sign changes of the second's distance over a sampling of the first.
    fn sign_changes(a: &Quadric, b: &Quadric, region: &Region) -> usize {
        chords(&grid(a, region, 40))
            .iter()
            .filter(|(p, q)| b.distance(*p) * b.distance(*q) < 0.0)
            .count()
    }

    #[test]
    fn curves_lie_on_both_surfaces_and_none_is_missed() {
        let mut rng = Rng(0x3EE7_0001);
        let mut sampled = 0;
        let mut outcomes = std::collections::BTreeMap::new();
        for case in 0..400 {
            let (p, u) = (rng.point(0.5), rng.unit());
            let region = ball(p + rng.point(0.2), 2.0);
            let r = rng.range(0.2, 0.8);
            let flip = |rng: &mut Rng| if rng.next() < 0.5 { u } else { -u };
            let across = |rng: &mut Rng| u.cross(rng.unit()).normalize();
            let cyl = Quadric::Cylinder {
                point: p,
                axis: u,
                radius: r,
            };
            let cone = Quadric::Cone {
                point: p + u * rng.range(-0.3, 0.3),
                axis: flip(&mut rng),
                radius: rng.range(0.2, 0.6),
                slope: 10f64.powf(rng.range(-1.0, 0.7)),
            };
            let sphere = Quadric::Sphere {
                center: p + u * rng.range(-0.4, 0.4),
                radius: rng.range(0.3, 1.0),
            };
            let square = Quadric::Plane {
                point: p + u * rng.range(-0.3, 0.3) + across(&mut rng) * 0.3,
                normal: flip(&mut rng),
            };
            let oblique = Quadric::Plane {
                point: p + rng.point(0.2),
                normal: (u * rng.range(0.1, 0.9) + across(&mut rng)).normalize(),
            };
            let along = Quadric::Plane {
                point: p + across(&mut rng) * rng.range(0.0, 0.9 * r),
                normal: across(&mut rng),
            };
            let beside = Quadric::Cylinder {
                point: p + across(&mut rng) * rng.range(0.1, 1.2) + u * rng.range(-1.0, 1.0),
                axis: flip(&mut rng),
                radius: rng.range(0.3, 0.8),
            };
            let other_sphere = Quadric::Sphere {
                center: p + rng.point(0.6),
                radius: rng.range(0.3, 1.0),
            };
            let coaxial_cone = Quadric::Cone {
                point: p + u * rng.range(-0.5, 0.5),
                axis: flip(&mut rng),
                radius: rng.range(0.2, 0.6),
                slope: 10f64.powf(rng.range(-1.0, 0.7)),
            };
            let plane2 = Quadric::Plane {
                point: p + rng.point(0.3),
                normal: rng.unit(),
            };
            let pairs: Vec<(&str, Quadric, Quadric)> = vec![
                ("planes", oblique, plane2),
                ("plane across cylinder", oblique, cyl),
                ("plane along cylinder", along, cyl),
                ("plane square to cylinder", square, cyl),
                ("plane square to cone", square, cone),
                ("plane and sphere", plane2, other_sphere),
                ("cylinders side by side", cyl, beside),
                ("coaxial cylinder and cone", cyl, cone),
                ("coaxial cones", cone, coaxial_cone),
                ("sphere on a cylinder's axis", cyl, sphere),
                ("sphere on a cone's axis", cone, sphere),
                ("spheres", sphere, other_sphere),
            ];
            for (name, a, b) in pairs {
                let m = meet(&a, &b, &region, EPS);
                let outcome = match &m {
                    Meet::Curves(_) => {
                        sampled += check(&a, &b, &region, &format!("{name} {case}"));
                        // The pair the other way round meets in the same
                        // curves (two of a kind run along them either way).
                        if std::mem::discriminant(&a) != std::mem::discriminant(&b) {
                            assert_eq!(meet(&b, &a, &region, EPS), m);
                        }
                        "curves"
                    }
                    Meet::Apart => {
                        // Nothing crosses where either is sampled.
                        assert_eq!(sign_changes(&a, &b, &region), 0, "{name} {case}");
                        assert_eq!(sign_changes(&b, &a, &region), 0, "{name} {case}");
                        "apart"
                    }
                    other => panic!("{name} {case}: {other:?}"),
                };
                *outcomes.entry((name, outcome)).or_insert(0) += 1;
            }
        }
        assert!(sampled > 100_000, "{sampled}");
        // Both outcomes are exercised where both occur.
        for name in [
            "plane along cylinder",
            "plane and sphere",
            "spheres",
            "coaxial cones",
        ] {
            assert!(
                outcomes.get(&(name, "curves")).copied().unwrap_or(0) > 20,
                "{outcomes:?}"
            );
        }
    }

    #[test]
    fn touching_coincident_apart_and_unsupported_pairs() {
        let region = ball(DVec3::new(0.1, 0.2, 0.3), 1.5);
        let (o, u) = (DVec3::new(0.1, -0.1, 0.2), DVec3::new(1.0, 2.0, 2.0) / 3.0);
        let x = u.any_orthonormal_vector();
        let cyl = |r: f64, at: DVec3| Quadric::Cylinder {
            point: at,
            axis: u,
            radius: r,
        };
        let plane = |at: DVec3, n: DVec3| Quadric::Plane {
            point: at,
            normal: n,
        };
        let sphere = |c: DVec3, r: f64| Quadric::Sphere {
            center: c,
            radius: r,
        };
        let cone = Quadric::Cone {
            point: o,
            axis: u,
            radius: 0.3,
            slope: 0.5,
        };
        let apex = o - u * 0.6;
        let e = EPS;
        let cases: Vec<(&str, Quadric, Quadric, Meet)> = vec![
            (
                "plane on a cylinder",
                plane(o + x * (0.5 + 0.5 * e), x),
                cyl(0.5, o),
                Meet::Touch,
            ),
            (
                "plane just off a cylinder",
                plane(o + x * (0.5 + 2.0 * e), x),
                cyl(0.5, o),
                Meet::Apart,
            ),
            (
                "cylinders side by side",
                cyl(0.3, o),
                cyl(0.2, o + x * (0.5 - 0.5 * e)),
                Meet::Touch,
            ),
            (
                "cylinder inside one",
                cyl(0.5, o),
                cyl(0.2, o + x * (0.3 + 0.5 * e)),
                Meet::Touch,
            ),
            (
                "concentric cylinders",
                cyl(0.5, o),
                cyl(0.5 + 2.0 * e, o),
                Meet::Apart,
            ),
            (
                "one cylinder",
                cyl(0.5, o),
                cyl(0.5 + 0.5 * e, o + x * (e / 64.0)),
                Meet::Same,
            ),
            (
                "plane on a sphere",
                plane(o + x * (0.4 - 0.5 * e), x),
                sphere(o, 0.4),
                Meet::Touch,
            ),
            (
                "spheres touching",
                sphere(o, 0.4),
                sphere(o + x * (0.7 + 0.5 * e), 0.3),
                Meet::Touch,
            ),
            (
                "sphere inside one",
                sphere(o, 0.4),
                sphere(o + x * (0.1 - 0.5 * e), 0.3),
                Meet::Touch,
            ),
            (
                "one sphere",
                sphere(o, 0.4),
                sphere(o + x * (e / 64.0), 0.4 + 0.5 * e),
                Meet::Same,
            ),
            (
                "sphere in a cylinder",
                cyl(0.4, o),
                sphere(o + u * 0.3, 0.4 + 0.5 * e),
                Meet::Touch,
            ),
            (
                "plane through a cone's apex",
                plane(apex, u),
                cone,
                Meet::Touch,
            ),
            (
                "plane behind a cone's apex",
                plane(apex - u * 0.1, -u),
                cone,
                Meet::Apart,
            ),
            (
                "one plane",
                plane(o, u),
                plane(o + x + u * (0.5 * e), -u),
                Meet::Same,
            ),
            (
                "parallel planes",
                plane(o, u),
                plane(o + x + u * (2.0 * e), u),
                Meet::Apart,
            ),
            (
                "one cone",
                cone,
                Quadric::Cone {
                    point: o + u * 0.2,
                    axis: u,
                    radius: 0.4 + 0.5 * e,
                    slope: 0.5,
                },
                Meet::Same,
            ),
            (
                // A sphere centred on the axis whose radius is its distance
                // from the cone's generators.
                "sphere in a cone",
                cone,
                sphere(o + u * 0.4, 0.5 / (1.25f64).sqrt() + 0.5 * e),
                Meet::Touch,
            ),
            (
                "cylinders at an angle",
                cyl(0.3, o),
                Quadric::Cylinder {
                    point: o,
                    axis: x,
                    radius: 0.2,
                },
                Meet::Unsupported(
                    "two cylinders at an angle meet in a quartic curve (design note 0012 §10.3)",
                ),
            ),
            (
                "plane along a cone's axis through its apex",
                plane(apex, x),
                cone,
                Meet::Unsupported(
                    "a plane through a cone's apex meets it in two rulings (design note 0012 §14)",
                ),
            ),
            (
                "plane steeper than the rulings through the apex",
                plane(apex + x * (0.5 * e), (u + x * 0.3).normalize()),
                cone,
                Meet::Touch,
            ),
            (
                // Normal to the ruling toward x, square to the axis's plane
                // through it: cos α x − sin α u with tan α = 0.5.
                "plane tangent along a ruling",
                plane(apex - u * (0.5 * e), (x * 2.0 - u).normalize()),
                cone,
                Meet::Touch,
            ),
            (
                "plane steeper than the rulings behind the apex",
                plane(apex - u * 0.05, (u + x * 0.3).normalize()),
                cone,
                Meet::Apart,
            ),
            (
                "sphere off a cylinder's axis",
                cyl(0.3, o),
                sphere(o + x * 0.2, 0.3),
                Meet::Unsupported(
                    "these surfaces meet in a quartic curve (design note 0012 §10.3)",
                ),
            ),
            (
                "cone beside a cylinder",
                cyl(0.3, o + x * 0.5),
                cone,
                Meet::Unsupported(
                    "these surfaces meet in a quartic curve (design note 0012 §10.3)",
                ),
            ),
        ];
        for (name, a, b, want) in cases {
            assert_eq!(meet(&a, &b, &region, e), want, "{name}");
            assert_eq!(
                meet(&b, &a, &region, e),
                want,
                "{name}, the other way round"
            );
        }
        // Just inside the touching distance, the surfaces cross.
        let m = meet(&plane(o + x * (0.5 - 2.0 * e), x), &cyl(0.5, o), &region, e);
        assert!(matches!(&m, Meet::Curves(t) if t.len() == 2), "{m:?}");
    }

    #[test]
    fn near_special_configurations_snap_within_a_sixteenth() {
        let mut rng = Rng(0x3EE7_0002);
        let region = ball(DVec3::ZERO, 1.0);
        // Points of the curves within the region, within `limit` of both.
        let within = |a: &Quadric, b: &Quadric, m: &Meet, limit: f64, case: &str| {
            let Meet::Curves(tracks) = m else {
                panic!("{m:?}");
            };
            for t in tracks {
                for i in 0..=200 {
                    let p = if t.closed() {
                        t.point(TAU * i as f64 / 200.0)
                    } else {
                        let c = t.param(region.center);
                        t.point(c - 1.0 + 2.0 * i as f64 / 200.0)
                    };
                    if (p - region.center).length() > region.radius {
                        continue;
                    }
                    let (da, db) = (a.distance(p), b.distance(p));
                    assert!(da.abs() <= limit && db.abs() <= limit, "{case}: {da} {db}");
                }
            }
        };
        for _ in 0..200 {
            let u = rng.unit();
            let w = u.cross(rng.unit()).normalize();
            let r = rng.range(0.2, 0.8);
            let at = rng.point(0.1);
            let cyl = Quadric::Cylinder {
                point: at,
                axis: u,
                radius: r,
            };
            // A plane nearly square to the axis: a circle, or an ellipse
            // once the tilt matters.
            for (sine, circle) in [(EPS / 32.0 / r, true), (4.0 * EPS / r, false)] {
                let n = (u + w * sine).normalize();
                let p = Quadric::Plane {
                    point: rng.point(0.1),
                    normal: n,
                };
                let m = meet(&p, &cyl, &region, EPS);
                let Meet::Curves(t) = &m else { panic!("{m:?}") };
                assert_eq!(matches!(t[0], Track::Circle { .. }), circle);
                within(
                    &p,
                    &cyl,
                    &m,
                    if circle { EPS / 16.0 } else { 1e-13 },
                    "square",
                );
            }
            // A plane nearly along the axis: two lines.
            let p = Quadric::Plane {
                point: at + w * rng.range(-0.9, 0.9) * r,
                normal: (w + u * (EPS / 32.0)).normalize(),
            };
            let m = meet(&p, &cyl, &region, EPS);
            within(&p, &cyl, &m, EPS / 16.0, "along");
            // Cylinders with nearly parallel axes: two lines.
            let r2 = rng.range(0.3, 0.8);
            let other = Quadric::Cylinder {
                point: at + w * ((r - r2).abs() + r + r2) * 0.5,
                axis: (u + w * (EPS / 32.0)).normalize(),
                radius: r2,
            };
            let m = meet(&cyl, &other, &region, EPS);
            within(&cyl, &other, &m, EPS / 16.0, "cylinders");
        }
    }

    #[test]
    fn parameters_and_elliptical_edges_are_exact() {
        let mut rng = Rng(0x3EE7_0003);
        for _ in 0..200 {
            let u = rng.unit();
            let n = (u * rng.range(0.2, 0.95) + u.cross(rng.unit()).normalize()).normalize();
            let (plane, cyl) = (
                Quadric::Plane {
                    point: rng.point(0.3),
                    normal: n,
                },
                Quadric::Cylinder {
                    point: rng.point(0.3),
                    axis: u,
                    radius: rng.range(0.1, 1.0),
                },
            );
            let region = ball(DVec3::ZERO, 1.0);
            let Meet::Curves(t) = meet(&plane, &cyl, &region, EPS) else {
                panic!()
            };
            let e = t[0];
            let tracks = [
                e,
                Track::Circle {
                    center: rng.point(1.0),
                    x: u.any_orthonormal_vector(),
                    y: u.cross(u.any_orthonormal_vector()),
                    radius: 0.7,
                },
                Track::Line {
                    origin: rng.point(1.0),
                    dir: u,
                },
            ];
            for t in tracks {
                for i in 0..20 {
                    let s = -FRAC_PI_2 + TAU * i as f64 / 20.0;
                    let q = t.param(t.point(s));
                    assert!((t.point(q) - t.point(s)).length() < 1e-14 * 8.0);
                    if t.closed() {
                        assert!((q - s).abs() < 1e-13);
                    }
                }
            }
            let ta = rng.range(-PI, PI);
            let tb = ta + rng.range(0.1, TAU - 0.01);
            let (curve, t0, t1) = e.edge(ta, tb).unwrap();
            assert!((curve.point(t0) - e.point(ta)).length() < 1e-14 * 8.0);
            assert!((curve.point(t1) - e.point(tb)).length() < 1e-14 * 8.0);
            for i in 0..=100 {
                let p = curve.point(t0 + (t1 - t0) * i as f64 / 100.0);
                let mag = 1.0 + cyl.magnitude() + plane.magnitude();
                assert!(plane.distance(p).abs() < 1e-13 * mag);
                assert!(cyl.distance(p).abs() < 1e-13 * mag);
            }
        }
    }

    /// A cone about `u` through `p`, its axis either way, with a plane
    /// whose normal is at angle `beta` from the axis, and a region around
    /// them.
    fn cone_and_plane(rng: &mut Rng, beta: f64, radius: f64) -> (Quadric, Quadric, Region) {
        let (p, u) = (rng.point(0.5), rng.unit());
        let k = 10f64.powf(rng.range(-1.0, 0.7));
        let axis = if rng.next() < 0.5 { u } else { -u };
        let cone = Quadric::Cone {
            point: p + u * rng.range(-0.3, 0.3),
            axis,
            radius: rng.range(0.2, 0.6),
            slope: k,
        };
        let across = u.cross(rng.unit()).normalize();
        let plane = Quadric::Plane {
            point: p + rng.point(0.4),
            normal: axis * beta.cos() + across * beta.sin(),
        };
        let region = ball(p + rng.point(0.2), radius);
        (plane, cone, region)
    }

    /// The kind of a conic meeting, for counting.
    fn kind(t: &Track) -> &'static str {
        match *t {
            Track::Line { .. } | Track::Ray { .. } => "rulings",
            Track::Circle { .. } => "circle",
            Track::Ellipse { .. } => "ellipse",
            Track::Conic { c, .. } if c.abs() < 1e-6 => "parabola",
            Track::Conic { c, .. } if c > 0.0 => "half ellipse",
            Track::Conic { .. } => "hyperbola",
        }
    }

    #[test]
    fn planes_across_cones_meet_them_in_conics() {
        let mut rng = Rng(0x3EE7_0004);
        let mut kinds = std::collections::BTreeMap::new();
        let mut sampled = 0;
        for case in 0..900 {
            let pick = case % 6;
            let draw = rng.range(0.05, FRAC_PI_2);
            let jitter = [1e-12, -1e-12, 1e-8, -1e-8, 1e-4, -1e-4][rng.below(6)];
            // The cone cone_and_plane is about to draw, for its half-angle.
            let (_, cone, _) = cone_and_plane(&mut Rng(rng.0), 0.0, 2.0);
            let Quadric::Cone { slope, .. } = cone else {
                unreachable!()
            };
            let alpha = slope.atan();
            // The normal's angle from the axis: random; along a ruling's
            // direction (a parabola) and near it; square to the axis (a
            // plane along it: a hyperbola); between (hyperbolas).
            let beta = match pick {
                0 | 1 => draw,
                2 => FRAC_PI_2 - alpha,
                3 => FRAC_PI_2 - alpha + jitter,
                4 => FRAC_PI_2,
                _ => FRAC_PI_2 - alpha * draw / FRAC_PI_2,
            };
            let (plane, cone, region) = cone_and_plane(&mut rng, beta, 2.0);
            let name = format!("case {case}");
            let m = meet(&plane, &cone, &region, EPS);
            match &m {
                Meet::Curves(tracks) => {
                    sampled += check(&plane, &cone, &region, &name);
                    assert_eq!(meet(&cone, &plane, &region, EPS), m, "{name}");
                    for t in tracks {
                        *kinds.entry(kind(t)).or_insert(0) += 1;
                    }
                }
                Meet::Apart => {
                    // Nothing crosses within the region (beyond it, the
                    // curve can).
                    let within = |p: &DVec3| (*p - region.center).length() <= region.radius;
                    let crossing = |on: &Quadric, other: &Quadric| {
                        chords(&grid(on, &region, 80))
                            .into_iter()
                            .filter(|(p, q)| within(p) && within(q))
                            .filter(|(p, q)| other.distance(*p) * other.distance(*q) < 0.0)
                            .count()
                    };
                    let crossing = crossing(&plane, &cone) + crossing(&cone, &plane);
                    assert_eq!(crossing, 0, "{name}");
                    *kinds.entry("apart").or_insert(0) += 1;
                }
                other => panic!("{name}: {other:?}"),
            }
        }
        assert!(sampled > 100_000, "{sampled}");
        for k in ["ellipse", "half ellipse", "parabola", "hyperbola", "apart"] {
            assert!(kinds.get(k).copied().unwrap_or(0) >= 20, "{kinds:?}");
        }
    }

    #[test]
    fn conic_parameters_edges_and_boxes_are_exact() {
        let mut rng = Rng(0x3EE7_0005);
        let mut checked = 0;
        for case in 0..600 {
            let beta = match case % 3 {
                0 => rng.range(0.05, FRAC_PI_2),
                1 => FRAC_PI_2,
                _ => FRAC_PI_2 - rng.range(0.0, 1.4),
            };
            let (plane, cone, region) = cone_and_plane(&mut rng, beta, 2.0);
            let Meet::Curves(tracks) = meet(&plane, &cone, &region, EPS) else {
                continue;
            };
            let mag = 1.0 + plane.magnitude() + cone.magnitude();
            for t in tracks {
                let Track::Conic { y, .. } = t else { continue };
                checked += 1;
                let (lo, hi) = span(&t, &region);
                // Parameters: exact for points on the curve, and the foot
                // of a point off it.
                for i in 0..=20 {
                    let s = lo + (hi - lo) * i as f64 / 20.0;
                    let q = t.point(s);
                    if (q - region.center).length() > 2.0 * region.radius {
                        continue;
                    }
                    assert!((t.point(t.param(q)) - q).length() < 1e-13 * mag);
                    // Off the curve along the plane, square to it.
                    let off = y.cross(t.deriv(s)).cross(t.deriv(s)).normalize_or_zero();
                    let off = (off - plane_normal(&plane) * off.dot(plane_normal(&plane)))
                        .normalize_or_zero();
                    if off == DVec3::ZERO {
                        continue;
                    }
                    let (u, d) = t.nearest(q + off * 1e-6);
                    assert!((d - 1e-6).abs() < 1e-12, "case {case}: {d}");
                    assert!((t.point(u) - q).length() < 1e-9, "case {case}");
                }
                // An edge over a random stretch: its ends, and every point
                // on both surfaces.
                let ta = rng.range(lo, hi);
                let tb = rng.range(ta, hi);
                if (t.point(tb) - t.point(ta)).length() < 1e-3 {
                    continue;
                }
                let (curve, t0, t1) = t.edge(ta, tb).unwrap();
                let Curve3::Nurbs(n) = &curve else {
                    panic!("case {case}: {curve:?}")
                };
                assert!(n.weights().iter().all(|&w| w > 0.0), "case {case}");
                let reach = mag
                    + (t.point(ta) - t.point(0.0)).length()
                    + (t.point(tb) - t.point(0.0)).length();
                assert!((curve.point(t0) - t.point(ta)).length() < 1e-14 * 8.0 * reach);
                assert!((curve.point(t1) - t.point(tb)).length() < 1e-14 * 8.0 * reach);
                for i in 0..=100 {
                    let q = curve.point(t0 + (t1 - t0) * i as f64 / 100.0);
                    assert!(plane.distance(q).abs() < 1e-13 * reach, "case {case}");
                    assert!(cone.distance(q).abs() < 1e-13 * reach, "case {case}");
                }
                // The stretch's box holds it, and is reached on every side.
                let (blo, bhi) = t.bounds(ta, tb);
                let (mut slo, mut shi) = (DVec3::splat(f64::MAX), DVec3::splat(f64::MIN));
                for i in 0..=2000 {
                    let q = t.point(ta + (tb - ta) * i as f64 / 2000.0);
                    assert!(
                        q.cmpge(blo - 1e-12 * reach).all() && q.cmple(bhi + 1e-12 * reach).all()
                    );
                    slo = slo.min(q);
                    shi = shi.max(q);
                }
                let sag = (tb - ta) / 2000.0;
                let slack = 1e-12 * reach
                    + sag * sag * t.deriv(ta).length().max(t.deriv(tb).length())
                    + 1e-6 * (tb - ta);
                assert!(
                    (slo - blo).max_element() < slack && (bhi - shi).max_element() < slack,
                    "case {case}"
                );
            }
        }
        assert!(checked > 300, "{checked}");
    }

    /// A plane's unit normal.
    fn plane_normal(q: &Quadric) -> DVec3 {
        match *q {
            Quadric::Plane { normal, .. } => normal,
            _ => unreachable!(),
        }
    }

    /// The foot of a point near half an ellipse, near the ends of the
    /// half (on its minor axis), where the conic's parameter runs ever
    /// faster: the distance found is the point's distance from the curve
    /// (found by review: 360 times it, Newton's method in the parameter
    /// stopping short, 10⁻⁸ from an end).
    #[test]
    fn nearest_holds_at_the_ends_of_half_an_ellipse() {
        let mut rng = Rng(0x3EE7_0007);
        let mut cases: Vec<(f64, f64, f64, f64)> = vec![(0.158, 0.002, 10f64.powf(-11.9), -4.8e-7)];
        for _ in 0..4000 {
            let c = 10f64.powf(rng.range(-3.0, 0.0));
            let semi = 10f64.powf(rng.range(-3.0, -1.0));
            let gap = 10f64.powf(rng.range(-14.0, -1.0));
            let off = 10f64.powf(rng.range(-9.0, -6.0)) * if rng.next() < 0.5 { -1.0 } else { 1.0 };
            cases.push((c, semi, gap, off));
        }
        for (k, &(c, semi, gap, off)) in cases.iter().enumerate() {
            let track = Track::Conic {
                vertex: DVec3::ZERO,
                x: DVec3::X,
                y: DVec3::Y,
                p: semi,
                c,
            };
            let end = 1.0 / c.sqrt();
            for t in [end * (1.0 - gap), -end * (1.0 - gap)] {
                // Off the curve in its plane, square to it, by less than
                // its least radius of curvature (`p`): the foot is `t`.
                let normal = DVec3::Z.cross(track.deriv(t)).normalize();
                let q = track.point(t) + normal * off;
                let (u, d) = track.nearest(q);
                // Near the end a step of the parameter's last bits moves
                // the point `h` along the curve, which puts the nearest
                // representable foot up to `h² / 2d` farther.
                let h = 4.0 * f64::EPSILON * end * track.deriv(u).length();
                let allowed = 1e-13 * (1.0 + semi / c) + h * h / off.abs();
                assert!(
                    (d - off.abs()).abs() <= allowed,
                    "case {k}: c {c}, p {semi}, {gap} from the end: {d} for {off}"
                );
                assert!(
                    (track.point(u) - track.point(t)).length() <= 1e-6 * semi,
                    "case {k}"
                );
            }
        }
    }

    /// Conics near a parabola, in a region ten times the cone's size: the
    /// vertex form keeps every point on both surfaces to rounding wherever
    /// the configuration falls (no centre far away, no snap).
    #[test]
    fn conics_near_a_parabola_are_exact() {
        let mut rng = Rng(0x3EE7_0006);
        let mut kinds = std::collections::BTreeMap::new();
        for case in 0..300 {
            let (_, cone, _) = cone_and_plane(&mut Rng(rng.0), 0.0, 2.0);
            let Quadric::Cone { slope, .. } = cone else {
                unreachable!()
            };
            let jitter = [
                0.0, 1e-15, -1e-15, 1e-12, -1e-12, 1e-9, -1e-9, 1e-6, -1e-6, 1e-3, -1e-3,
            ][case % 11];
            let beta = FRAC_PI_2 - slope.atan() + jitter;
            let (plane, cone, region) = cone_and_plane(&mut rng, beta, 20.0);
            let m = meet(&plane, &cone, &region, EPS);
            let Meet::Curves(tracks) = &m else {
                *kinds.entry("apart").or_insert(0) += 1;
                assert_eq!(m, Meet::Apart, "case {case}");
                continue;
            };
            let mag =
                1.0 + plane.magnitude() + cone.magnitude() + region.center.abs().max_element();
            for t in tracks {
                *kinds.entry(kind(t)).or_insert(0) += 1;
                let (lo, hi) = span(t, &region);
                for i in 0..=400 {
                    let q = t.point(lo + (hi - lo) * i as f64 / 400.0);
                    if (q - region.center).length() > region.radius {
                        continue;
                    }
                    let scale = mag + (q - region.center).length();
                    let (da, db) = (plane.distance(q), cone.distance(q));
                    assert!(
                        da.abs() < 1e-13 * scale && db.abs() < 1e-13 * scale,
                        "case {case}: {da} {db} {t:?}"
                    );
                }
            }
        }
        // Every side of the parabola is reached: open conics of both
        // signs, and ellipses closed and halved.
        for (k, least) in [
            ("parabola", 100),
            ("hyperbola", 20),
            ("half ellipse", 5),
            ("ellipse", 5),
        ] {
            assert!(kinds.get(k).copied().unwrap_or(0) >= least, "{kinds:?}");
        }
    }

    /// A plane near the apex of a frustum far from it: within a sixteenth
    /// of the tolerance, two rulings within it of both surfaces; beyond, a
    /// hyperbola's branch on both to rounding. Also where the region's
    /// centre lies below the apex's height, beside a wide cone.
    #[test]
    fn planes_near_an_apex_beyond_the_region() {
        let mut rng = Rng(0x3EE7_0007);
        let (mut lines, mut below_cases, mut clipped) = (0, 0, 0);
        for case in 0..200 {
            let u = rng.unit();
            let k = 10f64.powf(rng.range(-1.0, 0.5));
            let across = u.cross(rng.unit()).normalize();
            let apex = rng.point(1.0);
            // The frustum's part where the region is: two to four units up
            // the axis from the apex.
            let h = rng.range(2.0, 4.0);
            let cone = Quadric::Cone {
                point: apex + u * h,
                axis: u,
                radius: k * h,
                slope: k,
            };
            // A plane through the apex nearer the axis's direction than the
            // rulings (its normal at β from the axis, β beyond π/2 − α): it
            // crosses the cone in two rulings, h (u − cot β v ± √(k² −
            // cot² β) across) from the apex, v = across × u.
            let beta = FRAC_PI_2 - rng.range(0.0, 0.8) * k.atan();
            let v = across.cross(u);
            let n = u * beta.cos() + v * beta.sin();
            let cot = beta.cos() / beta.sin();
            let ruling = u - v * cot + across * (k * k - cot * cot).sqrt();
            let mut regions = vec![ball(apex + ruling * h + rng.point(0.1), 0.5 * h.min(1.0))];
            // Beside a wide cone, a region whose centre lies below the
            // apex's height and which reaches the ruling at `foot`.
            let foot = apex + ruling.normalize() * h;
            let center = foot - u * ((foot - apex).dot(u) + 0.05 * h);
            let radius = (foot - center).length() + 0.05 * h;
            if (center - apex).length() > radius + 1e-3 {
                regions.push(ball(center, radius));
            }
            for (region, (off, ruled)) in regions.iter().flat_map(|r| {
                [
                    (0.0, true),
                    (EPS / 32.0, true),
                    (-EPS / 32.0, true),
                    (4.0 * EPS, false),
                    (1e-3, false),
                ]
                .map(|o| (r, o))
            }) {
                let below = (region.center - apex).dot(u) < 0.0;
                below_cases += usize::from(below);
                let plane = Quadric::Plane {
                    point: apex + n * off,
                    normal: n,
                };
                let m = meet(&plane, &cone, region, EPS);
                assert_eq!(meet(&cone, &plane, region, EPS), m);
                let Meet::Curves(tracks) = &m else {
                    // Rulings that miss the region, which the one below
                    // the apex reaches.
                    assert!(!(below && ruled), "case {case} {off}: {m:?}");
                    assert_eq!(
                        m,
                        Meet::Apart,
                        "case {case} {off} {plane:?} {cone:?} {region:?}"
                    );
                    continue;
                };
                if below && ruled {
                    assert!(
                        tracks.iter().any(|t| off_track(t, foot) <= EPS / 16.0),
                        "case {case} {off}: {tracks:?}"
                    );
                }
                for t in tracks {
                    assert_eq!(t.straight(), ruled, "case {case} {off}: {t:?}");
                    lines += usize::from(t.straight());
                    if let Track::Ray {
                        origin,
                        dir,
                        start,
                        from_apex,
                    } = *t
                    {
                        // Along a ruling, away from the apex.
                        let w = apex - origin;
                        assert!(
                            (w - dir * w.dot(dir)).length() <= EPS / 16.0
                                && w.dot(dir) <= start + EPS / 16.0,
                            "case {case} {off}: {t:?}"
                        );
                        let c = region.center - origin;
                        if start == f64::INFINITY {
                            // Farther than CLIP tolerances from the region.
                            assert!(
                                (c - dir * c.dot(dir)).length() > region.radius + CLIP * EPS,
                                "case {case} {off}: {t:?}"
                            );
                            continue;
                        }
                        if from_apex {
                            assert!(
                                (t.point(start) - apex).length() <= EPS / 16.0,
                                "case {case} {off}: {t:?}"
                            );
                        } else {
                            // Where it comes within CLIP tolerances of the
                            // region, nearer its centre from there on.
                            assert!(
                                (t.point(start) - region.center).length()
                                    >= region.radius + (CLIP - 1e-3) * EPS
                                    && start <= c.dot(dir) + 1e-12,
                                "case {case} {off}: {t:?}"
                            );
                            clipped += 1;
                        }
                    }
                    let (lo, hi) = span(t, region);
                    for i in 0..=100 {
                        let q = t.point(lo + (hi - lo) * i as f64 / 100.0);
                        if (q - region.center).length() > region.radius {
                            continue;
                        }
                        let limit = if ruled { EPS / 16.0 } else { 1e-12 * (1.0 + h) };
                        assert!(plane.distance(q).abs() <= limit, "case {case} {off}");
                        assert!(cone.distance(q).abs() <= limit, "case {case} {off}");
                    }
                }
            }
        }
        assert!(lines > 400, "{lines}");
        assert!(below_cases > 100, "{below_cases}");
        assert!(clipped > 400, "{clipped}");
    }
}
