//! Where two surfaces of step 2 meet, in closed form (design note 0012
//! §10.3): lines, circles and ellipses.
//!
//! Two surfaces of revolution about one axis (a plane across the axis, a
//! cylinder, a cone, a sphere with its centre on the axis) meet in circles
//! about that axis, where their meridians (their sections by a half-plane
//! bounded by the axis: lines and circles) meet. Otherwise two planes meet
//! in a line, a plane and a cylinder in two lines along the axis or an
//! ellipse, and two cylinders with parallel axes in two lines. Other pairs
//! meet in curves of degree four, or (a plane obliquely across a cone) in
//! conics that come with cones (§10.7, item 3); they are refused as
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

use std::f64::consts::FRAC_PI_2;

use glam::{DVec2, DVec3};
use opendatum_nurbs::NurbsCurve;

use super::quadric::Quadric;
use crate::geom::Curve3;
use crate::KernelError;

/// The part of space where two faces are: a ball holding both faces'
/// common bounding box.
#[derive(Clone, Copy, Debug, PartialEq)]
pub(super) struct Region {
    pub center: DVec3,
    pub radius: f64,
}

/// A curve where two surfaces meet.
#[derive(Clone, Copy, Debug, PartialEq)]
pub(super) enum Track {
    /// `origin + t · dir`, `dir` a unit vector.
    Line { origin: DVec3, dir: DVec3 },
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
    /// Their curves are not handled in closed form yet.
    Unsupported(&'static str),
}

impl Track {
    /// The point at `t`.
    pub fn point(&self, t: f64) -> DVec3 {
        match *self {
            Track::Line { origin, dir } => origin + dir * t,
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
        }
    }

    /// The first derivative at `t`.
    pub fn deriv(&self, t: f64) -> DVec3 {
        match *self {
            Track::Line { dir, .. } => dir,
            Track::Circle { x, y, radius, .. } => {
                let (s, c) = t.sin_cos();
                (y * c - x * s) * radius
            }
            Track::Ellipse { x, y, a, b, .. } => {
                let (s, c) = t.sin_cos();
                y * (b * c) - x * (a * s)
            }
        }
    }

    /// The second derivative at `t`.
    pub fn second(&self, t: f64) -> DVec3 {
        match *self {
            Track::Line { .. } => DVec3::ZERO,
            Track::Circle { center, .. } | Track::Ellipse { center, .. } => center - self.point(t),
        }
    }

    /// The parameter of a point on the curve (within rounding of it); for
    /// a closed curve, in `[−π/2, 3π/2)` like the kernel's circles.
    pub fn param(&self, p: DVec3) -> f64 {
        match *self {
            Track::Line { origin, dir } => (p - origin).dot(dir),
            Track::Circle { center, x, y, .. } => {
                let d = p - center;
                wrap(d.dot(y).atan2(d.dot(x)))
            }
            Track::Ellipse { center, x, y, a, b } => {
                let d = p - center;
                wrap((d.dot(y) / b).atan2(d.dot(x) / a))
            }
        }
    }

    /// True for a closed curve (a circle or an ellipse), of period 2π.
    pub fn closed(&self) -> bool {
        !matches!(self, Track::Line { .. })
    }

    /// The edge geometry of the curve from `ta` to `tb` (`ta < tb`, less
    /// than a turn apart on a closed curve), with its parameter range:
    /// a line or circle as itself, an elliptical arc as an exact rational
    /// quadratic B-spline in pieces of at most a quarter turn (Piegl &
    /// Tiller A7.1 under the affine map from the unit circle, as the
    /// sketch's ellipses).
    pub fn edge(&self, ta: f64, tb: f64) -> Result<(Curve3, f64, f64), KernelError> {
        Ok(match *self {
            Track::Line { origin, dir } => (Curve3::Line { origin, dir }, ta, tb),
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
        })
    }
}

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
                return Meet::Unsupported("two cylinders at an angle meet in a quartic curve");
            }
            parallel_cylinders(a1, u1, r1, a2, u2, r2, region, eps)
        }
        (Quadric::Plane { .. }, Quadric::Cone { .. }) => {
            Meet::Unsupported("a plane across a cone off its axis meets it in a conic")
        }
        _ => Meet::Unsupported("these surfaces meet in a quartic curve"),
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
            let reach = match other {
                Quadric::Cone {
                    point: c,
                    axis: v,
                    radius,
                    slope,
                } => {
                    // The cone's radius where the plane crosses the axis.
                    let h = (point - c).dot(normal) / normal.dot(v);
                    (radius + slope * h).abs()
                }
                Quadric::Cylinder { radius, .. } => radius,
                _ => return None,
            };
            (sine * reach <= snap).then_some(axis)
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

    /// The distance from `p` to the curve, for points near it.
    fn off_track(t: &Track, p: DVec3) -> f64 {
        match *t {
            Track::Line { origin, dir } => {
                let w = p - origin;
                (w - dir * w.dot(dir)).length()
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
        for t in tracks {
            let (lo, hi) = if t.closed() {
                (0.0, TAU)
            } else {
                let c = t.param(region.center);
                (c - region.radius, c + region.radius)
            };
            for i in 0..=64 {
                let p = t.point(lo + (hi - lo) * i as f64 / 64.0);
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
            let region = Region {
                center: p + rng.point(0.2),
                radius: 2.0,
            };
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
        let region = Region {
            center: DVec3::new(0.1, 0.2, 0.3),
            radius: 1.5,
        };
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
                Meet::Unsupported("two cylinders at an angle meet in a quartic curve"),
            ),
            (
                "plane across a cone",
                plane(o, (u + x).normalize()),
                cone,
                Meet::Unsupported("a plane across a cone off its axis meets it in a conic"),
            ),
            (
                "sphere off a cylinder's axis",
                cyl(0.3, o),
                sphere(o + x * 0.2, 0.3),
                Meet::Unsupported("these surfaces meet in a quartic curve"),
            ),
            (
                "cone beside a cylinder",
                cyl(0.3, o + x * 0.5),
                cone,
                Meet::Unsupported("these surfaces meet in a quartic curve"),
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
        let region = Region {
            center: DVec3::ZERO,
            radius: 1.0,
        };
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
            let region = Region {
                center: DVec3::ZERO,
                radius: 1.0,
            };
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
}
