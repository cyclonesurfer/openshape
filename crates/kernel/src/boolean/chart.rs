//! Charts of faces (design note 0012 §10.4 and §12.1): plane coordinates
//! on a face's surface, the traces of curves in them with bounds on how
//! far each strays from its chords, and the outlines of faces, for
//! deciding whether a point of the surface lies in a face.
//!
//! A plane's chart is its frame; a cylinder's is the cylinder unrolled
//! about the middle of the face (which spans at most half a turn, ADR
//! 0008), so the face never reaches the cut of the unrolling; a cone's is
//! the cone developed about the face's middle ruling, anchored at the
//! face. These three are isometries of the surface onto the plane. A
//! sphere has none: its chart is the stereographic projection from the
//! antipode of the middle of the half of the sphere that holds the face,
//! which is conformal and takes circles to circles, with a scale between 1
//! and 2 over that half. Every chart is
//! oriented so that a face's outer loop, counter-clockwise about its
//! outward normal, is counter-clockwise in the chart.
//!
//! Lengths in a chart are compared with lengths in space by two bounds,
//! [`Chart::reach`] and [`Chart::space`]; no chart makes a chord in space
//! longer than its image, so a band measured in a chart bounds the
//! distance in space too.
//!
//! A curve's course in a chart is followed by chords. The distance from a
//! course to the chord between two of its points at parameters `t0` and
//! `t1` is at most `M (t1 − t0)² / 8`, `M` a bound on the second derivative
//! of the course (the error of linear interpolation); the second derivative
//! of a curve's image in an isometric chart of a developable surface is the
//! part of its acceleration along the surface, so `M` is at most the size
//! of the curve's acceleration, and less for the curves that unroll to
//! lines (rulings, and circles square to a cylinder's axis, which unroll
//! straight) and for circles about a cone's axis. A circle's image in a
//! sphere's chart is a circle, whose chords are bounded by their sagittas.
//! A NURBS curve in a plane is bounded instead by its Bézier pieces'
//! control points, whose convex hull holds the piece.

use std::f64::consts::{FRAC_PI_2, FRAC_PI_8, PI};

use glam::{DVec2, DVec3};

use super::meet::{self, Track};
use super::planar::{point_segment, Bend};
use super::quadric::{self, Quadric};
use super::shape::{self, Face, Shape};
use crate::bvh::Bvh;
use crate::geom::Curve3;
use crate::math::Frame;
use crate::KernelError;

/// A chart of a face's surface.
#[derive(Clone, Copy, Debug, PartialEq)]
pub(super) enum Chart {
    /// A plane's frame: `(x, y)`, `z` the face's outward normal.
    Plane(Frame),
    /// A cylinder unrolled: `(radius · θ, h)`, `θ ∈ (−π, π]` the angle about
    /// the unit `axis` from `x` toward `y` (both square to the axis), `h`
    /// the position along the axis from `point`.
    Cylinder {
        point: DVec3,
        axis: DVec3,
        x: DVec3,
        y: DVec3,
        radius: f64,
    },
    /// A cone developed (design note 0012 §12.1): a point at slant distance
    /// `s` from the apex and angle θ about the unit `axis` (from `x` toward
    /// `y`, both square to the axis) goes to polar coordinates `(s, θ sin
    /// α)` about the apex's image, shifted so that the anchor, the point
    /// at `radius` from `point` toward `x` (slant distance `s0 = radius /
    /// sin α` from the apex), is the origin: `(σ cos ψ − 2 s0 sin²(ψ/2),
    /// (s0 + σ) sin ψ)` with `ψ = θ sin α`, `σ = s − s0`.
    Cone {
        point: DVec3,
        axis: DVec3,
        x: DVec3,
        y: DVec3,
        radius: f64,
        sin: f64,
        cos: f64,
    },
    /// A sphere projected stereographically from the antipode of the unit
    /// direction `c` onto the plane tangent at `c`: the point at angle φ
    /// from `c` goes to `2 radius tan(φ/2)` from the origin, toward its
    /// direction in the plane, in coordinates along `x` and `y` (square to
    /// `c` and each other).
    Sphere {
        center: DVec3,
        radius: f64,
        c: DVec3,
        x: DVec3,
        y: DVec3,
    },
}

impl Chart {
    /// The chart point of `p` (a point on or near the surface).
    pub fn map(&self, p: DVec3) -> DVec2 {
        match *self {
            Chart::Plane(f) => {
                let d = p - f.origin;
                DVec2::new(d.dot(f.x), d.dot(f.y))
            }
            Chart::Cylinder {
                point,
                axis,
                x,
                y,
                radius,
            } => {
                let w = p - point;
                DVec2::new(radius * w.dot(y).atan2(w.dot(x)), w.dot(axis))
            }
            Chart::Cone {
                point,
                axis,
                x,
                y,
                radius,
                sin,
                cos,
            } => {
                let w = p - point;
                let h = w.dot(axis);
                let across = w - axis * h;
                // The slant position from the anchor's parallel, from small
                // numbers (the foot on the ruling for a point off the cone).
                let sigma = (across.length() - radius) * sin + h * cos;
                let psi = across.dot(y).atan2(across.dot(x)) * sin;
                let s0 = radius / sin;
                let half = (0.5 * psi).sin();
                DVec2::new(
                    sigma * psi.cos() - 2.0 * s0 * half * half,
                    (s0 + sigma) * psi.sin(),
                )
            }
            Chart::Sphere {
                center,
                radius,
                c,
                x,
                y,
            } => {
                let u = (p - center).normalize();
                let k = 2.0 * radius / (1.0 + u.dot(c));
                DVec2::new(u.dot(x) * k, u.dot(y) * k)
            }
        }
    }

    /// The point of the surface at chart point `q`.
    pub fn unmap(&self, q: DVec2) -> DVec3 {
        match *self {
            Chart::Plane(f) => f.point2(q),
            Chart::Cylinder {
                point,
                axis,
                x,
                y,
                radius,
            } => {
                let (s, c) = (q.x / radius).sin_cos();
                point + (x * c + y * s) * radius + axis * q.y
            }
            Chart::Cone {
                point,
                axis,
                x,
                y,
                radius,
                sin,
                cos,
            } => {
                let s0 = radius / sin;
                let (px, py) = (q.x + s0, q.y);
                let s = px.hypot(py);
                // σ = s − s0 without cancelling when the apex is far.
                let sigma = if s + s0 > 0.0 {
                    (q.x * (q.x + 2.0 * s0) + py * py) / (s + s0)
                } else {
                    0.0
                };
                let theta = py.atan2(px) / sin;
                let (st, ct) = theta.sin_cos();
                point + axis * (sigma * cos) + (x * ct + y * st) * (radius + sigma * sin)
            }
            Chart::Sphere {
                center,
                radius,
                c,
                x,
                y,
            } => {
                let w = q / (2.0 * radius);
                let w2 = w.length_squared();
                center
                    + (x * (2.0 * w.x) + y * (2.0 * w.y) + c * (1.0 - w2)) * (radius / (1.0 + w2))
            }
        }
    }

    /// The chart of the face turned over (its outward normal reversed), in
    /// which its loops run counter-clockwise again: a plane's frame with
    /// `y` and `z` reversed, a cylinder unrolled the other way.
    pub fn reversed(&self) -> Chart {
        match *self {
            Chart::Plane(f) => Chart::Plane(Frame {
                origin: f.origin,
                x: f.x,
                y: -f.y,
                z: -f.z,
            }),
            Chart::Cylinder {
                point,
                axis,
                x,
                y,
                radius,
            } => Chart::Cylinder {
                point,
                axis,
                x,
                y: -y,
                radius,
            },
            Chart::Cone {
                point,
                axis,
                x,
                y,
                radius,
                sin,
                cos,
            } => Chart::Cone {
                point,
                axis,
                x,
                y: -y,
                radius,
                sin,
                cos,
            },
            Chart::Sphere {
                center,
                radius,
                c,
                x,
                y,
            } => Chart::Sphere {
                center,
                radius,
                c,
                x,
                y: -y,
            },
        }
    }

    /// A distance in space within which no point of the face `d` or more
    /// from chart point `q` in the chart lies (both points within half a
    /// turn of the chart's middle): `d` itself on a plane; on a cylinder the
    /// chord of the arc `d` long; on a cone `d` times the least ratio of a
    /// chord to its image, `sin α / sin(π sin α / 2)`; on a sphere from the
    /// chordal distance, `d / √(λ(q) λ)` with `λ` the greatest scale within
    /// `d` of `q`.
    pub fn space(&self, q: DVec2, d: f64) -> f64 {
        match *self {
            Chart::Plane(_) => d,
            Chart::Cylinder { radius, .. } => {
                if d >= PI * radius {
                    2.0 * radius
                } else {
                    2.0 * radius * (d / (2.0 * radius)).sin()
                }
            }
            Chart::Cone { sin, .. } => d * sin / (FRAC_PI_2 * sin).sin(),
            Chart::Sphere { radius, .. } => {
                // The chord between chart points q and x is |q − x| /
                // √(λ(q) λ(x)), which grows with |q − x| when λ(x) is bounded
                // by its value |q| + |q − x| from the origin.
                let k = 1.0 / (4.0 * radius * radius);
                let far = q.length() + d;
                d / ((1.0 + q.length_squared() * k) * (1.0 + far * far * k)).sqrt()
            }
        }
    }

    /// A distance in the chart from chart point `q` beyond which every point
    /// of the face is farther than `d` from `q`'s point in space: `d` on a
    /// plane; on a cylinder the arc of the chord `d`; on a cone `d` over the
    /// ratio of a chord to its image for the angles about the axis that
    /// points that near can differ by; on a sphere from the chordal
    /// distance and the scale at the angles from the chart's middle that
    /// points that near can have.
    pub fn reach(&self, q: DVec2, d: f64) -> f64 {
        if d <= 0.0 {
            return d;
        }
        match *self {
            Chart::Plane(_) => d,
            Chart::Cylinder { radius, .. } => {
                let z = d / (2.0 * radius);
                if z >= 1.0 {
                    FRAC_PI_2 * d
                } else {
                    2.0 * radius * z.asin()
                }
            }
            Chart::Cone { radius, sin, .. } => {
                // Every chord is at least `least` times its image, so a
                // point within `d` in space is within `far` in the chart, and
                // seen from the apex it is within asin(far / s) of `q`.
                let least = sin / (FRAC_PI_2 * sin).sin();
                let far = d / least;
                let s = (q.x + radius / sin).hypot(q.y);
                if s <= 2.0 * far {
                    return far;
                }
                // Half the angle about the axis between the two points.
                let half = (far / s).asin() / (2.0 * sin);
                if half >= FRAC_PI_2 {
                    return far;
                }
                let ratio = sin * half.sin() / (sin * half).sin();
                (d / ratio).min(far)
            }
            Chart::Sphere { radius, .. } => {
                // The angle from the chart's middle of `q`'s point, and the
                // most that of a point within `d` of it can be.
                let phi = 2.0 * (q.length() / (2.0 * radius)).atan();
                let reach = phi + 2.0 * (d / (2.0 * radius)).min(1.0).asin();
                if reach >= PI {
                    return f64::INFINITY;
                }
                let scale = |a: f64| 2.0 / (1.0 + a.cos());
                d * (scale(phi) * scale(reach)).sqrt()
            }
        }
    }
}

/// The chart of face `f` of `shape`: its plane's frame, its cylinder
/// unrolled or its cone developed about the direction of the middle of its
/// boundary, or its sphere projected from the antipode of the middle of
/// the half of it that holds the face ([`Face::middle`]).
pub(super) fn chart_of(shape: &Shape, f: &Face) -> Result<Chart, KernelError> {
    // The mean of the directions from an axis to points along the
    // boundary: the face spans at most half a turn about it, and points
    // inside its edges weigh its middle.
    let middle = |point: DVec3, axis: DVec3| -> Result<DVec3, KernelError> {
        let mut sum = DVec3::ZERO;
        for l in &f.loops {
            for &(e, _) in l {
                let (curve, t0, t1) = shape.curve(e);
                for k in 0..=4 {
                    let w = curve.point(t0 + (t1 - t0) * k as f64 / 4.0) - point;
                    sum += (w - axis * w.dot(axis)).normalize_or_zero();
                }
            }
        }
        let x = sum - axis * sum.dot(axis);
        if x.length() <= 1e-6 * sum.length().max(1.0) {
            return Err(KernelError::Unsupported(
                "booleans with a face going around its axis (design note 0012 §10)".into(),
            ));
        }
        Ok(x.normalize())
    };
    match shape.surfaces[f.surface] {
        Quadric::Plane { .. } => Ok(Chart::Plane(f.frame)),
        Quadric::Cylinder {
            point,
            axis,
            radius,
        } => {
            let x = middle(point, axis)?;
            let y = axis.cross(x);
            let y = if f.outward { y } else { -y };
            Ok(Chart::Cylinder {
                point,
                axis,
                x,
                y,
                radius,
            })
        }
        Quadric::Cone {
            point,
            axis,
            radius,
            slope,
        } => {
            let x = middle(point, axis)?;
            // Radially out along the ruling, then about the axis: about the
            // outward normal (out from the axis and back toward the apex)
            // that turns from `x` toward `x × axis`.
            let y = x.cross(axis);
            let y = if f.outward { y } else { -y };
            let cos = 1.0 / (1.0 + slope * slope).sqrt();
            Ok(Chart::Cone {
                point,
                axis,
                x,
                y,
                radius,
                sin: slope * cos,
                cos,
            })
        }
        Quadric::Sphere { center, radius } => {
            let c = f.middle;
            let x = c.any_orthonormal_vector();
            let y = c.cross(x);
            let y = if f.outward { y } else { -y };
            Ok(Chart::Sphere {
                center,
                radius,
                c,
                x,
                y,
            })
        }
    }
}

/// The most first chords of one course (see [`Trace::division`]).
const MOST_CHORDS: f64 = 4096.0;

/// A bound on how far a track's image in a chart strays from its chords.
#[derive(Clone, Copy)]
enum Bound {
    /// The image's second derivative per unit parameter squared is at
    /// most this everywhere (zero for a straight image).
    Steady(f64),
    /// The image is an arc of a circle of this curvature, followed so that
    /// a stretch of a first chord turns less than a quarter turn: a chord's
    /// band is the sagitta of its arc.
    Arc(f64),
    /// A conic track (design note 0012 §14), in a chart that keeps lengths:
    /// the image's second derivative is at most `p / w³`, `w = √(1 − c t²)`
    /// least over each stretch.
    Conic { p: f64, c: f64 },
}

/// How a course is followed.
#[derive(Clone, Copy)]
enum Path<'a> {
    /// A track (a circle or an ellipse, by angle, or a conic), with a
    /// bound on its image.
    Track(Track, Bound),
    /// A NURBS curve in a plane's chart.
    Nurbs(&'a Curve3),
}

/// A curved segment's course in a chart.
pub(super) struct Trace<'a> {
    chart: Chart,
    path: Path<'a>,
    /// A bound on the rounding of chart points, added to every band.
    slack: f64,
}

impl Trace<'_> {
    fn point(&self, t: f64) -> DVec3 {
        match self.path {
            Path::Track(track, _) => track.point(t),
            Path::Nurbs(curve) => curve.point(t),
        }
    }

    /// Parameters dividing `[t0, t1]` into first chords: at most an
    /// eighth of a half turn each along a track (on a sphere, few enough
    /// that each chord's image turns at most a quarter turn; on a conic, of
    /// the turning of its tangent), two per Bézier piece of a NURBS curve.
    pub fn division(&self, t0: f64, t1: f64) -> Result<Vec<f64>, KernelError> {
        let mut out = vec![t0];
        match self.path {
            Path::Track(_, Bound::Conic { c, .. }) => {
                // A conic turns less than a half turn: at most eight.
                let (f0, f1) = (meet::conic_turn(c, t0), meet::conic_turn(c, t1));
                let n = ((f1 - f0).abs() / FRAC_PI_8).ceil().clamp(1.0, 8.0) as usize;
                for k in 1..n {
                    out.push(meet::conic_at_turn(c, f0 + (f1 - f0) * k as f64 / n as f64));
                }
            }
            Path::Track(track, bound) => {
                let span = (t1 - t0).abs();
                let n = match bound {
                    Bound::Steady(0.0) => 1.0,
                    Bound::Steady(_) | Bound::Conic { .. } => (span / FRAC_PI_8).ceil(),
                    Bound::Arc(kappa) => (span / FRAC_PI_8)
                        .ceil()
                        .max((kappa * self.chart_length(&track, t0, t1) / FRAC_PI_2).ceil()),
                };
                // A course in a face turns at most a few times (a sphere's
                // face lies within a half of it, where the chart's scale is
                // at most 2): more chords than this is a course that runs
                // toward the chart's cut, not one in the face.
                if n.is_nan() || n > MOST_CHORDS {
                    return Err(KernelError::InvalidResult(
                        "boolean: a curve was traced toward the cut of a face's chart".into(),
                    ));
                }
                let n = n.max(1.0) as usize;
                for k in 1..n {
                    out.push(t0 + (t1 - t0) * k as f64 / n as f64);
                }
            }
            Path::Nurbs(curve) => {
                let (lo, hi) = if t0 <= t1 { (t0, t1) } else { (t1, t0) };
                let mut ks: Vec<f64> = Vec::new();
                if let Curve3::Nurbs(n) = curve {
                    for seg in n.bezier_segments() {
                        for t in [seg.a, 0.5 * (seg.a + seg.b)] {
                            if t > lo && t < hi {
                                ks.push(t);
                            }
                        }
                    }
                }
                ks.push(0.5 * (lo + hi));
                ks.sort_by(f64::total_cmp);
                ks.dedup();
                if t0 > t1 {
                    ks.reverse();
                }
                out.extend(ks);
            }
        }
        out.push(t1);
        Ok(out)
    }

    /// A bound on the length in a sphere's chart of the circle `track`
    /// between `t0` and `t1`: its length in space times the greatest scale
    /// of the chart along it, `2 / (1 + cos φ)` at its greatest angle φ from
    /// the chart's middle, found in closed form (the cosine of φ along the
    /// circle is a sinusoid in its angle).
    fn chart_length(&self, track: &Track, t0: f64, t1: f64) -> f64 {
        let (
            Chart::Sphere {
                center,
                radius: big,
                c,
                ..
            },
            Track::Circle {
                center: o,
                x,
                y,
                radius,
            },
        ) = (self.chart, *track)
        else {
            return f64::INFINITY;
        };
        // cos φ(t) = ((o − center) · c + radius (cos t x + sin t y) · c) / big.
        let least = shape::least_of_wave(x.dot(c) * radius, y.dot(c) * radius, t0, t1);
        let cos = ((o - center).dot(c) + least) / big;
        if cos <= -1.0 + 1e-12 {
            return f64::INFINITY;
        }
        radius * (t1 - t0).abs() * 2.0 / (1.0 + cos)
    }
}

impl Bend for Trace<'_> {
    fn at(&self, t: f64) -> DVec2 {
        self.chart.map(self.point(t))
    }

    fn band(&self, t0: f64, t1: f64) -> f64 {
        let b = match self.path {
            Path::Track(_, Bound::Steady(m)) => m * (t1 - t0) * (t1 - t0) / 8.0,
            Path::Track(track, Bound::Conic { p, c }) => {
                // w is least at the greatest |t| on half an ellipse, else at
                // the |t| nearest zero; the rounding of chart points grows
                // with the stretch's reach from the vertex.
                let (lo, hi) = (t0.min(t1), t0.max(t1));
                let t = if c > 0.0 {
                    lo.abs().max(hi.abs())
                } else if lo <= 0.0 && hi >= 0.0 {
                    0.0
                } else {
                    lo.abs().min(hi.abs())
                };
                let w = meet::conic_w(c, t);
                let reach = (track.point(t0) - track.point(0.0))
                    .length()
                    .max((track.point(t1) - track.point(0.0)).length());
                p / (w * w * w) * (t1 - t0) * (t1 - t0) / 8.0 + 1e-13 * reach
            }
            Path::Track(_, Bound::Arc(kappa)) => {
                // The sagitta of an arc of curvature κ less than a half turn
                // over a chord of length 2h: h² κ / (1 + √(1 − (κ h)²)).
                let h = 0.5 * (self.at(t1) - self.at(t0)).length();
                let x = (kappa * h).min(1.0);
                h * x / (1.0 + (1.0 - x * x).sqrt())
            }
            Path::Nurbs(curve) => {
                let (p, q) = (self.at(t0), self.at(t1));
                quadric::hull_points(curve, t0, t1).map_or(f64::INFINITY, |pieces| {
                    pieces
                        .iter()
                        .flatten()
                        .map(|&c| point_segment(self.chart.map(c), p, q))
                        .fold(0.0, f64::max)
                })
            }
        };
        b + self.slack
    }
}

/// The course of a circle, an ellipse or a NURBS curve in `chart`, or
/// `None` for a course that is straight there (a line; a ruling or a
/// circle square to the axis on a cylinder; a ruling on a cone).
pub(super) fn trace<'a>(
    chart: &Chart,
    path: TraceOf<'a>,
) -> Result<Option<Trace<'a>>, KernelError> {
    let track =
        match path {
            TraceOf::Track(t) => t,
            TraceOf::Curve(c) => match c {
                Curve3::Line { .. } => return Ok(None),
                Curve3::Circle {
                    center,
                    x,
                    y,
                    radius,
                } => Track::Circle {
                    center: *center,
                    x: *x,
                    y: *y,
                    radius: *radius,
                },
                Curve3::Nurbs(_) => return match chart {
                    Chart::Plane(_) => Ok(Some(Trace {
                        chart: *chart,
                        path: Path::Nurbs(c),
                        slack: 1e-15 * (1.0 + c.magnitude()),
                    })),
                    _ => Err(KernelError::Unsupported(
                        "booleans with a free-form edge on a curved face (design note 0012 §10)"
                            .into(),
                    )),
                },
            },
        };
    let off_surface = || {
        KernelError::InvalidResult(
            "boolean: a curve that cannot lie on a sphere was traced on one".into(),
        )
    };
    let (bound, size) = match (chart, track) {
        (_, Track::Line { .. } | Track::Ray { .. }) => match chart {
            Chart::Sphere { .. } => return Err(off_surface()),
            _ => return Ok(None),
        },
        (Chart::Plane(_), Track::Circle { center, radius, .. }) => {
            (Bound::Steady(radius), center.abs().max_element() + radius)
        }
        (Chart::Plane(_), Track::Ellipse { center, a, b, .. }) => (
            Bound::Steady(a.max(b)),
            center.abs().max_element() + a.max(b),
        ),
        (
            Chart::Cylinder { axis, .. },
            Track::Circle {
                center,
                x,
                y,
                radius,
            },
        ) => {
            let size = center.abs().max_element() + radius;
            if x.cross(y).cross(*axis).length() <= 1e-9 {
                return Ok(None);
            }
            (Bound::Steady(radius), size)
        }
        (
            Chart::Cylinder { axis, .. },
            Track::Ellipse {
                center, x, y, a, b, ..
            },
        ) => {
            // On the cylinder the ellipse turns once about the axis at a
            // steady rate, so its image is `(r θ, h)` with θ linear in the
            // parameter: only its height bends, by at most its amplitude.
            let amp = ((a * x.dot(*axis)).powi(2) + (b * y.dot(*axis)).powi(2)).sqrt();
            (Bound::Steady(amp), center.abs().max_element() + a.max(b))
        }
        (
            Chart::Cone { axis, sin, .. },
            Track::Circle {
                center,
                x,
                y,
                radius,
            },
        ) => {
            // A circle about the axis develops into a circle about the
            // apex's image of radius `r / sin α`, followed at speed `r`: its
            // image's acceleration is `r sin α`. Any other circle is bounded
            // by its acceleration in space.
            let about = x.cross(y).cross(*axis).length() <= 1e-9;
            let m = if about { radius * sin } else { radius };
            (Bound::Steady(m), center.abs().max_element() + radius)
        }
        (Chart::Cone { .. }, Track::Ellipse { center, a, b, .. }) => (
            Bound::Steady(a.max(b)),
            center.abs().max_element() + a.max(b),
        ),
        (
            Chart::Sphere {
                center: o,
                radius: big,
                c,
                ..
            },
            Track::Circle {
                center,
                x,
                y,
                radius,
            },
        ) => {
            // The image of a circle of the sphere is a circle of curvature
            // |n · (C − S)| / (2 R r), `S` the projection's point.
            let s = *o - *c * *big;
            let kappa = x.cross(y).dot(center - s).abs() / (2.0 * big * radius);
            (Bound::Arc(kappa), o.abs().max_element() + 4.0 * big)
        }
        (Chart::Sphere { .. }, Track::Ellipse { .. }) => return Err(off_surface()),
        (Chart::Plane(_) | Chart::Cone { .. }, Track::Conic { vertex, p, c, .. }) => {
            (Bound::Conic { p, c }, vertex.abs().max_element())
        }
        (Chart::Cylinder { .. } | Chart::Sphere { .. }, Track::Conic { .. }) => {
            return Err(KernelError::InvalidResult(
                "boolean: a conic was traced on a surface it cannot lie on".into(),
            ))
        }
    };
    Ok(Some(Trace {
        chart: *chart,
        path: Path::Track(track, bound),
        slack: 1e-13 * (1.0 + size),
    }))
}

/// What a course follows.
#[derive(Clone, Copy)]
pub(super) enum TraceOf<'a> {
    Track(Track),
    Curve(&'a Curve3),
}

/// A chord of an outline.
#[derive(Clone, Copy, Debug)]
struct OChord {
    p: DVec2,
    q: DVec2,
    band: f64,
    /// The part of `band` from the course's bend, which halving narrows;
    /// the rest is how far the chord's ends are from the course's.
    sag: f64,
    /// For a chord of a curved edge: its course and its parameters at `p`
    /// and `q`.
    course: Option<(usize, f64, f64)>,
}

/// A face's boundary in its chart, for point-in-face tests.
pub(super) struct Outline<'a> {
    pub chart: Chart,
    chords: Vec<OChord>,
    traces: Vec<Trace<'a>>,
    tree: Bvh,
    /// The widest band of a chord.
    reach: f64,
}

impl<'a> Outline<'a> {
    pub fn new(shape: &'a Shape, f: &Face) -> Result<Outline<'a>, KernelError> {
        let chart = chart_of(shape, f)?;
        let mut chords = Vec::new();
        let mut traces: Vec<Trace<'a>> = Vec::new();
        for l in &f.loops {
            for &(e, _) in l {
                let (curve, t0, t1) = shape.curve(e);
                let (va, vb) = shape.edges[e];
                let (pa, pb) = (chart.map(shape.points[va]), chart.map(shape.points[vb]));
                match trace(&chart, TraceOf::Curve(curve))? {
                    None => chords.push(OChord {
                        p: pa,
                        q: pb,
                        band: 0.0,
                        sag: 0.0,
                        course: None,
                    }),
                    Some(c) => {
                        let ts = c.division(t0, t1)?;
                        let n = ts.len() - 1;
                        let at: Vec<DVec2> = (0..=n)
                            .map(|k| {
                                if k == 0 {
                                    pa
                                } else if k == n {
                                    pb
                                } else {
                                    c.at(ts[k])
                                }
                            })
                            .collect();
                        for k in 0..n {
                            let moved = (c.at(ts[k]) - at[k])
                                .length()
                                .max((c.at(ts[k + 1]) - at[k + 1]).length());
                            let sag = c.band(ts[k], ts[k + 1]);
                            chords.push(OChord {
                                p: at[k],
                                q: at[k + 1],
                                band: sag + moved,
                                sag,
                                course: Some((traces.len(), ts[k], ts[k + 1])),
                            });
                        }
                        traces.push(c);
                    }
                }
            }
        }
        let reach = chords.iter().fold(0.0f64, |m, c| m.max(c.band));
        let tree = Bvh::new(
            chords
                .iter()
                .map(|c| {
                    let r = DVec2::splat(c.band);
                    (
                        (c.p.min(c.q) - r).extend(0.0),
                        (c.p.max(c.q) + r).extend(0.0),
                    )
                })
                .collect(),
        );
        Ok(Outline {
            chart,
            chords,
            traces,
            tree,
            reach,
        })
    }

    /// Whether the chart point `q` is in the face: `Some(true)` inside,
    /// `Some(false)` outside, `None` within `margin` of its boundary in
    /// space (within [`Chart::reach`] of it in the chart). Chords near `q`
    /// are halved, for this question only, until their sags are within
    /// `fine` or their bands clear of `q`; the answer is the parity of a ray from `q`
    /// across the chords, which is the parity across the boundary for a
    /// point outside every band.
    pub fn locate(&self, q: DVec2, margin: f64, fine: f64) -> Option<bool> {
        let margin = self.chart.reach(q, margin);
        let r = self.reach + margin;
        let lo = DVec3::new(q.x - r, q.y - r, 0.0);
        let hi = DVec3::new(f64::INFINITY, q.y + r, 0.0);
        let mut near = Vec::new();
        self.tree.overlapping(lo, hi, &mut near);
        let mut inside = false;
        let mut stack: Vec<OChord> = Vec::new();
        for i in near {
            stack.push(self.chords[i]);
            while let Some(c) = stack.pop() {
                let d = point_segment(q, c.p, c.q);
                if d <= c.band + margin {
                    match c.course {
                        Some((k, t0, t1)) if c.sag > fine => {
                            let trace = &self.traces[k];
                            let tm = 0.5 * (t0 + t1);
                            let m = trace.at(tm);
                            for (p, a, b, s) in [(c.p, t0, tm, m), (m, tm, t1, c.q)] {
                                let moved =
                                    (trace.at(a) - p).length().max((trace.at(b) - s).length());
                                let sag = trace.band(a, b);
                                stack.push(OChord {
                                    p,
                                    q: s,
                                    band: sag + moved,
                                    sag,
                                    course: Some((k, a, b)),
                                });
                            }
                            continue;
                        }
                        _ => {
                            if d <= margin.max(c.band) {
                                return None;
                            }
                        }
                    }
                }
                if (c.p.y > q.y) != (c.q.y > q.y) {
                    let x = c.p.x + (q.y - c.p.y) / (c.q.y - c.p.y) * (c.q.x - c.p.x);
                    if x > q.x {
                        inside = !inside;
                    }
                }
            }
        }
        Some(inside)
    }
}

#[cfg(test)]
mod tests {
    use super::super::testing::Rng;
    use super::*;

    /// A random unit vector square to the unit `axis`.
    fn across(rng: &mut Rng, axis: DVec3) -> DVec3 {
        loop {
            let v = rng.unit();
            let w = v - axis * v.dot(axis);
            if w.length() > 0.1 {
                return w.normalize();
            }
        }
    }

    /// A random cone's chart, oriented for an outward face, and its
    /// surface.
    fn cone(rng: &mut Rng) -> (Chart, Quadric) {
        let point = rng.point(2.0);
        let axis = rng.unit();
        let x = across(rng, axis);
        let radius = rng.range(0.1, 2.0);
        // Half-angles from slight tapers to wide cones.
        let alpha = [0.01, 0.2, 0.7, 1.3][rng.below(4)] * rng.range(0.8, 1.2);
        let (sin, cos) = alpha.sin_cos();
        let chart = Chart::Cone {
            point,
            axis,
            x,
            y: x.cross(axis),
            radius,
            sin,
            cos,
        };
        let surface = Quadric::Cone {
            point,
            axis,
            radius,
            slope: sin / cos,
        };
        (chart, surface)
    }

    /// The point of a cone's chart at slant `sigma` from its reference
    /// circle and angle `theta` about its axis from the chart's middle.
    fn on_cone(chart: &Chart, sigma: f64, theta: f64) -> DVec3 {
        let Chart::Cone {
            point,
            axis,
            x,
            y,
            radius,
            sin,
            cos,
        } = *chart
        else {
            unreachable!()
        };
        point + axis * (sigma * cos) + (x * theta.cos() + y * theta.sin()) * (radius + sigma * sin)
    }

    /// A random point of a cone's face: within a quarter turn of the
    /// chart's middle, on the nappe, within a few radii of the reference.
    fn cone_point(rng: &mut Rng, chart: &Chart) -> DVec3 {
        let Chart::Cone { radius, sin, .. } = *chart else {
            unreachable!()
        };
        let sigma = rng.range((-radius / sin).max(-3.0), 3.0);
        on_cone(chart, sigma, rng.range(-FRAC_PI_2, FRAC_PI_2))
    }

    fn sphere(rng: &mut Rng) -> (Chart, Quadric) {
        let center = rng.point(2.0);
        let radius = rng.range(0.1, 3.0);
        let c = rng.unit();
        let x = across(rng, c);
        let chart = Chart::Sphere {
            center,
            radius,
            c,
            x,
            y: c.cross(x),
        };
        (chart, Quadric::Sphere { center, radius })
    }

    /// A random point of a sphere's face: within 120° of the chart's
    /// middle (the most a face may reach past its half, `shape::HALF`).
    fn sphere_point(rng: &mut Rng, chart: &Chart) -> DVec3 {
        let Chart::Sphere {
            center, radius, c, ..
        } = *chart
        else {
            unreachable!()
        };
        loop {
            let u = rng.unit();
            if u.dot(c) >= -0.5 {
                return center + u * radius;
            }
        }
    }

    fn cylinder(rng: &mut Rng) -> (Chart, Quadric) {
        let point = rng.point(2.0);
        let axis = rng.unit();
        let x = across(rng, axis);
        let radius = rng.range(0.1, 2.0);
        let chart = Chart::Cylinder {
            point,
            axis,
            x,
            y: axis.cross(x),
            radius,
        };
        (
            chart,
            Quadric::Cylinder {
                point,
                axis,
                radius,
            },
        )
    }

    fn cylinder_point(rng: &mut Rng, chart: &Chart) -> DVec3 {
        let Chart::Cylinder { x, y, .. } = *chart else {
            unreachable!()
        };
        let theta = rng.range(-FRAC_PI_2, FRAC_PI_2);
        let q = DVec2::new(0.0, rng.range(-3.0, 3.0));
        let p = chart.unmap(q);
        let Chart::Cylinder {
            point,
            axis,
            radius,
            ..
        } = *chart
        else {
            unreachable!()
        };
        let h = (p - point).dot(axis);
        point + axis * h + (x * theta.cos() + y * theta.sin()) * radius
    }

    type Maker = fn(&mut Rng) -> (Chart, Quadric);
    type Sampler = fn(&mut Rng, &Chart) -> DVec3;
    const KINDS: [(&str, Maker, Sampler); 3] = [
        ("cylinder", cylinder, cylinder_point),
        ("cone", cone, cone_point),
        ("sphere", sphere, sphere_point),
    ];

    /// Points go to chart points and back; the cone's chart keeps lengths
    /// and the sphere's scales them by `1 + |q|² / (4 R²)` the same way in
    /// every direction.
    #[test]
    fn charts_round_trip_and_keep_or_scale_lengths() {
        let mut rng = Rng(0xC4A2_7001);
        for (name, make, sample) in KINDS {
            for _ in 0..200 {
                let (chart, surface) = make(&mut rng);
                let size = surface.magnitude();
                for _ in 0..20 {
                    let p = sample(&mut rng, &chart);
                    let q = chart.map(p);
                    let back = chart.unmap(q);
                    assert!(
                        (back - p).length() <= 1e-12 * (1.0 + size),
                        "{name}: {p} {back}"
                    );
                    assert!(
                        (chart.map(back) - q).length() <= 1e-12 * (1.0 + size),
                        "{name}"
                    );
                    // Two steps along the surface, square to each other in
                    // the chart.
                    let h = 1e-5;
                    let ex = chart.unmap(q + DVec2::X * h) - p;
                    let ey = chart.unmap(q + DVec2::Y * h) - p;
                    let scale = match chart {
                        Chart::Sphere { radius, .. } => {
                            1.0 + q.length_squared() / (4.0 * radius * radius)
                        }
                        _ => 1.0,
                    };
                    for e in [ex, ey] {
                        let ratio = h / e.length();
                        assert!(
                            (ratio - scale).abs() <= 1e-4 * scale,
                            "{name}: {ratio} vs {scale}"
                        );
                    }
                    assert!(ex.dot(ey).abs() <= 1e-4 * h * h, "{name}: not conformal");
                    // Outward: the chart's x and y turn about the surface's
                    // outward normal; reversed, about the inward one.
                    let n = surface.gradient(p);
                    assert!(ex.cross(ey).dot(n) > 0.0, "{name}: chart faces inward");
                    let r = chart.reversed();
                    let qr = r.map(p);
                    let fx = r.unmap(qr + DVec2::X * h) - p;
                    let fy = r.unmap(qr + DVec2::Y * h) - p;
                    assert!(
                        fx.cross(fy).dot(n) < 0.0,
                        "{name}: reversed chart faces outward"
                    );
                }
            }
        }
    }

    /// `space` and `reach` bound distances in space by distances in the
    /// chart and back, for pairs of points of faces.
    #[test]
    fn reach_and_space_bound_distances_in_space() {
        let mut rng = Rng(0xC4A2_7002);
        let mut tight = [0.0f64; 3];
        for (k, (name, make, sample)) in KINDS.into_iter().enumerate() {
            for _ in 0..400 {
                let (chart, _) = make(&mut rng);
                for _ in 0..25 {
                    let (p, s) = (sample(&mut rng, &chart), sample(&mut rng, &chart));
                    // Near pairs too: a point toward the other.
                    let s = if rng.next() < 0.5 {
                        chart.unmap(
                            chart.map(p) + (chart.map(s) - chart.map(p)) * rng.range(1e-6, 0.1),
                        )
                    } else {
                        s
                    };
                    let (q, r) = (chart.map(p), chart.map(s));
                    let (d, space) = ((r - q).length(), (s - p).length());
                    if d == 0.0 {
                        continue;
                    }
                    let low = chart.space(q, d);
                    assert!(
                        space >= low * (1.0 - 1e-9),
                        "{name}: space {space} below {low}"
                    );
                    let reach = chart.reach(q, space);
                    assert!(
                        d <= reach * (1.0 + 1e-9),
                        "{name}: chart {d} beyond reach {reach}"
                    );
                    tight[k] = tight[k].max(d / reach);
                }
            }
        }
        // The bounds are not wildly loose: some pairs come near them.
        for (k, t) in tight.iter().enumerate() {
            assert!(*t > 0.5, "{}: reach loose ({t})", KINDS[k].0);
        }
    }

    /// The largest distance of `n` samples of the course from its chord.
    fn strays(trace: &Trace, t0: f64, t1: f64) -> f64 {
        let (a, b) = (trace.at(t0), trace.at(t1));
        (0..=200)
            .map(|i| point_segment(trace.at(t0 + (t1 - t0) * i as f64 / 200.0), a, b))
            .fold(0.0, f64::max)
    }

    /// Circles about a cone's axis and circles of a sphere stay within
    /// their chords' bands, over every first chord and over halves and
    /// quarters of them.
    #[test]
    fn bands_hold_courses_of_circles() {
        let mut rng = Rng(0xC4A2_7003);
        for _ in 0..300 {
            let (chart, _) = if rng.next() < 0.5 {
                cone(&mut rng)
            } else {
                sphere(&mut rng)
            };
            let track = match chart {
                Chart::Cone {
                    point,
                    axis,
                    x,
                    y,
                    radius,
                    sin,
                    cos,
                } => {
                    let sigma = rng.range((-radius / sin).max(-3.0) * 0.9, 3.0);
                    Track::Circle {
                        center: point + axis * (sigma * cos),
                        x,
                        y,
                        radius: radius + sigma * sin,
                    }
                }
                Chart::Sphere {
                    center, radius, c, ..
                } => {
                    // A plane across the sphere, its circle within the
                    // face's reach of the chart's middle.
                    let n = loop {
                        let n = rng.unit();
                        if n.dot(c) > 0.3 {
                            break n;
                        }
                    };
                    let off = radius * rng.range(-0.2, 0.95);
                    let x = across(&mut rng, n);
                    Track::Circle {
                        center: center + n * off,
                        x,
                        y: n.cross(x),
                        radius: (radius * radius - off * off).sqrt(),
                    }
                }
                _ => unreachable!(),
            };
            let trace = trace(&chart, TraceOf::Track(track)).unwrap().unwrap();
            // A stretch of the circle within the face.
            let t0 = rng.range(-1.0, 0.0);
            let t1 = t0 + rng.range(0.1, 2.0);
            let ts = trace.division(t0, t1).unwrap();
            for w in ts.windows(2) {
                let (a, b) = (w[0], w[1]);
                for (u, v) in [
                    (a, b),
                    (a, 0.5 * (a + b)),
                    (0.5 * (a + b), b),
                    (a, 0.75 * a + 0.25 * b),
                ] {
                    let band = trace.band(u, v);
                    let off = strays(&trace, u, v);
                    assert!(
                        off <= band,
                        "{chart:?} {track:?}: strays {off} beyond band {band}"
                    );
                }
            }
        }
    }
}
