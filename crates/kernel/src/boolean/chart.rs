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
//! control points, whose convex hull holds the piece; on a cylinder or a
//! cone, by the Bernstein coefficients of its derivatives (§16.2).

use std::f64::consts::{FRAC_PI_2, FRAC_PI_8, PI};

use glam::{DVec2, DVec3};
use opendatum_nurbs::BezierSegment;

use super::meet::{self, Track};
use super::planar::{halvable, point_segment, Bend, Sag};
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

    /// Size of the numbers [`Chart::map`] subtracts and scales by (the
    /// plane's origin, the surface's point plus its radius, a sphere's
    /// centre plus four radii, as far as a face's chart reaches): with the
    /// size of the point mapped, this sets the rounding of its chart point,
    /// as [`Curve3::magnitude`] does a curve's points'.
    pub fn magnitude(&self) -> f64 {
        match *self {
            Chart::Plane(f) => f.origin.abs().max_element(),
            Chart::Cylinder { point, radius, .. } | Chart::Cone { point, radius, .. } => {
                point.abs().max_element() + radius
            }
            Chart::Sphere { center, radius, .. } => center.abs().max_element() + 4.0 * radius,
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

    /// The cylinder or the cone a chart develops onto the plane.
    fn developed(&self) -> Option<Quadric> {
        match *self {
            Chart::Cylinder {
                point,
                axis,
                radius,
                ..
            } => Some(Quadric::Cylinder {
                point,
                axis,
                radius,
            }),
            Chart::Cone {
                point,
                axis,
                radius,
                sin,
                cos,
                ..
            } => Some(Quadric::Cone {
                point,
                axis,
                radius,
                slope: sin / cos,
            }),
            Chart::Plane(_) | Chart::Sphere { .. } => None,
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
    /// A NURBS curve: in a plane's chart bounded by its Bézier pieces'
    /// control points, in a cylinder's or a cone's by their derivatives
    /// (design note 0012 §16.2).
    Nurbs(&'a Curve3),
}

/// A curved segment's course in a chart.
pub(super) struct Trace<'a> {
    chart: Chart,
    path: Path<'a>,
    /// For a NURBS curve, its Bézier segments; empty for a track.
    segments: Vec<BezierSegment>,
    /// A bound on the rounding of chart points and of the bounds' terms,
    /// added to every band.
    slack: f64,
    /// For a NURBS curve in a plane's chart, the allowance in its bend and
    /// the bend's floor ([`Sag::floor`]): a bound on what the rounding of
    /// its clipped control points, of its chord's ends and of their
    /// distances adds to the distances computed, however short the stretch
    /// ([`hull_rounding`]). Zero otherwise.
    rounding: f64,
}

/// The allowance for the rounding of a conic's chart points across its
/// course, in units of `ε` times the stretch's reach from the vertex: at
/// least the first-order bounds (7.5 in a plane's chart, 27.2 in a
/// cone's; design note 0012 §16.2), and 16 times the most measured.
const CONIC_ROUNDING: f64 = 32.0;

/// A bound, to first order, on what rounding adds to the distances of a
/// NURBS curve's clipped control points from a chord in the plane chart
/// with origin `origin` ([`Trace::bend`]), however short the stretch
/// (design note 0012 §16.2). With `R` the control points' largest distance
/// from the world's origin, `D` their diameter, `d` their largest distance
/// from the chart's origin, `p` the degree and `ρ` the ratio of the largest
/// weight to the least, in units of ε:
///
/// - the clipped pieces' points, two runs of de Casteljau's `p` steps:
///   the products and sums in each step round by ε of the terms, which
///   the hull bounds by `R` times the weight (the weight cancels when the
///   point is divided by it), and `1 − s` by ε/2 of the coefficient, which
///   moves the point by at most that part of `D`: `(4p + ½) R + p D`;
/// - the clip's parameters, rounded before `1 − s` is taken from them:
///   moving a split by `δ` moves a point by up to `(1 + ρ) D δ`, for at
///   most `2p (1 + ρ) D`;
/// - splitting the curve into its pieces: `2p R + 2p D` the same way, and
///   (sketched, not derived) `1.5p (1 + ρ) D` for its parameters;
/// - the chord's ends, evaluated: `(p + 3/2) R + 3p D`;
/// - mapping the points to the chart and their distances from the chord:
///   `4.5 √2 d`.
///
/// In all `(7p + 2) R + (6p + 3.5p (1 + ρ)) D + 6.4 d`; measured at most
/// 0.47 of it over ten million stretches (degrees 2 to 5, weights' ratios
/// to 10⁴, from 1 m to 10¹⁰ m from the origin, down to 10⁻¹² of a knot
/// span; `hull_rounding`).
fn hull_rounding(n: &opendatum_nurbs::NurbsCurve, origin: DVec3) -> f64 {
    let points = n.control_points();
    let weights = n.weights();
    let (least, most) = weights
        .iter()
        .fold((f64::INFINITY, 0.0f64), |(l, m), &w| (l.min(w), m.max(w)));
    let (lo, hi) = points.iter().fold(
        (DVec3::splat(f64::INFINITY), DVec3::splat(f64::NEG_INFINITY)),
        |(lo, hi), q| (lo.min(*q), hi.max(*q)),
    );
    // The box's diagonal bounds the diameter.
    let diameter = (hi - lo).length();
    let (r, d) = points.iter().fold((0.0f64, 0.0f64), |(r, d), q| {
        (r.max(q.length()), d.max((*q - origin).length()))
    });
    let p = n.degree() as f64;
    let rho = most / least;
    f64::EPSILON * ((7.0 * p + 2.0) * r + (6.0 * p + 3.5 * p * (1.0 + rho)) * diameter + 6.4 * d)
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
    /// the turning of its tangent). A NURBS curve's are two per Bézier
    /// piece, so that no chord spans a joint, halved while a chord's bend
    /// is more than a quarter of its length, up to [`MOST_CHORDS`]
    /// (design note 0012 §16.2).
    pub fn division(&self, t0: f64, t1: f64) -> Result<Vec<f64>, KernelError> {
        let mut out = vec![t0];
        match self.path {
            Path::Track(track, bound) => {
                if let Bound::Conic { c, .. } = bound {
                    // A conic turns less than a half turn: at most eight.
                    let (f0, f1) = (meet::conic_turn(c, t0), meet::conic_turn(c, t1));
                    let n = ((f1 - f0).abs() / FRAC_PI_8).ceil().clamp(1.0, 8.0) as usize;
                    for k in 1..n {
                        out.push(meet::conic_at_turn(c, f0 + (f1 - f0) * k as f64 / n as f64));
                    }
                } else {
                    let span = (t1 - t0).abs();
                    let n = match bound {
                        Bound::Steady(0.0) => 1.0,
                        Bound::Arc(kappa) => (span / FRAC_PI_8)
                            .ceil()
                            .max((kappa * self.chart_length(&track, t0, t1) / FRAC_PI_2).ceil()),
                        _ => (span / FRAC_PI_8).ceil(),
                    };
                    // A course in a face turns at most a few times (a
                    // sphere's face lies within a half of it, where the
                    // chart's scale is at most 2): more chords than this is
                    // a course that runs toward the chart's cut, not one in
                    // the face.
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
            }
            Path::Nurbs(_) => {
                let (lo, hi) = if t0 <= t1 { (t0, t1) } else { (t1, t0) };
                let mut ks: Vec<f64> = Vec::new();
                for (a, b, _) in quadric::spans(&self.segments, lo, hi) {
                    if a > lo {
                        ks.push(a);
                    }
                    ks.push(0.5 * (a + b));
                }
                ks.retain(|&t| t > lo && t < hi);
                ks.dedup();
                if t0 > t1 {
                    ks.reverse();
                }
                out.extend(ks);
                out.push(t1);
                return self.refined(out);
            }
        }
        out.push(t1);
        Ok(out)
    }

    /// The parameters `first` with every chord between them whose bend is
    /// more than a quarter of its length halved, in rounds (each halving
    /// every such chord of the round before, so that the chords are spent
    /// evenly along the course), until there are [`MOST_CHORDS`] or the
    /// bend is down to its rounding (four times its [`Sag::floor`]): the
    /// outline's widest band bounds the searches of point location. A
    /// chord whose bend is still unbounded then (a course that runs to a
    /// cone's apex or a cylinder's axis, where no face's edge can) is
    /// refused.
    fn refined(&self, first: Vec<f64>) -> Result<Vec<f64>, KernelError> {
        // Each chord with whether it is settled.
        let mut chords: Vec<(f64, f64, bool)> =
            first.windows(2).map(|w| (w[0], w[1], false)).collect();
        let mut count = chords.len();
        loop {
            let mut next = Vec::with_capacity(2 * chords.len());
            let mut halved = false;
            for &(a, b, settled) in &chords {
                if settled {
                    next.push((a, b, true));
                    continue;
                }
                let bend = self.bend(a, b);
                let m = 0.5 * (a + b);
                let len = (self.at(b) - self.at(a)).length();
                let wide = bend.bound.is_nan()
                    || len.is_nan()
                    || bend.bound > (0.25 * len).max(4.0 * bend.floor);
                if wide && (count as f64) < MOST_CHORDS && m != a && m != b {
                    next.push((a, m, false));
                    next.push((m, b, false));
                    count += 1;
                    halved = true;
                    continue;
                }
                if !bend.bound.is_finite() {
                    return Err(KernelError::InvalidResult(
                        "boolean: a curved edge's course in its face's chart could not be \
                         bounded (it runs toward the surface's axis or apex)"
                            .into(),
                    ));
                }
                next.push((a, b, true));
            }
            chords = next;
            if !halved {
                break;
            }
        }
        let mut out = Vec::with_capacity(chords.len() + 1);
        out.push(first[0]);
        out.extend(chords.iter().map(|&(_, b, _)| b));
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

impl Trace<'_> {
    /// A bound on the distance from the course on `[t0, t1]` (either
    /// order) to the chord between its ends, but for the rounding of chart
    /// points ([`Trace::slack`]), with the part of it that halving does not
    /// narrow ([`Sag::floor`]).
    fn bend(&self, t0: f64, t1: f64) -> Sag {
        let exact = |bound| Sag {
            bound,
            floor: 0.0,
            settles: false,
        };
        match self.path {
            Path::Track(_, Bound::Steady(m)) => exact(m * (t1 - t0) * (t1 - t0) / 8.0),
            Path::Track(track, Bound::Conic { p, c }) => {
                // w is least at the greatest |t| on half an ellipse, else at
                // the |t| nearest zero. The rounding of chart points grows
                // with the stretch's reach from the vertex (the part that
                // grows with the vertex's and the chart's coordinates is the
                // slack's): across the course, by at most 7.5 ε reach in a
                // plane's chart and 27.2 ε reach in a cone's, to first order
                // (design note 0012 §16.2); measured at most 1.9 ε reach
                // where the reach is a metre or more (1.4 million points,
                // 2026-10-15; `rounding::conic_courses_lie_within_their_bands`
                // checks the bands). The allowance is 32 ε reach, reach the
                // farthest the stretch's ends are from the vertex (the
                // distance grows with |t|). The other term is a product of
                // the stretch's parameters, rounded only in proportion to
                // it: the bound settles on the allowance ([`Sag::settles`]).
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
                let floor = CONIC_ROUNDING * f64::EPSILON * reach;
                Sag {
                    bound: p / (w * w * w) * (t1 - t0) * (t1 - t0) / 8.0 + floor,
                    floor,
                    settles: true,
                }
            }
            Path::Track(_, Bound::Arc(kappa)) => {
                // The sagitta of an arc of curvature κ less than a half turn
                // over a chord of length 2h: h² κ / (1 + √(1 − (κ h)²)).
                let h = 0.5 * (self.at(t1) - self.at(t0)).length();
                let x = (kappa * h).min(1.0);
                exact(h * x / (1.0 + (1.0 - x * x).sqrt()))
            }
            Path::Nurbs(_) => match self.chart.developed() {
                None => {
                    // The convex hull of each Bézier piece's control points
                    // holds the piece (weights are positive), and the chart
                    // is affine: the hull's farthest point from the chord
                    // bounds the course's. The pieces are the trace's own,
                    // split from the curve once ([`quadric::hull_points`]
                    // split it again on every call).
                    // The distances are computed from rounded points; the
                    // allowance [`Trace::rounding`] bounds what that adds.
                    let (p, q) = (self.at(t0), self.at(t1));
                    let hull = quadric::spans(&self.segments, t0.min(t1), t0.max(t1))
                        .flat_map(|(_, _, pw)| pw)
                        .map(|w| point_segment(self.chart.map(w.truncate() / w.w), p, q))
                        .fold(0.0, f64::max);
                    Sag {
                        bound: hull + self.rounding,
                        floor: self.rounding,
                        settles: false,
                    }
                }
                Some(surface) => self.developed_bend(&surface, t0, t1),
            },
        }
    }

    /// [`Trace::bend`] of a NURBS curve in a cylinder's or a cone's chart
    /// (design note 0012 §16.2). On each Bézier piece's stretch the course
    /// is within an eighth of its second derivative's bound of the chord of
    /// that stretch (the error of linear interpolation); the chord of each
    /// stretch is within the largest distance of its ends from the whole
    /// chord, at the same parameters, of it (their difference is affine
    /// in the parameter). The floor is the pieces' ([`quadric::developed_bend`]);
    /// the joints' distances narrow with the chord to the rounding of chart
    /// points, a few units in the last place of their coordinates, well
    /// within the margin between the three floors the bound comes down to
    /// and the four [`Sag::wide`] stops at.
    fn developed_bend(&self, surface: &Quadric, t0: f64, t1: f64) -> Sag {
        let (lo, hi) = if t0 <= t1 { (t0, t1) } else { (t1, t0) };
        if lo == hi {
            return Sag::ZERO;
        }
        let (p, q) = (self.at(t0), self.at(t1));
        let mut bend = Sag::ZERO;
        let mut joints = 0.0f64;
        for (a, _, pw) in quadric::spans(&self.segments, lo, hi) {
            let Some((m, floor)) = quadric::developed_bend(surface, &pw) else {
                return Sag {
                    bound: f64::INFINITY,
                    floor: 0.0,
                    settles: false,
                };
            };
            bend.bound = bend.bound.max(m / 8.0);
            bend.floor = bend.floor.max(floor / 8.0);
            if a > lo {
                let on_chord = p + (q - p) * ((a - t0) / (t1 - t0));
                joints = joints.max((self.at(a) - on_chord).length());
            }
        }
        bend.bound += joints;
        bend
    }
}

impl Bend for Trace<'_> {
    fn at(&self, t: f64) -> DVec2 {
        self.chart.map(self.point(t))
    }

    fn sag(&self, t0: f64, t1: f64) -> Sag {
        self.bend(t0, t1)
    }

    fn slack(&self) -> f64 {
        self.slack
    }
}

/// The course of a circle, an ellipse or a NURBS curve in `chart`, or
/// `None` for a course that is straight there (a line; a ruling or a
/// circle square to the axis on a cylinder; a ruling on a cone).
pub(super) fn trace<'a>(
    chart: &Chart,
    path: TraceOf<'a>,
) -> Result<Option<Trace<'a>>, KernelError> {
    let track = match path {
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
            Curve3::Nurbs(n) => {
                return match chart {
                    Chart::Sphere { .. } => Err(KernelError::Unsupported(
                        "booleans with a NURBS edge on a sphere, where no conic lies \
                             (design note 0012 §16.2)"
                            .into(),
                    )),
                    _ => Ok(Some(Trace {
                        chart: *chart,
                        path: Path::Nurbs(c),
                        segments: n.bezier_segments(),
                        slack: 1e-13 * (1.0 + c.magnitude() + chart.magnitude()),
                        rounding: match chart {
                            Chart::Plane(f) => hull_rounding(n, f.origin),
                            _ => 0.0,
                        },
                    })),
                };
            }
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
        segments: Vec::new(),
        slack: 1e-13 * (1.0 + size + chart.magnitude()),
        rounding: 0.0,
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
    /// The part of `band` from the course's bend, which halving narrows
    /// down to its floor; the rest is the rounding allowance and how far
    /// the chord's ends are from the course's.
    sag: Sag,
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
                        sag: Sag::ZERO,
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
                            let sag = c.sag(ts[k], ts[k + 1]);
                            chords.push(OChord {
                                p: at[k],
                                q: at[k + 1],
                                band: sag.bound + c.slack + moved,
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
    /// `fine` (or held over it by their rounding, [`halvable`]) or their
    /// bands clear of `q`; the answer is the parity of a ray from `q`
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
                    // Halved while that can narrow its sag toward `fine`.
                    let halves = c
                        .course
                        .filter(|&(_, t0, t1)| halvable(c.sag, t0, t1, fine));
                    match halves {
                        Some((k, t0, t1)) => {
                            let trace = &self.traces[k];
                            let tm = 0.5 * (t0 + t1);
                            let m = trace.at(tm);
                            for (p, a, b, s) in [(c.p, t0, tm, m), (m, tm, t1, c.q)] {
                                let moved =
                                    (trace.at(a) - p).length().max((trace.at(b) - s).length());
                                let sag = trace.sag(a, b);
                                stack.push(OChord {
                                    p,
                                    q: s,
                                    band: sag.bound + trace.slack + moved,
                                    sag,
                                    course: Some((k, a, b)),
                                });
                            }
                            continue;
                        }
                        None => {
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
    use super::super::planar::SETTLED;
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

    /// The chart `chart` (a cylinder's or a cone's) turned about its axis to
    /// face `p`, and the angle about the axis of a point from its middle.
    fn facing(chart: &Chart, p: DVec3) -> Chart {
        match *chart {
            Chart::Cylinder {
                point,
                axis,
                radius,
                ..
            } => {
                let x = across_of(point, axis, p);
                Chart::Cylinder {
                    point,
                    axis,
                    x,
                    y: axis.cross(x),
                    radius,
                }
            }
            Chart::Cone {
                point,
                axis,
                radius,
                sin,
                cos,
                ..
            } => {
                let x = across_of(point, axis, p);
                Chart::Cone {
                    point,
                    axis,
                    x,
                    y: x.cross(axis),
                    radius,
                    sin,
                    cos,
                }
            }
            _ => unreachable!(),
        }
    }

    fn across_of(point: DVec3, axis: DVec3, p: DVec3) -> DVec3 {
        let w = p - point;
        (w - axis * w.dot(axis)).normalize()
    }

    /// The angle about a cylinder's or a cone's axis of `p` from the
    /// chart's middle.
    fn angle_in(chart: &Chart, p: DVec3) -> f64 {
        let (Chart::Cylinder {
            point, axis, x, y, ..
        }
        | Chart::Cone {
            point, axis, x, y, ..
        }) = *chart
        else {
            unreachable!()
        };
        let w = p - point;
        let a = w - axis * w.dot(axis);
        a.dot(y).atan2(a.dot(x))
    }

    /// A stretch of a plane's section of a random cylinder or cone (`kind`
    /// 0 or 1) through a point of a face of it, as the booleans write such
    /// an edge (a rational quadratic NURBS curve on `[0, 1]`), with the
    /// chart turned to face its middle; the stretch keeps within a little
    /// less than a quarter turn of it, as a face's edge keeps within half
    /// a turn. `None` when the plane misses or meets in rulings or a
    /// circle.
    fn rational_edge(rng: &mut Rng, kind: usize) -> Option<(Chart, Quadric, Curve3, Track)> {
        let (chart, surface) = if kind == 0 { cylinder(rng) } else { cone(rng) };
        let p = if kind == 0 {
            cylinder_point(rng, &chart)
        } else {
            cone_point(rng, &chart)
        };
        let chart = facing(&chart, p);
        let size = surface.magnitude() + 3.0;
        let plane = Quadric::Plane {
            point: p,
            normal: rng.unit(),
        };
        let region = meet::Region::of_box(p - DVec3::splat(size), p + DVec3::splat(size));
        let meet::Meet::Curves(tracks) = meet::meet(&plane, &surface, &region, 1e-7) else {
            return None;
        };
        let track = tracks
            .into_iter()
            .filter(|t| matches!(t, Track::Ellipse { .. } | Track::Conic { .. }))
            .min_by(|a, b| a.nearest(p).1.total_cmp(&b.nearest(p).1))?;
        let (tp, d) = track.nearest(p);
        if d > 1e-9 * size {
            return None;
        }
        let (lo_end, hi_end) = track.ends();
        let reach = match track {
            Track::Conic { p: semi, .. } => 4.0 * size / semi,
            _ => 1.5,
        };
        let (mut a, mut b) = (rng.range(0.0, reach), rng.range(0.0, reach));
        for _ in 0..60 {
            let (lo, hi) = ((tp - a).max(0.999 * lo_end), (tp + b).min(0.999 * hi_end));
            let inside = (0..=64).all(|k| {
                let q = track.point(lo + (hi - lo) * k as f64 / 64.0);
                angle_in(&chart, q).abs() < 0.45 * PI
                    && match surface {
                        Quadric::Cone {
                            point,
                            axis,
                            radius,
                            slope,
                        } => radius + slope * (q - point).dot(axis) > 1e-3 * radius,
                        _ => true,
                    }
            });
            if inside && hi - lo > 1e-6 {
                let (curve, _, _) = track.edge(lo, hi).ok()?;
                return Some((chart, surface, curve, track));
            }
            a *= 0.7;
            b *= 0.7;
        }
        None
    }

    /// The curve with its control points moved by up to `e` and its
    /// weights scaled by up to `1 ± e` (relative to `size`), off the
    /// surface.
    fn pushed(rng: &mut Rng, curve: &Curve3, e: f64, size: f64) -> Curve3 {
        let Curve3::Nurbs(n) = curve else {
            unreachable!()
        };
        let points = n
            .control_points()
            .iter()
            .map(|p| *p + rng.point(e * size))
            .collect();
        let weights = n
            .weights()
            .iter()
            .map(|w| w * (1.0 + rng.range(-e, e)))
            .collect();
        Curve3::Nurbs(
            opendatum_nurbs::NurbsCurve::new(
                n.degree(),
                n.knots().as_slice().to_vec(),
                points,
                Some(weights),
            )
            .unwrap(),
        )
    }

    /// Plane sections of cylinders and cones written as rational NURBS
    /// curves, as the booleans write a conic edge, and the same pushed off
    /// the surface: over their first chords, halves and quarters of them,
    /// and random stretches either way, their courses in the face's chart
    /// stay within the bands (design note 0012 §16.2), and the bands are
    /// not far looser than they need be on the surface.
    #[test]
    fn bands_hold_courses_of_rational_edges() {
        let mut rng = Rng(0xC4A2_7004);
        let mut kinds = [0usize; 2];
        let mut conics = 0;
        let mut chords = 0;
        let mut loosest = 0.0f64;
        for case in 0..1200 {
            let kind = case % 2;
            let Some((chart, surface, curve, track)) = rational_edge(&mut rng, kind) else {
                continue;
            };
            kinds[kind] += 1;
            conics += usize::from(matches!(track, Track::Conic { c, .. } if c <= 0.0));
            let size = surface.magnitude() + curve.magnitude();
            let e = [0.0, 1e-9, 1e-6, 1e-3][case / 2 % 4];
            let curve = if e > 0.0 {
                pushed(&mut rng, &curve, e, size)
            } else {
                curve
            };
            let trace = trace(&chart, TraceOf::Curve(&curve)).unwrap().unwrap();
            let (u, v) = (rng.next(), rng.next());
            for (t0, t1) in [
                (0.0, 1.0),
                (1.0, 0.0),
                (u.min(v), u.max(v)),
                (u.max(v), u.min(v)),
            ] {
                if t0 == t1 {
                    continue;
                }
                let ts = trace.division(t0, t1).unwrap();
                assert!(ts.len() >= 2 && ts[0] == t0 && ts[ts.len() - 1] == t1);
                for w in ts.windows(2) {
                    let (a, b) = (w[0], w[1]);
                    assert!((a < b) == (t0 < t1), "division out of order");
                    for (x, y) in [
                        (a, b),
                        (a, 0.5 * (a + b)),
                        (0.5 * (a + b), b),
                        (a, 0.75 * a + 0.25 * b),
                        (b, a),
                    ] {
                        let band = trace.sag(x, y).bound + trace.slack;
                        let off = strays(&trace, x, y);
                        assert!(
                            off <= band,
                            "case {case}: {track:?} pushed {e}: strays {off} beyond band {band}"
                        );
                        if e == 0.0 && off > 1e-9 * size {
                            loosest = loosest.max(band / off);
                        }
                        chords += 1;
                    }
                    // The division leaves no chord whose bend is more than a
                    // quarter of its length (else it hit the cap, or the
                    // bend is down to its rounding).
                    let len = (trace.at(b) - trace.at(a)).length();
                    let sag = trace.sag(a, b);
                    assert!(
                        sag.bound <= (0.25 * len).max(4.0 * sag.floor)
                            || ts.len() as f64 > MOST_CHORDS,
                        "case {case}: a wide first chord"
                    );
                }
            }
        }
        assert!(kinds[0] >= 300 && kinds[1] >= 300, "{kinds:?}");
        assert!(conics >= 100, "{conics} parabolas and hyperbolas");
        assert!(chords > 10_000);
        // Measured on the surface (2026-10-15): the bands are 4.7 times
        // what the samples stray at the median, 140 at the 90th percentile,
        // and up to 7.4 × 10⁵ far along a hyperbola's branch, where the
        // speed of the curve's parameter changes most.
        assert!(loosest < 1e7, "bands loose by {loosest}");
    }

    /// Halving a stretch of a curved course narrows its sag to under four
    /// times its floor, the part its rounding makes (a conic track's, which
    /// settles on its floor, to within a thousandth of it, [`Sag::settles`]),
    /// before the parameter's last bit, so that halving toward a sixteenth
    /// of the tolerance stops ([`halvable`]) wherever the course lies:
    /// NURBS ellipses in planes and on cylinders, NURBS circles on cones,
    /// conic tracks (some stretches 10⁶ m and more from their vertices) and
    /// rational Bézier pieces as large as their distance from the origin,
    /// with weights' ratios up to 10⁴, from the origin to 10⁸ m from it,
    /// and the hyperbolas of cones cut a millimetre off their axes
    /// (weights' ratios in the hundreds) to 5·10⁵ m (meet refuses conics
    /// farther out). Without the floors the halving went on to the
    /// parameter's last bit far from the origin, the chords doubling
    /// (found by review). Where it stops, the sags are under the sixteenth
    /// near the origin: within 10² m of it, all but the Bézier pieces,
    /// which are at the origin (up to 2 m across), and the hyperbolas up to
    /// 10⁴ m from it (a floor in proportion to the weights' ratio held them
    /// over it from 10² m, found by the third review); and a conic track's
    /// wherever its floor is under it by a thousandth.
    #[test]
    fn sags_narrow_to_their_floors() {
        let fine = 1e-7 / 16.0;
        let mut rng = Rng(0xC4A2_7010);
        // Per kind, the most a sag is of its floor once its stretch is
        // tiny, and of the coordinates' size times ε.
        let mut most = [[0.0f64; 2]; 6];
        for case in 0..1200 {
            let far = [0.0, 1e2, 1e4, 1e5, 1e6, 1e7, 1e8][case % 7];
            let kind = (case / 7) % 6;
            // Whether the stretch is far from where its course is small
            // enough to resolve to the sixteenth (far from a conic's
            // vertex, or on a large unevenly weighted piece).
            let mut distant = false;
            let shift = rng.unit() * far;
            let (delta, c) = (rng.range(0.2, 2.5), rng.range(0.3, 3.0));
            // A rational quadratic arc of `delta` about `center`, the
            // image of a circle's by `map`, reweighted (1, w c, c²): the
            // same curve, followed unevenly.
            let arc = |center: DVec3, map: &dyn Fn(f64, f64) -> DVec3, from: f64| {
                let w = (0.5 * delta).cos();
                let mid = from + 0.5 * delta;
                let points = vec![
                    center + map(from.cos(), from.sin()),
                    center + map(mid.cos() / w, mid.sin() / w),
                    center + map((from + delta).cos(), (from + delta).sin()),
                ];
                Curve3::Nurbs(
                    opendatum_nurbs::NurbsCurve::new(
                        2,
                        vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
                        points,
                        Some(vec![1.0, w * c, c * c]),
                    )
                    .unwrap(),
                )
            };
            let curve;
            let (trace, t0, t1) = match kind {
                0 => {
                    let z = rng.unit();
                    let x = across(&mut rng, z);
                    let frame = Frame {
                        origin: rng.point(2.0) + shift,
                        x,
                        y: z.cross(x),
                        z,
                    };
                    let (a, b) = (rng.range(0.1, 2.0), rng.range(0.1, 2.0));
                    let map = |u: f64, v: f64| frame.x * (a * u) + frame.y * (b * v);
                    curve = arc(frame.origin, &map, rng.range(-PI, PI));
                    let chart = Chart::Plane(frame);
                    (trace(&chart, TraceOf::Curve(&curve)), 0.0, 1.0)
                }
                1 | 2 => {
                    let (chart, surface) = if kind == 1 {
                        cylinder(&mut rng)
                    } else {
                        cone(&mut rng)
                    };
                    let (chart, surface) = shifted(&chart, &surface, shift);
                    let (point, axis, x, y) = match chart {
                        Chart::Cylinder {
                            point, axis, x, y, ..
                        }
                        | Chart::Cone {
                            point, axis, x, y, ..
                        } => (point, axis, x, y),
                        _ => unreachable!(),
                    };
                    let from = -0.5 * delta;
                    curve = match surface {
                        // A section by a plane through the axis's point,
                        // leaning k.
                        Quadric::Cylinder { radius, .. } => {
                            let k = rng.range(-2.0, 2.0);
                            let map = |u: f64, v: f64| {
                                x * (radius * u) + y * (radius * v) + axis * (k * radius * u)
                            };
                            arc(point, &map, from)
                        }
                        // A circle about the axis, at height h.
                        Quadric::Cone { radius, slope, .. } => {
                            let h = rng.range(-0.5, 0.5) * radius / slope.max(0.1);
                            let r = radius + slope * h;
                            let map = |u: f64, v: f64| x * (r * u) + y * (r * v);
                            arc(point + axis * h, &map, from)
                        }
                        _ => unreachable!(),
                    };
                    (trace(&chart, TraceOf::Curve(&curve)), 0.0, 1.0)
                }
                3 => {
                    let z = rng.unit();
                    let x = across(&mut rng, z);
                    let frame = Frame {
                        origin: rng.point(2.0) + shift,
                        x,
                        y: z.cross(x),
                        z,
                    };
                    let p = rng.range(0.01, 1.0);
                    let cc = [-2.0, -0.3, 0.0, 0.5][rng.below(4)];
                    let track = Track::Conic {
                        vertex: frame.origin + frame.x * rng.range(-1.0, 1.0),
                        x: frame.x,
                        y: frame.y,
                        p,
                        c: cc,
                    };
                    let chart = Chart::Plane(frame);
                    let (t0, t1) = if cc <= 0.0 && rng.below(2) == 0 {
                        // A stretch 1 m to 10⁸ m from the vertex, where
                        // the allowance for rounding that grows with the
                        // reach is over the sixteenth from about 10⁶ m.
                        distant = true;
                        let t = if cc < 0.0 {
                            10f64.powf(rng.range(2.0, 8.0))
                        } else {
                            10f64.powf(rng.range(1.0, 4.5))
                        };
                        let (a, b) = (t, t * (1.0 + rng.range(0.01, 0.5)));
                        if rng.below(2) == 0 {
                            (a, b)
                        } else {
                            (-b, -a)
                        }
                    } else {
                        let reach = if cc > 0.0 { 0.9 / cc.sqrt() } else { 3.0 };
                        let t0 = rng.range(-reach, 0.0);
                        (t0, t0 + rng.range(0.1, reach))
                    };
                    (trace(&chart, TraceOf::Track(track)), t0, t1)
                }
                4 => {
                    // Up to 5·10⁵ m from the origin (meet refuses conics
                    // whose vertices lie 9·10⁵ m out, §14).
                    let shift = shift * (5e5 / far.max(5e5));
                    let (c, frame, t0, t1) = off_axis_hyperbola(&mut rng, shift);
                    curve = c;
                    let chart = Chart::Plane(frame);
                    (trace(&chart, TraceOf::Curve(&curve)), t0, t1)
                }
                _ => {
                    // A rational Bézier piece of degree 2 to 5 as large
                    // as its distance from the origin, its weights'
                    // ratio up to 10⁴, on a span of the parameter
                    // anywhere to ±10³ of length 10⁻³ to 10³.
                    distant = far > 0.0;
                    let z = rng.unit();
                    let x = across(&mut rng, z);
                    let frame = Frame {
                        origin: rng.point(2.0) + shift,
                        x,
                        y: z.cross(x),
                        z,
                    };
                    let p = 2 + rng.below(4);
                    let size = (1.0 + far) * rng.range(0.1, 2.0);
                    let rho = 10f64.powf(rng.range(0.0, 4.0));
                    let points: Vec<DVec3> = (0..=p)
                        .map(|_| {
                            frame.point2(DVec2::new(rng.range(-size, size), rng.range(-size, size)))
                        })
                        .collect();
                    let mut weights: Vec<f64> = (0..=p)
                        .map(|_| 10f64.powf(rng.range(0.0, rho.log10())))
                        .collect();
                    weights[rng.below(p + 1)] = rho;
                    let k0 = rng.range(-1e3, 1e3);
                    let len = 10f64.powf(rng.range(-3.0, 3.0));
                    let mut knots = vec![k0; p + 1];
                    knots.extend(vec![k0 + len; p + 1]);
                    curve = Curve3::Nurbs(
                        opendatum_nurbs::NurbsCurve::new(p, knots, points, Some(weights)).unwrap(),
                    );
                    let chart = Chart::Plane(frame);
                    (trace(&chart, TraceOf::Curve(&curve)), k0, k0 + len)
                }
            };
            let trace = trace.unwrap().unwrap();
            let size = 1.0 + far;
            let u = t0 + (t1 - t0) * rng.range(0.2, 0.8);
            let (mut lo, mut hi) = (t0, t1);
            let mut stopped = None;
            for k in 0..80 {
                let s = trace.sag(lo, hi);
                if stopped.is_none() && !halvable(s, lo, hi, fine) {
                    assert!(
                        !s.wide(fine),
                        "case {case}: {far} kind {kind}: halved to the parameter's last bit, \
                         sag {s:?}"
                    );
                    stopped = Some((k, s));
                }
                if hi - lo <= 1e-9 * (t1 - t0) {
                    let k = kind;
                    if s.floor > 0.0 {
                        most[k][0] = most[k][0].max(s.bound / s.floor);
                    }
                    most[k][1] = most[k][1].max(s.bound / (f64::EPSILON * size));
                }
                let m = 0.5 * (lo + hi);
                if m == lo || m == hi {
                    break;
                }
                if u < m {
                    hi = m;
                } else {
                    lo = m;
                }
            }
            let (k, s) = stopped.expect("halving stops");
            assert!(k <= 50, "case {case}: {far} kind {kind}: {k} halvings");
            let near = (far <= 1e2 || kind == 4 && far <= 1e4) && !distant;
            if near || s.settles && s.floor * (1.0 + SETTLED) < fine {
                assert!(
                    s.bound <= fine,
                    "case {case}: {far} kind {kind}: stopped at {s:?}"
                );
            }
        }
        eprintln!("bound over floor, over size ε, at tiny stretches: {most:?}");
        for [floor, _] in most {
            assert!(floor <= 3.0, "{most:?}");
        }
    }

    /// A conic track's sag settles on its floor, the allowance for rounding
    /// that grows with the stretch's reach from the vertex, so halving
    /// brings a chord within a sixteenth of the tolerance wherever the
    /// floor is under it by a thousandth ([`Sag::settles`]): parabolas
    /// 260 km to 890 km from their vertices, with floors from 0.3 to 1.01
    /// of the sixteenth; the last is held over it, to be refused. Stopped
    /// at four floors, as other bounds are, chords whose floors were over a
    /// quarter of the sixteenth (from 220 km) could stop over it, and from
    /// 0.58 of it (500 km) always did, refusing booleans that more halving
    /// resolves (found by review).
    #[test]
    fn conic_sags_settle_on_their_floors() {
        let fine = 1e-7 / 16.0;
        let frame = Frame {
            origin: DVec3::new(0.5, -0.25, 0.0),
            x: DVec3::X,
            y: DVec3::Y,
            z: DVec3::Z,
        };
        let track = Track::Conic {
            vertex: DVec3::ZERO,
            x: DVec3::X,
            y: DVec3::Y,
            p: 1.0,
            c: 0.0,
        };
        let tr = trace(&Chart::Plane(frame), TraceOf::Track(track))
            .unwrap()
            .unwrap();
        for ratio in [0.3, 0.6, 0.8, 0.95, 0.999, 1.01] {
            // The point at `reach` from the vertex: (t²/2, t), so t² = 2
            // (√(1 + reach²) − 1).
            let reach = ratio * fine / (CONIC_ROUNDING * f64::EPSILON);
            let t1 = (2.0 * ((1.0 + reach * reach).sqrt() - 1.0)).sqrt();
            let (mut t0, mut halvings) = (0.98 * t1, 0);
            let stop = loop {
                let s = tr.sag(t0, t1);
                if !halvable(s, t0, t1, fine) {
                    break s;
                }
                t0 = 0.5 * (t0 + t1);
                halvings += 1;
            };
            assert!(stop.settles, "{stop:?}");
            assert!(
                (stop.floor / (ratio * fine) - 1.0).abs() < 1e-6,
                "{ratio}: {stop:?}"
            );
            assert!(!stop.wide(fine) && halvings <= 50, "{ratio}: {halvings}");
            if ratio < 1.0 {
                assert!(stop.bound <= fine, "{ratio}: stopped at {stop:?}");
            } else {
                assert!(
                    stop.bound <= stop.floor * (1.0 + SETTLED),
                    "{ratio}: {stop:?}"
                );
            }
        }
    }

    /// The hyperbola of a cone (radius 1 cm to 10 cm at its point, slope
    /// 0.1 to 3) cut by a plane parallel to its axis a millimetre off it,
    /// at `shift` (at most 5·10⁵ m from the origin), as a NURBS curve (its weights' ratio in the hundreds),
    /// with the plane's frame and the curve's parameter range.
    fn off_axis_hyperbola(rng: &mut Rng, shift: DVec3) -> (Curve3, Frame, f64, f64) {
        loop {
            let axis = rng.unit();
            let radius = 10f64.powf(rng.range(-2.0, -1.0));
            let cone = Quadric::Cone {
                point: shift,
                axis,
                radius,
                slope: 10f64.powf(rng.range(-1.0, 0.5)),
            };
            let x = across(rng, axis);
            let point = shift + x * 1e-3;
            let plane = Quadric::Plane { point, normal: x };
            let r = DVec3::splat(4.0 * radius);
            let m = meet::meet(
                &plane,
                &cone,
                &meet::Region::of_box(shift - r, shift + r),
                1e-7,
            );
            let meet::Meet::Curves(tracks) = m else {
                continue;
            };
            for track in tracks.iter().filter(|t| matches!(t, Track::Conic { .. })) {
                let (lo, hi) = track.ends();
                let reach = 10f64.powf(rng.range(-1.0, 4.0));
                let (lo, hi) = (lo.max(-reach), hi.min(reach));
                if lo.partial_cmp(&hi) != Some(std::cmp::Ordering::Less) {
                    continue;
                }
                if let Ok((curve @ Curve3::Nurbs(_), t0, t1)) = track.edge(lo, hi) {
                    let frame = Frame::from_normal(point, x, x.any_orthonormal_vector());
                    return (curve, frame, t0, t1);
                }
            }
        }
    }

    /// `chart` and `surface` (a cylinder's or a cone's) moved by `shift`.
    fn shifted(chart: &Chart, surface: &Quadric, shift: DVec3) -> (Chart, Quadric) {
        let chart = match *chart {
            Chart::Cylinder {
                point,
                axis,
                x,
                y,
                radius,
            } => Chart::Cylinder {
                point: point + shift,
                axis,
                x,
                y,
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
                point: point + shift,
                axis,
                x,
                y,
                radius,
                sin,
                cos,
            },
            _ => unreachable!(),
        };
        let surface = match *surface {
            Quadric::Cylinder {
                point,
                axis,
                radius,
            } => Quadric::Cylinder {
                point: point + shift,
                axis,
                radius,
            },
            Quadric::Cone {
                point,
                axis,
                radius,
                slope,
            } => Quadric::Cone {
                point: point + shift,
                axis,
                radius,
                slope,
            },
            _ => unreachable!(),
        };
        (chart, surface)
    }

    /// Straight NURBS segments off a cylinder or a cone, unevenly
    /// parameterised: in space they do not bend, so their courses in the
    /// chart bend only through the chart's map off the surface, which the
    /// bands' terms in the distance from the surface bound (design note
    /// 0012 §16.2).
    #[test]
    fn bands_hold_courses_of_segments_off_the_surface() {
        let mut rng = Rng(0xC4A2_7006);
        let mut checked = 0;
        for case in 0..2000 {
            let kind = case % 2;
            let (chart, surface) = if kind == 0 {
                cylinder(&mut rng)
            } else {
                cone(&mut rng)
            };
            let p = if kind == 0 {
                cylinder_point(&mut rng, &chart)
            } else {
                cone_point(&mut rng, &chart)
            };
            let chart = facing(&chart, p);
            let radius = match surface {
                Quadric::Cylinder { radius, .. } => radius,
                Quadric::Cone {
                    point,
                    axis,
                    radius,
                    slope,
                } => radius + slope * (p - point).dot(axis),
                _ => unreachable!(),
            };
            // Off the surface along its normal by up to a third of the
            // radius there, along a random direction.
            let n = surface.gradient(p).normalize();
            let start = p + n * (radius * rng.range(-0.3, 0.3));
            let end = start + rng.unit() * (radius * rng.range(0.01, 1.0));
            let ok = (0..=64).all(|k| {
                let q = start + (end - start) * (k as f64 / 64.0);
                angle_in(&chart, q).abs() < 0.45 * PI && surface.distance(q).abs() < 0.5 * radius
            });
            if !ok {
                continue;
            }
            let w = rng.range(0.3, 3.0);
            let curve = Curve3::Nurbs(
                opendatum_nurbs::NurbsCurve::new(
                    2,
                    vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
                    vec![start, (start + end) * 0.5, end],
                    Some(vec![1.0, w, 1.0]),
                )
                .unwrap(),
            );
            let trace = trace(&chart, TraceOf::Curve(&curve)).unwrap().unwrap();
            let Ok(ts) = trace.division(0.0, 1.0) else {
                continue;
            };
            for win in ts.windows(2) {
                let (a, b) = (win[0], win[1]);
                for (x, y) in [(a, b), (a, 0.5 * (a + b)), (b, 0.25 * a + 0.75 * b)] {
                    let band = trace.sag(x, y).bound + trace.slack;
                    let off = strays(&trace, x, y);
                    assert!(off <= band, "case {case}: strays {off} beyond band {band}");
                }
            }
            checked += 1;
        }
        assert!(checked > 1000, "{checked}");
    }

    /// A NURBS edge on a sphere's face is refused, naming the sphere.
    #[test]
    fn rational_edges_on_spheres_are_refused() {
        let mut rng = Rng(0xC4A2_7005);
        let (chart, _) = sphere(&mut rng);
        let curve = Curve3::Nurbs(
            opendatum_nurbs::NurbsCurve::new(
                2,
                vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
                vec![DVec3::X, DVec3::new(1.0, 1.0, 0.0), DVec3::Y],
                Some(vec![1.0, std::f64::consts::FRAC_1_SQRT_2, 1.0]),
            )
            .unwrap(),
        );
        match trace(&chart, TraceOf::Curve(&curve)) {
            Err(KernelError::Unsupported(m)) => assert!(m.contains("sphere"), "{m}"),
            Err(e) => panic!("{e:?}"),
            Ok(_) => panic!("traced"),
        }
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
                    let band = trace.sag(u, v).bound + trace.slack;
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

#[cfg(test)]
#[path = "chart_rounding.rs"]
mod rounding;
