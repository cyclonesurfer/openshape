//! Charts of faces (design note 0012 §10.4): plane coordinates on a face's
//! surface in which lengths are lengths along the surface, the traces of
//! curves in them with bounds on how far each strays from its chords, and
//! the outlines of faces, for deciding whether a point of the surface lies
//! in a face.
//!
//! A plane's chart is its frame; a cylinder's is the cylinder unrolled
//! about the middle of the face (which spans at most half a turn, ADR
//! 0008), so the face never reaches the cut of the unrolling. Both are
//! isometries of the surface onto the plane, and both are oriented so that
//! a face's outer loop, counter-clockwise about its outward normal, is
//! counter-clockwise in the chart.
//!
//! A curve's course in a chart is followed by chords. The distance from a
//! course to the chord between two of its points at parameters `t0` and
//! `t1` is at most `M (t1 − t0)² / 8`, `M` a bound on the second derivative
//! of the course (the error of linear interpolation); the second derivative
//! of a curve's image in an isometric chart of a developable surface is the
//! part of its acceleration along the surface, so `M` is at most the size
//! of the curve's acceleration, and less for the curves that unroll to
//! lines (rulings, and circles square to a cylinder's axis, which unroll
//! straight). A NURBS curve in a plane is bounded instead by its Bézier
//! pieces' control points, whose convex hull holds the piece.

use std::f64::consts::{FRAC_PI_8, PI};

use glam::{DVec2, DVec3};

use super::meet::Track;
use super::planar::{point_segment, Bend};
use super::quadric::{self, Quadric};
use super::shape::{Face, Shape};
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
        }
    }

    /// A lower bound on the distance in space between two points of the
    /// surface `d` apart in the chart (both within half a turn of the
    /// chart's middle): `d` itself on a plane; on a cylinder the chord of
    /// the arc, the shortest way two points `d` apart along the surface can
    /// be.
    pub fn space(&self, d: f64) -> f64 {
        match *self {
            Chart::Plane(_) => d,
            Chart::Cylinder { radius, .. } => {
                if d >= PI * radius {
                    2.0 * radius
                } else {
                    2.0 * radius * (d / (2.0 * radius)).sin()
                }
            }
        }
    }
}

/// The chart of face `f` of `shape`: its plane's frame, or its cylinder
/// unrolled about the direction of the middle of its boundary.
pub(super) fn chart_of(shape: &Shape, f: &Face) -> Result<Chart, KernelError> {
    match shape.surfaces[f.surface] {
        Quadric::Plane { .. } => Ok(Chart::Plane(f.frame)),
        Quadric::Cylinder {
            point,
            axis,
            radius,
        } => {
            // The mean of the directions from the axis to points along the
            // boundary: the face spans at most half a turn, and points
            // inside its edges weigh its middle.
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
                    "booleans with a face going around its cylinder (design note 0012 §10)".into(),
                ));
            }
            let x = x.normalize();
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
        _ => Err(KernelError::Unsupported(
            "booleans of cones and spheres (design note 0012 §10.7, item 3)".into(),
        )),
    }
}

/// How a course is followed.
#[derive(Clone, Copy)]
enum Path<'a> {
    /// A track (a circle or an ellipse, by angle), with a bound on the
    /// second derivative of its image per unit angle squared.
    Track(Track, f64),
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
    /// eighth of a half turn each along a track, two per Bézier piece of a
    /// NURBS curve.
    pub fn division(&self, t0: f64, t1: f64) -> Vec<f64> {
        let mut out = vec![t0];
        match self.path {
            Path::Track(_, m) => {
                let n = if m == 0.0 {
                    1
                } else {
                    ((t1 - t0).abs() / FRAC_PI_8).ceil().max(1.0) as usize
                };
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
        out
    }
}

impl Bend for Trace<'_> {
    fn at(&self, t: f64) -> DVec2 {
        self.chart.map(self.point(t))
    }

    fn band(&self, t0: f64, t1: f64) -> f64 {
        let b = match self.path {
            Path::Track(_, m) => m * (t1 - t0) * (t1 - t0) / 8.0,
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
/// circle square to the axis on a cylinder).
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
    let (m, size) = match (chart, track) {
        (_, Track::Line { .. }) => return Ok(None),
        (Chart::Plane(_), Track::Circle { center, radius, .. }) => {
            (radius, center.abs().max_element() + radius)
        }
        (Chart::Plane(_), Track::Ellipse { center, a, b, .. }) => {
            (a.max(b), center.abs().max_element() + a.max(b))
        }
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
            (radius, size)
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
            (amp, center.abs().max_element() + a.max(b))
        }
    };
    Ok(Some(Trace {
        chart: *chart,
        path: Path::Track(track, m),
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
                        let ts = c.division(t0, t1);
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
                            chords.push(OChord {
                                p: at[k],
                                q: at[k + 1],
                                band: c.band(ts[k], ts[k + 1]) + moved,
                                sag: c.band(ts[k], ts[k + 1]),
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
    /// `Some(false)` outside, `None` within `margin` of its boundary.
    /// Chords near `q` are halved, for this question only, until their
    /// sags are within `fine` or their bands clear of `q`; the answer is the parity of
    /// a ray from `q` across the chords, which is the parity across the
    /// boundary for a point outside every band.
    pub fn locate(&self, q: DVec2, margin: f64, fine: f64) -> Option<bool> {
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
                                stack.push(OChord {
                                    p,
                                    q: s,
                                    band: trace.band(a, b) + moved,
                                    sag: trace.band(a, b),
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
