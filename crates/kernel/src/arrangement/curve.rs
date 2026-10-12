//! Validated internal form of the input curves: exact evaluation in the
//! module's curve parameters, rational Bézier decomposition, and closest
//! points.

use std::f64::consts::{FRAC_PI_2, TAU};

use glam::{DVec2, DVec3};
use openshape_nurbs::NurbsCurve;

use super::bezier::{self, Bez};
use crate::{KernelError, ProfileCurve};

/// Highest supported spline degree.
const MAX_DEGREE: usize = 24;
/// Most control points accepted per spline.
const MAX_POINTS: usize = 10_000;
/// Coordinates must stay within this multiple of the linear tolerance, so
/// that the tolerance stays well above the floating-point resolution.
pub(super) const RANGE_IN_TOLERANCES: f64 = 1e13;

/// One polynomial piece (knot span) of a spline as a rational Bézier curve.
#[derive(Clone, Debug)]
pub(super) struct Span {
    /// Knot span start.
    pub t0: f64,
    /// Knot span end.
    pub t1: f64,
    /// The span's Bézier form on `u = (t − t0) / (t1 − t0)`.
    pub bez: Bez,
    /// Bounds of the control points.
    pub lo: DVec2,
    /// Bounds of the control points.
    pub hi: DVec2,
}

/// Exact geometry of a curve.
#[derive(Clone, Debug)]
pub(super) enum Geom {
    /// `a + t (b − a)`, `t ∈ [0, 1]`.
    Line { a: DVec2, b: DVec2 },
    /// `c + r (cos φ, sin φ)` with `φ = a0 + sgn · t`.
    Arc { c: DVec2, r: f64, a0: f64, sgn: f64 },
    /// B-spline, decomposed into its knot spans.
    Spline { spans: Vec<Span> },
}

/// A validated input curve.
#[derive(Clone, Debug)]
pub(super) struct Crv {
    /// Geometry.
    pub geom: Geom,
    /// Parameter domain start.
    pub lo: f64,
    /// Parameter domain end.
    pub hi: f64,
    /// True for a full circle, whose parameter wraps.
    pub periodic: bool,
    /// Bounding box.
    pub bmin: DVec2,
    /// Bounding box.
    pub bmax: DVec2,
}

/// A rational Bézier piece of a curve together with the curve parameters of
/// its ends.
#[derive(Clone, Debug)]
pub(super) struct BPiece {
    /// Geometry.
    pub bez: Bez,
    /// Curve parameter at `u = 0`.
    pub t0: f64,
    /// Curve parameter at `u = 1`.
    pub t1: f64,
    /// True if the curve parameter is affine in `u`; otherwise (arcs) it is
    /// recovered from the point.
    pub linear: bool,
}

fn bad(i: usize, m: impl std::fmt::Display) -> KernelError {
    KernelError::InvalidInput(format!("curve {i}: {m}"))
}

/// Control points of the Bézier form of knot span `k` (`u_k < u_{k+1}`) by
/// blossoming: the `i`-th is the blossom `f(u_k^{p−i}, u_{k+1}^i)`, computed
/// with de Boor's recurrence (Piegl & Tiller §5.6; Ramshaw, "Blossoms are
/// polar forms", 1989).
fn span_bezier(u: &[f64], pw: &[DVec3], p: usize, k: usize) -> Vec<DVec3> {
    let (a, b) = (u[k], u[k + 1]);
    (0..=p)
        .map(|i| {
            let mut d: Vec<DVec3> = pw[k - p..=k].to_vec();
            for r in 1..=p {
                let x = if r <= p - i { a } else { b };
                for j in (r..=p).rev() {
                    let lo = u[k - p + j];
                    let hi = u[k + j - r + 1];
                    let alpha = if hi > lo { (x - lo) / (hi - lo) } else { 0.0 };
                    d[j] = d[j - 1] * (1.0 - alpha) + d[j] * alpha;
                }
            }
            d[p]
        })
        .collect()
}

/// True if the control polygon's edge directions fit in an open half-plane,
/// so the curve is strictly monotone along some direction and cannot
/// intersect itself (variation diminishing property, Farin §5.3).
pub(super) fn is_monotone(hull: &[DVec2]) -> bool {
    let mut ang: Vec<f64> = hull
        .windows(2)
        .map(|w| w[1] - w[0])
        .filter(|e| e.length_squared() > 0.0)
        .map(|e| e.y.atan2(e.x))
        .collect();
    if ang.len() <= 1 {
        return true;
    }
    ang.sort_by(f64::total_cmp);
    let mut gap = ang[0] + TAU - ang[ang.len() - 1];
    for w in ang.windows(2) {
        gap = gap.max(w[1] - w[0]);
    }
    gap > std::f64::consts::PI + 1e-9
}

impl Crv {
    /// Validates input curve `i` and builds its internal form.
    pub fn new(i: usize, c: &ProfileCurve, tol: f64) -> Result<Crv, KernelError> {
        let limit = tol * RANGE_IN_TOLERANCES;
        let in_range = |p: DVec2| p.is_finite() && p.abs().max_element() <= limit;
        match c {
            ProfileCurve::Line { start, end } => {
                if !(in_range(*start) && in_range(*end)) {
                    return Err(bad(i, "line coordinates are not finite or out of range"));
                }
                if (*end - *start).length() <= tol {
                    return Err(bad(i, "line is shorter than the model tolerance"));
                }
                Ok(Crv {
                    geom: Geom::Line { a: *start, b: *end },
                    lo: 0.0,
                    hi: 1.0,
                    periodic: false,
                    bmin: start.min(*end),
                    bmax: start.max(*end),
                })
            }
            ProfileCurve::Arc {
                center,
                radius,
                start_angle,
                sweep_angle,
            } => {
                if !(center.is_finite()
                    && radius.is_finite()
                    && start_angle.is_finite()
                    && sweep_angle.is_finite())
                {
                    return Err(bad(i, "arc has non-finite parameters"));
                }
                if !in_range(*center) || radius.abs() + center.abs().max_element() > limit {
                    return Err(bad(i, "arc is out of the coordinate range"));
                }
                if *radius <= tol {
                    return Err(bad(i, "arc radius must exceed the model tolerance"));
                }
                if sweep_angle.abs() > TAU + 1e-12 {
                    return Err(bad(i, "arc sweep exceeds a full turn"));
                }
                if sweep_angle.abs() * radius <= tol {
                    return Err(bad(i, "arc is shorter than the model tolerance"));
                }
                let len = sweep_angle.abs().min(TAU);
                let mut crv = Crv {
                    geom: Geom::Arc {
                        c: *center,
                        r: *radius,
                        a0: start_angle.rem_euclid(TAU),
                        sgn: sweep_angle.signum(),
                    },
                    lo: 0.0,
                    hi: len,
                    periodic: false,
                    bmin: DVec2::ZERO,
                    bmax: DVec2::ZERO,
                };
                crv.arc_bounds();
                Ok(crv)
            }
            ProfileCurve::Circle { center, radius } => {
                if !(center.is_finite() && radius.is_finite()) {
                    return Err(bad(i, "circle has non-finite parameters"));
                }
                if !in_range(*center) || radius.abs() + center.abs().max_element() > limit {
                    return Err(bad(i, "circle is out of the coordinate range"));
                }
                if *radius <= tol {
                    return Err(bad(i, "circle radius must exceed the model tolerance"));
                }
                Ok(Crv {
                    geom: Geom::Arc {
                        c: *center,
                        r: *radius,
                        a0: 0.0,
                        sgn: 1.0,
                    },
                    lo: 0.0,
                    hi: TAU,
                    periodic: true,
                    bmin: *center - DVec2::splat(*radius),
                    bmax: *center + DVec2::splat(*radius),
                })
            }
            ProfileCurve::Spline {
                degree,
                knots,
                points,
                weights,
            } => Self::spline(i, *degree, knots, points, weights.as_deref(), tol),
        }
    }

    fn spline(
        i: usize,
        degree: usize,
        knots: &[f64],
        points: &[DVec2],
        weights: Option<&[f64]>,
        tol: f64,
    ) -> Result<Crv, KernelError> {
        let limit = tol * RANGE_IN_TOLERANCES;
        if degree > MAX_DEGREE {
            return Err(bad(i, format!("spline degree above {MAX_DEGREE}")));
        }
        if points.len() > MAX_POINTS {
            return Err(bad(i, format!("spline has more than {MAX_POINTS} points")));
        }
        if points
            .iter()
            .any(|p| !(p.is_finite() && p.abs().max_element() <= limit))
        {
            return Err(bad(i, "spline points are not finite or out of range"));
        }
        if let Some(w) = weights {
            let (mut wmin, mut wmax) = (f64::INFINITY, 0.0f64);
            for &x in w {
                if !(x.is_finite() && x > 0.0) {
                    return Err(bad(i, "spline weights must be positive"));
                }
                wmin = wmin.min(x);
                wmax = wmax.max(x);
            }
            if wmax > wmin * 1e12 {
                return Err(bad(
                    i,
                    "spline weights span more than 12 orders of magnitude",
                ));
            }
        }
        let p3: Vec<DVec3> = points.iter().map(|p| p.extend(0.0)).collect();
        let n = NurbsCurve::new(degree, knots.to_vec(), p3, weights.map(|w| w.to_vec()))
            .map_err(|e| bad(i, format!("spline: {e}")))?;
        let (a, b) = n.domain();
        if (a.abs() + b.abs()) > 1e6 * (b - a) {
            return Err(bad(i, "spline knot values are too large for its domain"));
        }
        let u = n.knots().as_slice();
        let pw: Vec<DVec3> = n
            .homogeneous_points()
            .iter()
            .map(|h| DVec3::new(h.x, h.y, h.w))
            .collect();
        let mut spans = Vec::new();
        for k in degree..pw.len() {
            let (t0, t1) = (u[k], u[k + 1]);
            if t1 <= t0 {
                continue;
            }
            if (b - a) > 1e14 * (t1 - t0) {
                return Err(bad(i, "spline knot spacing is below numerical resolution"));
            }
            let bez = Bez {
                pts: span_bezier(u, &pw, degree, k),
            };
            if bez
                .pts
                .iter()
                .any(|h| !(h.is_finite() && h.z > 0.0 && bezier::project(*h).is_finite()))
            {
                return Err(bad(i, "spline evaluation is not finite"));
            }
            let (lo, hi) = bezier::bounds(&bez.hull());
            spans.push(Span {
                t0,
                t1,
                bez,
                lo,
                hi,
            });
        }
        if spans.is_empty() {
            return Err(bad(i, "spline has an empty domain"));
        }
        let mut bmin = DVec2::splat(f64::INFINITY);
        let mut bmax = DVec2::splat(f64::NEG_INFINITY);
        for s in &spans {
            bmin = bmin.min(s.lo);
            bmax = bmax.max(s.hi);
        }
        let crv = Crv {
            geom: Geom::Spline { spans },
            lo: a,
            hi: b,
            periodic: false,
            bmin,
            bmax,
        };
        // Reject curves that never leave a tolerance-sized neighbourhood.
        let p0 = crv.point(a);
        let mut far = false;
        if let Geom::Spline { spans } = &crv.geom {
            'outer: for s in spans {
                for j in 0..=8 {
                    if (s.bez.point(j as f64 / 8.0) - p0).length() > tol {
                        far = true;
                        break 'outer;
                    }
                }
            }
        }
        if !far {
            return Err(bad(i, "spline is shorter than the model tolerance"));
        }
        Ok(crv)
    }

    fn arc_bounds(&mut self) {
        let (c, r) = match self.geom {
            Geom::Arc { c, r, .. } => (c, r),
            _ => return,
        };
        let mut pts = vec![self.point(self.lo), self.point(self.hi)];
        for k in 0..4 {
            let phi = k as f64 * FRAC_PI_2;
            let t = self.arc_param(c + r * DVec2::from_angle(phi));
            if t <= self.hi {
                pts.push(c + r * DVec2::from_angle(phi));
            }
        }
        let (lo, hi) = bezier::bounds(&pts);
        self.bmin = lo;
        self.bmax = hi;
    }

    /// Size of the coordinates the curve is evaluated from (its end
    /// points, its centre plus radius, or its control points' bounds): its
    /// points carry rounding errors of about this times the machine
    /// epsilon.
    pub fn magnitude(&self) -> f64 {
        match &self.geom {
            Geom::Line { a, b } => a.abs().max(b.abs()).max_element(),
            Geom::Arc { c, r, .. } => c.abs().max_element() + r,
            Geom::Spline { spans } => spans
                .iter()
                .map(|s| s.lo.abs().max(s.hi.abs()).max_element())
                .fold(0.0, f64::max),
        }
    }

    /// Spline spans, if this is a spline.
    pub fn spans(&self) -> &[Span] {
        match &self.geom {
            Geom::Spline { spans } => spans,
            _ => &[],
        }
    }

    /// Index of the span containing `t`: the one to the right of a knot when
    /// `right`, else the one to the left.
    fn span_at(spans: &[Span], t: f64, right: bool) -> usize {
        let k = if right {
            spans.partition_point(|s| s.t1 <= t)
        } else {
            spans.partition_point(|s| s.t1 < t)
        };
        k.min(spans.len() - 1)
    }

    /// Point at parameter `t`.
    pub fn point(&self, t: f64) -> DVec2 {
        match &self.geom {
            Geom::Line { a, b } => {
                if t == 1.0 {
                    *b
                } else {
                    *a + (*b - *a) * t
                }
            }
            Geom::Arc { c, r, a0, sgn } => *c + *r * DVec2::from_angle(a0 + sgn * t),
            Geom::Spline { spans } => {
                let s = &spans[Self::span_at(spans, t, true)];
                let u = ((t - s.t0) / (s.t1 - s.t0)).clamp(0.0, 1.0);
                if u == 0.0 {
                    s.bez.start()
                } else if u == 1.0 {
                    s.bez.end()
                } else {
                    s.bez.point(u)
                }
            }
        }
    }

    /// Point and first two derivatives with respect to `t`; at a spline knot
    /// the one-sided derivatives of the span on the given side.
    pub fn eval(&self, t: f64, right: bool) -> [DVec2; 3] {
        match &self.geom {
            Geom::Line { a, b } => [self.point(t), *b - *a, DVec2::ZERO],
            Geom::Arc { c, r, a0, sgn } => {
                let e = DVec2::from_angle(a0 + sgn * t);
                [*c + *r * e, e.perp() * (*r * *sgn), -e * *r]
            }
            Geom::Spline { spans } => {
                let s = &spans[Self::span_at(spans, t, right)];
                let h = s.t1 - s.t0;
                let u = ((t - s.t0) / h).clamp(0.0, 1.0);
                let [p, d1, d2] = s.bez.derivs(u);
                [p, d1 / h, d2 / (h * h)]
            }
        }
    }

    /// Wraps (periodic) or clamps a parameter into the domain.
    pub fn clamp(&self, t: f64) -> f64 {
        if self.periodic {
            t.rem_euclid(TAU)
        } else {
            t.clamp(self.lo, self.hi)
        }
    }

    /// Parameter (in `[0, 2π)`) of the point of an arc's circle nearest to `p`.
    pub fn arc_param(&self, p: DVec2) -> f64 {
        match self.geom {
            Geom::Arc { c, a0, sgn, .. } => {
                let v = p - c;
                (sgn * (v.y.atan2(v.x) - a0)).rem_euclid(TAU)
            }
            _ => 0.0,
        }
    }

    /// Like [`Crv::arc_param`], unwrapped to the representative nearest `near`.
    pub fn arc_param_near(&self, p: DVec2, near: f64) -> f64 {
        let t = self.arc_param(p);
        t + TAU * ((near - t) / TAU).round()
    }

    /// Closest point of the piece `[t0, t1]` to `p`: `(parameter, distance)`.
    pub fn closest(&self, p: DVec2, t0: f64, t1: f64) -> (f64, f64) {
        let dist = |t: f64| (self.point(t) - p).length();
        let mut best = (t0, dist(t0));
        let consider = |t: f64, best: &mut (f64, f64)| {
            let d = dist(t);
            if d < best.1 {
                *best = (t, d);
            }
        };
        consider(t1, &mut best);
        match &self.geom {
            Geom::Line { a, b } => {
                let d = *b - *a;
                let t = ((p - *a).dot(d) / d.length_squared()).clamp(t0, t1);
                consider(t, &mut best);
            }
            Geom::Arc { c, .. } => {
                if (p - *c).length_squared() > 0.0 {
                    let tau = self.arc_param(p);
                    let t = t0 + (tau - t0).rem_euclid(TAU);
                    if t <= t1 {
                        consider(t, &mut best);
                    }
                }
            }
            Geom::Spline { spans } => {
                for s in spans {
                    if s.t1 <= t0 || s.t0 >= t1 {
                        continue;
                    }
                    let q = p.clamp(s.lo, s.hi);
                    if (q - p).length() >= best.1 {
                        continue;
                    }
                    let h = s.t1 - s.t0;
                    let ua = ((t0.max(s.t0) - s.t0) / h).clamp(0.0, 1.0);
                    let ub = ((t1.min(s.t1) - s.t0) / h).clamp(0.0, 1.0);
                    if let Some((u, d)) = closest_on_bez(&s.bez, p, ua, ub) {
                        if d < best.1 {
                            best = ((s.t0 + u * h).clamp(t0, t1), d);
                        }
                    }
                }
            }
        }
        best
    }

    /// Rational Bézier pieces covering `[t0, t1]` (`t0 < t1`), each turning
    /// by at most a quarter circle for arcs.
    pub fn pieces(&self, t0: f64, t1: f64) -> Vec<BPiece> {
        match &self.geom {
            Geom::Line { .. } => vec![BPiece {
                bez: Bez::line(self.point(t0), self.point(t1)),
                t0,
                t1,
                linear: true,
            }],
            Geom::Arc { c, r, a0, sgn } => {
                let n = ((t1 - t0) / FRAC_PI_2 - 1e-9).ceil().max(1.0) as usize;
                (0..n)
                    .map(|k| {
                        let ta = t0 + (t1 - t0) * k as f64 / n as f64;
                        let tb = if k + 1 == n {
                            t1
                        } else {
                            t0 + (t1 - t0) * (k + 1) as f64 / n as f64
                        };
                        BPiece {
                            bez: Bez::arc(*c, *r, a0 + sgn * ta, a0 + sgn * tb),
                            t0: ta,
                            t1: tb,
                            linear: false,
                        }
                    })
                    .collect()
            }
            Geom::Spline { spans } => spans
                .iter()
                .filter(|s| s.t1 > t0 && s.t0 < t1)
                .map(|s| {
                    let a = t0.max(s.t0);
                    let b = t1.min(s.t1);
                    let h = s.t1 - s.t0;
                    let ua = ((a - s.t0) / h).clamp(0.0, 1.0);
                    let ub = ((b - s.t0) / h).clamp(0.0, 1.0);
                    let bez = if ua <= 0.0 && ub >= 1.0 {
                        s.bez.clone()
                    } else {
                        s.bez.sub(ua, ub)
                    };
                    BPiece {
                        bez,
                        t0: a,
                        t1: b,
                        linear: true,
                    }
                })
                .filter(|p| p.t1 > p.t0)
                .collect(),
        }
    }

    /// Pieces of the whole curve that are each monotone along some direction
    /// (so free of self-intersections), or smaller than `tol`.
    pub fn monotone_pieces(&self, tol: f64) -> Vec<BPiece> {
        let mut out = Vec::new();
        for p in self.pieces(self.lo, self.hi) {
            self.split_monotone(p, tol, 0, &mut out);
        }
        out
    }

    fn split_monotone(&self, p: BPiece, tol: f64, depth: u32, out: &mut Vec<BPiece>) {
        let hull = p.bez.hull();
        let (lo, hi) = bezier::bounds(&hull);
        if depth >= 24 || (hi - lo).length() <= tol || is_monotone(&hull) {
            out.push(p);
            return;
        }
        let (l, r) = p.bez.split(0.5);
        let tm = p.param_at(self, 0.5);
        self.split_monotone(
            BPiece {
                bez: l,
                t0: p.t0,
                t1: tm,
                linear: p.linear,
            },
            tol,
            depth + 1,
            out,
        );
        self.split_monotone(
            BPiece {
                bez: r,
                t0: tm,
                t1: p.t1,
                linear: p.linear,
            },
            tol,
            depth + 1,
            out,
        );
    }

    /// Interior knots of a spline where the curve turns by more than a right
    /// angle (sharp corners and reversals).
    pub fn corners(&self) -> Vec<f64> {
        let spans = self.spans();
        let first_dir = |h: &[DVec2]| {
            h.windows(2)
                .map(|w| w[1] - w[0])
                .find(|e| e.length_squared() > 0.0)
        };
        let mut out = Vec::new();
        for w in spans.windows(2) {
            if w[0].t1 != w[1].t0 {
                continue;
            }
            let mut left = w[0].bez.hull();
            left.reverse();
            let right = w[1].bez.hull();
            if let (Some(l), Some(r)) = (first_dir(&left), first_dir(&right)) {
                // `l` points backwards along the curve.
                if (-l).dot(r) < 0.0 {
                    out.push(w[1].t0);
                }
            }
        }
        out
    }
}

impl BPiece {
    /// Curve parameter at local parameter `u`.
    pub fn param_at(&self, crv: &Crv, u: f64) -> f64 {
        if self.linear {
            if u <= 0.0 {
                self.t0
            } else if u >= 1.0 {
                self.t1
            } else {
                self.t0 + (self.t1 - self.t0) * u
            }
        } else if u <= 0.0 {
            self.t0
        } else if u >= 1.0 {
            self.t1
        } else {
            crv.arc_param_near(self.bez.point(u), 0.5 * (self.t0 + self.t1))
        }
    }
}

/// Closest point of a Bézier curve's piece `[ua, ub]` to `p`, by sampling and
/// safeguarded Newton iteration on `(C(u) − p) · C'(u) = 0` (Piegl & Tiller
/// §6.1).
fn closest_on_bez(bez: &Bez, p: DVec2, ua: f64, ub: f64) -> Option<(f64, f64)> {
    if ub < ua {
        return None;
    }
    let m = 4 * (bez.degree() + 1) + 4;
    let us: Vec<f64> = (0..=m)
        .map(|k| ua + (ub - ua) * k as f64 / m as f64)
        .collect();
    let ds: Vec<f64> = us.iter().map(|&u| (bez.point(u) - p).length()).collect();
    let mut best = (us[0], ds[0]);
    for k in 0..us.len() {
        if ds[k] < best.1 {
            best = (us[k], ds[k]);
        }
    }
    let mut cands: Vec<usize> = (0..us.len())
        .filter(|&k| (k == 0 || ds[k] <= ds[k - 1]) && (k + 1 == us.len() || ds[k] <= ds[k + 1]))
        .collect();
    cands.sort_by(|&a, &b| ds[a].total_cmp(&ds[b]));
    cands.truncate(4);
    for k in cands {
        let mut u = us[k];
        let mut d = ds[k];
        for _ in 0..40 {
            let [c, c1, c2] = bez.derivs(u);
            let r = c - p;
            let g = r.dot(c1);
            let gp = c1.length_squared() + r.dot(c2);
            let step = if gp > 0.0 {
                -g / gp
            } else if c1.length_squared() > 0.0 {
                -g / c1.length_squared()
            } else {
                break;
            };
            if !step.is_finite() {
                break;
            }
            let mut alpha = 1.0;
            let mut moved = false;
            for _ in 0..30 {
                let un = (u + alpha * step).clamp(ua, ub);
                if un == u {
                    break;
                }
                let dn = (bez.point(un) - p).length();
                if dn <= d {
                    moved = (un - u).abs() > 1e-16 * (1.0 + u.abs());
                    u = un;
                    d = dn;
                    break;
                }
                alpha *= 0.5;
            }
            if !moved {
                break;
            }
        }
        if d < best.1 {
            best = (u, d);
        }
    }
    Some(best)
}
