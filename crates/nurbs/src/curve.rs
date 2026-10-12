//! Non-uniform rational B-spline curves.

use std::f64::consts::{FRAC_PI_2, PI, TAU};

use glam::{DMat4, DVec3, DVec4};
use serde::{Deserialize, Serialize};

use crate::knots::KnotVector;
use crate::util::{
    aabb_distance_squared, binomial, bounding_box, dehomogenize, homogenize, integrate,
    local_minima, solve_banded, transform_homogeneous, unit,
};
use crate::NurbsError;

/// Maximum number of Newton iterations in point inversion.
const MAX_NEWTON_ITERATIONS: usize = 64;
/// Maximum number of sampled local minima per knot span refined by point
/// inversion.
const MAX_CANDIDATES_PER_SPAN: usize = 4;

/// A NURBS curve of degree `p >= 1`.
///
/// Control points are stored in homogeneous form `(w x, w y, w z, w)` with
/// `w > 0`. The knot vector may be clamped or unclamped; the curve is defined
/// on its domain `[u_p, u_{n+1}]` and parameters outside it are clamped.
///
/// Invariants (checked by every constructor and by deserialization):
/// `knots.len() == points.len() + degree + 1`, the domain is non-empty, interior
/// knots have multiplicity at most `degree`, all weights are positive and all
/// coordinates are finite.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(try_from = "CurveData", into = "CurveData")]
pub struct NurbsCurve {
    degree: usize,
    knots: KnotVector,
    points: Vec<DVec4>,
}

/// Unvalidated serialized form of a [`NurbsCurve`].
#[derive(Serialize, Deserialize)]
struct CurveData {
    degree: usize,
    knots: KnotVector,
    points: Vec<DVec4>,
}

impl TryFrom<CurveData> for NurbsCurve {
    type Error = NurbsError;

    fn try_from(data: CurveData) -> Result<Self, Self::Error> {
        Self::from_homogeneous(data.degree, data.knots.into_vec(), data.points)
    }
}

impl From<NurbsCurve> for CurveData {
    fn from(curve: NurbsCurve) -> Self {
        Self {
            degree: curve.degree,
            knots: curve.knots,
            points: curve.points,
        }
    }
}

/// Checks a list of Cartesian points and optional weights and lifts them to
/// homogeneous coordinates.
pub(crate) fn lift_points(
    points: &[DVec3],
    weights: Option<&[f64]>,
    index_offset: usize,
) -> Result<Vec<DVec4>, NurbsError> {
    if let Some(i) = points.iter().position(|p| !p.is_finite()) {
        return Err(NurbsError::NonFiniteControlPoint {
            index: index_offset + i,
        });
    }
    match weights {
        None => Ok(points.iter().map(|&p| homogenize(p, 1.0)).collect()),
        Some(weights) => {
            if weights.len() != points.len() {
                return Err(NurbsError::WeightCountMismatch {
                    expected: points.len(),
                    found: weights.len(),
                });
            }
            if let Some(i) = weights.iter().position(|w| !(w.is_finite() && *w > 0.0)) {
                return Err(NurbsError::InvalidWeight {
                    index: index_offset + i,
                    value: weights[i],
                });
            }
            Ok(points
                .iter()
                .zip(weights)
                .map(|(&p, &w)| homogenize(p, w))
                .collect())
        }
    }
}

/// Checks homogeneous control points: positive finite weights and finite
/// Cartesian coordinates.
pub(crate) fn check_homogeneous(points: &[DVec4]) -> Result<(), NurbsError> {
    for (index, pw) in points.iter().enumerate() {
        if !(pw.w.is_finite() && pw.w > 0.0) {
            return Err(NurbsError::InvalidWeight { index, value: pw.w });
        }
        if !pw.is_finite() || !dehomogenize(*pw).is_finite() {
            return Err(NurbsError::NonFiniteControlPoint { index });
        }
    }
    Ok(())
}

/// Knot vector and per-control-point data of a circular arc template
/// (Algorithm A7.1): for every control point the angle offset from the start,
/// the radius scale factor and the weight.
pub(crate) struct ArcTemplate {
    /// Clamped quadratic knot vector on `[0, 1]`.
    pub knots: Vec<f64>,
    /// `(angle offset, radius scale, weight)` per control point.
    pub frames: Vec<(f64, f64, f64)>,
}

/// Builds the arc template for a sweep in `(0, 2 pi]` (Algorithm A7.1): the
/// sweep is split into at most four equal pieces of at most 90 degrees, each a
/// rational quadratic Bezier arc whose middle weight is `cos(dtheta / 2)`.
pub(crate) fn arc_template(sweep: f64) -> ArcTemplate {
    let slack = 1.0 + 1e-12;
    let narcs: usize = if sweep <= FRAC_PI_2 * slack {
        1
    } else if sweep <= PI * slack {
        2
    } else if sweep <= 1.5 * PI * slack {
        3
    } else {
        4
    };
    let dtheta = sweep / narcs as f64;
    let w1 = (0.5 * dtheta).cos();
    let mut frames = Vec::with_capacity(2 * narcs + 1);
    frames.push((0.0, 1.0, 1.0));
    for i in 1..=narcs {
        let end = if i == narcs { sweep } else { dtheta * i as f64 };
        frames.push((end - 0.5 * dtheta, 1.0 / w1, w1));
        frames.push((end, 1.0, 1.0));
    }
    let mut knots = vec![0.0; 3];
    for i in 1..narcs {
        let k = i as f64 / narcs as f64;
        knots.push(k);
        knots.push(k);
    }
    knots.extend([1.0; 3]);
    ArcTemplate { knots, frames }
}

/// Checks an angle for an arc or revolution and normalizes it into a sweep in
/// `(0, 2 pi]`.
pub(crate) fn normalize_sweep(sweep: f64) -> Result<f64, NurbsError> {
    if !sweep.is_finite() {
        return Err(NurbsError::InvalidArgument("sweep angle is not finite"));
    }
    if sweep <= 1e-12 {
        return Err(NurbsError::InvalidArgument("sweep angle is zero"));
    }
    if sweep > TAU * (1.0 + 1e-12) {
        return Err(NurbsError::InvalidArgument(
            "sweep angle exceeds a full turn",
        ));
    }
    Ok(sweep.min(TAU))
}

/// Rational derivatives from homogeneous ones (Algorithm A4.2, RatCurveDerivs).
fn rational_derivatives(aders: &[DVec4]) -> Vec<DVec3> {
    let mut ck: Vec<DVec3> = Vec::with_capacity(aders.len());
    let w0 = aders[0].w;
    for k in 0..aders.len() {
        let mut v = aders[k].truncate();
        for i in 1..=k {
            v -= ck[k - i] * (binomial(k, i) * aders[i].w);
        }
        ck.push(v / w0);
    }
    ck
}

impl NurbsCurve {
    /// Builds a curve from Cartesian control points and optional weights
    /// (all weights default to 1, giving a non-rational B-spline).
    ///
    /// Requires `degree >= 1`, `points.len() >= degree + 1`,
    /// `knots.len() == points.len() + degree + 1`, a non-decreasing finite
    /// knot vector with non-empty domain and admissible multiplicities, and
    /// positive finite weights.
    pub fn new(
        degree: usize,
        knots: Vec<f64>,
        points: Vec<DVec3>,
        weights: Option<Vec<f64>>,
    ) -> Result<Self, NurbsError> {
        let pw = lift_points(&points, weights.as_deref(), 0)?;
        Self::from_homogeneous(degree, knots, pw)
    }

    /// Builds a curve directly from homogeneous control points
    /// `(w x, w y, w z, w)`, validating like [`NurbsCurve::new`].
    pub fn from_homogeneous(
        degree: usize,
        knots: Vec<f64>,
        points: Vec<DVec4>,
    ) -> Result<Self, NurbsError> {
        let knots = KnotVector::new(knots)?;
        knots.validate_for(degree, points.len())?;
        check_homogeneous(&points)?;
        Ok(Self {
            degree,
            knots,
            points,
        })
    }

    /// Straight line segment from `a` (parameter 0) to `b` (parameter 1) as a
    /// degree-1 curve.
    pub fn line(a: DVec3, b: DVec3) -> Result<Self, NurbsError> {
        Self::new(1, vec![0.0, 0.0, 1.0, 1.0], vec![a, b], None)
    }

    /// Internal constructor for results that are valid by construction.
    pub(crate) fn from_parts(degree: usize, knots: Vec<f64>, points: Vec<DVec4>) -> Self {
        let curve = Self {
            degree,
            knots: KnotVector::from_vec_unchecked(knots),
            points,
        };
        debug_assert!(curve.knots.validate_for(degree, curve.points.len()).is_ok());
        curve
    }

    /// Polynomial degree `p`.
    pub fn degree(&self) -> usize {
        self.degree
    }

    /// The knot vector.
    pub fn knots(&self) -> &KnotVector {
        &self.knots
    }

    /// Cartesian control points (homogeneous points divided by their weight).
    pub fn control_points(&self) -> Vec<DVec3> {
        self.points.iter().map(|&pw| dehomogenize(pw)).collect()
    }

    /// Control point weights.
    pub fn weights(&self) -> Vec<f64> {
        self.points.iter().map(|pw| pw.w).collect()
    }

    /// Homogeneous control points `(w x, w y, w z, w)`.
    pub fn homogeneous_points(&self) -> &[DVec4] {
        &self.points
    }

    /// `true` if the weights are not all equal (a curve with equal weights is
    /// a plain polynomial B-spline).
    pub fn is_rational(&self) -> bool {
        let w0 = self.points[0].w;
        self.points.iter().any(|pw| pw.w != w0)
    }

    /// Parametric domain `(u_p, u_{n+1})`.
    pub fn domain(&self) -> (f64, f64) {
        self.knots.domain(self.degree)
    }

    fn last_index(&self) -> usize {
        self.points.len() - 1
    }

    fn clamp_param(&self, t: f64) -> f64 {
        let (a, b) = self.domain();
        t.clamp(a, b)
    }

    /// Point in homogeneous space (Algorithm A4.1, before the projection).
    fn point_homogeneous(&self, t: f64) -> DVec4 {
        let p = self.degree;
        let u = self.clamp_param(t);
        let span = self.knots.find_span(p, self.last_index(), u);
        let basis = self.knots.basis_funs(span, u, p);
        basis.iter().enumerate().fold(DVec4::ZERO, |acc, (j, &b)| {
            acc + self.points[span - p + j] * b
        })
    }

    /// Point at parameter `t` (Algorithm A4.1, CurvePoint for rational curves).
    /// `t` is clamped to the domain.
    pub fn point(&self, t: f64) -> DVec3 {
        dehomogenize(self.point_homogeneous(t))
    }

    /// Derivatives of the homogeneous curve `C^w` up to order `d`
    /// (Algorithm A3.2, CurveDerivsAlg1), evaluating the polynomial piece of
    /// knot span `span` at `u` (also valid at, or slightly beyond, the span
    /// ends: it is the analytic continuation of that piece).
    ///
    /// When the weights active on the span are all equal, the derivatives of
    /// the weight function are exactly zero (the basis derivatives sum to
    /// zero); they are then stored as exact zeros so that A4.2 reduces to the
    /// polynomial case without picking up rounding noise.
    fn homogeneous_derivatives_in_span(&self, u: f64, d: usize, span: usize) -> Vec<DVec4> {
        let p = self.degree;
        let du = d.min(p);
        let nders = self.knots.ders_basis_funs(span, u, p, du);
        let active = &self.points[span - p..=span];
        let polynomial = active.iter().all(|pw| pw.w == active[0].w);
        let mut ck = vec![DVec4::ZERO; d + 1];
        for (k, ders_k) in nders.iter().enumerate().take(du + 1) {
            ck[k] = ders_k
                .iter()
                .zip(active)
                .fold(DVec4::ZERO, |acc, (&n, &pw)| acc + pw * n);
            if polynomial && k > 0 {
                ck[k].w = 0.0;
            }
        }
        ck
    }

    /// Point and derivatives up to order `n` at `t`: element `k` is
    /// `C^(k)(t)`, element 0 the point itself.
    ///
    /// Computed with Algorithm A3.2 on the homogeneous curve followed by
    /// Algorithm A4.2 (RatCurveDerivs). At knots the derivatives of the span
    /// to the right are returned (left span at the end of the domain). `t` is
    /// clamped to the domain.
    pub fn derivatives(&self, t: f64, n: usize) -> Vec<DVec3> {
        let u = self.clamp_param(t);
        let span = self.knots.find_span(self.degree, self.last_index(), u);
        rational_derivatives(&self.homogeneous_derivatives_in_span(u, n, span))
    }

    /// Like [`NurbsCurve::derivatives`] but always evaluates the polynomial
    /// piece of knot span `span` (one-sided derivatives at knots).
    fn derivatives_in_span(&self, u: f64, n: usize, span: usize) -> Vec<DVec3> {
        rational_derivatives(&self.homogeneous_derivatives_in_span(u, n, span))
    }

    /// Typical magnitude of coordinates, used to scale tolerances.
    fn coordinate_scale(&self) -> f64 {
        let (lo, hi) = bounding_box(&self.points);
        (hi - lo).length().max(lo.abs().max(hi.abs()).max_element())
    }

    /// Unit tangent at `t`.
    ///
    /// Where the first derivative vanishes (for example at a doubled control
    /// point) the direction of the first non-vanishing higher derivative is
    /// used, which is the one-sided limit of the tangent direction. Returns
    /// [`DVec3::ZERO`] only if every derivative up to order `max(degree, 2)`
    /// vanishes.
    pub fn tangent(&self, t: f64) -> DVec3 {
        let order = self.degree.max(2);
        let ders = self.derivatives(t, order);
        let (a, b) = self.domain();
        let width = b - a;
        let scale = self.coordinate_scale();
        for (k, d) in ders.iter().enumerate().skip(1) {
            let threshold = 1e-12 * scale / width.powi(k as i32);
            if d.length() > threshold {
                if let Some(u) = unit(*d) {
                    return u;
                }
            }
        }
        DVec3::ZERO
    }

    /// Inserts the knot `t` up to `times` times (Algorithm A5.1, CurveKnotIns).
    ///
    /// The curve is geometrically and parametrically unchanged. The number of
    /// insertions is capped so that the multiplicity of `t` does not exceed
    /// the degree; `t` outside `[u_p, u_{n+1})` leaves the curve unchanged.
    pub fn insert_knot(&self, t: f64, times: usize) -> NurbsCurve {
        let (a, b) = self.domain();
        if times == 0 || !(t >= a && t < b) {
            return self.clone();
        }
        let p = self.degree;
        let s = self.knots.multiplicity(t);
        if s >= p {
            return self.clone();
        }
        let r = times.min(p - s);
        let np = self.last_index();
        let k = self.knots.find_span(p, np, t);
        let up = self.knots.as_slice();
        let pw = &self.points;
        let mp = np + p + 1;
        let nq = np + r;

        let mut uq = vec![0.0; mp + r + 1];
        uq[..=k].copy_from_slice(&up[..=k]);
        uq[k + 1..=k + r].fill(t);
        uq[k + 1 + r..].copy_from_slice(&up[k + 1..]);

        let mut qw = vec![DVec4::ZERO; nq + 1];
        qw[..=k - p].copy_from_slice(&pw[..=k - p]);
        qw[k - s + r..].copy_from_slice(&pw[k - s..]);
        let mut rw: Vec<DVec4> = pw[k - p..=k - s].to_vec();
        let mut l = k - p;
        for j in 1..=r {
            l = k - p + j;
            for i in 0..=p - j - s {
                let alpha = (t - up[l + i]) / (up[i + k + 1] - up[l + i]);
                rw[i] = rw[i + 1] * alpha + rw[i] * (1.0 - alpha);
            }
            qw[l] = rw[0];
            qw[k + r - j - s] = rw[p - j - s];
        }
        if k - s > l + 1 {
            qw[l + 1..k - s].copy_from_slice(&rw[1..k - s - l]);
        }
        Self::from_parts(p, uq, qw)
    }

    /// Splits the curve at `t` into the pieces on `[a, t]` and `[t, b]`,
    /// each keeping the original parametrization.
    ///
    /// Implemented by inserting `t` until its multiplicity equals the degree
    /// (Algorithm A5.1) and distributing the control points. A `t` within a
    /// few rounding steps of an interior knot splits at that knot instead
    /// (the pieces then meet at that knot's parameter): splitting beside it
    /// would leave a knot span of width ~1e-17 at a piece's end, on which
    /// end derivatives and quadrature mean nothing. The split point moves
    /// along the curve by a rounding-level distance only.
    ///
    /// # Panics
    ///
    /// Panics if `t` is not strictly inside the domain; use
    /// [`NurbsCurve::try_split`] for a non-panicking variant.
    pub fn split(&self, t: f64) -> (NurbsCurve, NurbsCurve) {
        match self.try_split(t) {
            Ok(pieces) => pieces,
            Err(_) => panic!("split parameter {t} is not strictly inside the curve domain"),
        }
    }

    /// Like [`NurbsCurve::split`] but returns an error when `t` is not
    /// strictly inside the domain.
    pub fn try_split(&self, t: f64) -> Result<(NurbsCurve, NurbsCurve), NurbsError> {
        let (a, b) = self.domain();
        if !(t > a && t < b) {
            return Err(NurbsError::InvalidArgument(
                "split parameter must lie strictly inside the domain",
            ));
        }
        let p = self.degree;
        let near_knot = 8.0 * f64::EPSILON * a.abs().max(b.abs()).max(b - a);
        let t = self
            .knots
            .as_slice()
            .iter()
            .copied()
            .filter(|&k| k > a && k < b && (k - t).abs() <= near_knot)
            .min_by(|x, y| (x - t).abs().total_cmp(&(y - t).abs()))
            .unwrap_or(t);
        let s = self.knots.multiplicity(t);
        let refined = if s < p {
            self.insert_knot(t, p - s)
        } else {
            self.clone()
        };
        let u = refined.knots.as_slice();
        // First index of t; t occupies indices first..first + p.
        let first = u.partition_point(|&k| k < t);
        let mut left_knots = u[..first + p].to_vec();
        left_knots.push(t);
        let left_points = refined.points[..first].to_vec();
        let mut right_knots = Vec::with_capacity(u.len() - first + 1);
        right_knots.push(t);
        right_knots.extend_from_slice(&u[first..]);
        let right_points = refined.points[first - 1..].to_vec();
        Ok((
            Self::from_parts(p, left_knots, left_points),
            Self::from_parts(p, right_knots, right_points),
        ))
    }

    /// The same curve traversed in the opposite direction, on the same domain:
    /// `reversed().point(a + b - t) == point(t)`.
    pub fn reversed(&self) -> NurbsCurve {
        let (a, b) = self.domain();
        let sum = a + b;
        let knots = self
            .knots
            .as_slice()
            .iter()
            .rev()
            .map(|&k| {
                // Reflect, keeping the domain end points exact and the
                // sequence monotone despite rounding.
                if k == a {
                    b
                } else if k == b {
                    a
                } else {
                    let r = sum - k;
                    if k < a {
                        r.max(b)
                    } else if k > b {
                        r.min(a)
                    } else {
                        r.clamp(a, b)
                    }
                }
            })
            .collect();
        let points = self.points.iter().rev().copied().collect();
        Self::from_parts(self.degree, knots, points)
    }

    /// Applies the affine transformation `m` (its projective row is ignored).
    /// Weights and parametrization are unchanged, so the result is exact.
    /// `m` must be finite (asserted in debug builds).
    pub fn transformed(&self, m: &DMat4) -> NurbsCurve {
        debug_assert!(m.is_finite(), "transform matrix must be finite");
        Self {
            degree: self.degree,
            knots: self.knots.clone(),
            points: self
                .points
                .iter()
                .map(|&pw| transform_homogeneous(m, pw))
                .collect(),
        }
    }

    /// Axis-aligned bounding box `(min, max)` of the control points, which
    /// contains the curve by the convex hull property.
    pub fn bounding_box(&self) -> (DVec3, DVec3) {
        bounding_box(&self.points)
    }

    /// Exact rational quadratic circular arc (Algorithm A7.1, MakeNurbsCircle).
    ///
    /// The arc lies in the plane spanned by `x_axis` and `y_axis` through
    /// `center` and runs counter-clockwise (from `x_axis` towards `y_axis`)
    /// from `start_angle` to `end_angle` (radians). `x_axis` is normalized and
    /// `y_axis` is made orthonormal to it. As in A7.1, an `end_angle` below
    /// `start_angle` is increased by `2 pi`. The sweep must be in `(0, 2 pi]`;
    /// a full circle is closed exactly. The domain is `[0, 1]`, split into at
    /// most four Bezier pieces of at most 90 degrees.
    pub fn circular_arc(
        center: DVec3,
        x_axis: DVec3,
        y_axis: DVec3,
        radius: f64,
        start_angle: f64,
        end_angle: f64,
    ) -> Result<Self, NurbsError> {
        if !center.is_finite() {
            return Err(NurbsError::InvalidArgument("arc center is not finite"));
        }
        if !(radius.is_finite() && radius > 0.0) {
            return Err(NurbsError::InvalidArgument(
                "arc radius must be positive and finite",
            ));
        }
        if !(start_angle.is_finite() && end_angle.is_finite()) {
            return Err(NurbsError::InvalidArgument("arc angles must be finite"));
        }
        let x = unit(x_axis).ok_or(NurbsError::InvalidArgument("arc x axis is degenerate"))?;
        if !y_axis.is_finite() {
            return Err(NurbsError::InvalidArgument("arc y axis is not finite"));
        }
        let y_perp = y_axis - x * y_axis.dot(x);
        if y_perp.length() <= 1e-12 * y_axis.length() {
            return Err(NurbsError::InvalidArgument(
                "arc y axis is zero or parallel to the x axis",
            ));
        }
        let y = unit(y_perp).ok_or(NurbsError::InvalidArgument("arc y axis is degenerate"))?;
        let mut sweep = end_angle - start_angle;
        if sweep < 0.0 {
            sweep += TAU;
        }
        let sweep = normalize_sweep(sweep)?;
        let template = arc_template(sweep);
        let mut points: Vec<DVec4> = template
            .frames
            .iter()
            .map(|&(offset, scale, w)| {
                let angle = start_angle + offset;
                let p = center + (x * angle.cos() + y * angle.sin()) * (radius * scale);
                homogenize(p, w)
            })
            .collect();
        if sweep == TAU {
            let last = points.len() - 1;
            points[last] = points[0];
        }
        Ok(Self::from_parts(2, template.knots, points))
    }

    /// Global interpolation through `points` (Algorithm A9.1).
    ///
    /// Uses chord-length parameters on `[0, 1]`, knots by averaging (Eq. 9.8)
    /// and a clamped knot vector; the curve passes through every point, the
    /// `k`-th at its chord-length parameter. Needs at least `degree + 1`
    /// points, no two consecutive ones equal.
    pub fn interpolate(points: &[DVec3], degree: usize) -> Result<Self, NurbsError> {
        if degree == 0 {
            return Err(NurbsError::ZeroDegree);
        }
        if points.len() < degree + 1 {
            return Err(NurbsError::TooFewPoints {
                required: degree + 1,
                found: points.len(),
            });
        }
        if let Some(index) = points.iter().position(|q| !q.is_finite()) {
            return Err(NurbsError::NonFiniteControlPoint { index });
        }
        let p = degree;
        let n = points.len() - 1;
        let chords: Vec<f64> = points.windows(2).map(|w| (w[1] - w[0]).length()).collect();
        let total: f64 = chords.iter().sum();
        if !(total.is_finite() && total > 0.0) {
            return Err(NurbsError::InvalidArgument(
                "interpolation points must not all coincide",
            ));
        }
        if chords.iter().any(|&c| c <= 0.0) {
            return Err(NurbsError::InvalidArgument(
                "consecutive interpolation points must be distinct",
            ));
        }
        // Chord-length parameters (Eq. 9.5).
        let mut params = vec![0.0; n + 1];
        let mut acc = 0.0;
        for k in 1..n {
            acc += chords[k - 1];
            params[k] = acc / total;
        }
        params[n] = 1.0;
        // Knots by averaging (Eq. 9.8).
        let m = n + p + 1;
        let mut knots = vec![0.0; m + 1];
        for j in 1..=n - p {
            knots[j + p] = params[j..j + p].iter().sum::<f64>() / p as f64;
        }
        knots[m - p..].fill(1.0);
        let kv = KnotVector::new(knots)?;
        kv.validate_for(p, n + 1)?;
        // Banded collocation matrix A[k][span - p + j] = N_{span-p+j,p}(params[k]).
        let mut rows = Vec::with_capacity(n + 1);
        let mut lo = Vec::with_capacity(n + 1);
        for &t in &params {
            let span = kv.find_span(p, n, t);
            rows.push(kv.basis_funs(span, t, p));
            lo.push(span - p);
        }
        let ctrl = solve_banded(rows, &lo, points.to_vec()).ok_or(NurbsError::SingularSystem)?;
        let pw = lift_points(&ctrl, None, 0)?;
        Ok(Self {
            degree: p,
            knots: kv,
            points: pw,
        })
    }

    /// Closest point on the curve to `p`, as `(parameter, point)`.
    ///
    /// Global search followed by Newton point inversion:
    ///
    /// 1. Every non-empty knot span gets a lower bound on its distance to `p`
    ///    from the bounding box of its `p + 1` active control points (strong
    ///    convex hull property). Spans are visited nearest bound first and
    ///    skipped once the bound cannot beat the best point found so far.
    /// 2. A visited span is sampled densely, uniformly in the parameter and
    ///    then bisected where consecutive samples are far apart or the
    ///    tangent turns by more than 30 degrees between them, so strongly
    ///    non-uniform rational parametrizations and near-cusps are resolved.
    /// 3. Each sampled local minimum of the distance within the span (best
    ///    first, plateaus reduced to their ends) is refined by Newton point
    ///    inversion (Eqs. 6.3 and 6.4) on that span's polynomial piece, with
    ///    the parameter clamped to the span (Eq. 6.5) and a backtracking
    ///    safeguard. Minima at knots with tangent discontinuities are found
    ///    because each span is treated on its own.
    ///
    /// The result is never farther from `p` than the best sample.
    pub fn closest_point(&self, p: DVec3) -> (f64, DVec3) {
        let deg = self.degree;
        let breaks = self.knots.breakpoints(deg);
        let mut intervals: Vec<(f64, f64, usize, f64)> = breaks
            .windows(2)
            .map(|w| {
                let span = self
                    .knots
                    .find_span(deg, self.last_index(), 0.5 * (w[0] + w[1]));
                let (lo, hi) = bounding_box(&self.points[span - deg..=span]);
                (w[0], w[1], span, aabb_distance_squared(lo, hi, p))
            })
            .collect();
        intervals.sort_by(|x, y| x.3.total_cmp(&y.3));
        let (a, _) = self.domain();
        let start = self.point(a);
        let mut best = (a, start, start.distance_squared(p));
        for &(lo, hi, span, bound) in &intervals {
            if bound >= best.2 {
                break;
            }
            self.closest_in_span(p, lo, hi, span, &mut best);
        }
        (best.0, best.1)
    }

    /// Samples knot span `[lo, hi]` and refines its local minima, updating
    /// `best = (t, point, squared distance)`.
    fn closest_in_span(
        &self,
        target: DVec3,
        lo: f64,
        hi: f64,
        span: usize,
        best: &mut (f64, DVec3, f64),
    ) {
        let per_span = 4 * (self.degree + 1) + 8;
        // Point and unit tangent (if defined) of the span's piece at `t`.
        let sample = |t: f64| {
            let d = self.derivatives_in_span(t, 1, span);
            (d[0], unit(d[1]))
        };
        let mut ts: Vec<f64> = (0..per_span)
            .map(|k| lo + (hi - lo) * (k as f64 / per_span as f64))
            .chain(std::iter::once(hi))
            .collect();
        let mut pts: Vec<(DVec3, Option<DVec3>)> = ts.iter().map(|&t| sample(t)).collect();
        // Bisect where consecutive samples are far apart (strongly
        // non-uniform rational parametrizations) or where the tangent turns
        // sharply (near-cusps, where two distance minima can hide between
        // samples), so that samples resolve the geometry evenly.
        let chord_sum: f64 = pts.windows(2).map(|w| w[0].0.distance(w[1].0)).sum();
        let max_chord = 2.0 * chord_sum / per_span as f64;
        let min_cos_turn = (30.0f64).to_radians().cos();
        let needs_split = |a: &(DVec3, Option<DVec3>), b: &(DVec3, Option<DVec3>)| {
            let turn = match (a.1, b.1) {
                (Some(ta), Some(tb)) => ta.dot(tb) < min_cos_turn,
                _ => false,
            };
            turn || a.0.distance(b.0) > max_chord
        };
        let max_samples = 16 * per_span;
        for _ in 0..16 {
            if ts.len() >= max_samples {
                break;
            }
            let mut new_ts = Vec::with_capacity(2 * ts.len());
            let mut new_pts = Vec::with_capacity(2 * ts.len());
            let mut added = 0;
            for k in 0..ts.len() {
                new_ts.push(ts[k]);
                new_pts.push(pts[k]);
                if k + 1 < ts.len()
                    && ts.len() + added < max_samples
                    && needs_split(&pts[k], &pts[k + 1])
                {
                    let tm = 0.5 * (ts[k] + ts[k + 1]);
                    if tm > ts[k] && tm < ts[k + 1] {
                        new_ts.push(tm);
                        new_pts.push(sample(tm));
                        added += 1;
                    }
                }
            }
            ts = new_ts;
            pts = new_pts;
            if added == 0 {
                break;
            }
        }
        let pts: Vec<DVec3> = pts.into_iter().map(|(q, _)| q).collect();
        let dist2: Vec<f64> = pts.iter().map(|q| q.distance_squared(target)).collect();
        for (k, &d) in dist2.iter().enumerate() {
            if d < best.2 {
                *best = (ts[k], pts[k], d);
            }
        }
        let mut candidates = local_minima(&dist2);
        candidates.sort_by(|&i, &j| dist2[i].total_cmp(&dist2[j]));
        candidates.truncate(MAX_CANDIDATES_PER_SPAN);
        for c in candidates {
            let (t, q) = self.refine_closest(target, ts[c], lo, hi, span);
            let d = q.distance_squared(target);
            if d < best.2 {
                *best = (t, q, d);
            }
        }
    }

    /// Safeguarded Newton point inversion (Eqs. 6.3-6.5) starting at `t0`,
    /// restricted to `[lo, hi]` and evaluating the piece of knot span `span`.
    fn refine_closest(
        &self,
        target: DVec3,
        t0: f64,
        lo: f64,
        hi: f64,
        span: usize,
    ) -> (f64, DVec3) {
        let width = hi - lo;
        let mut t = t0.clamp(lo, hi);
        let mut ders = self.derivatives_in_span(t, 2, span);
        let mut r = ders[0] - target;
        let mut dist2 = r.length_squared();
        let mut f = ders[1].dot(r);
        for _ in 0..MAX_NEWTON_ITERATIONS {
            if dist2 == 0.0 || f == 0.0 {
                break;
            }
            let speed2 = ders[1].length_squared();
            let fp = ders[2].dot(r) + speed2;
            let step = if fp > 0.0 {
                -f / fp
            } else if speed2 > 0.0 {
                // Away from a minimum: Gauss-Newton step, a descent direction.
                -f / speed2
            } else {
                break;
            };
            // Rounding-level tie allowance so Newton can keep polishing the
            // parameter once the distance stops changing measurably.
            let scale = ders[0].abs().max(target.abs()).max_element();
            let slack = 1e-14 * (dist2 + dist2.sqrt() * scale);
            let mut alpha = 1.0;
            let mut accepted = None;
            for _ in 0..40 {
                let tn = (t + alpha * step).clamp(lo, hi);
                if tn == t {
                    break;
                }
                let dn = self.derivatives_in_span(tn, 2, span);
                let rn = dn[0] - target;
                let d2n = rn.length_squared();
                let fnew = dn[1].dot(rn);
                if d2n < dist2 || (d2n <= dist2 + slack && fnew.abs() < f.abs()) {
                    accepted = Some((tn, dn, rn, d2n, fnew));
                    break;
                }
                alpha *= 0.5;
            }
            let Some((tn, dn, rn, d2n, fnew)) = accepted else {
                break;
            };
            let moved = (tn - t).abs();
            t = tn;
            ders = dn;
            r = rn;
            dist2 = d2n;
            f = fnew;
            if moved <= 1e-15 * width.max(t.abs()) {
                break;
            }
        }
        (t, ders[0])
    }

    /// Length of the curve between parameters `t0` and `t1` (order does not
    /// matter; both are clamped to the domain).
    ///
    /// Integrates the speed `|C'(t)|` with globally adaptive 16-point
    /// Gauss-Legendre quadrature whose initial panels are the knot spans, on
    /// which the speed is smooth.
    pub fn arc_length(&self, t0: f64, t1: f64) -> f64 {
        let (lo, hi) = if t0 <= t1 { (t0, t1) } else { (t1, t0) };
        let lo = self.clamp_param(lo);
        let hi = self.clamp_param(hi);
        if lo.is_nan() || hi.is_nan() || hi <= lo {
            return 0.0;
        }
        let ctrl = self.control_points();
        let polygon: f64 = ctrl.windows(2).map(|w| (w[1] - w[0]).length()).sum();
        // Relative accuracy on the length, but never below the rounding noise
        // of derivatives computed from coordinates of this magnitude.
        let tol = 1e-14 * polygon + 64.0 * f64::EPSILON * self.coordinate_scale();
        let breaks = self.knots.breakpoints(self.degree);
        let intervals: Vec<(f64, f64)> = breaks
            .windows(2)
            .map(|w| (w[0].max(lo), w[1].min(hi)))
            .filter(|(a, b)| b > a)
            .collect();
        let mut speed = |t: f64| self.derivatives(t, 1)[1].length();
        integrate(&mut speed, &intervals, tol)
    }

    /// Total length of the curve over its domain.
    pub fn length(&self) -> f64 {
        let (a, b) = self.domain();
        self.arc_length(a, b)
    }

    /// Unit normal of the plane containing the curve, if all control points
    /// lie within `tol` of a common plane (the curve then lies in it by the
    /// convex hull property).
    ///
    /// Returns `None` if the curve is not planar, and also when no unique
    /// plane exists because the control points are collinear within `tol`.
    pub fn is_planar(&self, tol: f64) -> Option<DVec3> {
        let pts = self.control_points();
        let p0 = pts[0];
        let p1 = pts
            .iter()
            .copied()
            .max_by(|x, y| x.distance_squared(p0).total_cmp(&y.distance_squared(p0)))?;
        let dir = unit(p1 - p0)?;
        if p1.distance(p0) <= tol {
            return None;
        }
        let off_line = |q: DVec3| {
            let v = q - p0;
            (v - dir * v.dot(dir)).length()
        };
        let p2 = pts
            .iter()
            .copied()
            .max_by(|x, y| off_line(*x).total_cmp(&off_line(*y)))?;
        if off_line(p2) <= tol {
            return None;
        }
        let normal = unit((p1 - p0).cross(p2 - p0))?;
        pts.iter()
            .all(|q| (*q - p0).dot(normal).abs() <= tol)
            .then_some(normal)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn quarter_circle() -> NurbsCurve {
        let w = std::f64::consts::FRAC_1_SQRT_2;
        NurbsCurve::new(
            2,
            vec![0., 0., 0., 1., 1., 1.],
            vec![
                DVec3::new(1., 0., 0.),
                DVec3::new(1., 1., 0.),
                DVec3::new(0., 1., 0.),
            ],
            Some(vec![1.0, w, 1.0]),
        )
        .unwrap()
    }

    #[test]
    fn quarter_circle_points() {
        let c = quarter_circle();
        for i in 0..=20 {
            let t = i as f64 / 20.0;
            assert!((c.point(t).length() - 1.0).abs() < 1e-15);
        }
        assert_eq!(c.point(0.0), DVec3::X);
        assert!((c.point(1.0) - DVec3::Y).length() < 1e-15);
        assert!(c.is_rational());
    }

    #[test]
    fn insertion_book_style() {
        let c = NurbsCurve::new(
            3,
            vec![0., 0., 0., 0., 1., 2., 3., 3., 3., 3.],
            vec![
                DVec3::new(0., 0., 0.),
                DVec3::new(1., 2., 0.),
                DVec3::new(3., 3., 1.),
                DVec3::new(4., 0., 2.),
                DVec3::new(6., 1., 0.),
                DVec3::new(7., 3., 1.),
            ],
            Some(vec![1.0, 2.0, 0.5, 1.0, 3.0, 1.0]),
        )
        .unwrap();
        let c2 = c.insert_knot(1.5, 2);
        assert_eq!(c2.knots().multiplicity(1.5), 2);
        assert_eq!(c2.homogeneous_points().len(), 8);
        for i in 0..=60 {
            let t = 3.0 * i as f64 / 60.0;
            assert!((c.point(t) - c2.point(t)).length() < 1e-12);
        }
        // Inserting at an existing knot up to full multiplicity.
        let c3 = c.insert_knot(2.0, 10);
        assert_eq!(c3.knots().multiplicity(2.0), 3);
        // Clamped end: nothing to insert.
        assert_eq!(c.insert_knot(0.0, 1), c);
        assert_eq!(c.insert_knot(3.0, 1), c);
    }

    #[test]
    fn tangent_at_doubled_control_point() {
        let c = NurbsCurve::new(
            2,
            vec![0., 0., 0., 1., 1., 1.],
            vec![DVec3::ZERO, DVec3::ZERO, DVec3::new(1., 1., 0.)],
            None,
        )
        .unwrap();
        let t = c.tangent(0.0);
        assert!((t - DVec3::new(1., 1., 0.).normalize()).length() < 1e-12);
    }

    #[test]
    fn planarity() {
        let c = quarter_circle();
        let n = c.is_planar(1e-12).unwrap();
        assert!((n.z.abs() - 1.0).abs() < 1e-15);
        let line = NurbsCurve::line(DVec3::ZERO, DVec3::X).unwrap();
        assert!(line.is_planar(1e-9).is_none());
        let twisted = NurbsCurve::new(
            2,
            vec![0., 0., 0., 0.5, 1., 1., 1.],
            vec![
                DVec3::ZERO,
                DVec3::X,
                DVec3::new(1., 1., 0.),
                DVec3::new(1., 1., 1.),
            ],
            None,
        )
        .unwrap();
        assert!(twisted.is_planar(1e-6).is_none());
    }
}
