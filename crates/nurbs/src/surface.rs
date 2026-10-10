//! Non-uniform rational B-spline surfaces.

use std::f64::consts::TAU;

use glam::{DMat4, DVec2, DVec3, DVec4};
use serde::{Deserialize, Serialize};

use crate::curve::{arc_template, check_homogeneous, lift_points, normalize_sweep};
use crate::knots::KnotVector;
use crate::util::{
    aabb_distance_squared, binomial, bounding_box, dehomogenize, homogenize, transform_homogeneous,
    unit,
};
use crate::{NurbsCurve, NurbsError};

/// Maximum number of Newton iterations in point inversion.
const MAX_NEWTON_ITERATIONS: usize = 64;
/// Maximum number of sampled local minima per knot-span cell refined by
/// point inversion.
const MAX_CANDIDATES_PER_CELL: usize = 6;
/// Number of best grid samples per cell refined in addition to the grid's
/// local minima.
const EXTRA_BEST_SAMPLES: usize = 3;
/// Target number of samples per parameter direction (over all knot spans)
/// for closest-point search; each span gets between 4 and `2 * degree + 6`
/// before adaptive refinement.
const SAMPLES_PER_DIRECTION: usize = 64;

/// A tensor-product NURBS surface of degrees `(p, q) >= (1, 1)`.
///
/// The control net has `n_u x n_v` homogeneous points `(w x, w y, w z, w)`,
/// stored row-major: point `(i, j)` (the `i`-th in the `u` direction and the
/// `j`-th in the `v` direction) is at index `i * n_v + j`.
///
/// Invariants (checked by every constructor and by deserialization): each
/// direction satisfies the same knot conditions as a [`NurbsCurve`], all
/// weights are positive and all coordinates are finite.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(try_from = "SurfaceData", into = "SurfaceData")]
pub struct NurbsSurface {
    degree_u: usize,
    degree_v: usize,
    knots_u: KnotVector,
    knots_v: KnotVector,
    count_u: usize,
    count_v: usize,
    points: Vec<DVec4>,
}

/// Unvalidated serialized form of a [`NurbsSurface`].
#[derive(Serialize, Deserialize)]
struct SurfaceData {
    degree_u: usize,
    degree_v: usize,
    knots_u: KnotVector,
    knots_v: KnotVector,
    count_u: usize,
    count_v: usize,
    points: Vec<DVec4>,
}

impl TryFrom<SurfaceData> for NurbsSurface {
    type Error = NurbsError;

    fn try_from(data: SurfaceData) -> Result<Self, Self::Error> {
        Self::from_homogeneous_flat(
            data.degree_u,
            data.degree_v,
            data.knots_u.into_vec(),
            data.knots_v.into_vec(),
            data.count_u,
            data.count_v,
            data.points,
        )
    }
}

impl From<NurbsSurface> for SurfaceData {
    fn from(s: NurbsSurface) -> Self {
        Self {
            degree_u: s.degree_u,
            degree_v: s.degree_v,
            knots_u: s.knots_u,
            knots_v: s.knots_v,
            count_u: s.count_u,
            count_v: s.count_v,
            points: s.points,
        }
    }
}

/// Checks that all rows of a nested net have the same length and returns it.
fn row_length<T>(rows: &[Vec<T>]) -> Result<usize, NurbsError> {
    let expected = rows.first().map_or(0, Vec::len);
    if let Some((row, r)) = rows.iter().enumerate().find(|(_, r)| r.len() != expected) {
        return Err(NurbsError::RaggedControlNet {
            row,
            expected,
            found: r.len(),
        });
    }
    Ok(expected)
}

/// Non-empty knot-span intervals of one direction: `(start, end, span)`.
fn span_intervals(knots: &KnotVector, degree: usize, count: usize) -> Vec<(f64, f64, usize)> {
    knots
        .breakpoints(degree)
        .windows(2)
        .map(|w| {
            (
                w[0],
                w[1],
                knots.find_span(degree, count - 1, 0.5 * (w[0] + w[1])),
            )
        })
        .collect()
}

/// `count + 1` evenly spaced parameters from `lo` to exactly `hi`.
fn even_params(lo: f64, hi: f64, count: usize) -> Vec<f64> {
    (0..count)
        .map(|k| lo + (hi - lo) * (k as f64 / count as f64))
        .chain(std::iter::once(hi))
        .collect()
}

/// Rational surface derivatives from homogeneous ones (Algorithm A4.4,
/// RatSurfaceDerivs).
fn rational_surface_derivatives(aders: &[Vec<DVec4>]) -> Vec<Vec<DVec3>> {
    let d = aders.len() - 1;
    let w00 = aders[0][0].w;
    let mut skl: Vec<Vec<DVec3>> = (0..=d).map(|k| vec![DVec3::ZERO; d + 1 - k]).collect();
    for k in 0..=d {
        for l in 0..=d - k {
            let mut val = aders[k][l].truncate();
            for j in 1..=l {
                val -= skl[k][l - j] * (binomial(l, j) * aders[0][j].w);
            }
            for i in 1..=k {
                val -= skl[k - i][l] * (binomial(k, i) * aders[i][0].w);
                let mut v2 = DVec3::ZERO;
                for j in 1..=l {
                    v2 += skl[k - i][l - j] * (binomial(l, j) * aders[i][j].w);
                }
                val -= v2 * binomial(k, i);
            }
            skl[k][l] = val / w00;
        }
    }
    skl
}

/// Parameter box of one knot-span cell, with its span indices.
#[derive(Clone, Copy)]
struct Cell {
    u: (f64, f64),
    v: (f64, f64),
    uspan: usize,
    vspan: usize,
}

/// Solves the symmetric 2x2 system `[[a, b], [b, c]] x = -g`, or `None` if
/// the matrix is not positive definite.
fn solve_spd2(a: f64, b: f64, c: f64, g: DVec2) -> Option<DVec2> {
    let det = a * c - b * b;
    if a > 0.0 && c > 0.0 && det > 1e-14 * a * c {
        Some(DVec2::new(
            -(c * g.x - b * g.y) / det,
            -(a * g.y - b * g.x) / det,
        ))
    } else {
        None
    }
}

impl NurbsSurface {
    /// Builds a surface from a Cartesian control net `points[i][j]` (`i` along
    /// `u`, `j` along `v`) and optional weights of the same shape (default 1).
    ///
    /// Each direction must satisfy the curve conditions of
    /// [`NurbsCurve::new`]; rows must all have the same length.
    pub fn new(
        degree_u: usize,
        degree_v: usize,
        knots_u: Vec<f64>,
        knots_v: Vec<f64>,
        points: Vec<Vec<DVec3>>,
        weights: Option<Vec<Vec<f64>>>,
    ) -> Result<Self, NurbsError> {
        let count_u = points.len();
        let count_v = row_length(&points)?;
        if let Some(w) = &weights {
            if w.len() != count_u {
                return Err(NurbsError::WeightCountMismatch {
                    expected: count_u * count_v,
                    found: w.iter().map(Vec::len).sum(),
                });
            }
            if let Some((row, r)) = w.iter().enumerate().find(|(_, r)| r.len() != count_v) {
                return Err(NurbsError::RaggedControlNet {
                    row,
                    expected: count_v,
                    found: r.len(),
                });
            }
        }
        let mut flat = Vec::with_capacity(count_u * count_v);
        for (i, row) in points.iter().enumerate() {
            let w = weights.as_ref().map(|w| w[i].as_slice());
            flat.extend(lift_points(row, w, i * count_v)?);
        }
        Self::from_homogeneous_flat(degree_u, degree_v, knots_u, knots_v, count_u, count_v, flat)
    }

    /// Builds a surface from a homogeneous control net `points[i][j]` of
    /// `(w x, w y, w z, w)` values, validating like [`NurbsSurface::new`].
    pub fn from_homogeneous(
        degree_u: usize,
        degree_v: usize,
        knots_u: Vec<f64>,
        knots_v: Vec<f64>,
        points: Vec<Vec<DVec4>>,
    ) -> Result<Self, NurbsError> {
        let count_u = points.len();
        let count_v = row_length(&points)?;
        let flat = points.into_iter().flatten().collect();
        Self::from_homogeneous_flat(degree_u, degree_v, knots_u, knots_v, count_u, count_v, flat)
    }

    fn from_homogeneous_flat(
        degree_u: usize,
        degree_v: usize,
        knots_u: Vec<f64>,
        knots_v: Vec<f64>,
        count_u: usize,
        count_v: usize,
        points: Vec<DVec4>,
    ) -> Result<Self, NurbsError> {
        let knots_u = KnotVector::new(knots_u)?;
        let knots_v = KnotVector::new(knots_v)?;
        knots_u.validate_for(degree_u, count_u)?;
        knots_v.validate_for(degree_v, count_v)?;
        if points.len() != count_u * count_v {
            return Err(NurbsError::RaggedControlNet {
                row: count_u,
                expected: count_u * count_v,
                found: points.len(),
            });
        }
        check_homogeneous(&points)?;
        Ok(Self {
            degree_u,
            degree_v,
            knots_u,
            knots_v,
            count_u,
            count_v,
            points,
        })
    }

    /// Bilinear patch with `S(0,0) = p00`, `S(1,0) = p10`, `S(0,1) = p01` and
    /// `S(1,1) = p11` on the domain `[0, 1] x [0, 1]`.
    pub fn bilinear(p00: DVec3, p10: DVec3, p01: DVec3, p11: DVec3) -> Result<Self, NurbsError> {
        Self::new(
            1,
            1,
            vec![0.0, 0.0, 1.0, 1.0],
            vec![0.0, 0.0, 1.0, 1.0],
            vec![vec![p00, p01], vec![p10, p11]],
            None,
        )
    }

    /// Degree in the `u` direction.
    pub fn degree_u(&self) -> usize {
        self.degree_u
    }

    /// Degree in the `v` direction.
    pub fn degree_v(&self) -> usize {
        self.degree_v
    }

    /// Knot vector in the `u` direction.
    pub fn knots_u(&self) -> &KnotVector {
        &self.knots_u
    }

    /// Knot vector in the `v` direction.
    pub fn knots_v(&self) -> &KnotVector {
        &self.knots_v
    }

    /// Control net dimensions `(n_u, n_v)`.
    pub fn control_net_size(&self) -> (usize, usize) {
        (self.count_u, self.count_v)
    }

    /// Homogeneous control net, row-major (`(i, j)` at `i * n_v + j`).
    pub fn homogeneous_points(&self) -> &[DVec4] {
        &self.points
    }

    /// Cartesian control net `[i][j]`.
    pub fn control_points(&self) -> Vec<Vec<DVec3>> {
        self.points
            .chunks(self.count_v)
            .map(|row| row.iter().map(|&pw| dehomogenize(pw)).collect())
            .collect()
    }

    /// Weights of the control net `[i][j]`.
    pub fn weights(&self) -> Vec<Vec<f64>> {
        self.points
            .chunks(self.count_v)
            .map(|row| row.iter().map(|pw| pw.w).collect())
            .collect()
    }

    /// `true` if the weights are not all equal.
    pub fn is_rational(&self) -> bool {
        let w0 = self.points[0].w;
        self.points.iter().any(|pw| pw.w != w0)
    }

    /// Parametric domain `((u_min, u_max), (v_min, v_max))`.
    pub fn domain(&self) -> ((f64, f64), (f64, f64)) {
        (
            self.knots_u.domain(self.degree_u),
            self.knots_v.domain(self.degree_v),
        )
    }

    fn hp(&self, i: usize, j: usize) -> DVec4 {
        self.points[i * self.count_v + j]
    }

    fn clamp_params(&self, u: f64, v: f64) -> (f64, f64) {
        let ((ua, ub), (va, vb)) = self.domain();
        (u.clamp(ua, ub), v.clamp(va, vb))
    }

    /// Point at `(u, v)` (Algorithm A4.3, SurfacePoint for rational
    /// surfaces). Parameters are clamped to the domain.
    pub fn point(&self, u: f64, v: f64) -> DVec3 {
        let (u, v) = self.clamp_params(u, v);
        let (p, q) = (self.degree_u, self.degree_v);
        let uspan = self.knots_u.find_span(p, self.count_u - 1, u);
        let vspan = self.knots_v.find_span(q, self.count_v - 1, v);
        let nu = self.knots_u.basis_funs(uspan, u, p);
        let nv = self.knots_v.basis_funs(vspan, v, q);
        let mut sw = DVec4::ZERO;
        for (l, &bv) in nv.iter().enumerate() {
            let mut temp = DVec4::ZERO;
            for (k, &bu) in nu.iter().enumerate() {
                temp += self.hp(uspan - p + k, vspan - q + l) * bu;
            }
            sw += temp * bv;
        }
        dehomogenize(sw)
    }

    /// Partial derivatives of the homogeneous surface `S^w` up to total order
    /// `d` (Algorithm A3.6, SurfaceDerivsAlg1), evaluating the polynomial
    /// piece of the knot span cell `(uspan, vspan)`. Row `k` has `d + 1 - k`
    /// entries.
    ///
    /// When the weights active on the cell are all equal, the derivatives of
    /// the weight function are stored as exact zeros so that A4.4 reduces to
    /// the polynomial case without picking up rounding noise.
    fn homogeneous_derivatives_in_span(
        &self,
        u: f64,
        v: f64,
        d: usize,
        uspan: usize,
        vspan: usize,
    ) -> Vec<Vec<DVec4>> {
        let (p, q) = (self.degree_u, self.degree_v);
        let du = d.min(p);
        let dv = d.min(q);
        let nu = self.knots_u.ders_basis_funs(uspan, u, p, du);
        let nv = self.knots_v.ders_basis_funs(vspan, v, q, dv);
        let mut skl: Vec<Vec<DVec4>> = (0..=d).map(|k| vec![DVec4::ZERO; d + 1 - k]).collect();
        let mut temp = vec![DVec4::ZERO; q + 1];
        for k in 0..=du {
            for (s, t) in temp.iter_mut().enumerate() {
                *t = (0..=p).fold(DVec4::ZERO, |acc, r| {
                    acc + self.hp(uspan - p + r, vspan - q + s) * nu[k][r]
                });
            }
            let dd = (d - k).min(dv);
            for l in 0..=dd {
                skl[k][l] = temp
                    .iter()
                    .zip(&nv[l])
                    .fold(DVec4::ZERO, |acc, (&t, &n)| acc + t * n);
            }
        }
        let w0 = self.hp(uspan - p, vspan - q).w;
        let polynomial =
            (0..=p).all(|r| (0..=q).all(|s| self.hp(uspan - p + r, vspan - q + s).w == w0));
        if polynomial {
            for (k, row) in skl.iter_mut().enumerate() {
                for (l, value) in row.iter_mut().enumerate() {
                    if k + l > 0 {
                        value.w = 0.0;
                    }
                }
            }
        }
        skl
    }

    /// Partial derivatives up to total order `d` at `(u, v)`.
    ///
    /// Returns `skl` with `skl[k][l] = d^(k+l) S / du^k dv^l` for
    /// `k + l <= d` (row `k` has `d + 1 - k` entries), so `skl[0][0]` is the
    /// point, `skl[1][0]` is `S_u` and `skl[0][1]` is `S_v`. Computed with
    /// Algorithm A3.6 on the homogeneous surface followed by Algorithm A4.4
    /// (RatSurfaceDerivs). Parameters are clamped to the domain; on knot
    /// lines the derivatives of the span above are returned (below at the
    /// end of the domain).
    pub fn derivatives(&self, u: f64, v: f64, d: usize) -> Vec<Vec<DVec3>> {
        let (u, v) = self.clamp_params(u, v);
        let uspan = self.knots_u.find_span(self.degree_u, self.count_u - 1, u);
        let vspan = self.knots_v.find_span(self.degree_v, self.count_v - 1, v);
        rational_surface_derivatives(&self.homogeneous_derivatives_in_span(u, v, d, uspan, vspan))
    }

    /// Like [`NurbsSurface::derivatives`] but always evaluates the
    /// polynomial piece of the given cell (one-sided at knot lines).
    fn derivatives_in_cell(&self, u: f64, v: f64, d: usize, cell: &Cell) -> Vec<Vec<DVec3>> {
        rational_surface_derivatives(
            &self.homogeneous_derivatives_in_span(u, v, d, cell.uspan, cell.vspan),
        )
    }

    /// Typical magnitude of coordinates, used to scale tolerances.
    fn coordinate_scale(&self) -> f64 {
        let (lo, hi) = bounding_box(&self.points);
        (hi - lo).length().max(lo.abs().max(hi.abs()).max_element())
    }

    /// Unit normal `S_u x S_v / |S_u x S_v|` at `(u, v)`, or `None` where the
    /// surface is degenerate (a vanishing or negligible partial derivative,
    /// as at the pole of a surface of revolution, or parallel partials).
    pub fn normal(&self, u: f64, v: f64) -> Option<DVec3> {
        let skl = self.derivatives(u, v, 1);
        let (su, sv) = (skl[1][0], skl[0][1]);
        let ((ua, ub), (va, vb)) = self.domain();
        let scale = self.coordinate_scale();
        let (lu, lv) = (su.length(), sv.length());
        let min_u = 1e-12 * scale / (ub - ua);
        let min_v = 1e-12 * scale / (vb - va);
        let n = su.cross(sv);
        let usable = lu > min_u && lv > min_v && n.length() > 1e-12 * lu * lv;
        if usable {
            unit(n)
        } else {
            None
        }
    }

    /// Closest point on the surface to `p`, as `((u, v), point)`.
    ///
    /// Global search followed by Newton point inversion:
    ///
    /// 1. Every knot-span cell gets a lower bound on its distance to `p` from
    ///    the bounding box of its `(p + 1)(q + 1)` active control points
    ///    (strong convex hull property). Cells are visited nearest bound first
    ///    and skipped once the bound cannot beat the best point found so far.
    /// 2. A visited cell is sampled on a tensor-product grid whose rows and
    ///    columns are bisected where neighbouring samples are far apart, so
    ///    non-uniform rational parametrizations are covered evenly.
    /// 3. Each local minimum of the sampled distance, plus the few best
    ///    samples (best first, neighbours of equal distance on plateaus
    ///    skipped), is refined by two-dimensional Newton point inversion
    ///    (Eqs. 6.6-6.8) on that cell's polynomial piece, with the parameters
    ///    clamped to the cell (an active-set form of Eq. 6.8) and a
    ///    backtracking safeguard.
    ///
    /// The result is never farther from `p` than the best sample.
    pub fn closest_point(&self, p: DVec3) -> ((f64, f64), DVec3) {
        let (pu, pv) = (self.degree_u, self.degree_v);
        let iu = span_intervals(&self.knots_u, pu, self.count_u);
        let iv = span_intervals(&self.knots_v, pv, self.count_v);
        let mut cells: Vec<(Cell, f64)> = Vec::with_capacity(iu.len() * iv.len());
        for &(u0, u1, uspan) in &iu {
            for &(v0, v1, vspan) in &iv {
                let mut lo = DVec3::splat(f64::INFINITY);
                let mut hi = DVec3::splat(f64::NEG_INFINITY);
                for i in uspan - pu..=uspan {
                    for j in vspan - pv..=vspan {
                        let q = dehomogenize(self.hp(i, j));
                        lo = lo.min(q);
                        hi = hi.max(q);
                    }
                }
                let cell = Cell {
                    u: (u0, u1),
                    v: (v0, v1),
                    uspan,
                    vspan,
                };
                cells.push((cell, aabb_distance_squared(lo, hi, p)));
            }
        }
        cells.sort_by(|x, y| x.1.total_cmp(&y.1));
        let per_u = SAMPLES_PER_DIRECTION
            .div_ceil(iu.len())
            .clamp(4, 2 * pu + 6);
        let per_v = SAMPLES_PER_DIRECTION
            .div_ceil(iv.len())
            .clamp(4, 2 * pv + 6);
        let ((ua, _), (va, _)) = self.domain();
        let start = self.point(ua, va);
        let mut best = ((ua, va), start, start.distance_squared(p));
        for (cell, bound) in &cells {
            if *bound >= best.2 {
                break;
            }
            self.closest_in_cell(p, cell, per_u, per_v, &mut best);
        }
        (best.0, best.1)
    }

    /// Point of the polynomial piece of `cell` at `(u, v)`.
    fn point_in_cell(&self, u: f64, v: f64, cell: &Cell) -> DVec3 {
        let (p, q) = (self.degree_u, self.degree_v);
        let nu = self.knots_u.basis_funs(cell.uspan, u, p);
        let nv = self.knots_v.basis_funs(cell.vspan, v, q);
        let mut sw = DVec4::ZERO;
        for (l, &bv) in nv.iter().enumerate() {
            let mut temp = DVec4::ZERO;
            for (k, &bu) in nu.iter().enumerate() {
                temp += self.hp(cell.uspan - p + k, cell.vspan - q + l) * bu;
            }
            sw += temp * bv;
        }
        dehomogenize(sw)
    }

    /// Samples one cell and refines its local minima, updating
    /// `best = ((u, v), point, squared distance)`.
    fn closest_in_cell(
        &self,
        target: DVec3,
        cell: &Cell,
        per_u: usize,
        per_v: usize,
        best: &mut ((f64, f64), DVec3, f64),
    ) {
        let mut us = even_params(cell.u.0, cell.u.1, per_u);
        let mut vs = even_params(cell.v.0, cell.v.1, per_v);
        let eval_grid = |us: &[f64], vs: &[f64]| -> Vec<DVec3> {
            us.iter()
                .flat_map(|&u| vs.iter().map(move |&v| (u, v)))
                .map(|(u, v)| self.point_in_cell(u, v, cell))
                .collect()
        };
        let mut points = eval_grid(&us, &vs);
        // Bisect grid rows / columns whose neighbouring samples are far apart
        // in space (non-uniform rational parametrizations), keeping the grid
        // tensor-product; the limits keep the cost bounded.
        let (limit_u, limit_v) = (4 * per_u + 1, 4 * per_v + 1);
        let mut max_gap: Option<(f64, f64)> = None;
        for _ in 0..4 {
            let nv = vs.len();
            let gap_u = |i: usize| {
                (0..nv)
                    .map(|j| points[i * nv + j].distance(points[(i + 1) * nv + j]))
                    .fold(0.0, f64::max)
            };
            let gap_v = |j: usize| {
                (0..us.len())
                    .map(|i| points[i * nv + j].distance(points[i * nv + j + 1]))
                    .fold(0.0, f64::max)
            };
            let gaps_u: Vec<f64> = (0..us.len() - 1).map(gap_u).collect();
            let gaps_v: Vec<f64> = (0..nv - 1).map(gap_v).collect();
            let (lim_u, lim_v) = *max_gap.get_or_insert_with(|| {
                (
                    2.0 * gaps_u.iter().sum::<f64>() / gaps_u.len() as f64,
                    2.0 * gaps_v.iter().sum::<f64>() / gaps_v.len() as f64,
                )
            });
            let refine = |params: &[f64], gaps: &[f64], lim: f64, limit: usize| {
                let mut out = Vec::with_capacity(2 * params.len());
                for (k, &t) in params.iter().enumerate() {
                    out.push(t);
                    if k < gaps.len() && gaps[k] > lim && params.len() + out.len() - k < limit {
                        let tm = 0.5 * (t + params[k + 1]);
                        if tm > t && tm < params[k + 1] {
                            out.push(tm);
                        }
                    }
                }
                out
            };
            let new_us = refine(&us, &gaps_u, lim_u, limit_u);
            let new_vs = refine(&vs, &gaps_v, lim_v, limit_v);
            if new_us.len() == us.len() && new_vs.len() == vs.len() {
                break;
            }
            us = new_us;
            vs = new_vs;
            points = eval_grid(&us, &vs);
        }
        let (nu, nv) = (us.len(), vs.len());
        let grid: Vec<f64> = points.iter().map(|q| q.distance_squared(target)).collect();
        for (idx, &dist2) in grid.iter().enumerate() {
            if dist2 < best.2 {
                *best = ((us[idx / nv], vs[idx % nv]), points[idx], dist2);
            }
        }
        let d = |i: usize, j: usize| grid[i * nv + j];
        let neighbours = |i: usize, j: usize| {
            (i.saturating_sub(1)..(i + 2).min(nu))
                .flat_map(move |a| (j.saturating_sub(1)..(j + 2).min(nv)).map(move |b| (a, b)))
        };
        let mut minima: Vec<(usize, usize)> = (0..nu)
            .flat_map(|i| (0..nv).map(move |j| (i, j)))
            .filter(|&(i, j)| neighbours(i, j).all(|(a, b)| d(i, j) <= d(a, b)))
            .collect();
        // The few best samples are refined too: two minima of nearly equal
        // distance can share one sampled basin.
        let mut best_samples: Vec<usize> = (0..nu * nv).collect();
        best_samples.sort_by(|&a, &b| grid[a].total_cmp(&grid[b]));
        minima.extend(
            best_samples
                .into_iter()
                .take(EXTRA_BEST_SAMPLES)
                .map(|idx| (idx / nv, idx % nv)),
        );
        minima.sort_by(|x, y| d(x.0, x.1).total_cmp(&d(y.0, y.1)));
        minima.dedup();
        // Thin out plateaus: skip a minimum adjacent to an accepted one of
        // equal distance.
        let mut candidates: Vec<(usize, usize)> = Vec::new();
        for (i, j) in minima {
            if candidates.len() == MAX_CANDIDATES_PER_CELL {
                break;
            }
            let redundant = candidates
                .iter()
                .any(|&(a, b)| a.abs_diff(i) <= 1 && b.abs_diff(j) <= 1 && d(a, b) == d(i, j));
            if !redundant {
                candidates.push((i, j));
            }
        }
        for (i, j) in candidates {
            let (uv, q) = self.refine_closest(target, us[i], vs[j], cell);
            let dist2 = q.distance_squared(target);
            if dist2 < best.2 {
                *best = (uv, q, dist2);
            }
        }
    }

    /// Safeguarded two-dimensional Newton point inversion (Eqs. 6.6-6.8)
    /// from `(u0, v0)`, restricted to `cell`.
    fn refine_closest(&self, target: DVec3, u0: f64, v0: f64, cell: &Cell) -> ((f64, f64), DVec3) {
        let ((ua, ub), (va, vb)) = (cell.u, cell.v);
        let (mut u, mut v) = (u0.clamp(ua, ub), v0.clamp(va, vb));
        let mut skl = self.derivatives_in_cell(u, v, 2, cell);
        let mut r = skl[0][0] - target;
        let mut dist2 = r.length_squared();
        // Gradient of |r|^2 / 2 with components on active bounds removed.
        let projected_gradient = |u: f64, v: f64, g: DVec2| {
            let gu = if (u <= ua && g.x > 0.0) || (u >= ub && g.x < 0.0) {
                0.0
            } else {
                g.x
            };
            let gv = if (v <= va && g.y > 0.0) || (v >= vb && g.y < 0.0) {
                0.0
            } else {
                g.y
            };
            DVec2::new(gu, gv)
        };
        for _ in 0..MAX_NEWTON_ITERATIONS {
            if dist2 == 0.0 {
                break;
            }
            let (su, sv) = (skl[1][0], skl[0][1]);
            let g = DVec2::new(r.dot(su), r.dot(sv));
            let pg = projected_gradient(u, v, g);
            if pg == DVec2::ZERO {
                break;
            }
            let free_u = pg.x != 0.0;
            let free_v = pg.y != 0.0;
            // Newton matrix of Eq. 6.7 and its Gauss-Newton part.
            let guu = su.dot(su);
            let gvv = sv.dot(sv);
            let guv = su.dot(sv);
            let huu = guu + r.dot(skl[2][0]);
            let hvv = gvv + r.dot(skl[0][2]);
            let huv = guv + r.dot(skl[1][1]);
            let step = if free_u && free_v {
                solve_spd2(huu, huv, hvv, g).or_else(|| {
                    let lambda = 1e-10 * (guu + gvv) + f64::MIN_POSITIVE;
                    solve_spd2(guu + lambda, guv, gvv + lambda, g)
                })
            } else if free_u {
                let h = if huu > 0.0 { huu } else { guu };
                (h > 0.0).then(|| DVec2::new(-g.x / h, 0.0))
            } else {
                let h = if hvv > 0.0 { hvv } else { gvv };
                (h > 0.0).then(|| DVec2::new(0.0, -g.y / h))
            };
            let Some(step) = step else {
                break;
            };
            let scale = skl[0][0].abs().max(target.abs()).max_element();
            let slack = 1e-14 * (dist2 + dist2.sqrt() * scale);
            let mut alpha = 1.0;
            let mut accepted = None;
            for _ in 0..40 {
                let un = (u + alpha * step.x).clamp(ua, ub);
                let vn = (v + alpha * step.y).clamp(va, vb);
                if un == u && vn == v {
                    break;
                }
                let sn = self.derivatives_in_cell(un, vn, 2, cell);
                let rn = sn[0][0] - target;
                let d2n = rn.length_squared();
                let gn = projected_gradient(un, vn, DVec2::new(rn.dot(sn[1][0]), rn.dot(sn[0][1])));
                if d2n < dist2 || (d2n <= dist2 + slack && gn.length() < pg.length()) {
                    accepted = Some((un, vn, sn, rn, d2n));
                    break;
                }
                alpha *= 0.5;
            }
            let Some((un, vn, sn, rn, d2n)) = accepted else {
                break;
            };
            let moved_u = (un - u).abs();
            let moved_v = (vn - v).abs();
            u = un;
            v = vn;
            skl = sn;
            r = rn;
            dist2 = d2n;
            if moved_u <= 1e-15 * (ub - ua).max(u.abs())
                && moved_v <= 1e-15 * (vb - va).max(v.abs())
            {
                break;
            }
        }
        ((u, v), skl[0][0])
    }

    /// Surface of revolution (Algorithm A8.1, MakeRevolvedSurf).
    ///
    /// Rotates `profile` by `angle` radians about the axis through
    /// `axis_origin` along `axis_dir` (right-hand rule; a negative angle turns
    /// the other way). The `u` direction runs around the axis (degree 2,
    /// knots on `[0, 1]` as in [`NurbsCurve::circular_arc`]) and the `v`
    /// direction follows the profile (its degree and knots), so
    /// `S(0, v) = profile(v)`. `|angle|` must be in `(0, 2 pi]`.
    pub fn revolution(
        profile: &NurbsCurve,
        axis_origin: DVec3,
        axis_dir: DVec3,
        angle: f64,
    ) -> Result<Self, NurbsError> {
        if !axis_origin.is_finite() {
            return Err(NurbsError::InvalidArgument("axis origin is not finite"));
        }
        let mut axis = unit(axis_dir).ok_or(NurbsError::InvalidArgument(
            "revolution axis direction is degenerate",
        ))?;
        if !angle.is_finite() {
            return Err(NurbsError::InvalidArgument(
                "revolution angle is not finite",
            ));
        }
        let mut sweep = angle;
        if sweep < 0.0 {
            sweep = -sweep;
            axis = -axis;
        }
        let sweep = normalize_sweep(sweep)?;
        let template = arc_template(sweep);
        let profile_points = profile.homogeneous_points();
        // Per profile point: axis foot point, radial and tangential directions,
        // radius, and the Cartesian point itself.
        let frames: Vec<(DVec3, DVec3, DVec3, f64, DVec3, f64)> = profile_points
            .iter()
            .map(|&pw| {
                let p = dehomogenize(pw);
                let foot = axis_origin + axis * (p - axis_origin).dot(axis);
                let radial = p - foot;
                let radius = radial.length();
                let (x, y) = match unit(radial) {
                    Some(x) => (x, axis.cross(x)),
                    None => (DVec3::ZERO, DVec3::ZERO),
                };
                (foot, x, y, radius, p, pw.w)
            })
            .collect();
        let count_u = template.frames.len();
        let count_v = profile_points.len();
        let mut points = Vec::with_capacity(count_u * count_v);
        for &(offset, scale, wa) in &template.frames {
            let (s, c) = offset.sin_cos();
            for &(foot, x, y, radius, p, w) in &frames {
                let q = if offset == 0.0 {
                    p
                } else {
                    foot + (x * c + y * s) * (radius * scale)
                };
                points.push(homogenize(q, w * wa));
            }
        }
        if sweep == TAU {
            let (first, rest) = points.split_at_mut(count_v);
            let last_row = rest.len() - count_v;
            rest[last_row..].copy_from_slice(first);
        }
        Self::from_homogeneous_flat(
            2,
            profile.degree(),
            template.knots,
            profile.knots().as_slice().to_vec(),
            count_u,
            count_v,
            points,
        )
    }

    /// Extrusion of `curve` along `dir`: `S(u, v) = curve(u) + v * dir` with
    /// `u` on the curve domain and `v` in `[0, 1]` (degree 1).
    ///
    /// Fails if `dir` is zero or not finite.
    pub fn extrusion(curve: &NurbsCurve, dir: DVec3) -> Result<Self, NurbsError> {
        if !dir.is_finite() || dir == DVec3::ZERO {
            return Err(NurbsError::InvalidArgument(
                "extrusion direction must be finite and non-zero",
            ));
        }
        let mut points = Vec::with_capacity(curve.homogeneous_points().len() * 2);
        for &pw in curve.homogeneous_points() {
            points.push(pw);
            points.push(pw + (dir * pw.w).extend(0.0));
        }
        Self::from_homogeneous_flat(
            curve.degree(),
            1,
            curve.knots().as_slice().to_vec(),
            vec![0.0, 0.0, 1.0, 1.0],
            curve.homogeneous_points().len(),
            2,
            points,
        )
    }

    /// Isoparametric curve `v -> S(u, v)` at fixed `u` (clamped to the
    /// domain; NaN is treated as the domain start). It has the `v` degree and
    /// knots and is exact.
    pub fn iso_curve_u(&self, u: f64) -> NurbsCurve {
        let ((ua, ub), _) = self.domain();
        let u = if u.is_nan() { ua } else { u.clamp(ua, ub) };
        let p = self.degree_u;
        let span = self.knots_u.find_span(p, self.count_u - 1, u);
        let basis = self.knots_u.basis_funs(span, u, p);
        let points = (0..self.count_v)
            .map(|j| {
                basis.iter().enumerate().fold(DVec4::ZERO, |acc, (k, &b)| {
                    acc + self.hp(span - p + k, j) * b
                })
            })
            .collect();
        NurbsCurve::from_parts(self.degree_v, self.knots_v.as_slice().to_vec(), points)
    }

    /// Isoparametric curve `u -> S(u, v)` at fixed `v` (clamped to the
    /// domain; NaN is treated as the domain start). It has the `u` degree and
    /// knots and is exact.
    pub fn iso_curve_v(&self, v: f64) -> NurbsCurve {
        let (_, (va, vb)) = self.domain();
        let v = if v.is_nan() { va } else { v.clamp(va, vb) };
        let q = self.degree_v;
        let span = self.knots_v.find_span(q, self.count_v - 1, v);
        let basis = self.knots_v.basis_funs(span, v, q);
        let points = (0..self.count_u)
            .map(|i| {
                basis.iter().enumerate().fold(DVec4::ZERO, |acc, (l, &b)| {
                    acc + self.hp(i, span - q + l) * b
                })
            })
            .collect();
        NurbsCurve::from_parts(self.degree_u, self.knots_u.as_slice().to_vec(), points)
    }

    /// Applies the affine transformation `m` (its projective row is ignored).
    /// Weights and parametrization are unchanged, so the result is exact.
    /// `m` must be finite (asserted in debug builds).
    pub fn transformed(&self, m: &DMat4) -> NurbsSurface {
        debug_assert!(m.is_finite(), "transform matrix must be finite");
        let mut out = self.clone();
        for pw in &mut out.points {
            *pw = transform_homogeneous(m, *pw);
        }
        out
    }

    /// Axis-aligned bounding box `(min, max)` of the control net, which
    /// contains the surface by the convex hull property.
    pub fn bounding_box(&self) -> (DVec3, DVec3) {
        bounding_box(&self.points)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bilinear_patch() {
        let s = NurbsSurface::bilinear(DVec3::ZERO, DVec3::X, DVec3::Y, DVec3::new(1.0, 1.0, 1.0))
            .unwrap();
        assert_eq!(s.point(0.0, 0.0), DVec3::ZERO);
        assert_eq!(s.point(1.0, 0.0), DVec3::X);
        assert_eq!(s.point(0.0, 1.0), DVec3::Y);
        assert_eq!(s.point(1.0, 1.0), DVec3::ONE);
        assert!((s.point(0.5, 0.5) - DVec3::new(0.5, 0.5, 0.25)).length() < 1e-15);
        let n = s.normal(0.0, 0.0).unwrap();
        assert!((n - DVec3::Z).length() < 1e-15);
    }

    #[test]
    fn sphere_pole_normal_is_none() {
        let profile = NurbsCurve::circular_arc(
            DVec3::ZERO,
            DVec3::Z,
            DVec3::X,
            2.0,
            0.0,
            std::f64::consts::PI,
        )
        .unwrap();
        let sphere = NurbsSurface::revolution(&profile, DVec3::ZERO, DVec3::Z, TAU).unwrap();
        assert!(sphere.normal(0.3, 0.0).is_none());
        assert!(sphere.normal(0.3, 1.0).is_none());
        let n = sphere.normal(0.3, 0.4).unwrap();
        let p = sphere.point(0.3, 0.4);
        assert!((p.length() - 2.0).abs() < 1e-14);
        assert!((n.dot(p / 2.0).abs() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn iso_curves_match_surface() {
        let s = NurbsSurface::new(
            2,
            1,
            vec![0., 0., 0., 1., 1., 1.],
            vec![0., 0., 2., 2.],
            vec![
                vec![DVec3::ZERO, DVec3::new(0., 2., 1.)],
                vec![DVec3::new(1., 0., 1.), DVec3::new(1., 2., 0.)],
                vec![DVec3::new(2., 0., 0.), DVec3::new(2., 2., 2.)],
            ],
            Some(vec![vec![1., 2.], vec![0.5, 1.], vec![1., 3.]]),
        )
        .unwrap();
        let cu = s.iso_curve_u(0.3);
        let cv = s.iso_curve_v(1.2);
        for i in 0..=10 {
            let t = i as f64 / 10.0;
            assert!((cu.point(2.0 * t) - s.point(0.3, 2.0 * t)).length() < 1e-14);
            assert!((cv.point(t) - s.point(t, 1.2)).length() < 1e-14);
        }
    }
}
