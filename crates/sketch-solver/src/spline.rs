//! Non-rational B-spline machinery for sketch splines (Piegl & Tiller,
//! *The NURBS Book*, 2nd ed.): knot span search (A2.1), basis functions and
//! their derivatives (A2.2, A2.3), clamped uniform knots, global
//! interpolation with knots by averaging (eq. 9.8, A9.1), knot insertion
//! (A5.1), and smooth closed (periodic) splines as exact clamped
//! equivalents (ADR 0016).
//!
//! A sketch spline is a fixed linear map from its referenced point entities
//! to its curve: `C(t) = Σ_k w_k(t) P_k`. For a control-point spline the
//! weights are the basis functions; for an interpolated spline they are the
//! basis functions composed with the inverse of the collocation matrix (the
//! hidden control points are `Q = A⁻¹ P`, design note 0005 §3). The solver
//! therefore needs no hidden control-point unknowns: the interpolation
//! equations are linear and nonsingular, and are eliminated exactly.

use glam::DVec2;
use openshape_linalg::{CompleteOrthogonal, Matrix};

/// A spline as a linear map from its referenced points to its curve.
#[derive(Clone, Debug)]
pub(crate) struct SplineMap {
    /// Degree `p`.
    pub degree: usize,
    /// Clamped knot vector on `[0, 1]`.
    pub knots: Vec<f64>,
    /// `inv[(j, k)]`: weight of referenced point `k` in control point `j`
    /// (an interpolated or a periodic spline). `None`: the referenced
    /// points are the control points.
    pub inv: Option<Matrix>,
    /// Number of control points of the clamped curve.
    pub count: usize,
    /// Number of referenced points (`count` unless periodic, where the
    /// clamped form has more control points than the sketch references).
    pub refs: usize,
    /// Closed smoothly (ADR 0016): the curve is periodic, its parameter
    /// wraps, and the last referenced point is the first again.
    pub periodic: bool,
    /// Interior knots of multiplicity ≥ degree, ascending: where the first
    /// derivative may jump (corners of a degree-1 spline).
    pub kinks: Vec<f64>,
}

/// Reasons a spline definition is unusable.
#[derive(Debug, Clone, PartialEq)]
pub(crate) enum SplineError {
    /// Wrong number of points for the degree, or bad degree.
    Shape(String),
    /// Parameters not strictly increasing from 0 to 1, or not finite.
    Params(String),
}

impl SplineMap {
    /// A clamped spline with uniform interior knots on `[0, 1]`.
    pub fn control(count: usize, degree: usize) -> Result<Self, SplineError> {
        if !(1..=5).contains(&degree) {
            return Err(SplineError::Shape(format!(
                "degree {degree} is outside 1..=5"
            )));
        }
        if count < degree + 1 {
            return Err(SplineError::Shape(format!(
                "{count} control points are fewer than degree + 1 = {}",
                degree + 1
            )));
        }
        let knots = clamped_uniform_knots(count, degree);
        Ok(SplineMap {
            degree,
            kinks: kinks(&knots, degree),
            knots,
            inv: None,
            count,
            refs: count,
            periodic: false,
        })
    }

    /// The interpolating spline through `count = params.len()` points at
    /// `params`: degree `min(3, count − 1)`, knots by averaging.
    pub fn interpolated(params: &[f64]) -> Result<Self, SplineError> {
        let count = params.len();
        if count < 2 {
            return Err(SplineError::Shape(format!(
                "an interpolated spline needs at least 2 points, got {count}"
            )));
        }
        check_params(params)?;
        let degree = 3.min(count - 1);
        let knots = averaged_knots(params, degree);
        // Collocation matrix A[i][j] = N_j(t_i).
        let mut a = Matrix::zeros(count, count);
        for (i, &t) in params.iter().enumerate() {
            let span = find_span(count - 1, degree, t, &knots);
            let n = basis_funs(span, t, degree, &knots);
            for (r, v) in n.iter().enumerate() {
                a[(i, span - degree + r)] = *v;
            }
        }
        let cod = CompleteOrthogonal::new(&a, 1e-13);
        if cod.rank() < count {
            return Err(SplineError::Params(
                "the interpolation system is singular".to_string(),
            ));
        }
        let mut inv = Matrix::zeros(count, count);
        let mut e = vec![0.0; count];
        for k in 0..count {
            e.iter_mut().for_each(|v| *v = 0.0);
            e[k] = 1.0;
            let x = cod.solve(&e);
            inv.col_mut(k).copy_from_slice(&x);
        }
        if !inv.is_finite() {
            return Err(SplineError::Params(
                "the interpolation system is singular".to_string(),
            ));
        }
        // Clamped knots make the first and last collocation rows unit rows,
        // so the end control points are the end through-points: exactly,
        // not to rounding (the curve meets what its ends are joined to).
        for k in 0..count {
            inv[(0, k)] = if k == 0 { 1.0 } else { 0.0 };
            inv[(count - 1, k)] = if k == count - 1 { 1.0 } else { 0.0 };
        }
        Ok(SplineMap {
            degree,
            kinks: kinks(&knots, degree),
            knots,
            inv: Some(inv),
            count,
            refs: count,
            periodic: false,
        })
    }

    /// A smooth closed control-point spline (ADR 0016) of degree `degree`
    /// on `refs` referenced points, the last being the first again: the
    /// uniform periodic B-spline on the wrapped polygon `P₀ … P_{n−1}, P₀ …
    /// P_{p−1}` (`n = refs − 1`; Piegl & Tiller §12.2), as the clamped
    /// spline on `n + p` control points that it equals on one period. The
    /// domain's end knots are raised to multiplicity `p` by Boehm insertion
    /// (P&T A5.1), which makes the curve's ends control points; the knots
    /// outside the domain then no longer act on it and are clamped.
    pub fn control_periodic(refs: usize, degree: usize) -> Result<Self, SplineError> {
        if !(2..=5).contains(&degree) {
            return Err(SplineError::Shape(format!(
                "a smooth closed spline has degree 2 to 5, not {degree}"
            )));
        }
        let n = refs.saturating_sub(1);
        if n < degree + 1 {
            return Err(SplineError::Shape(format!(
                "a smooth closed spline of degree {degree} needs at least {} distinct control points, got {n}",
                degree + 1
            )));
        }
        let p = degree;
        // The wrapped polygon, each point as its weights over the
        // referenced points, on the uniform knots 0, 1, …, n + 2p.
        let mut rows: Vec<Vec<f64>> = (0..n + p)
            .map(|j| {
                let mut w = vec![0.0; refs];
                w[j % n] = 1.0;
                w
            })
            .collect();
        let mut knots: Vec<f64> = (0..=n + 2 * p).map(|i| i as f64).collect();
        for u in [p as f64, (n + p) as f64] {
            for _ in 1..p {
                insert_knot(&mut knots, &mut rows, p, u);
            }
        }
        // C(p) is control point p − 1 and C(n + p) is control point
        // 2p + n − 2 (a knot of multiplicity p at index s ends the curve
        // on control point s − p).
        let mut rows: Vec<Vec<f64>> = rows.drain(p - 1..=2 * p + n - 2).collect();
        // The same combination of the same point at both ends: the curve
        // closes exactly, not to rounding.
        let last = rows.len() - 1;
        rows[last] = rows[0].clone();
        let controls = n + p;
        let mut inv = Matrix::zeros(controls, refs);
        for (j, row) in rows.iter().enumerate() {
            for (k, w) in row.iter().enumerate() {
                inv[(j, k)] = *w;
            }
        }
        let knots = clamped_uniform_knots(controls, p);
        match_end_derivatives(&mut inv, &knots, p);
        Ok(SplineMap {
            degree: p,
            kinks: kinks(&knots, p),
            knots,
            inv: Some(inv),
            count: controls,
            refs,
            periodic: true,
        })
    }

    /// A smooth closed through-point spline (ADR 0016) through `n =
    /// params.len() − 1` distinct points at `params` (first 0, last 1; the
    /// last referenced point is the first again): the periodic cubic
    /// spline with knots at the parameters (de Boor, *A Practical Guide to
    /// Splines*, rev. ed. 2001), as the clamped cubic on the knots `0⁴,
    /// t₁ … t_{n−1}, 1⁴`. Its `n + 3` control points satisfy the `n + 1`
    /// interpolation conditions and equal first and second derivatives at
    /// 0 and 1.
    pub fn interpolated_periodic(params: &[f64]) -> Result<Self, SplineError> {
        let refs = params.len();
        if refs < 4 {
            return Err(SplineError::Shape(format!(
                "a smooth closed spline needs at least 3 distinct points, got {}",
                refs.saturating_sub(1)
            )));
        }
        check_params(params)?;
        let p = 3;
        let n = refs - 1;
        let controls = n + p;
        let mut knots = vec![0.0; p + 1];
        knots.extend_from_slice(&params[1..n]);
        knots.extend(std::iter::repeat_n(1.0, p + 1));
        let mut a = Matrix::zeros(controls, controls);
        let row = |a: &mut Matrix, i: usize, t: f64, order: usize, sign: f64| {
            let span = find_span(controls - 1, p, t, &knots);
            let d = ders_basis_funs(span, t, p, order, &knots);
            for (r, v) in d[order].iter().enumerate() {
                a[(i, span - p + r)] += sign * v;
            }
        };
        for (i, &t) in params.iter().enumerate() {
            row(&mut a, i, t, 0, 1.0);
        }
        for order in 1..=2 {
            let i = n + order;
            row(&mut a, i, 0.0, order, 1.0);
            row(&mut a, i, 1.0, order, -1.0);
        }
        let cod = CompleteOrthogonal::new(&a, 1e-13);
        if cod.rank() < controls {
            return Err(SplineError::Params(
                "the interpolation system is singular".to_string(),
            ));
        }
        // Right-hand sides: interpolation row i takes referenced point i,
        // the closing row (i = n) the first point; the derivative rows are
        // zero.
        let mut inv = Matrix::zeros(controls, refs);
        let mut e = vec![0.0; controls];
        for k in 0..n {
            e.iter_mut().for_each(|v| *v = 0.0);
            e[k] = 1.0;
            if k == 0 {
                e[n] = 1.0;
            }
            let x = cod.solve(&e);
            inv.col_mut(k).copy_from_slice(&x);
        }
        if !inv.is_finite() {
            return Err(SplineError::Params(
                "the interpolation system is singular".to_string(),
            ));
        }
        // The ends are the first point, exactly.
        for k in 0..refs {
            let w = if k == 0 { 1.0 } else { 0.0 };
            inv[(0, k)] = w;
            inv[(controls - 1, k)] = w;
        }
        match_end_derivatives(&mut inv, &knots, p);
        Ok(SplineMap {
            degree: p,
            kinks: kinks(&knots, p),
            knots,
            inv: Some(inv),
            count: controls,
            refs,
            periodic: true,
        })
    }

    /// Weights `w_k(t)` and `w'_k(t)` of the referenced points (dense).
    /// `t` is clamped to `[0, 1]`.
    pub fn weights(&self, t: f64) -> (Vec<f64>, Vec<f64>) {
        let (w, d, _) = self.weights_ders(t, 1);
        (w, d)
    }

    /// Weights and their first and (if `ders ≥ 2`) second derivatives.
    /// `t` is clamped to `[0, 1]`, or wrapped into it when periodic.
    pub fn weights_ders(&self, t: f64, ders: usize) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
        let t = if self.periodic {
            t.rem_euclid(1.0)
        } else {
            t.clamp(0.0, 1.0)
        };
        let span = find_span(self.count - 1, self.degree, t, &self.knots);
        self.weights_in_span(span, t, ders)
    }

    /// At a kink `t`, the weight derivatives `w'_k` of the piece to the left
    /// of `t` (the right-hand piece's are those of [`Self::weights`]).
    /// `None` if `t` is not a kink.
    pub fn left_derivative(&self, t: f64) -> Option<Vec<f64>> {
        if !self.kinks.contains(&t) {
            return None;
        }
        let right = find_span(self.count - 1, self.degree, t, &self.knots);
        let mut left = right;
        while self.knots[left] >= t {
            left -= 1;
        }
        Some(self.weights_in_span(left, t, 1).1)
    }

    /// Weights and derivatives from the polynomial piece of knot span
    /// `span`, evaluated at `t`.
    fn weights_in_span(&self, span: usize, t: f64, ders: usize) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
        let p = self.degree;
        let d = ders_basis_funs(span, t, p, ders.min(p), &self.knots);
        let mut w = vec![0.0; self.refs];
        let mut w1 = vec![0.0; self.refs];
        let mut w2 = vec![0.0; self.refs];
        for r in 0..=p {
            let j = span - p + r;
            match &self.inv {
                None => {
                    w[j] += d[0][r];
                    if d.len() > 1 {
                        w1[j] += d[1][r];
                    }
                    if d.len() > 2 {
                        w2[j] += d[2][r];
                    }
                }
                Some(inv) => {
                    for k in 0..self.refs {
                        let c = inv[(j, k)];
                        w[k] += d[0][r] * c;
                        if d.len() > 1 {
                            w1[k] += d[1][r] * c;
                        }
                        if d.len() > 2 {
                            w2[k] += d[2][r] * c;
                        }
                    }
                }
            }
        }
        (w, w1, w2)
    }

    /// Weights of the referenced points in control point `j`.
    pub fn control_weights(&self, j: usize) -> Vec<f64> {
        match &self.inv {
            None => {
                let mut w = vec![0.0; self.count];
                w[j] = 1.0;
                w
            }
            Some(inv) => (0..self.refs).map(|k| inv[(j, k)]).collect(),
        }
    }

    /// The control points for the given referenced point positions.
    pub fn control_points(&self, pts: &[DVec2]) -> Vec<DVec2> {
        (0..self.count)
            .map(|j| {
                let w = self.control_weights(j);
                combine(&w, pts)
            })
            .collect()
    }
}

/// `Σ w_k p_k`, summed in index order.
pub(crate) fn combine(w: &[f64], pts: &[DVec2]) -> DVec2 {
    let mut x = 0.0;
    let mut y = 0.0;
    for (wk, pk) in w.iter().zip(pts) {
        x += wk * pk.x;
        y += wk * pk.y;
    }
    DVec2::new(x, y)
}

/// Makes a closed clamped spline's first derivative at its end equal the
/// one at its start exactly, not just to rounding (ADR 0016): the
/// second-to-last control row (`inv` maps referenced points to control
/// points, one row per control point; the last row is the first already)
/// is set from the first two, by `C′(0) = p (P₁ − P₀)/(u_{p+1} − u₁)` and
/// `C′(1) = p (P_m − P_{m−1})/(u_{m+p} − u_m)` (Piegl & Tiller eq. 3.7).
/// The periodic constructions make the two equal mathematically, so this
/// moves the row by rounding only; a display can then tell the seam is
/// smooth by comparing the derivatives exactly.
fn match_end_derivatives(inv: &mut Matrix, knots: &[f64], p: usize) {
    let m = inv.rows() - 1;
    let r = (knots[m + p] - knots[m]) / (knots[p + 1] - knots[1]);
    for k in 0..inv.cols() {
        inv[(m - 1, k)] = inv[(m, k)] - r * (inv[(1, k)] - inv[(0, k)]);
    }
}

/// Inserts knot `u` once into the B-spline of degree `p` on `knots` whose
/// control points are `rows` (each a vector of weights), in place (Boehm
/// insertion, Piegl & Tiller A5.1 with r = 1). `u` must lie inside the
/// knot range with multiplicity below `p`.
fn insert_knot(knots: &mut Vec<f64>, rows: &mut Vec<Vec<f64>>, p: usize, u: f64) {
    // The last span index k with knots[k] ≤ u, and u's multiplicity s.
    let k = knots.iter().rposition(|&x| x <= u).unwrap_or(0);
    let s = knots.iter().filter(|&&x| x == u).count();
    let mut out = Vec::with_capacity(rows.len() + 1);
    for i in 0..=rows.len() {
        if i + p <= k {
            out.push(rows[i].clone());
        } else if i + s <= k {
            let alpha = (u - knots[i]) / (knots[i + p] - knots[i]);
            let row = rows[i]
                .iter()
                .zip(&rows[i - 1])
                .map(|(a, b)| alpha * a + (1.0 - alpha) * b)
                .collect();
            out.push(row);
        } else {
            out.push(rows[i - 1].clone());
        }
    }
    knots.insert(k + 1, u);
    *rows = out;
}

/// Checks that interpolation parameters are finite, strictly increasing,
/// start at 0 and end at 1.
pub(crate) fn check_params(params: &[f64]) -> Result<(), SplineError> {
    if params.iter().any(|t| !t.is_finite()) {
        return Err(SplineError::Params("parameters must be finite".into()));
    }
    let (first, last) = (params[0], params[params.len() - 1]);
    if first != 0.0 || last != 1.0 {
        return Err(SplineError::Params(format!(
            "parameters must run from 0 to 1 (got {first} to {last})"
        )));
    }
    if params.windows(2).any(|w| w[1] <= w[0]) {
        return Err(SplineError::Params(
            "parameters must be strictly increasing".into(),
        ));
    }
    Ok(())
}

/// Interior knot values of multiplicity ≥ `p` (ascending, distinct).
fn kinks(knots: &[f64], p: usize) -> Vec<f64> {
    let mut out = Vec::new();
    let mut i = 0;
    while i < knots.len() {
        let mut j = i;
        while j + 1 < knots.len() && knots[j + 1] == knots[i] {
            j += 1;
        }
        let u = knots[i];
        if u > 0.0 && u < 1.0 && j + 1 - i >= p {
            out.push(u);
        }
        i = j + 1;
    }
    out
}

/// Clamped knots on `[0, 1]` with uniform interior knots for `count`
/// control points of degree `p` (the knots of a
/// [`Geometry::ControlSpline`](crate::Geometry::ControlSpline)); needs
/// `count > p`.
pub fn clamped_uniform_knots(count: usize, p: usize) -> Vec<f64> {
    let n = count - 1;
    let segments = n - p + 1;
    let mut u = Vec::with_capacity(n + p + 2);
    u.extend(std::iter::repeat_n(0.0, p + 1));
    for j in 1..segments {
        u.push(j as f64 / segments as f64);
    }
    u.extend(std::iter::repeat_n(1.0, p + 1));
    u
}

/// Knots by averaging (P&T eq. 9.8) for parameters `t_0 = 0 < … < t_n = 1`.
pub(crate) fn averaged_knots(t: &[f64], p: usize) -> Vec<f64> {
    let n = t.len() - 1;
    let mut u = Vec::with_capacity(n + p + 2);
    u.extend(std::iter::repeat_n(0.0, p + 1));
    for j in 1..=n - p {
        let mut s = 0.0;
        for ti in &t[j..j + p] {
            s += ti;
        }
        u.push(s / p as f64);
    }
    u.extend(std::iter::repeat_n(1.0, p + 1));
    u
}

/// Knot span index (P&T A2.1) for control points `0..=n`, degree `p`.
pub(crate) fn find_span(n: usize, p: usize, u: f64, knots: &[f64]) -> usize {
    if u >= knots[n + 1] {
        return n;
    }
    if u <= knots[p] {
        return p;
    }
    let (mut low, mut high) = (p, n + 1);
    let mut mid = (low + high) / 2;
    while u < knots[mid] || u >= knots[mid + 1] {
        if u < knots[mid] {
            high = mid;
        } else {
            low = mid;
        }
        mid = (low + high) / 2;
    }
    mid
}

/// Nonzero basis functions `N_{span−p..=span, p}(u)` (P&T A2.2).
pub(crate) fn basis_funs(span: usize, u: f64, p: usize, knots: &[f64]) -> Vec<f64> {
    ders_basis_funs(span, u, p, 0, knots).swap_remove(0)
}

/// Basis functions and their derivatives up to order `nd ≤ p`
/// (P&T A2.3): `ders[k][r]` is the k-th derivative of `N_{span−p+r, p}`.
pub(crate) fn ders_basis_funs(
    span: usize,
    u: f64,
    p: usize,
    nd: usize,
    knots: &[f64],
) -> Vec<Vec<f64>> {
    let mut ndu = vec![vec![0.0; p + 1]; p + 1];
    let mut left = vec![0.0; p + 1];
    let mut right = vec![0.0; p + 1];
    ndu[0][0] = 1.0;
    for j in 1..=p {
        left[j] = u - knots[span + 1 - j];
        right[j] = knots[span + j] - u;
        let mut saved = 0.0;
        for r in 0..j {
            ndu[j][r] = right[r + 1] + left[j - r];
            let temp = ndu[r][j - 1] / ndu[j][r];
            ndu[r][j] = saved + right[r + 1] * temp;
            saved = left[j - r] * temp;
        }
        ndu[j][j] = saved;
    }
    let mut ders = vec![vec![0.0; p + 1]; nd + 1];
    for j in 0..=p {
        ders[0][j] = ndu[j][p];
    }
    let mut a = vec![vec![0.0; p + 1]; 2];
    for r in 0..=p {
        let (mut s1, mut s2) = (0usize, 1usize);
        a[0][0] = 1.0;
        for k in 1..=nd {
            let mut d = 0.0;
            let rk = r as isize - k as isize;
            let pk = p - k;
            if r >= k {
                a[s2][0] = a[s1][0] / ndu[pk + 1][rk as usize];
                d = a[s2][0] * ndu[rk as usize][pk];
            }
            let j1 = if rk >= -1 { 1 } else { (-rk) as usize };
            let j2 = if (r as isize - 1) <= pk as isize {
                k - 1
            } else {
                p - r
            };
            for j in j1..=j2 {
                let idx = (rk + j as isize) as usize;
                a[s2][j] = (a[s1][j] - a[s1][j - 1]) / ndu[pk + 1][idx];
                d += a[s2][j] * ndu[idx][pk];
            }
            if r <= pk {
                a[s2][k] = -a[s1][k - 1] / ndu[pk + 1][r];
                d += a[s2][k] * ndu[r][pk];
            }
            ders[k][r] = d;
            std::mem::swap(&mut s1, &mut s2);
        }
    }
    let mut factor = p as f64;
    for (k, row) in ders.iter_mut().enumerate().skip(1) {
        for v in row.iter_mut() {
            *v *= factor;
        }
        factor *= (p - k) as f64;
    }
    ders
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn partition_of_unity_and_derivative_sum() {
        let m = SplineMap::control(7, 3).unwrap();
        for i in 0..=50 {
            let t = i as f64 / 50.0;
            let (w, d) = m.weights(t);
            let s: f64 = w.iter().sum();
            let ds: f64 = d.iter().sum();
            assert!((s - 1.0).abs() < 1e-14 && ds.abs() < 1e-12, "{t}");
        }
    }

    #[test]
    fn derivatives_match_finite_differences() {
        for (count, p) in [(4usize, 3usize), (7, 3), (6, 5), (3, 1), (5, 2)] {
            let m = SplineMap::control(count, p).unwrap();
            for i in 1..40 {
                let t = i as f64 / 40.0 + 0.003;
                let h = 1e-6;
                let (wp, _) = m.weights(t + h);
                let (wm, _) = m.weights(t - h);
                let (_, d, d2) = m.weights_ders(t, 2);
                let (_, dp, _) = m.weights_ders(t + h, 1);
                let (_, dm, _) = m.weights_ders(t - h, 1);
                for k in 0..count {
                    let fd = (wp[k] - wm[k]) / (2.0 * h);
                    assert!((fd - d[k]).abs() < 1e-6, "p={p} k={k} t={t}");
                    if p >= 2 {
                        let fd2 = (dp[k] - dm[k]) / (2.0 * h);
                        assert!((fd2 - d2[k]).abs() < 1e-4 * (1.0 + d2[k].abs()));
                    }
                }
            }
        }
    }

    #[test]
    fn kinks_and_one_sided_derivatives() {
        // Degree 1 through 4 points: corners at 1/3 and 2/3.
        let m = SplineMap::control(4, 1).unwrap();
        assert_eq!(m.kinks, vec![1.0 / 3.0, 2.0 / 3.0]);
        let t = m.kinks[0];
        let left = m.left_derivative(t).unwrap();
        let (_, right) = m.weights(t);
        let (_, before) = m.weights(t - 1e-3);
        let (_, after) = m.weights(t + 1e-3);
        let close = |a: &[f64], b: &[f64]| a.iter().zip(b).all(|(x, y)| (x - y).abs() < 1e-9);
        assert!(close(&left, &before));
        assert!(close(&right, &after));
        assert!(!close(&left, &right));
        assert!(m.left_derivative(0.5).is_none());
        // Smooth splines have none.
        assert!(SplineMap::control(7, 3).unwrap().kinks.is_empty());
        assert!(SplineMap::control(5, 2).unwrap().kinks.is_empty());
        assert!(SplineMap::interpolated(&[0.0, 0.3, 0.6, 1.0])
            .unwrap()
            .kinks
            .is_empty());
    }

    #[test]
    fn interpolation_passes_through_points() {
        let params = [0.0, 0.2, 0.45, 0.8, 1.0];
        let pts = [
            DVec2::new(0.0, 0.0),
            DVec2::new(1.0, 2.0),
            DVec2::new(2.0, -1.0),
            DVec2::new(4.0, 0.5),
            DVec2::new(5.0, 3.0),
        ];
        let m = SplineMap::interpolated(&params).unwrap();
        assert_eq!(m.degree, 3);
        for (t, p) in params.iter().zip(&pts) {
            let (w, _) = m.weights(*t);
            let c = combine(&w, &pts);
            assert!((c - *p).length() < 1e-13);
        }
        // Two points: the straight segment.
        let m2 = SplineMap::interpolated(&[0.0, 1.0]).unwrap();
        assert_eq!(m2.degree, 1);
        let (w, _) = m2.weights(0.25);
        assert!((w[0] - 0.75).abs() < 1e-15 && (w[1] - 0.25).abs() < 1e-15);
        // Three points: degree 2.
        assert_eq!(SplineMap::interpolated(&[0.0, 0.3, 1.0]).unwrap().degree, 2);
    }

    #[test]
    fn rejects_bad_definitions() {
        assert!(SplineMap::control(3, 3).is_err());
        assert!(SplineMap::control(3, 0).is_err());
        assert!(SplineMap::control(9, 6).is_err());
        assert!(SplineMap::interpolated(&[0.0]).is_err());
        assert!(SplineMap::interpolated(&[0.0, 0.5, 0.5, 1.0]).is_err());
        assert!(SplineMap::interpolated(&[0.1, 1.0]).is_err());
        assert!(SplineMap::interpolated(&[0.0, 0.9]).is_err());
        assert!(SplineMap::interpolated(&[0.0, f64::NAN, 1.0]).is_err());
    }

    /// The uniform periodic B-spline of degree `p` on the closed polygon
    /// `pts` at `t ∈ [0, 1]`, by de Boor's algorithm on the wrapped polygon
    /// and the integer knots (written independently of the map).
    fn periodic_point(pts: &[DVec2], p: usize, t: f64) -> DVec2 {
        let n = pts.len();
        let u = p as f64 + t * n as f64;
        let k = (u.floor() as usize).min(n + p - 1);
        let mut d: Vec<DVec2> = (0..=p).map(|j| pts[(j + k - p) % n]).collect();
        for r in 1..=p {
            for j in (r..=p).rev() {
                let i = (j + k - p) as f64;
                let alpha = (u - i) / (p + 1 - r) as f64;
                d[j] = d[j - 1] * (1.0 - alpha) + d[j] * alpha;
            }
        }
        d[p]
    }

    fn closed_polygon() -> Vec<DVec2> {
        vec![
            DVec2::new(0.0, 0.0),
            DVec2::new(3.0, -0.5),
            DVec2::new(4.0, 2.0),
            DVec2::new(1.5, 3.5),
            DVec2::new(-1.0, 2.0),
        ]
    }

    #[test]
    fn a_periodic_control_spline_is_the_wrapped_uniform_spline() {
        let poly = closed_polygon();
        let mut refs = poly.clone();
        refs.push(poly[0]);
        for p in 2..=4 {
            let m = SplineMap::control_periodic(refs.len(), p).unwrap();
            assert_eq!(m.count, poly.len() + p);
            assert!(m.kinks.is_empty());
            for i in 0..=200 {
                let t = i as f64 / 200.0;
                let (w, _) = m.weights(t);
                let c = combine(&w, &refs);
                let e = periodic_point(&poly, p, t);
                assert!((c - e).length() < 1e-13, "p={p} t={t}: {c} vs {e}");
            }
            // Closed exactly, and the parameter wraps.
            let ctrl = m.control_points(&refs);
            assert_eq!(ctrl[0], ctrl[ctrl.len() - 1]);
            let (a, _) = m.weights(0.25);
            let (b, _) = m.weights(1.25);
            let (c, _) = m.weights(-0.75);
            assert_eq!(a, b);
            assert_eq!(a, c);
        }
    }

    #[test]
    fn periodic_splines_are_smooth_through_the_seam() {
        let poly = closed_polygon();
        let mut refs = poly.clone();
        refs.push(poly[0]);
        let params = [0.0, 0.15, 0.4, 0.55, 0.8, 1.0];
        let maps = [
            SplineMap::control_periodic(refs.len(), 3).unwrap(),
            SplineMap::control_periodic(refs.len(), 2).unwrap(),
            SplineMap::interpolated_periodic(&params).unwrap(),
        ];
        for (k, m) in maps.iter().enumerate() {
            // The start, and the end of the period from its last piece (the
            // public evaluation wraps 1 to 0).
            let (w0, d0, dd0) = m.weights_ders(0.0, 2);
            let (w1, d1, dd1) = m.weights_in_span(m.count - 1, 1.0, 2);
            assert_eq!(combine(&w0, &refs), combine(&w1, &refs), "map {k}: closed");
            let (a, b) = (combine(&d0, &refs), combine(&d1, &refs));
            assert!(
                (a - b).length() < 1e-12 * a.length(),
                "map {k}: C′ {a} vs {b}"
            );
            // A quadratic is C1 at its knots, a cubic C2.
            if m.degree >= 3 {
                let (a2, b2) = (combine(&dd0, &refs), combine(&dd1, &refs));
                assert!(
                    (a2 - b2).length() < 1e-11 * a2.length(),
                    "map {k}: C″ {a2} vs {b2}"
                );
            }
            // A clamped spline on the same points has a corner there.
            let open = SplineMap::control(refs.len(), 3).unwrap();
            let (_, o0, _) = open.weights_ders(0.0, 1);
            let (_, o1, _) = open.weights_in_span(open.count - 1, 1.0, 1);
            let (oa, ob) = (combine(&o0, &refs), combine(&o1, &refs));
            assert!(oa.perp_dot(ob).abs() > 0.1 * oa.length() * ob.length());
        }
    }

    #[test]
    fn a_periodic_interpolated_spline_passes_through_its_points() {
        let poly = closed_polygon();
        let mut refs = poly.clone();
        refs.push(poly[0]);
        let params = [0.0, 0.15, 0.4, 0.55, 0.8, 1.0];
        let m = SplineMap::interpolated_periodic(&params).unwrap();
        assert_eq!((m.degree, m.count, m.refs), (3, 8, 6));
        for (t, q) in params.iter().zip(&refs) {
            let (w, _) = m.weights(*t);
            assert!((combine(&w, &refs) - *q).length() < 1e-13, "{t}");
        }
        let ctrl = m.control_points(&refs);
        assert_eq!(ctrl[0], poly[0]);
        assert_eq!(ctrl[ctrl.len() - 1], poly[0]);
        // The repeated point carries no weight of its own.
        for j in 0..m.count {
            assert_eq!(m.control_weights(j)[poly.len()], 0.0);
        }
    }

    #[test]
    fn periodic_definitions_are_checked() {
        assert!(SplineMap::control_periodic(4, 3).is_err()); // 3 points, degree 3
        assert!(SplineMap::control_periodic(4, 2).is_ok());
        assert!(SplineMap::control_periodic(5, 1).is_err());
        assert!(SplineMap::control_periodic(9, 6).is_err());
        assert!(SplineMap::interpolated_periodic(&[0.0, 0.5, 1.0]).is_err());
        assert!(SplineMap::interpolated_periodic(&[0.0, 0.3, 0.6, 1.0]).is_ok());
        assert!(SplineMap::interpolated_periodic(&[0.0, 0.6, 0.3, 1.0]).is_err());
    }

    #[test]
    fn knot_vectors() {
        assert_eq!(
            clamped_uniform_knots(5, 2),
            vec![0.0, 0.0, 0.0, 1.0 / 3.0, 2.0 / 3.0, 1.0, 1.0, 1.0]
        );
        let u = averaged_knots(&[0.0, 0.2, 0.45, 0.8, 1.0], 3);
        assert_eq!(u.len(), 9);
        assert!((u[4] - (0.2 + 0.45 + 0.8) / 3.0).abs() < 1e-15);
    }
}
