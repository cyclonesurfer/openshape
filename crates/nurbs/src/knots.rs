//! Knot vectors and B-spline basis functions (The NURBS Book, chapter 2).

use serde::{Deserialize, Serialize};

use crate::NurbsError;

/// A non-decreasing sequence of finite knot values.
///
/// The vector itself does not know the degree of the entity it belongs to;
/// methods that depend on it take the degree as an argument. Indices follow
/// The NURBS Book: for `n + 1` control points of degree `p` there are
/// `m + 1 = n + p + 2` knots `u_0 ..= u_m`, and the parametric domain is
/// `[u_p, u_{n+1}]`.
///
/// Deserialization validates the values exactly like [`KnotVector::new`].
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(try_from = "Vec<f64>", into = "Vec<f64>")]
pub struct KnotVector {
    knots: Vec<f64>,
}

impl KnotVector {
    /// Builds a knot vector after checking that it has at least two values,
    /// that every value is finite, that the sequence never decreases and that
    /// it is not constant.
    pub fn new(knots: Vec<f64>) -> Result<Self, NurbsError> {
        if knots.len() < 2 {
            return Err(NurbsError::TooFewKnots { found: knots.len() });
        }
        if let Some(index) = knots.iter().position(|k| !k.is_finite()) {
            return Err(NurbsError::NonFiniteKnot { index });
        }
        if let Some(index) = (1..knots.len()).find(|&i| knots[i] < knots[i - 1]) {
            return Err(NurbsError::DecreasingKnots { index });
        }
        if knots[0] == knots[knots.len() - 1] {
            return Err(NurbsError::EmptyDomain);
        }
        Ok(Self { knots })
    }

    /// Clamped (open) uniform knot vector on `[0, 1]` for `n_ctrl` control
    /// points of the given degree: `degree + 1` zeros, evenly spaced simple
    /// interior knots and `degree + 1` ones.
    pub fn clamped_uniform(degree: usize, n_ctrl: usize) -> Result<Self, NurbsError> {
        if degree == 0 {
            return Err(NurbsError::ZeroDegree);
        }
        if n_ctrl < degree + 1 {
            return Err(NurbsError::TooFewControlPoints {
                degree,
                required: degree + 1,
                found: n_ctrl,
            });
        }
        let interior = n_ctrl - degree - 1;
        let mut knots = Vec::with_capacity(n_ctrl + degree + 1);
        knots.extend(std::iter::repeat_n(0.0, degree + 1));
        knots.extend((1..=interior).map(|i| i as f64 / (interior + 1) as f64));
        knots.extend(std::iter::repeat_n(1.0, degree + 1));
        Ok(Self { knots })
    }

    /// The knot values.
    pub fn as_slice(&self) -> &[f64] {
        &self.knots
    }

    /// Number of knots (`m + 1`).
    pub fn len(&self) -> usize {
        self.knots.len()
    }

    /// Always `false`: a valid knot vector has at least two values.
    pub fn is_empty(&self) -> bool {
        self.knots.is_empty()
    }

    /// Consumes the knot vector and returns its values.
    pub fn into_vec(self) -> Vec<f64> {
        self.knots
    }

    /// Parametric domain `[u_p, u_{m-p}]` for the given degree.
    ///
    /// # Panics
    ///
    /// Panics if the vector has fewer than `2 * degree + 2` knots.
    pub fn domain(&self, degree: usize) -> (f64, f64) {
        let m = self.knots.len() - 1;
        (self.knots[degree], self.knots[m - degree])
    }

    /// Number of knots exactly equal to `u`.
    pub fn multiplicity(&self, u: f64) -> usize {
        // The sequence is sorted, so the equal values form one contiguous run.
        let start = self.knots.partition_point(|&k| k < u);
        let end = self.knots.partition_point(|&k| k <= u);
        end - start
    }

    /// Knot span index `i` with `u_i <= u < u_{i+1}` (Algorithm A2.1, FindSpan).
    ///
    /// `n` is the index of the last control point (the number of control
    /// points minus one), as in the book. Parameters outside the domain are clamped to it,
    /// `u = u_{n+1}` maps to the last non-empty span, and the returned span is
    /// always non-empty, so the basis function routines never divide by zero.
    ///
    /// # Panics
    ///
    /// Panics if `n + degree + 1` is not a valid knot index.
    pub fn find_span(&self, degree: usize, n: usize, u: f64) -> usize {
        let knots = &self.knots;
        let p = degree;
        if u >= knots[n + 1] {
            // Special case of A2.1, extended to skip empty trailing spans.
            let mut span = n;
            while span > p && knots[span] >= knots[n + 1] {
                span -= 1;
            }
            return span;
        }
        let u = if u < knots[p] { knots[p] } else { u };
        // Binary search for the span; also lands on the last of several equal
        // knots so that the span is non-empty.
        let mut low = p;
        let mut high = n + 1;
        let mut mid = (low + high) / 2;
        while u < knots[mid] || u >= knots[mid + 1] {
            if u < knots[mid] {
                high = mid;
            } else {
                low = mid;
            }
            let next = (low + high) / 2;
            if next == mid {
                // Only reachable for NaN input; avoids looping forever.
                break;
            }
            mid = next;
        }
        mid
    }

    /// The `degree + 1` non-vanishing basis functions `N_{span-p,p}(u) ..=
    /// N_{span,p}(u)` (Algorithm A2.2, BasisFuns).
    ///
    /// # Panics
    ///
    /// Panics if `span + degree` is not a valid knot index or `span < degree`.
    pub fn basis_funs(&self, span: usize, u: f64, degree: usize) -> Vec<f64> {
        let knots = &self.knots;
        let p = degree;
        let mut basis = vec![0.0; p + 1];
        let mut left = vec![0.0; p + 1];
        let mut right = vec![0.0; p + 1];
        basis[0] = 1.0;
        for j in 1..=p {
            left[j] = u - knots[span + 1 - j];
            right[j] = knots[span + j] - u;
            let mut saved = 0.0;
            for r in 0..j {
                let temp = basis[r] / (right[r + 1] + left[j - r]);
                basis[r] = saved + right[r + 1] * temp;
                saved = left[j - r] * temp;
            }
            basis[j] = saved;
        }
        basis
    }

    /// Non-vanishing basis functions and their derivatives up to order `n`
    /// (Algorithm A2.3, DersBasisFuns).
    ///
    /// Returns `ders` with `ders[k][j]` the `k`-th derivative of
    /// `N_{span-p+j,p}` at `u`, for `k in 0..=n` and `j in 0..=degree`.
    /// Derivatives of order greater than `degree` are zero.
    ///
    /// # Panics
    ///
    /// Panics if `span + degree` is not a valid knot index or `span < degree`.
    pub fn ders_basis_funs(&self, span: usize, u: f64, degree: usize, n: usize) -> Vec<Vec<f64>> {
        let knots = &self.knots;
        let p = degree;
        let mut ders = vec![vec![0.0; p + 1]; n + 1];
        // ndu stores the basis functions (upper triangle) and knot differences
        // (lower triangle).
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
        for j in 0..=p {
            ders[0][j] = ndu[j][p];
        }
        let nd = n.min(p);
        if nd == 0 {
            return ders;
        }
        // Two alternating rows of coefficients a_{k,j}.
        let mut a = [vec![0.0; p + 1], vec![0.0; p + 1]];
        let pi = p as isize;
        for r in 0..=p {
            let ri = r as isize;
            let (mut s1, mut s2) = (0usize, 1usize);
            a[0][0] = 1.0;
            for k in 1..=nd {
                let ki = k as isize;
                let mut d = 0.0;
                let rk = ri - ki;
                let pk = pi - ki;
                if ri >= ki {
                    let v = a[s1][0] / ndu[(pk + 1) as usize][rk as usize];
                    a[s2][0] = v;
                    d = v * ndu[rk as usize][pk as usize];
                }
                let j1: isize = if rk >= -1 { 1 } else { -rk };
                let j2: isize = if ri - 1 <= pk { ki - 1 } else { pi - ri };
                let mut j = j1;
                while j <= j2 {
                    let ju = j as usize;
                    let v = (a[s1][ju] - a[s1][ju - 1]) / ndu[(pk + 1) as usize][(rk + j) as usize];
                    a[s2][ju] = v;
                    d += v * ndu[(rk + j) as usize][pk as usize];
                    j += 1;
                }
                if ri <= pk {
                    let v = -a[s1][k - 1] / ndu[(pk + 1) as usize][r];
                    a[s2][k] = v;
                    d += v * ndu[r][pk as usize];
                }
                ders[k][r] = d;
                std::mem::swap(&mut s1, &mut s2);
            }
        }
        // Multiply through by the factors p! / (p - k)!.
        let mut factor = p as f64;
        for (k, row) in ders.iter_mut().enumerate().take(nd + 1).skip(1) {
            for value in row.iter_mut() {
                *value *= factor;
            }
            factor *= (p - k) as f64;
        }
        ders
    }

    /// Wraps values that are already known to be valid (used for the results
    /// of algorithms such as knot insertion).
    pub(crate) fn from_vec_unchecked(knots: Vec<f64>) -> Self {
        debug_assert!(knots.windows(2).all(|w| w[0] <= w[1]));
        Self { knots }
    }

    /// Checks that this knot vector fits `n_ctrl` control points of the given
    /// degree: correct length, non-empty domain and admissible multiplicities.
    pub(crate) fn validate_for(&self, degree: usize, n_ctrl: usize) -> Result<(), NurbsError> {
        if degree == 0 {
            return Err(NurbsError::ZeroDegree);
        }
        if n_ctrl < degree + 1 {
            return Err(NurbsError::TooFewControlPoints {
                degree,
                required: degree + 1,
                found: n_ctrl,
            });
        }
        let expected = n_ctrl + degree + 1;
        if self.knots.len() != expected {
            return Err(NurbsError::KnotCountMismatch {
                degree,
                control_points: n_ctrl,
                expected,
                found: self.knots.len(),
            });
        }
        let (a, b) = self.domain(degree);
        if a >= b {
            return Err(NurbsError::EmptyDomain);
        }
        let mut i = 0;
        while i < self.knots.len() {
            let value = self.knots[i];
            let mut j = i;
            while j < self.knots.len() && self.knots[j] == value {
                j += 1;
            }
            let multiplicity = j - i;
            let max = if value > a && value < b {
                degree
            } else {
                degree + 1
            };
            if multiplicity > max {
                return Err(NurbsError::KnotMultiplicityTooHigh {
                    value,
                    multiplicity,
                    max,
                });
            }
            i = j;
        }
        Ok(())
    }

    /// Distinct knot values strictly inside `(a, b)` together with the domain
    /// ends, i.e. the breakpoints of the non-empty spans of the domain.
    pub(crate) fn breakpoints(&self, degree: usize) -> Vec<f64> {
        let (a, b) = self.domain(degree);
        let mut out = vec![a];
        for &k in &self.knots {
            if k > a && k < b && k > *out.last().unwrap_or(&a) {
                out.push(k);
            }
        }
        out.push(b);
        out
    }
}

impl TryFrom<Vec<f64>> for KnotVector {
    type Error = NurbsError;

    fn try_from(knots: Vec<f64>) -> Result<Self, Self::Error> {
        Self::new(knots)
    }
}

impl From<KnotVector> for Vec<f64> {
    fn from(knots: KnotVector) -> Self {
        knots.knots
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn book_example_2_3() {
        // Example 2.3 of The NURBS Book: U = {0,0,0,1,2,3,4,4,5,5,5}, p = 2, u = 5/2.
        let kv = KnotVector::new(vec![0., 0., 0., 1., 2., 3., 4., 4., 5., 5., 5.]).unwrap();
        let span = kv.find_span(2, 7, 2.5);
        assert_eq!(span, 4);
        let n = kv.basis_funs(span, 2.5, 2);
        assert!((n[0] - 1.0 / 8.0).abs() < 1e-15);
        assert!((n[1] - 6.0 / 8.0).abs() < 1e-15);
        assert!((n[2] - 1.0 / 8.0).abs() < 1e-15);
        let d = kv.ders_basis_funs(span, 2.5, 2, 2);
        // N'_{2,2}(5/2) = -1/2, N'_{3,2} = 0, N'_{4,2} = 1/2.
        assert!((d[1][0] + 0.5).abs() < 1e-15);
        assert!(d[1][1].abs() < 1e-15);
        assert!((d[1][2] - 0.5).abs() < 1e-15);
        // Second derivatives 1, -2, 1.
        assert!((d[2][0] - 1.0).abs() < 1e-14);
        assert!((d[2][1] + 2.0).abs() < 1e-14);
        assert!((d[2][2] - 1.0).abs() < 1e-14);
    }

    #[test]
    fn find_span_ends_and_repeated_knots() {
        let kv = KnotVector::new(vec![0., 0., 0., 1., 2., 3., 4., 4., 5., 5., 5.]).unwrap();
        assert_eq!(kv.find_span(2, 7, 0.0), 2);
        assert_eq!(kv.find_span(2, 7, -1.0), 2);
        assert_eq!(kv.find_span(2, 7, 5.0), 7);
        assert_eq!(kv.find_span(2, 7, 7.0), 7);
        assert_eq!(kv.find_span(2, 7, 4.0), 7);
        assert_eq!(kv.find_span(2, 7, 3.999), 5);
        assert_eq!(kv.find_span(2, 7, 1.0), 3);
        // NaN must not hang or panic.
        let _ = kv.find_span(2, 7, f64::NAN);
    }

    #[test]
    fn multiplicity_and_domain() {
        let kv = KnotVector::new(vec![0., 0., 0., 1., 2., 3., 4., 4., 5., 5., 5.]).unwrap();
        assert_eq!(kv.multiplicity(0.0), 3);
        assert_eq!(kv.multiplicity(4.0), 2);
        assert_eq!(kv.multiplicity(2.5), 0);
        assert_eq!(kv.domain(2), (0.0, 5.0));
        assert_eq!(kv.breakpoints(2), vec![0., 1., 2., 3., 4., 5.]);
    }

    #[test]
    fn clamped_uniform_shape() {
        let kv = KnotVector::clamped_uniform(3, 6).unwrap();
        assert_eq!(
            kv.as_slice(),
            &[0., 0., 0., 0., 1. / 3., 2. / 3., 1., 1., 1., 1.]
        );
        assert!(KnotVector::clamped_uniform(0, 3).is_err());
        assert!(KnotVector::clamped_uniform(3, 3).is_err());
    }

    #[test]
    fn validation() {
        assert!(KnotVector::new(vec![]).is_err());
        assert!(KnotVector::new(vec![1.0]).is_err());
        assert!(KnotVector::new(vec![0.0, 0.0]).is_err());
        assert!(KnotVector::new(vec![0.0, f64::NAN, 1.0]).is_err());
        assert_eq!(
            KnotVector::new(vec![0.0, 1.0, 0.5]),
            Err(NurbsError::DecreasingKnots { index: 2 })
        );
        let kv = KnotVector::new(vec![0., 0., 0., 0.5, 0.5, 0.5, 1., 1., 1.]).unwrap();
        assert!(matches!(
            kv.validate_for(2, 6),
            Err(NurbsError::KnotMultiplicityTooHigh { .. })
        ));
        let kv = KnotVector::new(vec![0., 0., 0., 0.5, 0.5, 1., 1., 1.]).unwrap();
        assert!(kv.validate_for(2, 5).is_ok());
        assert!(kv.validate_for(2, 4).is_err());
    }

    #[test]
    fn serde_validates() {
        let kv = KnotVector::new(vec![0., 0., 1., 1.]).unwrap();
        let json = serde_json::to_string(&kv).unwrap();
        assert_eq!(json, "[0.0,0.0,1.0,1.0]");
        let back: KnotVector = serde_json::from_str(&json).unwrap();
        assert_eq!(back, kv);
        assert!(serde_json::from_str::<KnotVector>("[1.0,0.0]").is_err());
    }
}
