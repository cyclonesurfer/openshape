//! Small numerical helpers shared by curves and surfaces.

use std::sync::OnceLock;

use glam::{DMat4, DVec3, DVec4};

/// Binomial coefficient `n choose k` as a float.
pub(crate) fn binomial(n: usize, k: usize) -> f64 {
    if k > n {
        return 0.0;
    }
    let k = k.min(n - k);
    let mut result = 1.0;
    for i in 0..k {
        result = result * (n - i) as f64 / (i + 1) as f64;
    }
    result.round()
}

/// Projects a homogeneous point `(w x, w y, w z, w)` to Cartesian space.
pub(crate) fn dehomogenize(pw: DVec4) -> DVec3 {
    pw.truncate() / pw.w
}

/// Lifts a Cartesian point with weight `w` to homogeneous space.
pub(crate) fn homogenize(p: DVec3, w: f64) -> DVec4 {
    (p * w).extend(w)
}

/// Applies the affine part of `m` to a homogeneous control point, keeping its
/// weight: `(A (w P) + w t, w)`. Exact for affine maps; the projective row of
/// `m` is ignored.
pub(crate) fn transform_homogeneous(m: &DMat4, pw: DVec4) -> DVec4 {
    let a = m.x_axis.truncate() * pw.x + m.y_axis.truncate() * pw.y + m.z_axis.truncate() * pw.z;
    (a + m.w_axis.truncate() * pw.w).extend(pw.w)
}

/// Axis-aligned bounding box of a set of homogeneous points (after projection).
pub(crate) fn bounding_box(points: &[DVec4]) -> (DVec3, DVec3) {
    let mut lo = DVec3::splat(f64::INFINITY);
    let mut hi = DVec3::splat(f64::NEG_INFINITY);
    for &pw in points {
        let p = dehomogenize(pw);
        lo = lo.min(p);
        hi = hi.max(p);
    }
    (lo, hi)
}

/// Squared distance from `p` to the axis-aligned box `[lo, hi]` (zero inside).
pub(crate) fn aabb_distance_squared(lo: DVec3, hi: DVec3, p: DVec3) -> f64 {
    (lo - p).max(p - hi).max(DVec3::ZERO).length_squared()
}

/// Number of Gauss-Legendre nodes used per quadrature panel.
const GAUSS_ORDER: usize = 16;

/// Gauss-Legendre nodes and weights on `[-1, 1]`, computed once by Newton
/// iteration on the Legendre polynomial `P_16`.
fn gauss_legendre() -> &'static [(f64, f64)] {
    static TABLE: OnceLock<Vec<(f64, f64)>> = OnceLock::new();
    TABLE.get_or_init(|| {
        let n = GAUSS_ORDER;
        let mut table = Vec::with_capacity(n);
        for i in 0..n.div_ceil(2) {
            let mut x = (std::f64::consts::PI * (i as f64 + 0.75) / (n as f64 + 0.5)).cos();
            let mut dp = 1.0;
            for _ in 0..100 {
                let (mut p1, mut p2) = (1.0, 0.0);
                for j in 1..=n {
                    let p3 = p2;
                    p2 = p1;
                    p1 = ((2 * j - 1) as f64 * x * p2 - (j - 1) as f64 * p3) / j as f64;
                }
                dp = n as f64 * (x * p1 - p2) / (x * x - 1.0);
                let dx = p1 / dp;
                x -= dx;
                if dx.abs() < 1e-16 {
                    break;
                }
            }
            let w = 2.0 / ((1.0 - x * x) * dp * dp);
            table.push((x, w));
            table.push((-x, w));
        }
        table
    })
}

/// Gauss-Legendre approximation of the integral of `f` over `[a, b]`.
fn gauss_panel(f: &mut impl FnMut(f64) -> f64, a: f64, b: f64) -> f64 {
    let half = 0.5 * (b - a);
    let mid = 0.5 * (a + b);
    gauss_legendre()
        .iter()
        .map(|&(x, w)| w * f(mid + half * x))
        .sum::<f64>()
        * half
}

/// Maximum number of panel bisections in [`integrate`].
const MAX_BISECTIONS: usize = 400;

/// One quadrature panel: its interval, its integral estimate (sum of the two
/// half-panel rules), the half values, and the error estimate.
struct Panel {
    a: f64,
    b: f64,
    left: f64,
    right: f64,
    err: f64,
}

impl Panel {
    fn new(f: &mut impl FnMut(f64) -> f64, a: f64, b: f64, whole: f64) -> Self {
        let m = 0.5 * (a + b);
        let left = gauss_panel(f, a, m);
        let right = gauss_panel(f, m, b);
        let err = (left + right - whole).abs();
        // NaN errors are treated as converged so they cannot attract work.
        let err = if err.is_nan() { 0.0 } else { err };
        Self {
            a,
            b,
            left,
            right,
            err,
        }
    }
}

/// Globally adaptive Gauss-Legendre quadrature of `f` over the union of the
/// given intervals (each of which should contain no kink of `f`).
///
/// Every interval starts as one panel whose error is estimated by comparing
/// the 16-point rule on the whole panel with the sum over its two halves. The
/// panel with the largest error is bisected until the total error estimate
/// drops below `tol` or a fixed budget of bisections is spent, so the cost is
/// bounded even when `tol` is below the attainable rounding noise.
pub(crate) fn integrate(f: &mut impl FnMut(f64) -> f64, intervals: &[(f64, f64)], tol: f64) -> f64 {
    let mut panels: Vec<Panel> = intervals
        .iter()
        .filter(|(a, b)| b > a)
        .map(|&(a, b)| {
            let whole = gauss_panel(f, a, b);
            Panel::new(f, a, b, whole)
        })
        .collect();
    for _ in 0..MAX_BISECTIONS {
        let total: f64 = panels.iter().map(|p| p.err).sum();
        if total <= tol {
            break;
        }
        let Some(worst) = panels
            .iter()
            .enumerate()
            .max_by(|x, y| x.1.err.total_cmp(&y.1.err))
            .map(|(i, _)| i)
        else {
            break;
        };
        let panel = panels.swap_remove(worst);
        let m = 0.5 * (panel.a + panel.b);
        if m <= panel.a || m >= panel.b {
            // Cannot bisect further; keep it as converged.
            panels.push(Panel { err: 0.0, ..panel });
            continue;
        }
        panels.push(Panel::new(f, panel.a, m, panel.left));
        panels.push(Panel::new(f, m, panel.b, panel.right));
    }
    panels.iter().map(|p| p.left + p.right).sum()
}

/// Solves `A x = b` for a banded matrix whose row `i` has its non-zeros in
/// columns `lo[i] ..= lo[i] + width - 1`, with `lo` non-decreasing and
/// `lo[i] <= i <= lo[i] + width - 1`.
///
/// Gaussian elimination without pivoting is stable for the totally positive
/// B-spline collocation matrices this is used for (de Boor). `rows[i][k]` is
/// the entry in column `lo[i] + k`. Each right-hand side is a [`DVec3`].
/// Returns `None` when a pivot vanishes.
pub(crate) fn solve_banded(
    mut rows: Vec<Vec<f64>>,
    lo: &[usize],
    mut rhs: Vec<DVec3>,
) -> Option<Vec<DVec3>> {
    let n = rows.len();
    for c in 0..n {
        let pivot = rows[c][c - lo[c]];
        if !pivot.is_finite() || pivot.abs() < 1e-14 {
            return None;
        }
        let pivot_end = lo[c] + rows[c].len();
        let mut r = c + 1;
        while r < n && lo[r] <= c {
            let factor = rows[r][c - lo[r]] / pivot;
            if factor != 0.0 {
                for col in c..pivot_end {
                    let v = rows[c][col - lo[c]];
                    rows[r][col - lo[r]] -= factor * v;
                }
                let b = rhs[c];
                rhs[r] -= b * factor;
            }
            r += 1;
        }
    }
    let mut x = vec![DVec3::ZERO; n];
    for c in (0..n).rev() {
        let mut acc = rhs[c];
        let end = (lo[c] + rows[c].len()).min(n);
        for (col, xc) in x.iter().enumerate().take(end).skip(c + 1) {
            acc -= *xc * rows[c][col - lo[c]];
        }
        x[c] = acc / rows[c][c - lo[c]];
    }
    Some(x)
}

/// Indices of the local minima of a sampled sequence (`d[i]` no larger than
/// its neighbours). A run of equal minima (a plateau) is reduced to its two
/// end indices so that plateaus cannot crowd out other candidates.
pub(crate) fn local_minima(d: &[f64]) -> Vec<usize> {
    let n = d.len();
    let is_min = |i: usize| (i == 0 || d[i] <= d[i - 1]) && (i + 1 == n || d[i] <= d[i + 1]);
    let mut out = Vec::new();
    let mut i = 0;
    while i < n {
        if is_min(i) {
            let start = i;
            while i + 1 < n && is_min(i + 1) && d[i + 1] == d[start] {
                i += 1;
            }
            out.push(start);
            if i != start {
                out.push(i);
            }
        }
        i += 1;
    }
    out
}

/// Unit vector along `v`, or `None` if `v` is zero or not finite.
pub(crate) fn unit(v: DVec3) -> Option<DVec3> {
    let len = v.length();
    (len.is_finite() && len > 0.0).then(|| v / len)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn local_minima_collapse_plateaus() {
        assert_eq!(
            local_minima(&[3.0, 1.0, 2.0, 0.5, 0.5, 0.5, 4.0]),
            vec![1, 3, 5]
        );
        assert_eq!(local_minima(&[1.0, 2.0, 3.0]), vec![0]);
        assert_eq!(local_minima(&[2.0, 2.0]), vec![0, 1]);
        let _ = local_minima(&[f64::NAN, 1.0]);
        assert!(local_minima(&[]).is_empty());
    }

    #[test]
    fn binomials() {
        assert_eq!(binomial(0, 0), 1.0);
        assert_eq!(binomial(5, 2), 10.0);
        assert_eq!(binomial(10, 5), 252.0);
        assert_eq!(binomial(3, 4), 0.0);
    }

    #[test]
    fn gauss_legendre_is_exact_for_polynomials() {
        let table = gauss_legendre();
        assert_eq!(table.len(), GAUSS_ORDER);
        let weight_sum: f64 = table.iter().map(|&(_, w)| w).sum();
        assert!((weight_sum - 2.0).abs() < 1e-14);
        // x^30 is integrated exactly by a 16-point rule.
        let value = gauss_panel(&mut |x: f64| x.powi(30), -1.0, 1.0);
        assert!((value - 2.0 / 31.0).abs() < 1e-14);
    }

    #[test]
    fn adaptive_integration() {
        let value = integrate(&mut |x: f64| x.sqrt(), &[(0.0, 1.0)], 1e-13);
        assert!((value - 2.0 / 3.0).abs() < 1e-12);
        let value = integrate(&mut |x: f64| x.sin(), &[(0.0, std::f64::consts::PI)], 1e-14);
        assert!((value - 2.0).abs() < 1e-14);
        // Kink between intervals, plus an unreachable tolerance: still finite work.
        let value = integrate(&mut |x: f64| x.abs(), &[(-1.0, 0.0), (0.0, 2.0)], 0.0);
        assert!((value - 2.5).abs() < 1e-14);
        let mut calls = 0;
        let value = integrate(
            &mut |x: f64| {
                calls += 1;
                (x * 1e6).sin().abs()
            },
            &[(0.0, 1.0)],
            0.0,
        );
        assert!(value.is_finite());
        assert!(calls <= 16 * (3 + 4 * MAX_BISECTIONS));
    }

    #[test]
    fn banded_solver() {
        // Tridiagonal system with known solution.
        let rows = vec![
            vec![2.0, 1.0],
            vec![1.0, 2.0, 1.0],
            vec![1.0, 2.0, 1.0],
            vec![1.0, 2.0],
        ];
        let lo = [0, 0, 1, 2];
        let x_true = [1.0, -2.0, 3.0, 0.5];
        let b: Vec<DVec3> = (0..4)
            .map(|i| {
                let s: f64 = rows[i]
                    .iter()
                    .enumerate()
                    .map(|(k, a)| a * x_true[lo[i] + k])
                    .sum();
                DVec3::splat(s)
            })
            .collect();
        let x = solve_banded(rows, &lo, b).unwrap();
        for i in 0..4 {
            assert!((x[i].x - x_true[i]).abs() < 1e-14);
        }
    }
}
