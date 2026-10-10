//! Householder reflectors `H = I − τ v vᵀ` with `v[0] = 1` (Golub & Van Loan,
//! *Matrix Computations*, 4th ed., §5.1.2; the conventions of LAPACK
//! `dlarfg`/`dlarf`).

use crate::vector::{axpy, dot, norm2};

/// Turns `x` into the data of a reflector `H` with `H x = (β, 0, …, 0)`:
/// on return `x[0] = β` and `x[1..]` holds `v[1..]`. Returns `τ`; `τ = 0`
/// means `H = I` (then `x` is unchanged). `|β| = ‖x‖₂`.
pub(crate) fn make_reflector(x: &mut [f64]) -> f64 {
    if x.len() <= 1 {
        return 0.0;
    }
    let alpha = x[0];
    let xnorm = norm2(&x[1..]);
    if xnorm == 0.0 {
        return 0.0;
    }
    let h = norm2(&[alpha, xnorm]);
    let beta = if alpha >= 0.0 { -h } else { h };
    let tau = (beta - alpha) / beta;
    // |alpha − beta| = |alpha| + h ≥ every |x[i]|, so the quotients are ≤ 1.
    let denom = alpha - beta;
    for v in &mut x[1..] {
        *v /= denom;
    }
    x[0] = beta;
    tau
}

/// Applies `H = I − τ v vᵀ` (with `v = (1, tail)`) to `y` in place.
/// `y.len()` must be `tail.len() + 1`.
pub(crate) fn apply_reflector(tau: f64, tail: &[f64], y: &mut [f64]) {
    if tau == 0.0 {
        return;
    }
    debug_assert_eq!(tail.len() + 1, y.len());
    let (y0, rest) = y.split_first_mut().expect("non-empty");
    let w = *y0 + dot(tail, rest);
    let s = tau * w;
    *y0 -= s;
    axpy(-s, tail, rest);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reflector_annihilates_tail() {
        let orig = [3.0, 1.0, -2.0, 0.5];
        let mut x = orig;
        let tau = make_reflector(&mut x);
        let beta = x[0];
        assert!((beta.abs() - norm2(&orig)).abs() < 1e-14);
        let mut y = orig;
        apply_reflector(tau, &x[1..], &mut y);
        assert!((y[0] - beta).abs() < 1e-14);
        for v in &y[1..] {
            assert!(v.abs() < 1e-14);
        }
        // H is an involution.
        apply_reflector(tau, &x[1..], &mut y);
        for (a, b) in y.iter().zip(&orig) {
            assert!((a - b).abs() < 1e-14);
        }
    }

    #[test]
    fn zero_tail_is_identity() {
        let mut x = [-2.0, 0.0, 0.0];
        assert_eq!(make_reflector(&mut x), 0.0);
        assert_eq!(x, [-2.0, 0.0, 0.0]);
        let mut one = [5.0];
        assert_eq!(make_reflector(&mut one), 0.0);
    }
}
