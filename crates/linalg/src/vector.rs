//! Vector kernels with a fixed operation order (ADR 0010).
//!
//! Reductions use four interleaved partial sums combined in a fixed pattern.
//! The order is part of the code, not chosen at run time, so the result is
//! the same on every run and every machine with IEEE-754 doubles. Rust never
//! contracts `a * b + c` into a fused multiply-add unless asked, and nothing
//! here asks.

/// Dot product `Σ a[i]·b[i]` over the common length, with the fixed order
/// `((s0 + s1) + (s2 + s3)) + tail`, where `sk` sums the terms with index
/// `≡ k (mod 4)` in increasing index order and `tail` sums the last
/// `len % 4` terms.
pub fn dot(a: &[f64], b: &[f64]) -> f64 {
    let n = a.len().min(b.len());
    let (a, b) = (&a[..n], &b[..n]);
    let mut s0 = 0.0;
    let mut s1 = 0.0;
    let mut s2 = 0.0;
    let mut s3 = 0.0;
    let mut ca = a.chunks_exact(4);
    let mut cb = b.chunks_exact(4);
    for (x, y) in (&mut ca).zip(&mut cb) {
        s0 += x[0] * y[0];
        s1 += x[1] * y[1];
        s2 += x[2] * y[2];
        s3 += x[3] * y[3];
    }
    let mut tail = 0.0;
    for (x, y) in ca.remainder().iter().zip(cb.remainder()) {
        tail += x * y;
    }
    ((s0 + s1) + (s2 + s3)) + tail
}

/// Euclidean norm `‖x‖₂`, safe from overflow and underflow: values are
/// scaled by the largest magnitude when their squares would leave the normal
/// range.
pub fn norm2(x: &[f64]) -> f64 {
    let mut big = 0.0f64;
    for v in x {
        let a = v.abs();
        if a > big {
            big = a;
        }
    }
    if big == 0.0 || !big.is_finite() {
        return big;
    }
    // Squares of values in [1e-150, 1e150] stay normal and finite.
    if (1e-150..=1e150).contains(&big) {
        return dot(x, x).sqrt();
    }
    // Scale by an exact power of two so that the largest magnitude lands in
    // the normal range; the scaling is exact and is undone at the end.
    let (scale, unscale) = if big < 1e-150 {
        (2f64.powi(600), 2f64.powi(-600))
    } else {
        (2f64.powi(-600), 2f64.powi(600))
    };
    let mut s0 = 0.0;
    let mut s1 = 0.0;
    let mut s2 = 0.0;
    let mut s3 = 0.0;
    let mut c = x.chunks_exact(4);
    for q in &mut c {
        let (a, b, d, e) = (q[0] * scale, q[1] * scale, q[2] * scale, q[3] * scale);
        s0 += a * a;
        s1 += b * b;
        s2 += d * d;
        s3 += e * e;
    }
    let mut tail = 0.0;
    for v in c.remainder() {
        let a = v * scale;
        tail += a * a;
    }
    (((s0 + s1) + (s2 + s3)) + tail).sqrt() * unscale
}

/// `y[i] += alpha · x[i]` over the common length.
pub fn axpy(alpha: f64, x: &[f64], y: &mut [f64]) {
    for (yi, xi) in y.iter_mut().zip(x) {
        *yi += alpha * xi;
    }
}

/// `‖x‖∞`, the largest magnitude (0 for an empty slice).
pub fn norm_inf(x: &[f64]) -> f64 {
    x.iter().fold(0.0f64, |m, v| m.max(v.abs()))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn dot_matches_naive_on_exact_values() {
        let a: Vec<f64> = (0..11).map(|i| i as f64).collect();
        let b: Vec<f64> = (0..11).map(|i| (2 * i + 1) as f64).collect();
        let naive: f64 = a.iter().zip(&b).map(|(x, y)| x * y).sum();
        assert_eq!(dot(&a, &b), naive);
    }

    #[test]
    fn norm2_handles_extreme_magnitudes() {
        assert!((norm2(&[3e200, 4e200]) / 5e200 - 1.0).abs() < 1e-15);
        assert!((norm2(&[3e-200, 4e-200]) / 5e-200 - 1.0).abs() < 1e-15);
        assert_eq!(norm2(&[]), 0.0);
        assert_eq!(norm2(&[0.0, -0.0]), 0.0);
        let tiny = f64::MIN_POSITIVE / 8.0; // subnormal
        assert!((norm2(&[tiny, 0.0]) / tiny - 1.0).abs() < 1e-15);
        assert!((norm2(&[1e300, 1e300]) / (1e300 * 2f64.sqrt()) - 1.0).abs() < 1e-15);
        assert!((norm2(&[1.0, 2.0, 2.0]) - 3.0).abs() < 1e-15);
    }

    #[test]
    fn axpy_and_inf_norm() {
        let mut y = vec![1.0, 2.0];
        axpy(2.0, &[3.0, -4.0], &mut y);
        assert_eq!(y, vec![7.0, -6.0]);
        assert_eq!(norm_inf(&y), 7.0);
    }
}
