//! Ellipse geometry (design note 0008): `E(θ) = c + cos θ·a + sin θ·b`,
//! with `θ` the eccentric anomaly and `a`, `b` the semi-axis vectors (or,
//! while a solve is under way, any two conjugate semi-diameters).

use std::f64::consts::TAU;

use glam::DVec2;

/// `E(θ)`.
pub fn point(c: DVec2, a: DVec2, b: DVec2, theta: f64) -> DVec2 {
    let (s, co) = theta.sin_cos();
    c + a * co + b * s
}

/// `E′(θ) = −sin θ·a + cos θ·b`.
pub fn tangent(a: DVec2, b: DVec2, theta: f64) -> DVec2 {
    let (s, co) = theta.sin_cos();
    b * co - a * s
}

/// The eccentric anomaly of `p` seen from the centre: exact for a point
/// on the ellipse when `a ⟂ b`, and the parameter of the point where the
/// ray from the centre through `p`'s affine preimage meets it otherwise.
/// Zero when `p` is the centre or an axis has no length.
pub fn param_of(c: DVec2, a: DVec2, b: DVec2, p: DVec2) -> f64 {
    let (a2, b2) = (a.length_squared(), b.length_squared());
    if !(a2 > 0.0 && b2 > 0.0) {
        return 0.0;
    }
    let w = p - c;
    let (x, y) = (w.dot(a) / a2, w.dot(b) / b2);
    if x == 0.0 && y == 0.0 {
        0.0
    } else {
        y.atan2(x)
    }
}

/// Samples of the parameter used to find every local minimum of the
/// distance before refining (an ellipse has at most two).
const SAMPLES: usize = 64;

/// The parameter in `[0, 2π)` of the point of the ellipse closest to `p`.
///
/// The squared distance is sampled at [`SAMPLES`] parameters; each local
/// minimum of the samples is refined by golden-section search between its
/// neighbours (Kiefer 1953), then by safeguarded Newton on
/// `(E(θ) − p)·E′(θ) = 0` (steps must shrink that residual and stay in
/// the bracket); the nearest wins (the first on ties). Works
/// for any two non-parallel semi-diameters, so it can be used mid-solve.
pub fn closest_param(c: DVec2, a: DVec2, b: DVec2, p: DVec2) -> f64 {
    let d2 = |t: f64| (point(c, a, b, t) - p).length_squared();
    let step = TAU / SAMPLES as f64;
    let ts: Vec<f64> = (0..SAMPLES).map(|i| i as f64 * step).collect();
    let d: Vec<f64> = ts.iter().map(|&t| d2(t)).collect();
    let mut best = (f64::INFINITY, 0.0);
    for i in 0..SAMPLES {
        let left = d[(i + SAMPLES - 1) % SAMPLES];
        let right = d[(i + 1) % SAMPLES];
        if !(d[i] <= left && d[i] <= right) || (d[i] == left && d[i] == right) {
            continue;
        }
        let (t, dt) = refine(&d2, c, a, b, p, (ts[i], d[i]), ts[i] - step, ts[i] + step);
        if dt < best.0 {
            best = (dt, t);
        }
    }
    if best.0.is_infinite() {
        // Every sample equal (a circle about `p`): any parameter will do.
        return 0.0;
    }
    best.1.rem_euclid(TAU)
}

/// Refines a sampled local minimum `(θ, distance²)` within `[lo, hi]`.
#[allow(clippy::too_many_arguments)]
fn refine(
    d2: &dyn Fn(f64) -> f64,
    c: DVec2,
    a: DVec2,
    b: DVec2,
    p: DVec2,
    start: (f64, f64),
    lo: f64,
    hi: f64,
) -> (f64, f64) {
    let (mut best_t, mut best_d) = start;
    let (start_lo, start_hi) = (lo, hi);
    let golden = 0.5 * (5f64.sqrt() - 1.0);
    let (mut lo, mut hi) = (lo, hi);
    let (mut x1, mut x2) = (hi - golden * (hi - lo), lo + golden * (hi - lo));
    let (mut f1, mut f2) = (d2(x1), d2(x2));
    for _ in 0..40 {
        if f1 <= f2 {
            hi = x2;
            (x2, f2) = (x1, f1);
            x1 = hi - golden * (hi - lo);
            f1 = d2(x1);
        } else {
            lo = x1;
            (x1, f1) = (x2, f2);
            x2 = lo + golden * (hi - lo);
            f2 = d2(x2);
        }
    }
    for (t, d) in [(x1, f1), (x2, f2)] {
        if d < best_d {
            (best_t, best_d) = (t, d);
        }
    }
    // Newton polish on g(θ) = (E − p)·E′, g′ = ‖E′‖² + (E − p)·E″ with
    // E″ = −(E − c). The distance itself is flat to rounding near its
    // minimum (golden section stops at about √ε in θ), so steps are judged
    // by |g| instead, and kept within the sampled bracket.
    let g = |t: f64| {
        let e = point(c, a, b, t);
        let e1 = tangent(a, b, t);
        ((e - p).dot(e1), e1.length_squared() - (e - p).dot(e - c))
    };
    let (mut gr, mut gp) = g(best_t);
    for _ in 0..16 {
        if gp.is_nan() || gp <= 0.0 || gr == 0.0 {
            break;
        }
        let t = best_t - gr / gp;
        if !(t >= start_lo && t <= start_hi) || t == best_t {
            break;
        }
        let (gr2, gp2) = g(t);
        if gr2.is_nan() || gr2.abs() >= gr.abs() {
            break;
        }
        (best_t, gr, gp) = (t, gr2, gp2);
        best_d = d2(t);
    }
    (best_t, best_d)
}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;

    #[test]
    fn param_of_inverts_point() {
        let (c, a, b) = (
            DVec2::new(0.3, -0.2),
            DVec2::new(2.0, 1.0),
            DVec2::new(-0.5, 1.0),
        );
        for k in 0..12 {
            let t = -3.0 + 0.5 * k as f64;
            let back = param_of(c, a, b, point(c, a, b, t));
            let diff = (back - t).rem_euclid(TAU);
            assert!(diff.min(TAU - diff) < 1e-12, "{t} → {back}");
        }
    }

    proptest! {
        /// The closest point is a foot point (the offset is normal to the
        /// curve) and no farther than a dense sampling's nearest point.
        #[test]
        fn closest_param_is_the_nearest_foot_point(
            ax in 0.1f64..5.0, ratio in 0.05f64..1.0, rot in 0.0f64..TAU,
            px in -8.0f64..8.0, py in -8.0f64..8.0,
        ) {
            let c = DVec2::new(0.5, -1.0);
            let u = DVec2::from_angle(rot);
            let (a, b) = (u * ax, u.perp() * ax * ratio);
            let p = c + DVec2::new(px, py);
            let t = closest_param(c, a, b, p);
            let e = point(c, a, b, t);
            let d = (e - p).length();
            let dense = (0..20_000)
                .map(|i| (point(c, a, b, TAU * i as f64 / 20_000.0) - p).length())
                .fold(f64::INFINITY, f64::min);
            prop_assert!(d <= dense + 1e-12 * (1.0 + dense), "{d} vs {dense}");
            let tan = tangent(a, b, t);
            prop_assert!((e - p).dot(tan).abs() <= 1e-9 * tan.length() * (1.0 + d));
        }
    }
}
