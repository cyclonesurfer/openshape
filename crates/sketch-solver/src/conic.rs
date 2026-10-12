//! Conic arcs (design note 0009 §3): the rational quadratic Bézier
//! `C(t) = ((1 − t)² P₀ + 2t(1 − t) w P₁ + t² P₂) / D(t)` with
//! `D(t) = (1 − t)² + 2t(1 − t) w + t²`, `t ∈ [0, 1]`, and middle weight
//! `w = ρ/(1 − ρ)` (Lee 1987; Piegl & Tiller §7.2–7.3). The solver's shape
//! unknown is `u = ln w`, the logit of ρ (§4.1).

use glam::DVec2;

/// The solver keeps a conic's shape unknown `u` within `±SHAPE_MAX`
/// (`w` within about `e^±40`): far beyond the valid range of ρ (`|u| ≤
/// 13.8`), so that a solve heading for a flat or angular conic can be
/// seen and reported as a collapse, yet finite.
pub(crate) const SHAPE_MAX: f64 = 40.0;

/// The middle weight `w = ρ/(1 − ρ)` of shape value `ρ`.
pub fn weight(rho: f64) -> f64 {
    rho / (1.0 - rho)
}

/// The shape unknown `u = ln(ρ/(1 − ρ))` of shape value `ρ`.
pub fn logit(rho: f64) -> f64 {
    rho.ln() - (-rho).ln_1p()
}

/// The shape value `ρ = 1/(1 + e⁻ᵘ)` of shape unknown `u`.
pub fn rho_of(u: f64) -> f64 {
    1.0 / (1.0 + (-u).exp())
}

/// The rational basis `(R₀, R₁, R₂)` at `t` for middle weight `w`
/// (`C = Σ Rₖ Pₖ`, `Σ Rₖ = 1`).
pub fn basis(w: f64, t: f64) -> [f64; 3] {
    let s = 1.0 - t;
    let (b0, b1, b2) = (s * s, 2.0 * t * s * w, t * t);
    let d = b0 + b1 + b2;
    [b0 / d, b1 / d, b2 / d]
}

/// `C(t)`.
pub fn point(p: [DVec2; 3], w: f64, t: f64) -> DVec2 {
    let [r0, r1, r2] = basis(w, t);
    p[0] * r0 + p[1] * r1 + p[2] * r2
}

/// `C′(t) = (N′ − C·D′)/D`, with `N` the numerator of `C`.
pub fn tangent(p: [DVec2; 3], w: f64, t: f64) -> DVec2 {
    let s = 1.0 - t;
    let d = s * s + 2.0 * t * s * w + t * t;
    let c = point(p, w, t);
    let n1 = p[0] * (-2.0 * s) + p[1] * (2.0 * (1.0 - 2.0 * t) * w) + p[2] * (2.0 * t);
    let d1 = -2.0 * s + 2.0 * (1.0 - 2.0 * t) * w + 2.0 * t;
    (n1 - c * d1) / d
}

/// `C″(t) = (N″ − 2C′D′ − C·D″)/D`.
fn second(p: [DVec2; 3], w: f64, t: f64) -> DVec2 {
    let s = 1.0 - t;
    let d = s * s + 2.0 * t * s * w + t * t;
    let c = point(p, w, t);
    let c1 = tangent(p, w, t);
    let d1 = -2.0 * s + 2.0 * (1.0 - 2.0 * t) * w + 2.0 * t;
    let n2 = (p[0] - p[1] * (2.0 * w) + p[2]) * 2.0;
    let d2 = 4.0 - 4.0 * w;
    (n2 - c1 * (2.0 * d1) - c * d2) / d
}

/// Uniform samples of the parameter to start from.
const START: usize = 32;
/// Most samples [`samples`] takes.
const MAX_SAMPLES: usize = 2048;

/// Parameters of points along the arc spaced at most about 1/64 of the
/// sampled length apart: uniform in `t`, then intervals whose chord is
/// longer are halved. Near ρ = 1 the arc runs along its legs in a sliver
/// of `t` at each end, which uniform samples alone would step over.
pub fn samples(p: [DVec2; 3], w: f64) -> Vec<f64> {
    let mut ts: Vec<f64> = (0..=START).map(|i| i as f64 / START as f64).collect();
    let length = |ts: &[f64]| -> f64 {
        ts.windows(2)
            .map(|x| (point(p, w, x[1]) - point(p, w, x[0])).length())
            .sum()
    };
    let limit = length(&ts) / 64.0;
    loop {
        let mut out = Vec::with_capacity(ts.len() * 2);
        let mut split = false;
        for x in ts.windows(2) {
            out.push(x[0]);
            let chord = (point(p, w, x[1]) - point(p, w, x[0])).length();
            let mid = 0.5 * (x[0] + x[1]);
            if chord > limit && mid > x[0] && mid < x[1] && ts.len() + out.len() < MAX_SAMPLES {
                out.push(mid);
                split = true;
            }
        }
        out.push(1.0);
        ts = out;
        if !split {
            return ts;
        }
    }
}

/// The parameter in `[0, 1]` of the point of the arc closest to `q`.
///
/// The squared distance is taken at [`samples`]; each local minimum of the
/// samples (an end counts when its neighbour is farther) is refined by
/// golden-section search between its neighbours (Kiefer 1953), then by
/// safeguarded Newton on `(C(t) − q)·C′(t) = 0` (steps must shrink that
/// residual and stay in the bracket); the nearest wins (the first on
/// ties).
pub fn closest_param(p: [DVec2; 3], w: f64, q: DVec2) -> f64 {
    let d2 = |t: f64| (point(p, w, t) - q).length_squared();
    let ts = samples(p, w);
    let d: Vec<f64> = ts.iter().map(|&t| d2(t)).collect();
    let n = ts.len();
    let mut best = (f64::INFINITY, 0.0);
    for i in 0..n {
        let left = if i > 0 { d[i - 1] } else { f64::INFINITY };
        let right = if i + 1 < n { d[i + 1] } else { f64::INFINITY };
        if !(d[i] <= left && d[i] <= right) || (d[i] == left && d[i] == right) {
            continue;
        }
        let (lo, hi) = (ts[i.saturating_sub(1)], ts[(i + 1).min(n - 1)]);
        let (t, dt) = refine(&d2, p, w, q, (ts[i], d[i]), lo, hi);
        if dt < best.0 {
            best = (dt, t);
        }
    }
    if best.0.is_infinite() {
        return 0.0;
    }
    best.1
}

/// Refines a sampled local minimum `(t, distance²)` within `[lo, hi]`.
fn refine(
    d2: &dyn Fn(f64) -> f64,
    p: [DVec2; 3],
    w: f64,
    q: DVec2,
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
    for _ in 0..48 {
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
    // Newton polish on g(t) = (C − q)·C′, g′ = ‖C′‖² + (C − q)·C″, judged
    // by |g| (the distance is flat to rounding near its minimum).
    let g = |t: f64| {
        let e = point(p, w, t) - q;
        let c1 = tangent(p, w, t);
        (e.dot(c1), c1.length_squared() + e.dot(second(p, w, t)))
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
    // A stationary point beyond an end of the arc: the end is the nearest
    // point (golden section stops a rounding-level short of it, where the
    // distance no longer changes).
    let target = best_t - gr / gp;
    if gp > 0.0 && gr != 0.0 {
        for end in [0.0, 1.0] {
            let bracketed = if end == 0.0 {
                start_lo == 0.0
            } else {
                start_hi == 1.0
            };
            let beyond = if end == 0.0 {
                target <= 0.0
            } else {
                target >= 1.0
            };
            if bracketed && beyond {
                return (end, d2(end));
            }
        }
    }
    (best_t, best_d)
}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;

    /// Lee's implicit equation of the conic in barycentric coordinates of
    /// its control triangle: `τ₁² = 4w²τ₀τ₂`; relative residual.
    fn implicit(p: [DVec2; 3], w: f64, x: DVec2) -> f64 {
        let area = (p[1] - p[0]).perp_dot(p[2] - p[0]);
        let tau0 = (p[1] - x).perp_dot(p[2] - x) / area;
        let tau1 = (p[2] - x).perp_dot(p[0] - x) / area;
        let tau2 = (p[0] - x).perp_dot(p[1] - x) / area;
        (tau1 * tau1 - 4.0 * w * w * tau0 * tau2) / (1.0 + tau1 * tau1)
    }

    #[test]
    fn logit_and_rho_are_inverse() {
        for rho in [1e-6, 0.01, 0.3, 0.5, 0.7, 0.99, 1.0 - 1e-6] {
            let back = rho_of(logit(rho));
            assert!((back - rho).abs() <= 4.0 * f64::EPSILON * rho.max(1.0 - rho));
            assert!((weight(rho) - logit(rho).exp()).abs() <= 1e-9 * weight(rho));
        }
        assert_eq!(logit(0.5), 0.0);
    }

    #[test]
    fn the_shoulder_point_is_rho_along_the_median() {
        let p = [
            DVec2::new(0.0, 0.0),
            DVec2::new(1.0, 2.0),
            DVec2::new(3.0, 0.5),
        ];
        for rho in [0.1, 0.5, 0.8] {
            let m = (p[0] + p[2]) * 0.5;
            let s = point(p, weight(rho), 0.5);
            assert!((s - (m + (p[1] - m) * rho)).length() < 1e-14);
        }
    }

    proptest! {
        /// Points are on the conic (Lee's implicit equation) and the
        /// derivatives match central differences.
        #[test]
        fn points_and_derivatives(
            x1 in -2.0f64..2.0, y1 in 0.3f64..3.0, x2 in 1.0f64..4.0, y2 in -1.0f64..1.0,
            rho in 0.02f64..0.98, t in 0.0f64..1.0,
        ) {
            let p = [DVec2::ZERO, DVec2::new(x1, y1), DVec2::new(x2, y2)];
            prop_assume!((p[1] - p[0]).perp_dot(p[2] - p[0]).abs() > 0.1);
            let w = weight(rho);
            prop_assert!(implicit(p, w, point(p, w, t)).abs() < 1e-11);
            let h = 1e-6;
            let (lo, hi) = ((t - h).max(0.0), (t + h).min(1.0));
            let fd = (point(p, w, hi) - point(p, w, lo)) / (hi - lo);
            prop_assert!((fd - tangent(p, w, t)).length() < 1e-5 * (1.0 + fd.length()));
            let fd2 = (tangent(p, w, hi) - tangent(p, w, lo)) / (hi - lo);
            prop_assert!((fd2 - second(p, w, t)).length() < 1e-4 * (1.0 + fd2.length()));
        }

        /// The closest point is no farther than a dense sampling's nearest
        /// point, and is a foot point unless it is an end.
        #[test]
        fn closest_param_is_the_nearest_point(
            x1 in -2.0f64..2.0, y1 in 0.3f64..3.0, x2 in 1.0f64..4.0, y2 in -1.0f64..1.0,
            u in -13.0f64..13.0, qx in -3.0f64..5.0, qy in -2.0f64..4.0,
        ) {
            let p = [DVec2::ZERO, DVec2::new(x1, y1), DVec2::new(x2, y2)];
            prop_assume!((p[1] - p[0]).perp_dot(p[2] - p[0]).abs() > 0.1);
            let w = u.exp();
            let q = DVec2::new(qx, qy);
            let t = closest_param(p, w, q);
            prop_assert!((0.0..=1.0).contains(&t));
            let d = (point(p, w, t) - q).length();
            let dense = samples(p, w)
                .windows(2)
                .flat_map(|x| (0..64).map(move |k| x[0] + (x[1] - x[0]) * k as f64 / 64.0))
                .chain([1.0])
                .map(|s| (point(p, w, s) - q).length())
                .fold(f64::INFINITY, f64::min);
            prop_assert!(d <= dense + 1e-12 * (1.0 + dense), "{d} vs {dense}");
            if t > 0.0 && t < 1.0 {
                let tan = tangent(p, w, t);
                prop_assert!(
                    (point(p, w, t) - q).dot(tan).abs() <= 1e-8 * tan.length() * (1.0 + d)
                );
            }
        }
    }

    #[test]
    fn samples_follow_the_legs_of_a_nearly_angular_conic() {
        let p = [DVec2::ZERO, DVec2::new(0.0, 1.0), DVec2::new(1.0, 1.0)];
        let w = weight(1.0 - 1e-6);
        let ts = samples(p, w);
        let pts: Vec<DVec2> = ts.iter().map(|&t| point(p, w, t)).collect();
        let longest = pts
            .windows(2)
            .map(|x| (x[1] - x[0]).length())
            .fold(0.0, f64::max);
        assert!(longest < 0.04, "{longest}");
        // A point beside the middle of the first leg finds it.
        let t = closest_param(p, w, DVec2::new(0.1, 0.5));
        assert!((point(p, w, t) - DVec2::new(0.0, 0.5)).length() < 1e-5);
    }
}
