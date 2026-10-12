//! Shared proptest strategies and helpers for the integration tests.

#![allow(dead_code)]

use glam::DVec3;
use opendatum_nurbs::{NurbsCurve, NurbsSurface};
use proptest::prelude::*;

/// Clamped knot vector for `n_ctrl` control points of degree `p` whose spans
/// have the given relative lengths, mapped onto `[start, start + width]`.
pub fn clamped_knots(p: usize, n_ctrl: usize, spans: &[f64], start: f64, width: f64) -> Vec<f64> {
    assert_eq!(spans.len(), n_ctrl - p);
    let total: f64 = spans.iter().sum();
    let mut knots = vec![start; p + 1];
    let mut acc = 0.0;
    for s in &spans[..spans.len() - 1] {
        acc += s;
        knots.push(start + width * acc / total);
    }
    knots.extend(std::iter::repeat_n(start + width, p + 1));
    knots
}

fn arb_point() -> impl Strategy<Value = DVec3> {
    (-10.0..10.0f64, -10.0..10.0f64, -10.0..10.0f64).prop_map(|(x, y, z)| DVec3::new(x, y, z))
}

/// Random rational curve: degree 1..=5, up to 6 extra control points,
/// span lengths within a factor of 5, weights in `[0.2, 5]`, coordinates in
/// `[-10, 10]` and a shifted, scaled domain.
pub fn arb_curve() -> impl Strategy<Value = NurbsCurve> {
    (1usize..=5)
        .prop_flat_map(|p| (Just(p), (p + 1)..=(p + 6)))
        .prop_flat_map(|(p, n)| {
            (
                Just(p),
                prop::collection::vec(0.2..1.0f64, n - p),
                prop::collection::vec(arb_point(), n),
                prop::collection::vec(0.2..5.0f64, n),
                -5.0..5.0f64,
                0.5..10.0f64,
            )
        })
        .prop_map(|(p, spans, points, weights, start, width)| {
            let knots = clamped_knots(p, points.len(), &spans, start, width);
            NurbsCurve::new(p, knots, points, Some(weights)).expect("valid random curve")
        })
}

/// Random rational surface with a perturbed-grid control net (a "height
/// field" with in-plane jitter), degrees 1..=3 in each direction.
pub fn arb_surface() -> impl Strategy<Value = NurbsSurface> {
    ((1usize..=3), (1usize..=3))
        .prop_flat_map(|(p, q)| (Just(p), Just(q), (p + 1)..=(p + 3), (q + 1)..=(q + 3)))
        .prop_flat_map(|(p, q, nu, nv)| {
            (
                Just(p),
                Just(q),
                prop::collection::vec(0.3..1.0f64, nu - p),
                prop::collection::vec(0.3..1.0f64, nv - q),
                prop::collection::vec(
                    (-0.3..0.3f64, -0.3..0.3f64, -2.0..2.0f64, 0.3..3.0f64),
                    nu * nv,
                ),
                Just((nu, nv)),
            )
        })
        .prop_map(|(p, q, su, sv, jitter, (nu, nv))| {
            let ku = clamped_knots(p, nu, &su, 0.0, 1.0);
            let kv = clamped_knots(q, nv, &sv, -1.0, 3.0);
            let mut points = Vec::new();
            let mut weights = Vec::new();
            for i in 0..nu {
                let mut row = Vec::new();
                let mut wrow = Vec::new();
                for j in 0..nv {
                    let (dx, dy, z, w) = jitter[i * nv + j];
                    row.push(DVec3::new(
                        4.0 * i as f64 / (nu - 1) as f64 + dx,
                        4.0 * j as f64 / (nv - 1) as f64 + dy,
                        z,
                    ));
                    wrow.push(w);
                }
                points.push(row);
                weights.push(wrow);
            }
            NurbsSurface::new(p, q, ku, kv, points, Some(weights)).expect("valid random surface")
        })
}

/// Parameters strictly inside each non-empty knot span of the domain, at the
/// given relative positions.
pub fn interior_params(knots: &[f64], p: usize, rel: &[f64]) -> Vec<f64> {
    let m = knots.len() - 1;
    let (a, b) = (knots[p], knots[m - p]);
    let mut breaks = vec![a];
    for &k in knots {
        if k > *breaks.last().unwrap() && k < b {
            breaks.push(k);
        }
    }
    breaks.push(b);
    let mut out = Vec::new();
    for w in breaks.windows(2) {
        for &r in rel {
            out.push(w[0] + (w[1] - w[0]) * r);
        }
    }
    out
}

/// Width of the knot span containing `t` (domain interior).
pub fn span_width(knots: &[f64], t: f64) -> f64 {
    let lo = knots
        .iter()
        .copied()
        .filter(|&k| k <= t)
        .fold(f64::MIN, f64::max);
    let hi = knots
        .iter()
        .copied()
        .filter(|&k| k > t)
        .fold(f64::MAX, f64::min);
    hi - lo
}

/// Relative closeness check with a floor of 1 on the scale.
pub fn close(a: DVec3, b: DVec3, rel: f64) -> bool {
    (a - b).length() <= rel * a.length().max(b.length()).max(1.0)
}
