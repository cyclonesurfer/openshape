//! Exact measures on closed loops of curve pieces: signed area (Green's
//! theorem), winding numbers (ray crossing on rational Bézier pieces with
//! exact orientation tests), distance to the boundary, and interior points.

use std::f64::consts::FRAC_PI_4;

use glam::DVec2;

use super::bezier::{self, gauss_legendre, Bez};
use super::curve::{Crv, Geom};

/// A piece `[t0, t1]` of curve `curve` traversed forwards (`fwd`) or
/// backwards.
#[derive(Clone, Copy, Debug, PartialEq)]
pub(super) struct Side {
    /// Index into the curve list.
    pub curve: usize,
    /// Parameter interval start.
    pub t0: f64,
    /// Parameter interval end (`t0 < t1`).
    pub t1: f64,
    /// Traversed from `t0` to `t1`.
    pub fwd: bool,
}

impl Side {
    /// First point in the direction of travel.
    pub fn start(&self, crvs: &[Crv]) -> DVec2 {
        crvs[self.curve].point(if self.fwd { self.t0 } else { self.t1 })
    }

    /// Last point in the direction of travel.
    pub fn end(&self, crvs: &[Crv]) -> DVec2 {
        crvs[self.curve].point(if self.fwd { self.t1 } else { self.t0 })
    }
}

/// Exact orientation of `c` relative to the directed line `ab` (Shewchuk).
pub(super) fn orient(a: DVec2, b: DVec2, c: DVec2) -> f64 {
    robust::orient2d(
        robust::Coord { x: a.x, y: a.y },
        robust::Coord { x: b.x, y: b.y },
        robust::Coord { x: c.x, y: c.y },
    )
}

/// `∫ (x dy − y dx) / 2` over the curve piece `[t0, t1]`, forwards, with
/// coordinates relative to `o`. Lines and arcs in closed form; splines by
/// Gauss–Legendre quadrature per knot span (exact for polynomial spans,
/// adaptive for rational ones).
pub(super) fn piece_area(c: &Crv, t0: f64, t1: f64, o: DVec2) -> f64 {
    match &c.geom {
        Geom::Line { .. } => 0.5 * (c.point(t0) - o).perp_dot(c.point(t1) - o),
        Geom::Arc { c: ctr, r, sgn, .. } => {
            let p0 = c.point(t0);
            let p1 = c.point(t1);
            0.5 * ((*ctr - o).perp_dot(p1 - p0) + r * r * sgn * (t1 - t0))
        }
        Geom::Spline { .. } => c
            .pieces(t0, t1)
            .iter()
            .map(|p| bez_area(&p.bez, o))
            .sum(),
    }
}

fn bez_area(b: &Bez, o: DVec2) -> f64 {
    let f = |u: f64| {
        let [p, d1, _] = b.derivs(u);
        0.5 * (p - o).perp_dot(d1)
    };
    let gauss = |n: usize, a: f64, bb: f64| -> f64 {
        let h = 0.5 * (bb - a);
        gauss_legendre(n)
            .iter()
            .map(|&(x, w)| w * f(a + h * (x + 1.0)))
            .sum::<f64>()
            * h
    };
    if b.is_polynomial() {
        return gauss(b.degree() + 1, 0.0, 1.0);
    }
    let n = 2 * b.degree() + 8;
    fn adapt(g: &dyn Fn(usize, f64, f64) -> f64, n: usize, a: f64, b: f64, whole: f64, depth: u32) -> f64 {
        let m = 0.5 * (a + b);
        let l = g(n, a, m);
        let r = g(n, m, b);
        let scale = l.abs() + r.abs();
        if depth >= 24 || (l + r - whole).abs() <= 1e-15 * scale.max(f64::MIN_POSITIVE) {
            l + r
        } else {
            adapt(g, n, a, m, l, depth + 1) + adapt(g, n, m, b, r, depth + 1)
        }
    }
    let whole = gauss(n, 0.0, 1.0);
    adapt(&gauss, n, 0.0, 1.0, whole, 0)
}

/// Signed area of a closed loop (positive when counter-clockwise), with
/// straight connectors closing any rounding gaps between consecutive pieces.
pub(super) fn loop_area(crvs: &[Crv], lp: &[Side]) -> f64 {
    if lp.is_empty() {
        return 0.0;
    }
    let o = lp[0].start(crvs);
    let mut a = 0.0;
    for (k, s) in lp.iter().enumerate() {
        let c = &crvs[s.curve];
        let v = piece_area(c, s.t0, s.t1, o);
        a += if s.fwd { v } else { -v };
        let e = s.end(crvs);
        let n = lp[(k + 1) % lp.len()].start(crvs);
        a += 0.5 * (e - o).perp_dot(n - o);
    }
    a
}

/// The loop as a chain of oriented rational Bézier pieces.
pub(super) fn loop_bez(crvs: &[Crv], lp: &[Side]) -> Vec<Bez> {
    let mut out = Vec::new();
    for s in lp {
        let ps = crvs[s.curve].pieces(s.t0, s.t1);
        if s.fwd {
            out.extend(ps.into_iter().map(|p| p.bez));
        } else {
            out.extend(ps.into_iter().rev().map(|p| p.bez.reversed()));
        }
    }
    out
}

/// Net signed crossings of segment `ab` with the ray from `p` towards +x
/// (half-open rule: "above" means `y > p.y`). `None` if `p` is on `ab`.
fn seg_crossings(a: DVec2, b: DVec2, p: DVec2) -> Option<i32> {
    let ua = a.y > p.y;
    let ub = b.y > p.y;
    if ua == ub {
        return Some(0);
    }
    let o = orient(a, b, p);
    if o == 0.0 {
        return None;
    }
    if ub {
        Some(i32::from(o > 0.0))
    } else {
        Some(-i32::from(o < 0.0))
    }
}

/// Net signed crossings of a rational Bézier piece with the ray from `p`
/// towards +x, by subdivision until the convex hull is decided. `None` if
/// `p` is (numerically) on the curve.
fn bez_crossings(b: &Bez, p: DVec2, depth: u32) -> Option<i32> {
    if b.degree() == 1 {
        return seg_crossings(b.start(), b.end(), p);
    }
    let hull = b.hull();
    if hull.iter().all(|q| q.y > p.y) || hull.iter().all(|q| q.y <= p.y) {
        return Some(0);
    }
    if hull.iter().all(|q| q.x < p.x) {
        return Some(0);
    }
    if hull.iter().all(|q| q.x > p.x) {
        return Some(i32::from(b.end().y > p.y) - i32::from(b.start().y > p.y));
    }
    if depth >= 90 {
        return None;
    }
    let (l, r) = b.split(0.5);
    Some(bez_crossings(&l, p, depth + 1)? + bez_crossings(&r, p, depth + 1)?)
}

/// Winding number of a closed chain of Bézier pieces around `p`; straight
/// connectors close any gaps. `None` if `p` lies on the chain.
pub(super) fn chain_winding(chain: &[Bez], p: DVec2) -> Option<i32> {
    let mut w = 0;
    for (k, b) in chain.iter().enumerate() {
        w += bez_crossings(b, p, 0)?;
        let e = b.end();
        let n = chain[(k + 1) % chain.len()].start();
        if e != n {
            w += seg_crossings(e, n, p)?;
        }
    }
    Some(w)
}

/// Winding number of a loop around `p` (`None` if `p` is on it).
pub(super) fn winding(crvs: &[Crv], lp: &[Side], p: DVec2) -> Option<i32> {
    chain_winding(&loop_bez(crvs, lp), p)
}

/// Distance from `p` to a loop.
pub(super) fn loop_distance(crvs: &[Crv], lp: &[Side], p: DVec2) -> f64 {
    let mut d = f64::INFINITY;
    for (k, s) in lp.iter().enumerate() {
        d = d.min(crvs[s.curve].closest(p, s.t0, s.t1).1);
        let e = s.end(crvs);
        let n = lp[(k + 1) % lp.len()].start(crvs);
        if e != n {
            d = d.min(bezier::point_segment_distance(p, e, n));
        }
    }
    d
}

/// Bounds of a loop (of its Bézier control points).
pub(super) fn loop_bounds(crvs: &[Crv], lp: &[Side]) -> (DVec2, DVec2) {
    let mut lo = DVec2::splat(f64::INFINITY);
    let mut hi = DVec2::splat(f64::NEG_INFINITY);
    for b in loop_bez(crvs, lp) {
        let (l, h) = bezier::bounds(&b.hull());
        lo = lo.min(l);
        hi = hi.max(h);
    }
    (lo, hi)
}

/// A region given by its loops (outer first, then holes, all oriented with
/// the region on the left), prepared for point queries.
pub(super) struct RegionShape<'a> {
    crvs: &'a [Crv],
    loops: Vec<Vec<Side>>,
    chains: Vec<Vec<Bez>>,
    /// Bounds of the outer loop.
    pub lo: DVec2,
    /// Bounds of the outer loop.
    pub hi: DVec2,
}

impl<'a> RegionShape<'a> {
    /// Prepares the loops of a region.
    pub fn new(crvs: &'a [Crv], loops: Vec<Vec<Side>>) -> RegionShape<'a> {
        let chains: Vec<Vec<Bez>> = loops.iter().map(|l| loop_bez(crvs, l)).collect();
        let (lo, hi) = loops
            .first()
            .map(|l| loop_bounds(crvs, l))
            .unwrap_or((DVec2::ZERO, DVec2::ZERO));
        RegionShape {
            crvs,
            loops,
            chains,
            lo,
            hi,
        }
    }

    /// Sum of the loops' winding numbers around `p`: 1 inside the region,
    /// 0 outside. `None` if `p` is on the boundary.
    pub fn winding(&self, p: DVec2) -> Option<i32> {
        if p.x < self.lo.x || p.y < self.lo.y || p.x > self.hi.x || p.y > self.hi.y {
            return Some(0);
        }
        let mut w = 0;
        for c in &self.chains {
            w += chain_winding(c, p)?;
        }
        Some(w)
    }

    /// Distance from `p` to the region's boundary.
    pub fn clearance(&self, p: DVec2) -> f64 {
        self.loops
            .iter()
            .map(|l| loop_distance(self.crvs, l, p))
            .fold(f64::INFINITY, f64::min)
    }

    /// True if `p` is inside and farther than `tol` from the boundary.
    pub fn contains(&self, p: DVec2, tol: f64) -> bool {
        self.winding(p) == Some(1) && self.clearance(p) > tol
    }

    /// A point inside the region, as far from the boundary as this search
    /// finds, with its clearance (design note 0006, §4.8). Candidates are the
    /// midpoints of the inside intervals of horizontal and vertical scan
    /// lines across a fine polygonal approximation; each candidate is checked
    /// with the exact winding number and its exact clearance, and the best is
    /// improved by a pattern search on the clearance.
    pub fn interior_point(&self, tol: f64) -> Option<(DVec2, f64)> {
        let polys: Vec<Vec<DVec2>> = self
            .chains
            .iter()
            .map(|chain| {
                let mut pts = Vec::new();
                for b in chain {
                    let n = if b.degree() == 1 { 1 } else { 24 };
                    for k in 0..n {
                        pts.push(b.point(k as f64 / n as f64));
                    }
                }
                pts
            })
            .collect();
        let mut best: Option<(DVec2, f64)> = None;
        let consider = |p: DVec2, best: &mut Option<(DVec2, f64)>| {
            if !p.is_finite() || self.winding(p) != Some(1) {
                return;
            }
            let d = self.clearance(p);
            if best.is_none_or(|(_, bd)| d > bd) {
                *best = Some((p, d));
            }
        };
        for k in [16usize, 64, 256] {
            for axis in 0..2 {
                let o = 1 - axis;
                for i in 0..k {
                    let c = self.lo[axis] + (i as f64 + 0.5) / k as f64 * (self.hi[axis] - self.lo[axis]);
                    let mut xs = Vec::new();
                    for pl in &polys {
                        for j in 0..pl.len() {
                            let a = pl[j];
                            let b = pl[(j + 1) % pl.len()];
                            if (a[axis] > c) != (b[axis] > c) {
                                let f = (c - a[axis]) / (b[axis] - a[axis]);
                                xs.push(a[o] + f * (b[o] - a[o]));
                            }
                        }
                    }
                    xs.sort_by(f64::total_cmp);
                    for pair in xs.chunks_exact(2) {
                        let m = 0.5 * (pair[0] + pair[1]);
                        let mut p = DVec2::ZERO;
                        p[axis] = c;
                        p[o] = m;
                        consider(p, &mut best);
                    }
                }
            }
            if best.is_some_and(|(_, d)| d >= 4.0 * tol) {
                break;
            }
        }
        let (mut p, mut d) = best?;
        let mut step = 0.5 * d;
        for _ in 0..60 {
            let mut moved = false;
            for k in 0..8 {
                let q = p + step * DVec2::from_angle(k as f64 * FRAC_PI_4);
                if self.winding(q) != Some(1) {
                    continue;
                }
                let dq = self.clearance(q);
                if dq > d {
                    p = q;
                    d = dq;
                    moved = true;
                    break;
                }
            }
            if !moved {
                step *= 0.5;
                if step < 0.01 * tol || step < 1e-6 * d {
                    break;
                }
            }
        }
        Some((p, d))
    }
}
