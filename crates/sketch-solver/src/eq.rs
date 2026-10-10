//! Residual equations (design note 0005 §4) with analytic gradients.
//!
//! Unknowns live in a flat slot vector: point `p` has slots `2p` (x) and
//! `2p + 1` (y); circle radii and hidden spline parameters follow. Every
//! equation evaluates its residual and, on request, pushes `(slot, ∂r/∂slot)`
//! pairs (a slot may appear several times; the caller sums them).
//! Equations that divide by a length report degeneracy instead of producing
//! a non-finite value.

use crate::spline::SplineMap;

/// `Σ c_k P_k` over point indices: a 2D vector built linearly from points.
#[derive(Clone, Debug, PartialEq)]
pub(crate) struct VecExpr(pub Vec<(usize, f64)>);

impl VecExpr {
    /// `P_a − P_b`.
    pub fn diff(a: usize, b: usize) -> Self {
        VecExpr(vec![(a, 1.0), (b, -1.0)])
    }

    /// The vector's value.
    pub fn eval(&self, v: &[f64]) -> (f64, f64) {
        let mut x = 0.0;
        let mut y = 0.0;
        for &(p, c) in &self.0 {
            x += c * v[2 * p];
            y += c * v[2 * p + 1];
        }
        (x, y)
    }

    /// Pushes the gradient of `f(V)` given `(∂f/∂Vx, ∂f/∂Vy)`.
    pub fn push_grad(&self, gx: f64, gy: f64, out: &mut Vec<(usize, f64)>) {
        for &(p, c) in &self.0 {
            out.push((2 * p, c * gx));
            out.push((2 * p + 1, c * gy));
        }
    }

    /// The points referenced.
    pub fn points(&self) -> impl Iterator<Item = usize> + '_ {
        self.0.iter().map(|&(p, _)| p)
    }
}

/// A scalar length term of a right-hand side.
#[derive(Clone, Debug, PartialEq)]
pub(crate) enum RTerm {
    /// A constant (scaled).
    Const(f64),
    /// A slot (a circle radius).
    Slot(usize),
    /// The length of a vector (an arc radius `‖start − centre‖`).
    Norm(VecExpr),
}

/// `Σ cᵢ termᵢ`, or its absolute value.
#[derive(Clone, Debug, PartialEq)]
pub(crate) struct Rhs {
    /// Coefficients and terms.
    pub terms: Vec<(f64, RTerm)>,
    /// Take the absolute value of the sum.
    pub abs: bool,
}

impl Rhs {
    /// The constant `c`.
    pub fn constant(c: f64) -> Self {
        Rhs {
            terms: vec![(1.0, RTerm::Const(c))],
            abs: false,
        }
    }

    /// Zero.
    pub fn zero() -> Self {
        Rhs {
            terms: Vec::new(),
            abs: false,
        }
    }

    /// A single term with coefficient `c`.
    pub fn term(c: f64, t: RTerm) -> Self {
        Rhs {
            terms: vec![(c, t)],
            abs: false,
        }
    }

    fn raw(&self, v: &[f64]) -> f64 {
        let mut s = 0.0;
        for (c, t) in &self.terms {
            s += c * match t {
                RTerm::Const(k) => *k,
                RTerm::Slot(i) => v[*i],
                RTerm::Norm(e) => {
                    let (x, y) = e.eval(v);
                    norm(x, y)
                }
            };
        }
        s
    }

    /// The value.
    pub fn eval(&self, v: &[f64]) -> f64 {
        let s = self.raw(v);
        if self.abs {
            s.abs()
        } else {
            s
        }
    }

    /// Pushes `factor · ∂value/∂slot`.
    pub fn push_grad(&self, v: &[f64], factor: f64, out: &mut Vec<(usize, f64)>) {
        let f = if self.abs && self.raw(v) < 0.0 {
            -factor
        } else {
            factor
        };
        for (c, t) in &self.terms {
            match t {
                RTerm::Const(_) => {}
                RTerm::Slot(i) => out.push((*i, f * c)),
                RTerm::Norm(e) => {
                    let (x, y) = e.eval(v);
                    let (ux, uy) = unit(x, y);
                    e.push_grad(f * c * ux, f * c * uy, out);
                }
            }
        }
    }

    /// Slots and points referenced: `(slots, points)`.
    pub fn refs(&self, slots: &mut Vec<usize>, points: &mut Vec<usize>) {
        for (_, t) in &self.terms {
            match t {
                RTerm::Const(_) => {}
                RTerm::Slot(i) => slots.push(*i),
                RTerm::Norm(e) => points.extend(e.points()),
            }
        }
    }
}

/// How a pair of directions is compared, for [`EqKind::DirPair`].
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) enum PairMode {
    /// `cross(a, b)/(‖a‖‖b‖)`: parallel (or both normals parallel).
    Cross,
    /// `dot(a, b)/(‖a‖‖b‖)`: perpendicular.
    Dot,
    /// `sin(θ − v)` with θ the signed angle from `a` to `b`.
    AngleSin {
        /// The angle value `v`.
        value: f64,
    },
}

/// The kinds of residual equations.
#[derive(Clone, Debug, PartialEq)]
pub(crate) enum EqKind {
    /// `Σ c·slot − rhs`.
    Lin {
        /// Slots and coefficients.
        terms: Vec<(usize, f64)>,
        /// Constant.
        rhs: f64,
    },
    /// `|Σ c·slot| − rhs`.
    AbsLin {
        /// Slots and coefficients.
        terms: Vec<(usize, f64)>,
        /// Constant.
        rhs: f64,
    },
    /// `k‖v‖ − rhs`.
    Norm {
        /// Factor.
        k: f64,
        /// Vector.
        v: VecExpr,
        /// Right-hand side.
        rhs: Rhs,
    },
    /// `cross(dir, v)/‖dir‖ − rhs` (absolute value of the first term if
    /// `abs`): signed distance of a point from a line.
    SignedDist {
        /// Line direction.
        dir: VecExpr,
        /// Point minus a point of the line.
        v: VecExpr,
        /// Unsigned distance.
        abs: bool,
        /// Right-hand side.
        rhs: Rhs,
    },
    /// `dot(dir, v)/‖dir‖`.
    Proj {
        /// Direction.
        dir: VecExpr,
        /// Vector.
        v: VecExpr,
    },
    /// Normalized comparison of two directions.
    DirPair {
        /// First direction.
        a: VecExpr,
        /// Second direction.
        b: VecExpr,
        /// Comparison.
        mode: PairMode,
    },
    /// `C(t)_axis − p_axis` for a point on a spline.
    SplinePoint {
        /// Index into the problem's splines.
        spline: usize,
        /// Slot of the hidden parameter `t`.
        t: usize,
        /// The point.
        p: usize,
        /// 0 for x, 1 for y.
        axis: usize,
    },
}

/// A spline as used by equations: its linear map and its points.
#[derive(Clone, Debug)]
pub(crate) struct SplineDef {
    /// The map from points to curve.
    pub map: SplineMap,
    /// Point indices, in the spline's order.
    pub points: Vec<usize>,
}

/// Result of evaluating an equation.
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct Eval {
    /// Residual.
    pub r: f64,
    /// For direction comparisons, the shorter of the two direction lengths
    /// (the lever that turns an angle error into a position error);
    /// infinity otherwise.
    pub lever: f64,
}

/// `‖(x, y)‖` without overflow for the magnitudes the solver uses.
#[inline]
pub(crate) fn norm(x: f64, y: f64) -> f64 {
    (x * x + y * y).sqrt()
}

/// The unit vector, or `(1, 0)` for the zero vector (a fixed, deterministic
/// choice of sub-gradient for `‖v‖` at `v = 0`).
#[inline]
pub(crate) fn unit(x: f64, y: f64) -> (f64, f64) {
    let n = norm(x, y);
    if n > 0.0 {
        (x / n, y / n)
    } else {
        (1.0, 0.0)
    }
}

impl EqKind {
    /// Evaluates the residual at slot values `v`, pushing the gradient into
    /// `grad` if given. Returns `None` if a direction it divides by is
    /// shorter than `degen`.
    pub fn eval(
        &self,
        v: &[f64],
        splines: &[SplineDef],
        degen: f64,
        grad: Option<&mut Vec<(usize, f64)>>,
    ) -> Option<Eval> {
        let inf = f64::INFINITY;
        match self {
            EqKind::Lin { terms, rhs } => {
                let mut s = 0.0;
                for &(i, c) in terms {
                    s += c * v[i];
                }
                if let Some(g) = grad {
                    g.extend(terms.iter().copied());
                }
                Some(Eval {
                    r: s - rhs,
                    lever: inf,
                })
            }
            EqKind::AbsLin { terms, rhs } => {
                let mut s = 0.0;
                for &(i, c) in terms {
                    s += c * v[i];
                }
                let sg = if s >= 0.0 { 1.0 } else { -1.0 };
                if let Some(g) = grad {
                    g.extend(terms.iter().map(|&(i, c)| (i, sg * c)));
                }
                Some(Eval {
                    r: s.abs() - rhs,
                    lever: inf,
                })
            }
            EqKind::Norm { k, v: e, rhs } => {
                let (x, y) = e.eval(v);
                let n = norm(x, y);
                if let Some(g) = grad {
                    let (ux, uy) = unit(x, y);
                    e.push_grad(k * ux, k * uy, g);
                    rhs.push_grad(v, -1.0, g);
                }
                Some(Eval {
                    r: k * n - rhs.eval(v),
                    lever: inf,
                })
            }
            EqKind::SignedDist { dir, v: w, abs, rhs } => {
                let (dx, dy) = dir.eval(v);
                let n = norm(dx, dy);
                if n.is_nan() || n <= degen {
                    return None;
                }
                let (wx, wy) = w.eval(v);
                let c = dx * wy - dy * wx;
                let s = c / n;
                let sg = if *abs && s < 0.0 { -1.0 } else { 1.0 };
                if let Some(g) = grad {
                    let n3 = n * n * n;
                    w.push_grad(sg * (-dy / n), sg * (dx / n), g);
                    dir.push_grad(
                        sg * (wy / n - c * dx / n3),
                        sg * (-wx / n - c * dy / n3),
                        g,
                    );
                    rhs.push_grad(v, -1.0, g);
                }
                Some(Eval {
                    r: sg * s - rhs.eval(v),
                    lever: inf,
                })
            }
            EqKind::Proj { dir, v: w } => {
                let (dx, dy) = dir.eval(v);
                let n = norm(dx, dy);
                if n.is_nan() || n <= degen {
                    return None;
                }
                let (wx, wy) = w.eval(v);
                let s = (dx * wx + dy * wy) / n;
                if let Some(g) = grad {
                    w.push_grad(dx / n, dy / n, g);
                    let n2 = n * n;
                    dir.push_grad(wx / n - s * dx / n2, wy / n - s * dy / n2, g);
                }
                Some(Eval { r: s, lever: inf })
            }
            EqKind::DirPair { a, b, mode } => {
                let (ax, ay) = a.eval(v);
                let (bx, by) = b.eval(v);
                let na = norm(ax, ay);
                let nb = norm(bx, by);
                if na.is_nan() || nb.is_nan() || na <= degen || nb <= degen {
                    return None;
                }
                let nn = na * nb;
                let cr = ax * by - ay * bx;
                let dt = ax * bx + ay * by;
                let (cc, cd) = match mode {
                    PairMode::Cross => (1.0, 0.0),
                    PairMode::Dot => (0.0, 1.0),
                    PairMode::AngleSin { value } => (value.cos(), -value.sin()),
                };
                // r = (cc·cross + cd·dot)/(‖a‖‖b‖).
                let r = (cc * cr + cd * dt) / nn;
                if let Some(g) = grad {
                    let na2 = na * na;
                    let nb2 = nb * nb;
                    // ∂cross/∂a = (by, −bx), ∂dot/∂a = (bx, by); the
                    // normalization contributes −r·a/‖a‖².
                    let gax = (cc * by + cd * bx) / nn - r * ax / na2;
                    let gay = (-cc * bx + cd * by) / nn - r * ay / na2;
                    let gbx = (-cc * ay + cd * ax) / nn - r * bx / nb2;
                    let gby = (cc * ax + cd * ay) / nn - r * by / nb2;
                    a.push_grad(gax, gay, g);
                    b.push_grad(gbx, gby, g);
                }
                Some(Eval {
                    r,
                    lever: na.min(nb),
                })
            }
            EqKind::SplinePoint {
                spline,
                t,
                p,
                axis,
            } => {
                let def = &splines[*spline];
                let (w, dw) = def.map.weights(v[*t]);
                let mut c = 0.0;
                let mut dc = 0.0;
                for (k, &q) in def.points.iter().enumerate() {
                    c += w[k] * v[2 * q + axis];
                    dc += dw[k] * v[2 * q + axis];
                }
                if let Some(g) = grad {
                    g.push((*t, dc));
                    for (k, &q) in def.points.iter().enumerate() {
                        if w[k] != 0.0 {
                            g.push((2 * q + axis, w[k]));
                        }
                    }
                    g.push((2 * p + axis, -1.0));
                }
                Some(Eval {
                    r: c - v[2 * p + axis],
                    lever: inf,
                })
            }
        }
    }

    /// For an angle dimension, `cos(θ − v)` (positive on the intended
    /// branch); `None` for other equations or degenerate directions.
    pub fn angle_cos(&self, v: &[f64], degen: f64) -> Option<f64> {
        if let EqKind::DirPair {
            a,
            b,
            mode: PairMode::AngleSin { value },
        } = self
        {
            let (ax, ay) = a.eval(v);
            let (bx, by) = b.eval(v);
            let na = norm(ax, ay);
            let nb = norm(bx, by);
            if na <= degen || nb <= degen {
                return None;
            }
            let cr = ax * by - ay * bx;
            let dt = ax * bx + ay * by;
            Some((dt * value.cos() + cr * value.sin()) / (na * nb))
        } else {
            None
        }
    }

    /// The signed angle `θ ∈ (−π, π]` from `a` to `b` of an angle equation.
    pub fn angle_theta(&self, v: &[f64]) -> Option<f64> {
        if let EqKind::DirPair {
            a,
            b,
            mode: PairMode::AngleSin { .. },
        } = self
        {
            let (ax, ay) = a.eval(v);
            let (bx, by) = b.eval(v);
            Some(signed_angle(ax * by - ay * bx, ax * bx + ay * by))
        } else {
            None
        }
    }

    /// Slots this equation depends on (with repetitions).
    pub fn slots(&self, splines: &[SplineDef]) -> Vec<usize> {
        let mut slots = Vec::new();
        let mut points = Vec::new();
        match self {
            EqKind::Lin { terms, .. } | EqKind::AbsLin { terms, .. } => {
                slots.extend(terms.iter().map(|&(i, _)| i));
            }
            EqKind::Norm { v, rhs, .. } => {
                points.extend(v.points());
                rhs.refs(&mut slots, &mut points);
            }
            EqKind::SignedDist { dir, v, rhs, .. } => {
                points.extend(dir.points());
                points.extend(v.points());
                rhs.refs(&mut slots, &mut points);
            }
            EqKind::Proj { dir, v } => {
                points.extend(dir.points());
                points.extend(v.points());
            }
            EqKind::DirPair { a, b, .. } => {
                points.extend(a.points());
                points.extend(b.points());
            }
            EqKind::SplinePoint {
                spline, t, p, axis, ..
            } => {
                slots.push(*t);
                slots.push(2 * p + axis);
                for &q in &splines[*spline].points {
                    slots.push(2 * q + axis);
                }
            }
        }
        for p in points {
            slots.push(2 * p);
            slots.push(2 * p + 1);
        }
        slots
    }
}

/// `atan2(cross, dot)` mapped to `(−π, π]`.
pub(crate) fn signed_angle(cross: f64, dot: f64) -> f64 {
    let a = cross.atan2(dot);
    if a <= -std::f64::consts::PI {
        std::f64::consts::PI
    } else {
        a
    }
}

#[cfg(test)]
mod tests {
    //! Finite-difference checks of every equation's gradient (design note
    //! 0005 §8) at random configurations.

    use super::*;
    use proptest::prelude::*;

    /// Six points (12 slots), a radius slot (12), a radius slot (13), a
    /// spline parameter slot (14).
    const SLOTS: usize = 15;

    fn splines() -> Vec<SplineDef> {
        vec![
            SplineDef {
                map: SplineMap::control(4, 3).unwrap(),
                points: vec![1, 2, 3, 4],
            },
            SplineDef {
                map: SplineMap::interpolated(&[0.0, 0.3, 0.55, 1.0]).unwrap(),
                points: vec![2, 5, 3, 1],
            },
            SplineDef {
                map: SplineMap::interpolated(&[0.0, 0.4, 1.0]).unwrap(),
                points: vec![0, 2, 4],
            },
            SplineDef {
                map: SplineMap::control(7, 3).unwrap(),
                points: vec![0, 1, 2, 3, 4, 5, 1],
            },
            SplineDef {
                map: SplineMap::interpolated(&[0.0, 0.1, 0.35, 0.5, 0.8, 1.0]).unwrap(),
                points: vec![5, 4, 3, 2, 1, 0],
            },
        ]
    }

    /// One equation of every kind and variant the builder produces.
    fn all_equations() -> Vec<(&'static str, EqKind)> {
        let d01 = VecExpr::diff(1, 0);
        let d23 = VecExpr::diff(3, 2);
        let arc_r = RTerm::Norm(VecExpr::diff(4, 5));
        let mid = VecExpr(vec![(2, 0.5), (3, 0.5), (0, -1.0)]);
        let tangent_end = VecExpr(vec![(2, -0.7), (5, 1.3), (3, -0.4), (1, -0.2)]);
        vec![
            (
                "lin",
                EqKind::Lin {
                    terms: vec![(0, 1.0), (2, -0.5), (4, -0.5)],
                    rhs: 0.25,
                },
            ),
            (
                "abs lin",
                EqKind::AbsLin {
                    terms: vec![(1, 1.0), (3, -1.0)],
                    rhs: 0.3,
                },
            ),
            (
                "distance",
                EqKind::Norm {
                    k: 1.0,
                    v: d01.clone(),
                    rhs: Rhs::constant(0.7),
                },
            ),
            (
                "point on circle",
                EqKind::Norm {
                    k: 1.0,
                    v: VecExpr::diff(2, 3),
                    rhs: Rhs::term(1.0, RTerm::Slot(12)),
                },
            ),
            (
                "point on arc / equal lengths",
                EqKind::Norm {
                    k: 1.0,
                    v: VecExpr::diff(2, 3),
                    rhs: Rhs::term(1.0, arc_r.clone()),
                },
            ),
            (
                "diameter of arc",
                EqKind::Norm {
                    k: 2.0,
                    v: VecExpr::diff(4, 5),
                    rhs: Rhs::constant(1.1),
                },
            ),
            (
                "external circle tangency",
                EqKind::Norm {
                    k: 1.0,
                    v: VecExpr::diff(0, 3),
                    rhs: Rhs {
                        terms: vec![(1.0, RTerm::Slot(12)), (1.0, arc_r.clone())],
                        abs: false,
                    },
                },
            ),
            (
                "internal circle tangency",
                EqKind::Norm {
                    k: 1.0,
                    v: VecExpr::diff(0, 3),
                    rhs: Rhs {
                        terms: vec![(1.0, RTerm::Slot(12)), (-1.0, RTerm::Slot(13))],
                        abs: true,
                    },
                },
            ),
            (
                "point on line",
                EqKind::SignedDist {
                    dir: d01.clone(),
                    v: VecExpr::diff(2, 0),
                    abs: false,
                    rhs: Rhs::zero(),
                },
            ),
            (
                "point-line distance",
                EqKind::SignedDist {
                    dir: d01.clone(),
                    v: VecExpr::diff(2, 0),
                    abs: true,
                    rhs: Rhs::constant(0.2),
                },
            ),
            (
                "line-circle tangency",
                EqKind::SignedDist {
                    dir: d01.clone(),
                    v: VecExpr::diff(4, 0),
                    abs: false,
                    rhs: Rhs::term(-1.0, arc_r.clone()),
                },
            ),
            (
                "symmetric midpoint",
                EqKind::SignedDist {
                    dir: d01.clone(),
                    v: mid,
                    abs: false,
                    rhs: Rhs::zero(),
                },
            ),
            (
                "symmetric projection",
                EqKind::Proj {
                    dir: d01.clone(),
                    v: VecExpr::diff(3, 2),
                },
            ),
            (
                "parallel",
                EqKind::DirPair {
                    a: d01.clone(),
                    b: d23.clone(),
                    mode: PairMode::Cross,
                },
            ),
            (
                "perpendicular",
                EqKind::DirPair {
                    a: d01.clone(),
                    b: d23.clone(),
                    mode: PairMode::Dot,
                },
            ),
            (
                "angle",
                EqKind::DirPair {
                    a: d01.clone(),
                    b: d23.clone(),
                    mode: PairMode::AngleSin { value: 0.7 },
                },
            ),
            (
                "tangent at point, spline end vs line",
                EqKind::DirPair {
                    a: tangent_end,
                    b: d23,
                    mode: PairMode::Cross,
                },
            ),
            (
                "tangent at point, line vs arc normal",
                EqKind::DirPair {
                    a: d01,
                    b: VecExpr::diff(3, 5),
                    mode: PairMode::Dot,
                },
            ),
            (
                "point on control spline, x",
                EqKind::SplinePoint {
                    spline: 0,
                    t: 14,
                    p: 0,
                    axis: 0,
                },
            ),
            (
                "point on interpolated spline, y",
                EqKind::SplinePoint {
                    spline: 1,
                    t: 14,
                    p: 0,
                    axis: 1,
                },
            ),
            (
                "point on quadratic interpolated spline, x",
                EqKind::SplinePoint {
                    spline: 2,
                    t: 14,
                    p: 5,
                    axis: 0,
                },
            ),
            (
                "point on control spline with interior knots, y",
                EqKind::SplinePoint {
                    spline: 3,
                    t: 14,
                    p: 0,
                    axis: 1,
                },
            ),
            (
                "point on interpolated spline with interior knots, x",
                EqKind::SplinePoint {
                    spline: 4,
                    t: 14,
                    p: 2,
                    axis: 0,
                },
            ),
        ]
    }

    fn check_gradient(name: &str, eq: &EqKind, v: &[f64]) -> Result<(), TestCaseError> {
        // Directions shorter than this are skipped: the test is about the
        // formulas, not about conditioning near degenerate geometry.
        const DEGEN: f64 = 0.05;
        let sp = splines();
        let mut g = Vec::new();
        let Some(e) = eq.eval(v, &sp, DEGEN, Some(&mut g)) else {
            return Ok(());
        };
        prop_assert!(e.r.is_finite());
        let mut dense = vec![0.0; SLOTS];
        for (i, d) in g {
            dense[i] += d;
        }
        for (i, &analytic) in dense.iter().enumerate() {
            let h = 1e-6;
            let mut vp = v.to_vec();
            let mut vm = v.to_vec();
            vp[i] += h;
            vm[i] -= h;
            let (Some(rp), Some(rm)) = (
                eq.eval(&vp, &sp, DEGEN, None),
                eq.eval(&vm, &sp, DEGEN, None),
            ) else {
                continue;
            };
            let fd = (rp.r - rm.r) / (2.0 * h);
            let scale = 1.0 + analytic.abs().max(fd.abs());
            // An absolute value whose argument changes sign within ±h makes
            // the one-sided differences disagree: not a smooth point.
            let fwd = (rp.r - e.r) / h;
            let bwd = (e.r - rm.r) / h;
            if (fwd - bwd).abs() > 0.1 * scale {
                continue;
            }
            prop_assert!(
                (fd - analytic).abs() <= 2e-5 * scale,
                "{name}: slot {i}: analytic {analytic}, central difference {fd}"
            );
        }
        // Slots reported as referenced cover every nonzero derivative.
        let refs = eq.slots(&sp);
        for (i, d) in dense.iter().enumerate() {
            if *d != 0.0 {
                prop_assert!(refs.contains(&i), "{name}: slot {i} not listed");
            }
        }
        Ok(())
    }

    fn config() -> impl Strategy<Value = Vec<f64>> {
        (
            prop::collection::vec(-1.0f64..1.0, 12),
            0.2f64..1.0,
            0.2f64..1.0,
            0.02f64..0.98,
        )
            .prop_map(|(mut v, r1, r2, t)| {
                v.push(r1);
                v.push(r2);
                v.push(t);
                v
            })
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(512))]

        #[test]
        fn gradients_match_central_differences(v in config()) {
            for (name, eq) in all_equations() {
                // Skip samples where an absolute value sits at its kink.
                check_gradient(name, &eq, &v)?;
            }
        }
    }

    #[test]
    fn degenerate_directions_are_reported() {
        let v = vec![0.0; SLOTS];
        let sp = splines();
        for (name, eq) in all_equations() {
            let divides = matches!(
                eq,
                EqKind::SignedDist { .. } | EqKind::Proj { .. } | EqKind::DirPair { .. }
            );
            let out = eq.eval(&v, &sp, 1e-9, None);
            assert_eq!(out.is_none(), divides, "{name}");
            if let Some(e) = out {
                assert!(e.r.is_finite(), "{name}");
            }
        }
    }

    #[test]
    fn angle_branch_helpers() {
        // a = (1, 0), b = (0, 1): θ = π/2.
        let mut v = vec![0.0; SLOTS];
        v[2] = 1.0; // P1 = (1, 0)
        v[7] = 1.0; // P3 = (0, 1)
        let eq = EqKind::DirPair {
            a: VecExpr::diff(1, 0),
            b: VecExpr::diff(3, 2),
            mode: PairMode::AngleSin {
                value: std::f64::consts::FRAC_PI_2,
            },
        };
        assert!((eq.angle_theta(&v).unwrap() - std::f64::consts::FRAC_PI_2).abs() < 1e-15);
        assert!((eq.angle_cos(&v, 1e-9).unwrap() - 1.0).abs() < 1e-15);
        assert_eq!(signed_angle(0.0, -1.0), std::f64::consts::PI);
        assert_eq!(signed_angle(-0.0, -1.0), std::f64::consts::PI);
    }
}
